# Theorem 2: condition checks and evidence scope

This note corresponds to comment 03 of the revision memo. Theorem 2 transfers specified tracking and connectivity-estimation error budgets into sufficient tightening margins for the original ECBF constraints. Its regional hypotheses remain explicit conditions of the safety result.

Run `reproduce_paper('verify')` to regenerate the four-column tables and machine-readable report in the run's `verification/` folder. The committed copies are [theorem2_conditions_with_scope.csv](../verification/theorem2_conditions_with_scope.csv), [theorem2_audit.json](../verification/theorem2_audit.json), and the combined [theorem_checks_with_scope.csv](../verification/theorem_checks_with_scope.csv). The generator is [theorem2_condition_audit.m](../matlab/theorem2_condition_audit.m); full simulation mode also calls it when building the tables.

## Design inequalities

The declared physical domain D_R uses the position box [-3.1,3.1] × [-1.3,1.3], maximum pairwise distance 2.5, and maximum individual speed 0.5. The coefficient formulas give (a_int,a_obs,a_c) = (5,10.7689182372,2.69279673023). With r_u=0.061, delta_c=0.03 and tightening (0.5,0.7,0.4), the sufficient inequalities are:

| Design expression | Required bound | Value | Remaining margin |
| --- | --- | --- | --- |
| sqrt(2) a_int r_u | ≤0.5 | 0.431335136524 | 0.0686648634762 |
| a_obs r_u | ≤0.7 | 0.656904012471 | 0.0430959875294 |
| sqrt(5) a_c r_u + delta_c | ≤0.4 | 0.397297868842 | 0.00270213115808 |

These are analytic coefficient bounds and sufficient inequalities evaluated in double precision. The boundary-reactivation condition gives 4.4 ≤5.49624563 from the declared relative-speed bound 1.

## Tightened-trajectory observations

Diagnostics use accepted solver nodes, dense output and interpolation-based peak refinement over 0–1200 s. The standard and refined trajectories are recorded separately in the JSON report.

| Quantity | Required budget or check | Standard-run observation | Verification |
| --- | --- | --- | --- |
| Tracking error | r_u=0.061 | 0.0571576535274 | Mesh/dense/interpolant diagnostics along the tightened trajectory |
| Input-ball connectivity discrepancy | delta_c=0.03 | 0.00949398359918 | Maximize over the input ball at each audited trajectory state |
| Applied input norm | ‖u‖≤15 | 13.4800990276 (13.48010 rounded) | Maximum at audited trajectory states |
| Physical-domain residence | Position, diameter and speed bounds of D_R | Satisfied at checked states | Physical-state trajectory check |
| Reference-QP validity | All audited states valid | 138745/138745 | Offline pointwise numerical solves and KKT checks |

At a fixed trajectory state, the discrepancy is

`delta_c_tr(t) = max_{||u||<=15} |q_c(vartheta(t),u) - sum_i qhat_i_c(vartheta(t),u_i)| = |db| + 15 ||da||_2`.

The maximization is over inputs with the trajectory state fixed. It is different from maximizing over every state and input in Theta_R × U_R. The actual-input discrepancy is recorded separately and is smaller in this run; it is not used in place of the input-ball value.

The applied-input peak and the reference optimizer's input peak are also separate quantities. The latter is 13.4802921635 at audited states. Both recorded peaks are below 15. These observations do not certify input membership at all continuous times.

## Retained hypotheses

The following conditions have distinct quantifiers and are separate rows in the public table:

| Condition | Required statement | Current evidence status |
| --- | --- | --- |
| Regional uniform connectivity discrepancy | sup over Theta_R × U_R of the discrepancy ≤delta_c | Retained assumption; trajectory-state input-ball maxima do not establish it |
| Full slow-domain residence | vartheta(t)∈Theta_R on the maximal existence interval | Retained assumption; forward invariance of Theta_R has not been established |
| Region-wide QP feasibility | The decomposed regularized QP is feasible for every vartheta∈Theta_R | Retained assumption; pointwise trajectory solves do not establish it |
| Continuous-time uniform tracking bound | ‖u(t)−uhat_xi_star(vartheta(t))‖≤r_u on that interval | Retained assumption; recorded peaks do not establish the all-time bound |

D_R concerns physical positions and velocities. Theta_R also contains the decision-estimate components of the slow state. Observed physical residence in D_R does not prove full slow-state residence or forward invariance in Theta_R. Input membership u(t)∈U_R is another all-time premise; its recorded peak is classified as a trajectory diagnostic.

## KKT residuals

For the tightened standard run, the four recorded maximum reference-QP residuals are approximately:

| Coordinate | Maximum |
| --- | --- |
| Primal | 6.40349e-13 |
| Stationarity | 6.78521e-11 |
| Complementarity | 1.86351e-10 |
| Dual | 2.16416e-16 |

Their largest component is 1.86351e-10. They measure the quality of numerical reference solves in the documented row coordinates. The tracking error instead compares the dynamic filter input with the instantaneous regularized-QP optimizer. A small KKT residual establishes neither a small dynamic tracking error nor QP feasibility throughout Theta_R. Residual definitions and row normalization are given in [numerical_details.md](numerical_details.md).

The report therefore keeps `regionalUniformCertificate=false` and `continuousTimeProof=false` even when every design and recorded-trajectory check passes.
