# Numerical specification

This file supplements the numerical section archived in numerical_section.tex. The five archived EPS files are the actual files referenced by that section.

## Experiments and integration

| Experiment | Horizon | Relative tolerance | Absolute tolerance | Maximum step |
| --- | --- | --- | --- | --- |
| Tightened and untightened | 1200 s | 2e-9 | 1e-12 | 0.25 s |
| Refined main runs | 1200 s | 2e-10 | 1e-13 for 60 slow states; 1e-12 for 80 fast states | 0.125 s |
| Nominal HB | 1200 s | 2e-9 | 1e-12 | 0.25 s |
| Refined nominal HB | 1200 s | 2e-10 | 1e-13 | 0.125 s |
| Independent local | 200 s | 2e-10 | 1e-14 | 0.1 s |
| Refined local | 200 s | 2e-11 | 1e-15 | 0.05 s |

All use ode15s. Main saved outputs are spaced by 0.1 s; diagnostics combine those samples, accepted internal nodes and 0.01 s dense output. Communication-support events are localized using signed-distance crossings. Tracking and estimation peaks are refined on the saved numerical interpolants.

Fig. 2 displays 0–400 s. The reported main-run extrema, error budgets and reference-QP checks cover the complete 0–1200 s horizon.

The local output grid includes logarithmically spaced initial times to resolve the fast layer. Its slow error is the norm of col{y_est−1⊗y*,v}; the fast error is the norm of col{u−u_nom,zeta_0,lambda_0}, with the retained local rows used in the Appendix. The full projected flow is integrated in equilibrium-error coordinates to resolve the small initial perturbation. Fig. 3 displays slow-error decay and the fast transient above the numerical floor; the subsequent fast-error plateau near 10⁻¹⁵ reflects numerical precision. Its caption states this interpretation, and the inset retains t/epsilon as its time axis.

## Geometry, graph and physical residuals

Collision and obstacle barriers are squared-distance margins. The connectivity certificate is

h_c(y) = 5 det(L_chi^(r)(y))/512 − 0.04.

Communication support is defined by distance < d_c. The graph diagnostics record support crossings, component count, algebraic connectivity and the maximum-bottleneck spanning-tree weight. Near a support boundary, the minimum active weight can approach zero while a spanning tree remains well weighted.

The original physical ECBF residual is q = −h_ddot − K2 h_dot − K1 h; its safe direction is q≤0. Collision residuals apply to active rows, while geometric collision margins are evaluated for every pair.

The declared physical domain is the position box [−3.1,3.1]×[−1.3,1.3], intersected with pairwise diameter ≤2.5 and individual speed ≤0.5. It gives coefficient bounds (a_int,a_obs,a_c) = (5,10.7689182372,2.69279673023).

Collision/obstacle physical rows and margins use scale 0.06; connectivity rows and margins use scale 0.1. Mismatch terms are unscaled and cancel upon aggregation. The connectivity discrepancy on the declared input ball ||u||≤15 is |db|+15||da||; actual-input discrepancy is recorded separately.

### Scaled local rows

The local collision row in the manuscript is

$$g_{ij,\mathrm{sc}}^i=\chi_{ij}(y)\left[0.06\left(q_{ij,i}^{\mathrm{int}}+\frac12\Delta_{ij}^{\mathrm{int}}\right)+\zeta_{ij}^i-\zeta_{ij}^j\right]\leq0.$$

In [audit_core.m](../matlab/audit_core.m), the input coefficients and physical offset include `P.barrierScale=0.06`, while the two mismatch coefficients are `+chi_ij` and `−chi_ij`. The opposite local row reverses the mismatch difference, so summing the two rows cancels mismatch and gives `0.06 chi_ij (q_ij^int+Delta_ij^int)`. Obstacle rows have scale 0.06 and no pairwise mismatch difference; connectivity rows use scale 0.1 with their original mismatch allocation.

Scaling physical terms while leaving mismatch unchanged specifies the adopted regularized problem and primal–dual dynamics at fixed xi. It is not a uniform multiplication of the whole local inequality and need not preserve the regularized optimizer relative to a different physical scaling.

The offline reference starts from the same adopted constraints. [normalized_reference_rows.m](../matlab/normalized_reference_rows.m) cancels a strictly positive chi_ij from an active collision row before the solve, leaving the same physical scale 0.06 and mismatch coefficients ±1. Support is tested by distance < d_c; an exactly inactive pair is kept as a zero row. Subsequent positive whole-row conditioning in [audit_qp.m](../matlab/audit_qp.m) preserves this QP and maps multipliers to the stated reference-row coordinates. This conditioning is separate from the adopted physical-term scaling.

## QP references and KKT residuals

Offline references use quadprog with interior-point-convex, ConstraintTolerance=1e-11, OptimalityTolerance=1e-11 and MaxIterations=2000. Fixed-scale retries and verified active-set polishing solve the same QP. Positive collision-row weights are analytically cancelled before the reference solve; exact zero rows are retained for inactive pairs.

For min 0.5 wᵀHw+fᵀw subject to Aw+b≤0, the saved residuals are:

- Primal: max(0,max(Aw+b)).
- Stationarity: ||Hw+f+Aᵀlambda||∞.
- Complementarity: ||lambda.*(Aw+b)||∞.
- Dual: max(0,−min(lambda)).

Raw signed and normalized primal residuals are retained in the full data.

## Theorem checks

Theorem 1 gives tilde_l≈1.12161, c>185.085144 and gamma>56.689722. The selected values are c=200 and gamma=65.

Table I uses the connectivity minimum of the nominal flow, computed by `reproduce_paper('nominal')`. The independent RHS implements nominal_HB_stacked and is checked against the shared model equations at three deterministic states. It integrates the 50 decision estimates and 10 velocities, with instantaneous u=u_nom and the original initial estimates and velocities. The communication weights are evaluated from the agents' physical positions through audit_core('graph'). Diagnostics combine accepted solver nodes, 0.1 s outputs and 0.01 s dense output, with interpolation refinement for an interior sampled minimum and an event at lambda2=0.04. Both precision settings return a minimum of 1.58616280459112 at t=0. This is a finite-horizon numerical check of the nominal trajectory; the separate tightened-trajectory minimum remains in its own saved safety diagnostics.

The nominal command produces verification/nominal_hb.json and, under results/nominal and results/nominal_refined, summary.json, lambda2.csv and trajectory.mat. The MAT files include the numerical interpolants; these generated files stay in the run directory. The compact combined JSON is included in the public repository.

For Theorem 2, the three tightening design expressions are approximately 0.431335, 0.656904 and 0.397298, below 0.5, 0.7 and 0.4. Dense tightened-trajectory peaks are approximately 0.0571577 for tracking and 0.00949398 for connectivity discrepancy.

The applied-input peak at audited tightened-trajectory states is 13.4800990276 (13.48010 rounded), below the declared input-ball radius 15. The separate reference optimizer's input peak is 13.4802921635. These are trajectory diagnostics.

The regional uniform connectivity-discrepancy bound, full slow-domain residence, region-wide QP feasibility and continuous-time uniform tracking bound are distinct retained hypotheses. The first takes a supremum over Theta_R × U_R; the reported input-ball discrepancy instead maximizes over U_R with each trajectory state fixed. Physical residence in D_R does not certify forward invariance of the full estimate-containing slow domain Theta_R. Likewise, the 138745/138745 valid pointwise reference solves and their maximum KKT component of approximately 1.86e-10 do not certify regional feasibility or dynamic tracking accuracy.

The [Theorem 2 condition note](theorem2_conditions.md), verification/theorem2_conditions_with_scope.csv and verification/theorem2_audit.json provide these checks separately. `reproduce_paper('verify')` recomputes the report from model parameters and archived diagnostics and checks it against the committed record. The report explicitly retains regionalUniformCertificate=false and continuousTimeProof=false.

Theorem 3 uses explicitly constructed local radii from the Appendix. The selected epsilon is 5e-6 < 9.31743406e-6; initial slow/fast errors are approximately 4.45854e-9 and 3.29491e-7, below radii 5.94472e-9 and 6.58982e-7. These radii are conservative sufficient neighborhoods, not estimates of the largest region of attraction. All multipliers start at zero, including the two rows of locally omitted pair (3,4). The construction and sufficient inequalities are evaluated in double precision; the derivation is in theorem3_local_certificate.md.

The local slow perturbation is the normalized deterministic vector sin(k sqrt(2)), k=1,...,60, with norm 0.75 r_Theta1. The control starts at u_nom; the first retained mismatch coordinate is perturbed by 0.5 r_f, and all other fast mismatch coordinates and multipliers start at zero. Thus the nonzero initial fast error resolves the initial layer while remaining inside the checked local budget.

The stated slow/fast norm estimate concerns the retained errors e_vartheta and e_f. With zero omitted-row multipliers, the omitted mismatch coordinates satisfy `zeta_om(t)=exp(−xi t/epsilon) zeta_om(0)`. Their norm bound therefore includes `||zeta_om(0)||`; these coordinates are kept in the full 140-state ODE. In the archived local experiment their initial value is zero.

## Reproduction and archive

The quick redraw (`reproduce_paper`) reads the compact plot data and fills the final editable layout templates. Full mode (`reproduce_paper('full')`) reintegrates nominal HB and the six standard/refined filtered or local experiments, then regenerates dense diagnostics, local initializations and data-derived plots in a new run folder. The `nominal` mode independently integrates nominal HB; `verify` checks saved records, local inequalities and figure data and regenerates the scoped condition reports. The paper EPS/FIG files remain the source archive; regenerated exports can differ in file metadata and renderer encoding.

The v1.0.2 code-and-figures package includes the compact plotting data. Six unchanged full-data archives remain attached to the v1.0.0 reference-data release: selected, selected_refined, untightened, untightened_refined, local and local_refined. Each includes trajectory.mat, diagnostic tables and per-run summaries. The v1.0.2 file data/reference_archives.json binds their fixed v1.0.0 URLs, archive checksums and per-file checksums.
