# Safety-critical GNE seeking over state-dependent graphs

MATLAB code, data and editable figures for the manuscript **Distributed Safety-Critical Seeking of an Interior Generalized Nash Equilibrium over State-Dependent Graphs**, by Sichen Qian, Hongzhe Liu and Wenwu Yu.

This version reproduces the numerical section of the current manuscript: a five-agent planar game, a tightened-versus-untightened comparison with three obstacles, genuine communication-edge switching, and an independent local convergence experiment.

The [v1.0.1 package](https://github.com/SichenQian/safety-critical-gne/releases/tag/v1.0.1) is synchronized with the manuscript revised on 4 November 2026. [Manuscript alignment](docs/manuscript_alignment.md) records the numerical explanations and their implementation sources.

## Run

Use MATLAB R2024a. Optimization Toolbox is required for safety-filter simulations and QP diagnostics.

From the repository folder:

    reproduce_paper

This redraws all five paper figures from the included compact numerical data. Outputs are written to a new folder under runs/.

    reproduce_paper('full')     % Recompute nominal HB, six filtered/local runs and checks
    reproduce_paper('verify')   % Check parameters, initialization, saved diagnostics and figure data
    reproduce_paper('nominal')  % Recompute Table I's nominal HB connectivity; no figures

An output folder can be supplied as the second argument:

    reproduce_paper('full', fullfile(pwd, 'runs', 'my_experiment'))

Successful execution ends with REPRODUCTION_PASS. The full run simulates the 140-state projected primal–dual ODE. QP solutions provide initial conditions and offline reference diagnostics.

## Paper figures

The files below are the exact EPS figures referenced by the manuscript's includegraphics commands. Each has a matching editable FIG, vector PDF and PNG preview in [figures/](figures).

| Manuscript figure | File stem | Content |
| --- | --- | --- |
| Fig. 1 | trajectory | Tightened/untightened trajectories and geometric boundary views |
| Fig. 2(a) | safety | Collision and obstacle margins |
| Fig. 2(b) | connectivity | Matrix-tree certificate and algebraic connectivity |
| Fig. 2(c) | edge | Edge activation states and total active-edge count |
| Fig. 3 | exp | Slow/fast local error norms and the fast initial layer |

The final layouts are stored in presentation/. The plotting code fills those layouts from numerical data, preserving the manuscript's fonts, colors, labels and inset positions. The source correspondence is recorded in [docs/figure_sources.json](docs/figure_sources.json).

![Trajectories](figures/trajectory.png)
![Local convergence](figures/exp.png)

## Model and parameters

The agents have second-order dynamics and quadratic costs

J_i = ||y_i-p_i||²/2 + 0.03 Σ(j≠i)||y_i-y_j||²/2,

with B = 1.15 I − 0.03 11ᵀ and p = y*B. All agents start at rest with exact initial decision estimates. Initial and target positions, obstacle geometry and every gain are defined in [matlab/design_parameters.m](matlab/design_parameters.m) and [matlab/audit_core.m](matlab/audit_core.m).

| Parameter | Value |
| --- | --- |
| Communication radius / collision distance | 2.2 / 0.345 |
| Obstacle centers | (0.24,0), (1.72,1.1), (1.72,−1.1) |
| Obstacle radii | 0.625, 0.315, 0.315 |
| c / gamma / xi | 200 / 65 / 1.4 |
| epsilon / nu | 5×10⁻⁶ / 0.02 |
| K1 / K2 | 4 / 4.6 |
| Collision / obstacle / connectivity tightening | 0.5 / 0.7 / 0.4 |
| Tracking / connectivity-discrepancy budgets | 0.061 / 0.03 |
| Algebraic-connectivity lower bound | 0.04 |

The untightened run sets only the three tightening margins to zero and uses the same complete initial state. The tightened graph has 9→10→9 active edges: (2,5) activates at approximately 5.59002 s, and (3,4) deactivates at approximately 313.17490 s.

The local experiment uses initial conditions checked against the Appendix's sufficient inequalities, with u(0)=u_nom, a mismatch perturbation of norm 0.5 r_f, and all multipliers (including the locally omitted pair (3,4)) initialized at zero. The constructed radius r_Theta1≈5.94472×10⁻⁹ is a conservative sufficient neighborhood; it does not estimate the largest region of attraction.

Fig. 3 shows slow-error decay and the fast initial transient on t/epsilon. The subsequent fast-error plateau, near 10⁻¹⁵, reflects numerical precision. The retained error estimate and the separate decay zeta_om(t)=exp(−xi t/epsilon) zeta_om(0) follow the statement and proof of Theorem 3; [the local-certificate note](docs/theorem3_local_certificate.md) defines both sets of coordinates.

## Numerical data and checks

### Theorem 2 condition checks

The condition table distinguishes design arithmetic, observations along the tightened trajectory, and retained hypotheses. Table I states the respective nominal/tightened trajectory sources and the 0–1200 s window with numerical peak refinement. Its note retains regional residence, QP feasibility, input boundedness, and uniform estimation and tracking-error bounds as premises. The current manuscript hashes are recorded in [figure_sources.json](docs/figure_sources.json).

| Condition or quantity | Required bound | Numerical value | Verification |
| --- | --- | --- | --- |
| sqrt(2) a_int r_u | ≤0.5 | 0.431335 | Analytic coefficient bound evaluated in double precision |
| a_obs r_u | ≤0.7 | 0.656904 | Same |
| sqrt(5) a_c r_u + delta_c | ≤0.4 | 0.397298 | Same |
| Tightened-trajectory tracking peak | <r_u=0.061 | 0.0571577 | Mesh, dense output and interpolant diagnostics |
| Input-ball connectivity-discrepancy peak | <delta_c=0.03 | 0.00949398 | At each trajectory state maximize over ‖u‖≤15 |
| Regional uniform connectivity discrepancy | sup over Theta_R × U_R ≤delta_c | Retained assumption | Trajectory-state maxima do not establish this regional bound |
| Input-set membership | u(t)∈U_R, with ‖u‖≤15 | 13.48010 | Input-norm peak at audited trajectory states |
| Full slow-domain residence | vartheta(t)∈Theta_R | Retained assumption | Forward invariance of Theta_R not established |
| Region-wide QP feasibility | Feasible for every vartheta∈Theta_R | Retained assumption | Pointwise trajectory checks only |
| Continuous-time uniform tracking bound | ‖u(t)−uhat_xi_star(vartheta(t))‖≤r_u | Retained assumption | Numerical peaks do not establish the all-time bound |

The three remaining tightening margins are approximately 0.0686649, 0.0430960 and 0.00270213. Physical-state residence in D_R is checked separately from full slow-state residence in Theta_R. All 138745 audited tightened-trajectory states have valid reference QPs; the largest recorded KKT component is about 1.86×10⁻¹⁰. This residual measures reference-solver accuracy, separately from dynamic filter tracking error.

The regional discrepancy budget and its Lipschitz sufficient bound are distinguished in [the condition note](docs/theorem2_conditions.md#regional-discrepancy-budget): if E_R≤B_R, choosing delta_c≥B_R is sufficient. The numerical value delta_c=0.03 retains its stated regional-bound premise.

Fig. 2 displays the first 400 s. Error peaks, physical-domain checks and reference-QP diagnostics use the complete 0–1200 s runs.

`reproduce_paper('verify')` regenerates the scoped condition tables and JSON report in the run's verification/ folder. The public records are [theorem2_conditions_with_scope.csv](verification/theorem2_conditions_with_scope.csv), [theorem2_audit.json](verification/theorem2_audit.json) and the combined [theorem_checks_with_scope.csv](verification/theorem_checks_with_scope.csv). [Theorem 2 details](docs/theorem2_conditions.md) give the definitions, source quantities and verification scope.

### Nominal flow and reference data

Table I's nominal-flow connectivity minimum is reproduced independently by [matlab/run_nominal_hb.m](matlab/run_nominal_hb.m). It integrates the 60-state nominal HB equations with the same game, initial slow state and communication weights. Standard and refined runs give a minimum of 1.58616280459112 over 0–1200 s, reported as 1.58616 in the table. Parameters, solver settings and numerical checks are saved in [verification/nominal_hb.json](verification/nominal_hb.json).

The current code and scoped condition checks are fixed in [v1.0.1](https://github.com/SichenQian/safety-critical-gne/releases/tag/v1.0.1). The included data/paper_plot_data.mat is sufficient to redraw the five figures. The unchanged full trajectories, solver interpolants, dense diagnostic tables and refined-precision results remain in the [v1.0.0 reference-data release](https://github.com/SichenQian/safety-critical-gne/releases/tag/v1.0.0).

Download and verify all six full-data archives from MATLAB:

    reproduce_paper('download', fullfile(pwd, 'runs', 'reference_data'))
    reproduce_paper('verify', fullfile(pwd, 'runs', 'reference_data'))

The download checks archive and extracted-file SHA-256 values from [data/reference_archives.json](data/reference_archives.json). The same data can be generated directly with the full mode.

All experiments use ode15s. The main runs cover 1200 s with RelTol=2×10⁻⁹, AbsTol=10⁻¹² and MaxStep=0.25 s. The independent local run covers 200 s with 2×10⁻¹⁰, 10⁻¹⁴ and 0.1 s. Refined repeats use tighter tolerances and smaller maximum steps.

The checks cover gain inequalities, tightening arithmetic, local initialization, physical barriers, ECBF residuals, graph connectivity, QP/KKT residuals and precision repeats. [docs/numerical_details.md](docs/numerical_details.md) gives the solver settings, residual definitions, the scaled local constraint used by the controller and reference QP, and verification scope. [verification/theorem_checks_with_scope.csv](verification/theorem_checks_with_scope.csv) separates analytic design checks from trajectory diagnostics.

## Repository layout

| Folder | Contents |
| --- | --- |
| matlab/ | Dynamics, safety filter, certificates, diagnostics and paper plotting |
| data/ | Compact plotting data and full-data archive checksums |
| figures/ | Actual manuscript EPS files and matching FIG/PDF/PNG |
| presentation/ | Editable final-layout templates |
| reference/ | Per-experiment numerical summaries and initialization records |
| verification/ | Numerical tables, precision comparisons and validation logs |
| docs/ | Numerical specification, figure provenance and local certificate derivation |

This repository contains the current manuscript package. The entry point is reproduce_paper.m at the repository root.

## Citation

    @unpublished{qian_safety_critical_gne,
      author = {Sichen Qian and Hongzhe Liu and Wenwu Yu},
      title = {Distributed Safety-Critical Seeking of an Interior
               Generalized Nash Equilibrium over State-Dependent Graphs},
      year = {2026},
      note = {Manuscript},
      url = {https://github.com/SichenQian/safety-critical-gne}
    }

Code is distributed under the [MIT License](LICENSE).
