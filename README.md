# Distributed Safety-Critical Seeking of an Interior Generalized Nash Equilibrium over State-Dependent Graphs

MATLAB implementation and numerical reproducibility package for the paper:

**Distributed Safety-Critical Seeking of an Interior Generalized Nash Equilibrium over State-Dependent Graphs**  
Sichen Qian, Hongzhe Liu, and Wenwu Yu.

**Paper status:** Prepared for submission to *IEEE Transactions on Control of Network Systems* (IEEE TCNS).

## Numerical example

| MATLAB file | Example | Figures in the manuscript |
| --- | --- | --- |
| `reproduce_paper.m` | Five-agent second-order game with collision avoidance, three circular obstacles, and state-dependent communication | Figs. 1–3 |

The example implements a distributed heavy-ball equilibrium seeker with an
exponential control barrier function (ECBF) safety filter realized by projected
primal–dual dynamics. It compares the proposed tightened flow with an untightened
ablation that sets the collision, obstacle, and connectivity tightening margins
to zero while retaining the other algorithm parameters.

The MATLAB file is self-contained. The game, obstacle geometry, communication
weights, equilibrium target, controller parameters, and deterministic initial
conditions are defined in the file. Both variants start from the same positions,
zero velocities, and initial decision estimates; their fast variables are
initialized at the respective initial primal–dual QP solutions.

Each variant is integrated continuously over 0–1200 seconds. All proposed-method
panels use the same tightened trajectory, including the local error plot over
850–1050 seconds. Parameter, derivative, and initial-condition checks run before
integration, followed by sampled feasibility, tracking, and safety diagnostics.

## Requirements

- MATLAB; validated with **R2024a**.
- **Optimization Toolbox**, for `quadprog`.
- No external input data, saved trajectories, or additional user-supplied MATLAB files are required.

`quadprog` is used for fast-state initialization and offline diagnostics. The
online controller is simulated through its projected primal–dual differential
equations. Integration uses `ode15s` with `RelTol = 2e-8`, `AbsTol = 1e-11`,
and `MaxStep = 0.25`; trajectories and diagnostic quantities are sampled every
0.1 seconds.

## Usage

Open `reproduce_paper.m` in the MATLAB Editor and click **Run**, or execute the
following command from the repository folder:

```matlab
result = reproduce_paper;
```

The default outputs are saved relative to the MATLAB file's own location:

```text
results/
```

An alternative output folder can be selected by passing a directory:

```matlab
result = reproduce_paper(fullfile(pwd, 'my_results'));
```

Each invocation recomputes both trajectories and overwrites the corresponding
files in the selected output folder. Progress is printed during integration
and offline diagnostics. A successful run ends with `REPRODUCTION_PASS`.

## Outputs

The program saves seven figures in `.fig`, `.png`, `.pdf`, and `.eps` formats,
together with simulation data and numerical diagnostics:

- `results/experiment_data.mat`: full trajectories, initial states, model parameters, and diagnostic time series for both variants.
- `results/parameters.mat`: model parameters and computed analytical bounds.
- `results/plot_data.mat`: trajectory geometry, safety margins, and local error data used for plotting.
- `results/tightened_margins.csv` and `results/untightened_margins.csv`: collision and obstacle margins, the matrix-tree connectivity certificate, and the physical graph's algebraic connectivity.
- `results/same_run_error_tail.csv`: physical time, elapsed time, and aggregate, slow, and fast error norms over the local interval.
- `results/numerical_audit.json` and `results/figure_data_audit.json`: controller settings, sampled numerical checks, and figure correspondence.
- `results/figure_editability_check.json` and `results/runtime_validation.json`: native FIG edit-save-reopen checks, software versions, runtime, and completion status.

The program stops with an assertion error if a required check fails. Safety
diagnostics are evaluated on the saved samples; the analytical certificate
checks and their scope are recorded separately in the audit files. The
matrix-tree quantity is a sufficient connectivity certificate, and its negative
excursions in the untightened run are reported alongside the physical graph's
algebraic connectivity.

## Figure correspondence

All figures are exported to `results/` in `.fig`, `.png`, `.pdf`, and `.eps` formats:

| Figure | File stem | Content |
| --- | --- | --- |
| Fig. 1(a) | `my_trajectory` | Tightened planar trajectories and obstacle-clearance inset |
| Fig. 1(b) | `my_trajectory_2` | Untightened planar trajectories and obstacle-penetration inset |
| Fig. 2(a) | `inter_agent_avoidance` | Minimum inter-agent collision margin |
| Fig. 2(b) | `static_obstacle_avoidance` | Minimum obstacle-avoidance margin |
| Fig. 2(c) | `connectivity_preservation` | Matrix-tree connectivity certificate |
| Fig. 3(a) | `HB_flow_hebing` | Position-component convergence under the tightened flow |
| Fig. 3(b) | `exp_2` | Aggregate closed-loop error over 850–1050 seconds along the same tightened trajectory |

The manuscript uses `my_trajectory_1.eps` for Fig. 1(a); the program exports this
panel as `my_trajectory.eps`. Rename or copy this file when inserting the
generated figure into the manuscript. The horizontal axis in Fig. 3(b) shows
elapsed time from 850 seconds, so the plotted range is 0–200 seconds.

## Reproducibility checks

The standalone file was run from scratch with MATLAB R2024a and Optimization
Toolbox 24.1. Both 12001-sample, 140-component trajectories matched the stored
reference trajectories exactly in the validation run. All seven figures were
generated in four formats and passed native FIG edit-save-reopen checks.

`source_manifest.json` records the code's source hashes, and
`validation_summary.json` records the release validation. These files document
the release and are optional for running the MATLAB program.

## Full-precision verification and fixed version

The authoritative repository is
[`SichenQian/safety-critical-gne`](https://github.com/SichenQian/safety-critical-gne).
The annotated tag
[`TCNS-reproducibility-v1`](https://github.com/SichenQian/safety-critical-gne/releases/tag/TCNS-reproducibility-v1)
fixes the verified source and reference outputs:

```sh
git checkout TCNS-reproducibility-v1
git rev-parse HEAD
```

From MATLAB in that checkout, run:

```matlab
summary = verify_reproduction;
figures = verify_figure_data;
```

This entry recomputes both trajectories, runs the original checks, exports all
seven figures, and compares the numeric/logical arrays with the archived
`results/experiment_data.mat`. The default output is `reproduced_results/`,
preserving the archived reference. An alternative folder can be supplied as
`verify_reproduction(fullfile(pwd,'local_run'))`. This entry additionally uses
`compare_numeric_records.m` and `write_full_precision.m`; the original
`reproduce_paper.m` remains self-contained.

| Verification item | Published record in `verification/` |
| --- | --- |
| Actual command, environment, progress and completion | `reproduction_log.txt` |
| Full-precision parameters, gain thresholds and local bounds | `numerical_certificate.txt` |
| Complete 140-component initial states and state order | `initial_states.json` |
| Recomputed arrays versus the archived reference | `reference_comparison.json` |
| Fresh native figure data versus the archived figures | `figure_reference_comparison.json` |
| Actual environment and verification outcome | `verification_summary.json` |
| Manuscript parameters and figure-data correspondence | `paper_correspondence.json`, `paper_correspondence.md` |

Fresh runs write these numerical verification files into `reproduced_results/`.
The certificate uses 17 significant digits and records actual initial states,
selected gains, thresholds, and local-entry quantities. A failed assertion or
reference comparison terminates verification; success ends with
`VERIFICATION_PASS`. Elapsed runtime metadata is explicitly excluded from
array comparisons. The comparison records both its absolute/relative
tolerances and whether all numerical arrays are exactly identical.

`verification_manifest.json` associates the release with SHA-256 hashes of the
inspected manuscript TEX/PDF and figure sources, and Git blob hashes of the
released code and evidence. Manuscript files are not needed to run the code.
`source_manifest.json` retains the original single-file assembly provenance;
the tagged v1 release evidence is in `verification/`. The later tightening
diagnostic update has its own source hashes and evidence in
`verification/tightening/validation.json`.

### Scope of the numerical evidence

Both flows cover 0--1200 seconds at 0.1-second output spacing. They share the
physical initialization and controller settings, except for the stated
tightening margins; fast states start at their respective QP solutions.
The local error panel uses 850--1050 seconds from that same tightened run,
without a state reset or controller switch.

Safety, QP feasibility and tracking budgets are evaluated at saved samples.
Analytical local bounds are recorded separately. The records retain
`regionalUniformCertificate = false` and
`continuousTimeRegionalCertificate = false`: sampled maxima do not establish
a uniform regional certificate, and Corollary 1 is not invoked by this
experiment. A negative matrix-tree margin denotes loss of its sufficient
connectivity certificate; the physical graph's `lambda2` is saved separately.

Figure verification concerns numerical data and panel definitions. Fonts,
colors, crops and export metadata can differ from manuscript assembly; these
presentation differences are excluded from the numerical comparison.

## Theorem 2 tightening diagnostics

`diagnose_tightening_margins.m` adds **trajectory-based numerical diagnostics**
for Theorem 2 in `TCNS_QIAN_v1.pdf` (Section IV-A, p. 6). It reads the tightened
trajectory and produces a Command Window table, a text/JSON summary, and a MAT
record containing the optimizer inputs and diagnostic time series. It does
not integrate the system, change parameters or initial conditions, or create
or modify figures.

To diagnose the archived trajectory:

```matlab
D = diagnose_tightening_margins;
```

The default output is `reproduced_results/tightening_diagnostics/`. To reuse
the QPs from a new reproduction without solving them a second time:

```matlab
R = reproduce_paper(fullfile(pwd,'reproduced_results'));
D = diagnose_tightening_margins(R);
```

A saved MAT file and a separate output folder can also be supplied:

```matlab
D = diagnose_tightening_margins( ...
    fullfile(pwd,'reproduced_results','experiment_data.mat'), ...
    fullfile(pwd,'reproduced_results','tightening_diagnostics'));
```

Fresh `reproduce_paper` runs now retain the primal/dual solutions of the
offline QPs they already solve. For older data without these solutions, the
diagnostic calls the **original** `rows` and `solveScaledQP` functions through
a read-only helper interface in `reproduce_paper`. In either case it checks
the QP constraint and KKT residuals and computes
`norm(u - uStar)` for the regularized decomposable QP (9), including all
mismatch variables, the original `xi = 1.4`, and dissipativity rows. The
nominal input is never substituted for this optimizer.

The table uses physical, unscaled row units, matching the reported
`DeltaInt = 1.2`, `DeltaObs = 1.8`, and `DeltaC = 2.5`. The corresponding row
scales are 0.06, 0.06, and 0.1; the diagnostic also saves the scaled margins
and verifies the pairwise `1/2` and connectivity `1/N` tightening factors.
The matrix-tree barrier retains its existing normalization by 512.

| Theorem 2 quantity/condition | Diagnostic construction | Scope |
| --- | --- | --- |
| `r_u` | Maximum of `norm(u - uStar)` over all saved times | Sampled tracking error, not a continuous-time bound |
| `bar_a_int`, `bar_a_obs` | `2*regionDiameter` and twice the largest obstacle-to-box-corner distance | Analytical coefficient bounds on the existing declared physical region |
| `bar_a_c` | Existing `C.Gh`, bounding the full matrix-tree gradient and hence each local input block | Analytical global derivative bound for the implemented normalized barrier |
| `delta_c` | At each saved state, write `q_c - sum(qhat_c) = da*u + db` and compute `abs(db) + 15*norm(da)` | Exact supremum over the existing input ball at sampled states; not over all of `Theta_R` |
| QP feasibility and original ECBFs | Check the regularized QP/KKT residuals and the actual un-tightened physical ECBF rows at `u(t)` | Saved times only; collision rows are checked where their graph weight is positive |

The actual-input discrepancy is saved separately and is **not** used in place
of the stronger input-ball value. The original dynamic local-row residual is
also reported separately: a positive local-row residual during finite-speed
filter tracking is distinct from a violation of an original physical ECBF.

All three reported margins subtract the independently maximized sampled
tracking and input-ball discrepancy values, using the coefficient bounds
above. Passing these checks supports the margin mechanism along this run.
It does not establish uniform QP feasibility on `Theta_R`, a continuous-time
tracking/containment bound, the discrepancy bound on
`Theta_R x U_R^u`, or the a priori boundary-speed condition in Assumption 2.
Both regional-certificate flags remain `false`.

### Verified numerical values

MATLAB R2024a, 12,001 saved samples on `0:0.1:1200` seconds:

| Quantity (physical row units) | Value | Check |
| --- | ---: | --- |
| Sampled `r_u` | 0.15915799197842689 | From the regularized QP optimizer |
| Sampled-state/input-ball `delta_c` | 0.051054354215450876 | Input radius 15; used in `m_c` |
| Actual-input discrepancy maximum | 0.009725975511371665 | Reported separately |
| `m_int` | 0.074583045920201929 | Nonnegative |
| `m_obs` | 0.086040597583683320 | Nonnegative |
| `m_c` | 1.7900908262239106 | Nonnegative |

The coefficient bounds are `(5, 10.768918237223273, 1.8512977520357721)`.
The scaled margins are `(0.004474982755212116, 0.005162435855020999,
0.1790090826223911)`, so scaling does not change any sign. All 12,001 QPs
passed numerical feasibility/KKT checks. Their maximum constraint residual
was `3.776446132398208e-12` (numerical checking tolerance `1e-8`). The maximum
original collision, obstacle, and connectivity ECBF residuals, with the
convention `q <= 0`, were respectively
`(-1.1846159081172911, -1.7839027819218451, -2.4855712502043819)`.

The modified reproduction was run from scratch: all 3,648,948 original
numeric/logical entries matched the archived reference exactly, with maximum
absolute difference zero. The archived `results/` files are retained.
The new evidence, manuscript hash, and validation commands are recorded in
[`verification/tightening/`](verification/tightening/).

Suggested minimal text for Section V:

> Over the 12,001 saved samples, the regularized-QP tracking error was at most
> 0.15916, and the connectivity discrepancy maximized over the input ball of
> radius 15 at each sampled state was at most 0.05105. The resulting
> tightening margins were 0.07458, 0.08604, and 1.79009; all sampled QPs were
> numerically feasible, and the original ECBF residuals were nonpositive.
> These trajectory-based numerical diagnostics support the tightening
> mechanism but do not establish a uniform regional certificate for Theorem 2.

## Citation

If you use this code in your research, please cite the accompanying paper:

```bibtex
@unpublished{qian_safety_critical_gne,
  author = {Qian, Sichen and Liu, Hongzhe and Yu, Wenwu},
  title  = {Distributed Safety-Critical Seeking of an Interior Generalized {Nash} Equilibrium over State-Dependent Graphs},
  note   = {Manuscript prepared for submission to {IEEE} Transactions on Control of Network Systems}
}
```
