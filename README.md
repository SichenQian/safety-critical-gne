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

## Citation

If you use this code in your research, please cite the accompanying paper:

```bibtex
@unpublished{qian_safety_critical_gne,
  author = {Qian, Sichen and Liu, Hongzhe and Yu, Wenwu},
  title  = {Distributed Safety-Critical Seeking of an Interior Generalized {Nash} Equilibrium over State-Dependent Graphs},
  note   = {Manuscript prepared for submission to {IEEE} Transactions on Control of Network Systems}
}
```
