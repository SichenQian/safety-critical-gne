# Distributed Safety-Critical Seeking of an Interior Generalized Nash Equilibrium over State-Dependent Graphs

MATLAB reproducibility package for the TCNS manuscript by Sichen Qian, Hongzhe Liu, and Wenwu Yu.

## Current numerical example

The current case uses the same five-agent second-order game and tightened/untightened comparison, with **natural activation of communication edge (2,5)**. Each variant is one continuous 0–1200 s run. Fig. 1 has two trajectory panels, Fig. 2 three safety panels, and Fig. 3 position convergence plus the 850–1050 s local error segment from the same tightened run. There is no independent local experiment or scheduled adjacency switch.

The earlier case is preserved in Git at [`de46af7`](https://github.com/SichenQian/safety-critical-gne/tree/de46af77d4ecf07194040a6be1ba254b9ceab329). The existing `TCNS-reproducibility-v1` tag identifies the historical release, not the current topology case. See [TOPOLOGY_VALIDATION.md](TOPOLOGY_VALIDATION.md) for the parameter search, affected theoretical conditions, and evidence boundaries.

## Run and verify

Requires **MATLAB R2024a and Optimization Toolbox** (validated environment). No external user files or saved input trajectories are required.

```matlab
R = reproduce_paper;                    % recompute both variants and seven figures
D = diagnose_tightening_margins(R);      % reuse the already solved exact QPs
V = diagnose_topology_switch(R);         % edge events, connectivity, safety, local checks
```

`reproduce_paper` overwrites its selected output folder (`results/` by default). To preserve the archived results and compare a new run against them:

```matlab
summary = verify_reproduction;          % outputs to reproduced_results/
figures = verify_figure_data;           % compare all seven native figure data sets
```

The diagnostics also accept a saved `experiment_data.mat` path or no arguments (default: archived results). Their default output directories are under `reproduced_results/`. They do not integrate or change figures. The online controller is the projected primal–dual ODE; `quadprog` initializes the fast states and solves offline reference QPs only.

## Parameter changes from the earlier case

| Setting | Earlier case | Current case |
| --- | ---: | ---: |
| Initial vertical positions of agents 2 and 5 | -1.15, +1.15 | -1.605, +1.605 |
| Filter time scale `epsilon` | 2e-5 | 2e-7 |
| Dissipativity activation parameter `nu` | 0.02 | 0.1 |
| Connectivity physical-barrier normalization `connScale` | 0.1 | 1 |
| `ode15s` settings | default NDF; 2e-8 / 1e-11 / 0.25 | BDF, maximum order 2; 2e-9 / 1e-12 / 0.125 |
| Diagnostic geometry | diameter 2.5; vertical box [-1.3,1.3] | diameter 3.22; vertical box [-1.7,1.7] |
| Diagnostic input-ball radius | 15 | 35 |

The last three solver numbers are `RelTol / AbsTol / MaxStep`. The plant, game, targets, obstacles, `d_c=3.2`, safety distance 0.345, gains `c=40`, `gamma=45`, `K1=14.5`, `K2=8.8`, `xi=1.4`, and physical tightenings `(1.2,1.8,2.5)` remain unchanged. Both variants use the same settings except that the ablation sets all three tightenings to zero; their fast states start at their respective QP solutions.

The `connScale` change is a controller/QP parametrization change, not merely a reporting-unit change: it scales the physical input coefficients and affine/tightening terms while retaining the manuscript mismatch coefficients and their regularization. The QP and dynamic filter therefore must be recomputed, as done here. The physical collision/obstacle terms retain scale 0.06; the matrix-tree expression retains division by 512.

The declared geometry and input ball are **offline diagnostic sets**, not actuator limits or proved invariant regions. The old `r_u=0.165` and `delta_c=0.055` budgets are retained as explicitly historical metadata only. Current acceptance checks the computed inequalities directly; no replacement constants were selected to force a pass.

## Sampling and true edge events

There are **30,001 saved states per variant**, on the union of `0:0.1:1200` and:

```matlab
[0:1e-5:.02, .0201:1e-4:.1, .101:.001:2, 2.01:.01:150]
```

These are observations of the same integration. Dense initial sampling is necessary: a uniform 0.1 s grid can miss fast filter-tracking peaks. Nonterminal ODE events locate zeros of `||y_i-y_j||^2-d_c^2`; they do not reset states or change adjacency. Exactly zero inactive QP rows are excluded from dual-face event detection, while positive/subnormal weights are not thresholded.

The graph uses `chi = exp(-r^2/(d_c^2-r^2))` for `r<d_c` and zero otherwise. Mathematical support and floating-point nonzero support are reported separately because the exponential can underflow near the boundary.

| Quantity | Tightened | Untightened |
| --- | ---: | ---: |
| Edge count | 9 → 10 | 9 → 10 |
| Geometric activation time of (2,5), s | 0.0497658976656 | 0.218548545026 |
| Floating nonzero support activation bracket, s | (0.0552, 0.0553] | (0.261, 0.262] |
| Minimum physical `lambda2` | 1.93981681563 | 1.93981681563 |
| Minimum collision barrier | 0.0827527941966 | -4.24067958554e-07 |
| Minimum obstacle barrier | 0.124137019364 | -9.72712628944e-07 |
| Minimum matrix-tree barrier | 0.0287249940356 | 8.32903163772e-06 |
| Final full error norm at 1200 s | 1.35048015232e-09 | 1.16456245167e-09 |

All minima in this table are over saved states. Both graphs remain connected at all saved states. The ablation's collision/obstacle excursions are small and persist under tighter integration tolerances. Its connectivity certificate remains positive; it must not be described as disconnected or as violating all three safety quantities.

## Theorem 2 diagnostics

The tracking error is `norm(u-uStar)` for the **same regularized decomposable QP**, including the mismatch variables, `xi=1.4`, and dissipativity rows. Saved primal/dual reference solutions are reused; old data without these references are processed with the original row builder and QP solver. The nominal input is not substituted for `uStar`.

Coefficient bounds use `2*regionDiameter`, twice the maximum obstacle-to-box-corner distance, and the analytical global matrix-tree gradient bound `C.Gh`. They are `(6.44, 11.148524566058057, 1.8512977520357721)` in physical barrier units. At each saved state, the affine discrepancy `q_c-sum(qhat_c)=da*u+db` has the exact input-ball supremum `abs(db)+35*norm(da)`. The diagnostic also checks the actual pairwise `1/2` and connectivity `1/N` tightening factors and all scaling identities.

| Quantity, physical units | Value |
| --- | ---: |
| Maximum sampled `r_u` | 0.0207623188575 |
| `delta_c`: sampled states × input ball of radius 35 | 0.90388111748 |
| Actual-input discrepancy (reported separately) | 0.232718992275 |
| `m_int = DeltaInt-sqrt(2)*barAInt*r_u` | 1.01090644723 |
| `m_obs = DeltaObs-barAObs*r_u` | 1.56853077817 |
| `m_c = DeltaC-sqrt(N)*barAC*r_u-delta_c` | 1.51017061392 |

All three margins are positive. Their scaled values are `(0.0606543868338, 0.0941118466901, 1.51017061392)`. All **60,002 QPs** (both variants) pass numerical feasibility, stationarity, complementarity and dual-feasibility checks in the original assembled units at `1e-8`. A positive-row-scaling fallback improves conditioning near tiny communication weights; its solution is accepted only after the same original-unit KKT checks. It does not alter the QP objective or loosen validation tolerances.

For the tightened run, the maximum primal, stationarity, complementarity and dual violations are `(4.43425378376e-11, 1.00933263361e-10, 9.16825569543e-09, 3.75128949643e-09)`. The maximum original physical ECBF residuals (`q<=0`) are `(-1.19493346187, -1.79983697271, -2.49801459136)` for collision, obstacle and connectivity. Active collision rows are checked by the independent diagnostic; the main audit also evaluates physical pair residuals directly, including zero-weight pairs. The actual dynamic local-row residual is reported separately and can be positive during finite-speed tracking without violating the physical ECBFs.

These are **trajectory-based numerical diagnostics**, not a Theorem 2 certificate over the entire compact region or all continuous times. Uniform QP feasibility, a continuous-time tracking/containment bound, the connectivity discrepancy over all `Theta_R × U_R^u`, and the a priori sensing-boundary speed condition remain unverified. `regionalUniformCertificate` and `continuousTimeRegionalCertificate` remain `false`. Analytical local conditions are evaluated separately, and Corollary 1 is not invoked.

## Figures and evidence

| Manuscript panel | File stem in `results/` |
| --- | --- |
| Fig. 1(a), tightened planar trajectories | `my_trajectory` |
| Fig. 1(b), untightened trajectories | `my_trajectory_2` |
| Fig. 2(a), collision barrier | `inter_agent_avoidance` |
| Fig. 2(b), obstacle barrier | `static_obstacle_avoidance` |
| Fig. 2(c), matrix-tree barrier | `connectivity_preservation` |
| Fig. 3(a), position components | `HB_flow_hebing` |
| Fig. 3(b), same-run local error on [850,1050] | `exp_2` |

All seven panels are exported as editable `.fig`, vector `.pdf`/`.eps`, and `.png`. The existing panel structure is retained. Only the tiny penetration label's unit and the connectivity inset range were adapted to display the new data accurately. The full source manuscript is not required or distributed; [verification/section_v.tex](verification/section_v.tex) provides the revised numerical section.

The release was recomputed from scratch: **14,000,650 numerical/logical values match the archived current results exactly**, and all seven figures' numerical graphics data match exactly. This is reproduction of the new case, not equality with the historical trajectory. Native FIG edit-save-reopen checks pass. The precision comparison additionally checks the same case with tolerances 10 times looser and maximum step twice as large; edge-event times change by less than 4e-8 s and the tracking maximum by about 0.016%. It is a numerical convergence check, not a validated error enclosure.

| Record | Location |
| --- | --- |
| Full trajectories, QP references, model and audits | `results/experiment_data.mat` |
| All parameters and analytical local bounds | `verification/numerical_certificate.txt` |
| Fresh-run array and figure comparisons | `verification/reference_comparison.json`, `figure_reference_comparison.json` |
| Command Window output and MATLAB environment | `verification/reproduction_log.txt`, `verification_summary.json` |
| Topology events, sample minima, all numerical gates | `verification/topology/topology_summary.json` |
| Edge counts, physical lambda2 and matrix-tree curves | `verification/topology/*_topology.csv` |
| Tightening table and detailed diagnostic arrays | `verification/topology/tightening/` |
| Precision comparison | `verification/topology/precision_comparison.json` |
| Source/evidence hashes and manuscript correspondence | `source_manifest.json`, `verification_manifest.json` |

## Minimal Section V text

> Edge (2,5) activates naturally at 0.049766 s, increasing the number of communication edges from nine to ten. Over 30,001 saved states, the maximum regularized-QP tracking error is 0.0207623, and the connectivity discrepancy maximized over the input ball of radius 35 is 0.903881. The resulting tightening margins are 1.01091, 1.56853, and 1.51017. All sampled QPs pass feasibility/KKT checks at 1e-8, and the tightened physical ECBF residuals are nonpositive. These trajectory-based numerical diagnostics do not establish a uniform regional certificate for Theorem 2.

## Citation

```bibtex
@unpublished{qian_safety_critical_gne,
  author = {Qian, Sichen and Liu, Hongzhe and Yu, Wenwu},
  title = {Distributed Safety-Critical Seeking of an Interior Generalized {Nash} Equilibrium over State-Dependent Graphs},
  note = {Manuscript prepared for submission to {IEEE} Transactions on Control of Network Systems}
}
```
