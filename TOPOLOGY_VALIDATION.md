# Topology and tightening validation

Validated on 2026-10-24 with MATLAB R2024a / Optimization Toolbox 24.1. Scope: **trajectory-based numerical diagnostics** of the existing Section V experiment.

## Selection and preservation

The original case at `de46af77d4ecf07194040a6be1ba254b9ceab329` has ten edges throughout its saved trajectory. An independent working clone was used for all searches. The original local checkout and its archived figures/results remain unchanged. The current release uses the independently verified selected case; the old case remains retrievable by commit and the original v1 tag remains unchanged.

Before controller retuning was authorized, ten completed candidate configurations varied the radius or initial geometry (plus recomputed local proof constants where needed). None met both switching and the computed tightening inequalities: six retained ten edges, while four activated edge (2,5) but had large tracking errors and negative margins. Their complete numerical summary is [original_gain_candidates.csv](verification/topology/original_gain_candidates.csv). These failed cases were not promoted to paper figures. Prefix runs were parameter-screening runs of this same case. A slower-epsilon full attempt stalled numerically and was not classified as a completed feasible or infeasible trajectory.

After the permitted expansion of settings, the selected initial vertical coordinates of agents 2 and 5 are +/-1.605; epsilon is 2e-7, nu is 0.1, and connectivity physical-barrier normalization is 1. The communication radius remains 3.2. The same settings are used in both variants. These changes preserve the algorithm family and plant; they do change the regularized filter parametrization and trajectories. The full-row positive rescaling inside the offline QP solver is a separate, algebraically equivalent conditioning step.

The graph law, physical ECBF/mismatch row formulas, determinant derivatives, nominal seeker and ODE right-hand side are unchanged. A nonterminal observer detects the distance cutoff, and inactive identically zero constraints no longer cause false dual-face events. On the original case, adding this observer/zero-row handling produced bitwise identical trajectories; see [baseline_observer_validation.json](verification/topology/baseline_observer_validation.json).

## Numerical gates

| Requirement | Result | Interpretation |
| --- | --- | --- |
| Actual graph change | (2,5), 9 to 10 edges; tightened 0.0497658976656 s; untightened 0.218548545026 s | Distance-cutoff event, not scheduled adjacency |
| Physical graph connectivity | Minimum lambda2 1.93981681563 in both runs | Positive at every saved state |
| Tightened matrix-tree barrier | Minimum 0.0287249940356 | Positive sufficient certificate at samples |
| Tightened collision / obstacle | Minima 0.0827527941966, 0.124137019364 | Positive at samples |
| Exact regularized QPs | 30,001 per variant; all KKT checks <= 1e-8 | No nominal-input substitution or tolerance relaxation |
| Three tightening margins | 1.01090644723, 1.56853077817, 1.51017061392 | Independently maximized sampled tracking and input-ball discrepancy |
| Original tightened physical ECBFs | Maxima -1.19493346187, -1.79983697271, -2.49801459136 | Convention q <= 0 |
| Final aggregate error | 1.35048015232e-09, 1.16456245167e-09 | Both < 1e-6 at 1200 s |
| Same-run local segment | First admissible saved time 834.7; displayed [850,1050] | No state reset or separate local experiment |
| From-scratch reproduction | 14,000,650 values, maximum difference 0 | Current archive reproduced exactly |
| Figure correspondence | 7/7 native data sets exact | Existing 2/3/2 panel structure |

The untightened collision and obstacle minima are -4.24067958554e-07 and -9.72712628944e-07. The obstacle penetration in the highlighted inset is about 0.778 micrometres. Its matrix-tree minimum is **positive**, 8.32903163772e-06. The revised narrative reports these small excursions without claiming disconnection or large practical violations. This case supports natural support switching and the tightening mechanism, but the ablation's violation magnitude is weaker evidence than in the original case.

## Affected theoretical conditions

| Condition | Recomputed value | Status / boundary |
| --- | --- | --- |
| Strong monotonicity / Lipschitz constant | 1, 1.15 | Unchanged game |
| Nominal gains | c=40 > 37.0170288156; gamma=45 > 25.5680375357 | Analytical inequalities pass for lambdaProof=0.2 |
| Local fast/slow epsilon condition | epsilon=2e-7 < 4.24518763044e-05; min eig(Q)=0.636189779148 | Analytical local inequality passes |
| Target connectivity | lambda2=3.12057238925; h_c=0.317861115664 | Strict |
| Target reduced-QP slack | 0.001 | Strict, with the manuscript's locally omitted rows treated as specified |
| Local graph lower bound | 3.1153697883 > 3 on radius 0.0002 | Analytical local bound |
| Local ECBF/extended set bounds | 0.182173763523, 1.54047892, 0.317490856114 / 0.398653834385, 3.38086258584, 0.69666797142 | Positive |
| Retained local fast-tube slack | 0.0062989887831, 1.20495165753, 0.187892990153, 0.001 | All positive with current scaling and nu |
| Activation neighbourhood | rTheta=0.0002 <= nu/2=0.05 | Satisfied |
| Sensing-boundary inequality | 22.2201879955 >= 6.4 | Conditional on a priori relative speed <= 1 |
| Actual sampled speed | max agent speed 0.433936857661 < 0.5 | Does not prove the a priori boundary-speed assumption |
| Declared position region | diameter 3.22, x in [-3.1,3.1], y in [-1.7,1.7] | Contains saved physical states; not a proved invariant extended-state set |
| Diagnostic input ball | radius 35; maximum actual norm 32.7785857051 | Contains sampled actual and optimizer inputs; not a continuous-time proof |
| Coefficient bounds | 6.44, 11.1485245661, 1.85129775204 | Analytical bounds on the declared geometry / global matrix-tree derivative |
| QP feasibility on every point in Theta_R | Not established | Finite-sample feasibility only |
| Uniform tracking / discrepancy / containment | Not established | ru and delta are empirical state-sampled maxima |
| Corollary 1 entrance condition | Not invoked | Local-entry check is along this run |

The epsilon/local certificate is recomputed after both connScale and nu change; it is not copied from the original case. The original budgets 0.165 and 0.055 do not certify the new case, and are not used in acceptance. The reported ru is the observed maximum itself. Delta is evaluated over all inputs in the declared ball at each saved state, not just the applied input. Its value is 0.90388111748, compared with the actual-input-only discrepancy 0.232718992275.

## Precision and reproducibility

The released settings use ode15s BDF of maximum order 2, RelTol=2e-9, AbsTol=1e-12 and MaxStep=0.125. A full comparison used RelTol=2e-8, AbsTol=1e-11 and MaxStep=0.25 with all controller parameters and observation times unchanged. The tracking maximum changed by 0.01575369%; the tightened and untightened event times changed by 1.66214216904e-08 and 3.74030086903e-08 s. The signs of the three ablation minima were unchanged. The maximum physical-position differences over the entire trajectories were 2.42284179244e-06 and 3.69504063238e-06. Thus this comparison is not a certified enclosure for the very small penetrations.

To repeat the coarser-precision check on the **same** experiment without changing defaults or figures:

```matlab
loaded = load(fullfile('results','experiment_data.mat'),'result');
reference = loaded.result;
P = reference.P;
P.RelTol = 2e-8; P.AbsTol = 1e-11; P.maxStep = 0.25;
F = reproduce_paper([],'diagnostic_helpers');
coarse = struct('P',P,'C',F.certificates(P));
coarse.proposed = F.simulate(P,true);
coarse.auditProposed = F.audit(P,coarse.proposed,true,coarse.C);
coarse.untightened = F.simulate(P,false);
coarse.auditUntightened = F.audit(P,coarse.untightened,false,coarse.C);
V = diagnose_topology_switch(coarse,fullfile('reproduced_results','coarse_check'));
```

The geometric activation root and IEEE double-precision nonzero support differ near the exponential cutoff. For the tightened run, the floating support first changes in (0.0552,0.0553] s, compared with the geometric root at 0.0497658976656 s. There are 55 saved edge samples inside the geometric radius with an underflowed zero weight (43 in the ablation). Both definitions and their sample counts are archived; no ad hoc threshold is used to manufacture switching.

All saved data are ordinary floating-point computations. The result does not certify uniform regional feasibility, continuous-time safety between samples, or a global convergence basin. The numerical gates provide reproducible support for the revised example while keeping these theorem hypotheses separate.
