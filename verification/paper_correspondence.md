# Manuscript-to-code numerical correspondence

Audited on 2026-10-13 for `TCNS-reproducibility-v1`. The inspected manuscript is `generic-color.tex` / `generic-color.pdf`; the exact source hashes are in `paper_correspondence.json` and the root manifest.

The cost coefficients, initial/target positions, obstacle geometry, communication radius, safety distance, row scaling and controller gains match the inspected numerical section. Both 12001-by-140 state trajectories were recomputed. Their local-error panel comes from physical time 850--1050 in the same tightened run. The minimum tightened collision, obstacle and matrix-tree quantities are 0.082099993964505888, 0.12398039083157436 and 0.1716691868816303. The corresponding untightened minima are -0.00018601556003673692, -0.00036005806873201651 and -0.0022173313315357879. These are sampled diagnostics. Uniform regional feasibility/tracking is not asserted; analytical local bounds are reported separately.

The manuscript-side MAT records were compared against the archived repository records, and a fresh complete run was then compared against those same records. All compared numerical arrays are exactly identical in MATLAB R2024a; elapsed runtime and worker-count metadata are explicitly excluded. The fresh run ends with `VERIFICATION_PASS`.

| Manuscript panel | Active TEX figure source | Repository native figure | Numerical correspondence |
| --- | --- | --- | --- |
| Fig. 1(a) | `consensus_gain_c_20260925/my_trajectory_1.eps` | `results/my_trajectory.fig` | Exact scientific data |
| Fig. 1(b) | `consensus_gain_c_20260925/my_trajectory_2.eps` | `results/my_trajectory_2.fig` | Exact scientific data |
| Fig. 2(a) | `consensus_gain_c_20260925/inter_agent_avoidance.eps` | `results/inter_agent_avoidance.fig` | Exact scientific data |
| Fig. 2(b) | `consensus_gain_c_20260925/static_obstacle_avoidance.eps` | `results/static_obstacle_avoidance.fig` | Exact scientific data |
| Fig. 2(c) | `consensus_gain_c_20260925/connectivity_preservation.eps` | `results/connectivity_preservation.fig` | Exact scientific data |
| Fig. 3(a) | `consensus_gain_c_20260925/HB_flow_hebing.eps` | `results/HB_flow_hebing.fig` | Exact scientific data |
| Fig. 3(b) | `consensus_gain_c_20260925/exp_2.eps` | `results/exp_2.fig` | Exact scientific data |

Fresh figure data were independently compared by `verify_figure_data`: 7 generated figures passed, including all 7 manuscript panels. Native line/patch/errorbar data are checked; fonts, colors, axes layout and export metadata are presentation differences. The inspected manuscript EPS/PDF panels were visually checked. This release establishes data reproduction, not byte-identical EPS output.

The manuscript LaTeX source and PDF were left unchanged. This numerical check does not re-prove the manuscript's theorems; the explicit assumptions and claim boundaries continue to apply.
