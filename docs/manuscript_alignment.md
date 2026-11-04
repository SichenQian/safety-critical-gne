# Manuscript and implementation alignment

The archived [numerical section](numerical_section.tex) matches the current manuscript. Its TEX and PDF SHA-256 values are recorded in [figure_sources.json](figure_sources.json).

## Revision of 4 November 2026

| Manuscript location | Current statement | Reproduction source |
| --- | --- | --- |
| p. 6, before Theorem 2 | The regional discrepancy is bounded by a Lipschitz expression; a budget at least that expression is sufficient | [Budget definitions and scope](theorem2_conditions.md#regional-discrepancy-budget) |
| p. 7, numerical setup | Explicit local collision row with physical scale 0.06 and unscaled mismatch; reference QP uses the same adopted constraints | [Scaled-row coordinates](numerical_details.md#scaled-local-rows), audit_core.m and normalized_reference_rows.m |
| p. 8, Table I note | Nominal and tightened statistics use their respective 0–1200 s runs with numerical peak refinement; regional residence, QP feasibility, input boundedness, and uniform estimation/tracking-error bounds remain premises | [Scoped condition table](../verification/theorem2_conditions_with_scope.csv) |
| p. 8, local initialization | Appendix-derived slow/fast budgets and zero multiplier initialization | [Local construction](theorem3_local_certificate.md), run_local.m and reference/local/initial.json |
| pp. 8–9, Fig. 3 discussion/caption | Slow-error decay and fast transient reaching the numerical precision plateau | Compact data/paper_plot_data.mat and original figures/exp.fig |
| p. 8, reproduction link | tree/v1.0.2 is locked by an immutable release; the release snapshot records the full commit and asset hashes | [README](../README.md#run), reproduce_paper.m |
| p. 7, Theorem 3; p. 11, Appendix H | Retained-error estimate plus zeta_om(t)=exp(−xi t/epsilon) zeta_om(0) | [Retained and omitted coordinates](theorem3_local_certificate.md#4-实验初始化与结论范围) |

The manuscript also removes the visible first-page LOGO placeholder and remains 12 pages. The algorithm, parameters, archived figure files and six reference datasets are unchanged in this synchronization. The local radii are conservative sufficient neighborhoods, as specified in the construction; they do not estimate the largest region of attraction.

Table I includes the applied-input peak 13.480099 below the input-ball radius 15. The scoped Theorem 2 records distinguish this trajectory observation from the all-time input-membership premise and the four retained regional/tracking hypotheses; see [theorem2_conditions.md](theorem2_conditions.md).

The independent local experiment initializes $u(0)=u_{\mathrm{nom}}(\mathbf y(0))$ and $\|\zeta_0(0)\|=r_f/2$. The perturbation is applied to the first retained mismatch coordinate; all other mismatch coordinates and all multipliers start at zero, including the multipliers of the locally omitted pair $(3,4)$.

This initialization agrees with [run_local.m](../matlab/run_local.m), [the saved initial-state record](../reference/local/initial.json), and Table I of the manuscript. The initial fast error is approximately $3.29491\times10^{-7}$.

Table I reports the connectivity minimum of the nominal HB flow. The independent computation and precision repeat are recorded in [nominal_hb.json](../verification/nominal_hb.json) and reproduced with `reproduce_paper('nominal')`.

The R2 update changes only the manuscript reproduction URL. The supplied manuscript already includes the Appendix D correction; the implemented dissipativity row and all numerical results remain unchanged. See [release_policy.md](release_policy.md) for the retained earlier commits.
