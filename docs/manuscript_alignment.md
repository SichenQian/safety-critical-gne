# Manuscript and implementation alignment

The archived [numerical section](numerical_section.tex) matches the current manuscript. Its TEX and PDF SHA-256 values are recorded in [figure_sources.json](figure_sources.json).

Table I includes the applied-input peak 13.480099 below the input-ball radius 15. The scoped Theorem 2 records distinguish this trajectory observation from the all-time input-membership premise and the four retained regional/tracking hypotheses; see [theorem2_conditions.md](theorem2_conditions.md).

The independent local experiment initializes $u(0)=u_{\mathrm{nom}}(\mathbf y(0))$ and $\|\zeta_0(0)\|=r_f/2$. The perturbation is applied to the first retained mismatch coordinate; all other mismatch coordinates and all multipliers start at zero, including the multipliers of the locally omitted pair $(3,4)$.

This initialization agrees with [run_local.m](../matlab/run_local.m), [the saved initial-state record](../reference/local/initial.json), and Table I of the manuscript. The initial fast error is approximately $3.29491\times10^{-7}$.

Table I reports the connectivity minimum of the nominal HB flow. The independent computation and precision repeat are recorded in [nominal_hb.json](../verification/nominal_hb.json) and reproduced with `reproduce_paper('nominal')`.
