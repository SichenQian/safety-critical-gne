# Fig. 2 and Section V replacement fragments

The files in this directory are replacement snippets, not a modified manuscript.

1. Copy `connectivity_check/algebraic_connectivity.eps` into the manuscript's
   existing `consensus_gain_c_20260925/` figure directory.
2. Replace the Fig. 2 environment with `figure2.tex`. It retains panels (a)--(c)
   and adds the physical-graph algebraic connectivity as panel (d).
3. Replace only the Section V paragraph interpreting the safety margins with
   `sectionV_replacement.tex`. Preserve the neighboring QP tracking, KKT, and
   convergence paragraphs.

The four-panel figure was checked with the manuscript's IEEE class at
single-column width. The final float position depends on manuscript assembly.
The original manuscript and original figure files are not modified by these
scripts. The graph and spectral results are observations at saved times.
