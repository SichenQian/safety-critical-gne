# Release and citation policy

## Current snapshot: v1.0.2

The paper cites `https://github.com/SichenQian/safety-critical-gne/tree/v1.0.2`. This version is published with GitHub release immutability enabled: after publication, GitHub locks the associated tag and attached assets. The full commit SHA is written in the release notes and the attached `release_snapshot.json`; `SHA256SUMS.txt` records the asset hashes. GitHub also generates its release attestation.

The publication workflow is: prepare the repository and manuscript locally, create a new tag, upload all assets to a draft, check their hashes, then publish the immutable release. Subsequent changes use a new version. Published tags are not moved and attached files are not replaced.

Official behavior: [GitHub immutable releases](https://docs.github.com/en/code-security/concepts/supply-chain-security/immutable-releases).

## Earlier v1.0.1 review snapshots

The v1.0.1 label moved during earlier document synchronization. These complete commit identifiers retain the corresponding reviewed trees:

| Snapshot | Full commit |
| --- | --- |
| Initial scoped Theorem 2 audit | `6efd97554ec9c5893babef69f8d8224b0e39b67a` |
| V3 manuscript alignment | `b71ae33ef804d29b670959f2021bccae9531c9b7` |
| Current-source hash synchronization | `641c4378a2e07b33ae4913775ae2d56b3cf42b06` |

The existing v1.0.1 tag and release are preserved at the last row. They are historical snapshots; v1.0.2 is the citation target for the current manuscript.

## Numerical data

The model, parameters, original FIG files, compact plot data and six full reference datasets are unchanged from commit `641c4378a2e07b33ae4913775ae2d56b3cf42b06`. This R2 revision updates version identification and documents only; it does not rerun the ODEs.

The large reference archives remain attached to v1.0.0. Their exact bytes are bound by the unchanged versioned URLs and SHA-256 values in `data/reference_archives.json`. The new immutable code package includes that descriptor and the compact data needed to redraw the figures.
