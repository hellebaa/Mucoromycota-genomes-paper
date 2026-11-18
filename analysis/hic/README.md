# Hi-C reference annotations

The first shared release, **hic_reference_20261118_v1** (18 November 2026), contains
candidate-centromere and terminal-repeat annotations for three Mucoromycota genomes.
This is a versioned research snapshot for collaborative analysis and manuscript
preparation, not a final publication deposition or functional validation.

## Start here

1. Download the compact release (23 files, less than 100 kB).
2. Open the candidate BED and stringent telomere BED for the species of interest.
3. Consult the evidence tables before interpreting a region as a centromere or telomere.

```bash
# From the repository root, after configuring NeLS access (NELS_SETUP.md):
bash analysis/hic/download_data.sh
# Optional: download the three assemblies and original merged Hi-C matrices (~0.5 GB):
bash analysis/hic/download_data.sh --inputs
```

Downloads use the shared `scripts/download_files.sh` and verify every file's SHA256.
NeLS requires project access. Remote files mirror the repository-relative paths under
`Projects/NMBU_threeD_Sandve_2024/Mucor/`. Code and manifests are on GitHub; data and
results are not committed to Git. Public deposition with a persistent identifier
will be arranged for publication.

## Species and annotation scope

| Genome ID | Species | Candidates | Explicit unresolved/excluded contigs |
|---|---|---:|---|
| gzMucCirc1 | Mucor circinelloides | 10 | 5 and 10 unresolved |
| gzUmbVina1 | Umbelopsis vinacea | 12 | 13–16 have no retained candidates |
| ggMorZona1 | Mortierella zonata | 20 | 18 excluded |

The assemblies contain additional small contigs. An omitted contig has **no retained
annotation**, not demonstrated absence of a centromere. Complete contig lengths are
provided. Broad intervals are uncertainty envelopes, not molecular boundaries.

## Files on NeLS

Every release filename begins `hic_reference_20261118_v1__`. The shortened tree below
omits that prefix for readability; `<genome>` means one file for each of the three IDs.

```text
results/
├── centromeres/tables/                       # 12 files
│   ├── <genome>__candidates.bed               # 3: broad candidate envelopes
│   ├── <genome>__enrichment_targets.bed       # 3: windows used for contact analysis
│   ├── <genome>__chrom.sizes                  # 3: all contigs
│   ├── gzMucCirc1__shortlisted_motifs.bed     # 1: narrow sequence motif hits
│   ├── evidence.tsv                         # 1: 42 candidate records
│   └── unresolved.tsv                       # 1: explicit exclusions
├── assembly_structure/tables/                # 7 files
│   ├── <genome>__telomeres_all.bed            # 3: permissive terminal calls
│   ├── <genome>__telomeres_stringent.bed      # 3: stricter terminal calls
│   └── telomere_evidence.tsv                 # 1: all calls and filter decisions
└── hic/
    ├── tables/assemblies.tsv                 # 1: assembly/matrix paths and hashes
    └── logs/                                # 3
        ├── README.md                        # this guide
        ├── provenance.json                  # source identities and build record
        └── checksums.tsv                    # checksums of the other 22 files
```

The exact 23-file inventory, including the checksum table's own hash, is
[config/reference_v1_files.tsv](config/reference_v1_files.tsv).
The optional input inventory is [config/reference_v1_inputs.tsv](config/reference_v1_inputs.tsv).
Detailed enrichment tables, per-contig contact matrices, exploratory figures and
motif-discovery intermediates are **not** included in this first release.

## Coordinates and evidence

All BED intervals are **0-based, half-open**. BED score 0 is a placeholder, not a
confidence score; strand `.` means not applicable. Contig names match the assemblies
identified by SHA256 in the assembly table.

- **Candidate BEDs:** union of reviewed broad intervals and refined contact windows.
- **Enrichment-target BEDs:** windows used in the previous contact-enrichment analysis
  (nominally 96 kb, clipped at contig ends); preserve these for reproducing that analysis.
- **MucCirc1 motif BED:** ten shortlisted motif hits, separate from the published
  BLAST-based anchors and broad Hi-C regions. They are not interchangeable coordinates.
- **Evidence TSV:** source table, selection rationale, broad bounds, enrichment window,
  anchor, evidence type and qualifications. No numerical confidence score is inferred.

Candidates were reviewed using interchromosomal clustering and local Hi-C patterns,
then refined and evaluated using candidate-to-candidate contacts. MucCirc1 also has
prior published BLAST-based positions from Baalsrud et al. (G3,
[doi:10.1093/g3journal/jkaf207](https://doi.org/10.1093/g3journal/jkaf207)), using MucCirc3
centromere sequence evidence. UmbVina1 motif/LINE support varies among candidates;
MorZona has no established universal centromere motif. Same-data refinement and
contact enrichment are descriptive support, **not independent validation**.

MorZona terminal-repeat calls use **TTTTTTAGGG**, reverse complement CCCTAAAAAA.
MucCirc1 and UmbVina1 use the earlier multi-motif scan; motifs are recorded per call.
The stringent subset requires >=5 consecutive copies, <=100 bp from the contig end
and the expected strand orientation. These are descriptive filters, not a validated
classifier. Counts represent terminal-side calls, not chromosomes. There are 87
permissive calls and 85 stringent calls across the three assemblies.

## Matrix status and remaining work

The original merged matrices are already on NeLS in `data/HiC/mcool/`, with hashes
in the assembly table. Original per-library data are under
`data/HiC/Mucor_circinelloides/`, `Mortierella_zonata/` and `Umbelopsis_vinacea/`.
The optional download retrieves merged matrices, not all per-library files.

This snapshot **does not incorporate the pending short-cis-filtered/rebalanced
analysis**. Library review identified a large short-cis inward-facing component in
MorZona. Updated filtering, balancing and comparison of centromere versus telomere
clustering may lead to a new annotation version, especially for terminal candidates.
Cross-species synteny-based centromere predictions are outside this release.

## Reproduction and contributions

`python3 analysis/hic/scripts/build_reference_release.py` exports these annotations
from the existing curated source tables and genome FASTAs; it does not rediscover
centromeres. It checks coordinates and candidate counts, and refuses to change
existing frozen outputs. It uses the Python standard library. Source paths and
hashes are listed in [config/reference_v1_sources.tsv](config/reference_v1_sources.tsv).
Those upstream analytical tables are **not yet in this compact NeLS release**;
rebuilding from them requires the project's local source dataset. Downloading and
using the frozen BEDs requires none of those intermediates. The enrichment workflow
and its detailed inputs/results will be published in a separate step.

For revised annotations, use a new version prefix, document changed coordinates and
rationale, and retain this snapshot. Keep code in `analysis/<analysis>/`, inputs in
`data/`, and generated tables, figures and logs in the appropriate `results/<analysis>/`
folders. Commit code/documentation to GitHub, upload results to NeLS with checksums,
and record the Git commit and input identities. Do not overwrite the shared release.
