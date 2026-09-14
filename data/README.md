# Main analysis panel

The panel contains 43,797 cell–years at 2,021 cells in China during 2001–2024.
`panel/manifest.json` lists the 115 plain-CSV parts, row counts, and
SHA-256 checksums. `code/reproduce_main.py`, run from the repository root,
reads these parts by default, verifies their checksums and concatenated content,
and joins them in manifest order. No manual merge is needed.

The parts reconstruct the complete Dataset S1 CSV without changing its rows or
values. The complete `Dataset_S1_main_panel.csv.gz` is supplied separately with
the submission; it is not duplicated in Git. An optional `--panel` argument
accepts that complete CSV or CSV.gz instead of the parts.

`Dataset_S1_variable_dictionary.csv` and `Dataset_S1_metadata.json` describe
columns, units, and provenance. Frozen coefficients are in
`../reference_results/main_coefficients.csv`.

Figure and supplementary-table source data accompany the submission separately.
Original GOSIF, GLEAM, MODIS, ERA5-Land, CERN, and FLUXNET source datasets are not
redistributed here. Obtain them from the providers identified in the manuscript
and SI, subject to their registration, licensing, and use terms. This repository
contains the derived main analysis panel, not the original product archives.

Parts are below 900 KiB each for browser upload. Their concatenated CSV is
byte-identical to the complete Dataset S1. Candidate-location information
not needed by the estimator is retained with the submission materials.
