# Phenology and evaporation components

Core reproducibility files for **Total evaporation masks contrasting modeled water-loss responses to earlier spring activity**.

Yusupukadier Zimini, Qingjie Zhao, Liuyang Zhang, and Chao Yang.
Correspondence: yangchao161@szu.edu.cn.

## Reproduce the main estimates

The package contains the exact main panel: 43,797 cell–years at 2,021 sampled
cells in China during 2001–2024. 115 small plain-CSV parts reconstruct Dataset S1.

```bash
python -m pip install -r requirements-core.txt
python code/reproduce_main.py
```

The script verifies input checksums, fits the main common-panel model,
compares all 87 coefficients and confidence intervals with the frozen reference,
and checks component identities. It writes `results/reproduced_main.csv` and
`results/reproduced_main.verification.json`. Dependencies are pinned in
requirements-core.txt.

## Files and units

- `data/panel/`: main-panel parts and their manifest.
- `data/Dataset_S1_*`: variable dictionary and detailed metadata.
- `code/reproduce_main.py`: main estimator and verification routine.
- `reference_results/main_coefficients.csv`: frozen numerical comparison.

Positive A denotes earlier fluorescence onset. Coefficients are per one day;
the paper reports per ten days. Fraction coefficients are proportions, converted
to percentage points by multiplying by 100. E is GLEAM total evaporation, Et is
transpiration T, and En is the modeled non-transpiration remainder N = E − T.
Timing uses positive flux; water amounts retain signed rates. SI Methods provide
the model, window, sample, and uncertainty definitions.

This release covers the 87 main common-panel estimates. Figure and supplementary
table source data accompany the submission separately. The full local project
retains preprocessing, tower evaluation, sensitivity workflows, and plotting.
Those workflows and their caches are outside this core package. 

The manuscript is unpublished and has no assigned article DOI. The existing MIT
license applies to code; original data-provider terms continue to apply to
derived data. Sources are identified in the manuscript and SI.
