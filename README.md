# Sex Differences in Human Cortical Anatomy

Analysis code for:

> Lee, H. M. et al. **Regional sex differences in human cortical anatomy vary in their morphometric bases and overlap with sex chromosomal and gonadal influences.** *Nature Communications* 17, 7503 (2026). https://doi.org/10.1038/s41467-026-74274-8
## Overview

We mapped sex differences in regional cortical anatomy in three cohorts (HCP-YA, UK Biobank; combined N = 1,754). We then compared these maps with the anatomical effects of X- and Y-chromosome dosage and testicular hormone production from three clinical cohorts (combined N = 313).

## Dependencies

- MATLAB R2016b or later, with the Statistics and Machine Learning Toolbox and the Parallel Computing Toolbox
- [SurfStat](https://www.math.mcgill.ca/keith/surfstat/) (`SurfStatLinMod`, `SurfStatT`, `term`)
- [BrainSpace](https://github.com/MICA-MNI/BrainSpace) (`surface_to_graph`, `compute_mem`, `moran_randomization`)
- `fdr_bh` (Benjamini–Hochberg FDR) and `hexscatter` (MATLAB File Exchange)
- `daviolinplot` ([DataViz](https://github.com/povilaskarvelis/DataViz))

Some custom visualization functions (`loadColorMaps`, `ThreshColorMap`, `HyoSurfacePlotFunction`) are not included. They only affect figure rendering, not the statistics.

## Data

Participant-level data are not included in this repository.

- **HCP** (S1200 release): https://db.humanconnectome.org/, available to registered users who accept the data use terms
- **UK Biobank**: https://www.ukbiobank.ac.uk, available to approved researchers via application
- **NIH clinical cohorts** (sex chromosome aneuploidies, isolated gonadotropin-releasing hormone deficiency): not publicly available because of participant privacy; available from the corresponding author on reasonable request with a data use agreement

Group-level statistical maps are provided as source data on publication.

## License

Released under the MIT License (see `LICENSE`).

## Acknowledgments

Supported by the Intramural Research Program of the National Institute of Mental Health (1ZIAMH002949).

## Contact

Hyo Min Lee — hyo.lee@nih.gov
