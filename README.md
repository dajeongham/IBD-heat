# Cumulative high-temperature days and clinical prognosis in inflammatory bowel disease

Analysis code for a nationwide, population-based cohort study examining whether the
cumulative number of high-temperature days is associated with two markers of clinical
prognosis in inflammatory bowel disease (IBD): first biologic-agent initiation and
first IBD-related surgery.

This repository contains **code only**. The underlying National Health Insurance
Service (NHIS) claims data cannot be shared (see *Data availability* below), so the
scripts are provided for **transparency** — to document exactly how the exposures were
constructed and how the models were specified — rather than for end-to-end re-execution.

## Data availability

The NHIS claims data underlying this study cannot be shared, owing to data-sharing
restrictions imposed by the National Health Insurance Service of Korea. Daily
temperature records are from the Korea Meteorological Administration (KMA). No raw
data, identifiers, or individual-level values are included in this repository.

## Study design (summary)

- **Source population:** a 10% simple random sample of the NHIS-insured population.
- **Cohort:** adults with incident IBD, 2007–2023; a one-year landmark was applied
  to mitigate immortal-time bias.
- **Exposure:** district-specific warm-season (May–September) temperature thresholds
  at the 95th / 97.5th / 99th percentiles, **fixed** over the study period. The number
  of high-temperature days was counted per district-year, using the KMA station nearest
  to each district centroid (246 districts, 17 years). The **primary exposure** is the
  two-year cumulative number of high-temperature days at the 97.5th percentile
  (sum over t−1 and t−2), updated annually.
- **Outcomes:** first biologic-agent initiation; first IBD-related surgery (analyzed
  independently). All-cause death was treated as a censoring event.
- **Models:** time-varying (counting-process) Cox proportional-hazards models, reporting
  the hazard ratio per one-day increase in cumulative high-temperature days, adjusted for
  sex, baseline age, income level, region type, annual mean relative humidity, area
  deprivation index, and calendar year.

## Repository structure

```
IBD-heat/
├── 01_exposure/      # KMA station matching; fixed percentile thresholds; cumulative windows
│   └── 01_exposure.R
├── 02_cohort/        # incident-IBD cohort, 1-year landmark, counting-process dataset
│   └── 02_cohort.sas
├── 03_models/        # time-varying Cox models, subgroups, sensitivity (all parts)
│   └── 03_models.sas
├── 04_paf/           # population attributable fraction (event-weighted)
│   └── 04_paf.R
├── LICENSE
└── README.md
```

## Workflow

Run the scripts in numerical order; outputs of each step feed the next.

1. **`01_exposure/`** — build the district-year high-temperature-day table
   (thresholds, station matching, 2- and 3-year cumulative windows).
2. **`02_cohort/`** — assemble the analytic cohort, apply the landmark, and create the
   person-interval (counting-process) dataset with annually updated covariates.
3. **`03_models/`** — fit the primary models and the threshold/window grid, plus the
   subgroup, quartile, and sensitivity analyses (organized as PARTS A–E within the file).
4. **`04_paf/`** — compute the overall (event-weighted) and year-specific attributable
   fractions.

## Software

- **SAS** 9.4 (`PROC PHREG` with the counting-process / time-varying specification)
- **R** ≥ 4.2 (exposure construction, PAF, figures, mapping)

Required R packages are listed at the top of each script.

## Citation

If you use this code, please cite the associated article (citation to be added upon
publication) and this repository.

## License

Released under the MIT License (see `LICENSE`).

## Contact

Dajeong Ham — first author. Questions about the code may be directed via the issues
tab of this repository.
