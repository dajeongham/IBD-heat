/********************************************************************************
 02_cohort.sas
 Analytic cohort and counting-process (person-interval) dataset

 Builds the datasets used by 03_models.sas:
   - incident IBD adults in the NHIS 10% simple random sample, 2007-2023;
   - a one-year landmark to mitigate immortal-time bias;
   - annual follow-up intervals carrying the time-varying cumulative
     high-temperature days and time-varying covariates;
   - two outcomes analysed independently: first biologic-agent initiation
     and first IBD-related surgery; all-cause death is a censoring event.

 Inputs : NHIS claims (not shared) and the district-year exposure table
          from step 01 (district_year_exposure).
 Output : SEV_BIO_COUNTING, SEV_SURG_COUNTING (TSTART, TSTOP, status, covariates)

 Codes (defined in the article Supplement, not reproduced here):
   - IBD            : ICD-10 K50 (Crohn's disease), K51 (ulcerative colitis)
   - biologic agent : NHIS prescription/administration codes
   - IBD surgery    : NHIS procedure codes
 ********************************************************************************/

libname raw  "<PATH-TO-NHIS-CLAIMS>";     /* not shared */
libname work2 "<PATH-TO-WORKING-LIBRARY>";

/* ---- 1. Incident IBD cohort ------------------------------------------------
   First IBD diagnosis date (index), restricted to incident adult cases with
   no IBD record in a disease-free baseline window. The landmark date is one
   year after the index date; subjects with an outcome or death within the
   landmark window are excluded. */
data cohort;
  set raw.ibd_incident;                 /* one row per person */
  landmark_date = intnx('year', index_date, 1, 'same');
  start_age     = floor((landmark_date - birth_date) / 365.25);
  if start_age >= 18;
  base_year = year(landmark_date);      /* first year at risk */
  end_date  = min(event_date, death_date, '31DEC2023'd);
  if end_date <= landmark_date then delete;   /* landmark exclusion */
run;

/* ---- 2. Expand to annual person-intervals ----------------------------------
   Each subject contributes one row per follow-up year on the time-on-study
   scale (origin = landmark date). STD_YYYY indexes the calendar year of the
   interval and is used to attach the time-varying exposure and covariates. */
data person_year;
  set cohort;
  fu_years = ceil((end_date - landmark_date) / 365.25);
  do k = 1 to fu_years;
    tstart   = k - 1;
    tstop    = k;
    std_yyyy = base_year + (k - 1);     /* calendar year of this interval */
    output;
  end;
run;

/* ---- 3. Attach time-varying exposure and area-level covariates -------------
   district_year_exposure (from step 01) supplies, per district and calendar
   year: the cumulative high-temperature-day variables (HW*_2YR, HW*_3YR),
   annual mean humidity, and the area deprivation index (sum_z). */
proc sql;
  create table py_exp as
  select a.*, b.HW95_2YR, b.HW975_2YR, b.HW99_2YR,
              b.HW95_3YR, b.HW975_3YR, b.HW99_3YR,
              b.hum_2yr as MA1_hum, b.hum_3yr as MA2_hum, b.sum_z
  from person_year as a
  left join work2.district_year_exposure as b
    on a.district = b.district and a.std_yyyy = b.year;
quit;

/* ---- 4. Define the event indicator within each interval --------------------
   The event is marked in the interval containing the outcome year; later
   intervals for a subject who already had the event are removed. Separate
   status variables are built for each outcome. Deaths are handled as
   censoring (no event flag). */
%macro make_status(out=, eventyr=, status=);
  data &out;
    set py_exp;
    &status = 0;
    if &eventyr ne . and std_yyyy = year(&eventyr) then &status = 1;
    /* drop person-time after the event */
    if &eventyr ne . and std_yyyy > year(&eventyr) then delete;
  run;
%mend;

%make_status(out=SEV_BIO_COUNTING,  eventyr=bio_date,  status=BIO_STATUS);
%make_status(out=SEV_SURG_COUNTING, eventyr=surg_date, status=SURG_STATUS);

/* SEV_BIO_COUNTING and SEV_SURG_COUNTING are the counting-process datasets
   consumed by 03_models.sas. */

/* END */
