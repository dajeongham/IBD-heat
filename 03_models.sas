/********************************************************************************
 03_models.sas
 Time-varying (counting-process) Cox proportional-hazards models

 All models use PROC PHREG with the (TSTART, TSTOP) counting-process syntax and
 Efron's method for ties, on the datasets from step 02. The hazard ratio is
 expressed per one-day increase in the cumulative number of high-temperature
 days. The primary adjustment set matches the manuscript:

   exposure + SEX_TYPE + START_AGE + SES (income) + Urbanicity (region type)
            + MA1_hum (annual mean humidity) + sum_z (ADI) + STD_YYYY (calendar year)

 Outline
   PART A  Primary models (97.5th, 2-year cumulative) -> Table 2
   PART B  Threshold x window grid (Models 1-3)        -> Tables S6-S7
   PART C  Subgroups (sex, age, income, region type)   -> subgroup figure
   PART D  Exposure quartiles (p-trend and type-3)      -> Table S8
   PART E  Sensitivity analyses                         -> Table S9
 Reference for the exposure term: HazardRatio is reported per units=1 day.
 ********************************************************************************/

libname B "<PATH-TO-WORKING-LIBRARY>";   /* SEV_BIO_COUNTING, SEV_SURG_COUNTING */

/* =====================================================================
   PART A. Primary models (primary exposure: HW975_2YR)
   ===================================================================== */
%macro primary(ds=, status=, tag=);
  proc phreg data=&ds;
    class SEX_TYPE(ref='2') Urbanicity(ref='Urban') SES(ref='MISS') / param=ref;
    model (TSTART, TSTOP)*&status(0) =
          HW975_2YR SEX_TYPE START_AGE SES Urbanicity MA1_hum sum_z STD_YYYY
          / ties=efron risklimits;
    hazardratio HW975_2YR / units=1;
    title "&tag - primary (HW975_2YR), per 1-day";
  run;
%mend;
%primary(ds=B.SEV_BIO_COUNTING,  status=BIO_STATUS,  tag=Biologic initiation);
%primary(ds=B.SEV_SURG_COUNTING, status=SURG_STATUS, tag=IBD surgery);

/* =====================================================================
   PART B. Threshold x window grid (Models 1-3)
   Model 1 : exposure + sex + age + humidity
   Model 2 : Model 1 + region type + income + ADI + calendar year  (primary)
   Model 3 : Model 2 + smoking + alcohol + physical activity
   Humidity term matches the window: 2-year -> MA1_hum, 3-year -> MA2_hum.
   ===================================================================== */
%macro grid(ds=, status=, expo=, hum=, tag=);
  proc phreg data=&ds;
    model (TSTART, TSTOP)*&status(0) = &expo SEX_TYPE START_AGE &hum
          / ties=efron risklimits;
    hazardratio &expo / units=1;  title "&tag - &expo - Model 1"; run;

  proc phreg data=&ds;
    class SEX_TYPE(ref='2') Urbanicity(ref='Urban') SES(ref='MISS') / param=ref;
    model (TSTART, TSTOP)*&status(0) =
          &expo SEX_TYPE START_AGE Urbanicity SES sum_z STD_YYYY &hum
          / ties=efron risklimits;
    hazardratio &expo / units=1;  title "&tag - &expo - Model 2"; run;

  proc phreg data=&ds;
    class SEX_TYPE(ref='2') Urbanicity(ref='Urban') SES(ref='MISS')
          SMOKE_TV ALCOHOL_TV_C PA_TV_C / param=ref;
    model (TSTART, TSTOP)*&status(0) =
          &expo SEX_TYPE START_AGE Urbanicity SES sum_z STD_YYYY &hum
          SMOKE_TV ALCOHOL_TV_C PA_TV_C
          / ties=efron risklimits;
    hazardratio &expo / units=1;  title "&tag - &expo - Model 3"; run;
%mend;

/* Biologic-agent initiation */
%grid(ds=B.SEV_BIO_COUNTING, status=BIO_STATUS, expo=HW975_2YR, hum=MA1_hum, tag=BIO);
%grid(ds=B.SEV_BIO_COUNTING, status=BIO_STATUS, expo=HW975_3YR, hum=MA2_hum, tag=BIO);
%grid(ds=B.SEV_BIO_COUNTING, status=BIO_STATUS, expo=HW95_2YR,  hum=MA1_hum, tag=BIO);
%grid(ds=B.SEV_BIO_COUNTING, status=BIO_STATUS, expo=HW95_3YR,  hum=MA2_hum, tag=BIO);
%grid(ds=B.SEV_BIO_COUNTING, status=BIO_STATUS, expo=HW99_2YR,  hum=MA1_hum, tag=BIO);
%grid(ds=B.SEV_BIO_COUNTING, status=BIO_STATUS, expo=HW99_3YR,  hum=MA2_hum, tag=BIO);
/* IBD-related surgery */
%grid(ds=B.SEV_SURG_COUNTING, status=SURG_STATUS, expo=HW975_2YR, hum=MA1_hum, tag=SURG);
%grid(ds=B.SEV_SURG_COUNTING, status=SURG_STATUS, expo=HW975_3YR, hum=MA2_hum, tag=SURG);
%grid(ds=B.SEV_SURG_COUNTING, status=SURG_STATUS, expo=HW95_2YR,  hum=MA1_hum, tag=SURG);
%grid(ds=B.SEV_SURG_COUNTING, status=SURG_STATUS, expo=HW95_3YR,  hum=MA2_hum, tag=SURG);
%grid(ds=B.SEV_SURG_COUNTING, status=SURG_STATUS, expo=HW99_2YR,  hum=MA1_hum, tag=SURG);
%grid(ds=B.SEV_SURG_COUNTING, status=SURG_STATUS, expo=HW99_3YR,  hum=MA2_hum, tag=SURG);

/* =====================================================================
   PART C. Subgroups and multiplicative interaction
   Stratified HRs by sex, age group, income, and region type, with the
   p-interaction from a Wald type-3 test of the exposure x subgroup
   product term (a joint test across the subgroup levels).
   ===================================================================== */
%macro subgroup(var=);
  /* stratified estimates */
  proc sort data=B.SEV_BIO_COUNTING out=_s; by &var; run;
  proc phreg data=_s;
    by &var;
    class SEX_TYPE(ref='2') Urbanicity(ref='Urban') SES(ref='MISS') / param=ref;
    model (TSTART, TSTOP)*BIO_STATUS(0) =
          HW975_2YR SEX_TYPE START_AGE SES Urbanicity MA1_hum sum_z STD_YYYY
          / ties=efron risklimits;
    hazardratio HW975_2YR / units=1;
    title "Subgroup by &var";
  run;
  /* interaction (type-3 Wald test of the product term) */
  proc phreg data=B.SEV_BIO_COUNTING;
    class SEX_TYPE(ref='2') Urbanicity(ref='Urban') SES(ref='MISS') &var / param=ref;
    model (TSTART, TSTOP)*BIO_STATUS(0) =
          HW975_2YR SEX_TYPE START_AGE SES Urbanicity MA1_hum sum_z STD_YYYY
          &var HW975_2YR*&var
          / ties=efron risklimits type3(wald);
    title "Interaction: HW975_2YR x &var";
  run;
%mend;
%subgroup(var=SEX_TYPE);
%subgroup(var=AGE_GRP);
%subgroup(var=SES);
%subgroup(var=Urbanicity);

/* =====================================================================
   PART D. Exposure quartiles (Table S8)
   HW975_QCAT : quartiles of HW975_2YR (Q1 = reference).
   - p-trend  : quartile rank entered as a single ordinal score (1 df).
   - overall  : 3-df Wald type-3 test of the categorical quartile term.
   ===================================================================== */
proc rank data=B.SEV_BIO_COUNTING out=_q groups=4;
  var HW975_2YR; ranks HW975_Q;
run;
data _q; set _q; HW975_QCAT = HW975_Q + 1; run;  /* 1..4 */

/* per-quartile day distribution (cutpoints / medians for the table) */
proc means data=_q n min p50 max maxdec=0;
  class HW975_QCAT; var HW975_2YR;
  title "Quartile cutpoints and medians (days)";
run;

/* categorical quartiles: HRs + 3-df type-3 test */
proc phreg data=_q;
  class SEX_TYPE(ref='2') Urbanicity(ref='Urban') SES(ref='MISS')
        HW975_QCAT(ref='1') / param=ref;
  model (TSTART, TSTOP)*BIO_STATUS(0) =
        HW975_QCAT SEX_TYPE START_AGE SES Urbanicity MA1_hum sum_z STD_YYYY
        / ties=efron risklimits type3(wald);
  title "Quartiles (categorical) - HRs and 3-df type-3 test";
run;

/* trend test: quartile rank as a continuous score (1 df) */
proc phreg data=_q;
  class SEX_TYPE(ref='2') Urbanicity(ref='Urban') SES(ref='MISS') / param=ref;
  model (TSTART, TSTOP)*BIO_STATUS(0) =
        HW975_QCAT SEX_TYPE START_AGE SES Urbanicity MA1_hum sum_z STD_YYYY
        / ties=efron risklimits;
  title "Quartile trend (ordinal score) - p for trend";
run;

/* =====================================================================
   PART E. Sensitivity analyses (Table S9)
   ===================================================================== */
/* E1. + lifestyle covariates: see Model 3 in PART B. */

/* E2. Calendar-period robustness: pre/post-2016 indicator (POST2016) */
proc phreg data=B.SEV_BIO_COUNTING;
  class SEX_TYPE(ref='2') Urbanicity(ref='Urban') SES(ref='MISS') POST2016(ref='0') / param=ref;
  model (TSTART, TSTOP)*BIO_STATUS(0) =
        HW975_2YR SEX_TYPE START_AGE SES Urbanicity MA1_hum sum_z POST2016
        / ties=efron risklimits;
  hazardratio HW975_2YR / units=1;
  title "Sensitivity: pre/post-2016 indicator";
run;

/* E3. Restriction to 2007-2015 */
proc phreg data=B.SEV_BIO_COUNTING(where=(std_yyyy <= 2015));
  class SEX_TYPE(ref='2') Urbanicity(ref='Urban') SES(ref='MISS') / param=ref;
  model (TSTART, TSTOP)*BIO_STATUS(0) =
        HW975_2YR SEX_TYPE START_AGE SES Urbanicity MA1_hum sum_z STD_YYYY
        / ties=efron risklimits;
  hazardratio HW975_2YR / units=1;
  title "Sensitivity: restricted to 2007-2015";
run;

/* E4. Single-year lag structure (t-1, t-2, t-3) */
%macro lag(expo=);
  proc phreg data=B.SEV_BIO_COUNTING;
    class SEX_TYPE(ref='2') Urbanicity(ref='Urban') SES(ref='MISS') / param=ref;
    model (TSTART, TSTOP)*BIO_STATUS(0) =
          &expo SEX_TYPE START_AGE SES Urbanicity MA1_hum sum_z STD_YYYY
          / ties=efron risklimits;
    hazardratio &expo / units=1;  title "Lag: &expo"; run;
%mend;
%lag(expo=hd975_lag1);   /* previous year   */
%lag(expo=hd975_lag2);   /* two years prior */
%lag(expo=hd975_lag3);   /* three years prior */

/* Note: the landmark-period sensitivity (0 and 2 years) is run by rebuilding
   the cohort in step 02 with the alternative landmark and re-running PART A. */

/* END */
