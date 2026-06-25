# =====================================================================
# 01_exposure.R
# High-temperature exposure construction
#
# Builds the district-year exposure table used in the analysis:
#   - district-specific warm-season (May-September) temperature
#     thresholds at the 95th / 97.5th / 99th percentiles, fixed over
#     the study period (2007-2023);
#   - the annual number of high-temperature days per district;
#   - 1-year lags and the 2- and 3-year cumulative exposure windows;
#   - annual mean relative humidity and the area deprivation index (ADI).
#
# Input  : daily mean temperature and relative humidity by district
#          (Korea Meteorological Administration; not shared), each
#          district mapped to its nearest observation station by
#          district-centroid distance.
# Output : analysis-ready district-year exposure table.
#
# Note   : raw meteorological and NHIS files are not redistributable;
#          file paths below are placeholders.
# =====================================================================

library(dplyr)
library(lubridate)
library(zoo)
library(tidyr)
library(data.table)

# ---- 0. Inputs ------------------------------------------------------
# `daily` : district-day records with columns
#           district (NHIS district code), date, tem (mean temp, C),
#           hum (mean relative humidity, %).
# `adi`   : district-level area deprivation index (columns: district, sum_z).
# District codes are harmonised to a single vintage beforehand so that
# districts merged or renamed during 2007-2023 share one identifier.
daily <- readRDS("data/daily_temperature_humidity_by_district.rds")
adi   <- read.csv("data/area_deprivation_index.csv")

WARM_MONTHS <- 5:9   # warm season: May-September

# ---- 1. Date / season fields ---------------------------------------
daily <- daily %>%
  mutate(date  = as.Date(date),
         year  = year(date),
         month = month(date),
         warm  = as.integer(month %in% WARM_MONTHS))

# ---- 2. District-specific warm-season thresholds (fixed) -----------
# Percentiles are computed once over the full study period so that the
# threshold defining a "high-temperature day" is fixed within a
# district and does not drift across years.
thr <- daily %>%
  filter(warm == 1) %>%
  group_by(district) %>%
  summarise(thr_95  = quantile(tem, 0.95,  na.rm = TRUE),
            thr_975 = quantile(tem, 0.975, na.rm = TRUE),
            thr_99  = quantile(tem, 0.99,  na.rm = TRUE),
            .groups = "drop")

# ---- 3. Flag high-temperature days (warm season only) --------------
daily <- daily %>%
  left_join(thr, by = "district") %>%
  mutate(hd_95  = as.integer(warm == 1 & tem >= thr_95),
         hd_975 = as.integer(warm == 1 & tem >= thr_975),
         hd_99  = as.integer(warm == 1 & tem >= thr_99))

# ---- 4. Annual counts and annual mean humidity ---------------------
annual <- daily %>%
  group_by(district, year) %>%
  summarise(hd95  = sum(hd_95,  na.rm = TRUE),
            hd975 = sum(hd_975, na.rm = TRUE),
            hd99  = sum(hd_99,  na.rm = TRUE),
            mean_hum = mean(hum, na.rm = TRUE),
            .groups = "drop")

# ---- 5. Lags and cumulative exposure windows -----------------------
# Primary exposure: HW975_2YR (sum of high-temperature days at the
# 97.5th percentile over the two preceding years, t-1 + t-2).
exp_tab <- annual %>%
  arrange(district, year) %>%
  group_by(district) %>%
  mutate(across(c(hd95, hd975, hd99),
                list(lag1 = ~lag(.x, 1),
                     lag2 = ~lag(.x, 2),
                     lag3 = ~lag(.x, 3)),
                .names = "{.col}_{.fn}"),
         # annual mean humidity matched to each exposure window
         hum_2yr = rollmean(mean_hum, 1, fill = NA, align = "right"),
         hum_3yr = rollmean(mean_hum, 2, fill = NA, align = "right")) %>%
  ungroup() %>%
  mutate(HW95_2YR  = hd95_lag1  + hd95_lag2,
         HW975_2YR = hd975_lag1 + hd975_lag2,
         HW99_2YR  = hd99_lag1  + hd99_lag2,
         HW95_3YR  = hd95_lag1  + hd95_lag2  + hd95_lag3,
         HW975_3YR = hd975_lag1 + hd975_lag2 + hd975_lag3,
         HW99_3YR  = hd99_lag1  + hd99_lag2  + hd99_lag3)

# ---- 6. Merge area deprivation index -------------------------------
exp_tab <- exp_tab %>%
  left_join(adi, by = "district")

# ---- 7. Keep the analysis window and export ------------------------
exp_tab <- exp_tab %>%
  filter(year >= 2007 & year <= 2023) %>%
  select(district, year,
         hd95, hd975, hd99,
         hd975_lag1, hd975_lag2, hd975_lag3,
         HW95_2YR, HW975_2YR, HW99_2YR,
         HW95_3YR, HW975_3YR, HW99_3YR,
         mean_hum, hum_2yr, hum_3yr, sum_z)

fwrite(exp_tab, "output/district_year_exposure.txt", sep = "|", na = "")

# END
