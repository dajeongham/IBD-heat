# =====================================================================
# 04_paf.R
# Population attributable fraction (PAF) for biologic-agent initiation
#
# Quantifies the share of biologic-agent initiations attributable to
# high-temperature days, using the model coefficient from PART A of
# 03_models.sas (primary exposure: HW975_2YR).
#
#   AF_y = (RR_y - 1) / RR_y,   RR_y = exp(beta * mean_exposure_y)
#
# The overall PAF is the event-weighted sum of the year-specific
# attributable events divided by the total number of events, NOT the
# unweighted average of the year-specific AFs.
#
# Inputs : beta (log-HR per high-temperature day) and its SE from the
#          primary model; a year-level table with the mean exposure and
#          the observed number of events.
# Output : overall PAF and the year-specific AF series for the figure.
# =====================================================================

library(dplyr)

# ---- 0. Model coefficient (from 03_models.sas, primary model) ------
beta <- 0.00949      # log hazard ratio per high-temperature day
# exp(beta) ~ 1.010 (95% CI 1.004-1.015)

# ---- 1. Year-level inputs ------------------------------------------
# year       : calendar year
# mean_expo  : mean HW975_2YR among those at risk that year
# events     : observed biologic-agent initiations that year
yearly <- read.csv("data/yearly_exposure_events.csv")

# ---- 2. Year-specific attributable fraction and events -------------
yearly <- yearly %>%
  mutate(RR        = exp(beta * mean_expo),
         AF        = (RR - 1) / RR,            # year-specific AF
         attr_evt  = AF * events)              # attributable events

# ---- 3. Overall PAF = event-weighted sum ---------------------------
overall_paf <- sum(yearly$attr_evt) / sum(yearly$events)

cat(sprintf("Overall PAF: %.1f%% (%.0f of %.0f events)\n",
            100 * overall_paf, sum(yearly$attr_evt), sum(yearly$events)))

# year-specific series (for the burden-over-time figure)
yearly %>%
  transmute(year, mean_expo, events,
            AF_pct = round(100 * AF, 1),
            attr_events = round(attr_evt, 1)) %>%
  print(row.names = FALSE)

# END
