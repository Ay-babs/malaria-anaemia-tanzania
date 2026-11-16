#===============================================================================
# SAMPLE FLOW, MISSING DATA AND EXCLUSIONS
# Run this AFTER your main script has created these objects:
#   child_data, ml_an_df, mal_anemia, data_clean,
#   design_A, design_B, design_joint, model_A_final, model_B_final, model_joint
# (all counts are UNWEIGHTED numbers of children)
#===============================================================================
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

#-------------------------------------------------------------------------------
# 1. PARTICIPANT FLOW: how many children, and who was excluded at each step
#-------------------------------------------------------------------------------
n_all      <- nrow(child_data)                                   # all living children in the KR file
n_age      <- nrow(mal_anemia)                                   # aged 6-59 months
n_rdt      <- sum(!is.na(mal_anemia$malaria_child))              # valid RDT result
n_hb       <- sum(!is.na(mal_anemia$anaemia_child))              # valid haemoglobin
n_both     <- sum(!is.na(mal_anemia$malaria_child) &
                  !is.na(mal_anemia$anaemia_child))              # valid RDT AND Hb
n_analytic <- nrow(data_clean)                                   # complete on screening covariates

flow <- tibble(
  Step = c("Children in the KR recode (all living children)",
           "Aged 6-59 months (eligible for RDT / haemoglobin)",
           "With a valid RDT AND haemoglobin result",
           "With complete screening covariates (analytic sample)"),
  n = c(n_all, n_age, n_both, n_analytic),
  Excluded = c(NA, n_all - n_age, n_age - n_both, n_both - n_analytic),
  Reason = c("",
             "Aged under 6 months (not eligible for testing)",
             sprintf("No valid RDT (n = %d missing) or haemoglobin (n = %d missing) result",
                     n_age - n_rdt, n_age - n_hb),
             "Missing one or more screening covariates (see section 3)")
) %>%
  mutate(`% of previous step` = round(100 * n / lag(n), 1))
print(flow, width = Inf)

# Outcome groups within the analytic sample
group_tab <- data_clean %>%
  count(mal_anaemia_group) %>%
  mutate(pct = round(100 * n / sum(n), 1))
print(group_tab)

#-------------------------------------------------------------------------------
# 2. MISSINGNESS BY VARIABLE (among children aged 6-59 months)
#-------------------------------------------------------------------------------
miss_table <- function(df, vars) {
  vars <- intersect(vars, names(df))
  tibble(Variable = vars,
         n_missing = map_int(vars, ~ sum(is.na(df[[.x]]))),
         n_total   = nrow(df)) %>%
    mutate(pct_missing = round(100 * n_missing / n_total, 1)) %>%
    arrange(desc(n_missing))
}

all_vars <- c("malaria_child", "anaemia_child",
              "age_child", "gender_child", "residence", "birth_order",
              "underwt", "stunting", "wasting", "fever", "diarrhea", "cough",
              "deworming", "vitA", "iron_supp", "insurance",
              "age_mom", "marital_stats", "education", "anemia_mother",
              "zonal_region", "wealth_index", "hh_under5", "hh_size",
              "head_sex", "head_age", "water_source",
              "floor_material", "wall_material", "roof_material", "hh_owns_itn")

miss_all <- miss_table(mal_anemia, all_vars)
print(miss_all, n = Inf)

#-------------------------------------------------------------------------------
# 3. WHO WAS EXCLUDED AT THE ANALYTIC-SAMPLE STEP, AND WHY
#    (children with valid RDT + Hb but missing a screening covariate)
#-------------------------------------------------------------------------------
screen_vars <- c("stunting", "wasting", "underwt", "fever", "diarrhea", "cough",
                 "deworming", "anaemia_child", "malaria_child", "anemia_mother",
                 "zonal_region", "wealth_index", "hh_under5", "head_age")

valid_both <- mal_anemia %>%
  filter(!is.na(malaria_child), !is.na(anaemia_child))

excluded_step5 <- valid_both %>%
  filter(if_any(all_of(screen_vars), is.na))

cat("\nChildren excluded because of missing screening covariates:",
    nrow(excluded_step5), "\n")
print(miss_table(excluded_step5, screen_vars) %>% filter(n_missing > 0))
# (a child can be missing more than one variable, so these counts can sum to
#  more than the number excluded)

#-------------------------------------------------------------------------------
# 4. MODEL-SPECIFIC SAMPLES: eligible children, complete cases, children lost
#-------------------------------------------------------------------------------
model_spec <- list(
  "Model 1 (A vs D)" = list(
    groups = c("A_malaria_anaemia", "D_neither"),
    vars = c("age_child", "residence", "wealth_index", "education",
             "anemia_mother", "hh_owns_itn", "water_source", "hh_size",
             "stunting", "wasting")),
  "Model 2 (B vs D)" = list(
    groups = c("B_nonmalaria_anaemia", "D_neither"),
    vars = c("age_child", "residence", "zonal_region", "wealth_index",
             "education", "anemia_mother", "water_source", "hh_size",
             "stunting", "wasting", "fever", "deworming")),
  "Model 3 (A, B vs D)" = list(
    groups = c("A_malaria_anaemia", "B_nonmalaria_anaemia", "D_neither"),
    vars = c("age_child", "residence", "wealth_index", "education",
             "anemia_mother", "water_source", "hh_size", "stunting", "wasting"))
)

model_flow <- imap_dfr(model_spec, function(s, nm) {
  d <- data_clean %>% filter(mal_anaemia_group %in% s$groups)
  tibble(Model = nm,
         Eligible = nrow(d),
         Complete_cases = sum(complete.cases(d[, s$vars])),
         Lost = Eligible - Complete_cases,
         Pct_lost = round(100 * Lost / Eligible, 1))
})
print(model_flow)

# Which variables caused the loss in each model (eligible children only)
model_missing <- imap_dfr(model_spec, function(s, nm) {
  d <- data_clean %>% filter(mal_anaemia_group %in% s$groups)
  miss_table(d, s$vars) %>% mutate(Model = nm) %>% filter(n_missing > 0)
}) %>% select(Model, Variable, n_missing, n_total, pct_missing)
print(model_missing, n = Inf)

# Group C is never modelled
n_C <- sum(data_clean$mal_anaemia_group == "C_malaria_no_anaemia")
cat("\nGroup C (malaria without anaemia), not modelled:", n_C, "children\n")

# Cross-check against the fitted models (should match the 'Complete_cases' column)
cat("\nCross-check with the fitted models:\n")
if (exists("model_A_final")) cat("  Model 1 nobs:", nobs(model_A_final), "\n")
if (exists("model_B_final")) cat("  Model 2 nobs:", nobs(model_B_final), "\n")
if (exists("design_joint"))  cat("  Model 3 n   :", nrow(design_joint$variables), "\n")

#-------------------------------------------------------------------------------
# 5. SAVE THE TABLES
#-------------------------------------------------------------------------------
write.csv(flow,          "sample_flow.csv",          row.names = FALSE)
write.csv(group_tab,     "outcome_groups.csv",       row.names = FALSE)
write.csv(miss_all,      "missingness_by_variable.csv", row.names = FALSE)
write.csv(model_flow,    "model_sample_sizes.csv",   row.names = FALSE)
write.csv(model_missing, "model_missingness.csv",    row.names = FALSE)

#-------------------------------------------------------------------------------
# 6. READY-MADE SUMMARY SENTENCE (paste into the report after checking)
#-------------------------------------------------------------------------------
cat(sprintf(
  paste0("Of %s children in the Children's Recode, %s were aged 6-59 months and ",
         "%s had valid RDT and haemoglobin results. After excluding %s children with ",
         "missing screening covariates, %s children formed the analytic sample (%s). ",
         "The multivariable models used complete cases: %s of %s children in Model 1 ",
         "(%s lost), %s of %s in Model 2 (%s lost) and %s of %s in Model 3 (%s lost). ",
         "Group C (n = %s) was not modelled.\n"),
  format(n_all, big.mark = ","), format(n_age, big.mark = ","),
  format(n_both, big.mark = ","), format(n_both - n_analytic, big.mark = ","),
  format(n_analytic, big.mark = ","),
  paste(sprintf("%s: n = %s", group_tab$mal_anaemia_group, group_tab$n), collapse = "; "),
  format(model_flow$Complete_cases[1], big.mark = ","), format(model_flow$Eligible[1], big.mark = ","),
  paste0(model_flow$Pct_lost[1], "%"),
  format(model_flow$Complete_cases[2], big.mark = ","), format(model_flow$Eligible[2], big.mark = ","),
  paste0(model_flow$Pct_lost[2], "%"),
  format(model_flow$Complete_cases[3], big.mark = ","), format(model_flow$Eligible[3], big.mark = ","),
  paste0(model_flow$Pct_lost[3], "%"),
  n_C))
