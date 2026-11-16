#===============================================================================
# Malaria, Anaemia and the Double Burden among Children Under Five
# Tanzania 2022 DHS-MIS
# Author: Annabel Itaeli and Ayodele S. Babalola (revised)
# Revised: 28 Aug 2026
#
# This script:
#   0. Sets up correct file paths for the downloaded DHS recode bundle
#   1. Loads & merges the Children (KR) and Household Member (PR) recodes
#   2. Builds all analysis variables (with bug fixes flagged "# FIXED")
#   3. Builds the malaria/anaemia case-definition (Groups A-D) described in
#      the project's progress report:
#        A = malaria-related anaemia   (RDT+ & anaemic)
#        B = non-malaria-related anaemia (RDT- & anaemic)
#        C = malaria without anaemia   (RDT+ & not anaemic)
#        D = neither (reference group) (RDT- & not anaemic)
#   4. Produces descriptive tables (unweighted + survey-weighted)
#   5. Produces descriptive MAPS (cluster bubble maps + zone summary map)
#   6. Declares the complex survey design
#   7. Runs Rao-Scott chi-square tests
#   8. Runs bivariate modified Poisson regression (prevalence ratios)
#   9. Runs multivariable logistic regression for Model A and Model B
#  10. Compares the two models to identify SHARED vs SPECIFIC determinants
#===============================================================================

#-------------------------------------------------------------------------------
# 0. PACKAGES
#-------------------------------------------------------------------------------
# install.packages(c("expss","labelled","sjlabelled","survey","srvyr","dplyr",
#                     "tidyverse","haven","sf","ggplot2","broom","viridis",
#                     "geodata","scales"))

library(expss)      # tables that respect survey weights & labelled data
library(survey)     # complex survey design, svyglm, svychisq
library(srvyr)      # tidyverse-style wrapper around `survey`
library(dplyr)
library(tidyverse)
library(haven)      # read .DTA (Stata) files
library(labelled)
library(sjlabelled)
library(sf)         # read/handle the DHS GPS shapefile & polygons
library(ggplot2)
library(broom)       # tidy() model output, incl. svyglm
library(viridis)     # colour-blind-friendly map palettes
library(scales)

#-------------------------------------------------------------------------------
# 1. FILE PATHS  ----------------------------------------------------------------
#   Pattern confirmed by the analyst: each recode "XX" lives in its own folder
#   named  TZ<XX><phase>DT  (or, for the two geographic files, TZ<XX><phase>FL)
#   and the Stata file inside is named  TZ<XX><phase>FL.DTA
#
#   e.g. TZPR82FL.DTA  ->  .../TZ_2022_DHS_.../TZPR82DT/TZPR82FL.DTA
#-------------------------------------------------------------------------------

base_dir <- "C:/Users/babal/OneDrive/Desktop/PROJECT MAP ECR/DATA/TZ_2022_DHS_07202026_616_87810"


# Generic helper -- reuse this for ANY other recode you need later
# (default phase = "82"; TZFW81DT is the one exception that uses phase "81")
dhs_dta_path <- function(code, phase = "82", base = base_dir) {
  folder <- file.path(base, paste0("TZ", code, phase, "DT"))
  fname  <- paste0("TZ", code, phase, "FL.DTA")
  file.path(folder, fname)
}

# Quick reference for the folders you listed (standard DHS-8 recode codes;
# double-check FW/GC/GR/NR/SR against the "Recode Manual"/read-me that shipped
# with your download if you end up using them -- they are not used below):
#   BR = Births Recode            | CR = Couples' Recode
#   FW = Fieldworker/paradata     | GC = Geospatial Covariates (points, .DTA)
#   GE = GPS/Geographic Recode (cluster shapefile, .shp -- NOT .DTA)
#   GR = (verify)                 | HR = Household Recode
#   IR = Individual (women's) Recode | KR = Children's Recode  <-- used
#   MR = Men's Recode             | NR = (verify)
#   PR = Household Member Recode  <-- used (has hml35 malaria RDT & hc56 Hb)
#   SR = (verify)

kr_path <- dhs_dta_path("KR")   # Children's recode
pr_path <- dhs_dta_path("PR")   # Household member recode (malaria + Hb)

# The GE (GPS) recode is a shapefile bundle, not a .DTA -- folder has NO "DT" suffix
ge_folder <- file.path(base_dir, "TZGE81FL")

stopifnot(file.exists(kr_path), file.exists(pr_path))

#-------------------------------------------------------------------------------
# 2. LOAD & MERGE
#-------------------------------------------------------------------------------
child_data <- read_dta(kr_path)  # Children Recode (KR): one row per living child
hh_data    <- read_dta(pr_path)  # Household Member Recode (PR): malaria/Hb testing

dhs_TZ <- child_data %>%
  left_join(hh_data,
            by = c("v001" = "hv001",  # cluster number
                   "v002" = "hv002",  # household number
                   "b16"  = "hvidx")) # child's line number in household

#-------------------------------------------------------------------------------
# 3. SELECT VARIABLES OF INTEREST
#-------------------------------------------------------------------------------
ml_an_df <- dhs_TZ %>%
  select(caseid, v001, v002,           # IDs + cluster (needed for maps/design)
         v010, b19, b4, v102, hv024, bord,
         m19, hc70, hc71, hc72,
         v481, h22, h11, h31, h33, h42, h43,
         hc56,                          # child Hb (adjusted, g/dL x 10)
         ## Mother-related
         v012, v501, v106, v731, v457, v463a,
         ## Household-related
         hv012,                         # NEW: household size (for family-size var)
         v137, v151, v152, hv270, hv201,
         hv213, hv214, hv215,           # NEW: floor / wall / roof material
         hv227,                         # NEW: household owns >=1 mosquito net
         ## Malaria/community
         hv024, hml35, hml20,           # NEW: hml20 = # of ITNs in household
         hml12,                         # FIXED: needed for sleep_ITN below (was missing)
         ## Random-effect / design variables
         v005, v021, v022)

#-------------------------------------------------------------------------------
# 4. VARIABLE RECODING
#-------------------------------------------------------------------------------

# ---- Child age --------------------------------------------------------------
mal_anemia <- ml_an_df %>% filter(b19 > 5)   # keep 6-59 months (malaria/Hb-eligible)

mal_anemia <- mal_anemia %>% mutate(age_child = case_when(
  b19 >= 6  & b19 <= 11 ~ "Infant",
  b19 >= 12 & b19 <= 23 ~ "Year-1",
  b19 >= 24 & b19 <= 35 ~ "Year-2",
  b19 >= 36 & b19 <= 47 ~ "Year-3",
  b19 >= 48 & b19 <= 59 ~ "Year-4"))
mal_anemia$age_child <- factor(mal_anemia$age_child,
                                levels = c("Infant","Year-1","Year-2","Year-3","Year-4"))

# ---- Gender -------------------------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(gender_child = case_when(
  b4 == 1 ~ "male", b4 == 2 ~ "female"))

# ---- Residence -----------------------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(residence = case_when(
  v102 == 1 ~ "urban", v102 == 2 ~ "rural"))

# ---- Birth order -----------------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(birth_order = case_when(
  bord == 1 ~ "first", bord == 2 ~ "second", bord == 3 ~ "third",
  bord == 4 ~ "fourth", bord == 5 ~ "fifth", bord >= 6 ~ "Sixth&above"))

# ---- Stunting (HAZ) ---------------------------------------------------------
mal_anemia <- mal_anemia %>%
  mutate(haz = hc70 / 100,
         stunting = case_when(
           hc70 == 9998 ~ NA_character_,
           haz < -6 | haz > 6 ~ NA_character_,
           haz < -2 ~ "stunted",
           haz >= -2 ~ "normal",
           TRUE ~ NA_character_))

# ---- Wasting (WHZ) -----------------------------------------------------------
mal_anemia <- mal_anemia %>%
  mutate(whz = hc72 / 100,
         wasting = case_when(
           hc72 %in% c(9996, 9998, 9999) ~ NA_character_,
           whz < -5 | whz > 5 ~ NA_character_,
           whz < -2 ~ "wasted",
           whz >= -2 ~ "normal",
           TRUE ~ NA_character_))

# ---- Underweight (WAZ) -------------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(underwt = case_when(
  hc71 == 9998 ~ NA_character_,
  hc71 < -200 ~ "underweight",
  hc71 >= -200 & hc71 <= 200 ~ "normal",
  hc71 > 200 ~ "overweight",
  TRUE ~ NA_character_))

# ---- Fever ---------------------------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(fever = case_when(
  h22 == 0 ~ "no", h22 == 1 ~ "yes", h22 == 8 ~ NA_character_))

# ---- Diarrhoea -----------------------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(diarrhea = case_when(
  h11 == 0 ~ "no", h11 %in% c(1, 2) ~ "yes", h11 == 8 ~ NA_character_))

# ---- Cough ---------------------------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(cough = case_when(
  h31 == 0 ~ "no", h31 %in% c(1, 2) ~ "yes", h31 == 8 ~ NA_character_))

# ---- Deworming -------------------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(deworming = case_when(
  h43 == 0 ~ "no", h43 == 1 ~ "yes", h43 == 8 ~ NA_character_))

# ---- Vitamin A ---------------------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(vitA = case_when(
  h33 == 0 ~ "no", h33 >= 1 & h33 <= 3 ~ "yes", h33 == 8 ~ NA_character_))

# ---- Iron supplementation ---------------------------------------------------
# FIXED: original code referenced h33 (Vitamin A) here by mistake -> now uses h42
mal_anemia <- mal_anemia %>% mutate(iron_supp = case_when(
  h42 == 0 ~ "no", h42 == 1 ~ "yes", h42 == 8 ~ NA_character_))

# ---- Child anaemia status ----------------------------------------------------
# hc56 is stored in g/dL x 10 (e.g. 109 = 10.9 g/dL); WHO cut-off for 6-59 mo = 11.0 g/dL
mal_anemia <- mal_anemia %>% mutate(anaemia_child = case_when(
  hc56 %in% c(9994, 9995, 9997, 9998, 9999) ~ NA_character_,  # DHS missing/refused codes
  hc56 < 110 ~ "Anaemia",
  hc56 >= 110 ~ "Normal",
  TRUE ~ NA_character_))

# ---- Child malaria status (RDT) ----------------------------------------------
# FIXED: original "unique()" inspected hml32 (a variable that was never even
# selected) instead of hml35 (the RDT result actually used to build the label)
mal_anemia <- mal_anemia %>% mutate(malaria_child = case_when(
  hml35 == 0 ~ "negative", hml35 == 1 ~ "positive", hml35 >= 6 ~ NA_character_))

# ---- Health insurance --------------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(insurance = case_when(
  v481 == 0 ~ "no", v481 == 1 ~ "yes", TRUE ~ NA_character_))

# ---- Mother's age --------------------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(age_mom = case_when(
  v012 >= 15 & v012 <= 19 ~ "adolescent",
  v012 >= 20 & v012 <= 24 ~ "young adult",
  v012 >= 25 & v012 <= 29 ~ "late 20s",
  v012 >= 30 & v012 <= 34 ~ "30-34",
  v012 >= 35 & v012 <= 39 ~ "35-39",
  v012 >= 40 & v012 <= 44 ~ "40-44",
  v012 >= 45 & v012 <= 49 ~ "45-49"))

# ---- Marital status -----------------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(marital_stats = case_when(
  v501 == 0 ~ "Never married", v501 == 3 ~ "Widowed",
  v501 %in% c(1, 2) ~ "Married or living together",
  v501 %in% c(4, 5) ~ "Divorced or separated"))

# ---- Mother's education --------------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(education = case_when(
  v106 == 0 ~ "No education", v106 == 1 ~ "Primary", v106 == 2 ~ "Secondary",
  v106 == 3 ~ "More than sec", v106 == 8 ~ NA_character_))

# ---- Mother's anaemia status ---------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(anemia_mother = case_when(
  v457 == 1 ~ "Severe", v457 == 2 ~ "Moderate", v457 == 3 ~ "Mild",
  v457 == 4 ~ "Not anemic", TRUE ~ NA_character_))

# ---- Child slept under ITN last night ------------------------------------------
mal_anemia <- mal_anemia %>% mutate(sleep_ITN = case_when(
  hml12 == 0 ~ "No", hml12 == 1 ~ "Only ITN", hml12 == 2 ~ "Both untreated&ITN",
  hml12 == 3 ~ "Only untreated", TRUE ~ NA_character_))
# NB: hml12 was in the original variable list' comments but never select()-ed --
# if you need this variable, add hml12 to the select() step in Section 3 first.

# ---- Zonal region (DHS zone groupings) -----------------------------------------
mal_anemia <- mal_anemia %>% mutate(zonal_region = case_when(
  hv024 %in% c(17, 18, 19, 20, 24, 25) ~ "lake",
  hv024 %in% c(2, 3, 4) ~ "northern",
  hv024 %in% c(1, 13, 21) ~ "central",
  hv024 %in% c(14, 16) ~ "western",
  hv024 %in% c(10, 11, 22, 26) ~ "southern highlands",
  hv024 %in% c(5, 6, 7) ~ "eastern",
  hv024 %in% c(12, 15, 23) ~ "southwest",
  hv024 %in% c(8, 9) ~ "southern",
  hv024 %in% c(51, 52, 53, 54, 55) ~ "Zanzibar",
  TRUE ~ NA_character_))

# ---- Wealth index --------------------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(wealth_index = case_when(
  hv270 == 1 ~ "Poorest", hv270 == 2 ~ "Poorer", hv270 == 3 ~ "Middle",
  hv270 == 4 ~ "Richer", hv270 == 5 ~ "Richest", TRUE ~ NA_character_))
mal_anemia$wealth_index <- factor(mal_anemia$wealth_index,
                                   levels = c("Poorest","Poorer","Middle","Richer","Richest"))

# ---- Number of under-5s in household -------------------------------------------
mal_anemia <- mal_anemia %>% mutate(hh_under5 = case_when(
  v137 == 1 ~ "One", v137 == 2 ~ "Two", v137 >= 3 ~ "Three and more",
  TRUE ~ NA_character_))

# ---- NEW: Household size (family size) -----------------------------------------
mal_anemia <- mal_anemia %>% mutate(hh_size = case_when(
  hv012 < 5 ~ "<5",
  hv012 >= 5 & hv012 <= 10 ~ "5-10",
  hv012 > 10 ~ ">10",
  TRUE ~ NA_character_))
mal_anemia$hh_size <- factor(mal_anemia$hh_size, levels = c("<5","5-10",">10"))

# ---- Sex of household head -------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(head_sex = case_when(
  v151 == 1 ~ "male", v151 == 2 ~ "female"))

# ---- Age of household head -------------------------------------------------
mal_anemia <- mal_anemia %>% mutate(head_age = case_when(
  v152 == 98 ~ NA_character_,
  v152 >= 15 & v152 <= 29 ~ "Less than 30",
  v152 >= 30 & v152 <= 39 ~ "Thirties",
  v152 >= 40 & v152 <= 49 ~ "Forties",
  v152 >= 50 & v152 <= 59 ~ "Fifties",
  v152 >= 60 ~ "Sixty & more",
  TRUE ~ NA_character_))

# ---- Source of drinking water -----------------------------------------------
mal_anemia <- mal_anemia %>% mutate(water_source = case_when(
  hv201 %in% c(10,11,12,13,14,20,21,31,41,51,71,72,73) ~ "improved",
  hv201 %in% c(30,32,40,42,61,62,96) ~ "unimproved",
  TRUE ~ NA_character_))

# ---- NEW: Housing materials & ITN ownership ------------------------------------
# IMPORTANT: DHS value codes for floor/wall/roof differ slightly by country/round.
# Before trusting the classification below, inspect the actual codes/labels with:
#    sjlabelled::get_labels(hh_data$hv213, values = "p")   # floor
#    sjlabelled::get_labels(hh_data$hv214, values = "p")   # wall
#    sjlabelled::get_labels(hh_data$hv215, values = "p")   # roof
# and adjust the code lists below if Tanzania 2022 uses different numbering.
mal_anemia <- mal_anemia %>% mutate(
  floor_material = case_when(
    hv213 %in% c(11, 12) ~ "unimproved",                 # earth/sand, dung
    hv213 %in% c(96) ~ NA_character_,                    # other
    !is.na(hv213) ~ "improved",
    TRUE ~ NA_character_),
  wall_material = case_when(
    hv214 %in% c(11:17, 96) ~ "unimproved",               # no wall/cane/dirt/etc.
    !is.na(hv214) ~ "improved",
    TRUE ~ NA_character_),
  roof_material = case_when(
    hv215 %in% c(11:17, 96) ~ "unimproved",               # no roof/thatch/mud/etc.
    !is.na(hv215) ~ "improved",
    TRUE ~ NA_character_),
  hh_owns_itn = case_when(
    hv227 == 0 ~ "no", hv227 == 1 ~ "yes", TRUE ~ NA_character_)
)

## Inspect the columns with NA's
colSums(is.na(mal_anemia))

#-------------------------------------------------------------------------------
# 5. FINAL ANALYTIC SAMPLE (drop rows missing on key covariates)
#-------------------------------------------------------------------------------
data_clean <- mal_anemia %>%
  filter(
    !is.na(stunting), !is.na(wasting), !is.na(underwt),
    !is.na(fever), !is.na(diarrhea), !is.na(cough), !is.na(deworming),
    !is.na(anaemia_child), !is.na(malaria_child), !is.na(anemia_mother),
    !is.na(zonal_region), !is.na(wealth_index), !is.na(hh_under5),
    !is.na(head_age)
  )

nrow(data_clean)  # should match / exceed the n = 4,166 in the preliminary report

#===============================================================================
# 6. OUTCOME CONSTRUCTION: THE FOUR-GROUP MALARIA-ANAEMIA CLASSIFICATION
#    (this operationalises Section 4.2 of the progress report)
#===============================================================================
data_clean <- data_clean %>%
  mutate(
    mal_anaemia_group = case_when(
      malaria_child == "positive" & anaemia_child == "Anaemia" ~ "A_malaria_anaemia",
      malaria_child == "negative" & anaemia_child == "Anaemia" ~ "B_nonmalaria_anaemia",
      malaria_child == "positive" & anaemia_child == "Normal"  ~ "C_malaria_no_anaemia",
      malaria_child == "negative" & anaemia_child == "Normal"  ~ "D_neither"),
    mal_anaemia_group = factor(mal_anaemia_group,
                                levels = c("D_neither","A_malaria_anaemia",
                                           "B_nonmalaria_anaemia","C_malaria_no_anaemia")),
    # binary flags used throughout
    any_malaria  = if_else(malaria_child == "positive", 1L, 0L),
    any_anaemia  = if_else(anaemia_child == "Anaemia", 1L, 0L),
    double_burden = if_else(mal_anaemia_group == "A_malaria_anaemia", 1L, 0L)
  )

table(data_clean$mal_anaemia_group)
prop.table(table(data_clean$mal_anaemia_group)) * 100

#===============================================================================
# 7. DESCRIPTIVE STATISTICS
#===============================================================================

# ---- 7a. SURVEY-WEIGHTED descriptive summary (replaces the unweighted version) ----
options(survey.lonely.psu = "adjust")
mk_summary <- function(df, group_var) {
  df %>%
    mutate(wt = v005 / 1e6) %>%
    as_survey_design(ids = v021, strata = v022, weights = wt, nest = TRUE) %>%
    group_by({{ group_var }}) %>%
    summarise(
      n_unw     = unweighted(n()),
      w_total   = unweighted(sum(wt)),
      p_malaria = survey_mean(any_malaria,   vartype = "ci") * 100,
      p_anaemic = survey_mean(any_anaemia,   vartype = "ci") * 100,
      p_double  = survey_mean(double_burden, vartype = "ci") * 100) %>%
    mutate(n = n_unw, percent = 100 * w_total / sum(w_total)) %>%
    select(-n_unw, -w_total)
}

vars_to_summarise <- c("zonal_region","age_child","gender_child","residence",
                        "birth_order","underwt","stunting","wasting","fever",
                        "diarrhea","cough","deworming","age_mom","education",
                        "anemia_mother","marital_stats","wealth_index","hh_under5",
                        "hh_size","head_age","water_source","floor_material",
                        "wall_material","roof_material","hh_owns_itn")

desc_tables <- lapply(vars_to_summarise, function(v) mk_summary(data_clean, !!sym(v)))
names(desc_tables) <- vars_to_summarise
# e.g. view one: desc_tables[["wealth_index"]]

# ---- 7a-ii. Faceted descriptive plots, one facet-wrap plot per classification level ----
# Goal: 4 plots (Child / Maternal / Household / Community), each showing every
# variable in that level as its own facet panel, for the preliminary report.

# 1) Tag every covariate with its analytic level (used for Section 4.3 of the report)
var_classification <- tibble::tribble(
  ~variable,          ~level,
  "age_child",        "Child",
  "gender_child",      "Child",
  "birth_order",       "Child",
  "underwt",           "Child",
  "stunting",          "Child",
  "wasting",           "Child",
  "fever",             "Child",
  "diarrhea",          "Child",
  "cough",             "Child",
  "deworming",         "Child",
  "age_mom",           "Maternal",
  "education",         "Maternal",
  "anemia_mother",     "Maternal",
  "marital_stats",     "Maternal",
  "wealth_index",      "Household",
  "hh_under5",         "Household",
  "hh_size",           "Household",
  "head_age",          "Household",
  "water_source",      "Household",
  "floor_material",    "Household",
  "wall_material",     "Household",
  "roof_material",     "Household",
  "hh_owns_itn",       "Household",
  "zonal_region",      "Community",
  "residence",         "Community"
)
stopifnot(all(vars_to_summarise %in% var_classification$variable))  # catches un-tagged vars

# 2) Nicer facet-panel titles than the raw variable names
var_labels <- c(
  age_child = "Child's age", gender_child = "Sex of child",
  birth_order = "Birth order", underwt = "Underweight (WAZ)",
  stunting = "Stunting (HAZ)", wasting = "Wasting (WHZ)",
  fever = "Fever (2 wks)", diarrhea = "Diarrhoea (2 wks)",
  cough = "Cough (2 wks)", deworming = "Deworming",
  age_mom = "Mother's age", education = "Mother's education",
  anemia_mother = "Mother's anaemia status", marital_stats = "Marital status",
  wealth_index = "Wealth index", hh_under5 = "Under-5s in household",
  hh_size = "Household size", head_age = "Household head's age",
  water_source = "Drinking water source", floor_material = "Floor material",
  wall_material = "Wall material", roof_material = "Roof material",
  hh_owns_itn = "Household owns \u2265 1 net", zonal_region = "DHS zone",
  residence = "Residence")

# 3) Sensible (non-alphabetical) ordering for the ordinal/ordered variables
level_orders <- list(
  age_child     = c("Infant","Year-1","Year-2","Year-3","Year-4"),
  birth_order   = c("first","second","third","fourth","fifth","Sixth&above"),
  age_mom       = c("adolescent","young adult","late 20s","30-34","35-39","40-44","45-49"),
  education     = c("No education","Primary","Secondary","More than sec"),
  marital_stats = c("Never married","Married or living together",
                     "Divorced or separated","Widowed"),
  anemia_mother = c("Not anemic","Mild","Moderate","Severe"),
  wealth_index  = c("Poorest","Poorer","Middle","Richer","Richest"),
  hh_under5     = c("One","Two","Three and more"),
  hh_size       = c("<5","5-10",">10"),
  head_age      = c("Less than 30","Thirties","Forties","Fifties","Sixty & more")
)

# 4) Stack every table in `desc_tables` into one long data frame, tag level/labels,
#    and apply the ordering from (3) variable-by-variable
desc_long <- purrr::imap_dfr(desc_tables, function(df, var_name) {
  df %>%
    rename(category = 1) %>%                # 1st column = the grouping levels
    mutate(category = as.character(category), variable = var_name) %>%
    select(variable, category, n, percent, p_malaria, p_anaemic, p_double)
}) %>%
  left_join(var_classification, by = "variable") %>%
  mutate(variable_label = dplyr::recode(variable, !!!var_labels))

desc_long <- desc_long %>%
  split(.$variable) %>%
  purrr::map_dfr(function(df) {
    ord <- level_orders[[df$variable[1]]]
    df$category <- if (!is.null(ord)) factor(df$category, levels = ord) else factor(df$category)
    df
  })

# 5) Re-usable facet-wrap plotting function: one plot per LEVEL, one panel per VARIABLE
level_colours <- c(Child = "#2C7BB6", Maternal = "#7B3294",
                    Household = "#D7791F", Community = "#1B9E77")

plot_level_facet <- function(data, level_name, outcome = "p_anaemic",
                              outcome_label = "Anaemia prevalence (%)",
                              ncol = 3) {
  df <- data %>% filter(level == level_name)
  ggplot(df, aes(x = category, y = .data[[outcome]])) +
    geom_col(fill = level_colours[[level_name]], width = 0.7) +
    geom_text(aes(label = paste0(round(.data[[outcome]], 1), "%")),
              vjust = -0.4, size = 3.8) +
    facet_wrap(~ variable_label, scales = "free_x", ncol = ncol) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = paste0(level_name, "-level characteristics"),
         subtitle = outcome_label,
         x = NULL, y = outcome_label,
         caption = "Source: 2022 TDHS-MIS (KR & PR recodes)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))
}

levels_list <- c("Child", "Maternal", "Household", "Community")

# ---- Anaemia prevalence (primary outcome) -- 4 plots, one per level ------------
anaemia_facets <- purrr::map(levels_list, ~ plot_level_facet(
  desc_long, .x, outcome = "p_anaemic", outcome_label = "Anaemia prevalence (%)"))
names(anaemia_facets) <- levels_list
purrr::iwalk(anaemia_facets, ~ ggsave(paste0("facet_anaemia_", tolower(.y), ".png"),
                                       .x, width = 10, height = 7, dpi = 300))

# ---- Malaria prevalence -- 4 plots, one per level -------------------------------
malaria_facets <- purrr::map(levels_list, ~ plot_level_facet(
  desc_long, .x, outcome = "p_malaria", outcome_label = "Malaria prevalence (%)"))
names(malaria_facets) <- levels_list
purrr::iwalk(malaria_facets, ~ ggsave(paste0("facet_malaria_", tolower(.y), ".png"),
                                       .x, width = 10, height = 7, dpi = 300))

# ---- Double-burden (malaria-related anaemia) prevalence -- 4 plots, one per level ---
double_facets <- purrr::map(levels_list, ~ plot_level_facet(
  desc_long, .x, outcome = "p_double", outcome_label = "Double-burden prevalence (%)"))
names(double_facets) <- levels_list
purrr::iwalk(double_facets, ~ ggsave(paste0("facet_doubleburden_", tolower(.y), ".png"),
                                      .x, width = 10, height = 7, dpi = 300))

# ---- 6) COMBINED version: malaria + anaemia + double-burden as grouped bars ----
#    within each panel, one plot per level (4 plots total instead of 12) --
#    handy as a single compact "small multiples" summary figure per level.
desc_long_combined <- desc_long %>%
  select(variable, variable_label, level, category, p_malaria, p_anaemic, p_double) %>%
  pivot_longer(cols = c(p_malaria, p_anaemic, p_double),
               names_to = "outcome", values_to = "prevalence") %>%
  mutate(outcome = dplyr::recode(outcome,
                                  p_malaria = "Malaria",
                                  p_anaemic = "Anaemia",
                                  p_double  = "Double burden"),
         outcome = factor(outcome, levels = c("Malaria", "Anaemia", "Double burden")))

outcome_colours <- c(Malaria = "#2166AC", Anaemia = "#B2182B", "Double burden" = "#762A83")

plot_level_facet_combined <- function(data, level_name, ncol = 3) {
  df <- data %>% filter(level == level_name)
  ggplot(df, aes(x = category, y = prevalence, fill = outcome)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_text(aes(label = round(prevalence, 1)), position = position_dodge(width = 0.8),
              vjust = -0.3, size = 3.4, check_overlap = TRUE, show.legend = FALSE) +
    facet_wrap(~ variable_label, scales = "free_x", ncol = ncol) +
    scale_fill_manual(values = outcome_colours, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = paste0(level_name, "-level characteristics"),
         subtitle = "Survey-weighted prevalence (%), grouped per category",
         x = NULL, y = "Weighted prevalence (%)",
         caption = "Source: 2022 TDHS-MIS (KR & PR recodes)") +
    theme_minimal(base_size = 15) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          strip.text = element_text(face = "bold", size = 14),
          plot.title = element_text(face = "bold"),
          legend.position = "top", legend.text = element_text(size = 14))
}
combined_facets <- purrr::map(levels_list, ~ plot_level_facet_combined(desc_long_combined, .x))
names(combined_facets) <- levels_list
fig_heights <- c(Child = 9, Maternal = 6.5, Household = 9, Community = 5)
purrr::iwalk(combined_facets, ~ ggsave(paste0("facet_combined_", tolower(.y), ".png"),
                                       .x, width = 10, height = fig_heights[[.y]], dpi = 300))

# View any single plot in the RStudio Plots pane, e.g.:
print(anaemia_facets[["Household"]])
print(anaemia_facets[["Child"]])
print(combined_facets[["Child"]])   # all 3 outcomes side-by-side, Child-level vars
# 16 PNGs total are saved to the working directory:
#   - 12 single-outcome plots (4 levels x 3 outcomes: facet_<outcome>_<level>.png)
#   - 4 combined plots (one per level, all 3 outcomes grouped: facet_combined_<level>.png)
# ready to drop straight into the preliminary/descriptive-results section
# of the progress report.

# ---- 7b. Survey-weighted prevalence (the numbers that belong in the final report) ---
full_design <- data_clean %>%
  mutate(wt = v005 / 1e6) %>%
  as_survey_design(ids = v021, strata = v022, weights = wt, nest = TRUE)

weighted_overall <- full_design %>%
  summarise(
    malaria_prev = survey_mean(any_malaria, vartype = "ci") * 100,
    anaemia_prev = survey_mean(any_anaemia, vartype = "ci") * 100,
    double_burden_prev = survey_mean(double_burden, vartype = "ci") * 100)
print(weighted_overall)

weighted_by_zone <- full_design %>%
  group_by(zonal_region) %>%
  summarise(
    n = unweighted(n()),
    malaria_prev = survey_mean(any_malaria, vartype = "ci") * 100,
    anaemia_prev = survey_mean(any_anaemia, vartype = "ci") * 100,
    double_burden_prev = survey_mean(double_burden, vartype = "ci") * 100)
print(weighted_by_zone)

tab1 <- full_design %>% summarise(
  N = unweighted(n()),
  n_malaria = unweighted(sum(any_malaria)), n_anaemia = unweighted(sum(any_anaemia)),
  n_double  = unweighted(sum(double_burden)),
  malaria = survey_mean(any_malaria,   vartype = "ci", proportion = TRUE) * 100,
  anaemia = survey_mean(any_anaemia,   vartype = "ci", proportion = TRUE) * 100,
  double  = survey_mean(double_burden, vartype = "ci", proportion = TRUE) * 100)
print(as.data.frame(tab1))

tab2 <- full_design %>% group_by(mal_anaemia_group) %>%
  summarise(n = unweighted(n()), pct = survey_prop(vartype = "ci", proportion = TRUE) * 100)
print(as.data.frame(tab2))

full_design %>% filter(any_anaemia == 1) %>%
  summarise(malaria_among_anaemic = survey_mean(any_malaria, vartype = "ci", proportion = TRUE) * 100)
full_design %>% group_by(any_malaria) %>%
  summarise(anaemia_prev = survey_mean(any_anaemia, vartype = "ci", proportion = TRUE) * 100)


weighted_by_wealth <- full_design %>%
  group_by(wealth_index) %>%
  summarise(
    n = unweighted(n()),
    malaria_prev = survey_mean(any_malaria, vartype = "ci") * 100,
    anaemia_prev = survey_mean(any_anaemia, vartype = "ci") * 100)
print(weighted_by_wealth)

#===============================================================================
# 8. MAPS  (descriptive spatial analysis)
#===============================================================================

# ---- 8a. Load the DHS cluster GPS shapefile ------------------------------------
ge_shp_file <- list.files(ge_folder, pattern = "\\.shp$", full.names = TRUE)
stopifnot(length(ge_shp_file) == 1)   # fails loudly if the folder has 0 or >1 .shp

tz_clusters <- st_read(ge_shp_file)
# Typical DHS GE fields: DHSCLUST, LATNUM, LONGNUM, ADM1NAME, URBAN_RURA ...
# Clusters with missing GPS are coded LATNUM = LONGNUM = 0 -- drop these.
tz_clusters <- tz_clusters %>% filter(LATNUM != 0, LONGNUM != 0)

# ---- 8b. Cluster-level prevalence (unweighted; DHS clusters are ~equal-sized) --
cluster_summary <- data_clean %>%
  group_by(v001) %>%
  summarise(
    n_children      = n(),
    malaria_prev    = 100 * mean(malaria_child == "positive"),
    anaemia_prev    = 100 * mean(anaemia_child == "Anaemia"),
    double_burden_prev = 100 * mean(double_burden == 1),
    .groups = "drop")

map_data <- tz_clusters %>%
  left_join(cluster_summary, by = c("DHSCLUST" = "v001")) %>%
  filter(!is.na(n_children))   # keep only clusters with surveyed under-5 data

# ---- 8c. Bubble map: malaria prevalence by cluster -----------------------------
map_malaria <- ggplot(map_data) +
  geom_sf(aes(size = n_children, color = malaria_prev), alpha = 0.8) +
  scale_color_viridis_c(option = "C", name = "Malaria\nprevalence (%)") +
  scale_size_continuous(name = "Children\ntested", range = c(1, 6)) +
  labs(title = "Malaria prevalence among children 6-59 months",
       subtitle = "Cluster-level RDT positivity, Tanzania DHS-MIS 2022",
       caption = "Source: 2022 TDHS-MIS (PR & GE recodes)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right", panel.grid = element_blank())
print(map_malaria)
ggsave("map_malaria_prevalence.png", map_malaria, width = 8, height = 8, dpi = 300)

# ---- 8d. Bubble map: anaemia prevalence by cluster -----------------------------
map_anaemia <- ggplot(map_data) +
  geom_sf(aes(size = n_children, color = anaemia_prev), alpha = 0.8) +
  scale_color_viridis_c(option = "D", name = "Anaemia\nprevalence (%)") +
  scale_size_continuous(name = "Children\ntested", range = c(1, 6)) +
  labs(title = "Anaemia prevalence among children 6-59 months",
       subtitle = "Cluster-level Hb < 11 g/dL, Tanzania DHS-MIS 2022",
       caption = "Source: 2022 TDHS-MIS (PR & GE recodes)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right", panel.grid = element_blank())
print(map_anaemia)
ggsave("map_anaemia_prevalence.png", map_anaemia, width = 8, height = 8, dpi = 300)

# ---- 8e. Bubble map: double-burden (malaria-related anaemia) prevalence --------
map_double <- ggplot(map_data) +
  geom_sf(aes(size = n_children, color = double_burden_prev), alpha = 0.8) +
  scale_color_viridis_c(option = "B", name = "Double burden\nprevalence (%)") +
  scale_size_continuous(name = "Children\ntested", range = c(1, 6)) +
  labs(title = "Malaria-related anaemia (double burden)",
       subtitle = "Children RDT-positive AND anaemic, Tanzania DHS-MIS 2022",
       caption = "Source: 2022 TDHS-MIS (PR & GE recodes)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right", panel.grid = element_blank())
print(map_double)
ggsave("map_double_burden_prevalence.png", map_double, width = 8, height = 8, dpi = 300)

# ---- 8f. Zone-level bar chart (same outcome colours as Figures 5-8) ----
zone_bar <- weighted_by_zone %>%
  select(zonal_region, malaria_prev, anaemia_prev, double_burden_prev) %>%
  pivot_longer(-zonal_region, names_to = "indicator", values_to = "value") %>%
  mutate(indicator = factor(dplyr::recode(indicator, malaria_prev = "Malaria",
                                          anaemia_prev = "Anaemia", double_burden_prev = "Double burden"),
                            levels = c("Malaria", "Anaemia", "Double burden")))

ggplot(zone_bar, aes(x = reorder(zonal_region, value), y = value, fill = indicator)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = outcome_colours, name = NULL) +
  labs(x = NULL, y = "Weighted prevalence (%)", caption = "Source: 2022 TDHS-MIS") +
  theme_minimal(base_size = 16) + theme(legend.position = "top")

ggsave("zone_prevalence_barchart.png", width = 8, height = 6.5, dpi = 300)

# ---- 8g. Zone-level CHOROPLETH map (filled polygons, not just cluster points) ----
#
# The GE file only gives cluster POINTS -- there are no zone/region polygon
# boundaries anywhere in your downloaded DHS bundle, so a filled zone map
# needs an external administrative boundary shapefile. Two ways to get one:
#
#   OPTION A (most accurate, no name-matching needed): download Tanzania's
#   official "Subnational Boundaries" shapefile from the DHS Spatial Data
#   Repository (https://spatialdata.dhsprogram.com/boundaries/ -> Tanzania
#   -> 2022 -> download). It is built by DHS to match hv024 exactly. If you
#   have it, set `sdr_shp_path` to that .shp file and read it with st_read()
#   instead of the GADM block below, then skip straight to "DISSOLVE".
#
#   OPTION B (used below): auto-download Tanzania's 26 mainland + 5 Zanzibar
#   region (ADM1) boundaries from GADM via the `geodata` package. Needs
#   internet access on first run. The crosswalk below (hv024 code -> region
#   name -> zone) has been cross-checked against the exact code ranges
#   already used to build `zonal_region` in Section 4, so the ZONE
#   assignment is correct; the only remaining risk is that GADM's spelling
#   of a region name doesn't exactly match -- the diagnostic printout below
#   will flag any such mismatches so you can patch `manual_fix`.

# install.packages("geodata")   # only needed once (also installs `terra`)
library(geodata)
library(terra)

# ---- Validated crosswalk: hv024 code -> region name -> DHS zone ---------------
region_zone_lookup <- tibble::tribble(
  ~hv024, ~region_name,       ~zone,
   1,     "Dodoma",           "central",
   2,     "Arusha",           "northern",
   3,     "Kilimanjaro",      "northern",
   4,     "Tanga",            "northern",
   5,     "Morogoro",         "eastern",
   6,     "Pwani",            "eastern",
   7,     "Dar es Salaam",    "eastern",
   8,     "Lindi",            "southern",
   9,     "Mtwara",           "southern",
  10,     "Ruvuma",           "southern highlands",
  11,     "Iringa",           "southern highlands",
  12,     "Mbeya",            "southwest",
  13,     "Singida",          "central",
  14,     "Tabora",           "western",
  15,     "Rukwa",            "southwest",
  16,     "Kigoma",           "western",
  17,     "Shinyanga",        "lake",
  18,     "Kagera",           "lake",
  19,     "Mwanza",           "lake",
  20,     "Mara",             "lake",
  21,     "Manyara",          "central",
  22,     "Njombe",           "southern highlands",
  23,     "Katavi",           "southwest",
  24,     "Simiyu",           "lake",
  25,     "Geita",            "lake",
  26,     "Songwe",           "southern highlands",
  51,     "Kaskazini Unguja", "Zanzibar",
  52,     "Kusini Unguja",    "Zanzibar",
  53,     "Mjini Magharibi",  "Zanzibar",
  54,     "Kaskazini Pemba",  "Zanzibar",
  55,     "Kusini Pemba",     "Zanzibar"
)
# NB: if you ever edit the zonal_region case_when() in Section 4, update this
# table too so the map and the tables stay consistent.

# ---- Download Tanzania ADM1 (region) polygons ----------------------------------
tza_regions <- geodata::gadm(country = "TZA", level = 1, path = tempdir()) %>%
  sf::st_as_sf()

# ---- Diagnostic: check name-matching BEFORE trusting the map ------------------
norm_name <- function(x) tolower(trimws(gsub("[-_]", " ", x)))

gadm_names <- sort(unique(tza_regions$NAME_1))
dhs_names  <- sort(unique(region_zone_lookup$region_name))
cat("GADM region names (n =", length(gadm_names), "):\n"); print(gadm_names)
cat("\nDHS region names from the crosswalk (n =", length(dhs_names), "):\n"); print(dhs_names)
cat("\nGADM names with NO match in the DHS crosswalk -- fix these in `manual_fix`:\n")
unmatched <- gadm_names[!norm_name(gadm_names) %in% norm_name(dhs_names)]
print(unmatched)

# Populate this if the diagnostic above lists any mismatches, e.g.:
# e.g. manual_fix <- c("Zanzibar West" = "Mjini Magharibi", "Dar-Es-Salaam" = "Dar es Salaam")
# For this data, there was no mismatch. So we are good to go. Just run the empty manual_fix <- c()

manual_fix <- c()

tza_regions <- tza_regions %>%
  mutate(NAME_1_fixed = if (length(manual_fix) > 0) dplyr::recode(NAME_1, !!!manual_fix) else NAME_1,
         region_name_norm = norm_name(NAME_1_fixed))

region_zone_lookup <- region_zone_lookup %>%
  mutate(region_name_norm = norm_name(region_name))

# ---- Join zone labels onto the region polygons, then DISSOLVE to 9 zones ------
tza_zones <- tza_regions %>%
  left_join(region_zone_lookup %>% select(region_name_norm, zone),
            by = "region_name_norm")

if (any(is.na(tza_zones$zone))) {
  warning("Some GADM regions did not match a DHS zone -- inspect & fix via `manual_fix` above.")
  print(tza_zones %>% sf::st_drop_geometry() %>% filter(is.na(zone)) %>% distinct(NAME_1))
}

tza_zones <- tza_zones %>%
  filter(!is.na(zone)) %>%
  group_by(zone) %>%
  summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
  sf::st_as_sf()

# ---- Attach the zone-level prevalence already computed in Section 7b ----------
zone_prev <- weighted_by_zone %>%
  dplyr::select(zonal_region, malaria_prev, anaemia_prev, double_burden_prev)

zone_map_data <- tza_zones %>%
  left_join(zone_prev, by = c("zone" = "zonal_region"))

# ---- Re-usable choropleth plotting function ------------------------------------
plot_zone_choropleth <- function(data, fill_var, fill_label, palette_option) {
  ggplot(data) +
    geom_sf(aes(fill = .data[[fill_var]]), colour = "white", linewidth = 0.3) +
    geom_sf_text(aes(label = paste0(zone, "\n", round(.data[[fill_var]], 1), "%")),
                 size = 4, colour = "white", fontface = "bold") +
    scale_fill_viridis_c(option = palette_option, name = fill_label) +
    labs(subtitle = "Tanzania DHS-MIS 2022, children 6\u201359 months",
         caption = "Source: 2022 TDHS-MIS; boundaries: GADM v4.1 (TZA, level 1), dissolved to DHS zones") +
    theme_minimal(base_size = 14) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
}

choropleth_malaria <- plot_zone_choropleth(zone_map_data, "malaria_prev",
                        "Malaria\nprevalence (%)", "C") +
  labs(title = "Malaria prevalence by DHS zone")
print(choropleth_malaria)
ggsave("choropleth_malaria_by_zone.png", choropleth_malaria, width = 8, height = 8, dpi = 300)

choropleth_anaemia <- plot_zone_choropleth(zone_map_data, "anaemia_prev",
                        "Anaemia\nprevalence (%)", "D") +
  labs(title = "Anaemia prevalence by DHS zone")
print(choropleth_anaemia)
ggsave("choropleth_anaemia_by_zone.png", choropleth_anaemia, width = 8, height = 8, dpi = 300)

choropleth_double <- plot_zone_choropleth(zone_map_data, "double_burden_prev",
                        "Double burden\nprevalence (%)", "B") +
  labs(title = "Malaria-related anaemia (double burden) by DHS zone")
print(choropleth_double)
ggsave("choropleth_doubleburden_by_zone.png", choropleth_double, width = 8, height = 8, dpi = 300)

# ---- 8h. Cluster POINTS overlaid on the zone boundaries ------------------------
# Same cluster-level bubbles as Section 8c-8e, but now drawn on top of the
# dissolved DHS zone outlines (`tza_zones`, built in 8g) for geographic context.
plot_points_over_zones <- function(zones, points, fill_var, fill_label, palette_option) {
  ggplot() +
    geom_sf(data = zones, fill = "grey96", colour = "grey45", linewidth = 0.4) +
    geom_sf_text(data = zones, aes(label = zone), size = 4.2,
                 colour = "grey35", fontface = "italic") +
    geom_sf(data = points, aes(size = n_children, colour = .data[[fill_var]]), alpha = 0.85) +
    scale_colour_viridis_c(option = palette_option, name = fill_label) +
    scale_size_continuous(name = "Children\ntested", range = c(1, 6)) +
    labs(subtitle = "Cluster-level points within DHS zone boundaries, Tanzania DHS-MIS 2022",
         caption = "Source: 2022 TDHS-MIS (PR & GE recodes); zone boundaries: GADM v4.1, dissolved") +
    theme_minimal(base_size = 14) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), legend.position = "right")
}

points_malaria <- plot_points_over_zones(tza_zones, map_data, "malaria_prev",
                     "Malaria\nprevalence (%)", "C") +
  labs(title = "Malaria prevalence: cluster points within DHS zones")
print(points_malaria)
ggsave("points_over_zones_malaria.png", points_malaria, width = 9, height = 9, dpi = 300)

points_anaemia <- plot_points_over_zones(tza_zones, map_data, "anaemia_prev",
                     "Anaemia\nprevalence (%)", "D") +
  labs(title = "Anaemia prevalence: cluster points within DHS zones")
print(points_anaemia)
ggsave("points_over_zones_anaemia.png", points_anaemia, width = 9, height = 9, dpi = 300)

points_double <- plot_points_over_zones(tza_zones, map_data, "double_burden_prev",
                    "Double burden\nprevalence (%)", "B") +
  labs(title = "Double burden: cluster points within DHS zones")
print(points_double)
ggsave("points_over_zones_doubleburden.png", points_double, width = 9, height = 9, dpi = 300)

# ---- 8i. REGION-level choropleths (31 individual regions, not dissolved to zones) ----
# Uses the same `tza_regions` polygons and `region_zone_lookup` crosswalk from
# Section 8g, but skips the dissolve step so each of the 26 mainland + 5
# Zanzibar regions is shown individually -- a finer-grained view than the
# 9-zone maps above.
region_prev <- data_clean %>%
  left_join(region_zone_lookup %>% select(hv024, region_name) %>% distinct(), by = "hv024") %>%
  mutate(wt = v005 / 1e6) %>%
  as_survey_design(ids = v021, strata = v022, weights = wt, nest = TRUE) %>%
  group_by(region_name) %>%
  summarise(
    n = unweighted(n()),
    malaria_prev       = survey_mean(any_malaria)   * 100,
    anaemia_prev       = survey_mean(any_anaemia)   * 100,
    double_burden_prev = survey_mean(double_burden) * 100) %>%
  select(region_name, n, malaria_prev, anaemia_prev, double_burden_prev) %>%
  mutate(region_name_norm = norm_name(region_name))

region_map_data <- tza_regions %>%
  left_join(region_prev %>% select(-region_name), by = "region_name_norm")

if (any(is.na(region_map_data$n))) {
  warning("Some regions have no matched prevalence data -- check region_prev / region_name_norm.")
}

plot_region_choropleth <- function(data, fill_var, fill_label, palette_option) {
  ggplot(data) +
    geom_sf(aes(fill = .data[[fill_var]]), colour = "white", linewidth = 0.15) +
    scale_fill_viridis_c(option = palette_option, name = fill_label, na.value = "grey85") +
    labs(subtitle = "Tanzania DHS-MIS 2022, children 6\u201359 months \u2014 region level",
         caption = "Source: 2022 TDHS-MIS; boundaries: GADM v4.1 (TZA, level 1)") +
    theme_minimal(base_size = 14) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
}
# NB: region labels are omitted here (31 polygons, several very small Zanzibar
# ones) to avoid clutter -- pair this map with `print(region_prev)` for the
# exact per-region numbers, or add geom_sf_text(aes(label = region_name),
# size = 1.8) if you'd like labels anyway.

region_choro_malaria <- plot_region_choropleth(region_map_data, "malaria_prev",
                           "Malaria\nprevalence (%)", "C") +
  labs(title = "Malaria prevalence by region")
print(region_choro_malaria)
ggsave("choropleth_malaria_by_region.png", region_choro_malaria, width = 8, height = 8, dpi = 300)

region_choro_anaemia <- plot_region_choropleth(region_map_data, "anaemia_prev",
                           "Anaemia\nprevalence (%)", "D") +
  labs(title = "Anaemia prevalence by region")
print(region_choro_anaemia)
ggsave("choropleth_anaemia_by_region.png", region_choro_anaemia, width = 8, height = 8, dpi = 300)

region_choro_double <- plot_region_choropleth(region_map_data, "double_burden_prev",
                           "Double burden\nprevalence (%)", "B") +
  labs(title = "Double-burden prevalence by region")
print(region_choro_double)
ggsave("choropleth_doubleburden_by_region.png", region_choro_double, width = 8, height = 8, dpi = 300)


library(spdep)
nb <- knn2nb(knearneigh(sf::st_coordinates(map_data), k = 8, longlat = TRUE))
lw <- nb2listw(nb, style = "W")
for (v in c("malaria_prev", "anaemia_prev", "double_burden_prev")) { cat("\n", v, "\n"); print(moran.test(map_data[[v]], lw)) }

# ---- 8j. FACET-WRAP maps: Region | DHS zone | Survey cluster in ONE figure ----
lev3 <- c("Region", "DHS zone", "Survey cluster")


facet_map <- function(var, title_txt, hue) {
  poly_reg <- region_map_data %>%
    mutate(level = factor("Region", levels = lev3), value = .data[[var]]) %>%
    select(level, value)
  poly_zone <- zone_map_data %>%
    mutate(level = factor("DHS zone", levels = lev3), value = .data[[var]]) %>%
    select(level, value)
  zone_lab <- poly_zone
  sf::st_geometry(zone_lab) <- suppressWarnings(
    sf::st_point_on_surface(sf::st_geometry(poly_zone)))
  zone_out <- tza_zones %>%
    mutate(level = factor("Survey cluster", levels = lev3)) %>% select(level)
  pts <- map_data %>%
    sf::st_transform(sf::st_crs(tza_zones)) %>%
    mutate(level = factor("Survey cluster", levels = lev3), value = .data[[var]]) %>%
    select(level, value, n_children)
  cap <- as.numeric(quantile(pts$value, 0.95, na.rm = TRUE))  # colour scale cap
  
  ggplot() +
    geom_sf(data = poly_reg,  aes(fill = value), colour = "white", linewidth = 0.1) +
    geom_sf(data = poly_zone, aes(fill = value), colour = "white", linewidth = 0.3) +
    geom_sf(data = zone_out, fill = "grey95", colour = "grey55", linewidth = 0.3) +
    geom_sf(data = pts, aes(fill = value, size = n_children),
            shape = 21, colour = "grey25", stroke = 0.15, alpha = 0.9) +
    geom_sf_label(data = zone_lab,
                  aes(label = format(round(value, 1), nsmall = 1)),
                  size = 3.4, fill = scales::alpha("white", 0.75)) +
    scale_fill_gradient(low = "grey95", high = hue, limits = c(0, cap),
                        oob = scales::squish, na.value = "grey85",
                        name = "Prevalence (%)") +
    scale_size_continuous(range = c(1, 3), name = "Children\ntested") +
    facet_wrap(~ level, nrow = 1) +
    labs(title = title_txt,
         caption = paste0("Colour scale capped at the 95th percentile of cluster values. ",
                          "Region/zone: survey-weighted; cluster: unweighted. ",
                          "Source: 2022 TDHS-MIS; GADM v4.1.")) +
    theme_minimal(base_size = 13) +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 14),
          plot.title = element_text(face = "bold"),
          plot.caption = element_text(size = 9, hjust = 0),
          legend.position = "bottom", legend.key.width = unit(1.6, "cm"))
}

fm_anaemia <- facet_map("anaemia_prev", "Anaemia prevalence (any anaemia, Groups A + B)", "#B2182B")
fm_malaria <- facet_map("malaria_prev", "Malaria prevalence (RDT-positive, Groups A + C)", "#2166AC")
fm_double  <- facet_map("double_burden_prev", "Double burden (Group A: malaria-related anaemia)", "#762A83")
ggsave("Fig2_anaemia_region_zone_cluster.png", fm_anaemia, width = 9, height = 4.4, dpi = 300)
ggsave("Fig3_malaria_region_zone_cluster.png", fm_malaria, width = 9, height = 4.4, dpi = 300)
ggsave("Fig4_doubleburden_region_zone_cluster.png", fm_double, width = 9, height = 4.4, dpi = 300)


#===============================================================================
# 9. COMPLEX SURVEY DESIGN (already created above as `full_design`; repeated
#    here for clarity / in case Section 7b is skipped)
#===============================================================================
full_design <- data_clean %>%
  mutate(wt = v005 / 1e6) %>%
  as_survey_design(ids = v021, strata = v022, weights = wt, nest = TRUE)

#===============================================================================
# 10. RAO-SCOTT CHI-SQUARE TESTS
#     (bivariate association between each covariate and each outcome,
#      accounting for the complex survey design)
#===============================================================================
covariates <- c("age_child","gender_child","residence","zonal_region","birth_order",
                 "stunting","wasting","underwt","fever","diarrhea","cough","deworming",
                 "vitA","iron_supp","insurance","age_mom","marital_stats","education",
                 "anemia_mother","wealth_index","hh_under5","hh_size","head_sex",
                 "head_age","water_source","floor_material","wall_material",
                 "roof_material","hh_owns_itn")

run_rao_scott <- function(design, outcome, covariates) {
  purrr::map_dfr(covariates, function(v) {
    f <- as.formula(paste0("~", v, "+", outcome))
    test <- tryCatch(svychisq(f, design = design, statistic = "F"),
                      error = function(e) NULL)
    if (is.null(test)) return(tibble(variable = v, statistic = NA, df = NA, p_value = NA))
    tibble(variable = v,
           statistic = unname(test$statistic),
           df = paste(test$parameter, collapse = ", "),
           p_value = test$p.value)
  })
}

raoscott_malaria <- run_rao_scott(full_design, "malaria_child", covariates)
raoscott_anaemia <- run_rao_scott(full_design, "anaemia_child", covariates)
raoscott_group   <- run_rao_scott(full_design, "mal_anaemia_group", covariates)

print(raoscott_malaria)
print(raoscott_anaemia)
print(raoscott_group)

#===============================================================================
# 11. BIVARIATE MODIFIED POISSON REGRESSION (unadjusted prevalence ratios)
#     Modified Poisson = quasipoisson GLM w/ robust (sandwich) SEs, which the
#     `survey` package already supplies via svyglm on a survey design object.
#===============================================================================

# Two analytic subsets, created by SUBSETTING the full design (recommended
# over pre-filtering the data.frame, so PSU/strata info is preserved correctly)
design_A <- full_design %>%                      # Model A: malaria-related anaemia
  filter(mal_anaemia_group %in% c("A_malaria_anaemia", "D_neither")) %>%
  mutate(outcome_A = if_else(mal_anaemia_group == "A_malaria_anaemia", 1, 0))

design_B <- full_design %>%                      # Model B: non-malaria anaemia
  filter(mal_anaemia_group %in% c("B_nonmalaria_anaemia", "D_neither")) %>%
  mutate(outcome_B = if_else(mal_anaemia_group == "B_nonmalaria_anaemia", 1, 0))

bivariate_pr <- function(design, outcome, covariates) {
  purrr::map_dfr(covariates, function(v) {
    f <- as.formula(paste0(outcome, " ~ ", v))
    m <- tryCatch(svyglm(f, design = design, family = quasipoisson()),
                  error = function(e) NULL)
    if (is.null(m)) return(tibble(variable = v))
    broom::tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term != "(Intercept)") %>%
      mutate(variable = v) %>%
      select(variable, term, PR = estimate, conf.low, conf.high, p.value)
  })
}

pr_modelA <- bivariate_pr(design_A, "outcome_A", covariates)   # malaria-related anaemia
pr_modelB <- bivariate_pr(design_B, "outcome_B", covariates)   # non-malaria anaemia

print(pr_modelA, n = 100)
print(pr_modelB, n = 100)

#===============================================================================
# 12. MULTIVARIABLE LOGISTIC REGRESSION
#     Adjust for variables significant at p < 0.2 in the bivariate stage
#     (adjust the two vectors below once you have inspected pr_modelA/pr_modelB)
#I am still trying to build and fit this model. So it this part of the code is still a work in progress.
#===============================================================================

# --- edit these two lists based on your bivariate (Section 11) results ---------
final_vars_A <- c("age_child","residence","zonal_region","wealth_index",
                   "education","anemia_mother","hh_owns_itn","water_source",
                   "hh_size","stunting","wasting")

final_vars_B <- c("age_child","residence","zonal_region","wealth_index",
                   "education","anemia_mother","water_source","hh_size",
                   "stunting","wasting","fever","deworming")

formula_A <- as.formula(paste("outcome_A ~", paste(final_vars_A, collapse = " + ")))
formula_B <- as.formula(paste("outcome_B ~", paste(final_vars_B, collapse = " + ")))

model_A <- svyglm(formula_A, design = design_A, family = quasibinomial())
model_B <- svyglm(formula_B, design = design_B, family = quasibinomial())

summary(model_A)
summary(model_B)

tidy_A <- broom::tidy(model_A, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(model = "A: Malaria-related anaemia")

tidy_B <- broom::tidy(model_B, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(model = "B: Non-malaria-related anaemia")

results_table <- bind_rows(tidy_A, tidy_B) %>%
  rename(AOR = estimate) %>%
  select(model, term, AOR, conf.low, conf.high, p.value) %>%
  mutate(across(c(AOR, conf.low, conf.high), ~ round(.x, 2)),
         p.value = round(p.value, 3))
print(results_table, n = 100)
# write_csv(results_table, "AOR_results_ModelA_ModelB.csv")

#===============================================================================
# 13. OVERLAP ANALYSIS
#     Classify each covariate's TERMS as shared, malaria-anaemia-specific,
#     non-malaria-anaemia-specific, or divergent (per progress-report Section 4.4)
#Work in progress
#===============================================================================
sig_A <- tidy_A %>% filter(p.value < 0.05) %>% pull(term)
sig_B <- tidy_B %>% filter(p.value < 0.05) %>% pull(term)

overlap_table <- tibble(term = union(sig_A, sig_B)) %>%
  mutate(
    significant_in_A = term %in% sig_A,
    significant_in_B = term %in% sig_B,
    category = case_when(
      significant_in_A & significant_in_B ~ "Shared risk factor",
      significant_in_A & !significant_in_B ~ "Malaria-anaemia specific",
      !significant_in_A & significant_in_B ~ "Non-malaria-anaemia specific",
      TRUE ~ "Other")
  ) %>%
  left_join(tidy_A %>% select(term, AOR_A = estimate), by = "term") %>%
  left_join(tidy_B %>% select(term, AOR_B = estimate), by = "term")

print(overlap_table, n = 100)
# write_csv(overlap_table, "overlap_determinants_table.csv")

cat("\n--- Summary of overlap analysis ---\n")
print(table(overlap_table$category))

#===============================================================================
# END OF SCRIPT
# Next: paste the weighted prevalence tables, the three maps, results_table and
# overlap_table into the Results section of the progress/final report.
#===============================================================================





##### New model nfittings
# ============================================================
# MODEL A: MALARIA-RELATED ANAEMIA
# Excluding zonal_region because of separation and instability
# ============================================================

library(survey)
library(dplyr)
library(broom)
library(tibble)

# Fit the survey-weighted multivariable logistic regression
model_A_final <- survey::svyglm(
  outcome_A ~ age_child +
    residence +
    wealth_index +
    education +
    anemia_mother +
    hh_owns_itn +
    water_source +
    hh_size +
    stunting +
    wasting,
  design = design_A,
  family = quasibinomial()
)

# Display model summary
summary(model_A_final)

# Check analytic sample size
nobs(model_A_final)



# Extract estimates and survey-design-based confidence intervals
coef_A <- coef(model_A_final)
ci_A <- confint(model_A_final)
coef_summary_A <- summary(model_A_final)$coefficients

results_A <- tibble(
  term = names(coef_A),
  AOR = exp(coef_A),
  Lower_95_CI = exp(ci_A[, 1]),
  Upper_95_CI = exp(ci_A[, 2]),
  p_value = coef_summary_A[, 4]
) %>%
  # Remove the intercept from the publication table
  filter(term != "(Intercept)") %>%
  mutate(
    AOR = round(AOR, 2),
    Lower_95_CI = round(Lower_95_CI, 2),
    Upper_95_CI = round(Upper_95_CI, 2),
    p_value = signif(p_value, 3)
  )

# Display the results
print(as.data.frame(results_A), row.names = FALSE)

results_table_A <- results_A %>%
  mutate(
    `AOR (95% CI)` = sprintf(
      "%.2f (%.2f–%.2f)",
      AOR,
      Lower_95_CI,
      Upper_95_CI
    ),
    `p-value` = ifelse(
      p_value < 0.001,
      "<0.001",
      sprintf("%.3f", p_value)
    )
  ) %>%
  select(
    Predictor = term,
    `AOR (95% CI)`,
    `p-value`
  )

print(as.data.frame(results_table_A), row.names = FALSE)

# Install once if necessary:
# install.packages("openxlsx")

openxlsx::write.xlsx(
  results_table_A,
  file = "Model_A_Final_Adjusted_Results.xlsx",
  overwrite = TRUE
)





library(survey)
library(dplyr)
library(tibble)

# ------------------------------------------------------------
# MODEL B: NON-MALARIA-RELATED ANAEMIA
# Initial model retaining zonal_region
# ------------------------------------------------------------

model_B_region <- survey::svyglm(
  outcome_B ~ age_child +
    residence +
    zonal_region +
    wealth_index +
    education +
    anemia_mother +
    water_source +
    hh_size +
    stunting +
    wasting +
    fever +
    deworming,
  design = design_B,
  family = quasibinomial()
)

summary(model_B_region)

# Check regional separation using ordinary logistic regression.
# This is a diagnostic, not the final survey-weighted analysis.
sep_B <- glm(
  outcome_B ~ zonal_region +
    age_child +
    residence +
    wealth_index +
    education +
    anemia_mother +
    water_source +
    hh_size +
    stunting +
    wasting +
    fever +
    deworming,
  data = design_B$variables,
  family = binomial(),
  method = detectseparation::detect_separation
)

print(sep_B)



# ------------------------------------------------------------
# MODEL B: FINAL CANDIDATE WITHOUT REGIONAL FIXED EFFECTS
# ------------------------------------------------------------

model_B_final <- survey::svyglm(
  outcome_B ~ age_child +
    residence +
    wealth_index +
    education +
    anemia_mother +
    water_source +
    hh_size +
    stunting +
    wasting +
    fever +
    deworming,
  design = design_B,
  family = quasibinomial()
)

summary(model_B_final)

# Analytic sample size
nobs(model_B_final)



# ------------------------------------------------------------
# EXTRACT MODEL B RESULTS
# ------------------------------------------------------------

coef_B <- coef(model_B_final)
ci_B <- confint(model_B_final)
summary_B <- summary(model_B_final)$coefficients

results_B <- tibble(
  term = names(coef_B),
  AOR = exp(coef_B),
  Lower_95_CI = exp(ci_B[, 1]),
  Upper_95_CI = exp(ci_B[, 2]),
  p_value = summary_B[, 4]
) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    AOR = round(AOR, 2),
    Lower_95_CI = round(Lower_95_CI, 2),
    Upper_95_CI = round(Upper_95_CI, 2),
    p_value = signif(p_value, 3)
  )

print(as.data.frame(results_B), row.names = FALSE)




results_table_B <- results_B %>%
  mutate(
    `AOR (95% CI)` = sprintf(
      "%.2f (%.2f–%.2f)",
      AOR, Lower_95_CI, Upper_95_CI
    ),
    `p-value` = ifelse(
      p_value < 0.001,
      "<0.001",
      sprintf("%.3f", p_value)
    )
  ) %>%
  select(
    Predictor = term,
    `AOR (95% CI)`,
    `p-value`
  )

print(as.data.frame(results_table_B), row.names = FALSE)

# Export to Excel
openxlsx::write.xlsx(
  results_table_B,
  file = "Model_B_Final_Adjusted_Results.xlsx",
  overwrite = TRUE
)



# ------------------------------------------------------------
# HARMONISED MODELS FOR COMPARING DETERMINANTS
# Identical predictors and no regional fixed effects
# ------------------------------------------------------------

shared_formula <- outcome_A ~ age_child +
  residence +
  wealth_index +
  education +
  anemia_mother +
  water_source +
  hh_size +
  stunting +
  wasting

# Model A: malaria-related anaemia vs neither condition
model_A_overlap <- survey::svyglm(
  shared_formula,
  design = design_A,
  family = quasibinomial()
)

# Model B: non-malaria-related anaemia vs neither condition
# Update the outcome in the formula
shared_formula_B <- update(
  shared_formula,
  outcome_B ~ .
)

model_B_overlap <- survey::svyglm(
  shared_formula_B,
  design = design_B,
  family = quasibinomial()
)

summary(model_A_overlap)
summary(model_B_overlap)

# Compare sample sizes
nobs(model_A_overlap)
nobs(model_B_overlap)



# ------------------------------------------------------------
# FUNCTION TO EXTRACT RESULTS FROM SURVEY-WEIGHTED MODELS
# ------------------------------------------------------------

extract_svy_results <- function(model, model_name) {
  
  estimates <- coef(model)
  intervals <- confint(model)
  coefficient_table <- summary(model)$coefficients
  
  tibble(
    term = names(estimates),
    model = model_name,
    AOR = exp(estimates),
    Lower_95_CI = exp(intervals[, 1]),
    Upper_95_CI = exp(intervals[, 2]),
    p_value = coefficient_table[, 4]
  ) %>%
    filter(term != "(Intercept)")
}

overlap_A <- extract_svy_results(
  model_A_overlap,
  "Model A: Malaria-related anaemia"
)

overlap_B <- extract_svy_results(
  model_B_overlap,
  "Model B: Non-malaria-related anaemia"
)

# Combine the estimates into one comparison table
overlap_results <- bind_rows(overlap_A, overlap_B)

print(as.data.frame(overlap_results), row.names = FALSE)





####Fit the final candidate Model B with region retained
model_B_final <- survey::svyglm(
  outcome_B ~ age_child +
    residence +
    zonal_region +
    wealth_index +
    education +
    anemia_mother +
    water_source +
    hh_size +
    stunting +
    wasting +
    fever +
    deworming,
  design = design_B,
  family = quasibinomial()
)

summary(model_B_final)

# Check sample size
nobs(model_B_final)


coef_B <- coef(model_B_final)
ci_B <- confint(model_B_final)
summary_B <- summary(model_B_final)$coefficients

results_B <- tibble::tibble(
  term = names(coef_B),
  AOR = exp(coef_B),
  Lower_95_CI = exp(ci_B[, 1]),
  Upper_95_CI = exp(ci_B[, 2]),
  p_value = summary_B[, 4]
) %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(
    AOR = round(AOR, 2),
    Lower_95_CI = round(Lower_95_CI, 2),
    Upper_95_CI = round(Upper_95_CI, 2),
    p_value = signif(p_value, 3)
  )

print(as.data.frame(results_B), row.names = FALSE)


# Extract estimates from Model B
coef_B <- coef(model_B_final)
ci_B <- confint(model_B_final)
summary_B <- summary(model_B_final)$coefficients

results_B <- tibble::tibble(
  term = names(coef_B),
  AOR = exp(coef_B),
  Lower_95_CI = exp(ci_B[, 1]),
  Upper_95_CI = exp(ci_B[, 2]),
  p_value = summary_B[, 4]
) %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(
    AOR = round(AOR, 2),
    Lower_95_CI = round(Lower_95_CI, 2),
    Upper_95_CI = round(Upper_95_CI, 2),
    p_value = signif(p_value, 3)
  )

print(as.data.frame(results_B), row.names = FALSE)

# Publication-style table
results_table_B <- results_B %>%
  dplyr::mutate(
    `AOR (95% CI)` = sprintf(
      "%.2f (%.2f–%.2f)",
      AOR, Lower_95_CI, Upper_95_CI
    ),
    `p-value` = ifelse(
      p_value < 0.001,
      "<0.001",
      sprintf("%.3f", p_value)
    )
  ) %>%
  dplyr::select(
    Predictor = term,
    `AOR (95% CI)`,
    `p-value`
  )

print(as.data.frame(results_table_B), row.names = FALSE)

# Export
openxlsx::write.xlsx(
  results_table_B,
  "Model_B_Final_Adjusted_Results.xlsx",
  overwrite = TRUE
)


# Shared predictors for the overlap analysis
model_A_overlap <- survey::svyglm(
  outcome_A ~ age_child +
    residence +
    wealth_index +
    education +
    anemia_mother +
    water_source +
    hh_size +
    stunting +
    wasting,
  design = design_A,
  family = quasibinomial()
)

model_B_overlap <- survey::svyglm(
  outcome_B ~ age_child +
    residence +
    wealth_index +
    education +
    anemia_mother +
    water_source +
    hh_size +
    stunting +
    wasting,
  design = design_B,
  family = quasibinomial()
)

# Extract both models
extract_svy_results <- function(model, model_name) {
  b <- coef(model)
  ci <- confint(model)
  sm <- summary(model)$coefficients
  
  tibble::tibble(
    term = names(b),
    model = model_name,
    AOR = exp(b),
    Lower_95_CI = exp(ci[, 1]),
    Upper_95_CI = exp(ci[, 2]),
    p_value = sm[, 4]
  ) %>%
    dplyr::filter(term != "(Intercept)")
}

overlap_A <- extract_svy_results(
  model_A_overlap, "Malaria-related anaemia"
)

overlap_B <- extract_svy_results(
  model_B_overlap, "Non-malaria-related anaemia"
)

# Match equivalent predictors
overlap_comparison <- overlap_A %>%
  dplyr::select(
    term,
    AOR_A = AOR,
    Lower_A = Lower_95_CI,
    Upper_A = Upper_95_CI,
    p_A = p_value
  ) %>%
  dplyr::inner_join(
    overlap_B %>%
      dplyr::select(
        term,
        AOR_B = AOR,
        Lower_B = Lower_95_CI,
        Upper_B = Upper_95_CI,
        p_B = p_value
      ),
    by = "term"
  ) %>%
  dplyr::mutate(
    Direction_A = ifelse(AOR_A > 1, "Higher odds", "Lower odds"),
    Direction_B = ifelse(AOR_B > 1, "Higher odds", "Lower odds"),
    Same_direction = Direction_A == Direction_B,
    Evidence_A = Lower_A > 1 | Upper_A < 1,
    Evidence_B = Lower_B > 1 | Upper_B < 1,
    Pattern = dplyr::case_when(
      Same_direction & Evidence_A & Evidence_B ~
        "Same direction; CI excludes 1 in both",
      Same_direction & !Evidence_A & !Evidence_B ~
        "Same direction; CI includes 1 in both",
      Same_direction ~
        "Same direction; strength of evidence differs",
      TRUE ~ "Direction differs"
    )
  ) %>%
  dplyr::mutate(
    dplyr::across(
      c(AOR_A, Lower_A, Upper_A, AOR_B, Lower_B, Upper_B),
      ~ round(.x, 2)
    ),
    dplyr::across(c(p_A, p_B), ~ signif(.x, 3))
  )

print(as.data.frame(overlap_comparison), row.names = FALSE)

openxlsx::write.xlsx(
  overlap_comparison,
  "Model_A_Model_B_Overlap_Analysis.xlsx",
  overwrite = TRUE
)






######## Joint modelling to establish statistical difference in factors driving the two forms of anaemia

# Install once if required
# install.packages(c("svyVGAM", "VGAM"))

library(survey)
library(svyVGAM)
library(VGAM)
library(dplyr)
library(tibble)


# Keep only the three outcome groups
design_joint <- subset(
  full_design,
  mal_anaemia_group %in% c(
    "D_neither",
    "A_malaria_anaemia",
    "B_nonmalaria_anaemia"
  )
)

# Set the outcome reference explicitly:
# D_neither is the reference outcome
design_joint <- update(
  design_joint,
  outcome_joint = factor(
    as.character(mal_anaemia_group),
    levels = c(
      "D_neither",
      "A_malaria_anaemia",
      "B_nonmalaria_anaemia"
    )
  )
)

# Use the shared predictors from the harmonised models.
# Region is omitted because of Model A's separation problem.
joint_formula <- outcome_joint ~
  age_child +
  residence +
  wealth_index +
  education +
  anemia_mother +
  water_source +
  hh_size +
  stunting +
  wasting

# Use a common complete-case sample for all predictors
required_vars <- c(
  "outcome_joint", "age_child", "residence",
  "wealth_index", "education", "anemia_mother",
  "water_source", "hh_size", "stunting", "wasting"
)

design_joint <- subset(
  design_joint,
  complete.cases(
    design_joint$variables[, required_vars, drop = FALSE]
  )
)

# Fit the survey-weighted multinomial logistic regression
model_joint <- svyVGAM::svy_vglm(
  joint_formula,
  design = design_joint,
  family = VGAM::multinomial(refLevel = 1)
)

summary(model_joint)

# Check the common analytic sample
nrow(design_joint$variables)

# Inspect coefficient names before constructing contrasts
names(coef(model_joint))




# ------------------------------------------------------------
# INDIVIDUAL COEFFICIENT CONTRASTS
# Tests whether A and B coefficients differ
# ------------------------------------------------------------

b <- coef(model_joint)
V <- vcov(model_joint)
coef_names <- names(b)

# Confirm expected coefficient naming convention
print(coef_names)

# Identify predictors represented in both outcome equations
terms_A <- sub(":1$", "", coef_names[grepl(":1$", coef_names)])
terms_B <- sub(":2$", "", coef_names[grepl(":2$", coef_names)])

common_terms <- intersect(terms_A, terms_B)
common_terms <- setdiff(common_terms, "(Intercept)")

# Survey design degrees of freedom
df_design <- survey::degf(design_joint)

contrast_results <- lapply(common_terms, function(term) {
  
  name_A <- paste0(term, ":1")
  name_B <- paste0(term, ":2")
  
  # Difference between the two log-odds coefficients
  L <- setNames(rep(0, length(b)), coef_names)
  L[name_A] <- 1
  L[name_B] <- -1
  
  difference <- as.numeric(sum(L * b))
  variance_difference <- as.numeric(t(L) %*% V %*% L)
  se_difference <- sqrt(variance_difference)
  
  t_stat <- difference / se_difference
  p_value <- 2 * stats::pt(
    abs(t_stat),
    df = df_design,
    lower.tail = FALSE
  )
  
  critical_t <- stats::qt(0.975, df = df_design)
  
  tibble::tibble(
    term = term,
    log_OR_difference = difference,
    ratio_of_AORs = exp(difference),
    Lower_95_CI = exp(
      difference - critical_t * se_difference
    ),
    Upper_95_CI = exp(
      difference + critical_t * se_difference
    ),
    t_statistic = t_stat,
    df = df_design,
    p_value = p_value
  )
})

contrast_results <- dplyr::bind_rows(contrast_results) %>%
  dplyr::mutate(
    dplyr::across(
      c(
        log_OR_difference, ratio_of_AORs,
        Lower_95_CI, Upper_95_CI, t_statistic
      ),
      ~ round(.x, 3)
    ),
    p_value = signif(p_value, 3)
  )

print(as.data.frame(contrast_results), row.names = FALSE)



# ------------------------------------------------------------
# OMNIBUS TESTS: DO ASSOCIATIONS DIFFER BY OUTCOME?
# ------------------------------------------------------------

predictors_to_test <- c(
  "age_child",
  "anemia_mother",
  "residence",
  "water_source"
)

omnibus_results <- lapply(
  predictors_to_test,
  function(variable) {
    
    # Find coefficient terms belonging to this predictor
    predictor_terms <- common_terms[
      startsWith(common_terms, variable)
    ]
    
    if (length(predictor_terms) == 0) {
      return(NULL)
    }
    
    # Build one contrast row per category
    L <- matrix(
      0,
      nrow = length(predictor_terms),
      ncol = length(b),
      dimnames = list(predictor_terms, coef_names)
    )
    
    for (i in seq_along(predictor_terms)) {
      term <- predictor_terms[i]
      L[i, paste0(term, ":1")] <- 1
      L[i, paste0(term, ":2")] <- -1
    }
    
    difference_vector <- as.vector(L %*% b)
    contrast_vcov <- L %*% V %*% t(L)
    
    q <- nrow(L)
    
    wald_chisq <- as.numeric(
      t(difference_vector) %*%
        solve(contrast_vcov, difference_vector)
    )
    
    # Design-based approximate Wald F test
    F_stat <- wald_chisq / q
    p_value <- stats::pf(
      F_stat,
      df1 = q,
      df2 = df_design,
      lower.tail = FALSE
    )
    
    tibble::tibble(
      predictor = variable,
      contrasts_tested = q,
      Wald_F = F_stat,
      numerator_df = q,
      denominator_df = df_design,
      p_value = p_value
    )
  }
)

omnibus_results <- dplyr::bind_rows(omnibus_results) %>%
  dplyr::mutate(
    Wald_F = round(Wald_F, 3),
    p_value = signif(p_value, 3)
  )

print(as.data.frame(omnibus_results), row.names = FALSE)


openxlsx::write.xlsx(
  list(
    "Individual contrasts" = contrast_results,
    "Omnibus tests" = omnibus_results
  ),
  file = "Formal_Outcome_Overlap_Tests.xlsx",
  overwrite = TRUE
)




######## Some quality checks
# Packages
library(survey)
library(svyVGAM)
library(VGAM)

# 1. Inspect the model class and variance-covariance matrix
class(model_joint)

beta_joint <- coef(model_joint)
V_joint <- vcov(model_joint)

length(beta_joint)
dim(V_joint)

# Confirm matching coefficient names and matrix dimensions
stopifnot(
  nrow(V_joint) == length(beta_joint),
  ncol(V_joint) == length(beta_joint),
  all(names(beta_joint) %in% rownames(V_joint)),
  all(names(beta_joint) %in% colnames(V_joint))
)

# 2. Check the survey design degrees of freedom
survey::degf(design_joint)

# 3. Inspect the model's stored survey design, if available
names(model_joint)



###############################
library(survey)

beta <- coef(model_joint)
V <- vcov(model_joint)
df_design <- survey::degf(design_joint)

# Inspect the exact coefficient names first
print(names(beta))

# Helper: test whether coefficients differ between outcome equations
compare_outcome_coefficients <- function(term_names_A, term_names_B,
                                         beta, V, df_design) {
  
  stopifnot(
    length(term_names_A) == length(term_names_B),
    all(term_names_A %in% names(beta)),
    all(term_names_B %in% names(beta))
  )
  
  # Contrast matrix: beta_A - beta_B
  q <- length(term_names_A)
  L <- matrix(0, nrow = q, ncol = length(beta),
              dimnames = list(NULL, names(beta)))
  
  for (i in seq_len(q)) {
    L[i, term_names_A[i]] <-  1
    L[i, term_names_B[i]] <- -1
  }
  
  difference <- as.vector(L %*% beta)
  V_difference <- L %*% V %*% t(L)
  
  # Wald statistic
  W <- as.numeric(
    t(difference) %*%
      solve(V_difference, difference)
  )
  
  F_stat <- W / q
  p_value <- pf(F_stat, df1 = q, df2 = df_design,
                lower.tail = FALSE)
  
  data.frame(
    numerator_df = q,
    denominator_df = df_design,
    Wald_chisq = W,
    F_statistic = F_stat,
    p_value = p_value
  )
}



###############################################

compare_outcome_coefficients(
  term_names_A = "age_childYear-2:1",
  term_names_B = "age_childYear-2:2",
  beta = beta,
  V = V,
  df_design = df_design
)



# Compare coefficients between outcome equations:
# Equation 1 = malaria-related anaemia vs neither
# Equation 2 = non-malaria-related anaemia vs neither

compare_predictor <- function(prefix, beta, V, df_design) {
  
  coefficient_names <- names(beta)
  
  # Identify coefficients belonging to the predictor
  # and exclude the intercept
  names_1 <- coefficient_names[
    grepl(paste0("^", prefix), coefficient_names) &
      grepl(":1$", coefficient_names)
  ]
  
  names_2 <- sub(":1$", ":2", names_1)
  
  if (length(names_1) == 0) {
    stop("No matching coefficients found for: ", prefix)
  }
  
  if (!all(names_2 %in% coefficient_names)) {
    stop("Matching coefficients are missing for: ", prefix)
  }
  
  q <- length(names_1)
  
  # Construct the contrast matrix
  L <- matrix(
    0,
    nrow = q,
    ncol = length(beta),
    dimnames = list(NULL, coefficient_names)
  )
  
  for (i in seq_len(q)) {
    L[i, names_1[i]] <- 1
    L[i, names_2[i]] <- -1
  }
  
  # Differences between outcome-specific coefficients
  difference <- as.vector(L %*% beta)
  
  # Covariance matrix of the differences
  V_difference <- L %*% V %*% t(L)
  
  # Wald statistic
  W <- as.numeric(
    crossprod(difference, solve(V_difference, difference))
  )
  
  # Conventional design-df F approximation
  F_stat <- W / q
  p_value <- pf(
    F_stat,
    df1 = q,
    df2 = df_design,
    lower.tail = FALSE
  )
  
  data.frame(
    predictor = prefix,
    numerator_df = q,
    denominator_df = df_design,
    Wald_chisq = W,
    F_statistic = F_stat,
    p_value = p_value
  )
}

# Run the omnibus comparisons
results <- do.call(
  rbind,
  lapply(
    c("age_child", "residence", "water_source", "anemia_mother"),
    compare_predictor,
    beta = beta,
    V = V,
    df_design = df_design
  )
)

results$p_adjusted_Holm <- p.adjust(results$p_value, method = "holm")

print(results)


##############################################
# Compare the variance stored in the svyVGAM model
# with the covariance matrix returned by vcov()
dim(model_joint$var)
dim(model_joint$naive.var)

all.equal(
  unname(vcov(model_joint)),
  unname(model_joint$var)
)

# Inspect the fitted survey design stored in the model
class(model_joint$design)
survey::degf(model_joint$design)





# Compare design-based and naive covariance matrices
all.equal(
  unname(model_joint$var),
  unname(model_joint$naive.var)
)

# Compare their diagonal variances
variance_check <- data.frame(
  coefficient = names(coef(model_joint)),
  design_variance = diag(model_joint$var),
  naive_variance = diag(model_joint$naive.var)
)

head(variance_check, 10)


results$p_chisq <- pchisq(
  results$Wald_chisq,
  df = results$numerator_df,
  lower.tail = FALSE
)

results[, c(
  "predictor",
  "Wald_chisq",
  "numerator_df",
  "denominator_df",
  "p_value",
  "p_chisq"
)]
