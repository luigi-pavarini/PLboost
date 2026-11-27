# ==============================================================================
# build_df_lgd.r
# Builds data/df_lgd_all.RData (all states) from the Lending Club file
# accepted_2007_to_2018Q4.csv (Kaggle, not redistributed: see ../DATA.md).
# The analysis scripts then keep the Vermont loans (n = 320, p = 348).
#
# Input : $PLBOOST_OUT/data/accepted_2007_to_2018Q4.csv  (or pass its path)
# Output: $PLBOOST_OUT/data/df_lgd_all.RData
# Needs about 6 GB of RAM to read the 1.6 GB csv.
# ==============================================================================

library(data.table)
library(mltools)

# The one-hot column order depends on the collation; the master thesis data
# were built with en_US.UTF-8
if (!nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", "en_US.UTF-8"))) &&
  !nzchar(suppressWarnings(
    Sys.setlocale("LC_COLLATE", "English_United States.1252")
  ))) {
  warning(
    "could not set an en_US collation; ",
    "column order may differ from the master thesis"
  )
}

data_path <- file.path(Sys.getenv("PLBOOST_OUT", unset = "."), "data")
dir.create(data_path, showWarnings = FALSE, recursive = TRUE)
csv_in <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(csv_in)) csv_in <- file.path(data_path, "accepted_2007_to_2018Q4.csv")
if (!file.exists(csv_in)) {
  stop("Lending Club csv not found: ", csv_in, "  (see DATA.md)")
}

# ==============================================================================
# DEFAULTED LOANS AND SELECTED COLUMNS
# ==============================================================================

df <- fread(csv_in)
df_default <- df[loan_status %in% c("Charged Off", "Default")]

# recoveries is only known after default: used for the lgd, then dropped
cols_lgd_construction <- c("loan_amnt", "recoveries")

cols_feat_engineering <- c(
  "issue_d", "earliest_cr_line",
  "sec_app_earliest_cr_line", "last_pymnt_d"
)

cols_categorical <- c(
  "term", "grade", "sub_grade", "emp_length",
  "home_ownership", "verification_status", "purpose",
  "addr_state", "initial_list_status", "application_type",
  "disbursement_method", "pymnt_plan",
  "emp_title", "verification_status_joint"
)

cols_numeric_loan <- c(
  "loan_amnt", "funded_amnt", "funded_amnt_inv",
  "int_rate", "installment"
)

cols_numeric_borrower <- c(
  "annual_inc", "dti",
  "fico_range_low", "fico_range_high",
  "revol_bal", "revol_util",
  "delinq_2yrs", "inq_last_6mths",
  "open_acc", "pub_rec", "total_acc",
  "acc_now_delinq", "delinq_amnt",
  "mths_since_last_delinq", "mths_since_last_record",
  "collections_12_mths_ex_med",
  "mths_since_last_major_derog"
)

cols_numeric_bureau <- c(
  "tot_coll_amt", "tot_cur_bal",
  "total_rev_hi_lim", "avg_cur_bal",
  "bc_open_to_buy", "bc_util",
  "all_util", "il_util",
  "mort_acc", "pub_rec_bankruptcies", "tax_liens",
  "tot_hi_cred_lim", "total_bal_ex_mort",
  "total_bc_limit", "total_il_high_credit_limit",
  "acc_open_past_24mths",
  "pct_tl_nvr_dlq", "percent_bc_gt_75",
  "num_accts_ever_120_pd",
  "num_actv_bc_tl", "num_actv_rev_tl",
  "num_bc_sats", "num_bc_tl", "num_il_tl",
  "num_op_rev_tl", "num_rev_accts",
  "num_rev_tl_bal_gt_0", "num_sats",
  "num_tl_90g_dpd_24m", "num_tl_op_past_12m",
  "mo_sin_old_il_acct", "mo_sin_old_rev_tl_op",
  "mo_sin_rcnt_rev_tl_op", "mo_sin_rcnt_tl",
  "mths_since_recent_bc", "mths_since_recent_inq",
  "inq_fi", "inq_last_12m", "total_cu_tl",
  "max_bal_bc",
  "open_acc_6m", "open_act_il",
  "open_il_12m", "open_il_24m",
  "open_rv_12m", "open_rv_24m",
  "total_bal_il", "chargeoff_within_12_mths"
)

cols_numeric_joint <- c(
  "annual_inc_joint", "dti_joint", "revol_bal_joint",
  "sec_app_fico_range_low", "sec_app_fico_range_high",
  "sec_app_inq_last_6mths", "sec_app_mort_acc",
  "sec_app_open_acc", "sec_app_revol_util",
  "sec_app_open_act_il", "sec_app_num_rev_accts",
  "sec_app_chargeoff_within_12_mths",
  "sec_app_collections_12_mths_ex_med"
)

cols_numeric_other <- c("last_pymnt_amnt")

cols_select <- unique(c(
  cols_lgd_construction,
  cols_feat_engineering,
  cols_categorical,
  cols_numeric_loan,
  cols_numeric_borrower,
  cols_numeric_bureau,
  cols_numeric_joint,
  cols_numeric_other
))
cols_select <- cols_select[cols_select %in% names(df_default)]
df_selected <- df_default[, ..cols_select]

# ==============================================================================
# MISSING VALUES
# ==============================================================================

# Numeric columns with more than 10% NA are dropped, then rows with any NA in
# the remaining numeric columns. Categorical columns keep their NA and "".
num_cols_selected <- names(df_selected)[sapply(df_selected, is.numeric)]
na_frac <- colMeans(is.na(df_selected[, ..num_cols_selected]))
high_na_cols <- names(na_frac[na_frac > 0.10])
df_selected[, (high_na_cols) := NULL]

num_cols_clean <- names(df_selected)[sapply(df_selected, is.numeric)]
df_selected <- df_selected[complete.cases(df_selected[, ..num_cols_clean])]

cat("Rows with annual_inc = 0:", sum(df_selected$annual_inc == 0), "\n")
df_selected <- df_selected[annual_inc > 0]

# ==============================================================================
# LGD
# ==============================================================================

df_selected[, lgd := (loan_amnt - recoveries) / loan_amnt]
summary(df_selected$lgd)

df_lgd <- df_selected[lgd > 0 & lgd <= 1]
summary(df_lgd$lgd)
df_lgd[, recoveries := NULL]

# ==============================================================================
# FEATURE ENGINEERING
# ==============================================================================

df_lgd[, credit_age := as.numeric(difftime(
  as.Date(paste0("01-", issue_d), format = "%d-%b-%Y"),
  as.Date(paste0("01-", earliest_cr_line), format = "%d-%b-%Y"),
  units = "days"
)) / 365.25]

# 0 means no co-borrower
df_lgd[, sec_app_earliest_cr_line := as.character(sec_app_earliest_cr_line)]
df_lgd[, credit_age_coborrower := fifelse(
  sec_app_earliest_cr_line == "", 0,
  as.numeric(difftime(
    as.Date(paste0("01-", issue_d), format = "%d-%b-%Y"),
    as.Date(paste0("01-", sec_app_earliest_cr_line), format = "%d-%b-%Y"),
    units = "days"
  )) / 365.25
)]

df_lgd[, days_to_last_pymnt := as.numeric(difftime(
  as.Date(paste0("01-", last_pymnt_d), format = "%d-%b-%Y"),
  as.Date(paste0("01-", issue_d), format = "%d-%b-%Y"),
  units = "days"
))]

df_lgd[, c(
  "issue_d", "earliest_cr_line",
  "sec_app_earliest_cr_line", "last_pymnt_d"
) := NULL]

cat(
  "Rows dropped due to missing days_to_last_pymnt:",
  sum(is.na(df_lgd$days_to_last_pymnt)),
  "\n"
)
df_lgd <- df_lgd[!is.na(days_to_last_pymnt)]

# numeric columns that are zero in more than 95% of the rows
num_cols_check <- names(df_lgd)[sapply(df_lgd, is.numeric)]
num_cols_check <- num_cols_check[num_cols_check != "lgd"]
zero_frac <- sapply(
  num_cols_check,
  function(col) mean(df_lgd[[col]] == 0, na.rm = TRUE)
)
near_zero <- names(zero_frac[zero_frac > 0.95])
print(near_zero)
df_lgd[, (near_zero) := NULL]

df_lgd[, inter_loan_intrate := loan_amnt * int_rate]
df_lgd[, inter_loan_income := loan_amnt / annual_inc]
df_lgd[, inter_install_income := installment / annual_inc]
df_lgd[, inter_dti_intrate := dti * int_rate]
df_lgd[, inter_fico_intrate := fico_range_low * int_rate]
df_lgd[, inter_revol_income := revol_bal / annual_inc]
df_lgd[, inter_creditage_pubrec := credit_age * pub_rec]
df_lgd[, inter_inq_dti := inq_last_6mths * dti]
df_lgd[, inter_days_loan := days_to_last_pymnt * loan_amnt]

# ==============================================================================
# ONE-HOT ENCODING
# ==============================================================================

df_lgd[, term := as.factor(term)]
df_lgd <- one_hot(df_lgd, cols = "term", dropUnusedLevels = TRUE)

df_lgd[, grade := as.factor(grade)]
df_lgd <- one_hot(df_lgd, cols = "grade", dropUnusedLevels = TRUE)

df_lgd[, sub_grade := as.factor(sub_grade)]
df_lgd <- one_hot(df_lgd, cols = "sub_grade", dropUnusedLevels = TRUE)

df_lgd[, emp_length := as.factor(emp_length)]
df_lgd <- one_hot(df_lgd, cols = "emp_length", dropUnusedLevels = TRUE)

df_lgd[, home_ownership := as.factor(home_ownership)]
df_lgd <- one_hot(df_lgd, cols = "home_ownership", dropUnusedLevels = TRUE)

df_lgd[, verification_status := as.factor(verification_status)]
df_lgd <- one_hot(df_lgd, cols = "verification_status", dropUnusedLevels = TRUE)

df_lgd[, purpose := as.factor(purpose)]
df_lgd <- one_hot(df_lgd, cols = "purpose", dropUnusedLevels = TRUE)

state_counts <- sort(table(df_lgd$addr_state), decreasing = TRUE)
df_lgd[, addr_state := as.factor(addr_state)]
df_lgd <- one_hot(df_lgd, cols = "addr_state", dropUnusedLevels = TRUE)

df_lgd[, initial_list_status := as.factor(initial_list_status)]
df_lgd <- one_hot(df_lgd, cols = "initial_list_status", dropUnusedLevels = TRUE)

df_lgd[, application_type := as.factor(application_type)]
df_lgd <- one_hot(df_lgd, cols = "application_type", dropUnusedLevels = TRUE)

df_lgd[, disbursement_method := as.factor(disbursement_method)]
df_lgd <- one_hot(df_lgd, cols = "disbursement_method", dropUnusedLevels = TRUE)

df_lgd[, pymnt_plan := as.factor(pymnt_plan)]
df_lgd <- one_hot(df_lgd, cols = "pymnt_plan", dropUnusedLevels = TRUE)

# top 200 job titles; the rest, NA and "" go to the reference "Other"
df_lgd[, emp_title := as.character(emp_title)]
df_lgd[is.na(emp_title) | emp_title == "", emp_title := "Other"]
top_titles <- names(sort(table(df_lgd$emp_title), decreasing = TRUE)[1:200])
df_lgd[, emp_title := fifelse(emp_title %in% top_titles, emp_title, "Other")]
df_lgd[, emp_title := relevel(as.factor(emp_title), ref = "Other")]
df_lgd <- one_hot(df_lgd, cols = "emp_title", dropUnusedLevels = TRUE)
df_lgd[, emp_title_Other := NULL]

# "" means no co-borrower and is the reference "None"
df_lgd[, verification_status_joint := as.character(verification_status_joint)]
df_lgd[, verification_status_joint := fifelse(
  verification_status_joint == "", "None", verification_status_joint
)]
df_lgd[, verification_status_joint := relevel(
  as.factor(verification_status_joint),
  ref = "None"
)]
df_lgd <- one_hot(
  df_lgd,
  cols = "verification_status_joint",
  dropUnusedLevels = TRUE
)
df_lgd[, verification_status_joint_None := NULL]

setnames(df_lgd,
  old = names(df_lgd),
  new = gsub(" ", "_", names(df_lgd))
)

# ==============================================================================
# CHECK AND SAVE
# ==============================================================================

nrow(df_lgd)
ncol(df_lgd)
anyNA(df_lgd)
summary(df_lgd$lgd)
mean(df_lgd$lgd == 1)
print(state_counts)

save(df_lgd, file = file.path(data_path, "df_lgd_all.RData"))
cat(
  "Saved:", file.path(data_path, "df_lgd_all.RData"),
  " rows:", nrow(df_lgd),
  " cols:", ncol(df_lgd),
  "\n"
)
