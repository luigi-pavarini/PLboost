library(data.table)
library(mltools)

# Needs to be settled!!!
base_path <- "..."

# ── LOAD ──────────────────────────────────────────────────────────────────────
df <- fread(file.path(base_path, "accepted_2007_to_2018Q4.csv"))

# ── STEP 1: FILTER DEFAULTS ───────────────────────────────────────────────────
# loan_status: keep only loans where the borrower stopped paying
df_default <- df[loan_status %in% c("Charged Off", "Default")]

# ── STEP 2: SELECT COLUMNS ────────────────────────────────────────────────────
# We explicitly select only the columns we need upfront, grouped by purpose.
# Excluded variables are either administrative (id, url), leakage (post-default info),
# or temporal (dates already used to engineer features).
# Leakage variables are those that only exist because of or after the default event —
# including them would mean the model is learning from the outcome itself.

# Used only to compute lgd — recoveries dropped after, loan_amnt kept as predictor
# loan_amnt  : total amount the borrower requested and was approved for (in USD)
# recoveries : amount the lender clawed back after charge-off via collections
#              LEAKAGE: this is only known after the default has occurred and
#              the recovery process has concluded — it is part of what defines LGD
cols_lgd_construction <- c("loan_amnt", "recoveries")

# Used only to engineer credit_age, credit_age_coborrower, days_to_last_pymnt
# All four are dropped once the engineered features are created
# issue_d                 : month and year the loan was funded
# earliest_cr_line        : month and year the borrower opened their first credit account
# sec_app_earliest_cr_line: same as earliest_cr_line but for the co-borrower
# last_pymnt_d            : month and year the borrower made their last payment before default
#                           not leakage — the timing of last payment is known at default time
#                           and tells us how far into the loan term the borrower survived
cols_feat_engineering <- c("issue_d", "earliest_cr_line",
                           "sec_app_earliest_cr_line", "last_pymnt_d")

# Categorical variables for one-hot encoding (12 standard + 2 special = 14 total)
# term                    : loan repayment period — "36 months" or "60 months"
# grade                   : Lending Club risk grade from A (safest) to G (riskiest)
# sub_grade               : finer subdivision of grade, from A1 to G5 (35 levels)
# emp_length              : how long the borrower has been employed at current job
# home_ownership          : housing situation — RENT, OWN, MORTGAGE, OTHER
# verification_status     : whether Lending Club verified the borrower's income
# purpose                 : borrower's stated reason for the loan (debt consolidation, car, etc.)
# addr_state              : US state where the borrower resides — 51 levels
# initial_list_status     : whether the loan was initially offered as fractional or whole
# application_type        : individual application or joint application with a co-borrower
# disbursement_method     : how the borrower received the funds — Cash or DirectPay
# pymnt_plan              : whether a special payment plan was arranged with the borrower
# emp_title               : free-text job title entered by the borrower (special handling)
# verification_status_joint: whether the co-borrower's income was verified (special handling)
#                            empty string "" means no co-borrower — handled as "None" category
cols_categorical <- c("term", "grade", "sub_grade", "emp_length",
                      "home_ownership", "verification_status", "purpose",
                      "addr_state", "initial_list_status", "application_type",
                      "disbursement_method", "pymnt_plan",
                      "emp_title", "verification_status_joint")

# Core loan characteristics — all set at origination, zero leakage risk
# loan_amnt      : total amount the borrower requested (USD)
# funded_amnt    : amount actually funded, which can be less than requested
# funded_amnt_inv: portion of funded amount committed by investors
# int_rate       : annual interest rate assigned to the loan (%)
# installment    : fixed monthly payment amount the borrower owes (USD)
cols_numeric_loan <- c("loan_amnt", "funded_amnt", "funded_amnt_inv",
                       "int_rate", "installment")

# Borrower financial profile — income, debt burden, credit behaviour at origination
# annual_inc            : self-reported annual income of the borrower (USD)
# dti                   : debt-to-income ratio — monthly debt payments / monthly gross income (%)
#                         higher dti means the borrower is already financially stretched
# fico_range_low        : lower bound of the borrower's FICO credit score range at origination
# fico_range_high       : upper bound of the borrower's FICO credit score range at origination
# revol_bal             : total revolving credit balance outstanding (USD)
# revol_util            : revolving line utilisation rate — balance / credit limit (%)
# delinq_2yrs           : number of delinquencies (30+ days late) in the past 2 years
# inq_last_6mths        : number of credit inquiries in the last 6 months
#                         high inquiries suggest the borrower is actively seeking more credit
# open_acc              : number of currently open credit lines
# pub_rec               : number of derogatory public records (bankruptcies, liens, judgements)
# total_acc             : total number of credit lines ever opened
# acc_now_delinq        : number of accounts currently in delinquency
# delinq_amnt           : total dollar amount currently past due on delinquent accounts (USD)
# mths_since_last_delinq: months since the most recent delinquency (NA if never delinquent)
# mths_since_last_record: months since the most recent public record (NA if none on file)
# collections_12_mths_ex_med : number of collections in last 12 months excluding medical
# mths_since_last_major_derog: months since most recent 90-day or worse derogatory rating
cols_numeric_borrower <- c("annual_inc", "dti",
                           "fico_range_low", "fico_range_high",
                           "revol_bal", "revol_util",
                           "delinq_2yrs", "inq_last_6mths",
                           "open_acc", "pub_rec", "total_acc",
                           "acc_now_delinq", "delinq_amnt",
                           "mths_since_last_delinq", "mths_since_last_record",
                           "collections_12_mths_ex_med",
                           "mths_since_last_major_derog")

# Extended credit bureau variables — deeper view of credit portfolio composition
# tot_coll_amt           : total amount ever owed on collections accounts (USD)
# tot_cur_bal            : total current balance across all accounts (USD)
# total_rev_hi_lim       : total revolving high credit / credit limit (USD)
# avg_cur_bal            : average current balance across all open accounts (USD)
# bc_open_to_buy         : total open-to-buy on bankcard accounts (USD) — available credit
# bc_util                : ratio of bankcard balances to credit limits (%)
# all_util               : balance-to-limit ratio across all credit lines (%)
# il_util                : ratio of balance to high credit on instalment accounts (%)
# mort_acc               : number of mortgage accounts
# pub_rec_bankruptcies   : number of public record bankruptcies
# tax_liens              : number of tax liens on record
# tot_hi_cred_lim        : total high credit / credit limit across all accounts (USD)
# total_bal_ex_mort      : total credit balance excluding mortgage (USD)
# total_bc_limit         : total bankcard credit limit (USD)
# total_il_high_credit_limit : total instalment high credit / credit limit (USD)
# acc_open_past_24mths   : number of accounts opened in the past 24 months
# pct_tl_nvr_dlq         : percentage of trades that have never been delinquent (%)
# percent_bc_gt_75       : percentage of bankcard accounts with utilisation > 75% (%)
# num_accts_ever_120_pd  : number of accounts ever 120+ days past due
# num_actv_bc_tl         : number of currently active bankcard tradelines
# num_actv_rev_tl        : number of currently active revolving tradelines
# num_bc_sats            : number of satisfactory bankcard accounts
# num_bc_tl              : number of bankcard accounts
# num_il_tl              : number of instalment accounts
# num_op_rev_tl          : number of open revolving accounts
# num_rev_accts          : number of revolving accounts
# num_rev_tl_bal_gt_0    : number of revolving tradelines with balance > 0
# num_sats               : number of satisfactory accounts
# num_tl_90g_dpd_24m     : number of accounts 90+ days past due in last 24 months
# num_tl_op_past_12m     : number of accounts opened in the past 12 months
# mo_sin_old_il_acct     : months since oldest instalment account opened
# mo_sin_old_rev_tl_op   : months since oldest revolving account opened
# mo_sin_rcnt_rev_tl_op  : months since most recent revolving account opened
# mo_sin_rcnt_tl         : months since most recent account opened
# mths_since_recent_bc   : months since most recent bankcard account opened
# mths_since_recent_inq  : months since most recent credit inquiry
# inq_fi                 : number of personal finance inquiries
# inq_last_12m           : number of credit inquiries in last 12 months
# total_cu_tl            : number of finance trades (credit union tradelines)
# max_bal_bc             : maximum current balance owed on all bankcard accounts (USD)
# open_acc_6m            : number of new accounts opened in last 6 months
# open_act_il            : number of currently active instalment trades
# open_il_12m            : number of instalment accounts opened in past 12 months
# open_il_24m            : number of instalment accounts opened in past 24 months
# open_rv_12m            : number of revolving accounts opened in past 12 months
# open_rv_24m            : number of revolving accounts opened in past 24 months
# total_bal_il           : total current balance of all instalment accounts (USD)
# chargeoff_within_12_mths : number of charge-offs within last 12 months
cols_numeric_bureau <- c("tot_coll_amt", "tot_cur_bal",
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
                         "total_bal_il", "chargeoff_within_12_mths")

# Co-borrower variables — only populated for joint applications
# For individual applications these are NA or empty, handled during encoding
# annual_inc_joint              : combined annual income of both borrowers (USD)
# dti_joint                     : combined debt-to-income ratio for both borrowers (%)
# revol_bal_joint               : combined revolving balance for both borrowers (USD)
# sec_app_fico_range_low        : lower bound of co-borrower's FICO score range
# sec_app_fico_range_high       : upper bound of co-borrower's FICO score range
# sec_app_inq_last_6mths        : number of credit inquiries for co-borrower in last 6 months
# sec_app_mort_acc              : number of mortgage accounts for co-borrower
# sec_app_open_acc              : number of open credit lines for co-borrower
# sec_app_revol_util            : revolving utilisation rate for co-borrower (%)
# sec_app_open_act_il           : number of active instalment trades for co-borrower
# sec_app_num_rev_accts         : number of revolving accounts for co-borrower
# sec_app_chargeoff_within_12_mths     : co-borrower charge-offs within last 12 months
# sec_app_collections_12_mths_ex_med   : co-borrower collections in last 12 months excl. medical
cols_numeric_joint <- c("annual_inc_joint", "dti_joint", "revol_bal_joint",
                        "sec_app_fico_range_low", "sec_app_fico_range_high",
                        "sec_app_inq_last_6mths", "sec_app_mort_acc",
                        "sec_app_open_acc", "sec_app_revol_util",
                        "sec_app_open_act_il", "sec_app_num_rev_accts",
                        "sec_app_chargeoff_within_12_mths",
                        "sec_app_collections_12_mths_ex_med")

# last_pymnt_amnt: the dollar amount of the borrower's last payment before default (USD)
#                  not leakage — it reflects the borrower's behaviour just before stopping,
#                  e.g. a very small last payment may signal the borrower was already struggling
cols_numeric_other <- c("last_pymnt_amnt")

# Combine all groups — unique() handles loan_amnt appearing in both
# cols_lgd_construction and cols_numeric_loan
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

# Safety net: only keep columns that actually exist in this version of the dataset
cols_select <- cols_select[cols_select %in% names(df_default)]
df_selected <- df_default[, ..cols_select]

# ── STEP 3: DROP HIGH-NA NUMERIC COLUMNS, THEN DROP NA ROWS ─────────────────
# Text/categorical columns are never touched here — empty strings or NAs in those
# carry meaning (e.g. verification_status_joint = "" means no co-borrower)
# and are handled explicitly during one-hot encoding in Step 8.
#
# For numeric columns we apply a two-stage strategy:
#   Stage A — drop columns with >10% NA: these have structural missingness,
#             meaning the variable is simply not collected for most borrowers
#             (e.g. all sec_app_* joint variables are NA for 98% of individual loans,
#             il_util and all_util are missing for ~55% of borrowers).
#             Keeping them would force us to drop the majority of rows.
#   Stage B — drop rows where any remaining numeric column has NA: after removing
#             the high-NA columns, the survivors have at most ~7% NA (the bureau
#             block like tot_coll_amt, num_* etc.), so row loss is manageable.
num_cols_selected <- names(df_selected)[sapply(df_selected, is.numeric)]

# Stage A: identify and drop numeric columns with more than 10% NA
na_frac          <- colMeans(is.na(df_selected[, ..num_cols_selected]))
high_na_cols     <- names(na_frac[na_frac > 0.10])
df_selected[, (high_na_cols) := NULL]

# Stage B: drop rows where any remaining numeric column still has NA
num_cols_clean   <- names(df_selected)[sapply(df_selected, is.numeric)]
df_selected      <- df_selected[complete.cases(df_selected[, ..num_cols_clean])]

# Stage C: remove implausible annual_inc = 0
# An approved Lending Club loan with zero reported income is not credible —
# these are miscoded missing values, not genuine zero-income borrowers.
# Removing here, before feature engineering, prevents Inf/NaN in the
# income-ratio interaction terms constructed in Step 7.
cat("Rows with annual_inc = 0:", sum(df_selected$annual_inc == 0), "\n")
df_selected <- df_selected[annual_inc > 0]


# ── STEP 4: LGD CONSTRUCTION ──────────────────────────────────────────────────
# lgd = (loan_amnt - recoveries) / loan_amnt
# lgd = 1   means nothing was recovered — total loss
# lgd = 0.3 means 70% was recovered, only 30% lost
df_selected[, lgd := (loan_amnt - recoveries) / loan_amnt]

hist(df_selected$lgd, breaks = 50, main = "LGD before cleaning")
summary(df_selected$lgd)

# Keep only (0, 1]
# lgd <= 0: lender recovered more than original amount — not a real loss
# lgd >  1: loss exceeded loan amount (collection costs) — excluded
df_lgd <- df_selected[lgd > 0 & lgd <= 1]

hist(df_lgd$lgd, breaks = 50, main = "LGD after cleaning to (0,1]")
summary(df_lgd$lgd)

# Drop recoveries — used for LGD construction only, not a predictor
df_lgd[, recoveries := NULL]

# ── STEP 5: FEATURE ENGINEERING ───────────────────────────────────────────────

# credit_age: years of credit history at loan issue date
# earliest_cr_line to issue_d, leap-year adjusted (÷365.25)
df_lgd[, credit_age := as.numeric(difftime(
  as.Date(paste0("01-", issue_d),          format = "%d-%b-%Y"),
  as.Date(paste0("01-", earliest_cr_line), format = "%d-%b-%Y"),
  units = "days")) / 365.25]

# credit_age_coborrower: same for co-borrower
# sec_app_earliest_cr_line = "" means individual application → 0
# 0 signals "no co-borrower", not "zero years of history"
df_lgd[, sec_app_earliest_cr_line := as.character(sec_app_earliest_cr_line)]
df_lgd[, credit_age_coborrower := fifelse(
  sec_app_earliest_cr_line == "", 0,
  as.numeric(difftime(
    as.Date(paste0("01-", issue_d),                  format = "%d-%b-%Y"),
    as.Date(paste0("01-", sec_app_earliest_cr_line), format = "%d-%b-%Y"),
    units = "days")) / 365.25)]

# days_to_last_pymnt: days between loan issue and last payment before default
# measures how quickly the borrower stopped paying — known at default, not leakage
df_lgd[, days_to_last_pymnt := as.numeric(difftime(
  as.Date(paste0("01-", last_pymnt_d), format = "%d-%b-%Y"),
  as.Date(paste0("01-", issue_d),      format = "%d-%b-%Y"),
  units = "days"))]

# Drop date columns — feature engineering complete
df_lgd[, c("issue_d", "earliest_cr_line",
           "sec_app_earliest_cr_line", "last_pymnt_d") := NULL]

# Drop rows where days_to_last_pymnt is NA
# These correspond to defaulted loans with no last payment date recorded (last_pymnt_d = "")
# Only ~1,747 rows (<0.1% of data) — negligible loss, genuinely ambiguous records
# inter_days_loan inherits the same NAs so both are fixed by this single drop
cat("Rows dropped due to missing days_to_last_pymnt:", sum(is.na(df_lgd$days_to_last_pymnt)), "\n")
df_lgd <- df_lgd[!is.na(days_to_last_pymnt)]

# ── STEP 6: DROP NEAR-ZERO COLUMNS ────────────────────────────────────────────
# Numeric columns where >95% of values are exactly zero carry almost no signal
# and inflate the feature space unnecessarily
num_cols_check <- names(df_lgd)[sapply(df_lgd, is.numeric)]
num_cols_check <- num_cols_check[num_cols_check != "lgd"]

zero_frac  <- sapply(num_cols_check, function(col) mean(df_lgd[[col]] == 0, na.rm = TRUE))
near_zero  <- names(zero_frac[zero_frac > 0.95])

print(near_zero)

df_lgd[, (near_zero) := NULL]

# ── STEP 7: INTERACTION TERMS ─────────────────────────────────────────────────

# total cost burden: large loan at high rate is far riskier than either alone
df_lgd[, inter_loan_intrate     := loan_amnt      * int_rate]

# loan-to-income ratio: how large is the loan relative to annual earnings
df_lgd[, inter_loan_income      := loan_amnt      / annual_inc]

# monthly payment burden: share of income consumed by this loan installment
df_lgd[, inter_install_income   := installment    / annual_inc]

# indebted borrower at high interest: double pressure on cash flow
df_lgd[, inter_dti_intrate      := dti            * int_rate]

# poor credit quality + expensive borrowing: compounding risk signal
df_lgd[, inter_fico_intrate     := fico_range_low * int_rate]

# revolving debt burden relative to income
df_lgd[, inter_revol_income     := revol_bal      / annual_inc]

# long credit history but with derogatory public records
df_lgd[, inter_creditage_pubrec := credit_age     * pub_rec]

# actively seeking more credit while already heavily indebted
df_lgd[, inter_inq_dti          := inq_last_6mths * dti]

# quick default on large loan vs slow default on small loan
df_lgd[, inter_days_loan        := days_to_last_pymnt * loan_amnt]

# ── STEP 8: ONE-HOT ENCODING ──────────────────────────────────────────────────
# Each variable encoded separately to control the reference category.
# Reference = all zeros (no dummy column activated for that category).

# term: "36 months" or "60 months"
df_lgd[, term := as.factor(term)]
df_lgd <- one_hot(df_lgd, cols = "term", dropUnusedLevels = TRUE)

# grade: A (lowest risk) to G (highest risk)
df_lgd[, grade := as.factor(grade)]
df_lgd <- one_hot(df_lgd, cols = "grade", dropUnusedLevels = TRUE)

# sub_grade: A1 to G5 (up to 35 levels)
df_lgd[, sub_grade := as.factor(sub_grade)]
df_lgd <- one_hot(df_lgd, cols = "sub_grade", dropUnusedLevels = TRUE)

# emp_length: employment duration
df_lgd[, emp_length := as.factor(emp_length)]
df_lgd <- one_hot(df_lgd, cols = "emp_length", dropUnusedLevels = TRUE)

# home_ownership: RENT, OWN, MORTGAGE, OTHER
df_lgd[, home_ownership := as.factor(home_ownership)]
df_lgd <- one_hot(df_lgd, cols = "home_ownership", dropUnusedLevels = TRUE)

# verification_status: Not Verified, Source Verified, Verified
df_lgd[, verification_status := as.factor(verification_status)]
df_lgd <- one_hot(df_lgd, cols = "verification_status", dropUnusedLevels = TRUE)

# purpose: stated loan purpose (car, credit_card, debt_consolidation, etc.)
df_lgd[, purpose := as.factor(purpose)]
df_lgd <- one_hot(df_lgd, cols = "purpose", dropUnusedLevels = TRUE)

# addr_state: US state of residence (up to 51 levels)
# Snapshot state counts before encoding — original column is dropped by one_hot()
state_counts <- sort(table(df_lgd$addr_state), decreasing = TRUE)
df_lgd[, addr_state := as.factor(addr_state)]
df_lgd <- one_hot(df_lgd, cols = "addr_state", dropUnusedLevels = TRUE)

# initial_list_status: "f" = fractional, "w" = whole
df_lgd[, initial_list_status := as.factor(initial_list_status)]
df_lgd <- one_hot(df_lgd, cols = "initial_list_status", dropUnusedLevels = TRUE)

# application_type: Individual or Joint App
df_lgd[, application_type := as.factor(application_type)]
df_lgd <- one_hot(df_lgd, cols = "application_type", dropUnusedLevels = TRUE)

# disbursement_method: Cash or DirectPay
df_lgd[, disbursement_method := as.factor(disbursement_method)]
df_lgd <- one_hot(df_lgd, cols = "disbursement_method", dropUnusedLevels = TRUE)

# pymnt_plan: payment plan arranged — "y" or "n"
df_lgd[, pymnt_plan := as.factor(pymnt_plan)]
df_lgd <- one_hot(df_lgd, cols = "pymnt_plan", dropUnusedLevels = TRUE)

# emp_title: free-text job title — top 200 kept, rest → "Other" (reference = all zeros)
# NA and empty string both mean the borrower did not report a job title → "Other"
# NA must be replaced before the %in% check because NA %in% top_titles returns NA,
# not FALSE, which would propagate NAs and break relevel()
df_lgd[, emp_title := as.character(emp_title)]
df_lgd[is.na(emp_title) | emp_title == "", emp_title := "Other"]
top_titles <- names(sort(table(df_lgd$emp_title), decreasing = TRUE)[1:200])
df_lgd[, emp_title := fifelse(emp_title %in% top_titles, emp_title, "Other")]
df_lgd[, emp_title := relevel(as.factor(emp_title), ref = "Other")]
df_lgd <- one_hot(df_lgd, cols = "emp_title", dropUnusedLevels = TRUE)
df_lgd[, emp_title_Other := NULL]

# verification_status_joint: co-borrower income verification
# "" = no co-borrower → replaced with "None" as reference category
df_lgd[, verification_status_joint := as.character(verification_status_joint)]
df_lgd[, verification_status_joint := fifelse(
  verification_status_joint == "", "None", verification_status_joint)]
df_lgd[, verification_status_joint := relevel(
  as.factor(verification_status_joint), ref = "None")]
df_lgd <- one_hot(df_lgd, cols = "verification_status_joint", dropUnusedLevels = TRUE)
df_lgd[, verification_status_joint_None := NULL]

# ── CLEAN COLUMN NAMES ───────────────────────────────────────────────────────
# One-hot encoding produces column names with spaces (e.g. "term_36 months",
# "application_type_Joint App") because the original factor levels had spaces.
# We replace all spaces with underscores so every column name is safe to use
# in R formulas, data.table syntax, and downstream modelling functions.
setnames(df_lgd, old = names(df_lgd),
         new = gsub(" ", "_", names(df_lgd)))

# ── FINAL CHECK ───────────────────────────────────────────────────────────────
nrow(df_lgd)
ncol(df_lgd)
anyNA(df_lgd)
summary(df_lgd$lgd)
mean(df_lgd$lgd == 1)
hist(df_lgd$lgd, breaks = 50)
plot(density(df_lgd$lgd))

# Row count per state — computed before one-hot encoding in Step 8
print(state_counts)

# ── SAVE (ALL STATES) ────────────────────────────────────────────────────────
save(df_lgd, file = file.path(base_path, "df_lgd_all.RData"))
