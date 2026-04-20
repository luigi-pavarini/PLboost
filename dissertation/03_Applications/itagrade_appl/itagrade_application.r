library(data.table)
library(mltools)

# Needs to be settled!!!
base_path <- "..."

# ── LOAD ──────────────────────────────────────────────────────────────────────
url <- "https://zenodo.org/records/15423018/files/dataset.csv?download=1"
df  <- as.data.table(read.csv(url))

# ── STEP 1: BUILD RESPONSE y (GPA ratio) ─────────────────────────────────────
# y is each student's GPA divided by the maximum possible GPA.
# It is a continuous variable strictly bounded in (0, 1), which makes
# Beta regression the natural modelling framework.
# We apply a small epsilon correction to avoid exact boundary values,
# though in this dataset y_raw already lies in (0.6677, 0.9900).

gpa_col <- "GPA_ratio"

y_raw <- df[[gpa_col]]

eps <- 1e-4
y   <- pmin(pmax(y_raw, eps), 1 - eps)

df[, y := y]

# GPA is dropped as leakage — it is the raw fraction used to compute y
df[, GPA       := NULL]
df[, GPA_ratio := NULL]

# ── STEP 2: SELECT AND CLASSIFY COLUMNS ───────────────────────────────────────
# Variables described in Mai (2025), Section 5.4:
# Dataset: 171 university students from Milan and Rome.
#
# Demographic characteristics:
# Gender       : sex of the student
# Age          : age in years (continuous)
# City         : city of study — Milano, Roma, Other (3 levels)
# Major        : field of study — Economics, Engineering, Mathematics,
#                Medicine, and others (multiple levels)
#
# Study environment conditions (all binary):
# study_time        : period of study — Morning, Afternoon, Evening, Night
# study_location    : place of study — Home, Library, Café, Other
# enough_sleep      : whether the student sleeps enough
# natural_sun_exposure : whether the student has natural light exposure
# noisy_environment : whether the study environment is noisy
# heated_cooled     : whether heating/cooling is adequate
# ventilated        : whether the environment is well ventilated
# enough_desk_space : whether there is sufficient desk space
# often_distracted  : whether the student is frequently distracted
#                     KEY PREDICTOR — strong negative effect on GPA ratio
#                     (estimate ≈ −0.373 in the Horseshoe model)
# study_in_group    : whether the student studies in groups

pred_cols        <- setdiff(names(df), "y")
cols_categorical <- pred_cols[sapply(df[, ..pred_cols],
                                     function(x) is.character(x) | is.factor(x))]

# ── STEP 3: ONE-HOT ENCODING ──────────────────────────────────────────────────
# Each categorical variable encoded separately.
# Reference = all zeros (no dummy column activated for that category).

# City: Milano, Roma, Other
df[, City := as.factor(City)]
df <- one_hot(df, cols = "City", dropUnusedLevels = TRUE)

# Major: Economics, Engineering, Mathematics, Medicine, etc.
df[, Major := as.factor(Major)]
df <- one_hot(df, cols = "Major", dropUnusedLevels = TRUE)

# study_time: Morning, Afternoon, Evening, Night
if ("study_time" %in% names(df)) {
  df[, study_time := as.factor(study_time)]
  df <- one_hot(df, cols = "study_time", dropUnusedLevels = TRUE)
}

# study_location: Home, Library, Café, Other
if ("study_location" %in% names(df)) {
  df[, study_location := as.factor(study_location)]
  df <- one_hot(df, cols = "study_location", dropUnusedLevels = TRUE)
}

# Gender: encode any remaining character columns generically
remaining_cat <- names(df)[sapply(df, function(x) is.character(x) | is.factor(x))]
remaining_cat <- setdiff(remaining_cat, "y")
for (col in remaining_cat) {
  df[, (col) := as.factor(get(col))]
  df <- one_hot(df, cols = col, dropUnusedLevels = TRUE)
}

# ── CLEAN COLUMN NAMES ────────────────────────────────────────────────────────
setnames(df, old = names(df), new = gsub(" ", "_", names(df)))
setnames(df, old = names(df), new = gsub("[^A-Za-z0-9_]", "_", names(df)))

# ── FINAL CHECK ───────────────────────────────────────────────────────────────
nrow(df)
ncol(df)
anyNA(df)
summary(df$y)
hist(df$y, breaks = 30, main = "GPA ratio distribution", xlab = "y")
plot(density(df$y), main = "Density of GPA ratio")

# ── SAVE ──────────────────────────────────────────────────────────────────────
save(df, file = file.path(base_path, "df_itagrade.RData"))
colnames(df)
