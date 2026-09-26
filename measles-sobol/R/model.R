# model.R
# A toy measles transmission model with births, routine vaccination (two doses)
# and periodic mass vaccination campaigns.
#
# Design choices (kept deliberately simple):
#   - Discrete time, one step = 2 weeks (roughly one measles generation).
#   - Reed-Frost style infection: each susceptible escapes infection with
#     probability exp(-R0 * I / N) per step.
#   - Closed population of constant size: birth rate = death rate.
#   - No age structure. Campaigns act on the whole susceptible pool.
#   - Vectorised: every parameter set is simulated at the same time, so
#     tens of thousands of runs take seconds. This is what makes Sobol' cheap.

# ---- Input ranges -----------------------------------------------------------
# Each row: name, lower bound, upper bound, short description.
# Sources for each range are listed in README.md.
param_ranges <- data.frame(
  name  = c("R0", "birth_rate", "mcv1_cov", "mcv2_cov",
            "mcv1_age", "sia_interval", "sia_reach", "importation"),
  lower = c(12,   15,  0.60, 0.30,  9,  2, 0.20, 0.1),
  upper = c(18,   40,  0.95, 0.90, 12,  5, 0.80, 5.0),
  label = c("Basic reproduction number R0",
            "Crude birth rate (per 1,000 per year)",
            "MCV1 coverage (first dose)",
            "MCV2 coverage (among MCV1 recipients)",
            "Age at MCV1 (months)",
            "Years between campaigns",
            "Campaign reach (share of susceptible people reached)",
            "Imported infections (per million per 2 weeks)"),
  short = c("R0", "Birth rate", "MCV1 coverage", "MCV2 coverage",
            "Age at MCV1", "Years between campaigns", "Campaign reach",
            "Importation rate"),
  stringsAsFactors = FALSE
)

# Scale a matrix of U(0,1) draws to the parameter ranges.
scale_params <- function(U) {
  X <- as.data.frame(U)
  names(X) <- param_ranges$name
  for (j in seq_len(nrow(param_ranges))) {
    X[[j]] <- param_ranges$lower[j] +
      U[, j] * (param_ranges$upper[j] - param_ranges$lower[j])
  }
  X
}

# Efficacy of the first dose depends on age at vaccination: maternal
# antibodies interfere with the vaccine in younger infants. Linear
# interpolation between 77% at 9 months and 92% at 12 months: the median
# field effectiveness for single doses given at 9-11 and >=12 months in
# Uzicanin & Zimmerman (2011), J Infect Dis 204(S1):S133-48.
mcv1_efficacy <- function(age_months) {
  0.77 + (age_months - 9) / 3 * (0.92 - 0.77)
}
MCV2_EFFICACY <- 0.95   # among children the first dose did not protect
SIA_EFFICACY  <- 0.92   # campaign dose given to children mostly >= 12 months

# ---- Simulator --------------------------------------------------------------
# X: data frame of parameter sets (one row per run).
# Returns a list of outputs, one value per run, plus optional trajectories.
simulate_measles <- function(X, years = 30, burn_in = 10,
                             N = 1e6, keep_traj = FALSE) {
  steps_per_year <- 26
  n_steps <- years * steps_per_year
  n <- nrow(X)

  R0   <- X$R0
  mu   <- X$birth_rate / 1000 / steps_per_year      # per-step birth/death rate
  e1   <- mcv1_efficacy(X$mcv1_age)
  # Share of each birth cohort protected by routine vaccination:
  # protected by dose 1, or missed by dose 1 but protected by dose 2.
  v    <- X$mcv1_cov * e1 +
          X$mcv1_cov * X$mcv2_cov * (1 - e1) * MCV2_EFFICACY
  imp  <- X$importation * N / 1e6
  sia_every <- round(X$sia_interval * steps_per_year)
  sia_hit   <- X$sia_reach * SIA_EFFICACY

  # Start near the endemic equilibrium: S = N / R0.
  S <- N / R0
  I <- rep(N * 1e-4, n)
  R <- N - S - I

  record_from <- burn_in * steps_per_year
  cases       <- numeric(n)
  peak_reff   <- rep(0, n)    # highest R0 * S/N seen (worst moment)

  if (keep_traj) {
    traj_S <- matrix(NA_real_, n_steps, n)
    traj_I <- matrix(NA_real_, n_steps, n)
  }

  for (t in seq_len(n_steps)) {
    births  <- mu * N
    new_inf <- S * (1 - exp(-R0 * (I + imp) / N))

    S_next <- S - new_inf + births * (1 - v) - mu * S
    R_next <- R + I + births * v - mu * R
    I      <- new_inf
    S      <- S_next
    R      <- R_next

    # Mass campaign on schedule: move a share of susceptibles to immune.
    campaign <- (t %% sia_every) == 0
    moved    <- ifelse(campaign, S * sia_hit, 0)
    S <- S - moved
    R <- R + moved

    if (t > record_from) {
      cases     <- cases + new_inf
      # Measured before any campaign this step would have run, i.e. at the
      # point where susceptibles have built up the most.
      peak_reff <- pmax(peak_reff, R0 * (S + moved) / N)
    }
    if (keep_traj) {
      traj_S[t, ] <- S / N
      traj_I[t, ] <- I
    }
  }

  recorded_years <- years - burn_in
  out <- list(
    # Mean annual incidence per 100,000 over the recorded years.
    incidence = cases / recorded_years / N * 1e5,
    # Peak effective reproduction number, R0 * S/N. Above 1 means an
    # introduced case would, on average, start a growing outbreak.
    peak_reff = peak_reff
  )
  if (keep_traj) {
    out$traj_S <- traj_S
    out$traj_I <- traj_I
  }
  out
}
