# HYBRID SIM_SURVEY FUNCTION (updated June 12, 2025)
# Simulates trawl survey sampling by LENGTH (as in real surveys), while retaining true age.
# Modifications:
#     - simulate individuals based on abundance-at-age (sp_N)
#     - assign biological length via sim_length()
#     - sample LENGTHS, then subsample AGES from length bins
#     - minimize memory use with vectorized simulation loop
#     - preserve sim/year/set structure throughout

#' Closure for simulating logistic curve
#'
#' @description This closure is useful for simulating q inside the
#'              \code{\link{sim_survey}} function
#'
#' @param k      The steepness of the curve
#' @param x0     The x-value of the sigmoid's midpoint
#' @param plot   Plot relationship
#'
#' @return Returns a function for use in \code{\link{sim_survey}}.
#'
#' @examples
#' logistic_fun <- sim_logistic(k = 2, x0 = 3, plot = TRUE)
#' logistic_fun(x = 1:10)
#'
#' @export
#'

sim_logistic <- function(k = 2, x0 = 3, plot = FALSE) {
  function(x = NULL) {
    y <- 1 / (1 + exp(-k * (x - x0)))
    if (plot) plot(x, y, type = "b")
    y
  }
}

#' Function for deterministic distribution of I_at_length
#'
#' @description This function ...

selectivity_integral <- function(sim, l50=15, l95=20, cv = 0.1){

  ## Define logistic selectivity function based on l50 and l95
  q_length <- function(length) {
    k <- log(19) / (l95 - l50)  # Steepness
    x0 <- l50  # Midpoint
    1 / (1 + exp(-k * (length - x0)))
  }

  ## Define probability density function for length given age with sim$sim_length
  prob_len_given_age <- function(length, age){
    len_mean <- sim$sim_length(age) # return mean length at age
    len_sd <- len_mean * cv
    dnorm(length, mean=len_mean, sd=len_sd)
  }

  # Reconstruct length bins from the simulation’s length-age key
  lak <- sim$sim_length(age = sim$ages, length_age_key = TRUE)
  bin_midpoints <- as.numeric(rownames(lak))
  bin_width <- diff(bin_midpoints)[1]
  length_bins <- c(0, bin_midpoints + bin_width / 2)

  # Initialize output matrix
  years <- sim$years
  I_at_length_det <- matrix(0,
                            nrow = length(bin_midpoints),
                            ncol = length(years),
                            dimnames = list(length = as.character(bin_midpoints),
                                            year = as.character(years)))

  # Loop over bins and ages
  for (age in sim$ages) {
    for (i in seq_along(length_bins)[-length(length_bins)]) {
      bin_lower <- length_bins[i]
      bin_upper <- length_bins[i + 1]

      integral <- integrate(
        function(l) q_length(l) * prob_len_given_age(l, age),
        lower = bin_lower, upper = bin_upper
      )$value

      I_at_length_det[i, ] <- I_at_length_det[i, ] + sim$N[as.character(age), ] * integral
    }
  }

  return(I_at_length_det)
}


#' Round simulated population
#'
#' @param sim Simulation from \code{\link{sim_distribution}}
#'
#' @return Returns a rounded simulation object. Largely used as a helper in \code{\link{sim_survey}}.
#'
#' @export
#'

round_sim <- function(sim) {
  sim$sp_N$N <- round(sim$sp_N$N)
  N <- tapply(sim$sp_N$N, list(sim$sp_N$age, sim$sp_N$year), sum)
  N <- N[rownames(sim$N), colnames(sim$N)]
  dimnames(N) <- dimnames(sim$N)
  sim$N_at_length <- convert_N(N_at_age = N,
                               lak = sim$sim_length(age = sim$ages, length_age_key = TRUE))
  sim$N <- N
  sim$N0 <- N[, 1]
  sim$R <- N[1, ]
  sim
}



#' Simulate survey sets
#'
#' @param sim             Simulation object from \code{\link{sim_distribution}}
#' @param subset_cells    Logical expression indicating the elements (\code{x, y, depth, cell,
#'                        division, strat, year}) of the survey grid to keep (e.g., \code{cell
#'                        < 100})
#' @param n_sims          Number of simulations to produce
#' @param trawl_dim       Trawl width and distance (same units as grid)
#' @param min_sets        Minimum number of sets per strat
#' @param set_den         Set density (number of sets per grid unit squared)
#' @param resample_cells  Allow resampling of sampling units (grid cells)?
#'                        (Note: allowing resampling may create bias because
#'                        depletion is imposed at the cell level)
#'
#' @return Returns a data.table including details of each set location.
#'
#' @export
#'
#' @examples
#'
#'\donttest{
#' sim <- sim_abundance(ages = 1:5, years = 1:5) %>%
#'           sim_distribution(grid = make_grid(res = c(20, 20)))
#'
#' ## Multiple calls can be useful for defining a custom series of sets
#' standard_sets <- sim_sets(sim, year <= 2, set_den = 2 / 1000)
#' reduced_sets <- sim_sets(sim, year > 2 & !cell %in% 1:100, set_den = 1 / 1000)
#' sets <- rbind(standard_sets, reduced_sets)
#' sets$set <- seq(nrow(sets)) # Important - make sure set has a unique ID.
#'
#' survey <- sim_survey(sim, custom_sets = sets)
#'
#' plot_survey(survey, which_year = 3, which_sim = 1)
#' }
#'

sim_sets <- function(sim, subset_cells, n_sims = 1, trawl_dim = c(1.5, 0.02),
                     min_sets = 2, set_den = 2 / 1000,
                     resample_cells = FALSE) {

  strat_sets <- cell_sets <- NULL
  cells <- data.table(data.frame(sim$grid))

  ## Replicate cells data.table for each year in the simulation
  i <- rep(seq(nrow(cells)), times = length(sim$years))
  y <- rep(sim$years, each = nrow(cells))
  cells <- cells[i, ]
  cells$year <- y

  ## Replicate n_sims times
  i <- rep(seq(nrow(cells)), times = n_sims)
  s <- rep(seq(n_sims), each = nrow(cells))
  cells <- cells[i, ]
  cells$sim <- s

  ## Subset cells for sampling
  if (!missing(subset_cells)) {
    r <- eval(substitute(subset_cells), cells)
    cells <- cells[r,]
  }

  ## Strat area and sampling effort
  strat_det <- cells[, list(strat_cells = .N), by = c("sim", "year", "strat")]
  strat_det$tow_area <- prod(trawl_dim)
  strat_det$cell_area <- prod(stars::st_res(sim$grid))
  strat_det$strat_area <- strat_det$strat_cells * prod(stars::st_res(sim$grid))
  strat_det$strat_sets <- round(strat_det$strat_area * set_den) # set allocation
  strat_det$strat_sets[strat_det$strat_sets < min_sets] <- min_sets
  cells <- merge(cells, strat_det, by = c("sim", "year", "strat"))

  ## Simulate sets; randomly sample row id by group
  ind <- cells[, .I[sample(.N, size = unique(strat_sets), replace = resample_cells)],
               by = c("sim", "year", "strat")][[4]]
  sets <- cells[ind, ]
  sets[, cell_sets := .N, by = c("sim", "year", "cell")] # useful for identifying cells with more than one set (when resample_units = TRUE)
  sets$set <- seq(nrow(sets))
  sets

}


#' Simulate stratified-random survey
#'
#' @param sim                 Simulation from \code{\link{sim_distribution}}
#' @param n_sims              Number of surveys to simulate over the simulated population. Note: requesting
#'                            a large number of simulations may max out your RAM. Use
#'                            \code{\link{sim_survey_parallel}} if many simulations are required.
#' @param q                   Closure, such as \code{\link{sim_logistic}}, for simulating catchability at age
#'                            (returned values must be between 0 and 1)
#' @param trawl_dim           Trawl width and distance (same units as grid)
#' @param resample_cells      Allow resampling of sampling units (grid cells)? Setting to TRUE may introduce bias
#'                            because depletion is imposed at the cell level.
#' @param binom_error         Impose binomial error? Setting to FALSE may introduce bias in stratified estimates
#'                            at older ages because of more frequent rounding to zero.
#' @param min_sets            Minimum number of sets per strat
#' @param set_den             Set density (number of sets per grid unit squared). WARNING:
#'                            may return an error if \code{set_den} is high and
#'                            \code{resample_cells = FALSE} because the number of sets allocated may
#'                            exceed the number of cells in a strata.
#' @param lengths_cap         Maximum number of lengths measured per set
#' @param ages_cap            If \code{age_sampling = "stratified"}, this cap represents the maximum
#'                            number of ages to sample per length group (defined using the \code{age_length_group}
#'                            argument) per division or strat (defined using the \code{age_space_group} argument)
#'                            per year. If \code{age_sampling = "random"}, it is the maximum number of ages to sample
#'                            from measured fish per set.
#' @param age_sampling        Should age sampling be "stratified" (default) or "random"?
#' @param age_length_group    Numeric value indicating the size of the length bins for stratified
#'                            age sampling. Ignored if \code{age_sampling = "random"}.
#' @param age_space_group     Should age sampling occur at the "division" (default), "strat" or "set" spatial scale?
#'                            That is, age sampling can be spread across each "division", "strat" or "set"
#'                            in each year to a maximum number within each length bin (cap is defined using
#'                            the \code{age_cap} argument). Ignored if \code{age_sampling = "random"}.
#' @param custom_sets         Supply an object of the same structure as returned by \code{\link{sim_sets}} which
#'                            specifies a custom series of set locations to be sampled. Set locations are
#'                            automated if \code{custom_sets = NULL}.
#' @param light               Drop some objects from the output to keep object size low?
#'
#' @return A list including rounded population simulation, set locations and details
#' and sampling details. Note that that N = "true" population, I = population available
#' to the survey, n = number caught by survey.
#'
#' @examples
#'
#'\donttest{
#' sim <- sim_abundance(ages = 1:5, years = 1:5) %>%
#'            sim_distribution(grid = make_grid(res = c(20, 20))) %>%
#'            sim_survey(n_sims = 5, q = sim_logistic(k = 2, x0 = 3))
#' plot_survey(sim, which_year = 3, which_sim = 1)
#' }
#'
#' @export
#'

sim_survey_hybrid <- function(sim, n_sims = 1,
                              q = sim_logistic(),
                              q_length = sim_logistic(k=0.2, x0=20),
                              trawl_dim = c(1.5, 0.02),
                              resample_cells = FALSE,
                              binom_error = TRUE,
                              min_sets = 2,
                              set_den = 2 / 1000, lengths_cap = 500,
                              ages_cap = 10,
                              age_sampling = "stratified",
                              age_length_group = 1,
                              select_by_age = TRUE,
                              l50 = 15,
                              l95 = 20,
                              length_max = 120,
                              age_space_group = "division",
                              custom_sets = NULL,
                              light = TRUE) {

  n <- age <- id <- division <- strat <- N <- n_measured <- n_aged <- NULL

  ## Couple error traps
  if (!age_sampling %in% c("stratified", "random")) {
    stop('age_sampling must be either "stratified" or "random". Other options have yet to be implemented.')
  }
  if (age_sampling == "random" && ages_cap > lengths_cap) {
    stop('When age_sampling = "random", ages_cap cannot exceed lengths_cap.')
  }
  if (!age_space_group %in% c("division", "strat", "set")) {
    stop('age_space_group must be either "division", "strat" or "set". Other options have yet to be implemented.')
  }

  ## Round simulated population and calculate numbers available to survey
  sim <- round_sim(sim)

  ## Simulate sets conducted across survey grid
  if (is.null(custom_sets)) {
    sets <- sim_sets(sim, resample_cells = resample_cells, n_sims = n_sims,
                     trawl_dim = trawl_dim, set_den = set_den, min_sets = min_sets)
  } else {
    sets <- as.data.table(custom_sets)
    if (any(duplicated(sets$set))) {
      stop("When supplying 'custom_sets', please make sure each set has a unique number.")
    }
  }
  setkeyv(sets, c("sim", "year", "cell"))
  lak <- sim$sim_length(age = sim$ages, length_age_key = TRUE)

  I <- sim$N * q(replicate(length(sim$years), sim$ages))

  if (select_by_age ==TRUE) {

    ### SAMPLE BY AGE as in original sim_survey###

    ## Deterministic from convert_N
    I_at_length <- convert_N(N_at_age = I, lak = lak)
    sim$I_at_length <- I_at_length

    sim$sp_N$I <- sim$sp_N$N * q(sim$sp_N$age)

    ## Expand sp_N object n_sim times
    sp_I <- data.table(sim$sp_N[, c("cell", "age", "year", "N")])
    i <- rep(seq(nrow(sp_I)), times = n_sims)
    s <- rep(seq(n_sims), each = nrow(sp_I))
    sp_I <- sp_I[i, ]
    sp_I$sim <- s
    ## Subset population to surveyed cells and simulate portion caught by survey
    ## Introduce sampling error using rbinom
    ## (If more than one set is conducted in a cell, split population available to survey (I) amongst the sets)
    setdet <- merge(sets, sp_I, by = c("sim", "year", "cell"))

    if (binom_error) {
      setdet$n <- stats::rbinom(rep(1, nrow(setdet)), size = round(setdet$N / setdet$cell_sets),
                                prob = (setdet$tow_area / setdet$cell_area) * q(setdet$age))
    } else {
      setdet$n <- round((setdet$N / setdet$cell_sets) * ((setdet$tow_area / setdet$cell_area) * q(setdet$age)))
    }
    setkeyv(setdet, "set")
    setkeyv(sets, "set")
    rm(sp_I)

    ## Expand set catch to individuals and simulate length
    samp <- setdet[rep(seq(.N), n), list(set, age)]
    samp$id <- seq(nrow(samp))
    samp$length <- sim$sim_length(samp$age)

    ## Sample lengths
    measured <- samp[, list(id = id[sample(.N, ifelse(.N > lengths_cap, lengths_cap, .N),
                                           replace = FALSE)]), by = "set"]
    samp$measured <- samp$id %in% measured$id # tag lengths collected
    length_samp <- samp[samp$measured, ]
    rm(measured)

    ## Sample ages
    # length_samp$length_group <- group_lengths(length_samp$length, age_length_group)
    # length_samp <- merge(sets[, list(set, sim, year, division, strat)], length_samp,
    #                      by = c("set", "sim", "year"))
    length_samp <- merge(sets[, .(set, sim, year, division, strat)], length_samp, by = "set")
    length_samp$length_group <- group_lengths(length_samp$length, age_length_group)

    if (age_sampling == "stratified") {
      aged <- length_samp[, list(id = id[sample(.N, ifelse(.N > ages_cap, ages_cap, .N),
                                                replace = FALSE)]),
                          by = c("sim", "year", age_space_group, "length_group")]
    }
    if (age_sampling == "random") {
      aged <- length_samp[, list(id = id[sample(.N, ifelse(.N > ages_cap, ages_cap, .N),
                                                replace = FALSE)]),
                          by = c("set")]
    }
    samp$aged <- samp$id %in% aged$id # tag ages sampled
    rm(aged)
    rm(length_samp)

    ## Simplify samp object
    samp <- samp[, list(set, id, length, age, measured, aged)]
    if (light) samp$id <- NULL

    ## Summarise set catch and sampling
    if (!light) full_setdet <- setdet
    setdet <- merge(sets, setdet[, list(N = sum(N), n = sum(n)), by = "set"], by = "set")
    setdet <- merge(setdet,
                    samp[, list(n_measured = sum(measured), n_aged = sum(aged)), by = "set"],
                    by = "set", all.x = TRUE)
  } else {

    ### SAMPLE BY LENGTH ###

    ## Expand age-based abundance to individual fish, one row per fish
    sp_N <- as.data.table(sim$sp_N)[round(N) > 0]
    sp_N <- sp_N[rep(1:.N, times = n_sims)]
    sp_N[, sim := rep(seq_len(n_sims), each = .N / n_sims)]

    ## Rejoin spatial grid info
    grid_info <- as.data.table(sim$grid)

    ## Add cell_area from grid resolution
    cell_area_val <- prod(as.numeric(stars::st_res(sim$grid)))
    grid_info[, cell_area := cell_area_val]

    sp_N <- merge(sp_N, grid_info[, .(cell, cell_area)],
                  by = "cell", all.x = TRUE)

    ## Merge sp_N with sets to assign sampling info
    sp_N <- merge(sp_N, sets[, .(sim, year, cell, set, tow_area, cell_sets, x, y, division, strat)],
                  by = c("sim","year","cell"))

    ## Initialize I_at_length tally
    length_labels <- as.numeric(rownames(lak))
    length_bins <- c(0, unique(sort(length_labels + diff(length_labels)[1]/2)))  # define breaks

    # Initialize tallies as matrices (length bins × years) and (ages × years)
    I_at_length <- matrix(0,
                          nrow = length(length_labels),
                          ncol = length(unique(sp_N$year)))
    dimnames(I_at_length) <- list(length = as.character(length_labels),
                                  year = as.character(unique(sp_N$year)))

    # Sampling loop over each population cell
    samp_list <- lapply(1:nrow(sp_N), function(i) {
      row <- sp_N[i]
      N_fish <- round(row$N)
      if (N_fish == 0) return(NULL)

      # Simulate individual ages
      ages <- rep(row$age, N_fish)
      lengths <- sim$sim_length(ages)
      year <- as.character(row$year)

      # Bin by numeric lower bounds
      binned_lengths <- cut(lengths,
                            breaks = length_bins,
                            labels = as.character(length_labels),
                            include.lowest = TRUE,
                            right = FALSE)
      len_table <- table(binned_lengths)
      I_at_length[names(len_table), year] <<- I_at_length[names(len_table), year] + as.numeric(len_table)

      ## Simulate catch
      tow_ratio <- as.numeric(row$tow_area) / as.numeric(row$cell_area)
      catch_probs <- tow_ratio * q_length(lengths)
      caught <- rbinom(N_fish, 1, catch_probs)
      if (sum(caught) == 0) return(NULL)

      data.table(
        set = row$set,
        sim = row$sim,
        year = row$year,
        division = row$division,
        strat = row$strat,
        age = ages[caught == 1],
        length = lengths[caught == 1]
      )
    })

    samp <- rbindlist(samp_list)
    samp[, id := .I]
    samp[, length_group := group_lengths(length, age_length_group)]

    ## Subsample measured lengths
    measured <- samp[, if (.N > 0) .(id = sample(id, min(.N, lengths_cap)))
                     else .(id = integer(0)), by = set]
    samp[, measured := id %in% measured$id]

    length_samp <- samp[samp$measured == TRUE, ]
    rm(measured)

    ## Sample ages
    if (age_sampling == "stratified") {
      aged <- length_samp[, list(id = id[sample(.N, ifelse(.N > ages_cap, ages_cap, .N),
                                                replace = FALSE)]),
                          by = c("sim", "year", age_space_group, "length_group")]
    }
    if (age_sampling == "random") {
      aged <- length_samp[, list(id = id[sample(.N, ifelse(.N > ages_cap, ages_cap, .N),
                                                replace = FALSE)]),
                          by = c("set")]
    }
    ## Tag ages sampled
    samp[, aged := id %in% aged$id]
    samp <- samp[, list(set, id, length, age, measured, aged)]
    if (light) samp[, id := NULL]

    # sim$samp <- samp
    setdet <- merge(sets, samp[, .(n = .N, n_measured = sum(measured),
                                   n_aged = sum(aged)), by = "set"],
                    by = "set", all.x = TRUE)
    sim$I_at_length <- I_at_length

    sim$N <- tapply(round(sim$sp_N$N),
                    list(age = sim$sp_N$age, year = sim$sp_N$year),
                    sum, default = 0)

    sim$sp_N$I <- sim$sp_N$N * q(sim$sp_N$age)

    sim$I_at_length_det <- selectivity_integral(sim, l50=l50, l95=l95)
  }

  setdet$n_measured[is.na(setdet$n_measured)] <- 0
  setdet$n_aged[is.na(setdet$n_aged)] <- 0
  setdet$n[is.na(setdet$n)] <- 0


  sim$samp_totals <- setdet[, .(n_sets = .N,
                                n_caught = sum(n),
                                n_measured = sum(n_measured),
                                n_aged = sum(n_aged)
  ), by = .(sim, year)]

  sim$I <- I
  sim$setdet <- setdet
  sim$samp <- samp
  sim$sets <- sets
  rownames(sim$I_at_length) <- as.numeric(rownames(sim$I_at_length))

  return(sim)
}

#' Simulate stratified random surveys using parallel computation
#'
#' This function is a wrapper for \code{\link{sim_survey}} except it allows for
#' many more total iterations to be run than \code{\link{sim_survey}} before running
#' into RAM limitations. Unlike \code{\link{test_surveys}}, this function retains
#' the full details of the survey and it may therefore be more useful for testing
#' alternate approaches to a stratified analysis for obtaining survey indices.
#'
#' @param sim               Simulation from \code{\link{sim_distribution}}
#' @param n_sims            Number of times to simulate a survey over the simulated population.
#'                          Requesting a large number of simulations here may max out your RAM.
#' @param n_loops           Number of times to run the \code{\link{sim_survey}} function. Total
#'                          simulations run will be the product of \code{n_sims} and \code{n_loops}
#'                          arguments. Low numbers of \code{n_sims} and high numbers of \code{n_loops}
#'                          will be easier on RAM, but may be slower.
#' @param cores             Number of cores to use in parallel. More cores should speed up the process.
#' @param quiet             Print message on what to expect for duration?
#' @inheritDotParams sim_survey
#'
#' @details \code{\link{sim_survey}} is hard-wired here to be "light" to minimize object size.
#'
#' @return Returns an object of the same structure as \code{\link{sim_survey}}.
#'
#' @examples
#'
#' \donttest{
#' ## This call runs a total of 25 simulations of the same survey over
#' ## the same population (Note: total number of simulations are low to
#' ## decrease computation time for the example)
#' sim <- sim_abundance(ages = 1:20, years = 1:5) %>%
#'            sim_distribution(grid = make_grid(res = c(10, 10))) %>%
#'            sim_survey_parallel(n_sims = 5, n_loops = 5, cores = 1,
#'                                q = sim_logistic(k = 2, x0 = 3),
#'                                quiet = FALSE)
#' }
#'
#'
#' @export
#'

sim_survey_parallel <- function(sim, n_sims = 1, n_loops = 100,
                                cores = 1, quiet = FALSE, ...) {

  j <- loop <- new_set <- NULL

  start <- Sys.time()
  one_res <- sim_survey(sim, n_sims = n_sims, light = TRUE, ...)
  end <- Sys.time()
  elapsed <- end - start
  max_dur <- end + (elapsed * n_loops) - start

  if (!quiet) {
    message(paste("One run of sim_survey took ~",
                  round(elapsed), attr(elapsed, "units"),
                  "to run. It may take up to",
                  round(max_dur), attr(max_dur, "units"),
                  "to run all simulations."))
  }

  cl <- makeCluster(cores) # use parallel computation
  registerDoParallel(cl)
  loop_res <- foreach(j = seq(n_loops),
                      .packages = "SimSurvey") %dopar% {
                        res <- sim_survey(sim, n_sims = n_sims, light = TRUE, ...)
                        keep <- c("samp_totals", "setdet", "samp")
                        loop_res <- lapply(keep, function(nm) {
                          x <- res[[nm]]
                          x$loop <- j
                          x
                        })
                        names(loop_res) <- keep
                        loop_res
                      }
  stopCluster(cl) # stop parallel process

  ## Combine objects from loop
  samp_totals <- data.table::rbindlist(lapply(loop_res, `[[`, "samp_totals"))
  setdet <- data.table::rbindlist(lapply(loop_res, `[[`, "setdet"))
  samp <- data.table::rbindlist(lapply(loop_res, `[[`, "samp"))

  ## Fix numbering
  samp_totals$new_sim <- samp_totals$sim + (samp_totals$loop * n_sims - n_sims)
  setdet$new_sim <- setdet$sim + (setdet$loop * n_sims - n_sims)
  setdet$new_set <- seq.int(nrow(setdet))
  samp <- merge(samp, setdet[, list(set, loop, new_set)], by = c("set", "loop"),
                sort = FALSE)
  samp_totals$loop <- setdet$loop <- samp$loop <- NULL
  samp_totals$sim <- setdet$sim <- NULL
  setdet$set <- samp$set <- NULL
  setnames(samp_totals, "new_sim", "sim")
  setnames(setdet, "new_sim", "sim")
  setnames(setdet, "new_set", "set")
  setnames(samp, "new_set", "set")

  ## Add new stuff to main object
  sim$I <- one_res$I
  sim$I_at_length <- one_res$I_at_length
  sim$setdet <- setdet
  sim$samp <- samp
  sim$samp_totals <- samp_totals
  sim

}
