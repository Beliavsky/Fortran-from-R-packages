#' Simulate multiple scenarios of household lifetime finances
#' 
#' @inheritParams simulate_scenario
#' @param scenarios_parameters A `tibble` with column
#' `scenario_id` and nested column `events`.
#' Each scenario has defined one or more events in the tibbles 
#' that are stored in as a list in the `events` column. 
#' 
#' @returns A `tibble` with nested columns.
#' @examplesIf interactive()
#' older_member <- HouseholdMember$new(
#'   name       = "older",  
#'   birth_date = "1980-02-15",
#'   mode       = 80,
#'   dispersion = 10
#' )  
#' household <- Household$new()
#' household$add_member(older_member)  
#' 
#' household$expected_income <- list(
#'   "income" = c(
#'     "members$older$age <= 65 ~ 7000 * 12"
#'   )
#' )
#' household$expected_spending <- list(
#'   "spending" = c(
#'     "TRUE ~ 5000 * 12"
#'   )
#' )
#' 
#' portfolio <- create_portfolio_template() 
#' portfolio$accounts$taxable <- c(10000, 30000)
#' portfolio <- 
#'   portfolio |> 
#'   calc_effective_tax_rate(
#'     tax_rate_ltcg = 0.20, 
#'     tax_rate_ordinary_income = 0.40
#'   )
#' start_ages <- c(60, 65, 70)
#' scenarios_parameters <- 
#'   tibble::tibble(
#'     member    = "older",
#'     event      = "retirement",
#'     start_age = start_ages,
#'     years     = Inf,
#'     end_age   = Inf
#'    ) |> 
#'   dplyr::mutate(scenario_id = start_age) |> 
#'   tidyr::nest(events = -scenario_id)
#' 
#' scenarios_parameters
#' 
#' scenarios <- 
#'   simulate_scenarios(
#'     scenarios_parameters = scenarios_parameters,
#'     household            = household,
#'     portfolio            = portfolio,
#'     current_date         = "2020-07-15"
#'   )
#' scenarios$scenario_id |> unique()
#' @export
simulate_scenarios <- function(
  scenarios_parameters,
  household,
  portfolio,
  current_date        = get_current_date(),
  monte_carlo_samples = NULL,
  seeds               = NULL,
  auto_parallel       = FALSE,
  use_cache           = FALSE,
  debug               = FALSE,
  ...
) {

  scenario_id <- index <- NULL

  current_date <- lubridate::as_date(current_date)
    
  scenarios_ids <- unique(scenarios_parameters$scenario_id)

  cli::cli_h1(glue::glue(
    "Simulating {length(scenarios_ids)} scenarios"
  ))

  cli::cli_alert_info("Cache directory: {.file {(.pkg_env$cache_directory)}}")

  seeds <- 
    generate_random_seeds(
      monte_carlo_samples = monte_carlo_samples, 
      seeds = seeds
    )

  if (auto_parallel) {

    if (!is.null(monte_carlo_samples)) {

      n_workers <- future::availableCores()

      old_plan <- future::plan(future::multisession, workers = n_workers)

      cli::cli_alert_info(
        "Auto-parallelization enabled with {.field {n_workers}} workers."
      )
      
      on.exit(future::plan(old_plan))
    } 
  }

  if (use_cache) {

    cli::cli_alert_info("Caching is enabled!")

    memoised_functions       <- .pkg_env$memoised
    simulate_scenario <- memoised_functions$simulate_scenario
    
  } else {

    cli::cli_alert_warning(cli::col_yellow("Caching is NOT enabled."))

  }

  scenarios_ids |> 
    purrr::map(function(scenario_id) {

      cli::cli_h2(paste0(
        "Scenario ", 
        "{which(scenario_id == scenarios_ids)}/{length(scenarios_ids)} ",
        "({format_percent(which(scenario_id == scenarios_ids) / 
          length(scenarios_ids), 1)})"
      ))
      
      scenario_params <- 
        scenarios_parameters |> 
        dplyr::filter(scenario_id == !!scenario_id)
        
      scenario_events  <- scenario_params$events[[1]]
      household_cloned <- unserialize(serialize(household, NULL))

      for (i in NROW(scenario_events)) {
        household_cloned$get_members()[[scenario_events[i, ]$member]]$set_event(
          event     = scenario_events[i, ]$event,  
          start_age = scenario_events[i, ]$start_age, 
          end_age   = scenario_events[i, ]$end_age,
          years     = scenario_events[i, ]$years)
      }
        
      simulate_scenario(
        household           = household_cloned,
        portfolio           = portfolio,
        current_date        = current_date,
        scenario_id         = scenario_id,
        monte_carlo_samples = monte_carlo_samples,
        seeds               = seeds,
        use_cache           = use_cache,
        debug               = debug,
        ...
      ) 

    }) |> 
    purrr::list_rbind() |> 
    dplyr::select(
      scenario_id, 
      index, 
      dplyr::everything()
    )
}
