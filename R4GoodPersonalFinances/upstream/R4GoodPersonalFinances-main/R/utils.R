#' Printing currency values or percentages
#' 
#' Wrapper functions for printing nicely formatted values.
#' 
#' @seealso [scales::dollar()] 
#' 
#' @inheritParams scales::dollar
#' @inheritParams scales::percent
#' @param min_length A numeric. Minimum number of characters 
#' of the string with the formatted value.
#' @rdname format
#' 
#' @return A character. Formatted value.
#' 
#' @examples
#' format_currency(2345678, suffix = " PLN")
#' @export

format_currency <- function(
  x, 
  prefix   = "",
  suffix   = "",
  big.mark = ",",
  accuracy = NULL,
  min_length = NULL, 
  ...
) {
  
  if (is.list(x)) {
    return(
      purrr::map(
        x, 
        format_currency, 
        suffix   = suffix, 
        big.mark = big.mark,
        accuracy = accuracy,
        prefix   = prefix,
        min_length = min_length,
        ...
      )
    )
  }

  if (!is.numeric(x)) {
    return(x)
  }

  formatted_value <- 
    scales::dollar(
      x        = x, 
      prefix   = prefix, 
      suffix   = suffix,
      big.mark = big.mark,
      accuracy = accuracy,
      ...
    )

  if (is.null(min_length)) {
    return(formatted_value)
  }

  padding <- min_length - nchar(formatted_value)
  
  formatted_value <- 
    stringr::str_pad(
      formatted_value,
      width = nchar(formatted_value) + padding,
      side = "left"
    )
  
  formatted_value
}
    
#' @seealso [scales::percent()]
#' 
#' @inheritParams scales::percent
#' @inheritParams scales::dollar
#' @rdname format
#' 
#' @return A character. Formatted value.
#' 
#' @examples
#' format_percent(0.52366)
#' @export

format_percent <- function(x, 
                          accuracy = 0.1,
                          ...) {
  
  if (is.list(x)) {
    return(
      purrr::map(x, format_percent, accuracy = accuracy, ...)
    )
  }

  if (!is.numeric(x)) {
    return(x)
  }
  
  percents <- scales::percent(
    x = x,
    accuracy = accuracy,
    ...
  )

  names(percents) <- names(x)
  percents
}

#' Get current date
#' 
#' If `R4GPF.current_date` option is not set, the current system date is used.
#' 
#' @return A date.
#' @examples
#' get_current_date()
#' # Setting custom date using `R4GPF.current_date` option
#' options(R4GPF.current_date = as.Date("2023-01-01"))
#' get_current_date()
#' options(R4GPF.current_date = NULL) # Reset default date#' Working with cache
#' 
#' get_current_date()
#' @export
get_current_date <- function() {
  current_date <- getOption("R4GPF.current_date", default = Sys.Date())
  lubridate::as_date(current_date)
}

normalize <- function(x, min_val = 0, max_val = 1) {
  
  (x - base::min(x, na.rm = TRUE)) /
    (base::max(x, na.rm = TRUE) - base::min(x, na.rm = TRUE)) *
    (max_val - min_val) + min_val
}


generate_test_asset_returns <- function(n = 3) {
  expected_return <- capital_gains <- NULL

  if (n == 3) {

    test_asset_returns <- 
      dplyr::tribble(
        ~asset_class,          ~expected_return, ~standard_deviation,
        "DomesticStocks",      0.0472,           0.1588, 
        "InternationalStocks", 0.0504,           0.1718,
        "Bonds",               0.0275,           0.0562
      )

    test_asset_correlations <- dplyr::tribble(
      ~DomesticStocks, ~InternationalStocks, ~Bonds,
      1.00,            0.87,                 0.21,
      0.87,            1.00,                 0.37,
      0.21,            0.37,                 1.00
    )
    
  } else if (n == 2) {
    
    portfolio <- 
      dplyr::tribble(
        ~name,                     ~expected_return, ~standard_deviation, 
        "GlobalStockIndexFund",    0.0449,           0.15,                
        "InflationProtectedBonds", 0.02,             0,                  
      ) |> 
      dplyr::mutate(
        accounts = dplyr::tribble(
          ~taxable, ~taxadvantaged,
          200000,   50000,
          100000,   25000,
        ),
        weights = dplyr::tribble(
          ~human_capital, ~liabilities, 
          0.5,            0.5,          
          0.5,            0.5,          
        )
      )
      portfolio <- 
        portfolio |> 
        dplyr::mutate(
          correlations = {
            matrix <- diag(1, NROW(portfolio), NROW(portfolio)) 
            colnames(matrix) <- portfolio$name
            rownames(matrix) <- portfolio$name
            matrix
          },
          aftertax = dplyr::tibble(
            effective_tax_rate = rep(0.19, NROW(portfolio))
          )
        )
    
    test_asset_returns      <- portfolio
    test_asset_correlations <- portfolio$correlations

  } else if (n == 9) {

    test_asset_returns <- 
      dplyr::tribble(
        ~name,                  ~expected_return, ~standard_deviation, 
        "USLargeCapStocks",     0.0468,           0.1542,
        "USMidSmallCapStocks",  0.0501,           0.1795, 
        "GlobalDMxUSStocks",    0.0505,           0.1671, 
        "EmergingMarketStocks", 0.0540,           0.2142,
        "USBonds",              0.0269,           0.0379,
        "InflationLinkedBonds", 0.0288,           0.0581,
        "MuniBonds",            0.0190,           0.03138274,
        "GlobalBondExUS",       0.0329,           0.0833,
        "Cash",                 0.0250,           0.0055
      )
    
      test_asset_returns <- 
        test_asset_returns |> 
        dplyr::mutate(
          pretax = dplyr::tibble(
            capital_gains    = c(0.0349, 0.0387, 0.0336, 0.0388, rep(0, 5)),
            income           = expected_return - capital_gains,
            turnover         = c(0.3300, 0.3652, 0.1800, 0.3300, rep(1, 5)),
            cost_basis       = c(0.9364, 0.9393, 0.8750, 0.9301, rep(1, 5)),
            income_qualified = c(0.9762, 0.9032, 0.7998, 0.7387, rep(0, 5)),
            capital_gains_long_term = 
              c(0.9502, 0.9032, 0.8951, 0.9023, rep(0, 5))
          )
        )|> 
        dplyr::mutate(
          accounts = dplyr::tibble(
            taxable       = rep(1000, n), 
            taxadvantaged = rep(1000, n)
          ),
          weights = dplyr::tibble(
            human_capital = 1 / n, 
            liabilities   = 1 / n, 
          )
        ) 
    
      test_asset_correlations <- 
        diag(rep(1, length(test_asset_returns$expected_return)))
    
      test_asset_returns <- 
        test_asset_returns |> 
        dplyr::mutate(
          correlations = test_asset_correlations,
          effective_tax_rates = 0.19
        ) 
        
  }

  if (!is.null(test_asset_correlations)) {
    test_asset_correlations <- as.matrix(test_asset_correlations)
    rownames(test_asset_correlations) <- colnames(test_asset_correlations)
  }

  list(
    returns      = test_asset_returns,
    correlations = test_asset_correlations
  )

}

#' Working with cache
#' 
#' Get information about the cache
#' 
#' @rdname cache
#' @returns Invisibly returns the path to the cache directory 
#' or a list containing:
#' \item{path}{The path to the cache directory.} 
#' \item{files}{The number of files in the cache.}
#' @examplesIf interactive()
#' get_cache_info()
#' @export
get_cache_info <- function() {

  list(
    path  = .pkg_env$cache_directory,
    files = .pkg_env$cache$size()
  )
}

#' Working with cache
#' 
#' Reset the cache
#' 
#' @rdname cache
#' @examplesIf interactive()
#' reset_cache()
#' @export
reset_cache <- function() {

  .pkg_env$cache$reset()
  invisible(.pkg_env$cache_directory)
}

#' Working with cache
#' 
#' Set the cache directory
#' 
#' @param path The path to the cache directory.
#' Defaults to the '.cache' folder in the current working directory.
#' 
#' @rdname cache
#' @examplesIf interactive()
#' set_cache()
#' @export
set_cache <- function(
  path = file.path(getwd(), ".cache")
) {

  .pkg_env$cache_directory <- path

  .pkg_env$cache <- 
    cachem::cache_disk(
      dir        = .pkg_env$cache_directory,
      prune_rate = 100000
    )
  
  .pkg_env$memoised$simulate_single_scenario <- 
    memoise::memoise(
      f         = simulate_single_scenario,
      cache     = .pkg_env$cache,
      omit_args = c("debug", "verbose", "auto_parallel")
    )

  .pkg_env$memoised$simulate_scenario <- 
    memoise::memoise(
      f         = simulate_scenario,
      cache     = .pkg_env$cache,
      omit_args = c("debug", "verbose", "auto_parallel")
    )

  invisible(.pkg_env$cache_directory)
}

paste_labels <- function(breaks, scenario) {
  
  members_age <-
    names(scenario$members) |>
      purrr::map(function(member_name) {
    
    scenario[scenario$index %in% breaks,
    ]$members[[member_name]]$age |>
      unique()
  }) |>
    purrr::set_names(names(scenario$members))
  
  max_length <- max(purrr::map_int(members_age, length))
  
  ages <- purrr::map_chr(seq_len(max_length), function(i) {
    member_strings <- purrr::imap_chr(members_age, function(age_vec, name) {
      if (length(age_vec) >= i && !is.na(age_vec[i])) {
        paste0(abbreviate(name, minlength = 999), " (", age_vec[i], ")")
      } else {
        ""
      }
    })
    paste0(member_strings, collapse = "\n")
  })
  
  paste0(breaks, "\n", ages)
}

paste_year_index_axis_label <- function() {
  "Year index / household member (age)"
}

paste_scenario_id <- function(scenario) {
  
  if ("scenario_id" %in% names(scenario)) {
    glue::glue(
      "Scenario: <strong style='color: grey50;'>{
      unique(scenario$scenario_id)
      }</strong>. "
    )
  } else {
    ""
  }
}

  format_colored_names <- function(data, type_colors) {
    original_names <- as.character(data$type)
    colored_names <- glue::glue(
      "<span style='color: {type_colors[original_names]};'>**{original_names}**</span>"
    )
    return(colored_names)
  }

generate_random_seed_vector <- function(n) {

  stats::runif(n = n, min = 0, max = 1e8) |> 
    as.integer()
}

generate_random_seeds <- function(monte_carlo_samples, seeds) {

  if (!is.null(monte_carlo_samples)) {

    if (is.null(seeds)) {

      cli::cli_alert_info("Generating random seeds for Monte Carlo samples")
      seeds <- generate_random_seed_vector(n = monte_carlo_samples)

    } else if (length(seeds) == monte_carlo_samples) {

      cli::cli_alert_info("Using provided random seeds for Monte Carlo samples")
      
    } else if (length(seeds) == 1) {

      cli::cli_alert_info("Setting random seed to {.field {seeds}}")
      set.seed(seeds)

      cli::cli_alert_info("Generating random seeds for Monte Carlo samples")
      seeds <- generate_random_seed_vector(n = monte_carlo_samples)
    }
  }
  seeds
}
