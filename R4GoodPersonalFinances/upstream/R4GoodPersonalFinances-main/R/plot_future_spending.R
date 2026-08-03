#' Plot future spending structure over household life cycle
#' 
#' @description 
#' Plot future spending structure over household life cycle,
#' including discretionary and non-discretionary spending.
#' You can also plot discretionary and non-discretionary spending separately,
#' to see structure of non-discretionary spending and
#' possible levels of discretionary spending over time 
#' based on Monte Carlo simulations.
#' 
#' 
#' @inheritParams plot_future_income
#' @param type A character. Type of spending to plot:
#' discretionary, non-discretionary, or both (default).
#' @param discretionary_spending_position A character.
#' Position of discretionary spending in plot. 
#' Bottom is the default.
#' @returns A [ggplot2::ggplot()] object
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
#'     "members$older$age <= 65 ~ 9000 * 12"
#'   )
#' )
#' household$expected_spending <- list(
#'   "spending" = c(
#'     "members$older$age <= 65 ~ 5000 * 12",
#'     "TRUE ~ 4000 * 12"
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
#' 
#' scenario <- 
#'   simulate_scenario(
#'    household = household,
#'    portfolio = portfolio,
#'    # monte_carlo_samples = 100,
#'    current_date = "2020-07-15"
#'   )
#' 
#' plot_future_spending(scenario, "monthly")
#' plot_future_spending(
#'   scenario, 
#'   "monthly", 
#'   discretionary_spending_position = "top"
#' )
#' plot_future_spending(scenario, "monthly", "non-discretionary")
#' # If Monte Carlo samples are present: 
#' # plot_future_spending(scenario, "monthly", "discretionary")
#' @export
plot_future_spending <- function(
  scenario,
  period                          = c("yearly", "monthly"),
  type                            = c("both", "discretionary", "non-discretionary"),
  discretionary_spending_position = c("bottom", "top"),
  y_limits                        = c(NA, NA)
) {
  
  index <- discretionary_spending <- nondiscretionary_spending <- NULL
  
  period <- rlang::arg_match(period)
  type   <- rlang::arg_match(type)

  if (type == "both") {
    type <- c("discretionary", "non-discretionary")
  }

  discretionary_spending_position <- 
    rlang::arg_match(discretionary_spending_position)

  if (length(type) == 2) {
    return(
      plot_expected_spending(
        scenario, 
        period = period,
        discretionary_spending_position = discretionary_spending_position,
        y_limits = y_limits
      )
    )
  } 

  if (type == "discretionary") {
    return(
      plot_simulated_spending(
        scenario, 
        period   = period,
        y_limits = y_limits
      )
    )
  } 

  if (type == "non-discretionary") {

    return(
      plot_structure(
        scenario, 
        structure_of = "spending",
        period = period,
        y_limits = y_limits
      )
    )
  }
}

plot_simulated_spending <- function(
  scenario,
  period   = c("yearly", "monthly"),
  y_limits = c(NA, NA)
) {

  index <- discretionary_spending <- min_quantiles <- max_quantiles <- 
    quantile_group <- NULL
  
  min_alpha     <- 0.15
  period        <- rlang::arg_match(period)
  period_factor <- if (period == "yearly") 1 else 12

  
  if (length(unique(scenario$sample)) <= 1) {
    cli::cli_abort(
      call = NULL,
      "Plotting Monte Carlo samples requires more than one sample."
    )
  }

  quantile_min <- 
    scenario |>
    dplyr::filter(sample != 0) |>
    dplyr::group_by(index) |>
    dplyr::summarize(
      min_quantiles = list(
        stats::quantile(
          discretionary_spending, 
          probs = seq(0, 0.9, 0.1)
        ) |>
          stats::setNames(1:10)
      )
    ) |>
    tidyr::unnest_longer(
      min_quantiles,
      indices_to = "quantile_group",
      values_to  = "min"
    ) |> 
    dplyr::mutate(min = min / period_factor)
    
  quantile_max <- 
    scenario |>
    dplyr::filter(sample != 0) |>
    dplyr::group_by(index) |>
    dplyr::summarize(
      max_quantiles = list(
        stats::quantile(
          discretionary_spending, 
          probs = seq(0.1, 1, 0.1)
        ) |>
          stats::setNames(1:10)
      )
    ) |>
    tidyr::unnest_longer(
      max_quantiles,
      indices_to = "quantile_group",
      values_to  = "max"
    ) |> 
    dplyr::mutate(max = max / period_factor)

  quantile_data <- 
    dplyr::left_join(
      quantile_min, 
      quantile_max,
      by = c("index", "quantile_group")
    )  |> dplyr::filter(
      !quantile_group %in% c(1, 2, 9, 10)
    )   

  y_max <- max(abs(quantile_data$min), abs(quantile_data$max)) 
  if (y_max > 10000 * 3) {
    y_breaks_factor <- 10000
  } else if (y_max > 1000 * 3) {
    y_breaks_factor <- 1000
  } else {  
    y_breaks_factor <- 100
  }

  y_breaks <- 
    seq(
      round(
        min(quantile_data$min, y_limits[1], na.rm = TRUE) / y_breaks_factor
      ) * y_breaks_factor, 
      round(
        max(quantile_data$max, y_limits[2], na.rm = TRUE) / y_breaks_factor
      ) * y_breaks_factor, 
      by = y_breaks_factor
    )

  group_names <- sort(unique(quantile_data$quantile_group))
  n_groups    <- length(group_names)
  colors      <- PrettyCols::prettycols("Teals")

  plot <- 
    quantile_data |>
    ggplot2::ggplot(
      ggplot2::aes(x = index)
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin  = min, 
        ymax  = max, 
        group = quantile_group, 
        alpha = quantile_group
      ),
      fill = colors[3],
      show.legend = FALSE
    ) +
    ggplot2::scale_alpha_manual(
      values = c(
        seq(min_alpha, 1, length.out = n_groups/2) |> 
          stats::setNames(group_names[1:(n_groups/2)]),
        seq(1, min_alpha, length.out = n_groups/2) |> 
          stats::setNames(group_names[(n_groups/2 + 1):n_groups])
      )
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor.y = ggplot2::element_blank(),
      plot.caption       = ggtext::element_markdown(color = "grey60"),
      plot.subtitle      = ggtext::element_markdown(color = "grey60")
    ) +
    ggplot2::geom_line(
      data = dplyr::filter(scenario, sample == 0),
      ggplot2::aes(
        y = discretionary_spending / period_factor
      ),
      color     = PrettyCols::prettycols("Bold")[3],
      linetype  = "dashed",
      linewidth = 1
    ) +
    ggplot2::geom_line(
      data = 
        scenario |>
          dplyr::filter(sample != 0) |> 
          dplyr::group_by(index) |>
          dplyr::summarize(
            discretionary_spending = stats::median(discretionary_spending) 
          ),
      ggplot2::aes(
        y = discretionary_spending / period_factor
      ),
      color     = colors[1],
      linewidth = 1
    ) + 
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dotted"
    ) +
    ggplot2::labs(
      title = "Future Simulated Discretionary Spending",
      subtitle = glue::glue(paste0(
        paste_scenario_id(scenario),
        "Based on <strong>{max(scenario$sample)}</strong> Monte Carlo samples. ",
        "Median spending at year 0 is <strong>{format_currency(scenario$discretionary_spending[1] / period_factor)}</strong>."
      )),
      caption = "Yellow dashed line shows discretionary spending based on portfolio expected returns.<br>Solid teal line shows median of discretionary spending in Monte Carlo samples.<br>Teal bands show middle six decile groups of spending without top 2 and bottom 2 deciles.",
      x = "Year Index",
      y = glue::glue("Amount ({period})"),
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, max(scenario$index), by = 10),
      labels = function(breaks) paste_labels(breaks, scenario = scenario)
    ) +
    ggplot2::scale_y_continuous(
      labels = format_currency,
      breaks = y_breaks,
      # expand = c(0, NA)
    ) + 
    ggplot2::coord_cartesian(ylim = c(y_limits[1], y_limits[2])) 

  plot
}
