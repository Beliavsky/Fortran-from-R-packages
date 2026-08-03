#' Plot expected allocation over household life cycle
#' 
#' If multiple Monte Carlo samples are provided in the `scenario` argument,
#' the normalized median of expected allocation is plotted.
#' 
#' @param scenario A `tibble` with nested columns -
#' the result of [simulate_scenario()]. Data for a single scenario.
#' @param accounts A character. 
#' Plot allocation for specified types of accounts.
#' @returns A [ggplot2::ggplot()] object. 
#' @examplesIf interactive()
#' older_member <- HouseholdMember$new(
#'   name       = "older",  
#'   birth_date = "2000-02-15",
#'   mode       = 80,
#'   dispersion = 10
#' )  
#' household <- Household$new()
#' household$add_member(older_member)  
#' 
#' household$expected_income <- list(
#'   "income" = c(
#'     "members$older$age <= 65 ~ 10000 * 12"
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
#' portfolio$weights$human_capital <- c(0.2, 0.8)
#' portfolio$weights$liabilities <- c(0.1, 0.9)
#' portfolio <- 
#'   portfolio |> 
#'   calc_effective_tax_rate(
#'     tax_rate_ltcg            = 0.20, 
#'     tax_rate_ordinary_income = 0.40
#'   )
#' 
#' scenario <- 
#'   simulate_scenario(
#'    household    = household,
#'    portfolio    = portfolio,
#'    current_date = "2020-07-15"
#'   )
#' 
#' plot_expected_allocation(scenario)
#' @export
plot_expected_allocation <- function(
  scenario,
  accounts = c("all", "taxable", "taxadvantaged")
) {

  index <- asset <- allocation <- portfolio <- NULL

  accounts     <- rlang::arg_match(accounts)
  account_type <- ifelse(accounts == "all", "total", accounts)
  
  data_to_plot <- 
    scenario |> 
    dplyr::transmute(
      sample, 
      index, 
      portfolio = portfolio$allocation
    ) |> 
    tidyr::unnest_longer(portfolio) |>
    dplyr::transmute(
      sample,
      index, 
      asset = portfolio$asset, 
      taxable = portfolio$taxable,
      taxadvantaged = portfolio$taxadvantaged,
      total = portfolio$total,
    ) |> 
    dplyr::select(
      sample,
      index, 
      asset,
      allocation = dplyr::all_of(account_type)
    ) |> 
    dplyr::group_by(index, asset) |>
    dplyr::summarise(
      allocation = stats::median(allocation)
    ) |> 
    dplyr::group_by(index) |> 
    dplyr::mutate(allocation = allocation / sum(allocation)) 
  
  colors <- 
    grDevices::colorRampPalette(
      rev(PrettyCols::prettycols("Bold"))
    )(
      length(unique(data_to_plot$asset))
    )
  
  data_to_plot |> 
    ggplot2::ggplot(
    ggplot2::aes(
      x    = index, 
      y    = allocation, 
      fill = factor(asset, levels = unique(asset)))
      ) +
    ggplot2::geom_area() +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_x_continuous(
      breaks = function(limits) {
        max_index <- max(scenario$index)
        interval <- ifelse(max_index > 50, 10, ifelse(max_index > 20, 5, 2))
        seq(0, max_index, by = interval)
      },
       labels = function(breaks) paste_labels(breaks, scenario = scenario)
     ) +
    ggplot2::scale_y_continuous(labels = scales::percent) + 
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title     = ggplot2::element_blank(),
      plot.caption = 
        ggtext::element_markdown(
          color = "grey60", 
          size  = 10
        ),
      plot.subtitle = ggtext::element_markdown(color = "grey60")
    ) +
    ggplot2::labs(
      x = paste_year_index_axis_label(),
      y = "Allocation",
      title = glue::glue(paste0(
        "Optimal Asset Allocation Over Time",
        ifelse(
          account_type != "total",
          " (in {account_type} accounts)",
          ""
        )
      )), 
      subtitle = glue::glue(paste0(
        paste_scenario_id(scenario),
        "Based on expected returns",
        ifelse(
          max(scenario$sample) == 0, 
          ".",
          " and <strong>{max(scenario$sample)}</strong> Monte Carlo samples. "
        )
        
      ))
    )
}
