#' Plot optimal portfolio allocations
#' 
#' @description
#' The function plots current versus optimal 
#' portfolio allocations for each asset class 
#' and for taxable and tax-advantaged accounts.
#' 
#' @inheritParams simulate_scenario
#' @returns A [ggplot2::ggplot()] object.
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
#' portfolio$accounts$taxadvantaged <- c(0, 20000)
#' portfolio <- 
#'   portfolio |> 
#'   calc_effective_tax_rate(
#'     tax_rate_ltcg = 0.20, 
#'     tax_rate_ordinary_income = 0.40
#'   )
#' 
#' portfolio <- 
#'   calc_optimal_asset_allocation(
#'    household = household,
#'    portfolio = portfolio,
#'    current_date = "2020-07-15"
#'   )
#' 
#' plot_optimal_portfolio(portfolio)
#' @export
plot_optimal_portfolio <- function(
  portfolio
) {

  account <- name <- allocation <- total_allocation <- value <- total_value <- 
    portfolio_label <- NULL

  allocations <- portfolio$allocations

  accounts <- portfolio$accounts
 
  allocations_long <- 
    allocations |> 
    tidyr::pivot_longer(
      cols      = -"name",
      names_to  = "portfolio",
      values_to = "accounts"
    ) |> 
    tidyr::unnest(cols = "accounts") |>
    tidyr::pivot_longer(
      cols      = - c("name", "portfolio"),
      names_to  = "account",
      values_to = "allocation"
    ) |> 
    dplyr::filter(account != "total") |>
    dplyr::mutate(
      portfolio = factor(portfolio, levels = c("optimal", "current"))
    ) |> 
    dplyr::mutate(
      value = allocation * sum(accounts)
    )

  sums <- 
    allocations_long |>
    dplyr::group_by(name, portfolio) |>
    dplyr::summarise(total_allocation = sum(allocation), .groups = "drop") 

  total_values <- 
    allocations_long |>
    dplyr::group_by(portfolio) |>
    dplyr::summarise(total_value = sum(value), .groups = "drop")

  total_values_per_asset <- 
    allocations_long |>
    dplyr::group_by(name, portfolio) |>
    dplyr::summarise(total_value = sum(value), .groups = "drop")

  colors <- PrettyCols::prettycols("Bold")
  color_values <- c(
    "taxable"       = colors[3],
    "taxadvantaged" = colors[4]
  )

  allocations_with_labels <- 
    allocations_long |>
    dplyr::left_join(
      total_values_per_asset |>
        dplyr::mutate(
          portfolio_label = paste0(
            portfolio, "\n", 
            format_currency(total_value, accuracy = 1)
          )
        ) |>
        dplyr::select(name, portfolio, portfolio_label),
      by = c("name", "portfolio")
    ) |>
    dplyr::mutate(
      portfolio_label = factor(portfolio_label, 
                              levels = unique(portfolio_label[order(portfolio)]))
    )

  ggplot2::ggplot(
    allocations_with_labels,
    ggplot2::aes(
      x     = portfolio_label, 
      y     = allocation, 
      fill  = account
    )
  ) +
    ggplot2::scale_fill_manual(
      values = color_values
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::geom_col(
      position = ggplot2::position_stack(), 
      width = 0.9
    ) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ name, ncol = 1, scales = "free_y") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8)) + 
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(allocation < 0.001, "", format_percent(allocation))
      ),
      position      = ggplot2::position_stack(vjust = 0.5),
      size          = 3,
      check_overlap = TRUE
    ) + 
    ggplot2::geom_label(
      data = sums |>
        dplyr::left_join(
          total_values_per_asset |>
            dplyr::mutate(

              portfolio_label =
                paste0(
                  portfolio, "\n",
                  format_currency(total_value, accuracy = 1)
                )
            ) |>
            dplyr::select(name, portfolio, portfolio_label),
          by = c("name", "portfolio")
        ),
      ggplot2::aes(
        x = portfolio_label,
        y = total_allocation,
        label = format_percent(total_allocation)
      ),
      inherit.aes = FALSE,
      hjust = 0.5,
      size = 3
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = glue::glue("Current vs Optimal Portfolio Allocation"),
      x     = "Portfolio",
      y     = "Allocation",
    ) +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      plot.caption = 
        ggtext::element_markdown(
          color = "grey60", 
          size  = 10
        ),
      plot.subtitle = ggtext::element_markdown(color = "grey60"),
      axis.text.y = ggtext::element_markdown(
        color = c(
          "optimal" = "black", 
          "current" = "grey60"
        )
      )
    )
}
