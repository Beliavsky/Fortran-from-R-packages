#' Plotting risk adjusted returns
#' 
#' Plots the risk adjusted returns for portfolios
#' of various allocations to the risky asset.
#' 
#' @seealso 
#' * [How to Determine Our Optimal Asset Allocation?](https://www.r4good.academy/en/blog/optimal-asset-allocation/index.en.html#how-much-risk-is-enough)
#' * Haghani V., White J. (2023) "The Missing Billionaires: A Guide to Better Financial Decisions." ISBN:978-1-119-74791-8.
#' 
#' @inheritParams calc_risk_adjusted_return
#' @inheritParams calc_optimal_risky_asset_allocation
#' 
#' @param current_risky_asset_allocation A numeric. 
#' The current allocation to the risky asset.
#' For comparison with the optimal allocation.
#' 
#' @returns A [ggplot2::ggplot()] object.
#' 
#' @examples
#' plot_risk_adjusted_returns(
#'   safe_asset_return              = 0.02,
#'   risky_asset_return_mean        = 0.04,
#'   risky_asset_return_sd          = 0.15,
#'   risk_aversion                  = 2,
#'   current_risky_asset_allocation = 0.8
#' )
#' @export
#' 
plot_risk_adjusted_returns <- function(
  safe_asset_return,
  risky_asset_return_mean,
  risky_asset_return_sd,
  risk_aversion                  = 2,
  current_risky_asset_allocation = NULL
) {

  risky_asset_allocations <- seq(from = 0.00, to = 1.00, by = 0.001)

  risk_adjusted_returns <-
    calc_risk_adjusted_return(
      safe_asset_return       = safe_asset_return,
      risky_asset_return_mean = risky_asset_return_mean, 
      risky_asset_return_sd   = risky_asset_return_sd,
      risky_asset_allocation  = risky_asset_allocations,
      risk_aversion           = risk_aversion
    )
  
  optimal_risky_asset_allocation <-
    calc_optimal_risky_asset_allocation(
      safe_asset_return       = safe_asset_return,
      risky_asset_return_mean = risky_asset_return_mean,
      risky_asset_return_sd   = risky_asset_return_sd,
      risk_aversion           = risk_aversion
    )
    
  optimal_risk_adjusted_return <-
    calc_risk_adjusted_return(
      safe_asset_return       = safe_asset_return,
      risky_asset_return_mean = risky_asset_return_mean,
      risky_asset_return_sd   = risky_asset_return_sd,
      risky_asset_allocation  = optimal_risky_asset_allocation,
      risk_aversion           = risk_aversion
    )

  colours <-  PrettyCols::prettycols("Bold")
  current_risky_asset_allocation_colour <- colours[5]
  optimal_allocation_colour             <- colours[4]
  risky_asset_allocation_colour         <- colours[1]
  arrow_colour                          <- colours[2]
  safe_asset_return_color               <- "grey50"

  y_limit_min <- min(risk_adjusted_returns)
  y_limit_min <- ifelse(y_limit_min > 0, 0, y_limit_min)

  add_labels <- function(to_plot) {

    to_plot +
      ggplot2::geom_label(
        ggplot2::aes(
          x = 0,
          y = safe_asset_return,
          label = format_percent(safe_asset_return)
        ),
        color = safe_asset_return_color
      ) +
      ggplot2::geom_label(
        ggplot2::aes(
          x = optimal_risky_asset_allocation,
          y = y_limit_min,
          label = format_percent(optimal_risky_asset_allocation)
        ),
        color = optimal_allocation_colour
      ) +
      ggplot2::geom_label(
        ggplot2::aes(
          x = 0,
          y = max(risk_adjusted_returns),
          label = format_percent(max(risk_adjusted_returns))
        ),
        color = optimal_allocation_colour
      ) 
  }

  value_colour <- "grey40"

  to_plot <-
    ggplot2::ggplot() +
    ggplot2::geom_vline(
      xintercept = optimal_risky_asset_allocation,
      color      = optimal_allocation_colour,
      linetype   = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = max(risk_adjusted_returns),
      color      = optimal_allocation_colour,
      linetype   = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = safe_asset_return,
      color      = safe_asset_return_color,
      linetype   = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color      = safe_asset_return_color,
      linetype   = "dotted"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        x = risky_asset_allocations,
        y = risk_adjusted_returns
      ),
      color = risky_asset_allocation_colour,
      linewidth = 1
    ) +
    ggplot2::labs(
      title    = "Optimal Risky Asset Allocation",
      subtitle = glue::glue(paste0(
        "<span style='color: {optimal_allocation_colour};'>**Optimal**</span> ",
        ifelse(
          is.null(current_risky_asset_allocation), 
          "", 
          "vs <span style='color: {current_risky_asset_allocation_colour};'>**current**</span> "
        ),
        "risky asset allocations and corresponding ",
        "<span style='color: {risky_asset_allocation_colour};'>**risk adjusted returns**</span>"
      )),
      x = "Risky Asset Allocation",
      y = "Risk Adjusted Returns",
      caption = glue::glue(paste(
        "*Safe asset return*: <span style='color: {value_colour};'>**{format_percent(safe_asset_return)}**</span>;",
        "*Risky asset return*: <span style='color: {value_colour};'>**{format_percent(risky_asset_return_mean)}**</span>;",
        "*Risky asset sd*: <span style='color: {value_colour};'>**{format_percent(risky_asset_return_sd)}**</span>;",
        "*Risk aversion*: <span style='color: {value_colour};'>**{risk_aversion}**</span>."
      ))
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 1, by = 0.1),
      labels = scales::percent
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(
        y_limit_min, 
        max(risk_adjusted_returns), 
        by = 0.005
      ),
      limits = c(y_limit_min, NA),
      labels = scales::percent
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.caption = 
        ggtext::element_markdown(
          color = "grey60", 
          size  = 10
        ),
      plot.subtitle = ggtext::element_markdown(color = "grey60")
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = optimal_risky_asset_allocation,
        y = optimal_risk_adjusted_return
      ),
      color = optimal_allocation_colour,
      size  = 4)
  
  if (is.null(current_risky_asset_allocation)) {

    to_plot <- add_labels(to_plot)
    return(to_plot)
  }
  
  current_risk_adjusted_return <-
    calc_risk_adjusted_return(
      safe_asset_return         = safe_asset_return,
      risky_asset_return_mean   = risky_asset_return_mean,
      risky_asset_return_sd     = risky_asset_return_sd,
      risky_asset_allocation    = current_risky_asset_allocation,
      risk_aversion             = risk_aversion
    )
  
  arrow_length <- 9/10

  to_plot <- 
    to_plot +
      ggplot2::geom_hline(
        yintercept = current_risk_adjusted_return,
        color      = current_risky_asset_allocation_colour
      ) +
      ggplot2::geom_vline(
        xintercept = current_risky_asset_allocation,
        color      = current_risky_asset_allocation_colour
      ) +
      ggplot2::geom_label(
        ggplot2::aes(
          x     = 1,
          y     = current_risk_adjusted_return,
          label = format_percent(current_risk_adjusted_return)
        ),
        color = current_risky_asset_allocation_colour
      ) +
      ggplot2::geom_label(
        ggplot2::aes(
          x     = current_risky_asset_allocation,
          y     = y_limit_min + 0.002,
          label = format_percent(current_risky_asset_allocation)
        ),
        color = current_risky_asset_allocation_colour
      ) +
      ggplot2::geom_segment(
        ggplot2::aes(
          x    = current_risky_asset_allocation,
          y    = current_risk_adjusted_return,
          xend = 
            current_risky_asset_allocation - (
              current_risky_asset_allocation - optimal_risky_asset_allocation
            ) * arrow_length,
          yend = 
            current_risk_adjusted_return + (
              optimal_risk_adjusted_return - current_risk_adjusted_return
            ) * arrow_length
            ),
          arrow = ggplot2::arrow(type = "closed",
          length = ggplot2::unit(0.12, "inches")),
          color = arrow_colour
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          x = current_risky_asset_allocation,
          y = current_risk_adjusted_return
        ),
        color = current_risky_asset_allocation_colour,
        size = 4
      )
      
  add_labels(to_plot)
}
