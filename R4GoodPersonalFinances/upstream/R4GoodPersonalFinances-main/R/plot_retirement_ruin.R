#' Plotting retirement ruin
#' 
#' @param portfolio_return_mean A numeric. Mean of portfolio returns.
#' @param portfolio_return_sd A numeric. Standard deviation of portfolio returns.
#' @param age A numeric. Current age.
#' @param gompertz_mode A numeric. Gompertz mode.
#' @param gompertz_dispersion A numeric. Gompertz dispersion.
#' @param portfolio_value A numeric. Initial portfolio value.
#' @param monthly_spendings A numeric. Monthly spendings.
#' 
#' @return A [ggplot2::ggplot()] object showing the probability of 
#' retirement ruin for different monthly spending levels. 
#' If a specific 'monthly_spendings' value is provided,
#' it will be highlighted on the plot with annotations.
#' 
#' @examples
#' plot_retirement_ruin(
#'   portfolio_return_mean = 0.034,
#'   portfolio_return_sd   = 0.15,
#'   age                   = 65,
#'   gompertz_mode         = 88,
#'   gompertz_dispersion   = 10,
#'   portfolio_value       = 1000000,
#'   monthly_spendings     = 3000
#' )
#' @export

plot_retirement_ruin <- function(
  portfolio_return_mean,
  portfolio_return_sd,
  age,
  gompertz_mode,
  gompertz_dispersion,
  portfolio_value,
  monthly_spendings = NULL
) {
  
annotate_monthly_spendings <- FALSE
if (length(monthly_spendings) == 1) {
  annotate_monthly_spendings <- TRUE
  monthly_spendings_to_annotate <- monthly_spendings
}
  
if (is.null(monthly_spendings)) {
  monthly_spendings <- round(portfolio_value * 0.04 / 12 / 1000,  0) * 1000
}

if (length(monthly_spendings) == 1) {

  if (monthly_spendings >= 1000) scale <- 1000 else scale <- 100
  
  monthly_spendings <- 
    seq(
      from = round(monthly_spendings * 0.20 / scale,  0) * scale,
      to   = round(monthly_spendings * 2.50 / scale,  0) * scale,
      by   = scale / 2
    )
}

yearly_spendings <- 12 * monthly_spendings
spending_rate    <- yearly_spendings / portfolio_value

retirement_ruin <- 
  calc_retirement_ruin(
  yearly_spendings      = yearly_spendings,
  portfolio_value       = portfolio_value,
  portfolio_return_mean = portfolio_return_mean,
  portfolio_return_sd   = portfolio_return_sd,
  age                   = age,
  gompertz_mode         = gompertz_mode,
  gompertz_dispersion   = gompertz_dispersion
) 
  
value_colour <- "grey60"
  
the_plot <- 
  dplyr::tibble(
    monthly_spendings = monthly_spendings,
    retirement_ruin = retirement_ruin
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = monthly_spendings, 
      y = retirement_ruin,
      color = retirement_ruin
    )
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::scale_x_continuous(
    breaks = scales::breaks_extended(n = 10),
    labels = function(x) {
      paste0(
        format_currency(x / 1000, suffix = "k", accuracy = 0.1
      ),
        "<br><span style='color: grey60;'>(",
        format_currency(x / 1000 * 12, suffix = "k"),
        "/",
        format_percent(x * 12 / portfolio_value),
        ")</span>"
      )
    }
  ) + 
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = scales::percent,
    expand = c(0.005, 0.005)
  ) +
  ggplot2::theme_minimal() +
  ggplot2::guides(color = "none") +
  ggplot2::theme(
    axis.text.x = ggtext::element_markdown(),
    plot.caption = 
      ggtext::element_markdown(
        color = "grey60", 
        size  = 10
      ),
    plot.subtitle = ggtext::element_markdown(color = "grey60")
  ) +
  ggplot2::scale_color_gradientn(
    colors = c(
      PrettyCols::prettycols("Bold")[4], 
      rep(PrettyCols::prettycols("Bold")[5], 10)
    )
  ) + 
  ggplot2::labs(
    title = "Probability of Retirement Ruin",
    x = "Monthly spendings in thousands (yearly spendings / initial year withdrawal rate)",
    y = "Probability of Retirement Ruin",
    caption = glue::glue(paste(
      "*Mean of portfolio returns*: <span style='color: {value_colour};'>**{format_percent(portfolio_return_mean)}**</span>;",
      "*Standard deviation of portfolio returns*: <span style='color: {value_colour};'>**{format_percent(portfolio_return_sd)}**</span>.",
      "<br>",
      "*Gompertz mode*: <span style='color: {value_colour};'>**{round(gompertz_mode, 1)}**</span>;",
      "*Gompertz dispersion*: <span style='color: {value_colour};'>**{round(gompertz_dispersion, 1)}**</span>."
    )),
    subtitle = glue::glue(paste0(
      "*Current age*: <span style='color: {value_colour};'>**{age}**</span>;",
      " *Initial portfolio value*: <span style='color: {value_colour};'>**{format_currency(portfolio_value)}**</span>."
    ))
  ) 
      
  if (annotate_monthly_spendings) {

    retirement_ruin_to_annotate <- 
      calc_retirement_ruin(
        monthly_spendings     = monthly_spendings_to_annotate,
        portfolio_value       = portfolio_value,
        portfolio_return_mean = portfolio_return_mean,
        portfolio_return_sd   = portfolio_return_sd,
        age                   = age,
        gompertz_mode         = gompertz_mode,
        gompertz_dispersion   = gompertz_dispersion
      )

    the_plot <- 
      the_plot +
      ggplot2::geom_vline(
        xintercept = monthly_spendings_to_annotate,
        color      = PrettyCols::prettycols("Bold")[1],
        linetype   = "dashed"
      ) + 
      ggplot2::geom_hline(
        yintercept = retirement_ruin_to_annotate,
        color      = PrettyCols::prettycols("Bold")[1],
        linetype   = "dashed"
      ) +
      ggplot2::annotate(
        geom  = "label",
        x     = monthly_spendings_to_annotate,
        y     = max(retirement_ruin) * 0.95,
        label = paste0(
          format_currency(monthly_spendings_to_annotate, accuracy = 1),
          "\n(", 
          format_currency(
            monthly_spendings_to_annotate * 12 / 1000, 
            suffix = "k"
          ),
          " / ",
          format_percent(monthly_spendings_to_annotate * 12 / portfolio_value), 
          ")"
        ),
        color = PrettyCols::prettycols("Bold")[1]
      ) +
      ggplot2::annotate(
        geom  = "label",
        x     = min(monthly_spendings),
        y     = retirement_ruin_to_annotate,
        label = format_percent(retirement_ruin_to_annotate),
        color = PrettyCols::prettycols("Bold")[1]
      )
  }

  the_plot
}
