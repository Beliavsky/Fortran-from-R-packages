plot_expected_spending <- function(
  scenario, 
  period                          = c("yearly", "monthly"),
  discretionary_spending_position = c("bottom", "top"),
  y_limits                        = c(NA, NA)
) {
  
  index <- discretionary_spending <- nondiscretionary_spending <- spending <- 
    NULL
  
  stopifnot(
    length(unique(scenario$scenario_id)) == 1
  )

  type          <- c("discretionary", "non-discretionary")
  period        <- rlang::arg_match(period)
  period_factor <- if (period == "yearly") 1 else 12

  discretionary_spending_position <- 
    rlang::arg_match(discretionary_spending_position)
  
  if (discretionary_spending_position == "bottom") {
    type_levels <- c(
      "non-discretionary",
      "discretionary"
    )
  } else {
    type_levels <- c(
      "discretionary",
      "non-discretionary"
    )
  }

  colors <- PrettyCols::prettycols("Bold")
  
  data_to_plot <- 
    scenario |>
    dplyr::filter(sample == 0) |> 
    dplyr::select(
      sample,
      index, 
      discretionary_spending, 
      nondiscretionary_spending
    ) |>
    tidyr::pivot_longer(
      cols = c(
        discretionary_spending,
        nondiscretionary_spending
      ),
      names_to  = "type",
      values_to = "spending"
    ) |>
    dplyr::mutate(spending = spending / period_factor) |>
    dplyr::mutate(
      type = stringr::str_replace_all(type, "_spending", ""),
      type = dplyr::case_when(
        type == "discretionary"    ~ "discretionary", 
        type == "nondiscretionary" ~ "non-discretionary"
      ),
      type = factor(
        type, 
        levels = type_levels
      )
    ) 
  
  total_current_spending <- 
    data_to_plot |>
      dplyr::filter(index == min(index)) |> 
      dplyr::pull(spending) |> 
      sum()
  
  data_to_plot <-
    data_to_plot|> 
    dplyr::filter(type %in% !!type)

  current_year <- 
    data_to_plot |> 
    dplyr::filter(index == min(index)) 

  type_colors <- c(
    "discretionary"     = colors[4], 
    "non-discretionary" = colors[5]
  )

  current_year_spending <- current_year$spending
  names(current_year_spending) <- 
    glue::glue( 
      "<span style='color: {type_colors[as.character(current_year$type)]};'>**{as.character(current_year$type)}**</span> "
    )
    
  summarized_data <-
    data_to_plot |> 
    dplyr::group_by(type) |>
    dplyr::summarise(
      median_spending = stats::median(spending),
      max_spending    = max(spending),
      min_spending    = min(spending)
    )
  
 median_spending <- summarized_data$median_spending
 names(median_spending) <- format_colored_names(summarized_data, type_colors)
 
 max_spending <- summarized_data$max_spending
 names(max_spending) <- format_colored_names(summarized_data, type_colors)
 
 min_spending <- summarized_data$min_spending
 names(min_spending) <- format_colored_names(summarized_data, type_colors)

  max_y <- 
    ceiling(max(
      data_to_plot |> 
        dplyr::group_by(index) |>
        dplyr::summarise(spending = sum(spending)) |>
        dplyr::pull(spending) |> 
        max(na.rm = TRUE)
    ) / 1000) * 1000 
  
  min_y <- 
    ceiling(min(
      data_to_plot |> 
        dplyr::group_by(index, type) |>
        dplyr::summarise(spending = sum(spending)) |>
        dplyr::pull(spending) |> 
        min(na.rm = TRUE)
    ) / 1000) * 1000 

  min_length <- max(nchar(max_y)) + 2

  if (max_y > 10000) {
    y_breaks_factor <- 10000
  } else if (max_y > 1000) {
    y_breaks_factor <- 1000
  } else {  
    y_breaks_factor <- 100
  }
  
  y_breaks <- 
    seq(
      from = min(
        y_limits[1],
        min_y, 
        na.rm = TRUE
      ), 
      to = max(
        max_y,
        y_limits[2],
        na.rm = TRUE
      ), 
      by = y_breaks_factor
    )
  
  data_to_plot |> 
    ggplot2::ggplot(
      ggplot2::aes(
        x    = index, 
        y    = spending
        
      )
    ) +
    ggplot2::geom_area(
      ggplot2::aes(fill = type)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(title = "Spending type")
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "discretionary"     = colors[4], 
        "non-discretionary" = colors[5]
      )      
    ) +
    ggplot2::scale_y_continuous(
      labels = format_currency,
      breaks = y_breaks
    ) +
    ggplot2::labs(
      title = glue::glue("Expected Spending"),
      subtitle = glue::glue(paste0(
        paste_scenario_id(scenario),
        "Current spending: ",
        paste0(
          "<strong>", 
          glue::glue(
            "<span style='color: {type_colors[as.character(current_year$type)]};'>{format_currency(current_year_spending, accuracy = 1)}</span>"
          ),
          "</strong> ",
          collapse = " + "
        ),
        " = <strong>",
        format_currency(
          total_current_spending, 
          accuracy = 1,
          min_length = min_length
        ),
        "</strong>"
      )),
      caption = glue::glue(paste0(
        "Median spending: ",
        paste0(
          names(median_spending), 
          "=<strong>", 
          format_currency(
            median_spending, 
            accuracy = 1,
            min_length  = min_length
          ), 
          "</strong>",
          collapse = " & "
        ),
        "<br>",
        "Max spending: ",
        paste0(
          names(max_spending), 
          "=<strong>", 
          format_currency(
            max_spending, 
            accuracy = 1,
            min_length  = min_length
          ), 
          "</strong>",
          collapse = " & "
        ),
        "<br>",
        "Min spending: ",
        paste0(
          names(min_spending), 
          "=<strong>", 
          format_currency(
            min_spending, 
            accuracy = 1,
            min_length  = min_length
          ), 
          "</strong>",
          collapse = " & "
        ), 
        ""
      )),
      x = paste_year_index_axis_label(),
      y = glue::glue("Spending {period}"),
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, max(scenario$index), by = 10),
      labels = 
        function(breaks) paste_labels(breaks, scenario = scenario)
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
    ggplot2::coord_cartesian(ylim = c(y_limits[1], y_limits[2]))
}
