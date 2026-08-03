#' Plot life expectancy of household members
#' 
#' @description
#' Probability of dying at a given age is plotted for each member 
#' of a household.
#' Also for each member the life expectancy is shown
#' as dashed vertical line. 
#' 
#' 
#' @inheritParams simulate_scenario
#' @returns A `ggplot` object.
#' @examples
#' hm1 <- 
#'  HouseholdMember$new(
#'    name       = "member1",
#'    birth_date = "1955-01-01",
#'    mode       = 88,
#'    dispersion = 10.65
#'  )
#' hm2 <- 
#'  HouseholdMember$new(
#'    name       = "member2",
#'    birth_date = "1965-01-01",
#'    mode       = 91,
#'    dispersion = 8.88
#'  )
#' household <- Household$new()
#' household$add_member(hm1)
#' household$add_member(hm2)
#' 
#' plot_life_expectancy(household = household)
#' @export
plot_life_expectancy <- function(
  household
) {

  id <- density <- age <- life_expectancy <- NULL

  if ("Household" %in% class(household)) {

    members <- household$get_members()

    params <- purrr::map(members, function(member) {
      list(
        mode       = member$mode,
        dispersion = member$dispersion,
        life_expectancy = member$calc_life_expectancy()
      )
    })
  } else {
    params <- household
  }

  colors <- 
    grDevices::colorRampPalette(
      PrettyCols::prettycols("Bold")
    )(
      length(params)
    )

  min_age <- 0
  max_age <- 150

  dx <- 0.1
  ages <- seq(min_age, max_age, by = dx)
  
  densities <- 
    do.call(rbind, lapply(names(params), function(id) {
      
      m <- params[[id]]$mode
      d <- params[[id]]$dispersion
      data.frame(
        id      = id,
        age     = ages,
        density = calc_gompertz_pdf(ages, m, d)
      )
    }))
  densities$density[is.na(densities$density)] <- 0

  plot_range <- 
    densities |> 
    dplyr::group_by(id) |> 
    dplyr::filter(density >= 0.0005) |> 
    dplyr::reframe(age = range(age)) |> 
    dplyr::ungroup() |> 
    dplyr::summarise(min_age = min(age), max_age = max(age)) 

  if ("Household" %in% class(household)) {

    life_expectancies <- 
      params |> 
      dplyr::bind_rows(.id = "id") 

  } else {

    life_expectancies <- 
      do.call(
        rbind, 
        lapply(split(densities, densities$id), 
        function(df) {
          data.frame(
            id               = unique(df$id),
            life_expectancy  = sum(df$age * df$density) * dx
          )
        })
      )
    }

  max_density <- max(densities$density, na.rm = TRUE)
  label_y     <- max_density * 0.9
  
  ggplot2::ggplot(
    densities,
    ggplot2::aes(
      x     = age, 
      y     = density, 
      color = id
    )
  ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::geom_vline(
      data = life_expectancies,
      ggplot2::aes(
        xintercept = life_expectancy, 
        color = id,
        linetype = "Life Expectancy"
      ),
      linewidth = 0.8,
      show.legend = TRUE
    ) +
    ggplot2::geom_text(
      data = life_expectancies,
      ggplot2::aes(
        x     = life_expectancy, 
        label = round(life_expectancy, 1), 
        color = id
      ),
      y           = label_y,
      angle       = 90, 
      vjust       = -0.5,
      show.legend = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(
        floor(plot_range$min_age / 10) * 10, 
        ceiling(plot_range$max_age / 10) * 10, 
        by = 10
      ),
      limits = c(
        floor(plot_range$min_age / 10) * 10, 
        ceiling(plot_range$max_age / 10) * 10
      )
    ) +
    ggplot2::labs(
      x     = "Age",
      y     = "Density (probability of dying)",
      title = "Probability Distribution of Age of Death and Life Expectancy"
    ) +
    ggplot2::scale_linetype_manual(values = c("Life Expectancy" = "dashed")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.title     = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.caption = 
        ggtext::element_markdown(
          color = "grey60", 
          size  = 10
        ),
      plot.subtitle = ggtext::element_markdown(color = "grey60")
    )
}

calc_gompertz_pdf <- function(x, mode, dispersion) {
  eta   <- 1 / dispersion
  theta <- eta * exp(-mode * eta)
  theta * exp(eta * x) * exp(- (theta / eta) * (exp(eta * x) - 1))
}
