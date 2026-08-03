#' Plot survival of household members
#' 
#' @description 
#' Plot survival probabilities for each household members and
#' for the entire household when at least one member is alive.
#' The household joint survival probability is also
#' approximated by a Gompertz model.
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
#' hm3 <- 
#'  HouseholdMember$new(
#'    name       = "member3",
#'    birth_date = "1975-01-01",
#'    mode       = 88,
#'    dispersion = 7.77
#'  )
#' household <- Household$new()
#' household$add_member(hm1)
#' household$add_member(hm2)
#' household$add_member(hm3) 
#' 
#' plot_survival(
#'  household    = household, 
#'  current_date = "2020-01-01"
#' )
#' @export
plot_survival <- function(
  household, 
  current_date = get_current_date()
) {

  id <- value <- year <- name <- type <- NULL

  current_date   <- lubridate::as_date(current_date)
  members        <- household$get_members()
  members_params <- purrr::map(members, function(member) {
    list(
      name       = member$get_name(),
      age        = member$calc_age(current_date = current_date) |> round(0),
      mode       = member$mode,
      dispersion = member$dispersion
    )
  })

  data_to_plot <- 
    household$calc_survival(current_date = current_date)$data

  bold_colors     <- PrettyCols::prettycols("Bold")
  num_individual  <- length(members_params)
  required_colors <- num_individual + 2 
   
  all_colors <- bold_colors[1:required_colors]
  individual_names <- names(members_params)
  individual_colors <- all_colors[1:num_individual]
  names(individual_colors) <- individual_names

  fixed_colors <- c(
    individual_colors,
    "joint"    = all_colors[num_individual + 1],
    "gompertz" = all_colors[num_individual + 2]
  )

  cols <- c(names(members_params), "joint", "gompertz")

  data_to_plot_long <- 
    data_to_plot |> 
    tidyr::pivot_longer(
      cols = dplyr::all_of(cols),
      names_to = "name",
      values_to = "value"
    ) |> 
    dplyr::mutate(
      type = dplyr::case_when(
        name %in% c("joint", "gompertz") ~ "at_least_one",
        name %in% names(members_params) ~ "individual"
      )
    )
  
  data_to_plot_long|> 
    ggplot2::ggplot(
      ggplot2::aes(
        x     = year, 
        y     = value,
        color = name
      )
    ) + 
    ggplot2::geom_line(ggplot2::aes(linetype = type)) + 
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, by = 0.1),
      labels = scales::percent
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, household$get_lifespan(), by = 10),
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "Survival Probability of Household and Household Members",
      x = "Years from now",
      y = "Survival probability",
      color = NULL
    )  + 
    ggplot2::scale_color_manual(
      values = fixed_colors,
    ) +
    ggplot2::scale_linetype_manual(
      values = c("at_least_one" = "dashed", "individual" = "solid"),
      guide = "none"
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(
          linetype = c("dashed", "dashed", rep("solid", num_individual))
        )
      )
    )
}
