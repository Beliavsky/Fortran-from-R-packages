#' Calculating Gompertz survival probability
#' 
#' @param current_age Current age
#' @param target_age Target age
#' @param mode Mode of the Gompertz distribution
#' @param dispersion Dispersion of the Gompertz distribution
#' @param max_age Maximum age. Defaults to `NULL`.
#' 
#' @returns A numeric. The probability of survival from 'current_age' 
#' to 'target_age' based on the Gompertz distribution 
#' with the given parameters.
#' 
#' @examples
#' calc_gompertz_survival_probability(
#'   current_age = 65, 
#'   target_age  = 85, 
#'   mode        = 80, 
#'   dispersion  = 10
#' )
#' @export

calc_gompertz_survival_probability <- function(
  current_age, 
  target_age, 
  mode, 
  dispersion,
  max_age = NULL
) {

  if (is.null(max_age)) {
    return(calc_gompertz_surv_prob(
      current_age = current_age,
      target_age  = target_age,
      mode        = mode,
      dispersion  = dispersion
    ))
  }
  
  probabilities <- 
    (
      calc_gompertz_survival_probability(
        current_age = current_age, 
        target_age  = target_age, 
        mode        = mode, 
        dispersion  = dispersion
      ) - 
      calc_gompertz_survival_probability(
        current_age = current_age, 
        target_age  = max_age, 
        mode        = mode, 
        dispersion  = dispersion
      ) 
    ) / (
      1 - 
        calc_gompertz_survival_probability(
          current_age = current_age, 
          target_age  = max_age, 
          mode        = mode, 
          dispersion  = dispersion
        )
    )

  probabilities[target_age > max_age]  <- 0
  probabilities[is.nan(probabilities)] <- 0
  probabilities
}

calc_gompertz_surv_prob <- function(
  current_age, 
  target_age, 
  mode, 
  dispersion
) {

  exp(
    exp((current_age - mode) / dispersion) * 
    (1 - exp((target_age - current_age) / dispersion))
  )
}

#' Calculating Gompertz model parameters
#' 
#' @param mortality_rates A data frame 
#' with columns `mortality_rate` and `age`.
#' Usually the output of [read_hmd_life_tables()] function or 
#' filtered data from [life_tables] object.
#' @param current_age A numeric. Current age.
#' @param estimate_max_age A logical. Should the maximum age be estimated?
#' 
#' @returns A list containing:
#'   \item{data}{The input mortality rates data frame with additional columns like 'survival_rate' and 'probability_of_death'}
#'   \item{mode}{The mode of the Gompertz distribution}
#'   \item{dispersion}{The dispersion parameter of the Gompertz distribution}
#'   \item{current_age}{The current age parameter}
#'   \item{max_age}{The maximum age parameter}
#' 
#' @references Blanchet, David M., and Paul D. Kaplan. 2013. "Alpha, Beta, and Now... Gamma." Journal of Retirement 1 (2): 29-45. \doi{10.3905/jor.2013.1.2.029}.
#' 
#' @examples
#' mortality_rates <- 
#'   dplyr::filter(
#'     life_tables,
#'     country == "USA" & 
#'     sex     == "male" &
#'     year    == 2022
#'   )
#'   
#' calc_gompertz_parameters(
#'   mortality_rates = mortality_rates,
#'   current_age     = 65
#' )
#' @export

calc_gompertz_parameters <- function(
  mortality_rates,
  current_age,
  estimate_max_age = FALSE
) {

  age <- mortality_rate <- survival_rate <- probability_of_death <- NULL

  mortality_rates <- 
    mortality_rates |>
    dplyr::filter(age >= !!current_age) |> 
    dplyr::mutate(
      survival_rate = cumprod(c(1, 1 - mortality_rate[-1]))
    ) |> 
    dplyr::mutate(
      probability_of_death = 
        dplyr::lag(survival_rate, default = 1) - survival_rate
  )

  mode <- 
    mortality_rates |> 
    dplyr::filter(probability_of_death == max(probability_of_death)) |> 
    dplyr::pull(age)
  
  gompertz_objective_fun <- function(x) {

    dispersion <- x[1]
    if (estimate_max_age) max_age <- x[2] else max_age <- NULL

    sum((
      calc_gompertz_survival_probability(
        current_age = current_age, 
        target_age  = mortality_rates$age, 
        max_age     = max_age,
        mode        = mode, 
        dispersion  = dispersion
      ) - mortality_rates$survival_rate
    ) ^ 2)
  }

if (estimate_max_age) {
  par <- c(10, 100) 
  lower = -Inf
  upper = Inf
  method = "Nelder-Mead"
} else {
  par <- c(10)
  lower = 0
  upper = 100
  method = "Brent"
}

  results <- stats::optim(
    par    = par,
    lower  = lower, 
    upper  = upper,
    method = method,
    fn     = gompertz_objective_fun
  )

  dispersion <- results$par[1]
  max_age    <- results$par[2] 

  if (is.na(max_age) || is.null(max_age)) {
    max_age <- NULL 
  } 

  list(
    data        = mortality_rates,
    mode        = mode,
    dispersion  = dispersion,
    current_age = current_age,
    max_age     = max_age
  )
}

#' Plotting the results of Gompertz model calibration
#' 
#' @param params A list returned by [calc_gompertz_parameters()] function.
#' @param mode A numeric. The mode of the Gompertz model.
#' @param dispersion A numeric. The dispersion of the Gompertz model.
#' @param max_age A numeric. The maximum age of the Gompertz model.
#' 
#' @return A [ggplot2::ggplot()] object showing the comparison between
#' actual survival rates from life tables and the fitted Gompertz model.
#' 
#' @examples
#' mortality_rates <- 
#'   dplyr::filter(
#'     life_tables,
#'     country == "USA" & 
#'     sex     == "female" &
#'     year    == 2022
#'   )
#'   
#' params <- calc_gompertz_parameters(
#'   mortality_rates = mortality_rates,
#'   current_age     = 65
#' )
#' 
#' plot_gompertz_calibration(params = params)
#' @export

plot_gompertz_calibration <- function(
  params,
  mode,
  dispersion,
  max_age
) {

  age <- survival_rate <- survival_rate_gompertz <- NULL

  if (missing(max_age))    max_age    <- params$max_age
  if (missing(mode))       mode       <- params$mode
  if (missing(dispersion)) dispersion <- params$dispersion

  data_to_plot <- 
    params$data |> 
    dplyr::mutate(
      survival_rate_gompertz = 
        calc_gompertz_survival_probability(
          target_age  = age, 
          current_age = params$current_age, 
          max_age     = max_age,
          mode        = mode, 
          dispersion  = dispersion
        )
    ) 
  
  colours                    <- PrettyCols::prettycols("Bold")
  real_survival_rate_col     <- colours[4]
  gompertz_survival_rate_col <- colours[5]
  value_colour               <- "grey40"
  
  data_to_plot |> 
    ggplot2::ggplot(
      ggplot2::aes(x = age)
    ) + 
    ggplot2::geom_point(
      ggplot2::aes(y = survival_rate), 
      color = real_survival_rate_col,
      size = 2.5,
    ) + 
    ggplot2::geom_line(
      ggplot2::aes(y = survival_rate_gompertz),
      color = gompertz_survival_rate_col,
      linewidth = 1
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.minor = ggplot2::element_line(color = "grey90"),
      plot.caption = 
        ggtext::element_markdown(
          color = "grey60", 
          size  = 10
        ),
      plot.subtitle = ggtext::element_markdown(color = "grey60")
    ) +
    ggplot2::labs(
      title    = "Gompertz Model Calibration",
      subtitle = glue::glue("<span style='color: {real_survival_rate_col};'>**Life tables**</span> vs <span style='color: {gompertz_survival_rate_col};'>**Gompertz**</span> survival rates for 
      {unique(params$data$sex)} in {unique(params$data$country)} as of 
      {unique(params$data$year)}

      "),
      y        = "Survival rate", 
      x        = "Age",
      caption = glue::glue(paste(
        "*Current age*: <span style='color: {value_colour};'>**{params$current_age}**</span>;",
        "*Max age*: <span style='color: {value_colour};'>**{ifelse(is.null(max_age), 'NULL', round(max_age,1))}**</span>;",
        "*Mode*: <span style='color: {value_colour};'>**{mode}**</span>;",
        "*Dispersion*: <span style='color: {value_colour};'>**{round(dispersion,2)}**</span>.",
        ""
      ))
    )
}

#' Calculating the Gompertz model parameters for joint survival 
#' 
#' @param p1 A list with `age`, `mode` and `dispersion` parameters
#' for the first person (p1).
#' @param p2 A list with `age`, `mode` and `dispersion` parameters
#' for the second person (p2).
#' @param max_age A numeric. The maximum age for the Gompertz model.
#' 
#' @return A list containing:
#'   \item{data}{A data frame with survival rates for 'p1', 'p2', 'joint' survival, and the fitted Gompertz model}
#'   \item{mode}{The mode of the joint Gompertz distribution}
#'   \item{dispersion}{The dispersion parameter of the joint Gompertz distribution}
#'   
#' @examples
#' calc_gompertz_joint_parameters(
#'   p1 = list(
#'     age        = 65,
#'     mode       = 88,
#'     dispersion = 10.65
#'   ),
#'   p2 = list(
#'     age        = 60,
#'     mode       = 91,
#'     dispersion = 8.88
#'   ),
#'   max_age = 110
#' )
#' @export

calc_gompertz_joint_parameters <- function(
  p1 = list(
    age        = NULL,
    mode       = NULL,
    dispersion = NULL
  ),
  p2 = list(
    age        = NULL,
    mode       = NULL,
    dispersion = NULL
  ),
  max_age = 120
) {

  year <- NULL

  survival_rates <- 
    dplyr::tibble(
      year = 0:(max_age - min(p1$age, p2$age)),
      p1 = calc_gompertz_survival_probability(
        current_age = p1$age, 
        target_age  = p1$age + year, 
        mode        = p1$mode, 
        dispersion  = p1$dispersion
      ),
      p2 = calc_gompertz_survival_probability(
        current_age = p2$age, 
        target_age  = p2$age + year, 
        mode        = p2$mode, 
        dispersion  = p2$dispersion
      ),
      joint = p1 + p2 - p1 * p2
    ) 
  
  min_age <- min(p1$age, p2$age)

  objective_fun <- function(params) {

    mode       <- params[1]
    dispersion <- params[2]
    
    approx_surv <- 
      calc_gompertz_survival_probability(
        current_age = min_age, 
        target_age  = min_age + survival_rates$year,
        mode        = mode, 
        dispersion  = dispersion
      )
    
    actual_surv <- survival_rates$joint
    sum((approx_surv - actual_surv) ^ 2)
  }
  
  init_params <- c(
    mode       = mean(c(p1$mode, p2$mode)), 
    dispersion = mean(c(p1$dispersion, p2$dispersion))
  )
  
  params <- stats::optim(
    par = init_params, 
    fn  = objective_fun
  )

  mode       <- params$par[["mode"]]
  dispersion <- params$par[["dispersion"]]

  survival_rates <- 
    survival_rates |> 
    dplyr::mutate(
      gompertz = calc_gompertz_survival_probability(
        current_age = min_age,
        target_age  = min_age + year,
        mode        = mode,
        dispersion  = dispersion
      )
    )
    
  list(
    data       = survival_rates,
    mode       = mode,
    dispersion = dispersion
  )
} 

#' Plotting the results of Gompertz model calibration for joint survival
#' 
#' @param params A list returned by [calc_gompertz_joint_parameters()] function.
#' @param include_gompertz A logical. Should the Gompertz survival curve be included in the plot?
#' 
#' @return A [ggplot2::ggplot()] object showing the survival probabilities
#' for two individuals and their joint survival probability.
#' 
#' @examples
#' params <- calc_gompertz_joint_parameters(
#'   p1 = list(
#'     age        = 65,
#'     mode       = 88,
#'     dispersion = 10.65
#'   ),
#'   p2 = list(
#'     age        = 60,
#'     mode       = 91,
#'     dispersion = 8.88
#'   ),
#'   max_age = 110
#' )
#' 
#' plot_joint_survival(params = params, include_gompertz = TRUE)
#' @export

plot_joint_survival <- function(
  params, 
  include_gompertz = FALSE
) {

  year <- value <- name <- NULL

  fixed_colors <- c(
    "p1"  = PrettyCols::prettycols("Bold")[1],
    "p2"  = PrettyCols::prettycols("Bold")[4],
    "joint"  = PrettyCols::prettycols("Bold")[3],
    "gompertz"  = PrettyCols::prettycols("Bold")[5]
  )

  cols <- c("p1", "p2", "joint")
  if (include_gompertz) {
    cols <- c(cols,  "gompertz")
  }

  params$data |> 
    tidyr::pivot_longer(
      cols = dplyr::all_of(cols),
      names_to = "name",
      values_to = "value"
    ) |> 
    ggplot2::ggplot(
      ggplot2::aes(
        x = year, 
        y = value,
        color = name
      )
    ) + 
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, by = 0.1),
      labels = scales::percent
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      title = "Survival Probability Comparison",
      x = "Years from now",
      y = "Survival probability",
      color = NULL
    ) + 
    ggplot2::scale_color_manual(values = fixed_colors) 
}
