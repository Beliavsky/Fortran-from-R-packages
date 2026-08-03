#' Calculate Life Expectancy
#' 
#' @param current_age A numeric. Current age.
#' @param mode A numeric. Mode of the Gompertz distribution.
#' @param dispersion A numeric. Dispersion of the Gompertz distribution.
#' @param max_age A numeric. Maximum age. Defaults to 120.
#' 
#' @return A numeric. Total life expectancy in years.
#' 
#' @examples
#' calc_life_expectancy(
#'   current_age = 65, 
#'   mode        = 80, 
#'   dispersion  = 10
#' )
#' @export
calc_life_expectancy <- function(
  current_age,
  mode,
  dispersion,
  max_age = 120
) {

  target_ages <- (current_age + 1):max_age

  survival_probs <- 
    purrr::map_dbl(target_ages, function(target_age) {
      calc_gompertz_survival_probability(
        current_age = current_age,
        target_age  = target_age,
        mode        = mode,
        dispersion  = dispersion,
        max_age     = max_age
      )
    })

  life_expectancy <- sum(survival_probs)
  life_expectancy <- life_expectancy + 0.5
  life_expectancy <- life_expectancy + current_age
  life_expectancy
}

#' Calculate Gompertz mode for a given life expectancy
#' 
#' @param life_expectancy A numeric. Desired life expectancy.
#' @param current_age A numeric. Current age.
#' @param dispersion A numeric. Dispersion of the Gompertz distribution.
#' @param max_age A numeric. Maximum age. Defaults to 120.
#' 
#' @return A numeric. Mode of the Gompertz distribution.
#' 
#' @examples
#' calc_gompertz_mode(
#'  life_expectancy = 86,
#'   current_age    = 25,
#'   dispersion     = 8.88,
#'   max_age        = 115
#' )
#' @export
calc_gompertz_mode <- function(
  life_expectancy,
  current_age,
  dispersion,
  max_age = 120
) {

    mode_interval <- c(current_age + 1, 150)
    tolerance     <- 1e-5
  
    objective_function <- function(
      life_expectancy, 
      mode, 
      dispersion, 
      current_age, 
      max_age
    ) {

      calculated_life_expectancy <- 
        calc_life_expectancy(
          mode        = mode,
          current_age = current_age,
          dispersion  = dispersion,
          max_age     = max_age
        )
      
      return(calculated_life_expectancy - life_expectancy )
    }
  
    result <- stats::uniroot(
      f               = objective_function,
      interval        = mode_interval,
      current_age     = current_age,
      dispersion      = dispersion,
      life_expectancy = life_expectancy,
      max_age         = max_age,
      tol             = tolerance,
      extendInt       = "yes"
    )
  
  return(result$root)
}
