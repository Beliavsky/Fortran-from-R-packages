#' Get default Gompertz parameters
#' 
#' Calculates default Gompertz parameters for a given age, country and sex,
#' based on the package build-in HMD life tables.
#' 
#' @param age A numeric. The age of the individual.
#' @param country A character. The name of the country.
#' @param sex A character. The sex of the individual. 
#' 
#' @seealso [calc_gompertz_parameters()]
#' 
#' @returns A list containing:
#'   \item{mode}{The mode of the Gompertz distribution}
#'   \item{dispersion}{The dispersion parameter of the Gompertz distribution}
#'   \item{current_age}{The current age parameter}
#'   \item{max_age}{The maximum age parameter}
#' @examples
#' get_default_gompertz_parameters(
#'   age     = 65,
#'   country = "USA",
#'   sex     = "male"
#' )
#' @export
get_default_gompertz_parameters <- function(
  age, 
  country = unique(life_tables$country), 
  sex     = c("both", "male", "female")
) {

  life_tables <- life_tables 
  year <- NULL

  country <- rlang::arg_match(country)
  sex     <- rlang::arg_match(sex)

  mortality_rates <- 
    life_tables |> 
    dplyr::filter(
      country == !!country & 
      sex     == !!sex
    ) |>
    dplyr::filter(year == max(year)) 

  gompertz_parameters <-
    mortality_rates |> 
    calc_gompertz_parameters(
      current_age      = age,
      estimate_max_age = TRUE
    )
}
