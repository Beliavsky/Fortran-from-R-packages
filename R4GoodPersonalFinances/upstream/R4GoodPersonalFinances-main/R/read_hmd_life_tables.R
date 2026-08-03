#' Reading HMD life tables
#' 
#' @param path A character. Path to the folder with life tables.
#' @param files A character. Names of files with life tables.
#' 
#' @return A data frame containing mortality data with columns:
#'   \item{sex}{Character - sex ('male', 'female', or 'both')}
#'   \item{year}{Integer - the year of the data}
#'   \item{age}{Integer - age}
#'   \item{mortality_rate}{Numeric - mortality rate}
#'   \item{life_expectancy}{Numeric - life expectancy}
#' 
#' @references HMD. Human Mortality Database. Max Planck Institute for Demographic Research (Germany), University of California, Berkeley (USA), and French Institute for Demographic Studies (France). Available at www.mortality.org 
#' 
#' @examples
#' \dontrun{
#' # Download 'txt' files 
#' # ("mltper_1x1.txt", "fltper_1x1.txt", "bltper_1x1.txt") 
#' # for a given country to the working directory
#' # from https://www.mortality.org after registration.
#' 
#' read_hmd_life_tables(path = getwd())
#' }
#' @export
#' 
read_hmd_life_tables <- function(
  path  = getwd(),
  files = c(
    "mltper_1x1.txt", 
    "fltper_1x1.txt", 
    "bltper_1x1.txt"
  )
) {

  Year <- Age <- qx <- ex <- age <- NULL
  
  purrr::map(files, function(file) {

    if (!fs::file_exists(file.path(path, file))) {
      warning(call. = FALSE, glue::glue(
          "File {file} does NOT exists. Skipping..."
        ))
      return(NULL)
    }

    readr::read_table(
      file.path(path, file), 
      skip = 2,
      col_types = readr::cols(
        Year = readr::col_integer(),
        Age  = readr::col_character(),
        mx   = readr::col_double(),
        qx   = readr::col_double(),
        ax   = readr::col_double(),
        lx   = readr::col_double(),
        dx   = readr::col_double(),
        Lx   = readr::col_double(),
        Tx   = readr::col_double(),
        ex   = readr::col_double()
      )
    )  
  }) |> 
    purrr::set_names(files) |> 
    dplyr::bind_rows(.id = "file") |> 
    dplyr::transmute(
      sex = dplyr::case_when(
        file == "mltper_1x1.txt" ~ "male",
        file == "fltper_1x1.txt" ~ "female",
        file == "bltper_1x1.txt" ~ "both"
      ),
      year            = Year,
      age             = Age,
      mortality_rate  = qx,
      life_expectancy = ex
    ) |> 
    dplyr::mutate(
      age = stringr::str_replace(age, stringr::fixed("+"), "")
    ) |> 
    dplyr::mutate( age = as.integer(age)) 
}
