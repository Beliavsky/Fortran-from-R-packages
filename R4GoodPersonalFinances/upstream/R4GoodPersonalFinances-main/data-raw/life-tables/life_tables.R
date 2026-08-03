# HMD. Human Mortality Database. Max Planck Institute for Demographic Research (Germany), University of California, Berkeley (USA), and French Institute for Demographic Studies (France). Available at www.mortality.org (data downloaded on 2025-02-05).

countries <- c("Poland", "USA")

life_tables <- 
  purrr::map_dfr(countries, function(country) {
    read_hmd_life_tables(
      path = file.path("data-raw", "life-tables", tolower(country))
    ) |> 
    dplyr::mutate(country = !!country) |> 
    dplyr::select(country, dplyr::everything())
  })  

life_tables |> 
  dplyr::count(country)

usethis::use_data(
  life_tables, 
  overwrite = TRUE,
  internal  = FALSE
)
