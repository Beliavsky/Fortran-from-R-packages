apps_folders <- c(
  "risk-adjusted-returns",
  "purchasing-power",
  "retirement-ruin"
)

apps <- list()

for (app_folder in apps_folders) {

  print(app_folder)

  apps[[app_folder]] <- list()
  apps[[app_folder]]$ui <- 
    readLines(
      system.file(
        "apps", app_folder, "ui.R", 
        package = "R4GoodPersonalFinances"
      )
    )
  apps[[app_folder]]$server <- 
    readLines(
      system.file(
        "apps", app_folder, "server.R", 
        package = "R4GoodPersonalFinances"
      )
    )
}

apps$package_version <- utils::packageVersion("R4GoodPersonalFinances")
apps


# Table 1: "Annuity 2000 Basic Table" in Johansen, Robert J. 1998. “Annuity 2000 Mortality Tables.” Transactions, Society of Actuaries, pp. 264-290. http://www.soa.org/library/research/transactions-reports-of-mortality-moribidity-and-experience/1990-99/1995/january/TSR9510.pdf.

mortality_rates_males <-
  dplyr::tibble(
    age = 65:115,
    mortality_rate = c(
      10.993, 12.188, 13.572, 15.160, 16.946, 18.920, 21.071, 23.388, 25.871,
      28.552, 31.477, 34.686, 38.225, 42.132, 46.427, 51.128, 56.250, 61.809,
      67.826, 74.322, 81.326, 88.863, 96.958, 105.631, 114.858, 124.612,
      134.861, 145.575, 156.727, 168.290, 180.245, 192.565, 205.229, 218.683,
      233.371, 249.741, 268.237, 289.305, 313.391, 340.940, 372.398, 408.210,
      448.823, 494.681, 546.231, 603.917, 668.186, 739.483, 818.254, 904.945,
      1000.000
    )
  ) |> 
  dplyr::mutate(
    mortality_rate = mortality_rate / 1000,
    sex            = "male",
    country        = "USA",
    year           = 2000
  ) 
  
mortality_rates_females <-
  dplyr::tibble(
    age = 65:115,
    mortality_rate = c(
      7.017, 7.734, 8.491, 9.288, 10.163, 11.165, 12.339, 13.734, 15.391,
      17.326, 19.551, 22.075, 24.910, 28.074, 31.612, 35.580, 40.030, 45.017,
      50.600, 56.865, 63.907, 71.815, 80.682, 90.557, 101.307, 112.759, 
      124.733, 137.054, 149.552, 162.079, 174.492, 186.647, 198.403, 210.337,
      233.027, 237.051, 252.985, 271.406, 292.893, 318.023, 347.373, 381.520,
      421.042, 466.516, 518.520, 577.631, 644.427, 719.484, 803.380, 896.693,
      1000.000
    )
  ) |> 
  dplyr::mutate(
    mortality_rate = mortality_rate / 1000,
    sex            = "female",
    country        = "USA",
    year           = 2000
  ) 
  
test_mortality_rates <- 
  dplyr::bind_rows(
    mortality_rates_males,
    mortality_rates_females
  )




usethis::use_data(
  apps, 
  test_mortality_rates,
  overwrite = TRUE, 
  internal  = TRUE
)
