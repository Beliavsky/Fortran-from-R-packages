renv::status()
renv::clean()
renv::upgrade()
renv::update()
renv::snapshot()

devtools::dev_sitrep()
usethis::use_release_issue()
usethis::use_version()
# usethis::use_version(which = "patch")
# usethis::use_version(which = "minor")


rstudioapi::restartSession()
pak::local_install(upgrade = FALSE)
source("data-raw/internal.R")
R4GoodPersonalFinances:::apps$package_version

rstudioapi::restartSession()
pak::local_install(upgrade = FALSE)
R4GoodPersonalFinances:::apps$package_version
R4GoodPersonalFinances::run_app(shinylive = TRUE)
R4GoodPersonalFinances::run_app()
R4GoodPersonalFinances::run_app(which = "purchasing-power")
R4GoodPersonalFinances::run_app("retirement-ruin", shinylive = TRUE)
R4GoodPersonalFinances::run_app("purchasing-power", shinylive = TRUE)
R4GoodPersonalFinances::run_app("risk-adjusted-returns", shinylive = TRUE)

# badgecreatr::badge_lifecycle("experimental")
# badgecreatr::badge_projectstatus("active")


devtools::document()
devtools::build_readme()
pkgdown::build_site(preview = FALSE)
pkgdown::preview_site()

# roxygen2md::roxygen2md()

spelling::spell_check_package()
spelling::update_wordlist()

urlchecker::url_check()
# urlchecker::url_update()

# usethis::use_tidy_description()

# devtools::test_coverage()
withr::with_envvar(
  c(
    NOT_CRAN = "true",   # prevents skip_on_cran()
    CI = "false"         # prevents skip_on_ci()
  ),
  covr::package_coverage(quiet = FALSE, type = "tests")
)
# R4GoodPersonalFinances Coverage: 93.94%

# usethis::use_github_action()
# usethis::use_pkgdown_github_pages()


citation <- readLines("inst/CITATION")
citation
current_year <- lubridate::year(Sys.Date())
current_year <- 2025
citation[5] <- glue::glue("  year     = {current_year},")
citation[6] <- glue::glue("  note     = \"{sprintf(\"R package version %s, https://r4goodacademy.github.io/R4GoodPersonalFinances/\", as.character(packageVersion(\"R4GoodPersonalFinances\")))} \",")
citation
writeLines(citation, "inst/CITATION")
devtools::load_all()
citation('R4GoodPersonalFinances')

usethis::use_cran_comments()

devtools::check(remote = TRUE, manual = TRUE)

# usethis::use_release_issue(version = "1.0.0")


R4GoodPersonalFinances:::apps$package_version == packageVersion("R4GoodPersonalFinances")

devtools::dev_sitrep()


devtools::check_win_devel()

devtools::submit_cran()

usethis::use_github_release()
usethis::use_dev_version(push = TRUE)
