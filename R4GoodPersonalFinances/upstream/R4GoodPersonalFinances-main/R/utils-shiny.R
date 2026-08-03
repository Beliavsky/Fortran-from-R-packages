sidebar_footer <- function(
  package_version = utils::packageVersion("R4GoodPersonalFinances"),
  app_version     = apps$package_version
) {

  if (package_version == app_version) {
    version_number <- as.character(package_version)
  } else {
    version_number <- paste0(
      "(pkg version) ",
      as.character(package_version), 
      shiny::br(),
      "(app version) ",
      as.character(app_version)
    )
  }

  shiny::tagList(
    shiny::div(
      style = "text-align: center; font-size: 11px;",
      shiny::tags$a(
        "R4Good.Academy",
        href = "https://r4good.academy/",
        target = "_blank"
      ),
      shiny::br(),
      shiny::tags$a(
        "R4GoodPersonalFinances",
        href = "https://r4goodacademy.github.io/R4GoodPersonalFinances",
        target = "_blank"
      ),
      shiny::br(),
      shiny::HTML(version_number)
    )
  )

}
