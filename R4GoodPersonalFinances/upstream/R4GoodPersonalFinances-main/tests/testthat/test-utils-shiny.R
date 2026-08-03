test_that("sidebar footer returns current package version", {
  
  as.character(sidebar_footer()) |> 
    expect_match(
      utils::packageVersion("R4GoodPersonalFinances") |> 
        as.character()
    )
})

test_that("sidebar footer returns current app version", {
  
  as.character(sidebar_footer(
    app_version = "9.9.9"
  )) |> 
    expect_match("9.9.9")
})
