generate_household_timeline <- function(
  household, 
  current_date
) {
  index <- end_age <- start_age <- NULL

  current_date <- lubridate::as_date(current_date)
  max_lifespan <- household$get_lifespan(current_date = current_date)

  timeline <-
    dplyr::tibble(
      index         = as.integer(seq_len(max_lifespan + 1) - 1),
      years_left    = as.integer(max_lifespan - index),
      date          = current_date + lubridate::years(index),
      year          = as.integer(lubridate::year(date)),
      survival_prob = 
        household$calc_survival(current_date = current_date)$data$gompertz
    )
  
  members <- 
    household$get_members() |> 
    purrr::map(function(member) {
      
      member_specific <- 
        dplyr::tibble(
          age = member$calc_age(current_date = timeline$date) |>  round(0)
        ) 
      
      events <- member$get_events()
      if (length(events) > 0) {
        
        events <- 
          names(events) |> 
          purrr::map(function(event_name) {
            
            on <- 
              member_specific$age >= events[[event_name]]$start_age &
              member_specific$age <= events[[event_name]]$end_age

            dplyr::tibble(
              on        = on,
              off       = !on,
              start_age = events[[event_name]]$start_age,
              end_age   = events[[event_name]]$end_age,
              years     = end_age - start_age + 1
            ) 
          }) |> 
          purrr::set_names(names(events)) |> 
          dplyr::as_tibble() 

        member_specific <- 
          member_specific |>
          dplyr::mutate(events = events)
      }
      member_specific
    }) |> 
    dplyr::as_tibble()

  timeline |> 
    dplyr::mutate(members = members)
}
