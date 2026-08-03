generate_cashflow_streams <- function(
  timeline, 
  triggers
) {
  
  result <- purrr::map(triggers, function(trigger_formulas) {
    
    trigger_formulas <- c(trigger_formulas, "TRUE ~ 0")
    
    formulas <- purrr::map(trigger_formulas, rlang::parse_expr)

    streams <- vector("numeric", nrow(timeline))
    
    for (i in seq_len(nrow(timeline))) {
      
      eval_env <- rlang::env(.parent = parent.frame())
      
      eval_env$members <- timeline$members[i, ]
      
      eval_env$is_on <- function(member, event) {

        members_data <- get("members", envir = eval_env)
        if (!member %in% names(members_data)) {
          stop(paste("Member", member, "not found"))
        }
        
        member_data <- members_data[[member]]
        if (!"events" %in% names(member_data)) {
          stop(paste("No events found for member", member))
        }
        if (!event %in% names(member_data$events)) {
          stop(paste("Event", event, "not found for member", member))
        }
        
        member_data$events[[event]]$on
      }

      eval_env$is_not_on <- function(member, event) {
        !eval_env$is_on(member, event)
      }

      eval_env$age <- function(member) {
        
        members_data <- get("members", envir = eval_env)
        if (!member %in% names(members_data)) {
          stop(paste("Member", member, "not found"))
        }
        
        members_data[[member]]$age
      }
      
      row_result <- rlang::eval_tidy(
        rlang::expr(dplyr::case_when(!!!formulas)),
        env = eval_env
      )
      
      streams[i] <- row_result
    }
    
    return(streams)
  }) |>
  dplyr::as_tibble()
  
  return(result)
}
