#' Household class
#' 
#' @description
#' The `Household` class aggregates information about 
#' a household and its members.
#' 
#' @returns An object of class `Household`.
#' 
#' @examples
#' household <- Household$new()
#' household$risk_tolerance
#' household$consumption_impatience_preference
#' household$smooth_consumption_preference
#' @export 
Household <- R6::R6Class(
  classname = "Household",

  public = list(

    #' @description 
    #' Printing the household object
    #' @param current_date A date in the format "YYYY-MM-DD".
    print = function(current_date = get_current_date()) {

      cli::cli_h1("Household")

      cli::cli_bullets(c(
        "i" = "Current date: {.val {current_date}}",
        "i" = "Household lifespan: {.val {self$get_lifespan(current_date = current_date)}} years",
        "i" = "Consumption impatience preference: {.val {self$consumption_impatience_preference}}",
        "i" = "Smooth consumption preference: {.val {self$smooth_consumption_preference}}",
        "i" = "Risk tolerance: {.val {self$risk_tolerance}}"
      ))

      cli::cli_h2("Household members")
      
      self$get_members() |> 
        purrr::walk(function(member) {
          print(member)
        })
      
    },

    #' @description 
    #' Getting members of the household
    get_members = function() {
      private$.household_members
    },

    #' @description 
    #' Adding a member to the household
    #' It will fail if a member with the same name already exists.
    #' @param household_member A `HouseholdMember` object.
    add_member = function(household_member) {

      if (household_member$get_name() %in% 
          names(private$.household_members)) {
        
        cli::cli_abort(c(
          "Household member already exists:",
          "x" = "There is already a member named {.value {household_member$get_name()}}"
        ))
      }

      self$set_member(member = household_member)
    },

    #' @description 
    #' Setting a member of the household
    #' If a member already exists, it will be overwritten.
    #' @param member A `HouseholdMember` object.
    set_member = function(member) {
      private$.household_members[[member$get_name()]] <- member
    },

    #' @description 
    #' Setting an arbitrary lifespan of the household
    #' @param value A number of years.
    set_lifespan = function(value) {
      private$.lifespan <- value
    },

    #' @description 
    #' Getting a lifespan of the household
    #' If not set, it will be calculated based on the members' lifespans.
    #' @param current_date A date in the format "YYYY-MM-DD".
    get_lifespan = function(current_date = get_current_date()) {
      
      current_date <- lubridate::as_date(current_date)

      if (!is.null(private$.lifespan)) {
        return(private$.lifespan)
      }

      if (length(private$.household_members) == 0) {
        cli::cli_abort(c(
          "No members in the household:",
          "x" = "There are no members added to the household!"
        ))
      }

      private$.household_members |> 
        purrr::map_dbl(function(household_member) {
          household_member$get_lifespan(current_date = current_date)
        }) |> 
        max() |> 
        ceiling()
    },

    #' @description 
    #' Calculating a survival rate of the household
    #' based on its members' parameters of the Gompertz model.
    #' @param current_date A date in the format "YYYY-MM-DD".
    calc_survival = function(current_date = get_current_date()) {

      current_date       <- lubridate::as_date(current_date)
      household_lifespan <- self$get_lifespan(current_date = current_date)
      members            <- self$get_members()

      members_params <- purrr::map(members, function(member) {

        if (is.null(member$mode)) {
          cli::cli_abort(c(
            "Missing Gompertz mode parameter for member: {member$get_name()}"
          ))
        }

        if (is.null(member$dispersion)) {
          cli::cli_abort(c(
            "Missing Gompertz dispersion parameter for member: {
              member$get_name()
            }"
          ))
        }

        list(
          name       = member$get_name(),
          age        = member$calc_age(current_date = current_date) |> round(0),
          mode       = member$mode,
          dispersion = member$dispersion,
          max_age    = member$max_age
        )
      })

      survival_rates <- 
        dplyr::tibble(
          year = 0:(household_lifespan)
        ) 
      
      members_ages <- 
        members_params |> 
        purrr::map_dbl(function(x) x$age) |> 
        unname() 

      min_age <- min(members_ages)

      for (member in members_params) {

        survival_rates <- 
          survival_rates |>
          dplyr::mutate(
            !!member$name := 
              calc_gompertz_survival_probability(
                current_age = member$age, 
                target_age  = member$age + year, 
                mode        = member$mode, 
                dispersion  = member$dispersion,
                max_age     = member$max_age
              )
          )
      }

      survival_rates <- 
        survival_rates |>
        dplyr::mutate(
          joint = 
            1 - purrr::pmap_dbl(
              .l = dplyr::select(survival_rates, dplyr::all_of(names(members_params))),
              .f = function(...) {
                prod(1 - c(...))
              }
            )
        )

      objective_fun <- function(params) {

        mode       <- params[1]
        dispersion <- params[2]
        
        approx_surv <- 
          calc_gompertz_survival_probability(
            current_age = min_age, 
            target_age  = min_age + survival_rates$year,
            mode        = mode, 
            dispersion  = dispersion
          )
        
        actual_surv <- survival_rates$joint
        sum((approx_surv - actual_surv) ^ 2)
      }

      members_modes <- 
        purrr::map_dbl(members_params, function(x) x$mode)
      members_dispersions <- 
        purrr::map_dbl(members_params, function(x) x$dispersion)
      
      init_params <- c(
        mode       = mean(members_modes), 
        dispersion = mean(members_dispersions)
      )
      
      params <- stats::optim(
        par = init_params, 
        fn  = objective_fun
      )

      mode       <- params$par[["mode"]]
      dispersion <- params$par[["dispersion"]]

      survival_rates <- 
        survival_rates |> 
        dplyr::mutate(
          gompertz = calc_gompertz_survival_probability(
            current_age = min_age,
            target_age  = min_age + year,
            mode        = mode,
            dispersion  = dispersion
          )
        )
        
      list(
        data        = survival_rates,
        mode        = mode,
        dispersion  = dispersion,
        current_age = min_age
      )
    },

    #' @description 
    #' Calculating a minimum age of the household members.
    #' @param current_date A date in the format "YYYY-MM-DD".
    get_min_age = function(current_date = get_current_date()) {

      current_date <- lubridate::as_date(current_date)
      
      min_age <- 
        private$.household_members |>
        purrr::map_dbl(function(x) x$calc_age(current_date = current_date)) |>
        min()
      min_age
    }
  ),

  active = list(

    #' @field expected_income Set of rules that are used to 
    #' generate streams of expected income 
    expected_income = function(value) {
      if (missing(value)) {
        return(private$.expected_income)
      } 
      private$.expected_income <- value
    },

    #' @field expected_spending Set of rules that are used to
    #' generate streams of expected spending
    expected_spending = function(value) {
      if (missing(value)) {
        return(private$.expected_spending)
      } 
      private$.expected_spending <- value
    },

    #' @field risk_tolerance Risk tolerance of the household
    risk_tolerance = function(value) {
      if (missing(value)) {
        return(private$.risk_tolerance)
      } 
      private$.risk_tolerance <- value
    },

    #' @field consumption_impatience_preference Consumption 
    #' impatience preference of the household - 
    #' subjective discount rate (rho).
    #' Higher values indicate a stronger preference for consumption today 
    #' versus in the future.
    consumption_impatience_preference = function(value) {
      if (missing(value)) {
        return(private$.consumption_impatience_preference)
      }
      private$.consumption_impatience_preference <- value
    },

    #' @field smooth_consumption_preference Smooth consumption 
    #' preference of the household - 
    #' Elasticity of Intertemporal Substitution (EOIS) (eta).
    #' Higher values indicate more flexibility and a lower preference 
    #' for smooth consumption.
    # Usually between 0 (no flexibility) and 1 (high level of flexibility).
    smooth_consumption_preference = function(value) {
      if (missing(value)) {
        return(private$.smooth_consumption_preference)
      }
      private$.smooth_consumption_preference <- value
    }

  ),

  private = list(

    .household_members                 = list(),
    .expected_income                   = list(),
    .expected_spending                 = list(),
    .lifespan                          = NULL,
    .risk_tolerance                    = 0.5,
    .consumption_impatience_preference = 0.04,
    .smooth_consumption_preference     = 1,
    
    deep_clone = function(name, value) {

      if (name != ".household_members") return(value)
        
      purrr::map(value, function(item) {
        item$clone(deep = TRUE)
      })
    }

  )
)


