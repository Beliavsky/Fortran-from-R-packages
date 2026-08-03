#' @import R6
NULL

############# Class RoundingHelper ###########################################
#' Helper object to define rounding rules for the InsuranceContract,
#' InsuranceTarif and ProfitParticipation classes.
#'
#' @description The class \code{RoundingHelper} provides the code and settings
#' to define numeric rounding rules for premiums, reserves, benefits etc. of
#' a life insurance contract. By default, no rounding it applied.
#'
#' @param values Contract values calculated so far (in the \code{contract$Values}
#'      list) then this method is called by the contract object
#'
#' @param premiumCalculationTime The time when the premiums should be
#'        (re-)calculated according to the equivalence principle. A time 0
#'        means the initial premium calculation at contract closing, later
#'        premium calculation times can be used to re-calculate the new
#'        premium after a contract change (possibly including an existing reserve)
#'
#' @examples
#' # TODO
#' @export
RoundingHelper = R6Class(
  "RoundingHelper",

  ######################### PUBLIC METHODS ##################################
  public  = list(
    #' @field rounding The (named) list containing all declared rounding definitions
    rounding  = list(),

    #' @description Initialize the rounding settings
    #' @details Sets up the rounding helper by giving a list of named entries, specifying rounding accuracy for each particular value
    #'
    #' @param ... named entries specifying rounding accuracy
    #' @examples
    #' rounding = RoundingHelper$new(raw = 0, hundred = -2, accurate = 4)
    #' rounding$round("raw", c(1234.567891, 0.00012345, 1234))
    #' rounding$round("hundred", c(1234.567891, 0.00012345, 1234))
    #' rounding$round("accurate", c(1234.567891, 0.00012345, 1234))
    #' rounding$round("non-existing", c(1234.567891, 0.00012345, 1234))
    initialize = function(...) {
      self$rounding = list(...)
    },

    #' @description Round the given values using the pre-defined accuracy
    #' @details Rounds the given values using the accuracies defined in the
    #'          internal rounding list (set either via the 'initialize' function
    #'          or via a call to 'setRounding'. The accuracies are defined using
    #'          a 'spec' identifier, which allows to define different accuracies
    #'          for different uses
    #'
    #' @param spec the ID used for looking up the desired accuracy
    #' @param value the values to be rounded according to 'spec'
    #' @param ... currently unused
    #' @examples
    #' rounding = RoundingHelper$new(raw = 0, hundred = -2, accurate = 4)
    #' rounding$round("raw", c(1234.567891, 0.00012345, 1234))
    #' rounding$round("hundred", c(1234.567891, 0.00012345, 1234))
    #' rounding$round("accurate", c(1234.567891, 0.00012345, 1234))
    #' # If the given spec does not exist, no rounding it applied
    #' rounding$round("non-existing", c(1234.567891, 0.00012345, 1234))
    round = function(spec, value, ...) {
      if (is.character(spec)) {
        spec = self$getRounding(spec, ...)
      }
      if (is.function(spec)) {
        spec(value)
      } else if (is.numeric(spec)) {
        round(value, digits = spec)
      } else {
        value
      }
    },


    #' @description Define rounding accuracy for a certain identifier
    #' @details Configures the rounding helper for a given named entry,
    #'          specifying rounding accuracy for each particular value
    #'
    #' @param key the ID used for looking up the desired accuracy
    #' @param spec the rounding accuracy (number of digits)
    #' @param ... currently unused
    #' @examples
    #' rounding = RoundingHelper$new(raw = 0, hundred = -2, accurate = 4)
    #' rounding$round("raw", c(1234.567891, 0.00012345, 1234))
    #' rounding$round("hundred", c(1234.567891, 0.00012345, 1234))
    #' rounding$round("accurate", c(1234.567891, 0.00012345, 1234))
    #' # If the given spec does not exist, no rounding it applied
    #' rounding$round("non-existing", c(1234.567891, 0.00012345, 1234))
    #' # Add a new spec with different settings:
    #' rounding$setRounding("non-existing", 1)
    #' rounding$round("non-existing", c(1234.567891, 0.00012345, 1234))
    setRounding = function(key, spec, ...) {
      self$rounding[[key]] = spec
    },

    #' @description Extract rounding accuracy for a certain identifier
    #' @details Read out the rounding for a given named entry.
    #'
    #' @param key the ID used for looking up the desired accuracy
    #' @param ... currently unused
    #' @examples
    #' rounding = RoundingHelper$new(raw = 0, hundred = -2, accurate = 4)
    #' rounding$getRounding("hundred")
    getRounding = function(key, ...) {
      self$rounding[[key]]
    }
  )
)



