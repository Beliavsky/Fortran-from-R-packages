#' @include HelperFunctions.R InsuranceParameters.R ProfitParticipation.R
#'
#' @import MortalityTables
#' @import R6
#' @importFrom lubridate year month years days year<-
#' @importFrom objectProperties setSingleEnum
NULL

#' An enum specifying the main characteristics of the tarif.
#'
#' @description Possible values are:
#' \describe{
#'   \item{annuity}{Whole life or term annuity (periodic survival benefits)
#'        with flexible payouts (constand, increasing, decreasing, arbitrary,
#'        etc.)}
#'   \item{wholelife}{A whole or term life insurance with only death benefits.
#'        The benefit can be constant, increasing, decreasing, described by
#'        a function, etc.}
#'   \item{endowment}{An  endowment with death and survival benefits,
#'        potentially with different benefits.}
#'   \item{pureendowment}{A pure endowment with only a survival benefit at
#'        the end of the contract. Optionally, in case of death, all or part
#'        of the premiums paid may be refunded.}
#'   \item{terme-fix}{A terme-fix insurance with a fixed payout at the end
#'        of the contract, even if the insured dies before that time.
#'        Premiums are paid until death of the insured.}
#'   \item{dread-disease}{A dread-disease insurance, which pays in case of
#'        a severe illness (typically heart attacks, cancer, strokes, etc.),
#'        but not in case of death.}
#'   \item{endowment + dread-disease}{A combination of an endowment and a
#'        temporary dread-disease insurance. Benefits occur either on death,
#'        severe illness or survival, whichever comes first.}
#' }
#' @export
TariffTypeEnum = objectProperties::setSingleEnum(
  "TariffType",
  levels = c(
    "annuity",
    "wholelife",
    "endowment",
    "pureendowment",
    "terme-fix",
    "dread-disease",
    "endowment + dread-disease"
  ))
setValidity("TariffTypeSingleEnum", function(object) {
  if (length(object) != 1L) {
    "Only one tariff type can be given"
  } else if (!object %in% levels(object)) {
    paste("Tarif type '", object, "' does not exist. Valid tariff types are:",
          paste0("\n('", paste0(levels(object), collapse = "', '"),
                 "')"), sep = "")
  }
});


############# Class InsuranceTarif ###########################################
#' Base class for traditional Insurance Tarifs (with fixed guarantee, profit
#' sharing and no unit-linked component)
#'
#' @description The class \code{InsuranceTarif} provides the code and general
#' framework to implement contract-independent functionality of a life insurance
#' product.
#'
#' @details This is a base class for holding contract-independent values and
#' providing methods to calculate cash flows, premiums, etc. Objects of this
#' class do NOT contain contract-specific values like age, death probabilities,
#' premiums, reserves, etc. Rather, they are the calculation kernels that will
#' be called by the \code{\link{InsuranceContract}} objects to make the actual,
#' tariff-specific calculations.
#'
#' Most methods of this class are not meant to be called manually, but are supposed
#' to be called by the InsuranceContract object with contract-specific information.
#' The only methods that are typically used for defining an insurance tariff are
#' the constructor \href{#method-new}{\code{InsuranceTarif$new()}} and the cloning method
#' \href{#method-createModification}{\code{InsuranceTarif$createModification()}}.
#' All other methods should never be called manually.
#'
#' However, as overriding private methods is not possible in an R6 class, all the
#' methods need to be public to allow overriding them in derived classes.
#'
# # Parameters for the constructors
#' @param name The unique name / ID of the tariff
#' @param type An enum specifying the main characteristics of the tarif. See [TariffTypeEnum]
#' @param tarif The tariff's public name to be stored in the `tarif` field.
#' @param desc A short human-readable description to be stored in the `desc` field.
# # General parameters for (almost) all function
#' @param params Contract-specific, full set of parameters of the contract
#'      (merged parameters of the defaults, the tariff, the profit participation
#'      scheme and the contract)
#' @param values Contract values calculated so far (in the \code{contract$Values}
#'      list) then this method is called by the contract object
#'
#' @param premiumCalculationTime The time when the premiums should be
#'        (re-)calculated according to the equivalence principle. A time 0
#'        means the initial premium calculation at contract closing, later
#'        premium calculation times can be used to re-calculate the new
#'        premium after a contract change (possibly including an existing reserve)
#'
#' @import MortalityTables
#' @examples
#' # Define an insurance tariff for 10-year endowments, using a guaranteed interest
#' # rate of 1% and the Austrian population mortality table of the census 2011.
#' # Premiums are paid monthly in advance during the whole contract period.
#' MortalityTables::mortalityTables.load("Austria_Census")
#' # Cost structure:
#' #   - 4% up-front acquisition costs (of premium sum)
#' #   - 1% collection cost of each premium paid
#' #   - 1%o yearly administration cost (of the sum insured) as long as premiums are paid
#' #   - 2%o yearly administration cost for paid-up contracts
#' #   - 10 Euro yearly unit costs (as long as premiums are paid)
#' costs.endw = initializeCosts(alpha = 0.04, beta = 0.01, gamma = 0.001,
#'     gamma.paidUp = 0.002, gamma.premiumfree = 0.002, unitcosts = 10)
#'
#' endowment.AT1 = InsuranceTarif$new(
#'     name = "Endow AT 1%", type = "endowment", tarif = "Austrian Endowment",
#'     desc = "An endowment for Austrian insured with 1% interest and no profit participation",
#'     ContractPeriod = 10,
#'     i = 0.01, mortalityTable = mort.AT.census.2011.unisex,
#'     costs = costs.endw, premiumFrequency = 12)
#'
#' # The instantiation of the actual contract will provide the contract specific
#' # information and immediately calculate all further values:
#' ctr.end.AT1 = InsuranceContract$new(tarif = endowment.AT1,
#'     contractClosing = as.Date("2020-07-01"), age = 42)
#'
#' # All values for the contract are already calculated during construction and
#' # stored in the ctr.end.AT1$Values list:
#' ctr.end.AT1$Values$basicData
#' ctr.end.AT1$Values$transitionProbabilities
#' ctr.end.AT1$Values$cashFlowsCosts
#' ctr.end.AT1$Values$presentValues
#' ctr.end.AT1$Values$premiums
#' ctr.end.AT1$Values$reserves
#' ctr.end.AT1$Values$premiumComposition
#' # etc.
#' @export
InsuranceTarif = R6Class(
  "InsuranceTarif",

  ######################### PUBLIC METHODS ##################################
  public  = list(
    #' @field name The tariff's unique name. Will also be used as the key for exported data.
    name  = "Insurance Contract Type",
    #' @field tarif The tariff's public name (typically a product name), not necessarily unique.
    tarif = NULL,
    #' @field desc A short human-readable description of the tariff and its main features.
    desc  = NULL,
    #' @field tariffType An enum specifying the main characteristics of the tarif.
    #' Possible values are:
    #' \describe{
    #'   \item{annuity}{Whole life or term annuity (periodic survival benefits)
    #'        with flexible payouts (constand, increasing, decreasing, arbitrary,
    #'        etc.)}
    #'   \item{wholelife}{A whole or term life insurance with only death benefits.
    #'        The benefit can be constant, increasing, decreasing, described by
    #'        a function, etc.}
    #'   \item{endowment}{An  endowment with death and survival benefits,
    #'        potentially with different benefits.}
    #'   \item{pureendowment}{A pure endowment with only a survival benefit at
    #'        the end of the contract. Optionally, in case of death, all or part
    #'        of the premiums paid may be refunded.}
    #'   \item{terme-fix}{A terme-fix insurance with a fixed payout at the end
    #'        of the contract, even if the insured dies before that time.
    #'        Premiums are paid until death of the insured.}
    #'   \item{dread-disease}{A dread-disease insurance, which pays in case of
    #'        a severe illness (typically heart attacks, cancer, strokes, etc.),
    #'        but not in case of death.}
    #'   \item{endowment + dread-disease}{A combination of an endowment and a
    #'        temporary dread-disease insurance. Benefits occur either on death,
    #'        severe illness or survival, whichever comes first.}
    #' }
    tariffType = TariffTypeEnum("wholelife"),

    #' @field Parameters A data structure (nested list) containing all relevant
    #' parameters describing a contract, its underlying tariff, the profit
    #' participation scheme etc. See [InsuranceContract.ParameterStructure] for
    #' all fields.
    Parameters = InsuranceContract.ParameterStructure,

    #' @description Initialize a new tariff object
    #' @details The constructor function defines a tariff and generates the
    #' corresponding data structure, which can then be used with the [InsuranceContract]
    #' class to define an actual contract using the tariff.
    #'
    #' The arguments passed to this function will be stored inside the
    #' \code{Parameters} field of the class, inside one of the lists sublists.
    #' The parameters are stacked from different layers (higher levels override
    #' default values from lower layers):
    #'
    #' * InsuranceContract object (parameters passed directly to the individual
    #'     contract)
    #' * ProfitParticipation object (parameters for profit participation, passed
    #'     to the definition of the profit plan, which is used for the tarif
    #'     definition or the contract)
    #' * InsuranceTarif object (parameters passed to the definition of the tariff
    #'     that was used for the contract)
    #' * Defaults taken from [InsuranceContract.ParameterStructure]
    #'
    #' The general implementation of this parameter layering means that (a) a tariff
    #' can already provide default values for contracts (e.g. a default maturity,
    #' default sum insured, etc) and (b) individual contracts can override all
    #' parameters defined with the underlying tariff. In particular the latter
    #' feature has many use-cases in prototyping: E.g. when you have a tariff
    #' with a guaranteed interest rate of 1\% and a certain mortality table,
    #' one can immediately instantiate a contract with an updated interest rate
    #' or mortality table for comparison. There is no need to re-implement a
    #' tariff for such comparisons, as long as only parameters are changed.
    #'
    #' @param ... Parameters for the [InsuranceContract.ParameterStructure],
    #'            defining the characteristics of the tariff.
    #' @import MortalityTables
    #' @examples
    #' MortalityTables::mortalityTables.load("Austria_Annuities_AVOe2005R")
    #' tarif.male = InsuranceTarif$new(name = "Annuity Males", type = "annuity",
    #'     i = 0.01, mortalityTable = AVOe2005R.male)
    initialize = function(name = NULL, type = "wholelife", tarif = "Generic Tarif", desc = "Description of tarif", ...) {
      if (getOption('LIC.debug.Tarif.init', FALSE)) {
        browser();
      }
      if (!missing(name))           self$name = name;
      if (!missing(type)) {
        self$tariffType = TariffTypeEnum(type)
      }
      if (!missing(tarif))          self$tarif = tarif;
      if (!missing(desc))           self$desc = desc;

      # Set the passed arguments as tariff parameters
      self$Parameters = InsuranceContract.ParametersFill(self$Parameters, ...)

      # Use the profit participation's parameters as fallback for initialized parameters
      ppScheme = self$Parameters$ProfitParticipation$profitParticipationScheme;
      if (!is.null(ppScheme)) {
          self$Parameters$ProfitParticipation = InsuranceContract.ParametersFallback(self$Parameters$ProfitParticipation, ppScheme$Parameters)
      }

      # Fill all remaining uninitialized values with their defaults, except for profit participation params
      self$Parameters = InsuranceContract.ParametersFallback(self$Parameters, InsuranceContract.ParameterDefaults, ppParameters = FALSE);

      # Properly set up the rounding helper
      if (is.null(self$Parameters$Hooks$Rounding)) self$Parameters$Hooks$Rounding = list()
      if (is.list(self$Parameters$Hooks$Rounding)){
        self$Parameters$Hooks$Rounding = do.call(RoundingHelper$new, self$Parameters$Hooks$Rounding)
      }
    },

    #' @description create a copy of a tariff with certain parameters changed
    #' @details This method \code{createModification} returns a copy of the tariff
    #' with all given arguments changed in the tariff's `InsuranceTarif$Parameters`
    #' parameter list.
    #'
    #' As InsuranceTarif is a R6 class with reference logic, simply assigning
    #' the object to a new variable does not create a copy, but references the
    #' original tariff object. To create an actual copy, one needs to call this
    #' method, which first clones the whole object and then adjusts all parameters
    #' to the values passed to this method.
    #'
    #' @param tariffType An enum specifying the main characteristics of the tarif.
    #'       See [TariffTypeEnum]
    #' @param ... Parameters for the [InsuranceContract.ParameterStructure],
    #'            defining the characteristics of the tariff.
    #' @import MortalityTables
    #' @examples
    #' MortalityTables::mortalityTables.load("Austria_Annuities_AVOe2005R")
    #' tarif.male = InsuranceTarif$new(name = "Annuity Males", type = "annuity",
    #'     i = 0.01, mortalityTable = AVOe2005R.male)
    #' tarif.unisex = tarif.male$createModification(name = "Annuity unisex",
    #'     mortalityTable = AVOe2005R.unisex)
    createModification = function(name  = NULL, tarif = NULL, desc  = NULL, tariffType = NULL, ...) {
      if (getOption('LIC.debug.createModification', FALSE)) {
        browser();
      }
      cloned = self$clone();
      if (!missing(name))       cloned$name = name;
      if (!missing(tarif))      cloned$tarif = tarif;
      if (!missing(desc))       cloned$desc = desc;
      if (!missing(tariffType)) cloned$tariffType = tariffType;

      cloned$Parameters = InsuranceContract.ParametersFill(cloned$Parameters, ...);
      # Properly set up the rounding helper
      if (is.null(cloned$Parameters$Hooks$Rounding)) cloned$Parameters$Hooks$Rounding = list()
      if (is.list(cloned$Parameters$Hooks$Rounding)){
        cloned$Parameters$Hooks$Rounding = do.call(RoundingHelper$new, cloned$Parameters$Hooks$Rounding)
      }
      cloned
    },

    #' @description Retrieve the parameters for this tariff (can be overridden
    #' for each contract)
    #'
    #' @examples
    #' tarif.male = InsuranceTarif$new(name = "Annuity Males", type = "annuity",
    #'     i = 0.01, mortalityTable = AVOe2005R.male)
    #' tarif.male$getParameters()
    getParameters = function() {
      self$Parameters
    },

    #' @description Get some internal parameters cached (length of data.frames,
    #' policy periods cut at max.age, etc.)
    #'
    #' @details This methos is not meant to be called explicitly, but rather used
    #' by the InsuranceContract class. It returns a list of maturities and ages
    #' relevant for the contract-specific calculations
    #'
    #' @param ... currently unused
    getInternalValues = function(params, ...) {
      if (getOption('LIC.debug.getInternalValues', FALSE)) {
        browser();
      }
      age = params$ContractData$technicalAge
      maxAge = MortalityTables::getOmega(params$ActuarialBases$mortalityTable)
      policyTerm = min(maxAge + 1 - age, params$ContractData$policyPeriod)
      # Length of CF vectors / data.frames is always 1 more than the policy period, since there might be a survival payment at the end!
      list(
        l = policyTerm + 1, # Length of CF vectors (1 larger than policy period!)
        # maxAge is the last age with a given death probability, so it can be included in the policy term!
        # The policy must end AFTER that final year, not at the beginning!
        policyTerm = policyTerm,
        premiumTerm = min(policyTerm, params$ContractData$premiumPeriod)
      )
    },


    #' @description Calculate the contract-relevant age(s) given a certain
    #' parameter data structure (contract-specific values)
    #'
    #' @details This method is not meant to be called explicitly, but rather used
    #' by the InsuranceContract class. It returns the relevant ages during the
    #' whole contract period
    getAges = function(params) {
      if (getOption('LIC.debug.getAges', FALSE)) {
        browser();
      }
      ages = ages(params$ActuarialBases$mortalityTable, YOB = year(params$ContractData$birthDate));
      # Append one last age, just in case we need to pad the qx with one more entry:
      ages = c(ages, tail(ages, 1) + 1)
      age = params$ContractData$technicalAge;
      if (age > 0) {
        ages = ages[-age:-1];
      }
      ages
    },

    #' @description Calculate the transition probabilities from the contract-specific
    #'  parameters passed as \code{params} and the already-calculated contract
    #'  values \code{values}
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    getTransitionProbabilities = function(params, values) {
      if (getOption('LIC.debug.getTransitionProbabilities', FALSE)) {
        browser();
      }
      if (!is.null(values$transitionProbabilities)) {
        values$transitionProbabilities
      } else {
        age = params$ContractData$technicalAge;
        ages = self$getAges(params);
        qx = MortalityTables::deathProbabilities(params$ActuarialBases$mortalityTable, YOB = year(params$ContractData$birthDate), ageDifferences = params$ContractData$ageDifferences);
        if (age > 0) {
          qx    = qx[-age:-1];
        }
        qx = pad0(qx, length(ages), value = 1)
        if (!is.null(params$ActuarialBases$invalidityTable)) {
          ix = MortalityTables::deathProbabilities(params$ActuarialBases$invalidityTable, YOB = year(params$ContractData$birthDate), ageDifferences = params$ContractData$ageDifferences);
          if (age > 0) {
            ix    = ix[-age:-1];
          }
        } else {
          ix = rep(0, length(qx));
        }
        ix = pad0(ix, length(qx));
        # invalidity/disease does NOT end the contract if flag is set!
        if (params$ActuarialBases$invalidityEndsContract) {
          px = 1 - qx - ix
        } else {
          px = 1 - qx
        }
        rd = params$Hooks$Rounding
        df = data.frame(age = ages, qx = rd$round("qx", qx), ix = rd$round("ix", ix), px = rd$round("px", px), row.names = ages - age)
        df
      }
    },

    #' @description Obtain the cost structure from the cost parameter and the
    #' given paremeter set
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    #' The cost parameter can be either an array of costs (generated by [initializeCosts()])
    #' or a function with parameters \code{param} and \code{values}(=NULL) returning
    #' an array of the required dimensions. This function makes sures that the
    #' latter function is actually evaluated.
    #'
    #' @param params The parameters of the contract / tariff
    getCostValues = function(params) {
      if (getOption('LIC.debug.getCostValues', FALSE)) {
        browser();
      }
      costs = valueOrFunction(params$Costs, params = params, values = NULL)
      costs = applyHook(params$Hooks$adjustCosts, costs, params = params, values = NULL);
      baseCost = valueOrFunction(params$minCosts, params = params, values = NULL, costs = costs)
      baseCost = applyHook(params$Hooks$adjustMinCosts, baseCost, costs = costs, params = params, values = NULL);
      if (!is.null(baseCost)) {
        costWaiver = valueOrFunction(params$ContractData$costWaiver, params = params, values = NULL, costs = costs, minCosts = baseCost)
        if (is.numeric(costWaiver)) {
          costs = costs * (1 - costWaiver) + baseCost * costWaiver
        } else if (is.boolean(costWaiver)) {
          if (isTRUE(costWaiver)) {
            costs = baseCost
          }
        }
      }
      costs
    },

    #' @description Returns the unit premium cash flow for the whole contract period.
    #'   - For constant premiums it will be rep(1, premiumPeriod),
    #'   - for single premiums it will be c(1, 0, 0, ...),
    #'   - for increasing premiums it will be (1+increase)^(0:(premiumPeriod-1))
    #' and 0 after the premium period
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    #' @param len The desired length of the returned data frame (the number of contract periods desire)
    getPremiumCF = function(len, params, values) {
      if (getOption('LIC.debug.getPremiumCF', FALSE)) {
        browser();
      }
      premPeriod = min(params$ContractData$premiumPeriod, params$ContractData$policyPeriod, len);
      if (premPeriod <= 0) {
        return(rep(0, len))
      }
      if (is.null(params$ContractData$premiumIncrease)) {
        pad0(rep(1, premPeriod - 1), len)
      } else {
        inc = valueOrFunction(params$ContractData$premiumIncrease, premiumPeriod = premPeriod, params = params, values = values)
        if (is.vector(inc) && length(inc) > 1) {
          # If premiumIncrease is (or returns) a vector, treat it as
          # relative premium amounts, ie. c(1, 1.1, 1.2) means +10% of
          # the initial premium for the second and third year
          pad0(inc, len)
        } else {
          pad0(inc ^ (0:(premPeriod - 1)), len)
        }
      }
    },

    #' @description Returns the unit annuity cash flow (guaranteed and contingent) for
    #'     the whole annuity payment period (after potential deferral period)
    #'   - For constant annuity it will be rep(1, annuityPeriod),
    #'   - for increasing annuities it will be (1+increase)^(0:(premiumPeriod-1))
    #' and 0 after the premium period
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    #' @param len The desired length of the returned data frame (the number of contract periods desire)
    getAnnuityCF = function(len, params, values) {
      if (getOption('LIC.debug.getAnnuityCF', FALSE)) {
        browser();
      }
      annuityPeriod = min(params$ContractData$policyPeriod - params$ContractData$deferralPeriod, len);
      if (is.null(params$ContractData$annuityIncrease)) {
        pad0(rep(1, annuityPeriod), len);
      } else {
        inc = valueOrFunction(params$ContractData$annuityIncrease, annuityPeriod = annuityPeriod, params = params, values = values)
        if (is.vector(inc) && length(inc) > 1) {
          # If annuityIncrease is (or returns) a vector, treat it as
          # relative annuity amounts, ie. c(1, 1.1, 1.2) means +10% of
          # the initial annuity for the second and third year
          pad0(inc, len)
        } else {
          # a numeric value means constant yearly increases (multiplicative)
          pad0(inc ^ (0:annuityPeriod), len)
        }
      }
    },

    #' @description Returns the unit death cash flow for the whole protection
    #' period (after potential deferral period!)
    #'   - For constant death benefit it will be rep(1, policyPeriod),
    #'   - for linearly decreasing sum insured it will be (policyPeriod:0)/policyPeriod
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    #' @param len The desired length of the returned data frame (the number of contract periods desired)
    getDeathCF = function(len, params, values) {
      if (getOption('LIC.debug.getDeathCF', FALSE)) {
        browser();
      }
      period = params$ContractData$policyPeriod - params$ContractData$deferralPeriod;
      if (is.null(params$ContractData$deathBenefit)) {
        pad0(rep(1, period), len)
      } else {
        benefit = valueOrFunction(params$ContractData$deathBenefit, len = len, params = params, values = values)
        if (is.vector(benefit) && length(benefit) > 1) {
          # If deathBenefit is (or returns) a vector, treat it as
          # relative annuity amounts, ie. c(1, 1.1, 1.2) means +10% of
          # the initial annuity for the second and third year
          pad0(benefit, len)
        } else {
          # constant death benefit
          pad0(rep(benefit, period), len)
        }
      }
    },

    #' @description Returns the unit survival cash flow profile for the whole contract
    #' period (after potential deferral period!)
    #'   - a single numeric value indicates a single survival payment at the end of the contract
    #'   - a vector of numeric values indicates potentially multiple survival payments for the whole contract period (paddded with 0 to the full contract length if shorter)
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    #' @param len The desired length of the returned data frame (the number of contract periods desired)
    getSurvivalCF = function(len, params, values) {
      if (getOption('LIC.debug.getSurvivalCF', FALSE)) {
        browser();
      }
      benefit = valueOrFunction(params$ContractData$survivalBenefit, len = len, params = params, values = values)
      if (is.null(benefit)) {
        benefit = 1
      }

      if (is.vector(benefit) && length(benefit) == 1) {
        c(rep(0, len - 1), benefit)
      } else {
        # If survivalBenefit is (or returns) a vector, treat it as yearly
        # survival payments, pad it to the desired length
        pad0(benefit, len)
      }
    },

    #' @description Returns the basic (unit) cash flows associated with the type
    #' of insurance given in the InsuranceTarif's `tariffType` field
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    getBasicCashFlows = function(params, values) {
      if (getOption('LIC.debug.getBasicCashFlows', FALSE)) {
        browser();
      }
      deferralPeriod = params$ContractData$deferralPeriod;
      guaranteedPeriod = params$ContractData$guaranteedPeriod;

      zeroes = pad0(0, values$int$l)

      cf = data.frame(
        guaranteed = zeroes,
        survival = zeroes,
        death = zeroes,
        disease = zeroes,
        sumInsured = rep(1, values$int$l)
      );
      if (self$tariffType == "annuity") {
        annuityPeriod = values$int$policyTerm - deferralPeriod;
        annuityCF = self$getAnnuityCF(len = annuityPeriod, params = params, values = values)
        # guaranteed payments exist only with annuities (first n years of the payment)
        cf$guaranteed = pad0(
          c(
            rep(0, deferralPeriod),
            head(annuityCF, n = guaranteedPeriod)
          ), values$int$l);
        cf$survival = pad0(c(
          rep(0, deferralPeriod + guaranteedPeriod),
          if (guaranteedPeriod > 0) tail(annuityCF, n = -guaranteedPeriod) else annuityCF,
          0), values$int$l)

        # start current contract block after deferral period
#        cf$sumInsured = pad0(rep(#c(rep(0, deferralPeriod), annuityCF) # <= TODO!!!


      } else if (self$tariffType == "terme-fix") {
        # Begin of bock does not have any influence
        cf$guaranteed = c(rep(0, values$int$policyTerm), 1)

      } else if (self$tariffType == "dread-disease") {
        # potential Payments start after deferral period
        cf$disease = c(
          rep(0, deferralPeriod),
          rep(1, values$int$l - 1 - deferralPeriod),
          0)
      } else {
        # For endowments, use the death factor here in the basic death CF
        # to fix the relation of death to survival benefit
        deathCF = self$getDeathCF(values$int$l - 1 - deferralPeriod, params = params, values = values)

        if (self$tariffType == "endowment" || self$tariffType == "pureendowment" || self$tariffType == "endowment + dread-disease") {
          cf$survival = pad0(self$getSurvivalCF(len = values$int$l, params = params, values = values), values$int$l)
        }
        if (self$tariffType == "endowment" || self$tariffType == "wholelife" || self$tariffType == "endowment + dread-disease") {
          cf$death = c(rep(0, deferralPeriod), deathCF, 0)
          # cf$sumInsured = c(rep(0, deferralPeriod), deathCF, 1);
        }
        if (self$tariffType == "endowment + dread-disease") {
          cf$disease = c(
            rep(0, deferralPeriod),
            rep(1, values$int$l - 1 - deferralPeriod),
            0);
        }
      }
      cf
    },

    #' @description Returns the cash flows for the contract given the parameters
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    getCashFlows = function(params, values) {
      if (getOption('LIC.debug.getCashFlows', FALSE)) {
        browser();
      }
      age = params$ContractData$technicalAge;

      if (is.null(values$cashFlowsBasic)) {
        values$cashFlowsBasic = self$getBasicCashFlows(params, values);
      }
      cflen = values$int$l
      zeroes = pad0(0, cflen)
      ages = pad0(self$getAges(params), cflen);
      cf = data.frame(
        premiums_advance   = zeroes,
        premiums_arrears   = zeroes,
        additional_capital = zeroes,
        guaranteed_advance = zeroes,
        guaranteed_arrears = zeroes,
        survival_advance   = zeroes,
        survival_arrears   = zeroes,
        death_SumInsured   = zeroes,
        disease_SumInsured = zeroes,
        death_GrossPremium = zeroes,
        death_Refund_past  = zeroes,
        death_PremiumFree  = zeroes,
        row.names          = ages - age
      );

      cf$additional_capital = pad0(params$ContractData$initialCapital, cflen)

      # Premiums:
      if (!params$ContractState$premiumWaiver) {
        premiums = self$getPremiumCF(len = cflen, params = params, values = values)
        if (params$ContractData$premiumPayments == "in advance") {
          cf$premiums_advance = premiums;
        } else {
          cf$premiums_arrears = premiums;
        }
      }

      # Survival Benefits
      if (params$ContractData$benefitPayments == "in advance") {
        cf$guaranteed_advance = pad0(values$cashFlowsBasic$guaranteed, cflen);
        cf$survival_advance = pad0(values$cashFlowsBasic$survival, cflen);
      } else {
        cf$guaranteed_arrears = pad0(values$cashFlowsBasic$guaranteed, cflen);
        cf$survival_arrears = pad0(values$cashFlowsBasic$survival, cflen);
      }

      # Death Benefits
      cf$death_SumInsured = pad0(values$cashFlowsBasic$death, cflen);
      if ((!is.null(params$Features$absPremiumRefund)) && (params$Features$absPremiumRefund > 0)) {
        cf$death_SumInsured = cf$death_SumInsured + pad0(padLast(params$Features$absPremiumRefund, min(cflen - 1, params$ContractData$premiumRefundPeriod)), cflen);
      }
      cf$disease_SumInsured = pad0(values$cashFlowsBasic$disease, cflen);
      cf$death_PremiumFree = cf$death_SumInsured;
      # premium refund
      if (params$ContractData$premiumRefund != 0) {
        totalpremiumcf = cf$premiums_advance + pad0(c(0, cf$premiums_arrears), cflen);

        # death benefit for premium refund is the sum of all premiums so far, but only during the premium refund period, afterwards it's 0:
        cf$death_GrossPremium = pad0(pad0(Reduce("+", totalpremiumcf[0:params$ContractData$policyPeriod], accumulate = TRUE), params$ContractData$premiumRefundPeriod), cflen)
        cf$death_Refund_past = cf$death_GrossPremium
        cf$death_Refund_past[(cf$death_GrossPremium > 0)] = 1;
      }

      applyHook(params$Hooks$adjustCashFlows, cf, params, values)
    },

    #' @description Returns the cost cash flows of the contract given the contract
    #'  and tariff parameters
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    getCashFlowsCosts = function(params, values) {
      if (getOption('LIC.debug.getCashFlowsCosts', FALSE)) {
        browser();
      }
      dm = dim(params$Costs);
      dmnames = dimnames(params$Costs);

      cf = array(
        0,
        dim = list(values$int$l, dm[1], dm[2], 3),
        dimnames = list(0:(values$int$l - 1), dmnames[[1]], dmnames[[2]], c("survival", "guaranteed", "after.death"))
      );
      cf[1,,,"survival"] = cf[1,,,"survival"] + params$Costs[,,"once"]
      if (values$int$premiumTerm > 0){
        for (i in 1:values$int$premiumTerm) {
          cf[i,,,"survival"] = cf[i,,,"survival"] + params$Costs[,,"PremiumPeriod"];
        }
      }
      if (values$int$premiumTerm < values$int$policyTerm) {
        for (i in (values$int$premiumTerm + 1):values$int$policyTerm) {
          cf[i,,,"survival"] = cf[i,,,"survival"] + params$Costs[,,"PremiumFree"];
        }
      }
      for (i in 1:values$int$policyTerm) {
        cf[i,,,"survival"] = cf[i,,,"survival"] + params$Costs[,,"PolicyPeriod"];

        # Guaranteed cost charged (charged no matter if the person is still alive or death).
        # Used mainly for term-fix or premium waivers upton death
        cf[i,,,"guaranteed"] = params$Costs[,,"FullContract"]
        cf[i,,,"after.death"] = params$Costs[,,"AfterDeath"]
      }

      # Costs charged only for a specific time (i.e. acquisition costs / commissions)
      # There are several conventions to scale alpha costs over the commision period:
      # a) alpha cost once (i.e. not distributed), not scaled
      # b) uniformly over the period (i.e. sum of k equal commisions is given as cost)
      # c) by present value (i.e. present value of k equal commission is given as cost)
      if (params$Features$alphaCostsCommission == "sum") {
          params$Costs[,,"CommissionPeriod"] = params$Costs[,,"CommissionPeriod"] / params$Loadings$commissionPeriod
      } else if (params$Features$alphaCostsCommission == "presentvalue") {
          # Use yearly constant premiums in advance, irrespective of the actual
          # contract. This is a simplification, but the actual present values
          # are calculated later, so for now we just assume standard parameters!
          len = params$Loadings$commissionPeriod;
          q = self$getTransitionProbabilities(params, values);
          px = pad0(c(1,q$px), len); # by defualt, premiums are in advance, so first payment has 100% probability
          v = 1/(1 + params$ActuarialBases$i)^((1:len)-1)
          params$Costs[,,"CommissionPeriod"] = params$Costs[,,"CommissionPeriod"] / sum(cumprod(px)*v)
      } else if (params$Features$alphaCostsCommission == "actual") {
          # NOP, nothing to do
      } else {
          warning("unrecognized value given for commissionPeriod: ",params$Features$alphaCostsCommission )
      }
      for (i in 1:params$Loadings$commissionPeriod) {
          cf[i,,,"survival"] = cf[i,,,"survival"] + params$Costs[,,"CommissionPeriod"];
      }

      # After premiums are waived, use the gamma_nopremiums instead of gamma:
      if (params$ContractState$premiumWaiver) {
        cf[,"gamma",,"survival"] = cf[,"gamma_nopremiums",,"survival"];
      }

      # some values like sumInsured or gross premium might change over time,
      # so multiply them with the unit cash flows stored in values$cashFlows
      cf[,,"SumInsured",] = cf[,,"SumInsured",] * values$cashFlowsBasic$sumInsured

      applyHook(params$Hooks$adjustCashFlowsCosts, cf, params, values)
    },

    #' @description Returns the present values of the cash flows of the contract
    #' (cash flows already calculated and stored in the \code{cashFlows} data.frame)
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    #' @param cashFlows data.frame of cash flows calculated by a call to \href{#method-getCashFlows}{\code{InsuranceTarif$getCashFlows()}}
    presentValueCashFlows = function(params, values) {
      if (getOption('LIC.debug.presentValueCashFlows', FALSE)) {
        browser();
      }
      rd = params$Hooks$Rounding

      qq = self$getTransitionProbabilities(params, values);

      i = params$ActuarialBases$i;
      v = 1/(1 + i);
      benefitFreqCorr = correctionPaymentFrequency(
        i = i, m = params$ContractData$benefitFrequency,
        order = valueOrFunction(params$ActuarialBases$benefitFrequencyOrder, params = params, values = values));
      premiumFreqCorr = correctionPaymentFrequency(
        i = i, m = params$ContractData$premiumFrequency,
          order = valueOrFunction(params$ActuarialBases$premiumFrequencyOrder, params = params, values = values));

      pvf = PVfactory$new(qx = qq, m = params$ContractData$benefitFrequency, mCorrection = benefitFreqCorr, v = v);

      pvRefund = pvf$death(values$cashFlows$death_GrossPremium);
      pvRefundPast = pvf$death(values$cashFlows$death_Refund_past) *
        (values$cashFlows[,"death_GrossPremium"] - values$cashFlows[,"premiums_advance"]);

      pv = cbind(
        premiums = rd$round("PV Premiums", pvf$survival(values$cashFlows$premiums_advance, values$cashFlows$premiums_arrears,
          m = params$ContractData$premiumFrequency, mCorrection = premiumFreqCorr)),
        additional_capital = rd$round("PV AdditionalCapital", pvf$survival(advance = values$cashFlows$additional_capital)),
        guaranteed = rd$round("PV Guarantee", pvf$guaranteed(values$cashFlows$guaranteed_advance, values$cashFlows$guaranteed_arrears)),
        survival = rd$round("PV Survival", pvf$survival(values$cashFlows$survival_advance, values$cashFlows$survival_arrears)),
        death_SumInsured = rd$round("PV Death", pvf$death(values$cashFlows$death_SumInsured)),
        disease_SumInsured = rd$round("PV Disease", pvf$disease(values$cashFlows$disease_SumInsured)),
        death_GrossPremium = rd$round("PV Death PremiumRefund", pvRefund),
        death_Refund_past = rd$round("PV Death PremiumRefund Past", pvRefundPast),
        death_Refund_future = rd$round("PV Death PremiumRefund Future", pvRefund - pvRefundPast),
        death_PremiumFree = rd$round("PV Death PremiumFree", pvf$death(values$cashFlows$death_PremiumFree))
      );

      rownames(pv) <- pad0(rownames(qq), values$int$l);
      applyHook(hook = params$Hooks$adjustPresentValues, val = pv, params = params, values = values)
    },

    #' @description Calculates the present values of the cost cash flows of the
    #' contract (cost cash flows alreay calculated by \href{#method-getCashFlowsCosts}{\code{InsuranceTarif$getCashFlowsCosts()}}
    #' and stored in the \code{values} list
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    #' @param presentValues The present values of the insurance claims (without costs)
    presentValueCashFlowsCosts = function(params, values, presentValues) {
      if (getOption('LIC.debug.presentValueCashFlowsCosts', FALSE)) {
        browser();
      }
      len = values$int$l;
      q = self$getTransitionProbabilities(params, values);
      i = params$ActuarialBases$i
      v = 1/(1 + i)
      pvf = PVfactory$new(qx = q, v = v)

      costs = values$cashFlowsCosts;
      pvc = costs * 0;

      # survival cash flows (until death)
      if (any(costs[,,,"survival"] != 0)) {
        pvc[,,,"survival"] = pvf$survival(advance = costs[,,,"survival"]);
      }
      # guaranteed cash flows (even after death)
      if (any(costs[,,,"guaranteed"] != 0)) {
        pvc[,,,"guaranteed"] = pvf$guaranteed(advance = costs[,,,"guaranteed"]);
      }

      # Cash flows only after death
      if (any(costs[,,,"after.death"] != 0)) {
        pvc[,,,"after.death"] = pvf$afterDeath(advance = costs[,,,"after.death"]);
      }

      # Costs based on actual benefits need to be calculated separately
      # (need to be recalculated with the benefit CFs multiplied by the cost
      # rate for each type and each year). In most cases, the cost rates will be
      # constant, but in general we can not assume this => need to re-calculate
      if (any(values$cashFlowsCosts[,,"Benefits",] != 0)) {
        bFreq = params$ContractData$benefitFrequency
        benefitFreqCorr = correctionPaymentFrequency(
          i = i, m = bFreq,
          order = valueOrFunction(params$ActuarialBases$benefitFrequencyOrder, params = params, values = values));
        pvfben = PVfactory$new(qx = q, m = bFreq, mCorrection = benefitFreqCorr, v = v);

        cfCosts = values$cashFlowsCosts[,,"Benefits",];
        cf = values$cashFlows;

        # Guaranteed + Survival + Death cover + disease
        pvc[,,"Benefits",] =
            pvfben$guaranteed(cf$guaranteed_advance * cfCosts, cf$guaranteed_arrears * cfCosts) +
            pvfben$survival(cf$survival_advance * cfCosts, cf$survival_arrears * cfCosts) +
            pvfben$death(cf$death_SumInsured * cfCosts) +
            pvfben$disease(cf$disease_SumInsured * cfCosts);
      }

      rd = params$Hooks$Rounding;
      pvc = rd$round("PV Costs", pvc)
      applyHook(hook = params$Hooks$adjustPresentValuesCosts, val = pvc, params = params, values = values, presentValues = presentValues)
    },

    #' @description Calculate the cash flows in monetary terms of the insurance contract
    #' @details Once the premiums of the insurance contracts are calculated, all
    #' cash flows can also be expressed in absolute terms. This function
    #' calculates these time series in monetary terms, once the premiums
    #' are calculated by the previous functions of this class.
    #'
    #' This method is NOT to be called directly, but implicitly by the [InsuranceContract] object.
    getAbsCashFlows = function(params, values) {
      if (getOption('LIC.debug.getAbsCashFlows', FALSE)) {
        browser();
      }
      rd = params$Hooks$Rounding;

      # TODO: Set up a nice list with coefficients for each type of cashflow,
      # rather than multiplying each item manually (this also mitigates the risk
      # of forgetting a dimension, because then the dimensions would not match,
      # while here it's easy to overlook a multiplication)
      # Multiply each CF column by the corresponding basis
      #
      # All propSI cash flows are already set up with the correct multiple
      # of the sumInsured (in cashFlowsBasic) for non-constant sums insured.
      # So here, we don't need to multiply with  values$cashFlowsBasic$sumInsured!
      propGP = c("premiums_advance", "premiums_arrears");
      propSI = c("guaranteed_advance", "guaranteed_arrears",
                 "survival_advance", "survival_arrears", "death_SumInsured",
                 "death_PremiumFree", "disease_SumInsured");
      propPS = c("death_GrossPremium", "death_Refund_past");
      values$cashFlows[,propGP] = values$cashFlows[,propGP] * values$premiums[["gross"]];
      values$cashFlows[,propSI] = values$cashFlows[,propSI] * params$ContractData$sumInsured;
      values$cashFlows[,propPS] = values$cashFlows[,propPS] * values$premiums[["gross"]] * params$ContractData$premiumRefund;

      # Sum all death-related payments to "death"  and remove the death_GrossPremium column
      values$cashFlows[,"death_SumInsured"] = values$cashFlows[,"death_SumInsured"] + values$cashFlows[,"death_GrossPremium"]
      colnames(values$cashFlows)[colnames(values$cashFlows) == "death_SumInsured"] = "death";
      # cashFlows[,"death_GrossPremium"] = NULL;
      values$cashFlows = rd$round("CF absolute", values$cashFlows);


      # costs relative to sumInsured are already set up as the correct multiple
      # of the original SI, including the dynamic changes over time!
      values$cashFlowsCosts = values$cashFlowsCosts[,,"SumInsured",] * params$ContractData$sumInsured +
        values$cashFlowsCosts[,,"SumPremiums",] * values$unitPremiumSum * values$premiums[["gross"]] +
        values$cashFlowsCosts[,,"GrossPremium",] * values$premiums[["gross"]] +
        values$cashFlowsCosts[,,"NetPremium",] * values$premiums[["net"]] +
        # values$cashFlowsCosts[,,"Benefits",] * TODO!!!
        values$cashFlowsCosts[,,"Constant",];
      values$cashFlowsCosts = rd$round("CF costs absolute", values$cashFlowsCosts);

      # Handle survival CF differently, because we don't want ".survival" in the column names!
      cbind(values$cashFlows, values$cashFlowsCosts[,,"survival"], values$cashFlowsCosts[,,-1])
    },

    #' @description Calculate the absolute present value time series of the insurance contract
    #' @details Once the premiums of the insurance contracts are calculated, all
    #' present values can also be expressed in absolute terms. This function
    #' calculates these time series in monetary terms, once the premiums and the
    #'  unit-benefit present values are calculated by the previous functions of
    #'  this classe.
    #'
    #' This method is NOT to be called directly, but implicitly by the [InsuranceContract] object.
    getAbsPresentValues = function(params, values) {
      if (getOption('LIC.debug.getAbsPresentValues', FALSE)) {
        browser();
      }
      pv = values$presentValues;

      #pv[,"age"] = pv[,"premiums"];
      #colnames(pv)[colnames(pv)=="age"] = "premiums.unit";

      # Multiply each CF column by the corresponding basis
      pv[,"premiums"] = pv[,"premiums"] * values$premiums[["gross"]];
      pv[,c("guaranteed", "survival", "death_SumInsured", "disease_SumInsured", "death_PremiumFree")] =
        pv[,c("guaranteed", "survival", "death_SumInsured", "disease_SumInsured", "death_PremiumFree")] * params$ContractData$sumInsured;
      pv[,c("death_GrossPremium", "death_Refund_past", "death_Refund_future")] = pv[,c("death_GrossPremium", "death_Refund_past", "death_Refund_future")] * values$premiums[["gross"]] * params$ContractData$premiumRefund;
      pv[,c("benefits", "additional_capital", "benefitsAndRefund", "alpha", "Zillmer", "beta", "gamma", "gamma_nopremiums", "unitcosts")] =
        pv[,c("benefits", "additional_capital", "benefitsAndRefund", "alpha", "Zillmer", "beta", "gamma", "gamma_nopremiums", "unitcosts")] * params$ContractData$sumInsured;

      # Sum all death-related payments to "death"  and remove the death_SumInsured column
      pv[,"death_SumInsured"] = pv[,"death_SumInsured"] + pv[,"death_GrossPremium"]
      colnames(pv)[colnames(pv) == "death_SumInsured"] = "death";

      pv = params$Hooks$Rounding$round("PV absolute", pv);

      cbind("premiums.unit" = values$presentValues[,"premiums"], pv)
    },


    #' @description Calculate the absolute present value time series of the
    #' benefits of the insurance contract
    #' @details Once the premiums of the insurance contracts are calculated, all
    #' present values can also be expressed in absolute terms. This function
    #' calculates these time series of the benefits present values in monetary
    #' terms, once the premiums and the unit-benefit present values are calculated
    #'  by the previous functions of this classe.
    #'
    #' This method is NOT to be called directly, but implicitly by the [InsuranceContract] object.
    presentValueBenefits = function(params, values) {
      if (getOption('LIC.debug.presentValueBenefits', FALSE)) {
        browser();
      }
      # TODO: Here we don't use the securityLoading parameter => Shall it be used or are these values to be understood without additional security loading?
      benefits    = values$presentValues[,"survival"] +
                    values$presentValues[,"guaranteed"] +
                    values$presentValues[,"death_SumInsured"] +
                    values$presentValues[,"disease_SumInsured"];
      allBenefits = benefits +
          values$presentValues[,"death_GrossPremium"] * values$premiums[["unit.gross"]] * params$ContractData$premiumRefund;
      benefitsCosts = rowSums( # Sum over the fourth dimension, leave the first three intact
        values$presentValuesCosts[,,"SumInsured",] +
        values$presentValuesCosts[,,"SumPremiums",] * values$unitPremiumSum * values$premiums[["unit.gross"]] +
        values$presentValuesCosts[,,"GrossPremium",] * values$premiums[["unit.gross"]] +
        values$presentValuesCosts[,,"NetPremium",] * values$premiums[["unit.net"]] +
        values$presentValuesCosts[,,"Constant",] / params$ContractData$sumInsured,
        dims = 2)

      rd = params$Hooks$Rounding;
      benefits = rd$round("PV abs benefits", benefits)
      allBenefits = rd$round("PV abs allBenefits", allBenefits)
      benefitsCosts = rd$round("PV abs benefitsCosts", benefitsCosts)

      cbind(
        benefits = benefits,
        benefitsAndRefund = allBenefits,
        benefitsCosts)
    },

    #' @description Calculate the linear coefficients of the premium calculation formula for the insurance contract
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    #' When \code{getPremiumCoefficients} is called, the \code{values$premiums}
    #' array has NOT yet been filled! Instead, all premiums already calculated
    #' (and required for the premium coefficients) are passed in the \code{premiums}
    #' argument.
    #'
    #' @param type The premium that is supposed to be calculated ("gross", "Zillmer", "net")
    #' @param coeffBenefits (empty) data structure of the benefit coefficients.
    #'        The actual values have no meaning, this parameter is only used to
    #'        derive the required dimensions
    #' @param coeffCosts (empty) data structure of the cost coefficients. The
    #'        actual values have no meaning, this parameter is only used to
    #'        derive the required dimensions
    #' @param premiums The premium components that have already been calculated
    #'         (e.g. for net and Zillmer, the gross premium has already been
    #'         calculated to allow modelling the premium refund)
    getPremiumCoefficients = function(type = "gross", coeffBenefits, coeffCosts, premiums, params, values, premiumCalculationTime = values$int$premiumCalculationTime) {
      if (getOption('LIC.debug.getPremiumCoefficients', FALSE)) {
        browser();
      }
      # Merge a possibly passed loadings override with the defaults of this class:
      securityLoading = valueOrFunction(params$Loadings$security, params = params, values = values);
      t = as.character(premiumCalculationTime)

      coeff = list(
        "SumInsured" = list("benefits" = coeffBenefits*0, "costs" = coeffCosts*0),
        "Premium"    = list("benefits" = coeffBenefits*0, "costs" = coeffCosts*0)
      );

      coeff[["Premium"]][["benefits"]][["premiums"]]            = 1;
      if (!is.null(params$ContractData$sumInsured)) {
          coeff[["SumInsured"]][["benefits"]][["additional_capital"]]            = -1 / params$ContractData$sumInsured;
      } else {
          coeff[["SumInsured"]][["benefits"]][["additional_capital"]]            = 0;
      }

      # Costs proportional to NetPremium introduce a non-linearity, as the NP is not available when the gross premium is calculated
      # => the corresponding costs PV is included in the coefficient!
      coeff.benefits = (1 + securityLoading);
      if (type == "gross") {
          # TODO: How to include this into the Zillmer premium calculation?
          if (values$presentValues[[t,"premiums"]] != 0) {
            coeff.benefits = coeff.benefits * (1 + sum(values$presentValuesCosts[t, c("alpha", "beta", "gamma"), "NetPremium",]) / values$presentValues[[t,"premiums"]])
          }
      }
      coeff[["SumInsured"]][["benefits"]][["guaranteed"]]       = coeff.benefits;
      coeff[["SumInsured"]][["benefits"]][["survival"]]         = coeff.benefits;
      coeff[["SumInsured"]][["benefits"]][["death_SumInsured"]] = coeff.benefits;
      coeff[["SumInsured"]][["benefits"]][["disease_SumInsured"]] = coeff.benefits;

      # Premium refund is handled differently for gross and net premiums, because it is proportional to the gross premium
      if (type == "gross") {
        coeff[["Premium"]][["benefits"]][["death_GrossPremium"]] = -params$ContractData$premiumRefund * coeff.benefits;
      } else if (type == "net" || type == "Zillmer") {
        coeff[["SumInsured"]][["benefits"]][["death_GrossPremium"]] = premiums[["unit.gross"]] * params$ContractData$premiumRefund * (1 + securityLoading);
      }


      # coefficients for the costs

      if (type == "gross") {
        affected = c("alpha", "beta", "gamma")
        if (params$Features$unitcostsInGross) {
          affected = c(affected, "unitcosts")
        }
        coeff[["SumInsured"]][["costs"]][affected, "SumInsured",  ] = 1;
        coeff[["SumInsured"]][["costs"]][affected, "Benefits",  ] = 1;
        # TODO: How to handle beta costs proportional to Sum Insured
        coeff[["Premium"]]   [["costs"]][affected, "SumPremiums", ] = -values$unitPremiumSum;
        coeff[["Premium"]]   [["costs"]][affected, "GrossPremium",] = -1;
        if (!is.null(params$ContractData$sumInsured)) {
            coeff[["SumInsured"]][["costs"]][affected, "Constant",    ] = 1 / params$ContractData$sumInsured;
        } else {
            coeff[["SumInsured"]][["costs"]][affected, "Constant",    ] = 0;
        }
      } else if (type == "Zillmer") {
          # TODO: Include costs with basis NetPremium and fixed costs!
        affected = c("Zillmer")
        if (params$Features$betaGammaInZillmer) {
          affected = c(affected, "beta", "gamma")
        }
        if (params$Features$gammaInZillmer) {
          affected = c(affected, "gamma")
        }
        coeff[["SumInsured"]][["costs"]][affected,"SumInsured",  ] = 1;
        coeff[["SumInsured"]][["costs"]][affected,"SumPremiums", ] = values$unitPremiumSum * premiums[["unit.gross"]];
        coeff[["SumInsured"]][["costs"]][affected,"GrossPremium",] = premiums[["unit.gross"]];
      }

      applyHook(params$Hooks$adjustPremiumCoefficients, coeff, type = type, premiums = premiums, params = params, values = values, premiumCalculationTime = premiumCalculationTime)
    },

    #' @description Calculate the sumInsured of the InsuranceContract given the
    #' parameters and premiums given and teh , present values already calculated and
    #' stored in the \code{params} and \code{values} lists.
    #' @param calculationTime the time when the sumInsured should be recalculated from the given premium
    #'
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    sumInsuredCalculation = function(params, values, calculationTime = values$int$premiumCalculationTime) {
      if (getOption('sumInsuredCalculation', FALSE)) {
        browser();
      }
      loadings = params$Loadings;
      values$premiums = c(
        "unit.net" = 0, "unit.Zillmer" = 0, "unit.gross" = 0,
        "net" = 0, "Zillmer" = 0, "gross" = 0,
        "unitcost" = 0, "written_yearly" = 0,
        "written_beforetax" = 0, "tax" = 0, "written" = 0, "additional_capital" = 0);
      coefficients = list("gross" = c(), "Zillmer" = c(), "net" = c());

      # Get the present values of the premiums, claims and costs at time 'calculationTime' (where the premium is to be calculated)
      t = as.character(calculationTime)
      pv = values$presentValues[t,]
      pvCost = values$presentValuesCosts[t,,,]

      #======================================================================== =
      # Calculate sumInsured from Premium, if needed
      # ======================================================================= =
      sumInsured = 1
      params$ContractData$sumInsured = 1 # Temporarily set to 1!

      # Premium type can be given using a named array, e.g. premium = c(gross = 1000)
      premiumtype = names(params$ContractData$premium)
      if (is.null(premiumtype)) premiumtype = "written";
      premium = unname(params$ContractData$premium);
      # if ((premium == 0) || (params$ContractData$premiumPeriod == 0)) premiumtype = "noPremium";
      calculating = FALSE;


      # Calculate unit gross premium (sumInsured=1)
      values$premiums["additional_capital"] = values$cashFlows[t, "additional_capital"]
      coeff = self$getPremiumCoefficients("gross", pv * 0, pvCost * 0, premiums = values$premiums, params = params, values = values, premiumCalculationTime = calculationTime)
      #### BEGIN TWEAK_RK 21.1.2023
      # Since we don't know the sumInsured, the initial capital cannot be scaled to sumInsured=1 => handle differently!
      # Calculate pv of benefits without initialCapital:
      coeff.initialCapital = coeff[["SumInsured"]][["benefits"]][["additional_capital"]]
      coeff[["SumInsured"]][["benefits"]][["additional_capital"]] = 0
      #### END TWEAK_RK 21.1.2023

      enumerator  = sum(coeff[["SumInsured"]][["benefits"]] * pv) + sum(coeff[["SumInsured"]][["costs"]] * pvCost);
      denominator = sum(coeff[["Premium"   ]][["benefits"]] * pv) + sum(coeff[["Premium"   ]][["costs"]] * pvCost);
      if (is.na(denominator) || (denominator == 0)) {
        values$premiums[["unit.gross"]] = 0;
        denominator = 1
      } else {
        values$premiums[["unit.gross"]] = enumerator/denominator * (1 + loadings$ongoingAlphaGrossPremium);
      }

      # Calculate other premium components:
      # ATTENTION: This will not work if any of these depend on the absolute values of the premiums, or depend on net or Zillmer premium!
      tax           = valueOrFunction(loadings$tax,          params = params, values = values);
      unitCosts     = valueOrFunction(loadings$unitcosts,    params = params, values = values);
      noMedicalExam = valueOrFunction(loadings$noMedicalExam,params = params, values = values);
      noMedicalExam.relative = valueOrFunction(loadings$noMedicalExamRelative,params = params, values = values);
      extraRebate   = valueOrFunction(loadings$extraRebate,  params = params, values = values);
      sumRebate     = valueOrFunction(loadings$sumRebate,    params = params, values = values);
      premiumRebateRate = valueOrFunction(loadings$premiumRebate,params = params, values = values);
      premiumRebate = applyHook(params$Hooks$premiumRebateCalculation, premiumRebateRate, params = params, values = values);

      extraChargeGrossPremium = valueOrFunction(loadings$extraChargeGrossPremium, params = params, values = values);
      advanceProfitParticipation = 0;
      advanceProfitParticipationUnitCosts = 0;
      ppScheme      = params$ProfitParticipation$profitParticipationScheme;
      if (!is.null(ppScheme)) {
        advanceProfitParticipation = ppScheme$getAdvanceProfitParticipation(params = params, values = values)
        advanceProfitParticipationUnitCosts = ppScheme$getAdvanceProfitParticipationAfterUnitCosts(params = params, values = values)
      }
      if (is.null(advanceProfitParticipation)) advanceProfitParticipation = 0;
      if (is.null(advanceProfitParticipationUnitCosts)) advanceProfitParticipationUnitCosts = 0;

      partnerRebate = valueOrFunction(loadings$partnerRebate, params = params, values = values);

      # Start from the given premium to derive the sumInsured step-by-step:
      temp = premium
      #
      # Written premium after tax
      calculating = calculating | (premiumtype == "written");
      if (calculating) {
        temp = temp / (1 + tax);
      }
      # Written premium before tax
      calculating = calculating | (premiumtype == "written_beforetax");
      if (calculating) {
        temp = temp / (1 - premiumRebate - advanceProfitParticipationUnitCosts - partnerRebate);

        pv.unitcosts = sum(
          pvCost["unitcosts","SumInsured",] * sumInsured +
            pvCost["unitcosts","SumPremiums",] * values$unitPremiumSum * values$premiums[["gross"]] +
            pvCost["unitcosts","GrossPremium",] * values$premiums[["gross"]] +
            pvCost["unitcosts","NetPremium",] * values$premiums[["net"]] +
            pvCost["unitcosts","Constant",]
        )
        if (pv[["premiums"]] == 0) {
          premium.unitcosts = 0
        } else {
          premium.unitcosts = pv.unitcosts / pv[["premiums"]] + valueOrFunction(loadings$unitcosts, params = params, values = values);
        }
        if (!params$Features$unitcostsInGross) {
          temp = temp - premium.unitcosts;
        }
        temp = temp / (1 - advanceProfitParticipation)
      }
      calculating = calculating | (premiumtype == "gross");
      if (calculating) {
        # handle initialCapital here (coefficient for initialCapital is typically negative, as it reduces the premium!)
        temp = temp - coeff.initialCapital * values$premiums[["additional_capital"]] / denominator;
        temp = temp /
          (enumerator / denominator * (1 + noMedicalExam.relative + extraChargeGrossPremium) + noMedicalExam - sumRebate - extraRebate);
      }
      sumInsured = params$Hooks$Rounding$round("sumInsured", temp);

      sumInsured
    },

    #' @description Calculate the premiums of the InsuranceContract given the
    #' parameters, present values and premium cofficients already calculated and
    #' stored in the \code{params} and \code{values} lists.
    #'
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    premiumCalculation = function(params, values, premiumCalculationTime = values$int$premiumCalculationTime) {
      if (getOption('LIC.debug.premiumCalculation', FALSE)) {
        browser();
      }
      rd = params$Hooks$Rounding;

      loadings = params$Loadings;
      sumInsured = params$ContractData$sumInsured
      values$premiums = c(
        "unit.net" = 0, "unit.Zillmer" = 0, "unit.gross" = 0,
        "net" = 0, "Zillmer" = 0, "gross" = 0,
        "unitcost" = 0, "written_yearly" = 0,
        "written_beforetax" = 0, "tax" = 0, "written" = 0, "additional_capital" = 0);
      coefficients = list("gross" = c(), "Zillmer" = c(), "net" = c());

      # Get the present values of the premiums, claims and costs at time 'premiumCalculationTime' (where the premium is to be calculated)
      t = as.character(premiumCalculationTime)
      pv = values$presentValues[t,]
      pvCost = values$presentValuesCosts[t,,,]

      values$premiums["additional_capital"] = values$cashFlows[t, "additional_capital"]

      # If there are no premium payments, no need to calculate premium components
      # if (pv[["premiums"]] == 0) {
      #   return(list("premiums" = values$premiums, "coefficients" = coefficients, "sumInsured" = params$ContractData$sumInsured))
      # }



      #======================================================================== =
      # net, gross and Zillmer premiums are calculated from the present values using the coefficients on each present value as described in the formulas document
      # ======================================================================= =
      # GROSS Premium
      # ----------------------------------------------------------------------- -
      coeff = self$getPremiumCoefficients("gross", pv * 0, pvCost * 0, premiums = values$premiums, params = params, values = values, premiumCalculationTime = premiumCalculationTime)
      enumerator  = sum(coeff[["SumInsured"]][["benefits"]] * pv) + sum(coeff[["SumInsured"]][["costs"]] * pvCost);
      denominator = sum(coeff[["Premium"   ]][["benefits"]] * pv) + sum(coeff[["Premium"   ]][["costs"]] * pvCost);
      values$premiums[["unit.gross"]] = rd$round("Premium gross unit", enumerator/ifelse(denominator == 0, 1, denominator) * (1 + loadings$ongoingAlphaGrossPremium));
      values$premiums[["gross"]] = rd$round("Premium gross", values$premiums[["unit.gross"]] * sumInsured);
      coefficients[["gross"]] = coeff;

      # ======================================================================= =
      # NET Premium
      # ----------------------------------------------------------------------- -
      coeff = self$getPremiumCoefficients("net", pv*0, pvCost*0, premiums = values$premiums, params = params, values = values, premiumCalculationTime = premiumCalculationTime)
      enumerator  = sum(coeff[["SumInsured"]][["benefits"]] * pv) + sum(coeff[["SumInsured"]][["costs"]] * pvCost);
      denominator = sum(coeff[["Premium"   ]][["benefits"]] * pv) + sum(coeff[["Premium"   ]][["costs"]] * pvCost);
      values$premiums[["unit.net"]] = rd$round("Premium net unit", enumerator/ifelse(denominator == 0, 1, denominator));
      values$premiums[["net"]] = rd$round("Premium net", values$premiums[["unit.net"]] * sumInsured);
      coefficients[["net"]] = coeff;

      # ======================================================================= =
      # ZILLMER Premium
      # ----------------------------------------------------------------------- -
      coeff = self$getPremiumCoefficients("Zillmer", pv * 0, pvCost * 0, premiums = values$premiums, params = params, values = values, premiumCalculationTime = premiumCalculationTime);
      enumerator  = sum(coeff[["SumInsured"]][["benefits"]] * pv) + sum(coeff[["SumInsured"]][["costs"]] * pvCost);
      denominator = sum(coeff[["Premium"   ]][["benefits"]] * pv) + sum(coeff[["Premium"   ]][["costs"]] * pvCost);
      values$premiums[["unit.Zillmer"]] = rd$round("Premium Zillmer unit", enumerator/ifelse(denominator == 0, 1, denominator));
      values$premiums[["Zillmer"]] = rd$round("Premium Zillmer", values$premiums[["unit.Zillmer"]] * sumInsured);
      coefficients[["Zillmer"]] = coeff;


      # ======================================================================= =
      # Additional premium components (after gross premium)
      # ----------------------------------------------------------------------- -
      # The written premium is the gross premium with additional loadings, rebates, unit costs and taxes
      # TODO: Think through, how each of the components should / could be rounded
      tax           = valueOrFunction(loadings$tax,          params = params, values = values);
      unitCosts     = valueOrFunction(loadings$unitcosts,    params = params, values = values);
      noMedicalExam = valueOrFunction(loadings$noMedicalExam,params = params, values = values);
      noMedicalExam.relative = valueOrFunction(loadings$noMedicalExamRelative,params = params, values = values);
      extraRebate   = valueOrFunction(loadings$extraRebate,  params = params, values = values);
      sumRebate     = valueOrFunction(loadings$sumRebate,    params = params, values = values);
      premiumRebateRate = valueOrFunction(loadings$premiumRebate,params = params, values = values);
      premiumRebate = applyHook(params$Hooks$premiumRebateCalculation, premiumRebateRate, params = params, values = values);

      extraChargeGrossPremium = valueOrFunction(loadings$extraChargeGrossPremium, params = params, values = values);
      advanceProfitParticipation = 0;
      advanceProfitParticipationUnitCosts = 0;
      ppScheme      = params$ProfitParticipation$profitParticipationScheme;
      if (!is.null(ppScheme)) {
          advanceProfitParticipation = ppScheme$getAdvanceProfitParticipation(params = params, values = values)
          advanceProfitParticipationUnitCosts = ppScheme$getAdvanceProfitParticipationAfterUnitCosts(params = params, values = values)
      }
      if (is.null(advanceProfitParticipation)) advanceProfitParticipation = 0;
      if (is.null(advanceProfitParticipationUnitCosts)) advanceProfitParticipationUnitCosts = 0;

      partnerRebate = valueOrFunction(loadings$partnerRebate, params = params, values = values);

      pv.unitcosts = sum(
        pvCost["unitcosts","SumInsured",] * sumInsured +
        pvCost["unitcosts","SumPremiums",] * values$unitPremiumSum * values$premiums[["gross"]] +
        pvCost["unitcosts","GrossPremium",] * values$premiums[["gross"]] +
        pvCost["unitcosts","NetPremium",] * values$premiums[["net"]] +
        pvCost["unitcosts","Constant",]
      )
      premium.unitcosts = ifelse(pv[["premiums"]] == 0, 0, pv.unitcosts / pv[["premiums"]] + valueOrFunction(loadings$unitcosts, params = params, values = values));


      frequencyLoading = self$evaluateFrequencyLoading(loadings$premiumFrequencyLoading, params$ContractData$premiumFrequency, params = params, values = values)
      premiumBeforeTax = (values$premiums[["unit.gross"]]*(1 + noMedicalExam.relative + extraChargeGrossPremium) + noMedicalExam - sumRebate - extraRebate) * sumInsured * (1 - advanceProfitParticipation);
      if (!params$Features$unitcostsInGross) {
        premiumBeforeTax = premiumBeforeTax + premium.unitcosts;
      }
      premiumBeforeTax = premiumBeforeTax * (1 - premiumRebate - advanceProfitParticipationUnitCosts - partnerRebate);
      premiumBeforeTax.y = premiumBeforeTax * (1 + frequencyLoading);
      premiumBeforeTax = premiumBeforeTax.y / params$ContractData$premiumFrequency;
      values$premiums[["written_yearly"]] = rd$round("Premium written yearly", premiumBeforeTax.y * (1 + tax))
      values$premiums[["written_beforetax"]] = rd$round("Premium written beforeTax", premiumBeforeTax);
      values$premiums[["tax"]] = rd$round("Premium tax", premiumBeforeTax * tax);
      values$premiums[["written"]] = rd$round("Premium written", values$premiums[["written_beforetax"]] + values$premiums[["tax"]]);

      applyHook(
        params$Hooks$adjustPremiums,
        list("premiums" = values$premiums, "coefficients" = coefficients, "sumInsured" = params$ContractData$sumInsured),
        params = params, values = values
      )
    },

    #' @description Calculate the reserves of the InsuranceContract given the
    #' parameters, present values and premiums already calculated and stored in
    #' the \code{params} and \code{values} lists.
    #'
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    reserveCalculation = function(params, values) {
      if (getOption('LIC.debug.reserveCalculation', FALSE)) {
        browser();
      }
      rd = params$Hooks$Rounding;
      t = "0"
      securityFactor = (1 + valueOrFunction(params$Loadings$security, params = params, values = values));
      ppScheme      = params$ProfitParticipation$profitParticipationScheme;

      absPV = applyHook(params$Hooks$adjustPVForReserves, values$absPresentValues, params = params, values = values);

      # Net, Zillmer and Gross reserves
      resNet = rd$round("Res net", absPV[,"benefitsAndRefund"] * securityFactor - values$premiums[["net"]] * absPV[,"premiums.unit"]);
      BWZcorr = ifelse(absPV[t, "premiums"] == 0, 0,
                       absPV[t, "Zillmer"] / absPV[t, "premiums"]) * absPV[,"premiums"];
      resZ = rd$round("Res Zillmer", resNet - BWZcorr);

      resAdeq = absPV[,"benefitsAndRefund"] * securityFactor +
          absPV[,"alpha"] + absPV[,"beta"] + absPV[,"gamma"] -
        values$premiums[["gross"]] * absPV[,"premiums.unit"];
      if (params$Features$unitcostsInGross) {
          resAdeq = resAdeq + absPV[, "unitcosts"]
      }
      resAdeq = rd$round("Res adequate", resAdeq)

      resGamma = rd$round("Res gamma", absPV[,"gamma"] -
        ifelse(absPV[t, "premiums"] == 0, 0,
               absPV[t, "gamma"] / absPV[t, "premiums"]) * absPV[,"premiums"]);


      advanceProfitParticipation = 0;
      if (!is.null(ppScheme)) {
          advanceProfitParticipation = ppScheme$getAdvanceProfitParticipation(params = params, values = values)
      }

      # Alpha refund: Distribute alpha-costs to 5 years (or if shorter, the policy period), always starting at time 0:
      # If alphaRefunded==TRUE, don't refund a second time!
      if (params$ContractState$alphaRefunded) {
        alphaRefund = 0
      } else {
        r = min(params$ContractData$policyPeriod, params$Loadings$alphaRefundPeriod);
        ZillmerSoFar = Reduce("+", values$absCashFlows$Zillmer, accumulate = TRUE);
        ZillmerTotal = sum(values$absCashFlows$Zillmer);
        len = length(ZillmerSoFar);
        if (params$Features$alphaRefundLinear) {
          ZillmerVerteilungCoeff = pad0((0:r)/r, len, 1);
        } else {
          q = self$getTransitionProbabilities(params, values);
          pvf = PVfactory$new(qx = q, v = 1/(1 + params$ActuarialBases$i))
          # vector of all ä_{x+t, r-t}
          pvAlphaTmp = pvf$survival(advance = pad0(rep(1,r), len))
          ZillmerVerteilungCoeff = (1 - pvAlphaTmp/pvAlphaTmp[[1]]);
        }
        alphaRefund = ZillmerSoFar - ZillmerVerteilungCoeff * ZillmerTotal;
      }

      # Reduction Reserve: Reserve used for contract modifications:
      if (params$Features$zillmering) {
        resContractual = resZ + resGamma
        resReduction = resZ + alphaRefund;
      } else {
        resContractual = resAdeq + resGamma
        resReduction = resAdeq + alphaRefund;
      }
      resContractual = rd$round("Res contractual", resContractual)
      resReduction = rd$round("Res reduction", resReduction)

      resConversion = rd$round("Res conversion", resContractual * (1 - advanceProfitParticipation));
      if (params$Features$surrenderIncludesCostsReserves) {
        resReduction = resReduction + resGamma;
      }
      resReduction = pmax(0, resReduction) # V_{x,n}^{Rkf}

      # Collect all reserves to one large matrix
      res = cbind(
            "SumInsured"  = head0(rep(params$ContractData$sumInsured, values$int$l)),
            "net"         = resNet,
            "Zillmer"     = resZ,
            "adequate"    = resAdeq,
            "gamma"       = resGamma,
            "contractual" = resContractual,
            "conversion"  = resConversion,
            "alphaRefund" = alphaRefund,
            "reduction"   = resReduction
            #, "Reserve.premiumfree"=res.premiumfree, "Reserve.gamma.premiumfree"=res.gamma.premiumfree);
      );
      rownames(res) <- rownames(absPV);
      values$reserves = res;

      # The surrender value functions can have arbitrary form, so we store a function
      # here in the tarif and call that, passing the reduction reserve as
      # starting point, but also all reserves, cash flows, premiums and present values
      if (!params$ContractState$surrenderPenalty) {
          # No surrender penalty any more (has already been applied to the first contract change!)
          surrenderValue = resReduction;
      } else if (!is.null(params$ActuarialBases$surrenderValueCalculation)) {
          surrenderValue = params$ActuarialBases$surrenderValueCalculation(resReduction, params, values);
      } else {
          # by default, refund the full reduction reserve, except the advance profit participation, which is also included in the reserve, but not charged on the premium!
          advanceProfitParticipationUnitCosts = 0;
          ppScheme      = params$ProfitParticipation$profitParticipationScheme;
          if (!is.null(ppScheme)) {
              advanceProfitParticipationUnitCosts = ppScheme$getAdvanceProfitParticipationAfterUnitCosts(params = params, values = values)
          }
          partnerRebate = valueOrFunction(params$Loadings$partnerRebate, params = params, values = values);
          surrenderValue = resReduction * (1 - advanceProfitParticipationUnitCosts - partnerRebate);
      }
      surrenderValue = rd$round("Surrender Value", surrenderValue)


      # Calculate new sum insured after premium waiver
      if (!is.null(params$ActuarialBases$premiumWaiverValueCalculation)) {
        premiumfreeValue = params$ActuarialBases$premiumWaiverValueCalculation(resReduction, params, values);
      } else {
        premiumfreeValue = surrenderValue
      }
      Storno = 0; # TODO: Implement storno costs
      premiumfreePV = (absPV[, "benefits"] * securityFactor + absPV[, "gamma_nopremiums"]); # PV of future premium free claims + costs
      newSI = ifelse(premiumfreePV == 0, 0,
        (premiumfreeValue - absPV[,"death_Refund_past"] * securityFactor - c(Storno)) /
        premiumfreePV * params$ContractData$sumInsured);
      newSI = rd$round("Premiumfree SI", newSI)

      cbind(res,
            "PremiumsPaid" = Reduce("+", values$absCashFlows$premiums_advance, accumulate = TRUE),
            "Surrender" = surrenderValue,
            "PremiumFreeSumInsured" = newSI
      )
    },

    #' @description Calculate the (linear) interpolation factors for the balance
    #' sheet reserve (Dec. 31) between the yearly contract closing dates
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    #' @param method The method for the balance sheet interpolation (30/360, act/act, act/360, act/365 or a function)
    #' @param years how many years to calculate (for some usances, the factor
    #'      is different in leap years!)
    getBalanceSheetReserveFactor = function(method, params, years = 1) {
      if (getOption('LIC.debug.getBalanceSheetReserveFactor', FALSE)) {
        browser();
      }
      balanceDate = params$ActuarialBases$balanceSheetDate
      year(balanceDate) = year(params$ContractData$contractClosing);
      if (balanceDate < params$ContractData$contractClosing) {
        balanceDate = balanceDate + years(1);
      }

      # contractDates = params$ContractData$contractClosing + years(1:years);
      # balanceDates = balanceDate + years(1:years - 1);
      contractDates = seq(params$ContractData$contractClosing, length.out = years, by = "year")
      balanceDates = seq(balanceDate, length.out = years, by = "year")

      if (is.function(method)) {
        baf = method(params = params, contractDates = contractDates, balanceDates = balanceDates)
      } else if (method == "30/360") {
        baf = ((month(balanceDates + days(1)) - month(contractDates) - 1) %% 12 + 1) / 12
      } else if (method == "act/act") {
        baf = as.numeric((balanceDates + days(1)) - contractDates, units = "days" ) / as.numeric(balanceDates - (balanceDates - years(1)), units = "days")
      } else if (method == "act/360") {
        baf = pmin(as.numeric((balanceDates + days(1)) - contractDates, units = "days" ) / 360, 1)
      } else if (method == "act/365") {
        baf = pmin(as.numeric((balanceDates + days(1)) - contractDates, units = "days" ) / 365, 1)
      }
      data.frame(date = balanceDates, time = baf + (1:years) - 1, baf = baf)
    },

    #' @description Calculate the reserves for the balance sheets at Dec. 31 of each
    #'              year by interpolation from the contract values calculated for
    #'              the yearly reference date of the contract
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    reserveCalculationBalanceSheet = function(params, values) {
      if (getOption('LIC.debug.reserveCalculationBalanceSheet', FALSE)) {
        browser();
      }
      rd = params$Hooks$Rounding
      reserves = values$reserves;
      years = length(reserves[,"Zillmer"]);
      # Balance sheet reserves:
      factors = self$getBalanceSheetReserveFactor(method = params$ActuarialBases$balanceSheetMethod, params = params, years = years);
      baf = factors$baf
      factors$baf = NULL

      useUnearnedPremiums = valueOrFunction(params$Features$useUnearnedPremiums, params = params, values = values)
      resN_BS = (1 - baf) * (reserves[,"net"] + if (!useUnearnedPremiums) values$premiumComposition[,"net"] else 0) + baf * c(reserves[-1, "net"], 0)
      resZ_BS = (1 - baf) * (reserves[,"Zillmer"] + if (!useUnearnedPremiums) values$premiumComposition[,"Zillmer"] else 0) + baf * c(reserves[-1, "Zillmer"], 0)
      resGamma_BS = (1 - baf) * (reserves[,"gamma"] + if (!useUnearnedPremiums) values$premiumComposition[,"gamma"] else 0) + baf * c(reserves[-1, "gamma"], 0)
      resGross_BS = (1 - baf) * (reserves[,"adequate"] + if (!useUnearnedPremiums) sum(values$premiumComposition[,c("alpha", "beta", "gamma")]) else 0) + baf * c(reserves[-1, "adequate"], 0)
      if (params$Features$zillmering) {
        res_BS = resZ_BS + resGamma_BS;
      } else {
        res_BS = resGross_BS;
      }

      # Premium transfer / unearned premium:
      if (useUnearnedPremiums) {
        fact = valueOrFunction(params$ActuarialBases$unearnedPremiumsMethod, params = params, dates = factors$date)
        if (is.null(fact) || is.na(fact)) {
          freq = params$ContractData$premiumFrequency
          bm = month(params$ContractData$contractClosing)

                    fact = (bm - month(factors$date) + 12 - 1) %% (12/freq) * (freq/12)
        }
        # TODO: We have no vector of actual written premiums (implicit assumption
        # seems to be that the premium stays constant!). Once we have such a vector,
        # rewrite the following code
        unearnedPremiums = fact * values$cashFlows$premiums_advance * values$premiums[["written_beforetax"]] # TODO
        # If advance profit participation is granted, unearned premiums still apply to the whole gross premium without PP and partner rebate!
        ppScheme      = params$ProfitParticipation$profitParticipationScheme;
        if (!is.null(ppScheme)) {
          partnerRebate = valueOrFunction(params$Loadings$partnerRebate, params = params, values = values);
          advanceProfitParticipation = ppScheme$getAdvanceProfitParticipationAfterUnitCosts(params = params, values = values);
          unearnedPremiums = unearnedPremiums / (1 - partnerRebate - advanceProfitParticipation);
        }
      } else {
        # If reserves contain single-premium, no unearned premiums are shown in the balance sheet!
        unearnedPremiums = 0
      }

      # Collect all reserves to one large matrix
      res = cbind(factors,
                  "net"                   = pmax(resN_BS,0),
                  "Zillmer"               = pmax(resZ_BS,0),
                  "gamma"                 = pmax(resGamma_BS,0),
                  "gross"                 = pmax(resGross_BS,0),
                  "Balance Sheet Reserve" = pmax(res_BS,0),
                  "unearned Premiums"     = unearnedPremiums
      );
      rownames(res) <- rownames(reserves);
      rd$round("Balance Sheet", res)
    },

    #' @description Calculate the profit participation given the contract
    #' parameters and the already calculated reserves of the contract.
    #'
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    #' @param ... Additional parameters for the profit participation calculation, passed
    #'            through to the profit participation scheme's \href{../../LifeInsureR/html/ProfitParticipation.html#method-getProfitParticipation}{\code{ProfitParticipation$getProfitParticipation()}}
    calculateProfitParticipation = function(params, ...) {
      if (getOption('LIC.debug.calculateProfitParticipation', FALSE)) {
        browser();
      }
      ppScheme = params$ProfitParticipation$profitParticipationScheme;
        if (!is.null(ppScheme)) {
            ppScheme$getProfitParticipation(params = params, ...)
        }
    },

    #' @description Calculate the reserves after profit participation for the given profit scenario
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    #' @param profitScenario The ID of the profit scenario for which to calculate the reserves
    #' @param ... TODO
    reservesAfterProfit = function(profitScenario, params, values, ...) {
      if (getOption('LIC.debug.reservesAfterProfit', FALSE)) {
        browser();
      }
      # TODO
    },


    #' @description Return the time series of the basic contract
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    getBasicDataTimeseries = function(params, values) {
      if (getOption('LIC.debug.getBasicDataTimeseries', FALSE)) {
        browser();
      }
      res = cbind(
            "PremiumPayment" = values$premiumComposition[, "charged"] > 0,
            "SumInsured" = values$reserves[, "SumInsured"],
            "Premiums" = values$absCashFlows$premiums_advance + values$absCashFlows$premiums_arrears,
            "InterestRate" = rep(params$ActuarialBases$i, values$int$l),
            "PolicyDuration" = rep(values$int$policyTerm, values$int$l),
            "PremiumPeriod" = rep(values$int$premiumTerm, values$int$l)
        );
        rownames(res) = 0:(values$int$l-1);
        res
    },

    #' @description Calculate the time series of the premium decomposition of the contract
    #' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
    #'          All premiums, reserves and present values have already been calculated.
    premiumDecomposition = function(params, values) {
      # TODO: No rounding applied yet!
      if (getOption('LIC.debug.premiumDecomposition', FALSE)) {
        browser();
      }
      loadings   = params$Loadings;
      sumInsured = params$ContractData$sumInsured;
      premiums   = values$premiums;
      v          = 1/(1 + params$ActuarialBases$i);
      l          = dim(values$reserves)[[1]];
      ppScheme   = params$ProfitParticipation$profitParticipationScheme;
      t          = as.character(0) # Time for original premium calculation => TODO: Use values stored in ContractData?

      # TODO: This assumes all premiums are paid in advance!
      premium.gross    = values$absCashFlows[,"premiums_advance"];

      # First get the charges and rebates that are added to the gross premium to obtain the charged premium:

      # charge for no medical exam:
      extraChargeGrossPremium = valueOrFunction(loadings$extraChargeGrossPremium, params = params, values = values);
      noMedExam        = valueOrFunction(loadings$noMedicalExam,params = params, values = values);
      noMedExam.rel    = valueOrFunction(loadings$noMedicalExamRelative,params = params, values = values);
      withMedExam      = premium.gross * (1 + noMedExam.rel + extraChargeGrossPremium) + noMedExam * sumInsured;
      charge.noMedicalExam = withMedExam - premium.gross;

      # sum rebate:
      sumRebate        = valueOrFunction(loadings$sumRebate,    params = params, values = values);
      extraRebate      = valueOrFunction(loadings$extraRebate,  params = params, values = values);
      afterSumRebate   = withMedExam - (sumRebate + extraRebate) * sumInsured; # calculate the charge as the difference, because we want a vector!
      rebate.sum       = afterSumRebate - withMedExam;

      # advance profit participation has two parts, one before and one after unit costs. Part 1:
      advanceProfitParticipation = 0;
      if (!is.null(ppScheme)) {
          advanceProfitParticipation = ppScheme$getAdvanceProfitParticipation(params = params, values = values)
      }
      afterProfit      = afterSumRebate * (1 - advanceProfitParticipation);
      profits.advance  = afterProfit - afterSumRebate;

      # unit costs
      unitCosts        = premiums[["unitcost"]];
      # unit costs are only charged if a premium is paid, so exclude all times with premium==0!
      if (!params$Features$unitcostsInGross) {
          afterUnitCosts   = afterProfit + (afterProfit != 0)*unitCosts;
          unitcosts        = afterUnitCosts - afterProfit;
      } else {
          afterUnitCosts   = afterProfit;
          unitcosts        = 0;
      }

      # advance profit participation, Part 2:
      advanceProfitParticipationUnitCosts = 0;
      if (!is.null(ppScheme)) {
          advanceProfitParticipationUnitCosts = ppScheme$getAdvanceProfitParticipationAfterUnitCosts(params = params, values = values)
      }
      afterProfit      = afterUnitCosts * (1 - advanceProfitParticipationUnitCosts);
      profits.advance  = profits.advance + afterProfit - afterUnitCosts;

      # premium rebate
      premiumRebateRate = valueOrFunction(loadings$premiumRebate,params = params, values = values);
      premiumRebate = applyHook(params$Hooks$premiumRebateCalculation, premiumRebateRate, params = params, values = values);

      afterPremiumRebate = afterUnitCosts * (1 - premiumRebate);
      rebate.premium   = afterPremiumRebate - afterUnitCosts;

      # partner rebate
      partnerRebate    = valueOrFunction(loadings$partnerRebate,params = params, values = values);
      afterPartnerRebate = afterUnitCosts * (1 - partnerRebate);
      rebate.partner   = afterPartnerRebate - afterUnitCosts;

      # value after all rebates
      afterRebates     = afterProfit + rebate.premium + rebate.partner;

      # premium frequency loading
            frequencyLoading = self$evaluateFrequencyLoading(loadings$premiumFrequencyLoading, params$ContractData$premiumFrequency, params = params, values = values)

            afterFrequency   = afterRebates * (1 + frequencyLoading);
      charge.frequency = afterFrequency - afterRebates;

      # insurance tax
      taxRate          = valueOrFunction(loadings$tax,          params = params, values = values);
      premium.charged  = afterFrequency * (1 + taxRate);
      tax              = premium.charged - afterFrequency;


      # Gross premium = net + zillmeredAlpha + unzillmeredAlpha + beta + gamma premium
      unit.premiumCF   = if (premiums[["gross"]] == 0) { premium.gross * 0 } else { premium.gross / premiums[["gross"]] }
      if (values$absPresentValues[t, "premiums.unit"] == 0) {
        premium.gamma    = 0
        premium.beta     = 0
        premium.alpha    = 0
        premium.alpha.Zillmer = 0
      } else {
        premium.gamma    = unit.premiumCF * values$absPresentValues[t, "gamma"] / values$absPresentValues[t, "premiums.unit"];
        premium.beta     = unit.premiumCF * values$absPresentValues[t, "beta"]  / values$absPresentValues[t, "premiums.unit"];
        premium.alpha    = unit.premiumCF * values$absPresentValues[t, "alpha"] / values$absPresentValues[t, "premiums.unit"];
        premium.alpha.Zillmer = unit.premiumCF * values$absPresentValues[t, "Zillmer"] / values$absPresentValues[t, "premiums.unit"];
      }
      premium.Zillmer  = unit.premiumCF * premiums[["Zillmer"]];
      premium.alpha.noZ = premium.alpha - premium.alpha.Zillmer; # ungezillmerter Teil der Abschlusskosten

      premium.net       = unit.premiumCF * premiums[["net"]];

      securityLoading   = valueOrFunction(params$Loadings$security, params = params, values = values);
      qq = self$getTransitionProbabilities(params, values);
      premium.risk.actual   = v * (values$absCashFlows[,"death"] - c(values$reserves[,"net"][-1], 0)) * pad0(qq$qx, l);
      premium.risk.security = v * (values$absCashFlows[,"death"] * securityLoading) * pad0(qq$qx, l);
      premium.risk          = premium.risk.actual + premium.risk.security;

      premium.risk.disease.actual   = v * (values$absCashFlows[,"disease_SumInsured"] - c(values$reserves[,"net"][-1], 0)) * pad0(qq$ix, l);
      premium.risk.disease.security = v * (values$absCashFlows[,"disease_SumInsured"] * securityLoading) * pad0(qq$ix, l);
      premium.risk.disease          = premium.risk.disease.actual + premium.risk.disease.security;
      premium.savings       = getSavingsPremium(
          values$reserves[,"net"], v = v,
          survival_advance = values$absCashFlows[,"survival_advance"] + values$absCashFlows[,"guaranteed_advance"],
          survival_arrears = values$absCashFlows[,"survival_arrears"] + values$absCashFlows[,"guaranteed_arrears"]
      )

      premium.Zillmer.risk.actual   = v * (values$absCashFlows[,"death"] - c(values$reserves[,"contractual"][-1], 0)) * pad0(qq$qx, l);
      premium.Zillmer.risk.security = v * (values$absCashFlows[,"death"] * securityLoading) * pad0(qq$qx, l);
      premium.Zillmer.risk          = premium.Zillmer.risk.actual + premium.Zillmer.risk.security;
      premium.Zillmer.risk.disease.actual   = v * (values$absCashFlows[,"disease_SumInsured"] - c(values$reserves[,"contractual"][-1], 0)) * pad0(qq$ix, l);
      premium.Zillmer.risk.disease.security = v * (values$absCashFlows[,"disease_SumInsured"] * securityLoading) * pad0(qq$ix, l);
      premium.Zillmer.risk.disease          = premium.Zillmer.risk.disease.actual + premium.Zillmer.risk.disease.security;


      premium.Zillmer.savings  = getSavingsPremium(
          values$reserves[,"contractual"], v = v,
          survival_advance = values$absCashFlows[,"survival_advance"] + values$absCashFlows[,"guaranteed_advance"],
          survival_arrears = values$absCashFlows[,"survival_arrears"] + values$absCashFlows[,"guaranteed_arrears"]
      )
      premium.Zillmer.amortization = getSavingsPremium(
              pmin(0, values$reserves[,"contractual"]), v = v
      )
      premium.Zillmer.actsavings = getSavingsPremium(
              pmax(0, values$reserves[,"contractual"]), v = v,
              survival_advance = values$absCashFlows[,"survival_advance"] + values$absCashFlows[,"guaranteed_advance"],
              survival_arrears = values$absCashFlows[,"survival_arrears"] + values$absCashFlows[,"guaranteed_arrears"]
      )

      res = cbind(
        "charged"         = premium.charged,
        "tax"             = tax,
        "loading.frequency" = charge.frequency,
        "rebate.premium"  = rebate.premium,
        "rebate.partner"  = rebate.partner,
        "unitcosts"       = unitcosts,
        "profit.advance"  = profits.advance,
        "rebate.sum"      = rebate.sum,
        "charge.noMedicalExam" = charge.noMedicalExam,
        "gross"           = premium.gross,

        "gamma"           = premium.gamma,
        "beta"            = premium.beta,
        "alpha"           = premium.alpha,
        "alpha.noZillmer" = premium.alpha.noZ,
        "alpha.Zillmer"   = premium.alpha.Zillmer,
        "Zillmer"         = premium.Zillmer,

        "net"                   = premium.net,

        "risk"                          = premium.risk,
        "premium.risk.actual"           = premium.risk.actual,
        "premium.risk.security"         = premium.risk.security,
        "risk.disease"                  = premium.risk.disease,
        "premium.risk.disease.actual"   = premium.risk.disease.actual,
        "premium.risk.disease.security" = premium.risk.disease.security,
        "savings"                       = premium.savings,

        "Zillmer.risk"                  =  premium.Zillmer.risk,
        "Zillmer.risk.actual"           = premium.Zillmer.risk.actual,
        "Zillmer.risk.security"         = premium.Zillmer.risk.security,
        "Zillmer.risk.disease"          = premium.Zillmer.risk.disease,
        "Zillmer.risk.disease.actual"   = premium.Zillmer.risk.disease.actual,
        "Zillmer.risk.disease.security" = premium.Zillmer.risk.disease.security,

        "Zillmer.savings"               = premium.Zillmer.savings,
        "Zillmer.amortization"          = premium.Zillmer.amortization,
        "Zillmer.savings.real"          = premium.Zillmer.actsavings
      )
      rownames(res) <- rownames(premiums);
      res
    },


    #' @description Generic function to calculate future sums of the values
    #' @param cf The time series, for which future sums at all times are desired
    #' @param ... currently unused
    calculateFutureSums = function(cf, ...) {
      rcumsum = function(vec) rev(cumsum(rev(vec)))
      apply(cf, 2, rcumsum)
    },
    #' @description Calculate all present values for a given time series. The
    #' mortalities are taken from the contract's parameters.
    #' @param cf The time series, for which future present values at all
    #'      times are desired
    #' @param ... currently unused
    calculatePresentValues = function(cf, params, values) {
      len = dim(cf)[1];
      qq = self$getTransitionProbabilities(params, values);
      i = params$ActuarialBases$i;
      premFreq = params$ContractData$premiumFrequency;
      premFreqCorr = correctionPaymentFrequency(
        i = i, m = premFreq,
        order = valueOrFunction(params$ActuarialBases$premiumFrequencyOrder, params = params, values = values));

      pvf = PVfactory$new(qx = qq, m = premFreq, mCorrection = premFreqCorr, v = 1/(1 + i));
      pvf$survival(advance = cf)
    },

    #' @description Calculate the premium frequency loading, i.e. the surcharge
    #' on the premium for those cases where the premium is not paid yearly.
    #' Return values can be either a numeric value or a named list with all
    #' possible premium frequencies as keys.
    #' @param loading The premiumFrequencyLoading parameter of the Contract or Tariff to be evaluated
    #' @param frequency The premiumFrequency parameter of the contract
    evaluateFrequencyLoading = function(loading, frequency, params, values) {
      frequencyLoading = valueOrFunction(loading, frequency = frequency, params = params, values = values);
      if (is.null(frequencyLoading)) {
        0
      } else if (is.list(frequencyLoading)) {
        if (as.character(frequency) %in% names(frequencyLoading)) {
          frequencyLoading[[as.character(frequency)]]
        } else {
          warning("Unable to handle premium frequency ", frequency, " with the given loading ", frequencyLoading);
        }
      } else if (is.numeric(frequencyLoading)) {
        frequencyLoading
      } else {
        warning("premiumFrequencyLoading must be a number or a named list, given: ", frequencyLoading);
        0
      }
    },



    #' @field dummy Dummy field to allow commas after the previous method
    dummy = 0
  )
)

