#' @importFrom abind abind
#' @importFrom objectProperties setSingleEnum
#' @importFrom utils head tail
#' @importFrom methods new
#' @importFrom lubridate time_length interval
NULL



#' Enum to describe when a benefit or premium payment is due (in advance or in arrears)
#'
#' @details Currently, only two values are allowed;
#' \itemize{
#'     \item "in advance"
#'     \item "in arrears"
#' }
#'
#' @export
PaymentTimeEnum = objectProperties::setSingleEnum("PaymentTime", levels = c("in advance", "in arrears"));

#' Enum to describe possible sexes in an insurance contract or tariff.
#'
#' @details
#' Currently, the only possible values are:
#' * "unisex"
#' * "male"
#' * "female"
#'
#' @export
SexEnum = objectProperties::setSingleEnum("Sex", levels = c("unisex", "male", "female"));

#' Enum to define how much of a contract needs to be calculated automatically.
#'
#' @details
#' When an [InsuranceContract] object is created, all time series are immediately
#' calculated. However, sometimes, one only needs part of the values, so it
#' would be a waste of resources to calculate e.g. all future reserves and
#' profit participation, if only premiums are of interest.
#'
#' Possible values are:
#' * "all"
#' * "probabilities"
#' * "cashflows"
#' * "presentvalues"
#' * "premiums"
#' * "absvalues"
#' * "reserves"
#' * "premiumcomposition"
#' * "profitparticipation"
#' * "history"
#'
#' @export
CalculationEnum = objectProperties::setSingleEnum("Calculation",
    levels = c(
      "all",
      "probabilities",
      "cashflows",
      "presentvalues",
      "premiums",
      "absvalues",
      "reserves",
      "premiumcomposition",
      "profitparticipation",
      "history"
    )
)


#' Enum to define the different components of profit participation.
#'
#' @details
#' Profit participation schemes typically consist of different components,
#' which are calculated independently. Typical components are interest profit
#' to distribute investment gains to the customer, risk profit and expense profit
#' to return security margins in the biometric risk and the expenses to the customer
#' and sum profit, which aplies to contracts with higher sums insured, where
#' charged expenses are calculated from the sum insured, while the actual
#' expenses are more or less constant. Thus, high contracts are charged more,
#' which causes profits that are returned as sum profit.
#'
#' As a special case, part of the profits can be stored in a terminal bonus
#' reserve and only distributed on maturity (or potentially on death), but
#' not on surrender. Some (older) profit participation schemes add an independently
#' calculated bonus on maturity (e.g. twice the total profit assignment of the
#' last year) at maturity to give customers an additional incentive not to
#' surrender a contract.
#'
#' Possible values are (multiple can be given):
#' * "interest"
#' * "risk"
#' * "expense"
#' * "sum"
#' * "terminal"
#' * "TBF"
#'
#' @export
ProfitComponentsEnum = objectProperties::setMultipleEnum("ProfitComponents",
    levels = c(
      "advance",
      "interest",
      "risk",
      "expense",
      "sum",
      "terminal",
      "TBF"
    )
)

#' Determine whether a contract (given all parameters) is a single-premium contract or with regular premiums
#'
#' Single premium contracts are identified by the parameter \code{premiumPeriod = 1}.
#'
#' @param params The parameters of the contract.
#' @param values Unused by default (already calculated values of the contract)
#' @returns TRUE or FALSE, indicating a single-premium contract or not
#'
#' @export
isSinglePremiumContract = function(params, values) { params$ContractData$premiumPeriod <= 1 }


#' Determine whether a contract (given all parameters) is a contract with regular premiums
#'
#' Regular premium contracts are identified by the parameter \code{premiumPeriod > 1}.
#'
#' @param params The parameters of the contract.
#' @param values Unused by default (already calculated values of the contract)
#' @returns TRUE or FALSE, indicating a contract with regular premiums
#'
#' @export
isRegularPremiumContract = function(params, values) { params$ContractData$premiumPeriod > 1 }


#' Describes the death benefit of a linearly decreasing whole life insurance (after a possible deferall period)
#'
#' The death benefit will be the full sumInsured for the first year after the
#' deferral period and then decrease linearly to 0 at the end of the policyPeriod.
#' This can be used with the \code{deathBenefit} parameter for insurance
#' contracts, but should not be called directly.
#'
#' @param len The desired length of the Cash flow vector (can be shorter than
#'            the policyPeriod, if q_x=1 before the end of the contract, e.g.
#'            for life-long insurances)
#' @param params The full parameter set of the insurance contract (including
#'               all inherited values from the tariff and the profit participation)
#' @param values The values calculated from the insurance contract so far
#' @returns vector of linearly decreasing benefits over the whole insurance coverage period
#'
#' @export
deathBenefit.linearDecreasing = function(len, params, values) {
  protectionPeriod = params$ContractData$policyPeriod - params$ContractData$deferralPeriod;
  pad0((protectionPeriod:0) / protectionPeriod, l = len)
}


#' Default premium refund period: for deferred contracts the deferral period, otherwise the whole contract
#'
#' If a premium refund is set for the tariff, the default is the full contract
#' period, except for deferred contracts (typically deferred life annuities),
#' for which the deferral period is the refund period.
#'
#' @param params The full parameter set of the insurance contract (including
#'               all inherited values from the tariff and the profit participation)
#' @param values The values calculated from the insurance contract so far
#' @returns The length of the premium refund period (deferral period for deferred contracts, contract period otherwise)
#'
#' @export
premiumRefundPeriod.default = function(params, values) {
  ifelse(params$ContractData$deferralPeriod > 0, params$ContractData$deferralPeriod, params$ContractData$policyPeriod)
}


#' Describes the death benefit of a decreasing whole life insurance (after a possible deferall period)
#'
#' The death benefit will be the full sumInsured for the first year after the
#' deferral period and then decrease like an annuity to 0 at the end of the policyPeriod.
#' This can be used with the \code{deathBenefit} parameter for insurance
#' contracts, but should not be called directly.
#'
#' This function is a mere generator function, which takes the interest rate and
#' generates a function that describes a decreasing annuity.
#'
#' The generated function has the following parameters:
#' \describe{
#'     \item{len}{The desired length of the Cash flow vector (can be shorter than
#'            the policyPeriod, if q_x=1 before the end of the contract, e.g.
#'            for life-long insurances)}
#'     \item{params}{The full parameter set of the insurance contract (including
#'               all inherited values from the tariff and the profit participation)}
#'     \item{values}{The values calculated from the insurance contract so far}
#' }
#'
#'
#' @param interest The interest rate of the loan, which is underlying the insurance.
#' @returns vector of decreasing death benefits
#'
#'
#' @export
deathBenefit.annuityDecreasing = function(interest) {
  function(len, params, values) {
    protectionPeriod = params$ContractData$policyPeriod - params$ContractData$deferralPeriod;
    vk = 1/(1 + interest);
    if (interest == 0) {
      benefit = (protectionPeriod:0) / protectionPeriod
    } else {
      benefit = (vk ^ (protectionPeriod:0) - 1) / (vk ^ protectionPeriod - 1)
    }
    pad0(benefit, l = len)
  }
}


#' Calculate the age of the insured based on exact age at contract closing, rounded
#' to the nearest birthday.
#'
#' @param params The parameters of the contract.
#' @param values Unused by default (already calculated values of the contract)
#' @returns the resulting age
#'
#' @export
age.exactRounded = function(params, values) {
  round(time_length(
    interval(params$ContractData$birthDate, params$ContractData$contractClosing),
  "years"))
}

#' Calculate the age of the insured based on the difference of the bith year and
#' contract closing year.
#'
#' @param params The parameters of the contract.
#' @param values Unused by default (already calculated values of the contract)
#' @returns the resulting age
#'
#' @export
age.yearDifference = function(params, values) {
  year(params$ContractData$contractClosing) - year(params$ContractData$birthDate)
}


#' Defines a frequency charge (surcharge for monthly/quarterly/semiannual) premium payments
#'
#' @description Tariffs are typically calculated with yearly premium installments. When
#' premiums are paid more often then one a year (in advance), the insurance
#' receives part of the premium later (or not at all in case of death), so a
#' surcharge for premium payment frequencies higher than yearly is applied to
#' the  premium, typically in the form of a percentage of the premium.
#'
#' This function generates the internal data structure to define surcharges for
#' monthly, quarterly and semiannual premium payments. The given surcharges can
#' be either given as percentage points (e.g. 1.5 means 1.5% = 0.015) or as
#' fractions of 1 (i.e. 0.015 also means 1.5% surcharge). The heuristics applied
#' to distinguish percentage points and fractions is that all values larger than 0.1
#' are understood as percentage points and values 0.1 and lower are understood
#' as fractions of 1.
#' As a consequence, a frequency charge of 10% or more MUST be given as percentage points.
#'
#' Currently, the frequency charges are internally represented as a named list,
#' \code{list("1" = 0, "2" = 0.01, "4" = 0.02, "12" = 0.03)}, but that might
#' change in the future, so it is advised to use this function rather than
#' explicitly using the named list in your code.
#'
#' @param monthly Surcharge for monthly premium payments
#' @param quarterly Surcharge for quarterly premium payments
#' @param semiannually Surcharge for semi-annual premium payments
#' @param yearly Surcharge for yearly premium payments (optiona, default is no surcharge)
#' @returns a named list of surcharges for names 1 (yearly), 2 (semiannually), 4 (quarterly) and 12 (monthly) payments
#'
#' @export
freqCharge = function(monthly = 0, quarterly = 0, semiannually = 0, yearly = 0) {
  # Apply the heuristics to allow percentage points given
  if (monthly > 0.1) monthly = monthly / 100;
  if (quarterly > 0.1) quarterly = quarterly / 100;
  if (semiannually > 0.1) semiannually = semiannually / 100;
  if (yearly > 0.1) yearly = yearly / 100;

  # internal representation for now is a named list:
  list("1" = yearly, "2" = semiannually, "4" = quarterly, "12" = monthly)
}



mergeValues = function(starting, ending, t) {
  # if either starting or ending is missing, always use the other, irrespective of t:
  if (missing(ending) || is.null(ending)) {
    starting
  } else if (missing(starting) || is.null(starting)) {
    ending
  } else if (t == 0) {
    ending
  } else {
    rbind(starting[1:t,], ending[-1:-t,])
  }
}
mergeValues3D = function(starting, ending, t) {
  # if either starting or ending is missing, always use the other, irrespective of t:
  if (missing(ending) || is.null(ending)) {
    starting
  } else if (missing(starting) || is.null(starting)) {
    ending
  } else if (t == 0) {
    ending
  } else {
    abind::abind(starting[1:t,,,], ending[-1:-t,,,], along = 1)
  }
}

#' PVfactory (R6Class for present values with arbitrary dimensions)
#'
#' provides functions to calculate present values for survival, death, dread
#' disease, invalidity and guaranteed benefits in various dimensions
#' @param qx the vector of mortality / death probabilities
#' @param m the number of yearly payments
#' @param mCorrection the list of alpha/beta to correct present values for multiple payments per year
#' @param v the discount factor (1 / (1+i))
#' @param advance the payment cashflows in advance
#' @param arrears the payment cashflows in arrears
#' @param start the time index, where the PV calculation should be based
#' @param benefits the vector of death / invalidity / disease benefits
#' @param ... other parameters (currently not used, for future use)
#' @export
#'
PVfactory = R6Class(
  "PVfactory",

  ######################### PUBLIC METHODS ################################# #
  public  = list(
    #' @description Initialize the present value factory with defaults for qx, interest and multiple payments per year
    initialize = function(qx, m = 1, mCorrection = list(alpha = 1, beta = 0), v = 1) {
      private$qx = qx;
      private$m = m;
      private$mCorrection = mCorrection;
      private$v = v;
    },
    #' @description Present values of guaranteed benefits (paid in advance or arrears, possible multiple times per year)
    guaranteed = function(advance = NULL, arrears = NULL, start = 0, ..., m = private$m, mCorrection = private$mCorrection, v = private$v) {
      # General Note: Since the CF vectors can have an arbitrary number of
      # dimensions, we cannot directly access them via advance[1,..]. Rather,
      # we have to construct the `[` function manually as a quoted expression,
      # inserting the required number of dimensions and then evaluating that
      # expression. This makes this function a little harder to read, but the
      # performance should not take a hit and it is implemented in a very
      # general way.

      cfs = list(advance, arrears)
      # https://stackoverflow.com/a/16896422/920231
      deflt = cfs[!unlist(lapply(cfs, is.null))][[1]] * 0

      if (missing(advance)     || is.null(advance))     advance = deflt
      if (missing(arrears)     || is.null(arrears))     arrears = deflt

      l = max(unlist(lapply(cfs, function(cf) if(!is.null(dim(cf))) dim(cf)[[1]] else length(cf))))

      # TODO: Make sure all CF tensors have the same number of dimensions
      dms = if (is.null(dim(advance))) length(advance + 1) else dim(advance);

      # Resulting PV tensor has one timestep more than the CF tensors!
      dms[1] = dms[1] + 1
      dmNr = if (is.null(dim(advance))) 1 else length(dim(advance))

      # To be able to access the CF tensors in arbitrary dimensions, we
      # construct the [..] operator manually by quoting it and then inserting
      # arguments matching the number of dimensions
      Qadv     = Quote(advance[]    )[c(1,2,rep(3, dmNr))];
      Qarr     = Quote(arrears[]    )[c(1,2,rep(3, dmNr))];
      Qres     = Quote(res[])[c(1,2,rep(3, dmNr))]; # Access the correct number of dimensions
      QresAss  = Quote(res <- tmp)

      VL = function(quoted, time) {
        eval(quoted %>% `[<-`(3, time))
      }

      init = VL(Qadv, 1) * 0;

      # assuming advance and arrears have the same dimensions...
      # TODO: Pad to length l
      # advance = pad0(advance, l, value = init);
      # arrears = pad0(arrears, l, value = init);

      # TODO: Replace loop by better way (using Reduce?)
      res = array(0, dim = dms)

      # Starting value for the recursion:
      tmp = init
      QresAss[[2]] = Qres %>% `[<-`(3, dms[[1]])
      eval(QresAss)

      advcoeff = mCorrection$alpha - mCorrection$beta * (1 - v);
      arrcoeff = mCorrection$alpha - (mCorrection$beta + 1/m) * ( 1 - v);
      for (i in l:(start + 1)) {
        # coefficients for the payments (including corrections for payments during the year (using the alpha(m) and beta(m)):
        # The actual recursion:
        tmp = VL(Qadv, i) * advcoeff + VL(Qarr, i) * arrcoeff + v * VL(Qres, i+1);
        # Assign tmp to the slice res[i, ....]
        QresAss[[2]] = Qres %>% `[<-`(3, i)
        eval(QresAss)
      }
      VL(Qres, list(1:l))
    },

    #' @description Present values of survival benefits (paid in advance or arrears, possible multiple times per year)
    survival = function(advance = NULL, arrears = NULL, start = 0, ..., m = private$m, mCorrection = private$mCorrection, v = private$v) {
      # General Note: Since the CF vectors can have an arbitrary number of
      # dimensions, we cannot directly access them via advance[1,..]. Rather,
      # we have to construct the `[` function manually as a quoted expression,
      # inserting the required number of dimensions and then evaluating that
      # expression. This makes this function a little harder to read, but the
      # performance should not take a hit and it is implemented in a very
      # general way.

      cfs = list(advance, arrears)
      # https://stackoverflow.com/a/16896422/920231
      deflt = cfs[!unlist(lapply(cfs, is.null))][[1]] * 0

      if (missing(advance)     || is.null(advance))     advance = deflt
      if (missing(arrears)     || is.null(arrears))     arrears = deflt

      l = max(unlist(lapply(cfs, function(cf) if(!is.null(dim(cf))) dim(cf)[[1]] else length(cf))))

      # TODO: Make sure all CF tensors have the same number of dimensions
      dms = if (is.null(dim(advance))) length(advance) else dim(advance);

      # Resulting PV tensor has one timestep more than the CF tensors!
      dms[1] = dms[1] + 1
      dmNr = if (is.null(dim(advance))) 1 else length(dim(advance))

      # To be able to access the CF tensors in arbitrary dimensions, we
      # construct the [..] operator manually by quoting it and then inserting
      # arguments matching the number of dimensions
      Qadv     = Quote(advance[]    )[c(1,2,rep(3, dmNr))];
      Qarr     = Quote(arrears[]    )[c(1,2,rep(3, dmNr))];
      Qres     = Quote(res[])[c(1,2,rep(3, dmNr))]; # Access the correct number of dimensions
      QresAss  = Quote(res <- tmp)

      VL = function(quoted, time) {
        eval(quoted %>% `[<-`(3, time))
      }

      init = VL(Qadv, 1) * 0;

      # assuming advance and arrears have the same dimensions...
      p = pad0(private$qx$px, l, value=0);
      # TODO: Pad to length l
      # advance = pad0(advance, l, value = init);
      # arrears = pad0(arrears, l, value = init);

      # TODO: Replace loop by better way (using Reduce?)
      res = array(0, dim = dms)

      # Starting value for the recursion:
      tmp = init
      QresAss[[2]] = Qres %>% `[<-`(3, dms[[1]])
      eval(QresAss)

      advcoeff = mCorrection$alpha - mCorrection$beta * (1 - p * v);
      arrcoeff = mCorrection$alpha - (mCorrection$beta + 1/m) * (1 - p * v);
      for (i in l:(start + 1)) {
        # coefficients for the payments (including corrections for payments during the year (using the alpha(m) and beta(m)):
        # The actual recursion:
        tmp = VL(Qadv, i) * advcoeff[i] + VL(Qarr, i) * arrcoeff[i] + v * p[i] * VL(Qres, i+1);
        # Assign tmp to the slice res[i, ....]
        QresAss[[2]] = Qres %>% `[<-`(3, i)
        eval(QresAss)
      }
      VL(Qres, list(1:l))
    },

    #' @description Present values of death benefits
    death = function(benefits, start = 0, ..., v = private$v) {
      # General Note: Since the CF vectors can have an arbitrary number of
      # dimensions, we cannot directly access them via advance[1,..]. Rather,
      # we have to construct the `[` function manually as a quoted expression,
      # inserting the required number of dimensions and then evaluating that
      # expression. This makes this function a little harder to read, but the
      # performance should not take a hit and it is implemented in a very
      # general way.

      cfs = list(benefits)
      if (missing(benefits) || is.null(benefits)) return(0);

      l = max(unlist(lapply(cfs, function(cf) if(!is.null(dim(cf))) dim(cf)[[1]] else length(cf))))

      # TODO: Make sure all CF tensors have the same number of dimensions
      dms = if (is.null(dim(benefits))) length(benefits) else dim(benefits);

      # Resulting PV tensor has one timestep more than the CF tensors!
      dms[1] = dms[1] + 1
      dmNr = if (is.null(dim(benefits))) 1 else length(dim(benefits))

      # To be able to access the CF tensors in arbitrary dimensions, we
      # construct the [..] operator manually by quoting it and then inserting
      # arguments matching the number of dimensions
      Qben     = Quote(benefits[]    )[c(1,2,rep(3, dmNr))];
      Qres     = Quote(res[])[c(1,2,rep(3, dmNr))]; # Access the correct number of dimensions
      QresAss  = Quote(res <- tmp)

      VL = function(quoted, time) {
        eval(quoted %>% `[<-`(3, time))
      }

      init = VL(Qben, 1) * 0;

      # assuming advance and arrears have the same dimensions...
      p = pad0(private$qx$px, l, value = 0);
      q = pad0(private$qx$qx, l, value = 1);
      # TODO: Pad to length l
      # benefits = pad0(benefits, l, value = init);

      # TODO: Replace loop by better way (using Reduce?)
      res = array(0, dim = dms)

      # Starting value for the recursion:
      tmp = init
      QresAss[[2]] = Qres %>% `[<-`(3, dms[[1]])
      eval(QresAss)

      for (i in l:(start + 1)) {
        # coefficients for the payments (including corrections for payments during the year (using the alpha(m) and beta(m)):
        # The actual recursion:
        tmp = v * q[i] * VL(Qben, i) + v * p[i] * VL(Qres, i+1);
        # Assign tmp to the slice res[i, ....]
        QresAss[[2]] = Qres %>% `[<-`(3, i)
        eval(QresAss)
      }
      VL(Qres, list(1:l))
    },
    #' @description Present values of disease benefits
    disease = function(benefits, start = 0, ..., v = private$v) {
      # General Note: Since the CF vectors can have an arbitrary number of
      # dimensions, we cannot directly access them via advance[1,..]. Rather,
      # we have to construct the `[` function manually as a quoted expression,
      # inserting the required number of dimensions and then evaluating that
      # expression. This makes this function a little harder to read, but the
      # performance should not take a hit and it is implemented in a very
      # general way.

      cfs = list(benefits)
      if (missing(benefits) || is.null(benefits)) return(0);

      l = max(unlist(lapply(cfs, function(cf) if(!is.null(dim(cf))) dim(cf)[[1]] else length(cf))))

      # TODO: Make sure all CF tensors have the same number of dimensions
      dms = if (is.null(dim(benefits))) length(benefits) else dim(benefits);

      # Resulting PV tensor has one timestep more than the CF tensors!
      dms[1] = dms[1] + 1
      dmNr = if (is.null(dim(benefits))) 1 else length(dim(benefits))

      # To be able to access the CF tensors in arbitrary dimensions, we
      # construct the [..] operator manually by quoting it and then inserting
      # arguments matching the number of dimensions
      Qben     = Quote(benefits[]    )[c(1,2,rep(3, dmNr))];
      Qres     = Quote(res[])[c(1,2,rep(3, dmNr))]; # Access the correct number of dimensions
      QresAss  = Quote(res <- tmp)

      VL = function(quoted, time) {
        eval(quoted %>% `[<-`(3, time))
      }

      init = VL(Qben, 1) * 0;

      # assuming advance and arrears have the same dimensions...
      p = pad0(private$qx$px, l, value = 0);
      ix = pad0(private$qx$ix, l, value = 0);
      # TODO: Pad to length l
      # benefits = pad0(benefits, l, value = init);

      # TODO: Replace loop by better way (using Reduce?)
      res = array(0, dim = dms)

      # Starting value for the recursion:
      tmp = init
      QresAss[[2]] = Qres %>% `[<-`(3, dms[[1]])
      eval(QresAss)

      for (i in l:(start + 1)) {
        # coefficients for the payments (including corrections for payments during the year (using the alpha(m) and beta(m)):
        # The actual recursion:
        tmp = v * ix[i] * VL(Qben, i) + v * p[i] * VL(Qres, i+1);
        # Assign tmp to the slice res[i, ....]
        QresAss[[2]] = Qres %>% `[<-`(3, i)
        eval(QresAss)
      }
      VL(Qres, list(1:l))
    },
    # Cash flows only after death
    # This case is more complicated, as we have two possible states of
    # payments (present value conditional on active, but payments only when
    # dead => need to write the Thiele difference equations as a pair of
    # recursive equations rather than a single recursive formula...)
    #' @description Present values of guaranteed benefits after death (paid in advance or arrears, possible multiple times per year)
    afterDeath = function(advance = NULL, arrears = NULL, start = 0, ..., m = private$m, mCorrection = private$mCorrection, v = private$v) {
      # General Note: Since the CF vectors can have an arbitrary number of
      # dimensions, we cannot directly access them via advance[1,..]. Rather,
      # we have to construct the `[` function manually as a quoted expression,
      # inserting the required number of dimensions and then evaluating that
      # expression. This makes this function a little harder to read, but the
      # performance should not take a hit and it is implemented in a very
      # general way.

      cfs = list(advance, arrears)
      # https://stackoverflow.com/a/16896422/920231
      deflt = cfs[!unlist(lapply(cfs, is.null))][[1]] * 0

      if (missing(advance)     || is.null(advance))     advance = deflt
      if (missing(arrears)     || is.null(arrears))     arrears = deflt

      l = max(unlist(lapply(cfs, function(cf) if(!is.null(dim(cf))) dim(cf)[[1]] else length(cf))))

      # TODO: Make sure all CF tensors have the same number of dimensions
      dms = if (is.null(dim(advance))) length(advance) else dim(advance);

      # Resulting PV tensor has one timestep more than the CF tensors!
      dms[1] = dms[1] + 1
      dmNr = if (is.null(dim(advance))) 1 else length(dim(advance))

      # To be able to access the CF tensors in arbitrary dimensions, we
      # construct the [..] operator manually by quoting it and then inserting
      # arguments matching the number of dimensions
      Qadv     = Quote(advance[]    )[c(1,2,rep(3, dmNr))];
      Qarr     = Quote(arrears[]    )[c(1,2,rep(3, dmNr))];
      Qres     = Quote(res[])[c(1,2,rep(3, dmNr))]; # Access the correct number of dimensions
      QresAss  = Quote(res <- tmp)

      VL = function(quoted, time) {
        eval(quoted %>% `[<-`(3, time))
      }

      init = VL(Qadv, 1) * 0;

      # assuming advance and arrears have the same dimensions...
      p = pad0(private$qx$px, l, value=0);
      q = pad0(private$qx$qx, l, value=0);
      # TODO: Pad to length l
      # advance = pad0(advance, l, value = init);
      # arrears = pad0(arrears, l, value = init);

      # TODO: Replace loop by better way (using Reduce?)
      res = array(0, dim = dms)

      # Starting value for the recursion:
      prev = init;
      prev.dead = init;
      QresAss[[2]] = Qres %>% `[<-`(3, dms[[1]])
      eval(QresAss)

      advcoeff = mCorrection$alpha - mCorrection$beta * (1 - v);
      arrcoeff = mCorrection$alpha - (mCorrection$beta + 1/m) * (1 - v);
      for (i in l:(start + 1)) {
        # The actual recursion:
        tmp = p[i] * v * prev + q[i] * v * prev.dead;
        # Assign tmp to the slice res[i, ....]
        QresAss[[2]] = Qres %>% `[<-`(3, i)
        eval(QresAss)
        prev = tmp
        prev.dead = VL(Qadv, i) * advcoeff[i] + VL(Qarr, i) * arrcoeff[i] + v * prev.dead;
      }
      VL(Qres, list(1:l))
    }
  ),
  private = list(
    qx = data.frame(Alter = c(), qx = c(), ix = c(), px = c()),
    m = 1,
    mCorrection = list(alpha = 1, beta = 0),
    v = 1
  )
);



getSavingsPremium = function(reserves, v=1, survival_advance=c(0), survival_arrears=c(0)) {
  pad0(reserves[-1], length(reserves))*v - reserves + survival_advance + survival_arrears*v
}

correctionPaymentFrequency = function(i, m = 1, order = 0) {
  # 0th-order approximation
  alpha = 1;
  beta = 0;
  # negative orders mean that NO correction is done, e.g. because other means of
  # correction are used like an explicit premium frequency loading on the premium.
  if (order >= 0 ) beta = beta + (m - 1) / (2 * m);
  # For higher orders, simply add one term after the other!
  if (order >= 1)     beta = beta + (m ^ 2 - 1) / (6 * m ^ 2) * i; # S-Versicherung: *(1-i/2)
  # order 1.5 has a special term that should NOT be used for higher-order approximations!
  if (order == 1.5)   beta = beta + (1 - m ^ 2) / (12 * m ^ 2) * i ^ 2;

  if (order >= 2) {
    beta = beta + (1 - m ^ 2) / (24 * m ^ 2) * i ^ 2;
    alpha = alpha + (m ^ 2 - 1) / (12 * m ^ 2) * i ^ 2;
  }
  # Exact value
  if (order == Inf) {
    d = i / (1 + i);
    im = m * ((1 + i) ^ (1/m) - 1);
    dm = im / (1 + im/m);

    alpha = d*i / (dm*im);
    beta = (i - im) / (dm * im);
  }
  list(alpha = alpha, beta = beta);
}

#' Pad a vector with 0 to a desired length
#'
#' @param v the vector to pad with 0
#' @param l the desired (resulting) length of the vector
#' @param value the value to pad with (if padding is needed). Default to 0, but
#'     can be overridden to pad with any other value.
#' @param start the first \code{start} values are always set to 0 (default is 0,
#'     can be changed using the \code{value.start} argument),
#'     the vector \code{v} starts only after these leading zeroes. The number of
#'     leading zeroes counts towards the desired length
#' @param value.start the value to insert before the start index.
#'
#' @return returns the vector \code{v} padded to length \code{l} with value \code{value} (default 0).
#'
#' @examples
#' pad0(1:5, 7)   # Pad to length 7 with zeroes
#' pad0(1:5, 3)   # no padding, but cut at length 3
#'
#' # 3 leading zeroes, then the vector start (10 elements of vector, no additional padding needed):
#' pad0(1:10, 13, start = 3)
#'
#' # padding with value other than zero:
#' pad0(1:5, 7, value = "pad")
#' @export
pad0 = function(v, l, value = 0, start = 0, value.start = 0) {
  # 3 cases: desired length<=start => only 0/value.start
  #          desired length within start+v => cut v
  #          desired length longer than start+v => pad with 0/value
  if (l <= start) {
    rep(value.start, l)
  } else if (start <= l && l <= start + length(v)) {
    c(rep(value.start, start), v[0:(l - start)])
  } else {
    # Need padding
    c(rep(value.start, start), v, rep(value, l - length(v) - start))
  }
}



#' Set all entries of the given vector to 0 up until index 'start'
#'
#' @param v the vector to modify
#' @param start how many leading elements to zero out
#' @param value.start the value to insert before the start index.
#'
#' @return the vector \code{v} with the first \code{start} elements replaced by 0.
#'
#' @examples
#' head0(1:10, 3)
#' @export
head0 = function(v, start = 0, value.start = 0) {
  if (start == 0) {
    v
  } else {
    c(rep(value.start, start), tail(v, -start))
  }
}

#' Pad the vector \code{v} to length \code{l} by repeating the last entry of the
#' vector.
#'
#' This function is just a trivial wrapper around \code{pad0} and only calls [pad0()]
#' with the last element of the vector as padding value instead of the default 0.
#'
#' @param v the vector to pad by repeating the last element
#' @param ... arguments passed through to \code{pad0}
#' @returns vector padded to the correct length
#'
#' @examples
#' padLast(1:5, 7) # 5 is repeated twice
#' padLast(1:5, 3) # no padding needed
#'
#' @export
padLast = function(v, ...) {
    pad0(v, value = tail(v, n = 1), ...)
}

#' Replace all \code{NA} entries of a vector with the previous non-NA value
#'
#' Sometimes one has a vector with some gaps (\code{NA}) values, which cause
#' problems for several numeric functions. This function \code{fillNAgaps} fills
#' these missing values by inserting the last preceeding non-NA-value. Leading
#' NA values (at the start of the vector will not be modified). If the
#' argument \code{firstBack = TRUE}, leading \code{NA}-values are replaced by
#' the first non-NA value.
#' Trailing NAs are always replaced by the last previous NA-value.
#'
#' This code was taken from the R Cookbook:
#' http://www.cookbook-r.com/Manipulating_data/Filling_in_NAs_with_last_non-NA_value/
#' LICENSE (from that page): The R code is freely available for use without any restrictions.
#' In other words: you may reuse the R code for any purpose (and under any license).
#'
#' @param x The vector where NA-values should be filled by repeating the last preceeding non-NA value
#' @param firstBack if \code{TRUE}, leading NAs are replaced by the first non-NA
#'     value in the vector, otherwise leading NAs are left untouched.
#' @returns vector with NAs properly replaced
#'
#' @export
fillNAgaps <- function(x, firstBack=FALSE) {
    ## NA's in a vector or factor are replaced with last non-NA values
    ## If firstBack is TRUE, it will fill in leading NA's with the first
    ## non-NA value. If FALSE, it will not change leading NA's.

    # If it's a factor, store the level labels and convert to integer
    lvls <- NULL
    if (is.factor(x)) {
        lvls <- levels(x)
        x    <- as.integer(x)
    }

    goodIdx <- !is.na(x)

    # These are the non-NA values from x only
    # Add a leading NA or take the first good value, depending on firstBack
    if (firstBack)   goodVals <- c(x[goodIdx][1], x[goodIdx])
    else             goodVals <- c(NA,            x[goodIdx])

    # Fill the indices of the output vector with the indices pulled from
    # these offsets of goodVals. Add 1 to avoid indexing to zero.
    fillIdx <- cumsum(goodIdx) + 1

    x <- goodVals[fillIdx]

    # If it was originally a factor, convert it back
    if (!is.null(lvls)) {
        x <- factor(x, levels = seq_along(lvls), labels = lvls)
    }

    x
}


#' If \code{val} is a function, evaluate it, otherwise return \code{val}
#' @param val Function or value
#' @param ... Argument passed to \code{val} if it is a function
#' @returns the value, potentially with the given function applied
#' 
#' @examples
#' valueOrFunction(3) # returns 3
#' valueOrFunction(`+`, 1, 2) # also returns 3
#' A = `+`
#' valueOrFunction(A, 1, 2)
#' @export
valueOrFunction = function(val, ...) {
  if (is.function(val)) {
    val(...)
  } else {
    val
  }
}

#' If \code{hook} is a function, apply it to \code{val}, otherwise return \code{val} unchanged
#' @param hook (optional) function to apply to \code{val} and the other parameters
#' @param val The value to which the hook is applied (ifgiven)
#' @param ... optional parameters passed to the hook function (if it is a function)
#' @returns the value with the hook function applied (if a function) or unchanged otherwiese
#' 
#' @examples
#' applyHook(NULL, 3) # returns 3 unchanged
#' applyHook(function(x) 2*x, 3) # applies the function, returns 6
#' applyHook(`+`, 3, 1) # returns 4
#' @export
applyHook = function(hook, val, ...) {
  if (is.function(hook)) {
    hook(val, ...)
  } else if (is.null(hook)) {
    val
  } else {
    warning("Hook function", hook, "is neither a function nor NULL. Please provide a function or leave it empty!")
    val
  }
}



#' Overwrite all existing fields with default values given
#'
#' @description Overwrite all existing fields in the first argument with
#' values given in valuelist. Members of valuelist that are not yet in
#' fields are ignored. This allows a huge valuelist to be used to fill
#' fields in multiple lists with given structure.
#'
#' @param fields existing list
#' @param valuelist list of fields to replace in \code{fields}. Only keys that exist in \code{fields} are overwritten, no new fields are added to \code{fields}
#' @returns updated fields list
#'
#' @export
fillFields = function(fields, valuelist) {
  fieldsToInsert = intersect(names(fields), names(valuelist));
  fields[fieldsToInsert] = valuelist[fieldsToInsert]
  fields
}


#' Replace missing values in ields by default fallback values
#'
#' @description Replace all missing values in fields (either missing or NA) with
#' their corresponding values from fallback. Members in fallback that are missing
#' in fields are inserted
#' @param fields existing list
#' @param valuelist list of fields to replace in \code{fields}. Only keys that are missing in \code{fields} are added, no existing fields in \code{fields} are overwritten
#' @return fields with missing values filled from the fallback data
#' @export
fallbackFields = function(fields, valuelist) {
  keepFields = !sapply(fields, is.null);
  # We need to set all fields of valuelist, except those that are NOT NA in fields:
  useFields = setdiff(names(valuelist), names(fields[keepFields]))
  fields[useFields] = valuelist[useFields]
  fields
}

#' Calculate the rolling mean of length 2
#' @param x vector of values, for which the rolling mean is calculated
#' @examples
#' rollingmean(1:10)
#' @export
rollingmean = function(x) (tail(x, -1) + head(x, -1))/2


# Sum two or more vectors and correctly handle (i.e. ignore) NULL values given
plusNULL = function(v1, v2, ...) {
  if (missing(v2) && length(list(...)) == 0) {
    if (missing(v1) || is.null(v1)) {
      return(0)
    } else {
      return(v1)
    }
  }
  if (missing(v1) || is.null(v1)) {
    return(plusNULL(v2, ...));
  }
  if (missing(v2) || is.null(v2)) {
    return(plusNULL(v1, ...));
  } else {
    return(plusNULL(v1 + v2, ...))
  }
}


######################################################################=#
# Functions for handling sub-contract blocks                        ####

# Helper functions to prepend/append rows to the arrays and sum them up
padArray = function(arr = NULL, pad = 0, len = 0, value = 0) {
  padEnd = max(0, len - pad - NROW(arr)) # if len is too short, return an array containing at least the arr
  nrcols = ifelse(is.null(arr), 0, NCOL(arr))
  rbind(
    array(value, dim = c(pad, nrcols)) %>% `colnames<-`(colnames(arr)),
    arr,
    array(value, dim = c(padEnd, nrcols)) %>% `colnames<-`(colnames(arr))
  ) %>% `colnames<-`(colnames(arr))
}

sumPaddedArrays = function(arr1 = NULL, arr2 = NULL, pad1 = 0, pad2 = 0) {
  newlen = max(pad1 + NROW(arr1), pad2 + NROW(arr2))
  if (is.null(arr2)) {
    padArray(arr1, pad = pad1, len = newlen)
  } else if (is.null(arr1)) {
    padArray(arr2, pad = pad2, len = newlen)
  } else {
    # First prepend trailing zero rows according to pad1/pad2:
    arr1 = padArray(arr1, pad = pad1, len = newlen)
    arr2 = padArray(arr2, pad = pad2, len = newlen)

    # arr1 and arr2 now should have the same dimensions => sum them up
    arr1 + arr2
  }
}




