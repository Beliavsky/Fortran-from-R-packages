#' @include HelperFunctions.R InsuranceParameters.R
#'
#' @import dplyr
NULL

#' Helper functions for profit participation
#'
#' Various helper functions for the \code{ProfitParticipation} class that
#' provide the building blocks for the individual components of profit participation,
#' the rates and how the assigned profit is calculated.
#'
#' @name ProfitParticipationFunctions
#' @param rates data.frame of profit rates
#' @param params Contract-specific, full set of parameters of the contract
#'      (merged parameters of the defaults, the tariff, the profit participation
#'      scheme and the contract)
#' @param values Contract values calculated so far (guaranteed component of the
#'     insurance contract, including cash flows, premiums, reserves etc.).
#' @param ... Other values that might be used for the calculation (currently unused)
#'
#' @param base The profit calculation base, on which the rate is to be applied
#' @param rate The profit participation rate
#' @param waiting A possible waiting period
NULL

shiftBy = function(rate, n = 1) {
  nm = names(rate)
  res = c(rep(0, n), head(rate, -n))
  names(res) = nm
  res

}

##########################################################################m##
# Calculation bases for the various types of profit                      ####
##########################################################################m##

#' @describeIn ProfitParticipationFunctions
#' Basis for profit: NONE (i.e. always returns 0)
#' @export
PP.base.NULL = function(rates, params, values, ...) {
  rep(0, values$int$l)
};

#' @describeIn ProfitParticipationFunctions
#' Basis for profit: Previous Zillmer reserve (no administration cost reserve)
#' @export
PP.base.PreviousZillmerReserve = function(rates, params, values, ...) {
    shiftBy(values$reserves[,"Zillmer"], n = 1)
}

#' @describeIn ProfitParticipationFunctions
#' Basis for profit: Zillmer reserve (no administration cost reserve) at time t-2
#' @export
PP.base.ZillmerReserveT2 = function(rates, params, values, ...) {
    shiftBy(values$reserves[,"Zillmer"], n = 2)
}

#' @describeIn ProfitParticipationFunctions
#' Basis for profit: Contractual reserve (including administration costs) at time t
#' @export
PP.base.contractualReserve = function(rates, params, values, ...) {
  pmax(0, values$reserves[,"contractual"])
};

#' @describeIn ProfitParticipationFunctions
#' Basis for profit: Contractual reserve (including administration costs) at time t-1
#' @export
PP.base.previousContractualReserve = function(rates, params, values, ...) {
  shiftBy(values$reserves[,"contractual"], n = 1)
}

#' @describeIn ProfitParticipationFunctions
#' Basis for profit: Contractual reserve (including administration costs) averaged over t and t-1
#' @export
PP.base.meanContractualReserve = function(rates, params, values, ...) {
    # Rolling mean of the value for the current and previous year.
    pmax(0, rollingmean(c(0, values$reserves[,"contractual"])))
};

#' @describeIn ProfitParticipationFunctions
#' Basis for risk/mortality profit: Zillmer Risk Premium of the past year
#' @export
PP.base.ZillmerRiskPremium = function(rates, params, values, ...) {
    # The risk premium of t=0 is used to determine the risk profit at time
    # t=1, so shift the whole vector!
    shiftBy(values$premiumComposition[,"Zillmer.risk"], n = 1)
};

#' @describeIn ProfitParticipationFunctions
#' Basis for expense/sum profit: sum insured
#' @export
PP.base.sumInsured = function(rates, params, values, ...) {
  params$ContractData$sumInsured
};

#' @describeIn ProfitParticipationFunctions
#' Basis for Terminal Bonus Fund Assignment: total profit assignment of the year
#'
#' @param res the data.frame of reserves.
#' @export
PP.base.totalProfitAssignment = function(res, ...) {
  res[,"totalProfitAssignment"]
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Profit rates for the various types of profit
# Can / shall be overridden in child classes that use other schemes!
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#' @describeIn ProfitParticipationFunctions
#' Returns the array of interest profit rates (keyed by year)
#' @export
PP.rate.interestProfit = function(rates, ...) {
  rates$interestProfitRate
};

#' @describeIn ProfitParticipationFunctions
#' Returns the array of risk profit rates (keyed by year)
#' @export
PP.rate.riskProfit = function(rates, ...) {
  rates$mortalityProfitRate
};
#' @describeIn ProfitParticipationFunctions
#' Returns the array of expense profit rates (keyed by year)
#' @export
PP.rate.expenseProfit = function(rates, ...) {
  rates$expenseProfitRate
};

#' @describeIn ProfitParticipationFunctions
#' Returns the array of sum profit rates (keyed by year)
#' @export
PP.rate.sumProfit = function(rates, ...) {
  rates$sumProfitRate
};

#' @describeIn ProfitParticipationFunctions
#' Returns the array of terminal bonus rates (keyed by year)
#' @export
PP.rate.terminalBonus = function(rates, ...) {
  rates$terminalBonusRate
};

#' @describeIn ProfitParticipationFunctions
#' Returns the array of terminal bonus rates (keyed by year) as the terminal bonus fund ratio
#' @export
PP.rate.terminalBonusFund = function(rates, ...) {
  rates$terminalBonusFundRate
};


#' @describeIn ProfitParticipationFunctions
#' Rate for interest on past profits: total credited rate, but at least the guarantee
#' @export
PP.rate.interestProfitPlusGuarantee = function(rates, ...) {
  rates$interestProfitRate + rates$guaranteedInterest
};

#' @describeIn ProfitParticipationFunctions
#' Rate for interest on past profits: total creditedrate2, but at least the guarantee
#' @export
PP.rate.interestProfit2PlusGuarantee = function(rates, ...) {
  rates$interestProfitRate2 + rates$guaranteedInterest
};

#' @describeIn ProfitParticipationFunctions
#' Rate for interest on past profits: total interest rate
#' @export
PP.rate.totalInterest = function(rates, ...) {
  rates$totalInterest
};

#' @describeIn ProfitParticipationFunctions
#' Rate for interest on past profits: second total interest rate
#' @export
PP.rate.totalInterest2 = function(rates, ...) {
  rates$totalInterest2
};

#' @describeIn ProfitParticipationFunctions
#' Rate for interest on past profits: second interest profit rate (not including guaranteed interest), keyed by year
#' @export
PP.rate.interestProfit2 = function(rates, ...) {
  rates$interestProfitRate2
};


# TODO
#' @describeIn ProfitParticipationFunctions
#'
#' Calculate the terminal bonus reserve.
#' @param profits The data.frame of profits assigned
#' @param terminalBonus The terminal bonus calculated
#' @param terminalBonusAccount The terminal bonus account (like a bank account,
#'        where terminal bonuses are accrued, potentiall discounted from the maturity)
getTerminalBonusReserve = function(profits, rates, terminalBonus, terminalBonusAccount, params, values, ...) {
  n = length(terminalBonusAccount)
  terminalBonusAccount * 1/(1.07) ^ ((n - 1):0)
};


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculations of the assigned profit amounts, based on the bases and
# rates defined with the functions above.
# Can / shall be overridden in child classes that use other bases!
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#' @describeIn ProfitParticipationFunctions
#' Calculate profit by a simple rate applied on the basis (with an optional waiting vector of values 0 or 1)
#' @export
PP.calculate.RateOnBase = function(base, rate, waiting, rates, params, values, ...) {
  base * rate * waiting
};

#' @describeIn ProfitParticipationFunctions
#' Calculate profit by a simple rate applied on the basis (with an optional waiting vector of values 0 or 1), bound below by 0
#' @export
PP.calculate.RateOnBaseMin0 = function(base, rate, waiting, rates, params, values, ...) {
  pmax(0, base * rate * waiting)
};

#' @describeIn ProfitParticipationFunctions
#' Calculate profit by a rate + guaranteed interest applied on the basis (with an optional waiting vector of values 0 or 1)
#' @export
PP.calculate.RatePlusGuaranteeOnBase = function(base, rate, waiting, rates, params, values, ...) {
  base * (rate + rates$guaranteedInterest) * waiting
};

#' @describeIn ProfitParticipationFunctions
#' Calculate profit by a simple rate applied on the basis (with only (1-SGFFactor) put into profit participation, and an optional waiting vector of values 0 or 1)
#' @export
PP.calculate.RateOnBaseSGFFactor = function(base, rate, waiting, rates, params, values, ...) {
  base * rate * waiting * (1 - rates$terminalBonusFund)
};



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculations of the benefits, based on the assigned profit amounts and
# rates defined with the functions above.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#' @describeIn ProfitParticipationFunctions
#' Extract the given columns of the profit participation array of values and sum
#' them up. Columns that do not exist, because the profit scheme does not
#' provide the corresponding profit component will be silently ignored.
#' This allows generic benefit calculation functions to be written that do
#' not need to distinguish e.g. whether an old-style terminal bonus or a terminal
#' bonus fund is provided.
#'
#' This function is not meant to be called directly, but within a profit benefit
#' calculation function.
#'
#' @param profits The array of profit participation component values
#' @param cols The columns of the profit values array to be summed (columns given that do not exist in the profits array are ignired)
#' @export
sumProfits = function(profits, cols) {
  # extract the columns -- if they exist -- and sum them up:
  rowSums(
    profits[, intersect(cols, colnames(profits)), drop = FALSE]
  )
}

#' @describeIn ProfitParticipationFunctions
#' Calculate survival benefit as total profit amount plus the terminal bonus reserve
#' @export
PP.benefit.ProfitPlusTerminalBonusReserve = function(profits, ...) {
  sumProfits(profits, c("regularBonus", "TBF", "terminalBonusReserve"))
};

#' @describeIn ProfitParticipationFunctions
#' Calculate benefit as total profit accrued so far
#' @export
PP.benefit.Profit = function(profits, ...) {
  profits[,"regularBonus"]
};

#' @describeIn ProfitParticipationFunctions
#' Calculate accrued death benefit as total profit with (guaranteed) interest for one year
#' @export
PP.benefit.ProfitPlusGuaranteedInterest = function(profits, rates, ...) {
  profits[,"regularBonus"] * (1 + rates$guaranteedInterest)
};

#' @describeIn ProfitParticipationFunctions
#' Calculate accrued death benefit as total profit with total interest (interest on profit rate) for one year
#' @export
PP.benefit.ProfitPlusTotalInterest = function(profits, rates, params, values) {
  profits[,"regularBonus"] * (1 + profits[,"interestOnProfitRate"])
};

#' @describeIn ProfitParticipationFunctions
#' Calculate accrued benefit as total profit with total interest (interest on profit rate) for half a year
#' @export
PP.benefit.ProfitPlusHalfTotalInterest = function(profits, ...) {
  profits[,"regularBonus"] * (1 + profits[,"interestOnProfitRate"]/2)
};

#' @describeIn ProfitParticipationFunctions
#' Calculate death benefit as total profit with (guaranteed) interest for one year
#' @export
PP.benefit.ProfitPlusHalfGuaranteedInterest = function(profits, rates, ...) {
  profits[,"regularBonus"] * (1 + rates$guaranteedInterest/2)
};

#' @describeIn ProfitParticipationFunctions
#' Calculate accrued benefit as total profit with interest for one year (max of guarantee and total interest)
#' @export
PP.benefit.ProfitPlusInterestMinGuaranteeTotal = function(profits, rates, ...) {
  profits[,"regularBonus"] * (1 + pmax(rates$guaranteedInterest, rates$totalInterest))
};

#' @describeIn ProfitParticipationFunctions
#' Calculate accrued benefit as total profit with interest for half a year (max of guarantee and total interest)
#' @export
PP.benefit.ProfitPlusHalfInterestMinGuaranteeTotal = function(profits, rates, ...) {
  profits[,"regularBonus"] * (1 + pmax(rates$guaranteedInterest, rates$totalInterest)/2)
};

#' @describeIn ProfitParticipationFunctions
#' Calculate accrued benefit as regular profit, but used to cover initial Zillmerization
#' @export
PP.benefit.ProfitGuaranteeSupporting = function(profits, rates, params, values, ...) {
  pmax(0, values$reserves[,"contractual"] + profits[,"regularBonus"] - pmax(0, values$reserves[,"contractual"]))
};

#' @describeIn ProfitParticipationFunctions
#' Calculate benefit from terminal bonus as 1/n parts of the terminal bonus reserve during the last 5 years
#' @export
PP.benefit.TerminalBonus5YearsProRata = function(profits, params, ...) {
  n = params$ContractData$policyPeriod;
  sumProfits(profits, c("TBF", "terminalBonusReserve")) * (0:n)/n * ((0:n) >= max(10, n - 5))
};

#' @describeIn ProfitParticipationFunctions
#' Terminal bonus is only paid out during the last 5 years of the contract (but never during the first 10 years)
#' @export
PP.benefit.TerminalBonus5Years = function(profits, params, ...) {
  n = params$ContractData$policyPeriod;
  sumProfits(profits, c("TBF", "terminalBonusReserve")) * ((0:n) >= max(10, n - 5))
}

#' @describeIn ProfitParticipationFunctions
#' Calculate benefit from terminal bonus (full bonus), either old-style terminal bonus reserve or Terminal Bonus Fund (TBF)
#' @export
PP.benefit.TerminalBonus = function(profits, params, ...) {
  sumProfits(profits, c("TBF", "terminalBonusReserve"))
};

#' @describeIn ProfitParticipationFunctions
#' No benefit paid out
#' @export
PP.benefit.None = function(profits, ...) {
  0
};


"dummy"
