#' @include HelperFunctions.R InsuranceParameters.R ProfitParticipation_Functions.R
#'
#' @import R6
#' @import dplyr
NULL
# Shut up the package checker:
# if (getRversion() >= "2.15.1")  utils::globalVariables(c("profitClass"))

#' Filter the whole data.frame of profit rates for the given profit classes
#'
#' This is a rather trivial helper function, which just calls [dplyr::filter()].
#'
#' @param rates data.frame containing all profit rates for multiple profit classes
#' @param classes the profit classes, for which rates should be extracted
#'
#' @export
filterProfitRates = function(rates, classes) {
    rates %>%
        dplyr::filter(.data$profitClass %in% classes)
}


#' Base Class for Profit Participation Schemes
#'
#' Base class for Profit Participation schemes  (holding contract-independent values and
#' providing methods to calculate the profit participation values from the given
#' reserves).
#'
#' The profit participation object is typically not used directly, but rather
#' defined once and then passed on to an [InsuranceTarif] or [InsuranceContract]
#' object, where it will be used internally when profit participation is
#' calculated.
#'
#' This class provides the technical implementation of a profit plan for traditional
#' life insurance contracts with a guaranteed component (calculated before the
#' profit scheme comes into play) and a discretionary profit on top.
#'
#' @param params Contract-specific, full set of parameters of the contract
#'      (merged parameters of the defaults, the tariff, the profit participation
#'      scheme and the contract)
#' @param values Contract values calculated so far (guaranteed component of the
#'     insurance contract, including cash flows, premiums, reserves etc.).
#'
#' @export
ProfitParticipation = R6Class(
    "ProfitParticipation",
    public  = list(
        #' @field name The human-readable name of the profit plan.
        name  = "Name des Gewinnplans",
        #' @field Parameters Parameter template for profit-participation-specific
        #' parameters, i.e. the \code{ProfitParticipation} element of the
        #' [InsuranceContract.ParameterStructure] data structure.
        #'
        #' All elements defined in the profit scheme can be overriden per contract
        #' in the call to \code{[InsuranceContract]$new} or even in the explicit
        #' call to \ifelse{html}{\href{../../LifeInsureR/html/InsuranceContract.html#method-profitScenario}{\code{InsuranceContract$profitScenario()}}}{\code{InsuranceContract$profitScenario()()}}
        #' or \ifelse{html}{\href{../../LifeInsureR/html/InsuranceContract.html#method-addProfitScenario}{\code{InsuranceContract$addProfitScenario()}}}{\code{InsuranceContract$addProfitScenario()()}}.
        #'
        #'
        Parameters = InsuranceContract.ParameterStructure$ProfitParticipation,

        ########################################################################m#
        # Function blocks (modular) to determine bases, rates and calculation ##m#
        ########################################################################m#
        #' @field Functions list of functions defined to calculate the individual
        #' components. For each of the profit components
        #' \itemize{
        #'     \item interest profit
        #'     \item risk profit
        #'     \item expense profit
        #'     \item sum profit
        #'     \item terminal bonus
        #'     \item terminal bonus fund
        #' }
        #' a rate,  a profit base and a calculation function can be defined, by assigning one of the pre-defined
        #' [ProfitParticipationFunctions] or proving your own function with signature
        #' \code{function(rates, params, values, ...)}. Additionally, for each of the
        #' benefit types (survival, death, surrender, premium waiver) a function can
        #' be provided to calculate the benefit stemming from profit participation.
        Functions = list(

            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # Calculation bases for the various types of profit
            # Can / shall be overridden in child classes that use other bases!
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            getInterestProfitBase   = PP.base.meanContractualReserve,
            getRiskProfitBase       = PP.base.ZillmerRiskPremium,
            getExpenseProfitBase    = PP.base.sumInsured,
            getSumProfitBase        = PP.base.sumInsured,
            getTerminalBonusBase    = PP.base.sumInsured,
            getTerminalBonusFundBase = PP.base.totalProfitAssignment,

            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # Profit rates for the various types of profit
            # Can / shall be overridden in child classes that use other schemes!
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            getInterestProfitRate   = PP.rate.interestProfit,
            getRiskProfitRate       = PP.rate.riskProfit,
            getExpenseProfitRate    = PP.rate.expenseProfit,
            getSumProfitRate        = PP.rate.sumProfit,
            getTerminalBonusRate    = PP.rate.terminalBonus,
            getTerminalBonusFundRate= PP.rate.terminalBonusFund,


            getInterestOnProfits    = PP.rate.totalInterest,



            calculateInterestProfit = PP.calculate.RateOnBase,
            calculateRiskProfit     = PP.calculate.RateOnBase,
            calculateExpenseProfit  = PP.calculate.RateOnBase,
            calculateSumProfit      = PP.calculate.RateOnBase,

            calculateTerminalBonus  = PP.calculate.RateOnBase,
            getTerminalBonusReserve = function(profits, rates, terminalBonus, terminalBonusAccount, params, values) {
                n = length(terminalBonusAccount)
                terminalBonusAccount * 1/(1.07) ^ ((n - 1):0)
            },
            calculateTerminalBonusFund = PP.calculate.RateOnBase,

            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # Calculations of the assigned profit amounts, based on the bases and
            # rates defined with the functions above.
            # Can / shall be overridden in child classes that use other bases!
            #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            calculateSurvivalBenefit      = PP.benefit.ProfitPlusTerminalBonusReserve,
            calculateDeathBenefitAccrued  = PP.benefit.ProfitPlusGuaranteedInterest,
            calculateDeathBenefitTerminal = PP.benefit.TerminalBonus5YearsProRata,
            calculateSurrenderBenefitAccrued = PP.benefit.ProfitPlusHalfGuaranteedInterest,
            calculateSurrenderBenefitTerminal = PP.benefit.TerminalBonus5YearsProRata,
            calculatePremiumWaiverBenefitAccrued = PP.benefit.Profit,
            calculatePremiumWaiverBenefitTerminal = PP.benefit.TerminalBonus5YearsProRata,


            dummy = 0
        ),


        #' @description Create a new profit participation scheme
        #' @details This function is called when a new profit participation scheme
        #' is created with a call to \code{ProfitParticipation$new(...)}. Possible
        #' parameters to the \code{new}-Call are all parameters from the `ProfitParticipation` sublist of the
        #' [InsuranceContract.ParameterStructure] parameter
        #' structure (which are understood as template values that can be overridden
        #' per contract or even per profit participation scenario) and the components
        #' of the `Functions` field defining the functions to calculate the individual
        #' components of the profit participation (rates, calculation bases, calculation, benefits)
        #'
        #' @param name The name of the profit scheme (typicall the name of the profit plan and its version)
        #' @param ... profit participation parameters to be stored in the
        #' `Parameters field or calculation functions to be stored in the `Functions`
        #' field
        initialize = function(name = NULL, ...) {
            if (!missing(name))           self$name = name;
            self$setParameters(...);
            self$setFunctions(...);
            self$setFallbackParameters();
        },


        #' @description Store all passed parameters in the `Parameters` field
        #' @param ... any of the named fields defined in the `ProfitParticipation` sublist of the
        #'      [InsuranceContract.ParameterStructure]. All other arguments will be ignored
        setParameters = function(...) {
            self$Parameters = fillFields(self$Parameters, list(...));
        },

        #' @description Store all passed functions in the `Functions` field
        #' @param ... any of the functions defined in the `Functions` field. All other
        #'     arguments will be ignored
        setFunctions = function(...) {
            self$Functions = fillFields(self$Functions, list(...));
        },

        #' @description Fill all missing parameters with the default fall-back values
        setFallbackParameters = function() {
            self$Parameters = fallbackFields(self$Parameters, list(profitParticipationScheme = self));
            self$Parameters = fallbackFields(self$Parameters, InsuranceContract.ParameterDefaults$ProfitParticipation);
        },

        #' @description create a copy of a profit scheme with certain parameters changed
        #' @details This method \code{createModification} returns a copy of the profit scheme
        #' with all given arguments changed in the schmes's `Parameters`
        #' parameter list.
        #'
        #' As ProfitParticipation is a R6 class with reference logic, simply assigning
        #' the object to a new variable does not create a copy, but references the
        #' original profit scheme object. To create an actual copy, one needs to call this
        #' method, which first clones the whole object and then adjusts all parameters
        #' to the values passed to this method.
        #'
        #' @param name The new name for the cloned [ProfitParticipation] object
        #' @param ... Parameters for the [InsuranceContract.ParameterStructure],
        #'            defining the characteristics of the tariff.
        createModification = function(name  = NULL, ...) {
            cloned = self$clone();
            if (!missing(name)) cloned$name = name;
            cloned$Parameters = fillFields(cloned$Parameters, list(...));
            cloned$Functions  = fillFields(cloned$Functions, list(...));
            cloned
        },


        ##########################################################################m#
        # Advance Profit Participation                                          ####
        ##########################################################################m#

        #' @description Calculate and return the advance profit participation (to be
        #' applied on the actuarial gross premium)
        #'
        #' @details The [InsuranceContract]'s param structure [InsuranceContract.ParameterStructure]
        #'     contains the field \code{params$ProfitParticipation$advanceProfitParticipation},
        #'     which can either be numeric rate for advance profit participation, or
        #'     a function with signature \code{function(params, values, ...)} that
        #'     returns the advance profit participation rate when called with the
        #'     contract's parameters and the values calculated so far (cash flows and premiums)
        #' @return Return either one numerical value (constant for the whole premium payment period)
        #' of a vector of numerical values for the whole contract period
        #'
        #' @param ... optional parameters, to be passed to the advanceProfitParticipation
        #'     field of the parameter structure (if that is a function)
        getAdvanceProfitParticipation = function(params, values, ...) {
            valueOrFunction(params$ProfitParticipation$advanceProfitParticipation, params, values, ...)
        },


        #' @description Calculate and return the advance profit participation (to be
        #' applied after unit costs are added to the gross premium)
        #'
        #' @details The [InsuranceContract]'s param structure [InsuranceContract.ParameterStructure]
        #'     contains the field \code{params$ProfitParticipation$advanceProfitParticipationInclUnitCost},
        #'     which can either be numeric rate for advance profit participation, or
        #'     a function with signature \code{function(params, values, ...)} that
        #'     returns the advance profit participation rate when called with the
        #'     contract's parameters and the values calculated so far (cash flows and premiums)
        #' @return Return either one numerical value (constant for the whole premium payment period)
        #' of a vector of numerical values for the whole contract period
        #'
        #' @param ... optional parameters, to be passed to the advanceProfitParticipationInclUnitCost
        #'     field of the parameter structure (if that is a function)
        getAdvanceProfitParticipationAfterUnitCosts = function(params, values, ...) {
            valueOrFunction(params$ProfitParticipation$advanceProfitParticipationInclUnitCost, params, values, ...)
        },


        ##########################################################################m#
        # Traditional Profit participation:                                     ####
        #   - Interest profit                                                   ##m#
        #   - Risk profit                                                       ##m#
        #   - Expense profit                                                    ##m#
        #   - Sum profit                                                        ##m#
        #   - Terminal bonus                                                    ##m#
        ##########################################################################m#


        #' @description Set up the data.frame containing the profit participation rates
        #' @param ... additional parameters passed to the profit calculation functions
        #'     stored in the `Functions` field.
        setupRates = function(params, values, ...) {
            # 1) Profit scheme or contract provides general company-wide profit rates for some years:
            #       profitRates
            # 2) Contract can override individual rates (either for calendar years or contract years):
            #       guaranteedInterest, interestProfitRate, totalInterest, mortalityProfitRate,
            #       expenseProfitRate, sumProfitRate, terminalBonusRate, terminalBonusFundRate
            # 3) Explicit function arguments (either for calendar years or contract years).
            # 4) Any missing values will be taken from the last given year
            startYear = year(params$ContractData$contractClosing);
            policyPeriod = params$ContractData$policyPeriod;
            years = startYear:(startYear + policyPeriod);

            columns = c(
                "year",
                "guaranteedInterest", "interestProfitRate", "totalInterest", "interestProfitRate2", "totalInterest2",
                "mortalityProfitRate", "expenseProfitRate", "expenseProfitRate_premiumfree",
                "sumProfitRate",
                "terminalBonusRate", "terminalBonusFundRate")
            rates = data.frame(matrix(ncol = length(columns), nrow = length(years), dimnames = list(years = years, rates = columns)))
            rates$year = years;

            # 1) profitRates are general company-wide tables with all default rates
            # => use them as default, unless overridden
            defrates = params$ProfitParticipation$profitRates
            if (!is.null(defrates)) {
                profclass = params$ProfitParticipation$profitClass
                defrates = dplyr::filter(defrates, profitClass == profclass)
                rownames(defrates) = defrates$year

                defcols = intersect(columns, colnames(defrates))

                rates[rownames(defrates), defcols] = defrates[, defcols]
            }

            # 2) Add the explicit overrides per profit rate (from the contract)
            for (col in columns) {
                if (!is.null(params$ProfitParticipation[[col]])) {
                    rt = valueOrFunction(params$ProfitParticipation[[col]], params = params, values = values);
                    if (is.null(names(rt))) {
                        # numeric value or vector => assume values from the contract start year
                        # and fill ALL policy years with the given rates (last entry repeated)
                        rates[as.character(years), col] = padLast(rt, length(years));
                    } else {
                        # values with years assigned => use them only for the given year(s), don't override any other year
                        rates[names(rt), col] = rt;
                        rates[names(rt), "year"] = names(rt);
                    }
                }
            }

            # 3) Use explicit function param overrides per profit rate (from this function call)
            argcols = match.call(expand.dots = FALSE)$`...`;
            matching.cols = intersect(columns, names(argcols))
            for (col in matching.cols) {
                rt = eval(argcols[[col]]);
                if (!is.null(rt)) {
                    rt = valueOrFunction(rt, params = params, values = values);
                    if (is.null(names(rt))) {
                        # numeric value or vector => assume values from the contract start year
                        # and fill ALL policy years with the given rates (last entry repeated)
                        rates[as.character(years), col] = padLast(rt, length(years));
                    } else {
                        # values with years assigned => use them only for the given year(s), don't override any other year
                        rates[names(rt), col] = rt;
                        rates[names(rt), "year"] = names(rt);
                    }
                }
            }
            rownames(rates) = rates[, "year"]

            # 4) Fill all NAs with the last known value
            # First, make sure that all entries are in the correct order (sorted by year)
            rates
            newrates = rates %>% dplyr::arrange(year) %>% dplyr::mutate_all(fillNAgaps)
            rownames(newrates) = newrates$year

            # 5) Replace all NA values with 0, so we don't run into any problems later on
            newrates[is.na(newrates)] <- 0

            # TODO: Fix guaranteedInterest + interestProfitRate = totalInterest, where one of them might be missing!

            # Return only the policy years...
            self$adjustRates(newrates[as.character(years),], params = params, values = values)
        },

        #' @description Adjust the data.frame of profit participation rates after their setup
        #' @details This function provides an easy way to modify the whole set of
        #' profit rates after their initial setup. Possible applications are waiting
        #' periods, which can be implemented once for all rates rather than inside
        #' each individual calculation period.
        #' @param rates data.frame of profit paticipation rates
        adjustRates = function(rates, params, values) {
            rates[1,] = 0;
            rates
        },





        #' @description Calculation the full time series of profit participation for
        #' the given contract values
        #'
        #' @param calculateFrom The time from which to start calculating the profit
        #' participation. When a contract is changed at some time t (possibly even
        #' changing the profit scheme), all future profit participation needs to be
        #' re-calculated from that time on, without changing past profit participation.
        #' All values before \code{calculateFrom} will not be calculated.
        #' @param profitScenario profit participation values from a previous calculation
        #' (NULL if profit calculation is to be calculated from the contract inception).
        #' Values before \code{calculateFrom} will be used from this data.frame.
        #' @param ... additional parameters to be passed to \href{../../LifeInsureR/html/ProfitParticipation.html#method-setupRates}{\code{ProfitParticipation$setupRates()}}
        getProfitParticipation = function(calculateFrom = 0, profitScenario = NULL, params, values, ...) {
            waiting      = valueOrFunction(params$ProfitParticipation$waitingPeriod, params = params, values = values);
            if (is.numeric(waiting) && waiting > 0) {
                waitingFactor = c(rep(0, waiting + 1), rep(1, params$ContractData$policyPeriod - waiting));
        } else {
            waitingFactor = 1;
        }
        rates        = self$setupRates(params = params, values = values, ...)

        # Initialize all rates, bases and calc functions to NULL and then set
        # only those that are actually used in this profit scheme (all values
        # with NULL will silently be ignored in the cbind call)
        intBase      = riskBase     = expenseBase   = sumBase      = NULL;
        intRate      = riskRate     = expenseRate   = sumRate      = NULL;
        intProfit    = riskProfit   = expenseProfit = sumProfit   = NULL

        interestOnProfitRate = self$Functions$getInterestOnProfits(rates = rates, params = params, values = values);
        if ("interest" %in% params$ProfitParticipation$profitComponents) {
          intBase      = self$Functions$getInterestProfitBase(rates = rates, params = params, values = values);
          intRate      = self$Functions$getInterestProfitRate(rates = rates, params = params, values = values);
          intProfit     = self$Functions$calculateInterestProfit(base = intBase, rate = intRate, waiting = waitingFactor, rates = rates, params = params, values = values);
        }
        if ("risk" %in% params$ProfitParticipation$profitComponents) {
          riskBase     = self$Functions$getRiskProfitBase(rates = rates, params = params, values = values);
          riskRate     = self$Functions$getRiskProfitRate(rates = rates, params = params, values = values);
          riskProfit    = self$Functions$calculateRiskProfit(base = riskBase, rate = riskRate, waiting = waitingFactor, rates = rates, params = params, values = values);
        }
        if ("expense" %in% params$ProfitParticipation$profitComponents) {
          expenseBase  = self$Functions$getExpenseProfitBase(rates = rates, params = params, values = values);
          expenseRate  = self$Functions$getExpenseProfitRate(rates = rates, params = params, values = values);
          expenseProfit = self$Functions$calculateExpenseProfit(base = expenseBase, rate = expenseRate, waiting = waitingFactor, rates = rates, params = params, values = values);
        }
        if ("sum" %in% params$ProfitParticipation$profitComponents) {
          sumBase      = self$Functions$getSumProfitBase(rates = rates, params = params, values = values);
          sumRate      = self$Functions$getSumProfitRate(rates = rates, params = params, values = values);
          sumProfit     = self$Functions$calculateSumProfit(base = sumBase, rate = sumRate, waiting = waitingFactor, rates = rates, params = params, values = values);
        }

        res = cbind(
            # Profit Calculation Bases
            interestBase = c(intBase),
            riskBase = c(riskBase),
            expenseBase = c(expenseBase),
            sumBase = c(sumBase),

            # Profit Rates
            guaranteedInterest = c(rates$guaranteedInterest),
            totalInterest = c(rates$totalInterest),
            interestProfitRate = c(intRate),
            riskProfitRate = c(riskRate),
            expenseProfitRate = c(expenseRate),
            sumProfitRate = c(sumRate),
            interestOnProfitRate = c(interestOnProfitRate),

            # Profit components
            interestProfit = c(intProfit),
            riskProfit = c(riskProfit),
            expenseProfit = c(expenseProfit),
            sumProfit = c(sumProfit),
            componentsProfit = plusNULL(intProfit, riskProfit, expenseProfit, sumProfit),
            interestOnProfit = c(0),
            totalProfitAssignment = c(0),

            totalProfit = c(0)
        );

        # Use only newly calculated values starting at 'calculateFrom', but use
        # old values up to that moment (might be a contract change with a
        # completely different profit participation system before!)
        res = mergeValues(starting = profitScenario[,colnames(res)], ending = res, t = calculateFrom);

        if (calculateFrom > 0 && !is.null(profitScenario)) {
          prev = profitScenario[calculateFrom - 1, "totalProfit"]
        } else {
          prev = 0;
        }
        # TODO: turn the interest on profit into a calculator function!
        # res = self$Functions$calculateInterestOnProfit(base = sumBase, rate = sumRate, waiting = waitingFactor, rates = rates, params = params, values = values);
        for (i in (calculateFrom + 1):nrow(res)) {
            res[i,"interestOnProfit"] = res[i,"interestOnProfitRate"] * prev;
            res[i,"totalProfitAssignment"] = res[i, "componentsProfit"] + res[i,"interestOnProfit"];
            res[i,"totalProfit"] = prev + res[i,"totalProfitAssignment"];
            prev = res[i,"totalProfit"];
        }
        regularBonusAssignment = res[,"totalProfitAssignment"]

        ###########################################################################################################%#
        #### OLD Terminal bonus (not through terminal bonus fund, i.e. part of ongoing profits, but in addition) ####
        #### Terminal Bonus calculations (might depend on the individual profit assignments calculated above!
        ###########################################################################################################%#
        #### => TODO: Pass the current profit calculation inside the values!)
        if ("terminal" %in% params$ProfitParticipation$profitComponents) {
          terminalBase = self$Functions$getTerminalBonusBase(res, rates = rates, params = params, values = values);
          terminalRate = self$Functions$getTerminalBonusRate(res, rates = rates, params = params, values = values);
          terminalBonus = self$Functions$calculateTerminalBonus(res,
                  base = terminalBase, rate = terminalRate, calculateFrom = calculateFrom,
                  waiting = waitingFactor, rates = rates, params = params, values = values); # TODO: Add the AF(v) factor!

          if (calculateFrom == 0) {
            terminalBonusAccount = cumsum(terminalBonus); # TODO: Generalize! Not every scheme uses a cumulative account!
          } else {
            past = profitScenario[1:calculateFrom, "terminalBonusAccount"]
            # Preserve values up to calculateFrom, start from the last known value at calculateFrom-1 and sum all further contributions:
            terminalBonusAccount = c(head(past, -1), cumsum(c(tail(past,1), tail(terminalBonus, -calculateFrom))))
          }
          terminalBonusReserve = self$Functions$getTerminalBonusReserve(res, rates = rates, terminalBonus, terminalBonusAccount, params = params, values = values)

          resTerminal = cbind(
            terminalBase,
            terminalRate,
            terminalBonus,
            terminalBonusAccount,
            terminalBonusReserve
          )
          resTerminal = mergeValues(starting = profitScenario[,colnames(resTerminal)], ending = resTerminal, t = calculateFrom)
          # Add the terminal bonus values to the array:
          res = cbind(res, resTerminal)
        }

        ###########################################################################################################%#
        #### NEW Terminal bonus fund (part of regular profits, but not paid out on surrender, reserved as part of the free RfB) ####
        ###########################################################################################################%#
        if ("TBF" %in% params$ProfitParticipation$profitComponents) {
          TBFBase = self$Functions$getTerminalBonusFundBase(res, rates = rates, params = params, values = values);
          TBFRate = self$Functions$getTerminalBonusFundRate(res, rates = rates, params = params, values = values);
          TBFBonusAssignment = self$Functions$calculateTerminalBonusFund(res,
                base = TBFBase, rate = TBFRate, calculateFrom = calculateFrom,
                waiting = waitingFactor, rates = rates, params = params, values = values);
          regularBonusAssignment = res[,"totalProfitAssignment"] - TBFBonusAssignment

          # Calculate TBF and regular bonus as cumulative sum  of the assignments starting at t = calculateFrom plus the previous value!
          if (calculateFrom == 0) {
            TBF = cumsum(TBFBonusAssignment)
          } else {
            past = profitScenario[1:calculateFrom, "TBF"]
            # Preserve values up to calculateFrom, start from the last known value at calculateFrom-1 and sum all further contributions:
            TBF = c(head(past, -1), cumsum(c(tail(past,1), tail(TBFBonusAssignment, -calculateFrom))))
          }


          resTBF = cbind(
            TBFBase,
            TBFRate,
            TBFBonusAssignment,
            TBF
          )
          resTBF = mergeValues(starting = profitScenario[,colnames(resTBF)], ending = resTBF, t = calculateFrom)
          # Add the terminal bonus fund values to the array:
          res = cbind(res, resTBF)
        }

        ###########################################################################################################%#
        #### Regular bonus assignment / accrued regular bonus AFTER TBF                                          ####
        ###########################################################################################################%#
        # Calculate regular bonus (after potential TBF subtraction) as cumulative sum  of the assignments starting at t = calculateFrom plus the previous value!
        if (calculateFrom == 0) {
          regularBonus = cumsum(regularBonusAssignment)
        } else {
          past = profitScenario[1:calculateFrom, "regularBonus"]
          regularBonus = c(head(past, -1), cumsum(c(tail(past,1), tail(regularBonusAssignment, -calculateFrom))))
        }
        resRegular = cbind(regularBonusAssignment, regularBonus)
        resRegular = mergeValues(starting = profitScenario[,colnames(resRegular)], ending = resRegular, t = calculateFrom)
        res = cbind(res, resRegular)


        ###########################################################################################################%#
        #### BENEFITS                                                                                            ####
        ###########################################################################################################%#

        survival       = self$Functions$calculateSurvivalBenefit(res, rates = rates, params = params, values = values);

        deathAccrued   = self$Functions$calculateDeathBenefitAccrued(res, rates = rates, params = params, values = values);
        deathTerminalBonus = self$Functions$calculateDeathBenefitTerminal(res, rates = rates, params = params, values = values);

        surrenderAccrued  = self$Functions$calculateSurrenderBenefitAccrued(res, rates = rates, params = params, values = values);
        surrenderTerminalBonus = self$Functions$calculateSurrenderBenefitTerminal(res, rates = rates, params = params, values = values);

        premiumWaiverAccrued  = self$Functions$calculatePremiumWaiverBenefitAccrued(res, rates = rates, params = params, values = values);
        premiumWaiverTerminalBonus = self$Functions$calculatePremiumWaiverBenefitTerminal(res, rates = rates, params = params, values = values);

        resBenefit = cbind(
          survival = survival,

          deathAccrued = deathAccrued,
          deathTerminalBonus = deathTerminalBonus,
          death = deathAccrued + deathTerminalBonus,

          surrenderAccrued  = surrenderAccrued,
          surrenderTerminalBonus = surrenderTerminalBonus,
          surrender      = surrenderAccrued + surrenderTerminalBonus,

          premiumWaiverAccrued  = premiumWaiverAccrued,
          premiumWaiverTerminalBonus = premiumWaiverTerminalBonus,
          premiumWaiver  = premiumWaiverAccrued + premiumWaiverTerminalBonus
        )
        # Preserve values up to time t=calculateFrom of the old scenario values
        resBenefit = mergeValues(starting = profitScenario[,colnames(resBenefit)], ending = resBenefit, t = calculateFrom)

        res = cbind(res, resBenefit);
        res
    },



   #' @field dummy Dummy to allow commas in the previous method
    dummy = 0
  )
)
