#' @include HelperFunctions.R InsuranceParameters.R InsuranceTarif.R ProfitParticipation.R
#'
#' @import MortalityTables
#' @import R6
#' @import lubridate
#  Prevent spurious imports warnings in CRAN checks:
#' @importFrom kableExtra kable
#' @importFrom pander pander
#' @importFrom rmarkdown render
NULL



############ Class InsuranceContract ###########################################
#' Base Class for Insurance Contracts
#'
#' R6 class that models a complete, general insurance contract.
#' The corresponding tariff and the profit participation scheme, as well as
#' all other relevant contract parameters (if not defined by the tariff or
#' explicitly overridden by the contract) can be given in the constructor.
#'
#' # Usage
#'
#' The typical usage of this class is to simply call
#' \ifelse{html}{\href{#method-new}{\code{InsuranceContract$new()}}}{\code{InsuranceContract$new()()}}.
#'
#' All parameters from the [InsuranceContract.ParameterDefaults] can be passed
#' to the constructor of the class (i.e. the \ifelse{html}{\href{#method-new}{\code{InsuranceContract$new()}}}{\code{InsuranceContract$new()()}}-call).
#' Parameters not explicitly given, will be taken from the tariff or as a fall-back
#' mechanism from the [InsuranceContract.ParameterDefaults] defaults.
#'
#' Immediately upon construction, all premiums, reserves and cash flows for the
#' whole contract period are calculated and can be accessed via the \code{Values}
#' field of the object.
#'
#'
#' # Calculation approach: Valuation
#'
#' The calculation of all contract values is controlled by the function
#' \ifelse{html}{\href{#method-calculateContract}{\code{InsuranceContract$calculateContract()}}}{\code{InsuranceContract$calculateContract()()}} (using methods of the [InsuranceTarif]
#' object) and follows the following logic:
#'
#' 1. First the **contingent (unit) cash flows** and the **transition probbilities**
#' are determined.
#' 2. The **actuarial equivalence principle** states that at time of inception, the
#' (net and gross) premium must be determined in a way that the present value
#' of the future benefits and costs minus the present value of the future premiums
#' must be equal, i.e. in expectation the future premiums ove the whole lifetime
#' of the contract will exactly cover the benefits and costs. Similarly, at all
#' later time steps, the difference between these two present values needs to be
#' reserved (i.e. has already been paid by the customer by previous premiums).
#' 2. This allows the premiums to be calculated by first calculating the **present
#' values** for all of the **benefit and costs cash flow** vectors.
#' 3. The formulas
#' to calculate the gross, Zillmer and net **premiums** involve simple linear
#' combinations of these present values, so the **coefficients of these formulas**
#' are determined next.
#' 4. With the coefficients of the premium formulas calculated, all **premiums
#' can be calculated** (first the gross premium, because due to potential gross
#' premium refunds in case of death, the formula for the net premium requires
#' the gross premium, which the formula for the gross premium involves no other
#' type of premuim).
#' 5. With premiums determined, all unit cash flows and unit present values can
#' now be expressed in monetary terms / as **absolute cash flows** (i.e. the actual Euro-amount that flows
#' rather than a percentage).
#' 6. As described above, the difference between the present values of premiums
#' and present values of benefits and costs is defined as the required amount
#' of reserves, so the **reserves (net, gross, administration cost, balance sheet)**
#' and all values derived from them (i.e. surrender value, sum insured in case of
#' premium waiver, etc.) are calculated.
#' 7. The **decomposition of the premium** into parts dedicated to specific purposes
#' (tax, rebates, net premium, gross premium, Zillmer premium, cost components,
#' risk premium, savings premium, etc.) can be done once the reserves are
#' ready (since e.g. the savings premium is defined as the difference of
#' discounted reserves at times $t$ and $t+1$).
#' 8. If the contract has **(discretionary or obligatory) profit sharing**B mechanisms
#' included, the corresponding [ProfitParticipation] object can calculate that
#' profit sharing amounts, once all guaranteed values are calculated. This can
#' also be triggered manually (with custom profit sharing rates) by calling
#' the methods \ifelse{html}{\href{#method-profitScenario}{\code{InsuranceContract$profitScenario()}}}{\code{InsuranceContract$profitScenario()()}}]
#' or \ifelse{html}{\href{#method-addProfitScenario}{\code{InsuranceContract$addProfitScenario()}}}{\code{InsuranceContract$addProfitScenario()()}}.
#'
#'
#'
#'
#' # Calculation approach: Cash Flows
#'
#' An insurance contract is basically defined by the (unit) cash flows it produces:
#' \itemize{
#'   \item **Premium payments** (in advance or in arrears) at each timestep
#'   \item **Survival payments** at each timestep
#'   \item **Guaranteed payments** at each timestep
#'   \item **Death benefits** at each timestep
#'   \item **Disease benefits** at each timestep
#' }
#' Together with the transition probabilities (mortalityTable parameter)
#' the present values can be calculated, from which the premiums follow and
#' finally the reserves and a potential profit sharing.
#'
#' For example, a _**term life insurance with regular premiums**_ would have the following
#' cash flows:
#'
#' * premium cash flows: 1, 1, 1, 1, 1, ...
#' * survival cash flows: 0, 0, 0, 0, 0, ...
#' * guaranteed cash flows: 0, 0, 0, 0, 0, ...
#' * death benefit cash flows: 1, 1, 1, 1, 1, ...
#'
#' A _**single-premium term life insurance**_ would look similar, except for the premiums:
#'
#' * premium cash flows: 1, 0, 0, 0, 0, ...
#'
#' A _**pure endowment**_ has no death benefits, but a survival benefit of 1 at the
#' maturity of the contract:
#'
#' * premium cash flows: 1, 1, 1, 1, 1, ...
#' * survival cash flows: 0, 0, ..., 0, 1
#' * guaranteed cash flows: 0, 0, 0, 0, 0, ...
#' * death benefit cash flows: 0, 0, 0, 0, 0, ...
#'
#' An _**endowment**_ has also death benefits during the contract duration:
#'
#' * premium cash flows: 1, 1, 1, 1, 1, ...
#' * survival cash flows: 0, 0, ..., 0, 1
#' * guaranteed cash flows: 0, 0, 0, 0, 0, ...
#' * death benefit cash flows: 1, 1, 1, 1, 1, ...
#'
#' A _**(deferred) annuity**B_ has premium cash flows only during the deferral peroid
#' and only survival cash flows during the annuity payment phase. Often, in case
#' of death during the deferral period, all premiums paid are refunded as a death
#' benefit.:
#'
#' * premium cash flows: 1, 1, ...,  1, 0, 0, 0, ...
#' * survival cash flows: 0, 0, ..., 0, 1, 1, 1,...
#' * guaranteed cash flows: 0, 0, 0, 0, 0, ...
#' * death benefit cash flows: 1, 2, 3, 4, 5, ..., 0, 0, ...
#'
#' A _**terme-fix insurance**_ has a guaranteed payment at maturity, even if the insured
#' has already died. The premiums, however, are only paid until death (which is
#' not reflected in the contingent cash flows, but rather in the transition
#' probabilities):
#'
#' * premium cash flows: 1, 1, 1, 1, ...,  1
#' * survival cash flows: 0, 0, 0, 0, ..., 0
#' * guaranteed cash flows: 0, 0, 0, ..., 0, 1
#' * death benefit cash flows: 0, 0, 0, 0, ..., 0
#'

#'
#' @examples
#' # TODO
#'
#' @export
InsuranceContract = R6Class(
    "InsuranceContract",

    ######################### PUBLIC METHODS ##################################
    public = list(
        #' @field tarif
        #' The [InsuranceTarif] underlying this contract. The tarif is the abstract
        #' product description (i.e. defining the type of insurance, fixing tpyes
        #' of benefits, specifying costs, guaranteed interest rate, mortality tables,
        #' potential profit sharing mechanisms, etc.), while the contract holds
        #' the individual parts like age, sum insured, contract duration, premium
        #' payment frequency, etc.
        tarif = NULL,
        #' @field parent
        #' A pointer to the parent contract. Some contracts consist of multiple
        #' parts (e.g. a main savings contract with a dread-disease rider, or
        #' a contract with multiple dynamic increases). These are internally
        #' represented by one [InsuranceContract] object per contract part, plus
        #' one contract object combining them and deriving combined premiums,
        #' reserves and profit participation. The child contracts (i.e. the
        #' objects representing the individual parts) have a pointer to their
        #' parent, while the overall contract holds a list of all its child contract parts.
        parent = NULL,

        #' @field ContractParameters
        #' Insurance contract parameters explicitly specified in the contract
        #' (i.e. parameters that are NOT taken from the tariff of the defaults).
        ContractParameters = InsuranceContract.ParameterStructure, # Only values explicitly given for this contract, not including fallback values from the tariff
        #' @field Parameters
        #' Full set of insurance contract parameters applying to this contract.
        #' The set of parameters is a combination of explicitly given (contract-specific)
        #' values, parameters determined by the tariff and default values.
        Parameters = InsuranceContract.ParameterStructure,         # The whole parameter set, including values given by the tariff

        #' @field Values
        #' List of all contract values (cash flows, present values, premiums,
        #' reserves, premium decomposition, profit participation, etc.). These
        #' values will be calculated and filled when the contract is created
        #' and updated whenever the contract is changed.
        Values = InsuranceContract.Values,

        #' @field blocks
        #' For contracts with multiple contract parts: List of all tariff blocks
        #' (independently calculated [InsuranceContract] objects, that are combined
        #' to one contract, e.g. dynamic/sum increases). If this field is empty,
        #' this object describes a contract block (calculated as a stand-alone
        #' tariff), otherwise it will simply be the sum of its blocks (adjusted
        #' to span the same time periods)
        blocks = list(),

        #' @field history
        #' A list keeping track of all contract changes (including the whole
        #' contract state and its values before the change).
        history = list(),


        #### The code:

        #' @description  Create a new insurance contract (for the given tariff/product) and calculate all time series
        #'
        #' @details The \code{InsuranceContract$new()} function creates a new
        #' insurance contract for the given tariff, using the parameters passed
        #' to the function (and the defaults specified in the tariff).
        #'
        #' As soon as this function is called, the contract object calculates
        #' all time series (cash flows, premiums, reserves, profit participation)
        #' for the whole contract duration.
        #'
        #' The most important parameters that are typically passed to the
        #' constructor are:
        #' * \code{age} ... Age of the insured person (used to derive mortalities / transition probabilities)
        #' * \code{policyPeriod} ... Maturity of the policy (in years)
        #' * \code{premiumPeriod} ... How long premiums are paid (\code{premiumPeriod = 1}
        #'        for single-premium contracts, \code{premiumPeriod} equals
        #'        \code{policyPeriod} for regular premium payments for the whole
        #'        contract period, while other premium payment durations indicate
        #'        premium payments only for shorter periods than the whole contract
        #'        duration. Contract extensions without any premium payments are
        #'        indicated by \code{premiumPeriod}=0). Default is equal to
        #'        \code{policyPeriod}
        #' * \code{sumInsured} ... The sum insured (i.e. survival benefit for
        #'         endowments, death benefit for whole/term life insurances,
        #'         annuity payments for annuities)
        #' * \code{contractClosing} ... Date of the contract beginning (typically
        #'        created using something like \code{as.Date("2020-08-01")})
        #' * \code{YOB} ... Year of birth of the insured (for cohort mortality
        #'        tables). If not given, YOB is derived from \code{age} and
        #'        \code{contractClosing}.
        #' * \code{deferralPeriod} ... Deferral period for deferred annuities
        #'        (i.e. when annuity payments start at a future point in time).
        #'        Default is 0.
        #' * \code{premiumFrequency} ... How many premium payments per year are
        #'        made (e.g. 1 for yearly premiums, 4 for quarterly premiumd,
        #'        12 for monthly premium payments). Default is 1 (yearly premiums).
        #'
        #' While these are the most common and most important parameters, all
        #' parameters can be overwritten on a per-contract basis, even those
        #' that are defined by the tariff. For a full list and explanation of all
        #' parameters, see [InsuranceContract.ParameterDefaults].
        #'
        #' @param tarif The [InsuranceTarif] object describing the Tariff/Product
        #'        and providing defaults for the parameters.
        #' @param parent For contracts with multiple contract blocks (dynamic
        #'        increases, sum increases, riders), each child is created with
        #'        a pointer to its parent. NULL for single-block contracts or
        #'        for the overall-contract of a multi-block contract. This
        #'        parameter is used internally, but should not be used in
        #'        user-written code.
        #' @param calculate how much of the contract's time series need to be
        #'        calculated. See [CalculationEnum] for all possible values. This
        #'        is usefull to prevent calculation of e.g. reserves and profit
        #'        participation, when one only wants to create a grid of premiums.
        #' @param profitid The ID of the default profit participation scenario.
        #'        The default profit participation scenario uses the default
        #'        values passed, while further scenarios can be added by
        #'        \ifelse{html}{\href{#method-addProfitScenario}{\code{InsuranceContract$addProfitScenario()}}}{\code{InsuranceContract$addProfitScenario()()}}.
        #' @param ... Further parameters (age, sum insured, contract closing /
        #'        begin, premium payment details, etc.) of the contract, which
        #'        can also override parameters defined at the tariff-level.
        #'        Possible values are all sub-fields of the
        #'        [InsuranceContract.ParameterDefaults] data structure.
        #'
        initialize = function(tarif, parent = NULL, calculate = "all", profitid = "default", ...) {
            if (getOption('LIC.debug.init', FALSE)) {
                browser();
            }
            private$initParams = c(list(tarif = tarif, parent = parent, calculate = calculate, profitid = profitid), list(...))
            self$tarif = tarif;
            self$parent = parent;

            # TODO-block: If parent exists, use its parameters as fallback, too
            self$ContractParameters = InsuranceContract.ParametersFill(
                ...,
                premiumWaiver = FALSE,
                surrenderPenalty = TRUE,
                alphaRefunded = FALSE
            );
            # Set default values for required contract-specific data
            # First, take the tariff defaults, then the  ProfitParticipation
            # defaults, so a tariff can override the profit participation scheme
            self$Parameters = self$ContractParameters;
            self$Parameters = InsuranceContract.ParametersFallback(
                self$ContractParameters,
                self$tarif$getParameters()
            );

            ppScheme = self$Parameters$ProfitParticipation$profitParticipationScheme;
            if (!is.null(ppScheme)) {
                self$Parameters$ProfitParticipation = fallbackFields(
                    self$Parameters$ProfitParticipation,
                    ppScheme$Parameters);
                self$Parameters$ProfitParticipation$scenarios[[profitid]] = list()
            }

            private$consolidateContractData(tarif = tarif, ...);
            self$calculateContract(calculate = calculate);


            invisible(self)
        },

        #' @description Add the current state of the contract to the history list
        #'
        #' @details The \code{InsuranceContract$addHistorySnapshot()} function
        #' adds the current (or the explicitly given) state of the contract
        #' (parameters, calculated values, tarif, list of all contract blocks)
        #' to the history list (available in the \code{history} field of the
        #' contract, i.e. \code{InsuranceContract$history}).
        #'
        #' @param time the time described by the snapshot
        #' @param comment a comment to store together with the contract state
        #' @param type The type of action that caused a history snapshot to
        #'        be stored. Typical values are "Contract" to describe the initial
        #'        contract, "Premium Waiver" or "Dynamic Increase".
        #' @param params The set of params to be stored in the history snapshot
        #'        (default is \code{self$Parameters}, if not explicitly given)
        #' @param values The calculated time series of all contract values
        #'        calculated so far. Default is \code{self$Values}, if not
        #'        explicitly given
        #' @param tarif The underlying [InsuranceTarif] object describing the
        #'        Product/Tariff. Default is \code{self$tarif}, if not explicitly given.
        #' @param blocks The list of all contract children for contracts with
        #'        multiple insurance blocks (e.g. dynamic increases, riders, etc.)
        #'
        #' @examples
        #' # TODO
        addHistorySnapshot = function(time = 0, comment = "Initial contract values", type = "Contract", params = self$Parameters, values = self$Values, tarif = self$tarif, blocks = self$blocks) {
            self$history = rbind(
                self$history,
                list(
                    time = list(
                        "time"    = time,
                        "comment" = comment,
                        "type"    = type,
                        "params"  = params,
                        "values"  = values,
                        "tarif"   = tarif,
                        "blocks"  = blocks
                    )
                )
            );
            invisible(self)
        },

        #' @description Add a child contract block (e.g. a dynamic increase or a rider) to an insurance contract
        #'
        #' @details Contracts with multiple contract blocks (typically either
        #' contracts with dynamic increases, sum increases or protection riders)
        #' are constructed by instantiating the child block (e.g. a single
        #' dynamic increase or the rider) independently with its own (shorter)
        #' duration and then inserting it into the parent contract with this
        #' function at the given time.
        #'
        #' If no [InsuranceContract] object is passed as \code{block}, a copy
        #' of the parent is created with overriding parameters given in \code{...}.
        #'
        #' @param id The identifier of the child block to be inserted
        #' @param block The [InsuranceContract] object describing the child block.
        #'        If NULL (or not given at all), a copy of the parent will be
        #'        created.
        #' @param t Then the child block starts, relative to the parent block.
        #'        The child block is calculated independently (with time 0
        #'        describing its own start), so when aggregating all values from
        #'        the individual blocks to overall values for the whole contract,
        #'        the child's values need to be translated to the parent contracts's
        #'        time frame using this parameter
        #' @param comment The comment to use in the history snapshot.
        #' @param blockType The type of block to be added (e.g. Dynamics, Extension,
        #'        etc.). Can be any (short) string.
        #' @param ... parameters to be passed to \ifelse{html}{\href{#method-new}{\code{InsuranceContract$new()}}}{\code{InsuranceContract$new()()}} when
        #'        \code{block} is not given and a copy of the parent should be
        #'        created with overrides.
        #'
        #' @examples
        #' # TODO
        addBlock = function(id = NULL, block = NULL, t = block$Values$int$blockStart, comment = paste0("Additional block at time t=", t), blockType = "Dynamics", ...) {
            if (getOption('LIC.debug.addBlock', FALSE)) {
                browser();
            }
            if (missing(id) || is.null(id)) {
                id = paste0("block", length(self$blocks) + 1)
            }
            if (missing(block) || is.null(block) || !is(block, "InsuranceContract")) {
                # Create a block with the same tariff and parameters as the main contract, but allow overriding params with the ... arguments
                block = InsuranceContract$new(id = id, ...)
            } else {
                block$Parameters$ContractData$id = id
            }
            if (missing(t) || is.null(t)) {
                t = 0
            }
            # declare as child of 'self', store the time offset to the parent contract
            block$parent = self
            block$Parameters$ContractData$blockStart = t

            if (length(self$blocks) == 0) {
                main = self$clone()
                main$parent = self
                self$blocks[[main$Parameters$ContractData$id]] = main
                self$Parameters$ContractData$id = "Gesamt"
            }

            self$blocks[[id]] = block
            # recalculate the whole contract by consolidating values from each block
            self$consolidateBlocks(valuesFrom = t)

            self$addHistorySnapshot(time = t, comment = comment,
                                    type = blockType, params = self$Parameters, values = self$Values);

            invisible(self)
        },

        #' @description Add a dynamic increase with the same parameters as the main contract part
        #'
        #' @details This method adds a new contract block describing a dynamic
        #'          or sum increase (increasing the sum insured at a later time
        #'          $t$ than contract inception). This increase is modelled by a
        #'          separate [InsuranceContract] object with the sum difference
        #'          as its own \code{sumInsured}.
        #'
        #'          By default, all parameters are taken from the main contract,
        #'          with the maturity adjusted to match the original contract's
        #'          maturity.
        #'
        #'          The main contract holds all child blocks, controls their
        #'          valueation and aggregates all children's values to the
        #'          total values of the overall contract.
        #'
        #' @param t The time within the main contract when the sum increase happens.
        #'          The [InsuranceContract] object describing the dynamic increase
        #'          will still internally start at its own time 0, but the
        #'          aggregation by the main contract will correctly offset to
        #'          the time $t$ within the main contract.
        #' @param NewSumInsured The over-all new sum insured (sum of original
        #'          contract and all dynamica increaeses). The \code{sumInsured}
        #'          of the new dynamic increase block will be determined as the
        #'          difference of the old and new overall sum insured. Alternatively,
        #'          it can directly be given as the \code{SumInsuredDelta}
        #'          argument instead.
        #' @param SumInsuredDelta The sum insured of only the dynamic increase,
        #'          i.e. the sumInsured of the dynamic contract block only. The
        #'          overall sum insured will increase by this amount. Only one of
        #'          \code{NewSumInsured} and \code{SumInsuredDelta} is needed,
        #'          the other one will be calculated accordingly. If both are
        #'          given, the \code{SumInsuredDelta} will take precedence.
        #' @param id The identifier of the contract block describing the dynamic
        #'          increase. This is a free-form string that should be unique
        #'          within the list of child blocks. It will be displayed in the
        #'          'Excel' export feature and in the history snapshot list.
        #' @param ... Parameters to override in the dynamic block. By default,
        #'          all parameters of the main contract block will be used, but
        #'          they can be overridden per dynamic increase block.
        #'
        #' @examples
        #' # TODO
        addDynamics = function(t, NewSumInsured, SumInsuredDelta, id, ...) {
            if (getOption('LIC.debug.addDynamics', FALSE)) {
                browser();
            }

            # TODO: Override only the required parameters
            params = private$initParams
            if (is.null(params)) params = list()
            if (!is.null(params$age)) params$age = params$age + t
            if (!is.null(params$policyPeriod)) params$policyPeriod = params$policyPeriod - t
            if (!is.null(params$premiumPeriod)) params$premiumPeriod = max(0, params$premiumPeriod - t)
            if (!is.null(params$deferralPeriod)) params$deferralPeriod = max(0, params$deferralPeriod - t)
            if (!is.null(params$contractClosing)) params$contractClosing = params$contractClosing + years(t)
            params$initialCapital = NULL
            # TODO: Adjust non-constant parameters (e.g. profit rates or benefits given as vector) to the later start time

            # TODO: Generalize this to also allow specifying dynamic premium rather than sum insured
            if (!missing(SumInsuredDelta)) {
                SIdelta = SumInsuredDelta
                NewSumInsured = SIdelta + self$Values$reserves[t + 1, "SumInsured"]
            } else if (!missing(NewSumInsured)) {
                SIdelta = NewSumInsured - self$Values$reserves[t + 1, "SumInsured"]
            } else {
                warning("Neither NewSumInsured nor SumInsuredDelta are given. Unable to determine sum insured of the dynamic increase.")
                return(invisible(self))
            }
            params$sumInsured = SIdelta

            if (missing(id)) {
                # numbering dynamics: Use # of blocks (except main) as a simplification
                id = paste0("dyn", max(1, length(self$blocks)))
            }

            params$t = t
            params$id = id
            # Override with arguments explicitly given
            arguments = list(...)
            params[names(arguments)] = arguments[names(arguments)]
            params$comment = sprintf("Dynamic increase at time %d to sum %0.2f", t, NewSumInsured)
            params$blockType = "Dynamics";
            do.call(self$addBlock, params)
        },


        #' @description Add a contract extension after the contract has ended
        #' (existing reserve is used as initial capital of the follow-up contract).
        #'
        #' @details When a contract expires, this function adds a follow-up contract
        #' (with either the same or a different tariff), using the existing
        #' reserve as `additionalCapital` at inception.
        #' Technically, a child block using the new contract data of the extension
        #' is added to the original contract. The over-all contract values are then
        #' the sum of the original contract (providing values until expiration)
        #' and the extension (providing values after the extension).
        #'
        #'
        #' @param id The identifier of the child block to be inserted
        #' @param t The time of the extension (relative to the parent block),
        #'        by default contract expiration of the parent block.
        #'        The extension is calculated independently (with time 0
        #'        describing its own start), but using the existing reserve as
        #'        initialCapital and the parent's parameters as fall-back values.
        #' @param comment The comment to use in the history snapshot.
        #' @param ... Additional parameters to be passed to
        #'        \ifelse{html}{\href{#method-new}{\code{InsuranceContract$new}}}{\code{InsuranceContract$new}} to create the contract
        #'        extension object.
        #'
        #' @examples
        #' # TODO
        addExtension = function(id = NULL, t = NULL, comment = paste0("Contract extension at time t=", t), ...) {
            if (getOption('LIC.debug.addExtension', FALSE)) {
                browser();
            }
            if (missing(id) | is.null(id)) {
                # numbering extensions: Use nr. of blocks (except main) as a
                # simplification => numbering is withint all dynamics,
                # extensions, riders, etc.!
                id = paste0("dyn", max(1, length(self$blocks)))
            }
            if (missing(t) | is.null(t)) {
                # By default, use the parent's expiration, so the extension
                # is appended after the original contract has ended.
                t = self$Parameters$ContractData$policyPeriod
            }

            # TODO: Override only the required parameters
            params = private$initParams
            if (is.null(params)) params = list()
            if (!is.null(params$age)) params$age = params$age + t
            # Remaining premium period is kept, can be overwritten in the
            # arguments to this method. If premiumPeriod has already ended,
            # a premium-free extension is added by default
            if (!is.null(params$premiumPeriod)) params$premiumPeriod = max(0, params$premiumPeriod - t)
            if (!is.null(params$deferralPeriod)) params$deferralPeriod = max(0, params$deferralPeriod - t)
            if (!is.null(params$contractClosing)) params$contractClosing = params$contractClosing + years(t)
            # Use the existing reserve as initialCapital, reset premium parameter and sumInsured of the old contract
            params$initialCapital = self$Values$reserves[t + 1, "contractual"]
            params$sumInsured = NULL
            params$premium = NULL


            # TODO: Adjust non-constant parameters (e.g. profit rates or benefits given as vector) to the later start time


            params$t = t
            params$id = id
            # Override with arguments explicitly given
            arguments = list(...)
            params[names(arguments)] = arguments[names(arguments)]
            params$comment = comment
            params$blockType = "Contract Extension"

            # Two cases:
            # 1) extension with premium => either premium or sumInsured must be given
            # 2) premium-free extension => no premium, no sumInsured given => set premium to 0, calculate sumInsured

            noPremiums = (is.null(params$sumInsured) || (params$sumInsured == 0)) && (is.null(params$premium) || (params$premium == 0));
            if (noPremiums) {
                params$premium = 0;
            }
            do.call(self$addBlock, params)
        },

        #' @description Calculate all time series of the contract from the parameters
        #'
        #' @details This method calculates all contract values (potentially
        #'          starting from and preserving all values before a later time
        #'          \code{valuesFrom}). This function is not meant to be called
        #'          directly, but internally, whenever a contract is created or
        #'          modified.
        #'
        #'          There is, however, a legitimate case to call this function
        #'          when a contract was initially created with a value of
        #'          \code{calculate} other than "all", so not all values of the
        #'          contract were calculated. When one later needs more values
        #'          than were initially calculated, this function can be called.
        #'          However, any contract changes might need to be rolled back
        #'          and reapplied again afterwards. So even in this case it is
        #'          probably easier to create the contract object from scratch
        #'          again.
        #' @param calculate Which values to calculate. See [CalculationEnum]
        #' @param valuesFrom Calculate only values starting from this time step
        #'        on (all values before that time will be preserved). This is
        #'        required when a contract is changed significantly (potentially
        #'        even switching to a new tariff), so that the calculation bases
        #'        for previous periods are no longer available.
        #' @param premiumCalculationTime The time point when the premium should
        #'        be re-calculated (including existing reserves) based on the
        #'        actuarial equivalence principle. All reserves will be based on
        #'        these new premiums.
        #' @param preservePastPV Whether present value before the recalculation
        #'        time \code{valuesFrom} should be preserved or recalculated.
        #'        When they are recalculated, the present values are consistent
        #'        to the new cash flows over the whole contract period, but
        #'        they no longer represent the actual contract state at these
        #'        times. If values are not recalculated, the reserves at each
        #'        time step represent the proper state at that point in time.
        #' @param additionalCapital The capital that is added to the contract
        #'        (e.g. capital carried over from a previous contract) at the
        #'        premium calculation time.
        #' @param recalculatePremiums Whether the premiums should be recalculated
        #'        at time \code{premiumCalculationTime} at all.
        #' @param recalculatePremiumSum Whether to recalculate the overall premium
        #'        sum when the premium is recalculated.
        #' @param history_comment The comment for the history snapshot entry
        #' @param history_type The type (free-form string) to record in the history snapshot
        #'
        calculateContract = function(calculate = "all", valuesFrom = 0, premiumCalculationTime = 0, preservePastPV = TRUE, additionalCapital = 0, recalculatePremiums = TRUE, recalculatePremiumSum = TRUE, history_comment = NULL, history_type = "Contract") {
            if (getOption('LIC.debug.calculateContract', FALSE)) {
                browser();
            }
            if (!is.null(self$blocks)) {
                for (b in self$blocks) {
                    #
                    # .args = as.list(match.call()[-1])
                    # # correctly shift the valuesFrom by each block's blockStart parameter
                    # .args$valuesFrom = max(0, .args$valuesFrom - b$Parameters$ContractData$blockStart)
                    # do.call(b$calculateContract, .args)
                    #
                    b$calculateContract(
                        calculate = calculate,
                        valuesFrom = max(0, valuesFrom - b$Parameters$ContractData$blockStart),
                        premiumCalculationTime = max(0, premiumCalculationTime - b$Parameters$ContractData$blockStart),
                        preservePastPV = preservePastPV,
                        additionalCapital = additionalCapital,
                        recalculatePremiums = recalculatePremiums,
                        recalculatePremiumSum = recalculatePremiumSum,
                        history_comment = history_comment,
                        history_type = history_type)
                }
            }
            self$Values$int = private$determineInternalValues()
            self$Values$transitionProbabilities = mergeValues(
              starting = self$Values$transitionProbabilities,
              ending = private$determineTransitionProbabilities(),
              t = valuesFrom)
            if (calculate == "probabilities") return(invisible(self));

            self$Values$cashFlowsBasic = mergeValues(
                starting = self$Values$cashFlowsBasic,
                ending = private$determineCashFlowsBasic(),
                t = valuesFrom);
            self$Values$cashFlows = mergeValues(
                starting = self$Values$cashFlows,
                ending = private$determineCashFlows(),
                t = valuesFrom);

            if (additionalCapital > 0) {
                self$Values$cashFlows[as.character(premiumCalculationTime), "additional_capital"] = additionalCapital
            }

            if (recalculatePremiumSum) {
                # Premium waiver: Premium sum is not affected by premium waivers, i.e. everything depending on the premium sum uses the original premium sum!
                self$Values$unitPremiumSum = private$determinePremiumSum();
            }
            self$Values$cashFlowsCosts = mergeValues3D(
                starting = self$Values$cashFlowsCosts,
                ending = private$determineCashFlowsCosts(),
                t = valuesFrom);
            if (calculate == "cashflows") return(invisible(self));


            # Shall we re-calculate PV or preserve the old ones???
            pv = private$calculatePresentValues()
            pvCost = private$calculatePresentValuesCosts(presentValues = pv)
            oldPV = self$Values$presentValues
            if (preservePastPV) {
                # Preserve past present values, i.e. the PV represents the PV
                # with the knowledge of the past, even though the future CF
                # might have changed meanwhile, so the PV at time 0 is no
                # longer the PV of the current cash flows... The PV at time t
                # always represents the information available at time t, but no
                # future chagnes.
                # This is useful to preserver the PV information neede to
                # calculate the premiums from the past.
                if (!is.null(self$Values$presentValues)) {
                    self$Values$presentValues = self$Values$presentValues[,1:NCOL(pv)]
                }
                self$Values$presentValues = mergeValues(starting = self$Values$presentValues, ending = pv, t = valuesFrom)
                self$Values$presentValuesCosts = mergeValues3D(starting = self$Values$presentValuesCosts, ending = pvCost, t = valuesFrom)
            } else {
                # Recalculate present value for times before start, i.e. make all PV consistent with the current cash flows
                self$Values$presentValues = pv
                self$Values$presentValuesCosts = pvCost
            }
            if (calculate == "presentvalues") return(invisible(self));


            # If we have the premium given, determine the sumInsured from it
            # Since the cash flows depend on the sumInsured (in particular, )
            if (is.null(self$Parameters$ContractData$sumInsured)) {
                self$Parameters$ContractData$sumInsured = private$calculateSumInsured(calculationTime = premiumCalculationTime)
            }

            # the premiumCalculation function returns the premiums AND the cofficients,
            # so we have to extract the coefficients and store them in a separate variable
            if (recalculatePremiums) {
                res = private$calculatePremiums(premiumCalculationTime = premiumCalculationTime);
                self$Parameters$ContractData$sumInsured = res[["sumInsured"]]
                self$Values$premiumCoefficients = res[["coefficients"]];
                # TODO: Store premiums in a data.frame, including the time they are calculated???
                self$Values$premiums = res[["premiums"]]
                self$Values$int$premiumCalculationTime = premiumCalculationTime
            }
            if (calculate == "premiums") return(invisible(self));

            # Update the cash flows and present values with the values of the premium
            pvAllBenefits = private$calculatePresentValuesBenefits()
            if (preservePastPV) {
                self$Values$presentValues = mergeValues(starting = oldPV, ending = cbind(pv, pvAllBenefits), t = valuesFrom)
            } else {
                self$Values$presentValues = cbind(pv, pvAllBenefits)
            }

            self$Values$absCashFlows       = mergeValues(starting = self$Values$absCashFlows,       ending = private$calculateAbsCashFlows(), t = valuesFrom);
            self$Values$absPresentValues   = mergeValues(starting = self$Values$absPresentValues,   ending = private$calculateAbsPresentValues(), t = valuesFrom);
            if (calculate == "absvalues") return(invisible(self));

            self$Values$reserves           = mergeValues(starting = self$Values$reserves,           ending = private$calculateReserves(), t = valuesFrom);
            if (calculate == "reserves") return(invisible(self));
            self$Values$premiumComposition = mergeValues(starting = self$Values$premiumComposition, ending = private$premiumAnalysis(), t = valuesFrom);
            self$Values$premiumCompositionSums = mergeValues(starting = self$Values$premiumCompositionSums, ending = private$premiumCompositionSums(), t = valuesFrom);
            self$Values$premiumCompositionPV = mergeValues(starting = self$Values$premiumCompositionPV, ending = private$premiumCompositionPV(), t = valuesFrom);
            self$Values$basicData          = mergeValues(starting = self$Values$basicData,          ending = private$getBasicDataTimeseries(), t = valuesFrom);
            if (calculate == "premiumcomposition") return(invisible(self));

            self$Values$reservesBalanceSheet = mergeValues(starting = self$Values$reservesBalanceSheet,ending = private$calculateReservesBalanceSheet(), t = valuesFrom);
            if (calculate == "reservesbalancesheet") return(invisible(self));

            private$profitParticipation(calculateFrom = valuesFrom);
            if (calculate == "profitparticipation") return(invisible(self));

            self$addHistorySnapshot(
                time    = valuesFrom,
                comment = ifelse(is.null(history_comment),
                                 ifelse(valuesFrom == 0, "Initial contract values", paste("Contract recalculation at time", premiumCalculationTime)),
                                 history_comment),
                type    = history_type,
                params  = self$Parameters,
                values  = self$Values
            );
            if (calculate == "history") return(invisible(self));

            invisible(self)
        },

        #' @description Aggregate values from all child contract blocks (if any)
        #'
        #' @details This function is an internal function for contracts with
        #'        multiple child blocks (dynamic increases, sum increases, riders).
        #'        It takes the values from all child blocks and calculates the
        #'        overall values from all child blocks aggregated.
        #'
        #'        This function should not be called manually.
        #'
        #' @param valuesFrom The time from when to aggragate values. Values before
        #'        that time will be left unchanged.
        consolidateBlocks = function(valuesFrom = 0) {
            if (getOption('LIC.debug.consolidateBlocks', FALSE)) {
                browser();
            }
            # First, Re-calculate all children that have children on their own
            for (b in self$blocks) {
                if (length(b$blocks) > 0) {
                    # correctly shift the valuesFrom by each block's blockStart parameter
                    b$consolidateBlocks(valuesFrom = max(0, valuesFrom - b$Parameters$ContractData$blockStart))
                }
            }

            sumKeyedArrays = function(arr1 = NULL, arr2 = NULL) {
                if (is.null(arr2)) {
                    arr1
                } else if (is.null(arr1)) {
                    arr2
                } else {
                    bind_rows(arr1, arr2) %>%
                        group_by(date) %>%
                        summarise_all(list(sum))
                }
            }
            consolidateField = function(field, keyed = FALSE) {
                vals = NULL
                if (length(self$blocks) == 0) {
                    vals = self$Values[[field]]
                }
                for (b in self$blocks) {
                    if (keyed) {
                        # The rows of the two data.frames can be associated by the values of a certain column
                        vals = sumKeyedArrays(arr1 = vals, arr2 = b$Values[[field]])
                    } else {
                        # Simply pad the arrays and sum them up:
                        vals = sumPaddedArrays(arr1 = vals, arr2 = b$Values[[field]], pad2 = b$Parameters$ContractData$blockStart)
                    }
                }
                mergeValues(starting = self$Values[[field]],   ending = vals, t = valuesFrom);
            }


            # Some values do not make sense for consolidated values
            self$Values$cashFlowsBasic = NULL
            self$Values$cashFlows = NULL
            self$Values$cashFlowsCosts = NULL
            self$Values$presentValues = NULL
            self$Values$presentValuesCosts = NULL

            self$Values$unitPremiumSum = NULL
            self$Values$premiumCoefficients = NULL
            self$Values$premiums = NULL

            # self$Values$transitionProbabilities = consolidateField("transitionProbabilities")
            # self$Values$cashFlowsBasic       = consolidateField("cashFlowsBasic")
            # self$Values$cashFlows            = consolidateField("cashFlows")
            # self$Values$cashFlowsCosts       = consolidateField("cashFlowsCosts")
            self$Values$absCashFlows           = consolidateField("absCashFlows")
            self$Values$absPresentValues       = consolidateField("absPresentValues")
            self$Values$premiumComposition     = consolidateField("premiumComposition")
            self$Values$premiumCompositionSums = consolidateField("premiumCompositionSums")
            self$Values$premiumCompositionPV   = consolidateField("premiumCompositionPV")
            self$Values$reserves               = consolidateField("reserves")
            self$Values$reservesBalanceSheet   = consolidateField("reservesBalanceSheet", keyed = TRUE)
            # TODO: Basic Data cannot simply be summed, e.g. the interest rate!
            self$Values$basicData              = consolidateField("basicData")

            # Some fields can NOT be summed, but have to be left untouched.
            # Hard-code these to use the values from the main contract part:
            rows = nrow(self$Values$reservesBalanceSheet)
            colDt = rep(as.Date(NA), rows)
            colTime = rep(NA_real_, rows)
            colIntR = rep(NA_real_, rows)
            colDur = rep(NA_real_, rows)
            colPrem = rep(NA_real_, rows)
            polPeriod = self$Parameters$ContractData$policyPeriod
            for (b in self$blocks) {
                start = b$Parameters$ContractData$blockStart
                polPeriod = max(polPeriod, b$Parameters$ContractData$policyPeriod + start)
                colDt = coalesce(colDt, pad0(b$Values$reservesBalanceSheet[,"date"],     start = start, value = as.Date(NA), value.start = as.Date(NA), l = rows))
                colTime = coalesce(colTime, pad0(b$Values$reservesBalanceSheet[,"time"] + start, start = start, value = NA, value.start = NA, l = rows))

                colIntR = coalesce(colIntR, pad0(b$Values$basicData[,"InterestRate"],    start = start, value = NA, value.start = NA, l = rows))
                colDur  = coalesce(colDur,  pad0(b$Values$basicData[,"PolicyDuration"],  start = start, value = NA, value.start = NA, l = rows))
                colPrem = coalesce(colPrem, pad0(b$Values$basicData[,"PremiumPeriod"],   start = start, value = NA, value.start = NA, l = rows))
            }
            self$Parameters$ContractData$policyPeriod = polPeriod;
            self$Values$reservesBalanceSheet[,"date"] = colDt;
            self$Values$reservesBalanceSheet[,"time"] = colTime;
            self$Values$basicData[,"InterestRate"] = colIntR
            self$Values$basicData[,"PolicyDuration"] = colDur
            self$Values$basicData[,"PremiumPeriod"] = colPrem

            self$Values$int$l = rows

            invisible(self)
        },

        #' @description Stop premium payments and re-calculate sumInsured of the paid-up contract
        #'
        #' @details This function modifies the contract at time $t$ so that
        #'        no further premiums are paid (i.e. a paid-up contract) and the
        #'        \code{sumInsured} is adjusted according to the existing reserves.
        #'
        #' @param t Time of the premium waiver.
        #' @param ... Further parameters (currently unused)
        #'
        #' @examples
        #' # TODO
        premiumWaiver = function(t, ...) {
            if (getOption('LIC.debug.premiumWaiver', FALSE)) {
                browser();
            }
            if (length(self$blocks) > 0) {
                for (b in self$blocks) {
                    b$premiumWaiver(t - b$Parameters$ContractData$blockStart, ...)
                }
            } else {
                newSumInsured = self$Values$reserves[[toString(t), "PremiumFreeSumInsured"]];
                self$Parameters$ContractData$sumInsured = newSumInsured;
            }
            self$Parameters$ContractState$premiumWaiver = TRUE;
            self$Parameters$ContractState$surrenderPenalty = FALSE; # Surrender penalty has already been applied, don't apply a second time
            self$Parameters$ContractState$alphaRefunded = TRUE;     # Alpha cost (if applicable) have already been refunded partially, don't refund again
            # TODO: Extract current amount of premium refund and feed that into the calculateContract function...

            self$calculateContract(
                valuesFrom = t,
                premiumCalculationTime = t,
                preservePastPV = TRUE, recalculatePremiums = TRUE, recalculatePremiumSum = FALSE,
                history_comment = sprintf("Premium waiver at time %d", t), history_type = "PremiumWaiver")

            invisible(self)
        },



        #' @description Calculate one profit scenario and return all values
        #'
        #' @details This function calculates one profit scenario with the
        #'         provided profit participation parameters (all parameters
        #'         not given in the call are taken from their values of the
        #'         contract, profit participation scheme or tariff).
        #'
        #' @param ... Scenario-specific profit sharing parameters, overriding
        #'         the default values. Typically, adjusted profit rates are required
        #'         in a profitScenario.
        #' @return a data.frame holding all profit participation values (rates,
        #'         bases for the different profit types, profit allocations,
        #'         terminal bonus funds, profit in case of death/surrender/premium waiver)
        #'
        #' @examples
        #' # TODO
        profitScenario = function(...) {
            if (getOption('LIC.debug.profitScenario', FALSE)) {
                browser();
            }
            private$calculateProfitParticipation(...)
        },

        #' @description Calculate one profit scenario and store it in the contract
        #'
        #' @details This function calculates one profit scenario with the
        #'         provided profit participation parameters (all parameters
        #'         not given in the call are taken from their values of the
        #'         contract, profit participation scheme or tariff). The results
        #'         are stored in a list of profit scenarios inside the contract.
        #'
        #'         This function can be chained to calculate and add multiple
        #'         profit scenarios.
        #'
        #' @param id The unique ID of the profit scenario. Will be used as key
        #'         in the list of profit scenarios and printed out in the 'Excel'
        #'         export.
        #' @param ... Scenario-specific profit sharing parameters, overriding
        #'         the default values. Typically, adjusted profit rates are required
        #'         in a profitScenario.
        #'
        #' @examples
        #' # TODO
        addProfitScenario = function(id, ...) {
            if (getOption('LIC.debug.addProfitScenario', FALSE)) {
                browser();
            }
            .args = as.list(match.call()[-1])
            self$Parameters$ProfitParticipation$scenarios[[id]] = list(...)
            if (length(self$blocks) > 0) {
                vals = NULL
                for (b in self$blocks) {
                    # TODO: shift the profit rates by b$Parameters$ContractData$blockStart
                    do.call(b$addProfitScenario, .args)
                    vals = sumPaddedArrays(arr1 = vals, arr2 = b$Values$profitScenarios[[id]], pad2 = b$Parameters$ContractData$blockStart)
                    # TODO: consolidate reserves after profit!
                }
                # Consolidate all child blocks
                self$Values$profitScenarios[[id]] = vals
            } else {
                # profitParticipation will assign the values to Values$profitScenarios[[id]] and Values$reservesAfterProfit[[id]]
                # private$profitParticipation(id = id, ...)
                pp = private$calculateProfitParticipation(...)
                self$Values$profitScenarios[[id]] = pp
                self$Values$reservesAfterProfit[[id]] = private$calculateReservesAfterProfit(profitScenario = pp, ...)
            }

            invisible(self)
        },


        #' @field dummy.public
        #' dummy field to allow a trailing comma after the previous field/method
        dummy.public = NULL
    ),

    ######################### PRIVATE METHODS ##################################
    private = list(
        initParams = NULL,

        consolidateContractData = function(...) {
            args = list(...);
            # TODO-blocks

            if (getOption('LIC.debug.consolidateContractData', FALSE)) {
                browser();
            }
            # YOB is deprecated in favor of birthDate. If it is given, use January 1
            if (is.null(self$Parameters$ContractData$birthDate) && !is.null(self$Parameters$ContractData$YOB)) {
                self$Parameters$ContractData$birthDate = make_date(self$Parameters$ContractData$YOB, 1, 1)
            }

            # Calculate date/year of birth, age, contract closing etc. from each other
            # 1. Contract date (if not given) is NOW, unless age + YOB is given => Then year is derived as YOB+age
            if (is.null(self$Parameters$ContractData$contractClosing)) {
                # Default is contractClosing is NOW:
                self$Parameters$ContractData$contractClosing = Sys.Date()

                # However, if age and DOB / YOB is given, calculate contract closing from that:
                # age is given (and not a function that calculates age from DOB and Contract closing)
                if (!is.null(self$Parameters$ContractData$age) &&
                    !is.function(self$Parameters$ContractData$age)
                ) {
                    if (!is.null(self$Parameters$ContractData$birthDate)) {
                        ag = self$Parameters$ContractData$age
                        # Whole years are added as period (so the day stays the same), remaining fractions are added as dyears
                        self$Parameters$ContractData$contractClosing = as.Date(self$Parameters$ContractData$birthDate +
                            years(floor(self$Parameters$ContractData$age)) +
                            dyears(self$Parameters$ContractData$age %% 1))
                        # TODO: Always start at the nearest beginning of a month? Or leave the contract closing at any day?
                    }
                }
            }

            # 2. Current age: If age is given, use it
            if (!is.null(self$Parameters$ContractData$age)) {
                self$Parameters$ContractData$age = valueOrFunction(
                    self$Parameters$ContractData$age,
                    params = self$Parameters, values = self$Values);
            } else {
            # 3. Otherwise, either use the birth date to calculate the age
                if (!is.null(self$Parameters$ContractData$birthDate)) {
                    # TODO: Decide for variant 1 or 2...
                    # Variant 1: Exact age rounded to the nearest whole number
                    self$Parameters$ContractData$age = age.exactRounded(self$Parameters, self$Values)
                    # Variant 2: Year of contract closing - YOB
                    self$Parameters$ContractData$age = age.yearDifference(self$Parameters, self$Values)
                } else {
            # 4. Or use age=40 as default
                    self$Parameters$ContractData$age = 40
                    warning("InsuranceContract: Missing age, no information to derive age from YOB and contractClosing => Assuming default age 40. Tariff: ", self$tarif$name)
                }
            }
            if (is.null(self$Parameters$ContractData$birthDate)) {
                self$Parameters$ContractData$birthDate = as.Date(self$Parameters$ContractData$contractClosing -
                    years(floor(self$Parameters$ContractData$age)) -
                    dyears(self$Parameters$ContractData$age %% 1))
            }



            # Evaluate policy period, i.e. if a function is used, calculate its numeric value
            self$Parameters$ContractData$policyPeriod = valueOrFunction(
                self$Parameters$ContractData$policyPeriod,
                params = self$Parameters, values = self$Values);

            self$Parameters$ContractData$guaranteedPeriod = valueOrFunction(
                self$Parameters$ContractData$guaranteedPeriod,
                params = self$Parameters, values = self$Values);

            #### #
            # PREMIUM PAYMENT PERIOD (default: policyPeriod, can be given as function or numeric value)
            #### #
            if (is.null(self$Parameters$ContractData$premiumPeriod)) {
                self$Parameters$ContractData$premiumPeriod = self$Parameters$ContractData$policyPeriod
            }
            self$Parameters$ContractData$premiumPeriod = valueOrFunction(
                self$Parameters$ContractData$premiumPeriod,
                params = self$Parameters, values = self$Values);
            # premium period is at most contract duration!
            self$Parameters$ContractData$premiumPeriod =
                min(
                    self$Parameters$ContractData$premiumPeriod,
                    self$Parameters$ContractData$policyPeriod
                );

            self$Parameters$Loadings$commissionPeriod = valueOrFunction(
                self$Parameters$Loadings$commissionPeriod,
                params = self$Parameters, values = self$Values);
            self$Parameters$Loadings$commissionPeriod =
                min(
                    self$Parameters$Loadings$commissionPeriod,
                    self$Parameters$ContractData$policyPeriod
                )


            # Evaluate deferral period, i.e. if a function is used, calculate its numeric value from the other parameters
            self$Parameters$ContractData$deferralPeriod = valueOrFunction(
                self$Parameters$ContractData$deferralPeriod,
                params = self$Parameters, values = self$Values);
            self$Parameters$ContractData$deferralPeriod =
                min(
                    self$Parameters$ContractData$deferralPeriod,
                    self$Parameters$ContractData$policyPeriod
                )

            # Premium refund period (if applicable):
            self$Parameters$ContractData$premiumRefundPeriod = valueOrFunction(
                self$Parameters$ContractData$premiumRefundPeriod,
                params = self$Parameters, values = self$Values);

            #### #
            # AGES for multiple joint lives:
            #### #
            # For joint lives, some parameters can be given multiple times: age, sex
            # Collect all given values into one vector!

            # TODO: First person has birthDate handled properly, handle all other persons, too!
            age = unlist(args[names(args) == "age"], use.names = FALSE)
            if (!is.null(age) && length(age) > 1) {
                self$Parameters$ContractData$age = c(self$Parameters$ContractData$age[[1]], tail(age, -1));
                # TODO: Calculate ages for all other persons, too. Or rather, allow multiple birthDate values, too
            }
            sex = unlist(args[names(args) == "sex"], use.names = FALSE)
            if (!is.null(sex) && length(sex) > 1) {
                self$Parameters$ContractData$sex = sex;
            }
            if (is.null(self$Parameters$ContractData$ageDifferences)) {
                self$Parameters$ContractData$ageDifferences = diff(self$Parameters$ContractData$age);
            } else {
                self$Parameters$ContractData$ageDifferences = valueOrFunction(
                    self$Parameters$ContractData$ageDifferences,
                    params = self$Parameters, values = self$Values);
            }


            #### #
            # TECHNICAL AGE
            #### #
            # Calculate the technical age (e.g. female are made younger, contracts on joint lives, etc.)
            if (is.null(self$Parameters$ContractData$technicalAge)) {
                self$Parameters$ContractData$technicalAge = self$Parameters$ContractData$age[1]
            } else {
                self$Parameters$ContractData$technicalAge = valueOrFunction(
                    self$Parameters$ContractData$technicalAge,
                    params = self$Parameters, values = self$Values);
            }

            # Evaluate all possibly variable values (mortalityTable depending on sex, etc.)
            self$Parameters$ActuarialBases$mortalityTable = valueOrFunction(
                self$Parameters$ActuarialBases$mortalityTable,
                params = self$Parameters, values = self$Values)

            #### #
            # COSTS PARAMETERS: can be a function => evaluate it to get the real costs
            # This needs to be last, as the costs can depend on present values
            # => mortality table is needed
            #### #
            self$Parameters$Costs = private$evaluateCosts()


            ##### #
            # SUM INSURED: By default, sum insured is required and the premium
            # is calculated from the SI. However, it is also possible (not yet
            # fully implemented) to prescribe the premium and/or initial capital
            # to determine the SI. If neither is given, a warning should be printed
            # and a default of sumInsured=10000 should be applied
            if (is.null(self$Parameters$ContractData$sumInsured) &&
                is.null(self$Parameters$ContractData$premium) &&
                (is.null(self$Parameters$ContractData$initialCapital) ||
                 self$Parameters$ContractData$initialCapital == 0)) {
                self$Parameters$ContractData$sumInsured = 10000
                warning("InsuranceContract: Neither sumInsured nor premium nor initialCapital given! => Assuming default sum insured 10.000. Tariff: ", self$tarif$name)
            }

            invisible(self)
        },

        evaluateCosts = function() {
            self$tarif$getCostValues(params = self$Parameters)
        },

        determineInternalValues = function(...) {
            self$tarif$getInternalValues(params = self$Parameters, ...);
        },

        determineTransitionProbabilities = function(...) {
            self$tarif$getTransitionProbabilities(params = self$Parameters, values = self$Values, ...);
        },
        determineCashFlowsBasic = function(...) {
            self$tarif$getBasicCashFlows(params = self$Parameters, values = self$Values, ...);
        },
        determineCashFlows = function(...) {
            self$tarif$getCashFlows(params = self$Parameters, values = self$Values, ...);
        },
        determinePremiumSum = function(...) {
            sum(self$Values$cashFlows$premiums_advance + self$Values$cashFlows$premiums_arrears);
        },
        determineCashFlowsCosts = function(...) {
            self$tarif$getCashFlowsCosts(params = self$Parameters, values = self$Values, ...);
        },
        calculatePresentValues = function(...) {
            self$tarif$presentValueCashFlows(params = self$Parameters, values = self$Values, ...);
        },
        calculatePresentValuesCosts = function(...) {
            self$tarif$presentValueCashFlowsCosts(params = self$Parameters, values = self$Values, ...);
        },
        calculateSumInsured = function(...) {
            self$tarif$sumInsuredCalculation(params = self$Parameters, values = self$Values, ...)
        },
        calculatePremiums = function(...) {
            self$tarif$premiumCalculation(params = self$Parameters, values = self$Values, ...);
        },
        calculatePresentValuesBenefits = function(...) {
            self$tarif$presentValueBenefits(params = self$Parameters, values = self$Values, ...);
        },
        calculateAbsCashFlows = function(...) {
            self$tarif$getAbsCashFlows(params = self$Parameters, values = self$Values, ...);
        },
        calculateAbsPresentValues = function(...) {
            self$tarif$getAbsPresentValues(params = self$Parameters, values = self$Values, ...);
        },
        calculateReserves = function(...) {
            self$tarif$reserveCalculation(params = self$Parameters, values = self$Values, ...);
        },
        calculateReservesBalanceSheet = function(...) {
            self$tarif$reserveCalculationBalanceSheet(params = self$Parameters, values = self$Values, ...);
        },
        premiumAnalysis = function(...) {
            self$tarif$premiumDecomposition(params = self$Parameters, values = self$Values, ...);
        },
        premiumCompositionSums = function(...) {
            self$tarif$calculateFutureSums(cf = self$Values$premiumComposition, params = self$Parameters, values = self$Values, ...);
        },
        premiumCompositionPV = function(...) {
            self$tarif$calculatePresentValues(self$Values$premiumComposition, params = self$Parameters, values = self$Values, ...);
        },

        profitParticipation = function(...) {
            scens = self$Parameters$ProfitParticipation$scenarios
            lapply(seq_along(scens), function(x) {
                nm = names(scens)[x]
                scn = NULL
                if (!is.null(self$Values$profitScenarios) && length(self$Values$profitScenarios) >= x) {
                    scn = self$Values$profitScenarios[[x]]
                }
                pp = do.call(private$calculateProfitParticipation, c(list(profitScenario = scn, ...), scens[x]))
                res = do.call(private$calculateReservesAfterProfit, c(list(profitScenario = pp, ...), scens[x]))
                if (nm != "" && !is.null(nm)) {
                    self$Values$profitScenarios[[nm]] = pp
                    self$Values$reservesAfterProfit[[nm]] = res
                } else {
                    self$Values$profitScenarios[[x]] = pp
                    self$Values$reservesAfterProfit[[x]] = res
                }
            })

            # For convenience, return the profit participation table:
            # (self$Values$reservesAfterProfit was also changed, but is not returned!)
            self$Values$profitScenarios
        },

        calculateProfitParticipation = function(...) {
            self$tarif$calculateProfitParticipation(params = self$Parameters, values = self$Values, ...);
        },
        calculateReservesAfterProfit = function(profitScenario, ...) {
            self$tarif$reservesAfterProfit(profitScenario = profitScenario, params = self$Parameters, values = self$Values, ...);
        },


        getBasicDataTimeseries = function(...) {
            self$tarif$getBasicDataTimeseries(params = self$Parameters, values = self$Values, ...);
        },

        dummy.private = NULL
    )
)

