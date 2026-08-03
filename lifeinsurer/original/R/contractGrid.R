#' @title Create a grid of InsuranceContract objects or premiums with each axis representing one varying parameter
#'
#' @description The function \code{contractGrid} creates a (two- or multi-dimensional) grid
#' of InsuranceContract objects, where each axis represents one of the insurance
#' parameters varying as given in the \code{axes} argument (as a named list).
#'
#' @description The function \code{contractGridPremium} returns a grid of premiums as requested in
#' the \code{premium} parameter rather than the full InsuranceContract objects.
#' It is a convenience wrapper around \code{contractGrid} and is recommended if
#'  one is only interested in a grid of one particular value (typically some
#' kind of premium).
#' The function \code{contractGridPremium} can also be used on an existing
#' \code{contractGrid}-generated grid of contracts to extract grid of numerical
#' values of the specified premiums. If no contract grid is passed to
#' \code{contractGridPremium}, \code{contractGrid} will be called to create it.
#'
#'
#' @details The function \code{contractGrid} will return the full [InsuranceContract]
#' objects, so apply can later be used to extract premiums, reserves and other
#' values to display in a grid. For this feature, one can also use the convenience
#' function \code{contractGridPremium}.
#'
#' The \code{axes} list describing the parameters changing along the axes of the
#' resulting grid is internally expanded with [expand.grid()]. The resulting flat
#' list of parameter (together with the fixed parameters passed as \code{...})
#' is then passed to the \ifelse{html}{\href{../../LifeInsureR/html/InsuranceContract.html#method-new}{\code{InsuranceContract$new()}}}{\code{InsuranceContract$new()()}} call to create the corresponding
#' contract object.
#'
#' To create the human-readable row-/columnnames of the resulting array,
#' the function [makeContractGridDimname()] for each value of the axes, allowing
#' human-readable representations e.g. of a tariff or a mortality table as
#' the dimension label.
#'
#'
#'
#' @param axes List of paramters spanning the dimensions of the grid.
#' @param YOB optional year of bith. If missing, the \code{observationYear} and the contract's age
#' @param observationYear The observation year, for which the grid shall be calculated. If given, the YOB is calculated from it, otherwise the contract's YOB is used
#' @param ... In \code{contractGrid}: Additional parameters to be passed to \ifelse{html}{\href{../../LifeInsureR/html/InsuranceContract.html#method-new}{\code{InsuranceContract$new()}}}{\code{InsuranceContract$new()()}}; In \code{contractGridPremium}: Additional parameters to be passed to \code{contractGrid}.
#
# Params of the contractGridPreimium function:
#' @param contractGrid (optional) existing contract grid from which to derive
#' premiums. If not given, [contractGrid()] is called with all parameters, so
#' \code{...} should contain an \code{axes} argument in that case.
#' @param premium The type of premium to derive (key of the \code{contract$Values$premiums} list.
#' @param .fun The function to extract the desired premium from a contract
#'         object. By default it accesses the premium vector and extracts the
#'         type of premium given in the \code{premium} parameter. One can,
#'         however pass any other extractor function to access e.g. reserves,
#'         cash flows etc. at any desired time.
#' @returns multi-dimensional array of insurance contracts (dimensions defined by the axes argument)
#'
#' @rdname contractGrid
#'
#' @examples
#' library("MortalityTables")
#' mortalityTables.load("Austria_Annuities_AVOe2005R")
#' # A trivial deferred annuity tariff with no costs, premiums during whole
#' # deferral period, 30 years annuity payments:
#' tariff = InsuranceTarif$new(name="Test Annuity", type="annuity", tarif = "Annuity 1A",
#'     mortalityTable = AVOe2005R.unisex, i=0.01,
#'     deferralPeriod = function(params, ...) { params$ContractData$premiumPeriod },
#'     policyPeriod = function(params, ...) { params$ContractData$premiumPeriod + 30 }
#' )
#' contractGrid(
#'     axes = list(
#'         age = seq(20, 60, 10),
#'         premiumPeriod = seq(5,30, 5)
#'     ),
#'     tarif = tariff,
#'     sumInsured = 1000,
#'     contractClosing = as.Date("2023-11-01")
#' )
#'
#' @export
contractGrid = function(axes = list(age = seq(20, 60, 10), policyPeriod = seq(5, 35, 5)), YOB = NULL, observationYear = NULL, ...) {

    obsYear = observationYear;
    # Create all combinations of the variables given for the axes:
    gridByRow = expand.grid(axes, KEEP.OUT.ATTRS = FALSE);
    # Apply InsuranceContract$new to each combination (and add the additional arguments)
    vals = apply(gridByRow, 1, function(axisVals) {
        args = c(as.list(axisVals), list(...));
        if (!is.null(observationYear)) {
            args$YOB = obsYear - args$age;
        }
        do.call(InsuranceContract$new, args)
    })
    dimnames = makeContractGridDimnames(axes)
    array(vals, dim = sapply(axes, length), dimnames = dimnames)
}

# describeIn makeContractGridDimname Create a dimensional name for an [InsuranceTarif] object
#' @export
makeContractGridDimname.InsuranceTarif = function(value) { value$name }
# describeIn makeContractGridDimname Create a dimensional name for an R6 object (using its \code{name} field)
#' @export
makeContractGridDimname.R6 = function(value) { value$name }
# describeIn makeContractGridDimname Create a dimensional name for an \link[MortalityTables]{mortalityTable} object
#' @export
makeContractGridDimname.mortalityTable = function(value) { value@name }
# describeIn makeContractGridDimname Create a dimensional name for a numeric parameter value
#' @export
makeContractGridDimname.numeric = function(value) { value }
# describeIn makeContractGridDimname Create a dimensional name for a numeric parameter value
#' @export
makeContractGridDimname.double = function(value) { value }
# describeIn makeContractGridDimname Create a dimensional name for an object that can be directly used as a human-readable row-/columnname
#' @export
makeContractGridDimname.default = function(value) { value }
#' Create human-readable labels for the dimensions in a [contractGrid()]
#'
#' The function \code{makeContractGridDimname} generates a short, human-readable
#' dimension label for the entries along the axes of a [contractGrid()].
#' The default is to use the \code{value} unchanged as the row-/columnname, but
#' for some parameter values (like a [InsuranceTarif] or \link[MortalityTables]{mortalityTable})
#' a custom method of this function is needed to create the (short) human-readable
#' representation for the axes in the grid.
#'
#' @param value the value along the axis, for which a name should be generated
#' @returns The name of the entry in the dimnames of [contractGrid()]
#' @describeIn makeContractGridDimname Create a short, human-readable dimensional name for an object (default S3 method)
#' @examples
#' library(MortalityTables)
#' mortalityTables.load("Austria_Census")
#'
#' makeContractGridDimname(mort.AT.census.2011.unisex)
#'
#' makeContractGridDimnames(axes = list(
#'     age = seq(30,60,10),
#'     mortalityTable = c(mort.AT.census.2011.unisex, mort.AT.census.2011.male,
#'                        mort.AT.census.2011.female))
#' )
#' @export
makeContractGridDimname = function(value) { UseMethod("makeContractGridDimname", value) }

#' @description The function \code{makeContractGridDimnames} generate proper
#' dimnames for all entries of the axes of a [contractGrid()] by calling
#' \code{makeContractGridDimname} on each of the axes' values
#' @param axes the axes with all names, for which a name should be generated
#' @returns an array of dimnames derived from the axes definitions of [contractGrid()]
#' @describeIn makeContractGridDimname Generate proper dimnames for all entries of the axes of a [contractGrid()]
#' @export
makeContractGridDimnames = function(axes) {
    lapply(axes, function(axis) { lapply(axis, makeContractGridDimname); } )
}

# Create a grid of insurance premiums with each axes representing one varying parameter
#
# The function \code{contractGridPremium} creates a (two- or multi-dimensional) grid
#  of premiums for insurance contracts, where each axes represents one of the
# insurance parameters varying as given in the axes argument (as named list).
#
# This function will return the full InsuranceContract objects, so apply can
# later be used to extract premiums, reserves and other values to display in
# a grid.
#
# If one is only interested in a grid of one particular value (typically some
# kind of premium), the convenience function \code{contractGridPremium} is
# provided, which returns a grid of only the desired value.
#
#' @param contractGrid (optional) existing contract grid from which to derive
#' premiums. If not given, [contractGrid] is called with all parameters, so
#' \code{...} should contain an \code{axes} argument in that case.
#' @param premium The type of premium to derive (key of the \code{contract$Values$premiums} list.
#' @param .fun The function to extract the desired premium from a contract
#'         object. By default it accesses the premium vector and extracts the
#'         type of premium given in the \code{premium} parameter. One can,
#'         however pass any other extractor function to access e.g. reserves,
#'         cash flows etc. at any desired time.
#' @returns a array of premiums (or other contract-specific value) for the grid defined by the \code{axes} argument to \code{contractGrid}
#'
#' @rdname contractGrid
#'
#' @examples
#' library("MortalityTables")
#' mortalityTables.load("Austria_Annuities_AVOe2005R")
#' # A trivial deferred annuity tariff with no costs, premiums during whole
#' # deferral period, 30 years annuity payments:
#' tariff = InsuranceTarif$new(name="Test Annuity", type="annuity", tarif = "Annuity 1A",
#'     mortalityTable = AVOe2005R.unisex, i=0.01,
#'     deferralPeriod = function(params, ...) { params$ContractData$premiumPeriod },
#'     policyPeriod = function(params, ...) { params$ContractData$premiumPeriod + 30 }
#' )
#' contractGridPremium(
#'     axes = list(
#'         age = seq(20, 60, 10),
#'         premiumPeriod = seq(5,30, 5)
#'     ),
#'     tarif = tariff,
#'     sumInsured = 1000,
#'     contractClosing = as.Date("2023-11-01")
#' )
#'
#'
#' @export
contractGridPremium = function(contractGrid = NULL, premium="written", .fun = function(cntr) { cntr$Values$premiums[[premium]] }, ...) {
    if (missing(contractGrid) || is.null(contractGrid)) {
        contractGrid = contractGrid(...)
    }
    apply(contractGrid, 1:length(dim(contractGrid)), function(c) { .fun(c[[1]])})
}


