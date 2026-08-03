#' @include showVmGlgExamples.R exportInsuranceContract_xlsx.R InsuranceContract.R InsuranceTarif.R
#'
#' @import stringr
NULL

#' Export the example calculations of an insurance contract
#'
#' Export the given contract to 'Excel' (full history/timeseries of all cash
#' flows, reserves, premiums, etc.) and to a text file (sample calculation
#' required by the Austrian regulation).
#'
#' Three output files are generated:
#'   - \{DATE\}_\{TARIFF\}_Example.xlsx: Full history/timeseries
#'   - \{DATE\}_\{TARIFF\}_Example_PremiumWaiver_t10.xlsx: Full history/timeseries
#'                 after a premium waiver at the given time \code{prf}
#'   - \{DATE\}_\{TARIFF\}_Examples_VmGlg.txt: Example calculation required for the
#'                 Austrian regulation (LV-VMGLV)
#'
#' @param contract The \code{\link{InsuranceContract}} object to be exported
#' @param prf The time of the premium waiver
#' @param outdir The output directory (the file names are not configurable)
#' @param basename The base output filename (sans .xlsx). If missing, a name of
#'       the form 2020-08-01_TARIFNAME_EXTRANAME_RZ0.01_x35_YoB1977_LZ45_PrZ20_VS100000
#'       is used. If given, the main contract without modification will be
#'       exported to basename.xlsx, while the example with premium waiver will be
#'       exported to basename_PremiumWaiver_t10.xlsx and the text file containing
#'       the examples required by the LV-VMGV is exported to basename_VmGlg.txt.
#' @param extraname If basename is not given, this allows a suffix to distinguish
#'       multiple exports.
#' @param ... Further parameters (passed on to \code{\link{showVmGlgExamples}})
#' @returns None
#'
#' @examples
#' library("MortalityTables")
#' mortalityTables.load("Austria_Annuities_AVOe2005R")
#' # A trivial deferred annuity tariff with no costs:
#' tariff = InsuranceTarif$new(name="Test Annuity", type="annuity", tarif = "Annuity 1A",
#'     mortalityTable = AVOe2005R.unisex, i=0.01)
#' contract = InsuranceContract$new(
#'     tariff,
#'     age = 35, YOB = 1981,
#'     policyPeriod = 30, premiumPeriod = 15, deferralPeriod = 15,
#'     sumInsured = 1000,
#'     contractClosing = as.Date("2016-10-01")
#' );
#' \donttest{exportInsuranceContractExample(contract, prf = 10, outdir = tempdir())}
#'
#' @export
exportInsuranceContractExample = function(contract, prf = 10, outdir = ".", basename=NULL, extraname = NULL, ...) {
    if (!("InsuranceContract" %in% class(contract))) {
        stop("First argument to function showVmGlgExamples need to be an InsuranceContract object! ",
             "Given object is of class: ",
             paste(class(contract), collapse = ", "));
    }
    cleanname = contract$tarif$name;
    cleanname = str_replace_all(cleanname, " ", "_");
    cleanname = str_replace_all(cleanname, "[/:]", "-");
    if (missing(basename)) {
        basename = paste(outdir, "/", Sys.Date(), "_", cleanname, sep = "");
        if (!missing(extraname) && !is.null(extraname)) {
            basename = paste(basename, "_", extraname, sep = "")
        }
        params = contract$Parameters
        basename = paste(basename, "_RZ", sprintf("%.2f", params$ActuarialBases$i), "_x", params$ContractData$age, "_YoB", year(params$ContractData$birthDate), "_LZ", params$ContractData$policyPeriod, "_PrZ", params$ContractData$premiumPeriod, "_VS", params$ContractData$sumInsured, sep = "" )
    }

    filename = paste(basename, ".xlsx", sep = "");
    exportInsuranceContract.xlsx(contract, filename);

    contract.prf = contract$clone()
    contract.prf$premiumWaiver(t = prf)
    filename = paste(basename, "_PremiumWaiver_t", prf, ".xlsx", sep = "");
    exportInsuranceContract.xlsx(contract.prf, filename);

    filename = paste(basename, "_VmGlg.txt", sep = "")
    showVmGlgExamples(contract, prf = prf, ..., file = filename)
}
