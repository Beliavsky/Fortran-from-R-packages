#' @include InsuranceContract.R
#'
#' @import MortalityTables
#' @import scales
#' @importFrom methods is
NULL

# Internal helper function to calculate all values => will be used by
# showVmGlgExamples, testVmGlgExample and vmGlgExample.generateTest
calcVmGlgExample = function(contract, prf = 10, t = 10, t_prf = 12, ...) {

    if (is(contract, "InsuranceTarif")) {
        tariff = contract
        contract = InsuranceContract$new(tariff, ...)
    }
    has.prf = prf < contract$Parameters$ContractData$premiumPeriod;
    if (!is(contract, "InsuranceContract")) {
        stop("First argument to the functions showVmGlgExamples, testVmGlgExample
             and vmGlgExample.generateTest need to be an InsuranceContract or InsuranceTarif object! ",
             "Given object is of class: ",
             paste(class(contract), collapse = ", "));
    }
    if (has.prf) {
        contract.prf = contract$clone();
        contract.prf$premiumWaiver(t = prf)
    }

    vals = list(
        t = t,
        prf = prf,
        t_prf = t_prf,
        contractClosing = contract$Parameters$ContractData$contractClosing,

        TarifName = contract$tarif$name,
        TarifDesc = contract$tarif$desc,
        TarifType = as.character(contract$tarif$tariffType),
        VS = contract$Parameters$ContractData$sumInsured,
        Age = contract$Parameters$ContractData$age,
        policyPeriod = contract$Parameters$ContractData$policyPeriod,
        premiumPeriod = contract$Parameters$ContractData$premiumPeriod,
        deferralPeriod = contract$Parameters$ContractData$deferralPeriod,
        guaranteedPeriod = contract$Parameters$ContractData$guaranteedPeriod,
        Interest = contract$Parameters$ActuarialBases$i,
        MortalityTable = contract$Parameters$ActuarialBases$mortalityTable@name,

        net = contract$Values$premiums[["net"]],
        Zillmer = contract$Values$premiums[["Zillmer"]],
        gross = contract$Values$premiums[["gross"]],
        written = contract$Values$premiums[["written"]],
        savings = contract$Values$premiumComposition[[t + 1, "Zillmer.savings"]],
        risk = contract$Values$premiumComposition[[t + 1, "Zillmer.risk"]],
        ZillmerRes = contract$Values$reserves[[t + 1, "Zillmer"]],
        ZillmerRes.prf = if(has.prf) contract.prf$Values$reserves[[t_prf + 1, "Zillmer"]],
        VwKostenRes = contract$Values$reserves[[t + 1, "gamma"]],
        VwKostenRes.prf = if(has.prf) contract.prf$Values$reserves[[t_prf + 1, "gamma"]],
        Bilanzstichtag = contract$Values$reservesBalanceSheet[[t + 1, "time"]],
        Bilanzreserve = contract$Values$reservesBalanceSheet[[t + 1, "Balance Sheet Reserve"]],
        Praemienuebertrag = contract$Values$reservesBalanceSheet[[t + 1, "unearned Premiums"]],
        Rueckkaufsreserve = contract$Values$reserves[[t + 1, "reduction"]],
        Rueckkaufswert = contract$Values$reserves[[t + 1, "Surrender"]],
        Abschlusskostenruecktrag = contract$Values$reserves[[t + 1, "alphaRefund"]],
        Rueckkaufswert.prf = if(has.prf) contract.prf$Values$reserves[[t_prf + 1, "Surrender"]],
        VS.after_prf = if(has.prf) contract.prf$Values$reserves[[t_prf + 1, "PremiumFreeSumInsured"]],
        VS.prf = if(has.prf) contract$Values$reserves[[t + 1, "PremiumFreeSumInsured"]]

    );
    # manually remove all NULL elements (premium waiver not possible  )
    vals = Filter(Negate(is.null), vals)
    vals
}


#' Display insurance contract calculation example
#'
#' Display the values of the example calculation of the given insurance contract
#' as required by the Austrian regulation (LV-VMGLV, "LV
#' Versicherungsmathematische Grundlagen Verordnung").
#'
#' @param contract The insurance contract to calculate and show
#' @param t Time for which to show all values (except premium-free values)
#' @param prf Time of premium waiver (premium-free)
#' @param t_prf Time for which to show all values after the premium waiver
#' @param file If given, outputs all information to the file rather than the console
#' @param ... Further parameters for generating the contract for a tariff object
#' @returns Name of the output file (invisible) or empty if printed to console
#'
#' @examples
#' library(MortalityTables)
#' mortalityTables.load("Austria_Annuities_AVOe2005R")
#' # A trivial deferred annuity tariff with no costs:
#' tariff = InsuranceTarif$new(name="Test Annuity", type="annuity",
#'     mortalityTable = AVOe2005R.unisex, i=0.01)
#' contract = InsuranceContract$new(
#'     tariff,
#'     age = 35, YOB = 1981,
#'     policyPeriod = 30, premiumPeriod = 15, deferralPeriod = 15,
#'     sumInsured = 1000,
#'     contractClosing = as.Date("2016-10-01")
#' );
#' showVmGlgExamples(contract)
#'
#' # Optionally output to a file rather than the console:
#' \donttest{
#' showVmGlgExamples(contract, file =  tempfile("annuity-example", fileext = ".txt"))
#' }
#' @export
showVmGlgExamples = function(contract, prf = 10, t = 10, t_prf = 12, file = "", ...) {
    if (getOption('LIC.debug.showVmGlgExamples', FALSE)) {
        browser();
    }
    vals = calcVmGlgExample(contract, prf = prf, t = t, t_prf = t_prf, ...);

    has.prf = prf < contract$Parameters$ContractData$premiumPeriod;

    output = c(
        sprintf("Tarif: %s", vals$TarifName),
        vals$tarifDesc,
        sprintf("Typ: %s", vals$tariffType),
        sprintf("VS: %.2f", vals$VS),
        sprintf("Alter: %d, LZ: %d, Pr\u00e4mienzahlung: %d, Aufschubzeit: %d, Garantiezeit: %d",
                vals$age,
                vals$policyPeriod,
                vals$premiumPeriod,
                vals$deferralPeriod,
                vals$guaranteedPeriod),

        sprintf("Rechenzins: %s, Sterbetafel: %s",
                scales::percent(vals$Interest), vals$MortalityTable),
        "",

        "Pr\u00e4mien:",
        "========",
        sprintf("Nettopr\u00e4mie:         %8.2f", vals$net),
        sprintf("Zillmerpr\u00e4mie:       %8.2f", vals$Zillmer),
        sprintf("Bruttopr\u00e4mie:        %8.2f", vals$gross),
        sprintf("Vorgeschr.Pr.:       %8.2f", vals$written),
        "",

        sprintf("Sparpr\u00e4mie (t=%d):   %8.2f", vals$t, vals$savings),
        sprintf("Risikopr\u00e4mie (t=%d): %8.2f", vals$t, vals$risk),
        "",

        "Reserven:",
        "=========",
        sprintf("Pr\u00e4mienpflichtig (t=%d):             %8.2f", vals$t, vals$ZillmerRes),
        sprintf("Pr.Frei (mangels Zahlung) (tprf=%d): %8.2f", vals$t_prf, vals$ZillmerRes.prf),
        sprintf("VwKostenreserve (t=%d):              %8.2f", vals$t, vals$VwKostenRes),
        sprintf("VwKostenreserve (Pr.Frei) (tprf=%d): %8.2f", vals$t_prf, vals$VwKostenRes.prf),
        "",

        sprintf("Bilanzreserve (t=%.2f):             %8.2f",
                vals$Bilanzstichtag, vals$Bilanzreserve),
        sprintf("Pr\u00e4mien\u00fcbertrag (BM=%2d):             %8.2f",
                lubridate::month(vals$contractClosing), vals$Praemienuebertrag),
        "",

        "R\u00fcckkauf und Pr\u00e4mienfreistellung:",
        "=================================",
        sprintf("R\u00fcckkaufsreserve (t=%d):          %8.2f", vals$t, vals$Rueckkaufsreserve),
        sprintf("R\u00fcckkaufswert (t=%d):             %8.2f", vals$t, vals$Rueckkaufswert),
        sprintf("Abschlusskostenr\u00fccktrag (t=%d):   %8.2f", vals$t, vals$Abschlusskostenruecktrag),
        "",

        sprintf("R\u00fcckkaufswert (Prf.) (t=%d):      %8.2f (VS: %.2f)",
                vals$t_prf, vals$Rueckkaufswert, vals$VS.after_prf),
        "",

        sprintf("Pr\u00e4mienfreie VS (t=%d):           %8.2f",
                vals$t, vals$VS.prf)
    );
    output.str = paste(output, collapse = '\r\n')
    cat(output.str, file = file)
    invisible(file)
}


#' Perform unit tests of given standard values of the insurance contract example
#'
#' Check the values of the example calculation of the given insurance contract
#' as required by the Austrian regulation (LV-VMGLV, "LV
#' Versicherungsmathematische Grundlagen Verordnung").
#' Missing params not passed to the function call will be silently ignored and
#' not cause unit test failures.
#'
#' The easiest way to write unit-tests is using the function \code{vmGlgExample.generateTest}
#'
#' @param contract The insurance contract to calculate and check
#' @param t Time for which to check all values (except premium-free values)
#' @param prf Time of premium waiver (premium-free)
#' @param t_prf Time for which to check all values after the premium waiver
#' @param net,Zillmer,gross,written,savings,risk,ZillmerRes,ZillmerRes.prf,VwKostenRes,VwKostenRes.prf,Bilanzreserve,Praemienuebertrag,Rueckkaufsreserve,Rueckkaufswert,Abschlusskostenruecktrag,Rueckkaufswert.prf,VS.prf Values as printed out by showVmGlgExamples
#' @param absTolerance If non-NULL, will ignore small floating point differences. It uses same algorithm as all.equal()
#' @param ... Further parameters for generating the contract for a tariff object
#' @returns None
#'
#' @examples
#' \dontrun{
#' library(MortalityTables)
#' library(testthat)
#' mortalityTables.load("Austria_Annuities_AVOe2005R")
#'
#' test_that("Testtarif", {
#'     # A trivial deferred annuity tariff with no costs:
#'     library(MortalityTables)
#'     mortalityTables.load("Austria_Annuities")
#'     tariff = InsuranceTarif$new(name="Test Annuity", type="annuity",
#'         mortalityTable = AVOe2005R.unisex, i = 0.01)
#'     contract = InsuranceContract$new(
#'         tariff,
#'         age = 35, YOB = 1981,
#'         policyPeriod = 30, premiumPeriod = 15, deferralPeriod = 15,
#'         sumInsured = 1000,
#'         contractClosing = as.Date("2016-10-01")
#'     );
#'
#'     testVmGlgExample(
#'         contract, t = 10,
#'         net = 850.09, # NOT_CHECKED: Zillmer = 950.09,
#'         gross = 850.09,
#'         written = 884.09,
#'         savings = 857.09, risk = -7.00,
#'         ZillmerRes = 9011.40,
#'         ZillmerRes.prf = 9205.96,
#'         VwKostenRes = 0.00,
#'         VwKostenRes.prf = 0.00,
#'         Bilanzreserve = 9250.35,
#'         Praemienuebertrag = 637.567,
#'         Rueckkaufsreserve = 9011.40,
#'         Rueckkaufswert = 9011.40,
#'         Abschlusskostenruecktrag = 0.00,
#'         Rueckkaufswert.prf = 9205.96,
#'         VS.prf = 684.117
#'     )
#' })
#' }
#'
#' @export
testVmGlgExample = function(contract, prf = 10, t = 10, t_prf = 12, net, Zillmer, gross, written, savings, risk,
                            ZillmerRes, ZillmerRes.prf, VwKostenRes, VwKostenRes.prf,
                            Bilanzreserve, Praemienuebertrag,
                            Rueckkaufsreserve, Rueckkaufswert, Abschlusskostenruecktrag,
                            Rueckkaufswert.prf, VS.prf, absTolerance = 0.015,
                            ...
) {
    if (getOption('LIC.debug.testVmGlgExample', FALSE)) {
        browser();
    }
    vals = calcVmGlgExample(contract, prf = prf, t = t, t_prf = t_prf, ...);
    has.prf = prf < contract$Parameters$ContractData$premiumPeriod;

    if (!missing(net)) {
        eval(bquote(expect_equal(vals$net, .(net), tolerance = abs(absTolerance / .(net)))))
    }
    if (!missing(Zillmer)) {
        eval(bquote(expect_equal(vals$Zillmer, .(Zillmer), tolerance = abs(absTolerance / .(Zillmer)))))
    }
    if (!missing(gross)) {
        eval(bquote(expect_equal(vals$gross, .(gross), tolerance = abs(absTolerance / .(gross)))))
    }
    if (!missing(written)) {
        eval(bquote(expect_equal(vals$written, .(written), tolerance = abs(absTolerance / .(written)))))
    }
    if (!missing(savings)) {
        eval(bquote(expect_equal(vals$savings, .(savings), tolerance = abs(absTolerance / .(savings)))))
    }
    if (!missing(risk)) {
        eval(bquote(expect_equal(vals$risk, .(risk), tolerance = abs(absTolerance / .(risk)))))
    }

    if (!missing(ZillmerRes)) {
        eval(bquote(expect_equal(vals$ZillmerRes, .(ZillmerRes), tolerance = abs(absTolerance / .(ZillmerRes)))))
    }
    if (!missing(ZillmerRes.prf)) {
        eval(bquote(expect_equal(vals$ZillmerRes.prf, .(ZillmerRes.prf), tolerance = abs(absTolerance / .(ZillmerRes.prf)))))
    }
    if (!missing(VwKostenRes)) {
        eval(bquote(expect_equal(vals$VwKostenRes, .(VwKostenRes), tolerance = abs(absTolerance / .(VwKostenRes)))))
    }
    if (!missing(VwKostenRes.prf)) {
        eval(bquote(expect_equal(vals$VwKostenRes.prf, .(VwKostenRes.prf), tolerance = abs(absTolerance / .(VwKostenRes.prf)))))
    }


    if (!missing(Bilanzreserve)) {
        eval(bquote(expect_equal(vals$Bilanzreserve, .(Bilanzreserve), tolerance = abs(absTolerance / .(Bilanzreserve)))))
    }
    if (!missing(Praemienuebertrag)) {
        eval(bquote(expect_equal(vals$Praemienuebertrag, .(Praemienuebertrag), tolerance = abs(absTolerance / .(Praemienuebertrag)))))
    }

    if (!missing(Rueckkaufsreserve)) {
        eval(bquote(expect_equal(vals$Rueckkaufsreserve, .(Rueckkaufsreserve), tolerance = abs(absTolerance / .(Rueckkaufsreserve)))))
    }
    if (!missing(Rueckkaufswert)) {
        eval(bquote(expect_equal(vals$Rueckkaufswert, .(Rueckkaufswert), tolerance = abs(absTolerance / .(Rueckkaufswert)))))
    }
    if (!missing(Abschlusskostenruecktrag)) {
        eval(bquote(expect_equal(vals$Abschlusskostenruecktrag, .(Abschlusskostenruecktrag), tolerance = abs(absTolerance / .(Abschlusskostenruecktrag)))))
    }
    if (!missing(Rueckkaufswert.prf)) {
        eval(bquote(expect_equal(vals$Rueckkaufswert.prf, .(Rueckkaufswert.prf), tolerance = abs(absTolerance / .(Rueckkaufswert.prf)))))
    }
    if (!missing(VS.prf)) {
        eval(bquote(expect_equal(vals$VS.prf, .(VS.prf), tolerance = abs(absTolerance / .(VS.prf)))))
        # OR: contract$Values$reserves[t + 1, "PremiumFreeSumInsured"]
    }
}




#' Generate testthat output for unit-testing a tarif implementation
#'
#' This function calculates the required reference values for the given
#' insurance contract as required by the Austrian regulation (LV-VMGLV, "LV
#' Versicherungsmathematische Grundlagen Verordnung") and generates the
#' code for unit-testing the contract with these values. The code printed
#' can be directly copied into a unit test file.
#'
#' @param contract The insurance contract to calculate and generate unit-testing code.
#'                 If an InsuranceTarif object is given, a new contract with default
#'                 values is generated.
#' @param t Time for which to calculate all values (except premium-free values)
#' @param prf Time of premium waiver (premium-free)
#' @param t_prf Time for which to calculated all values after the premium waiver
#' @param ... Further parameters for generating the contract for a tariff object
#' @returns None (code is printed to the console, but no value returned)
#'
#' @examples
#' library(MortalityTables)
#' mortalityTables.load("Austria_Annuities_AVOe2005R")
#' # A trivial deferred annuity tariff with no costs:
#' tariff = InsuranceTarif$new(name="Test Annuity", type="annuity",
#'     mortalityTable = AVOe2005R.unisex, i=0.01)
#' vmGlgExample.generateTest(tariff,
#'     age = 35, YOB = 1981,
#'     policyPeriod = 30, premiumPeriod = 15, deferralPeriod = 15,
#'     sumInsured = 1000,
#'     contractClosing = as.Date("2016-10-01")
#' )
#'
#' @export
vmGlgExample.generateTest = function(contract, prf = 10, t = 10, t_prf = 12, ...) {
    if (getOption('LIC.debug.vmGlgExample.generateTest', FALSE)) {
        browser();
    }
    cntr = deparse(substitute(contract));

    vals = calcVmGlgExample(contract, prf = prf, t = t, t_prf = t_prf, ...);


    code = paste0("test_that(\"", vals$TarifName, "\", {\n");
    code = paste0(code, "\tcontract = InsuranceContract$new(\n\t\t", cntr, ",\n\t\t");
    arguments = sapply(substitute(list(...))[-1], deparse);
    code = paste0(code,
                  paste(names(arguments), arguments, sep = " = ", collapse = ",\n\t\t")
           );
    code = paste0(code, "\n\t);\n")
    code = paste0(code, "\t# showVmGlgExamples(contract, t = ", t, ", prf = ", prf, ", t_prf = ", t_prf, ");\n")
    code = paste0(code, "\t# exportInsuranceContract.xlsx(contract, filename = here(\"", cntr, ".xlsx\"));\n")
    code = paste0(code, "\t# openxlsx::openXL(here(\"", cntr, ".xlsx\"));\n\n")
    code = paste0(code, "\ttestVmGlgExample(\n\t\tcontract, \n\t\tt = ", t, ", prf = ", prf, ", t_prf = ", t_prf, ",\n")

    check.keys = c("net", "Zillmer", "gross", "written", "savings", "risk",
                   "ZillmerRes", "ZillmerRes.prf", "VwKostenRes", "VwKostenRes.prf",
                   "Bilanzreserve", "Praemienuebertrag",
                   "Rueckkaufsreserve", "Rueckkaufswert", "Abschlusskostenruecktrag",
                   "Rueckkaufswert.prf", "VS.prf");

    # Subsetting a list creates NULL entries for missing keys => filter them out
    # E.g. single-premium contracts do not have any premium-free values, so they are NULL.
    cmpvals = Filter(Negate(is.null), vals[check.keys]);
    check.str =  paste(
        names(cmpvals),
        sprintf("%.2f", cmpvals),
        sep = " = ", collapse = ", \n\t\t");
    code = paste0(code, "\t\t", check.str, ",\n\tabsTolerance = 0.01\n\t);\n");
    code = paste0(code, "})\n");
    cat(code);
}



