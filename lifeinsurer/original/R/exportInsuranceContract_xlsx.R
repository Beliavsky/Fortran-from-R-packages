#' @include HelperFunctions.R InsuranceContract.R InsuranceParameters.R InsuranceTarif.R ProfitParticipation.R
#'
#' @import openxlsx
#' @import MortalityTables
#' @import R6
#' @import tidyr
#' @importFrom rlang .data
NULL



################################################ #
# Helper Functions                            ####
################################################ #

addValuesWorksheet = function(wb, sheet, ...) {
  addWorksheet(wb, sheet, gridLines = FALSE, ...)
  # showGridLines(wb, sheet, showGridLines = FALSE)
  setColWidths(wb, sheet, cols = 1:50, widths = "auto", ignoreMergedCells = TRUE)
}


writeAgeQTable = function(wb, sheet, probs, crow = 1, ccol = 1, styles = list()) {
  writeData(wb, sheet, "Sterblichkeiten", startCol = ccol + 2, startRow = crow);
  addStyle(wb, sheet, style = styles$header, rows = crow, cols = ccol + 2, stack = TRUE);
  mergeCells(wb, sheet, rows = crow, cols = (ccol + 2):(ccol + 4))
  writeDataTable(wb, sheet, setInsuranceValuesLabels(probs),
                 startRow = crow + 1, startCol = ccol, colNames = TRUE, rowNames = TRUE,
                 tableStyle = "tableStyleMedium3", withFilter = FALSE, headerStyle = styles$tableHeader);
  # freezePane(wb, sheet, firstActiveRow = crow + 2, firstActiveCol = ccol + 2)
  addStyle(wb, sheet, style = styles$center, rows = (crow + 2):(crow + 1 + dim(probs)[[1]]), cols = ccol:(ccol + 1), gridExpand = TRUE, stack = TRUE);
  addStyle(wb, sheet, style = styles$qx, rows = (crow + 2):(crow + 1 + dim(probs)[[1]]), cols = (ccol + 2):(ccol + 3), gridExpand = TRUE, stack = TRUE);
  dim(probs)[[2]] + 2;
};

writeTableCaption = function(wb, sheet, caption, rows, cols, style = NULL) {
    r = min(rows);
    c = min(cols);
    writeData(wb, sheet, caption, startCol = c, startRow = r);
    if (!is.null(style)) {
        addStyle(wb, sheet, style = style, rows = r, cols = c, stack = TRUE);
    }
    mergeCells(wb, sheet, rows = rows, cols = cols);
}

writeValuesTable = function(wb, sheet, values, caption = NULL, crow = 1, ccol = 1, rowNames = FALSE, tableStyle = "tableStyleLight17", tableName = NULL, withFilter=FALSE, styles = list(), valueStyle = NULL) {
  nrrow = dim(values)[[1]];
  nrcol = dim(values)[[2]];
  addcol = if (rowNames) 1 else 0;
  ecol = ccol + addcol + nrcol - 1;
  if (!missing(caption)) {
      writeTableCaption(wb, sheet, caption, rows = crow, cols = (ccol + addcol):ecol, style = styles$header)
  }

  writeDataTable(wb, sheet, values, startRow = crow + 1, startCol = ccol, colNames = TRUE,
                 rowNames = rowNames, tableStyle = tableStyle,
                 tableName = tableName, withFilter = withFilter, headerStyle = styles$tableHeader)
  if (!missing(valueStyle)) {
    addStyle(wb, sheet, style = valueStyle, rows = (crow + 2):(crow + nrrow + 1), cols = (ccol + addcol):ecol, gridExpand = TRUE, stack = TRUE);
  }
  # width of table is the return value
  nrcol + addcol
};

writePremiumCoefficients = function(wb, sheet, values, tarif = NULL, type = "benefits", crow = crow, ccol = ccol) {
  writeData(wb, sheet, matrix(c(
    "Nettopr\u00e4mie", "", "Zillmerpr\u00e4mie", "", "Bruttopr\u00e4mie", "",
    "rel. zu VS", "rel. zu Pr\u00e4mie", "rel. zu VS", "rel. zu Pr\u00e4mie", "rel. zu VS", "rel. zu Pr\u00e4mie"), 6, 2
  ), startCol = ccol, startRow = crow, colNames = FALSE, borders = "rows", borderColour = "gray5", borderStyle = "thin");
  mergeCells(wb, sheet, cols = ccol, rows = crow:(crow + 1));
  mergeCells(wb, sheet, cols = ccol, rows = (crow + 2):(crow + 3));
  mergeCells(wb, sheet, cols = ccol, rows = (crow + 4):(crow + 5));
  addStyle(wb, sheet, style = createStyle(halign = "left", valign = "center",
                                        borderColour = "gray5", border = "LeftBottomTop",
                                        borderStyle = "thin"),
           rows = crow:(crow + 5), cols = ccol);
  addStyle(wb, sheet, style = createStyle(halign = "right", valign = "center",
                                        borderColour = "gray5", border = "RightBottomTop",
                                        borderStyle = "thin"),
           rows = crow:(crow + 5), cols = ccol + 1);

  # The first column of the benefits coefficients is for "age", which we want to remove
  mod = function(vals) { if (!is.null(vals)) as.data.frame(t(vals)) else NULL };
  if (type == "costs") {
    mod = function(vals) {
      if (is.null(vals)) return(NULL)
      vals = setInsuranceValuesLabels(vals);
      newvals = vals;
      dimn = dimnames(newvals);
      dim(newvals) = c(1, dim(vals));
      dimnames(newvals) = c(list("Coeff"), dimn);
      costValuesAsDF(newvals)
    };
  }
  coeff = rbind(mod(values[["net"]][["SumInsured"]][[type]]),
                mod(values[["net"]][["Premium"]][[type]]),
                mod(values[["Zillmer"]][["SumInsured"]][[type]]),
                mod(values[["Zillmer"]][["Premium"]][[type]]),
                mod(values[["gross"]][["SumInsured"]][[type]]),
                mod(values[["gross"]][["Premium"]][[type]]));

  writeData(wb, sheet, coeff, startCol = ccol + 2, startRow = crow, colNames = FALSE, borders = "rows", borderColour = "gray5", borderStyle = "thin");
  dim(coeff)[[2]]
}

labelsReplace = function(labels) {
  # TODO: Use recode here!

  # Prämienarten
  labels[labels == "unit.net"] = "Netto";
  labels[labels == "unit.Zillmer"] = "Zillmer";
  labels[labels == "unit.gross"] = "Brutto";
  labels[labels == "written_yearly"] = "Verrechnet";
  labels[labels == "written"] = "netto";
  labels[labels == "unitcost"] = "St\u00fcckkosten";
  labels[labels == "written_beforetax"] = "vor Steuer";

  # Kosten
  labels[labels == "alpha"] = "\u03b1";
  labels[labels == "Zillmer"] = "Zill.";
  labels[labels == "beta"] = "\u03b2";
  labels[labels == "gamma"] = "\u03b3";
  labels[labels == "gamma_nopremiums"] = "\u03b3 prf.";
  labels[labels == "unitcosts"] = "StkK";

  # Kosten-Basen
  labels[labels == "SumInsured"] = "VS";
  labels[labels == "SumPremiums"] = "PS";
  labels[labels == "GrossPremium"] = "BP";
  labels[labels == "NetPremium"] = "NP";
  labels[labels == "Constant"] = "";
  labels[labels == "Reserve"] = "Res.";

  # Cash Flows
  labels[labels == "premiums_advance"] = "Pr\u00e4m. vorsch.";
  labels[labels == "premiums_arrears"] = "Pr\u00e4m. nachsch.";
  labels[labels == "additional_capital"] = "Einschuss";
  labels[labels == "guaranteed_advance"] = "Gar. vorsch.";
  labels[labels == "guaranteed_arrears"] = "Gar. nachsch.";
  labels[labels == "survival_advance"] = "Erl. vorsch.";
  labels[labels == "survival_arrears"] = "Erl. nachsch.";

  # Barwerte
  labels[labels == "premiums"] = "Pr\u00e4m.";
  labels[labels == "guaranteed"] = "Gar.";
  labels[labels == "survival"] = "Erl.";
  labels[labels == "after.death"] = "tot";
  labels[labels == "death_SumInsured"] = "Abl. VS";
  labels[labels == "death_GrossPremium"] = "Abl. BP";
  labels[labels == "death"] = "Abl.";
  labels[labels == "disease_SumInsured"] = "Krkh.";
  labels[labels == "death_Refund_past"] = "PrRG(verg.)";
  labels[labels == "death_Refund_future"] = "PrRG(zuk.)";

  labels[labels == "death_PremiumFree"] = "Abl. prf";
  labels[labels == "benefits"] = "Abl.Lst.";
  labels[labels == "benefitsAndRefund"] = "Abl. + RG";

  labels[labels == "once"] = "einm."
  labels[labels == "PremiumPeriod"] = "PD"
  labels[labels == "PremiumFree"] = "Pr.Fr."
  labels[labels == "PolicyPeriod"] = "LZ"

  # R\u00fcckstellungen
  labels[labels == "adequate"] = "ausr.";
  labels[labels == "contractual"] = "vertragl.";
  labels[labels == "conversion"] = "Umrechn.";
  labels[labels == "alphaRefund"] = "\u03b1-R\u00fccktrag";
  labels[labels == "reduction"] = "Sparpr.f\u00fcr DK";
  labels[labels == "PremiumsPaid"] = "Pr.Summe";
  labels[labels == "Surrender"] = "R\u00fcckkauf";
  labels[labels == "PremiumFreeSumInsured"] = "Prf.VS";
  labels[labels == "Balance Sheet Reserve"] = "Bilanzreserve"

  # Pr\u00e4mienzerlegung
  labels[labels == "charged"] = "verrechnet"
  labels[labels == "tax"] = "VSt."
  labels[labels == "loading.frequency"] = "UJZ"
  labels[labels == "rebate.premium"] = "Pr\u00e4m.Rab."
  labels[labels == "rebate.partner"] = "Partn.Rab."
  labels[labels == "unitcosts"] = "StkK"
  labels[labels == "profit.advance"] = "Vw.GB"
  labels[labels == "rebate.sum"] = "Summenrab."
  labels[labels == "charge.noMedicalExam"] = "o.\u00e4rztl.U."
  labels[labels == "gross"] = "Brutto"
  labels[labels == "alpha.noZillmer"] = "\u03b1 (ungez.)";
  labels[labels == "alpha.Zillmer"] = "\u03b1 (gezill.)";
  labels[labels == "net"] = "Netto";
  labels[labels == "risk"] = "Risikopr.";
  labels[labels == "savings"] = "Sparpr.";
  labels[labels == "Zillmer.risk"] = "gez.Risikopr.";
  labels[labels == "Zillmer.savings"] = "gez.Sparpr.";
  labels[labels == "Zillmer.amortization"] = "gez.AK-Tilgung";
  labels[labels == "Zillmer.savings.real"] = "Sparpr.f\u00fcr DK";

  # Vertragseigenschaften
  labels[labels == "InterestRate"] = "i";
  labels[labels == "PolicyDuration"] = "LZ";
  labels[labels == "PremiumPayment"] = "Pr\u00e4mienzhlg.";
  labels[labels == "Premiums"] = "Pr\u00e4mien";
  labels[labels == "age"] = "Alter";
  labels[labels == "Sum insured"] = "Vers.summe";
  labels[labels == "Mortality table"] = "Sterbetafel";
  labels[labels == "i"] = "Garantiezins";
  labels[labels == "Age"] = "Alter";
  labels[labels == "Policy duration"] = "Laufzeit";
  labels[labels == "Premium period"] = "Pr\u00e4mienzahlung";
  labels[labels == "Deferral period"] = "Aufschub";
  labels[labels == "Guaranteed payments"] = "Garantiezeit";

  labels[labels == "time"] = "ZP t";
  labels[labels == "Comment"] = "Bemerkung";
  labels[labels == "type"] = "Art";


  labels
}

setInsuranceValuesLabels = function(vals) {
  dimnames(vals) = lapply(dimnames(vals), labelsReplace);
  vals
}

tableName = function(...) {
  gsub('[^A-Za-z0-9_]', '', paste0(...))
}


exportBlockID = function(wb, sheet, id, cols, rows, styles = c()) {
  writeData(wb, sheet, x = id, xy = c(cols[1], rows[1]))
  addStyle(wb, sheet, style = styles$blockID, rows = rows, cols = cols, stack = TRUE);
}

getContractBlockValues = function(contract) {
  values = data.frame(
      "ID"                  = contract$Parameters$ContractData$id,
      "Tariff"              = contract$tarif$tarif,
      "Start of Contract"   = contract$Parameters$ContractData$contractClosing,
      "Sum insured"         = contract$Parameters$ContractData$sumInsured,
      "Mortality table"     = contract$Parameters$ActuarialBases$mortalityTable@name,
      i                     = contract$Parameters$ActuarialBases$i,
      "Birth Date"          = contract$Parameters$ContractData$birthDate,
      "Age"                 = contract$Parameters$ContractData$age,
      "Technical Age"       = contract$Parameters$ContractData$technicalAge,
      "Policy duration"     = contract$Parameters$ContractData$policyPeriod,
      "Premium period"      = contract$Parameters$ContractData$premiumPeriod,
      "Deferral period"     = contract$Parameters$ContractData$deferralPeriod,
      "Guaranteed payments" = contract$Parameters$ContractData$guaranteed,
      stringsAsFactors = FALSE, check.names = FALSE
    )
  # Append all values from sub-blocks, one per line in the data.frame
  for (b in contract$blocks) {
    values = bind_rows(values, getContractBlockValues(b))
  }
  values
}

getContractBlockPremiums = function(contract) {
  values = NULL
  if (!is.null(contract$Values$premiums)){
    values = bind_cols(
      data.frame(
        "ID" = contract$Parameters$ContractData$id,
        stringsAsFactors = FALSE, check.names = FALSE
      ),
      data.frame(t(contract$Values$premiums), stringsAsFactors = FALSE, check.names = FALSE)
    )
  }
  # Append all values from sub-blocks, one per line in the data.frame
  for (b in contract$blocks) {
    values = bind_rows(values, getContractBlockPremiums(b))
  }
  values
}

#' Convert the multi-dimensional costs array to a data.frame for output to a file
#'
#' @description Convert the cost values array to a tx15 matrix
#' @details Not to be called directly, but implicitly by the [InsuranceContract] object.
#' Convert the array containing cost values like cashflows, present
#' values, etc. (objects of dimension tx5x3) to a matrix with dimensions (tx15)
#'
#' @param costValues Cost definition data structure
costValuesAsDF = function(costValues) {
  as.data.frame.table(costValues, responseName = "Value", stringsAsFactors = TRUE) %>%
    mutate(Var4 = recode(.data$Var4, "Erl." = "")) %>%
    arrange(.data$Var4, .data$Var2, .data$Var3, .data$Var1) %>%
    unite("costtype", "Var2", "Var3", "Var4", sep = " ") %>%
    pivot_wider(names_from = .data$costtype, values_from = .data$Value) %>%
    mutate(Var1 = NULL)
}



exportLoadingsTable = function(wb, sheet, contract, crow, ccol, styles = styles, seprows = 3, tariffs.handled = c()) {
  tarifname = contract$tarif$tarif
  if (!(tarifname %in% tariffs.handled)) {
    costtable = costsDisplayTable(contract$Parameters$Costs)
    # costtable = as.data.frame.table(setInsuranceValuesLabels(contract$Parameters$Costs) )
    # colnames(costtable) = c("Kostenart", "Basis", "Periode", "Kostensatz");
    # costtable = costtable[costtable[,"Kostensatz"] != 0.0000,]
    cap = sprintf("Kosten (Tarif %s)", tarifname)
    writeValuesTable(wb, sheet, costtable, crow = crow, ccol = 1, tableName = tableName("Kosten_", tarifname), styles = styles, caption = cap);
    # writeDataTable(wb, sheet, costtable, startCol = 1, startRow = crow + 1, colNames = TRUE, rowNames = FALSE,
    # tableStyle = "tableStyleMedium3", headerStyle = styles$tableHeader);
    addStyle(wb, sheet, style = styles$cost0, rows = (crow + 2):(crow + dim(costtable)[[1]] + 1), cols = 4, stack = TRUE);
    crow = crow + dim(costtable)[[1]] + 3;
    tariffs.handled = c(tariffs.handled, tarifname)
  }

  for (b in contract$blocks) {
    values = exportLoadingsTable(wb = wb, sheet = sheet, contract = b, crow = crow, ccol = ccol, styles = styles, seprows = seprows, tariffs.handled = tariffs.handled)
    crow = values$crow
    tariffs.handled = values$tariffs.handled
  }
  list(crow = crow, tariffs.handled = tariffs.handled)
}




#' @importFrom rlang .data
exportContractDataTable = function(wb, sheet, contract, ccol = 1, crow = 1, styles = c()) {
  contractValues = getContractBlockValues(contract)
  contractPremiums = getContractBlockPremiums(contract)
  # Some types of tables don't need the birth year -> leave it out rather than throwing an error on opening in Excel!
  # if (is.null(values["YOB"])) values["YOB"] = NULL;

  # General Contract data ####
  # TODO!
  crow = 1;
  writeData(wb, sheet, matrix(c(
    "Tarif:", contract$tarif$tarif,
    "Tarifname:", contract$tarif$name,
    "Description:", contract$tarif$desc
  ), 3, 2, byrow = TRUE), startCol = 1, startRow = 1, colNames = FALSE, rowNames = FALSE);
  mergeCells(wb, sheet, cols = 2:10, rows = 1);
  mergeCells(wb, sheet, cols = 2:10, rows = 2);
  mergeCells(wb, sheet, cols = 2:10, rows = 3);
  addStyle(wb, sheet, style = styles$wrap, rows = 3, cols = 2:11, stack = TRUE);
  addStyle(wb, sheet, style = createStyle(valign = "top"), rows = 1:3, cols = 1:11, gridExpand = TRUE, stack = TRUE);

  crow = crow + 4;
  # Values (parameters, premiums, etc.) of all blocks   ####
  additional_capital = contractPremiums$additional_capital
  if (is.null(contract$Values$premiums)) {
      # Contracts with multiple child blocks do not have their own premium structure => set additional capital to 0
      additional_capital = c(0, additional_capital)
  }
  tmp = contractValues %>%
    mutate(`Initial Capital` = additional_capital) %>%
    select(
      Vertragsteil = .data$ID, Beginn = .data$`Start of Contract`, Tarif = .data$Tariff, .data$`Sum insured`,
      .data$`Initial Capital`,
      .data$`Mortality table`, .data$i, .data$Age, .data$`Policy duration`, .data$`Premium period`,
      .data$`Deferral period`, .data$`Guaranteed payments`)
  writeValuesTable(wb, sheet, values = setInsuranceValuesLabels(tmp),
                   caption = "Basisdaten der Vertragsteile", crow = crow, ccol = 1,
                   tableName = "BlocksBasicData", styles = styles)
  # Second column is start date of contract, fourth is sum insured, sixth is guaranteed interest rate
  addStyle(wb, sheet, style = styles$currency0, rows = crow + 1 + (1:NROW(tmp)), cols = 4:5, gridExpand = TRUE, stack = TRUE);
  addStyle(wb, sheet, style = styles$cost0, rows = crow + 1 + (1:NROW(tmp)), cols = 7, gridExpand = TRUE, stack = TRUE);

  crow = crow + NROW(tmp) + 2 + 2 # 2 rows for caption/table header, 2 rows padding

  # Unit Premiums ####
  tmp = contractPremiums %>%
    select(Vertragsteil = .data$ID, .data$unit.net, .data$unit.Zillmer, .data$unit.gross)
  writeValuesTable(wb, sheet, values = setInsuranceValuesLabels(tmp),
                   caption = "Pr\u00e4miens\u00e4tze (auf VS 1)", crow = crow, ccol = 1,
                   tableName = "UnitPremiums", styles = styles, valueStyle = styles$unitpremiums)
  crow = crow + NROW(tmp) + 2 + 2 # 2 rows for caption/table header, 2 rows padding

  # Yearly Premiums ####
  tmp = contractPremiums %>%
    select(Vertragsteil = .data$ID, .data$net, .data$Zillmer, .data$gross, .data$written_yearly)
  writeValuesTable(wb, sheet, values = setInsuranceValuesLabels(tmp),
                   caption = "Jahrespr\u00e4mien", crow = crow, ccol = 1,
                   tableName = "YearlyPremiums", styles = styles, valueStyle = styles$currency0)
  crow = crow + NROW(tmp) + 2 + 2 # 2 rows for caption/table header, 2 rows padding

  # Written Premiums ####
  tmp = contractPremiums %>%
    select(Vertragsteil = .data$ID, .data$written, .data$unitcost, .data$written_beforetax, .data$tax)
  writeValuesTable(wb, sheet, values = setInsuranceValuesLabels(tmp),
                   caption = "Pr\u00e4mien (pro Zahlungsweise)", crow = crow, ccol = 1,
                   tableName = "WrittenPremiums", styles = styles, valueStyle = styles$currency0)
  crow = crow + NROW(tmp) + 2 + 2 # 2 rows for caption/table header, 2 rows padding


  # Cost structure #######
  crow.history = crow

  vals = exportLoadingsTable(wb = wb, sheet = sheet, contract = contract, crow = crow, ccol = 1, styles = styles, seprows = 3)
  crow = vals$crow

  # costtable = as.data.frame.table(setInsuranceValuesLabels(contract$Parameters$Costs) )
  # colnames(costtable) = c("Kostenart", "Basis", "Periode", "Kostensatz");
  # costtable = costtable[costtable[,"Kostensatz"] != 0.0000,]
  # writeValuesTable(wb, sheet, costtable, crow = crow, ccol = 1, tableName = "Kosten", styles = styles, caption = "Kosten");
  # # writeDataTable(wb, sheet, costtable, startCol = 1, startRow = crow + 1, colNames = TRUE, rowNames = FALSE,
  # # tableStyle = "tableStyleMedium3", headerStyle = styles$tableHeader);
  # addStyle(wb, sheet, style = styles$cost0, rows = (crow + 2):(crow + dim(costtable)[[1]] + 1), cols = 4, stack = TRUE);
  # crow = crow + dim(costtable)[[1]] + 3;

  # Contract history
  # time=t, comment=sprintf("Premium waiver at time %d", t), type = "PremiumWaiver"
  histtime = unlist(lapply(contract$history, function(xl) xl$time));
  histcomment = unlist(lapply(contract$history, function(xl) xl$comment));
  histtype = unlist(lapply(contract$history, function(xl) xl$type));
  writeValuesTable(wb, sheet, setInsuranceValuesLabels(data.frame(time = histtime, Comment = histcomment, type = histtype)),
                   crow = crow.history, ccol = 6, tableName = "Vertragshistorie", styles = styles,
                   caption = "Vertragshistorie");
  crow.history = crow.history + dim(histtime)[[1]] + 3;


}

exportBasicDataTable = function(wb, sheet, contract, ccol = 1, crow = 1, styles = c(), seprows = 5, freeze = TRUE) {
  id = contract$Parameters$ContractData$id
  nrrow = contract$Values$int$l

  blockid.row = crow
  crow = crow + 2
  endrow = (crow + 1 + nrrow)
  if (freeze) {
    freezePane(wb, sheet, firstActiveRow = crow + 2, firstActiveCol = ccol + 2)
  }

  qp = contract$Values$transitionProbabilities[1:nrrow,]; # extract the probabilities once, will be needed in
  cl = ccol

  tbl = qp[,"age", drop = FALSE];
  writeDataTable(wb, sheet, setInsuranceValuesLabels(tbl),
                 startRow = crow + 1, startCol = cl, colNames = TRUE, rowNames = TRUE,
                 tableStyle = "tableStyleMedium3", withFilter = FALSE, headerStyle = styles$tableHeader);
  addStyle(wb, sheet, style = styles$center, rows = (crow + 2):endrow, cols = cl:(cl + 1), gridExpand = TRUE, stack = TRUE);
  cl = cl + dim(tbl)[[2]] + 2;

  cl.table = cl - 1;
  cl = cl + writeValuesTable(
    wb, sheet, as.data.frame(setInsuranceValuesLabels(contract$Values$basicData)),
    crow = crow, ccol = cl, tableName = tableName("Grunddaten_", id), styles = styles,
    caption = "Vertragsgrunddaten im Zeitverlauf") + 1;

  # Change InterestRate column to percent format
  # Change premiumPayment column to single-digit column
  # Change period columnts to normal numbers
  cnames = colnames(contract$Values$basicData);

  r = (crow + 2):endrow;
  addStyle(wb, sheet, style = styles$rate, rows = r, cols = grep("^InterestRate$", cnames) + cl.table, gridExpand = TRUE, stack = TRUE);
  addStyle(wb, sheet, style = styles$digit, rows = r,
           cols = grep("^(PremiumPayment|PolicyDuration|PremiumPeriod)$", cnames) + cl.table,
           gridExpand = TRUE, stack = TRUE);
  addStyle(wb, sheet, style = styles$currency0, rows = r, cols = grep("^(SumInsured|Premiums)$", cnames) + cl.table, gridExpand = TRUE, stack = TRUE);

  exportBlockID(wb, sheet, id = id, rows = blockid.row, cols = ccol:cl, styles = styles)
  crow = endrow + seprows

  for (b in contract$blocks) {
    crow = exportBasicDataTable(
      wb = wb, sheet = sheet, contract = b,
      ccol = ccol, crow = crow, styles = styles, seprows = seprows, freeze = FALSE)
  }
  crow
}

exportReserveTable = function(wb, sheet, contract, ccol = 1, crow = 1, styles = c(), seprows = 5, freeze = TRUE) {
  id = contract$Parameters$ContractData$id
  nrrow = contract$Values$int$l

  blockid.row = crow
  crow = crow + 2

  if (freeze) {
    freezePane(wb, sheet, firstActiveRow = crow + 2, firstActiveCol = ccol + 2)
  }
  qp = contract$Values$transitionProbabilities[1:nrrow,]; # extract the probabilities once, will be needed in
  cl = ccol

  cl = cl + writeAgeQTable(wb, sheet, probs = qp, crow = crow, ccol = cl, styles = styles);
  cl = cl + writeValuesTable(
    wb, sheet, as.data.frame(setInsuranceValuesLabels(contract$Values$reserves)),
    crow = crow, ccol = cl, tableName = tableName("Reserves_", id), styles = styles,
    caption = "Reserven", valueStyle = styles$currency0) + 1;

  oldccol = cl
  cl = cl + writeValuesTable(
    wb, sheet, as.data.frame(setInsuranceValuesLabels(contract$Values$reservesBalanceSheet)),
    crow = crow, ccol = cl, tableName = tableName("Bilanzreserve_", id), styles = styles,
    caption = "Bilanzreserve", valueStyle = styles$currency0) + 1;

  endrow = (crow + 1 + nrrow)
  # First column of balance sheet reserves is the date, the second the fractional time within the contract
  addStyle(wb, sheet, style = createStyle(numFmt = "DATE"), cols = oldccol,
           rows = (crow + 2):endrow, gridExpand = TRUE, stack = TRUE);
  addStyle(wb, sheet, style = createStyle(numFmt = "0.0##"), cols = oldccol + 1,
           rows = (crow + 2):endrow, gridExpand = TRUE, stack = TRUE);

  exportBlockID(wb, sheet, id = id, rows = blockid.row, cols = ccol:cl, styles = styles)
  crow = endrow + seprows

  for (b in contract$blocks) {
    crow = exportReserveTable(
      wb = wb, sheet = sheet, contract = b,
      ccol = ccol, crow = crow, styles = styles, seprows = seprows, freeze = FALSE)
  }
  crow
}

exportProfitParticipationTable = function(wb, sheet, contract, ccol = 1, crow = 1, styles = c(), seprows = 5, freeze = TRUE) {
  id = contract$Parameters$ContractData$id
  nrrow = contract$Values$int$l
  blockid.row = crow
  crow = crow + 2

  if (freeze) {
    freezePane(wb, sheet, firstActiveRow = crow + 2, firstActiveCol = ccol + 2)
  }
  qp = contract$Values$transitionProbabilities[1:contract$Values$int$l,]; # extract the probabilities once, will be needed in

  for (s in names(contract$Values$profitScenarios)) {
    cl = ccol
    sc = contract$Values$profitScenarios[[s]]
    writeData(wb = wb, sheet = sheet, x = s, startRow = crow, startCol = ccol)
    addStyle(wb = wb, sheet = sheet, rows = crow, cols = ccol, style = styles$scenarioID, stack = TRUE)
    cl = cl + writeAgeQTable(wb, sheet, probs = qp, crow = crow, ccol = cl, styles = styles);
    ccol.table = cl - 1;
    cl = cl + writeValuesTable(wb, sheet, as.data.frame(setInsuranceValuesLabels(sc)),
                               crow = crow, ccol = cl, # tableName = tableName("ProfitParticipation_", id, "_", s),
                               styles = styles,
                               # caption = s,
                               valueStyle = styles$currency0) + 1;

    cnames = colnames(sc);
    # Make sure "terminalBonusRate" is NOT matched! Need to use a negative lookahead..
    baseCols = grep("^(?!terminal|TBF).*Base$", cnames, perl = TRUE);
    rateCols = grep("^(?!terminal|TBF).*(Interest|Rate)$", cnames, perl = TRUE);
    profitCols = grep(".*Profit$", cnames);
    terminalBonusCols = grep("^terminal.*", cnames);
    TBFCols = grep("^TBF.*", cnames);
    deathCols = grep("^death.*", cnames);
    surrenderCols = grep("^surrender.*", cnames);
    premiumWaiverCols = grep("^premiumWaiver.*", cnames);

    endrow = (crow + 1 + nrrow)

    # Rates are displayed in %:
    addStyle(wb, sheet, style = styles$rate, rows = (crow + 2):endrow, cols = rateCols + ccol.table, gridExpand = TRUE, stack = TRUE);

    # Add table headers for the various sections:
    if (length(baseCols) > 0) {
      writeTableCaption(wb, sheet, "Basisgr\u00f6\u00dfen", rows = crow, cols = baseCols + ccol.table, style = styles$header);
    }
    if (length(rateCols) > 0) {
      writeTableCaption(wb, sheet, "Gewinnbeteiligungss\u00e4tze", rows = crow, cols = rateCols + ccol.table, style = styles$header);
    }
    if (length(profitCols) > 0) {
      writeTableCaption(wb, sheet, "GB Zuweisungen", rows = crow, cols = profitCols + ccol.table, style = styles$header);
    }
    if (length(terminalBonusCols) > 0) {
      writeTableCaption(wb, sheet, "Schlussgewinn", rows = crow, cols = terminalBonusCols + ccol.table, style = styles$header);
    }
    if (length(TBFCols) > 0) {
      writeTableCaption(wb, sheet, "Schlussgewinnfonds", rows = crow, cols = TBFCols + ccol.table, style = styles$header);
    }
    if (length(deathCols) > 0) {
      writeTableCaption(wb, sheet, "Todesfallleistung", rows = crow, cols = deathCols + ccol.table, style = styles$header);
    }
    if (length(surrenderCols) > 0) {

            writeTableCaption(wb, sheet, "R\u00fcckkauf", rows = crow, cols = surrenderCols + ccol.table, style = styles$header);
    }
    if (length(premiumWaiverCols) > 0) {
      writeTableCaption(wb, sheet, "Pr\u00e4mienfreistellung", rows = crow, cols = premiumWaiverCols + ccol.table, style = styles$header);
    }

    exportBlockID(wb, sheet, id = id, rows = blockid.row, cols = ccol:cl, styles = styles)
    crow = endrow + seprows
  }

  for (b in contract$blocks) {
    crow = exportProfitParticipationTable(
      wb = wb, sheet = sheet, contract = b,
      ccol = ccol, crow = crow + seprows, styles = styles, seprows = seprows, freeze = FALSE)
  }
  crow
}

exportPremiumCompositionTable = function(wb, sheet, contract, ccol = 1, crow = 1, styles = c(), seprows = 5, freeze = TRUE) {
  id = contract$Parameters$ContractData$id
  nrrow = contract$Values$int$l

  blockid.row = crow
  crow = crow + 2

  if (freeze) {
    freezePane(wb, sheet, firstActiveRow = crow + 2, firstActiveCol = ccol + 2)
  }
  qp = contract$Values$transitionProbabilities[1:nrrow,]; # extract the probabilities once, will be needed in
  cl = ccol

  cl = cl + writeAgeQTable(wb, sheet, probs = qp, crow = crow, ccol = cl, styles = styles);
  cl = cl + writeValuesTable(
    wb, sheet, as.data.frame(setInsuranceValuesLabels(contract$Values$premiumComposition)),
    crow = crow, ccol = cl, tableName = tableName("Premium_Decomposition_", id), styles = styles,
    caption = "Pr\u00e4mienzerlegung", valueStyle = styles$currency0) + 1;

  crow = crow + nrrow + 4;

  cl = ccol
  cl = cl + writeAgeQTable(wb, sheet, probs = qp, crow = crow, ccol = cl, styles = styles);
  cl = cl + writeValuesTable(
    wb, sheet, as.data.frame(setInsuranceValuesLabels(contract$Values$premiumCompositionSums)),
    crow = crow, ccol = cl, tableName = tableName("Premium_DecompositionSums_", id), styles = styles,
    caption = "Pr\u00e4mienzerlegung (Summe zuk\u00fcnftiger Pr\u00e4mien)", valueStyle = styles$currency0) + 1;

  crow = crow + nrrow + 4;

  cl = ccol
  cl = cl + writeAgeQTable(wb, sheet, probs = qp, crow = crow, ccol = cl, styles = styles);
  cl = cl + writeValuesTable(
    wb, sheet, as.data.frame(setInsuranceValuesLabels(contract$Values$premiumCompositionPV)),
    crow = crow, ccol = cl, tableName = tableName("Premium_DecompositionPV_", id), styles = styles,
    caption = "Pr\u00e4mienzerlegung(Barwerte zuk\u00fcnftiger Pr\u00e4mien)", valueStyle = styles$currency0) + 1;

  endrow = (crow + 1 + nrrow)

  # Insert a separator line (with minimum height and dark background)
  addStyle(wb, sheet, style = styles$separator, rows = (endrow + 2), cols = ccol:cl, gridExpand = TRUE, stack = TRUE)


  exportBlockID(wb, sheet, id = id, rows = blockid.row, cols = ccol:cl, styles = styles)
  crow = endrow + 2 + seprows

  for (b in contract$blocks) {
    crow = exportPremiumCompositionTable(
      wb = wb, sheet = sheet, contract = b,
      ccol = ccol, crow = crow, styles = styles, seprows = seprows, freeze = FALSE)
  }
  crow
}

exportAbsPVTable = function(wb, sheet, contract, ccol = 1, crow = 1, styles = c(), seprows = 5, freeze = TRUE) {
  id = contract$Parameters$ContractData$id
  nrrow = contract$Values$int$l

  blockid.row = crow
  crow = crow + 2
  endrow = (crow + 1 + nrrow)

  if (freeze) {
    freezePane(wb, sheet, firstActiveRow = crow + 2, firstActiveCol = ccol + 2)
  }
  qp = contract$Values$transitionProbabilities[1:nrrow,]; # extract the probabilities once, will be needed in
  cl = ccol

  cl = cl + writeAgeQTable(wb, sheet, probs = qp, crow = crow, ccol = cl, styles = styles);
  cl = cl + writeValuesTable(
    wb, sheet, as.data.frame(setInsuranceValuesLabels(contract$Values$absPresentValues)),
    crow = crow, ccol = cl, tableName = tableName("PVabsolute_", id), styles = styles,
    caption = "abs. Leistungs- und Kostenbarwerte", valueStyle = styles$currency0) + 1;

  exportBlockID(wb, sheet, id = id, rows = blockid.row, cols = ccol:cl, styles = styles)
  crow = endrow + seprows

  for (b in contract$blocks) {
    crow = exportAbsPVTable(
      wb = wb, sheet = sheet, contract = b,
      ccol = ccol, crow = crow, styles = styles, seprows = seprows, freeze = FALSE)
  }
  crow
}

exportAbsCFTable = function(wb, sheet, contract, ccol = 1, crow = 1, styles = c(), seprows = 5, freeze = TRUE) {
  id = contract$Parameters$ContractData$id
  nrrow = contract$Values$int$l

  blockid.row = crow
  crow = crow + 2
  endrow = (crow + 1 + nrrow)

  if (freeze) {
    freezePane(wb, sheet, firstActiveRow = crow + 2, firstActiveCol = ccol + 2)
  }
  qp = contract$Values$transitionProbabilities[1:nrrow,]; # extract the probabilities once, will be needed in
  cl = ccol

  cl = cl + writeAgeQTable(wb, sheet, probs = qp, crow = crow, ccol = cl, styles = styles);
  cl = cl + writeValuesTable(
    wb, sheet, as.data.frame(setInsuranceValuesLabels(contract$Values$absCashFlows)),
    crow = crow, ccol = cl, tableName = tableName("CFabsolute_", id), styles = styles,
    caption = "abs. Leistungs- und Kostencashflows", valueStyle = styles$currency0) + 1;

  exportBlockID(wb, sheet, id = id, rows = blockid.row, cols = ccol:cl, styles = styles)
  crow = endrow + seprows

  for (b in contract$blocks) {
    crow = exportAbsCFTable(
      wb = wb, sheet = sheet, contract = b,
      ccol = ccol, crow = crow, styles = styles, seprows = seprows, freeze = FALSE)
  }
  crow
}

# Return list of all column indices witout any non-zero entries => these columns are to be hidden later
empty.col.indices = function(df) {
  unname(which(colSums(df != 0) == 0))
}

exportPVTable = function(wb, sheet, contract, ccol = 1, crow = 1, styles = c(), seprows = 5, freeze = TRUE) {
  empty.cols = NULL;
  if (!is.null(contract$Values$presentValues)) {
    id = contract$Parameters$ContractData$id
    nrrow = contract$Values$int$l

    blockid.row = crow
    crow = crow + 2
    if (freeze) {
      freezePane(wb, sheet, firstActiveRow = crow + 2 + 6, firstActiveCol = ccol + 2)
    }

    # Time the premium was last calculated (i.e. access the present values at that time rather than 0 in the formulas for the premium)
    tPrem = contract$Values$int$premiumCalculationTime

    qp = contract$Values$transitionProbabilities[1:nrrow,]; # extract the probabilities once, will be needed in
    costPV = costValuesAsDF(setInsuranceValuesLabels(contract$Values$presentValuesCosts))
    cl = ccol

    # We add six lines before the present values to show the coefficients for the premium calculation
    cl = cl + writeAgeQTable(wb, sheet, probs = qp, crow = crow + 6, ccol = cl, styles = styles);

    # Store the start/end columns of the coefficients, since we need them later in the formula for the premiums!
    w1 = writePremiumCoefficients(wb, sheet, contract$Values$premiumCoefficients, type = "benefits", crow = crow, ccol = cl - 2, tarif = contract$tarif);
    area.premiumcoeff = paste0(int2col(cl), "%d:", int2col(cl + w1 - 1), "%d");
    area.premiumvals  = paste0("$", int2col(cl), "$", crow + 6 + 2 + tPrem, ":$", int2col(cl + w1 - 1), "$", crow + 6 + 2 + tPrem);
    empty.cols = c(empty.cols, empty.col.indices(contract$Values$presentValues) - 1 + cl)
    cl = cl + writeValuesTable(wb, sheet, as.data.frame(setInsuranceValuesLabels(contract$Values$presentValues)),
                                  crow = crow + 6, ccol = cl, tableName = tableName("PresentValues_Benefits_", id), styles = styles,
                                  caption = "Leistungsbarwerte", valueStyle = styles$pv0) + 1;

    w2 = writePremiumCoefficients(wb, sheet, contract$Values$premiumCoefficients, type = "costs", crow = crow, ccol = cl - 2, tarif = contract$tarif);
    area.costcoeff = paste0(int2col(cl), "%d:", int2col(cl + w2 - 1), "%d");
    area.costvals  = paste0("$", int2col(cl), "$", crow + 6 + 2 + tPrem, ":$", int2col(cl + w2 - 1), "$", crow + 6 + 2 + tPrem);
    empty.cols = c(empty.cols, empty.col.indices(as.data.frame(costPV)) - 1 + cl)
    cl = cl + writeValuesTable(wb, sheet, as.data.frame(costPV),
                               crow = crow + 6, ccol = cl, tableName = tableName("PresentValues_Costs_", id), styles = styles,
                               caption = "Kostenbarwerte", valueStyle = styles$cost0) + 1;

    # Now print out the formulas for premium calculation into the columns 2 and 3:
    writeData(wb, sheet, as.data.frame(c("Nettopr\u00e4mie", contract$Values$premiums[["net"]],"Zillmerpr\u00e4mie", contract$Values$premiums[["Zillmer"]], "Bruttopr\u00e4mie", contract$Values$premiums[["gross"]])), startCol = ccol, startRow = crow, colNames = FALSE, borders = "rows");
    for (i in 0:5) {
      writeFormula(wb, sheet, paste0("SUMPRODUCT(", sprintf(area.premiumcoeff, crow + i, crow + i), ", ", area.premiumvals, ") + SUMPRODUCT(", sprintf(area.costcoeff, crow + i, crow + i), ", ", area.costvals, ")"), startCol = ccol + 2, startRow = crow + i);
      addStyle(wb, sheet, style = styles$pv0, rows = crow + i, cols = ccol + 2, stack = TRUE);
    }
    for (i in c(0,2,4)) {
      writeFormula(wb, sheet, paste0(int2col(3), crow + i, "/", int2col(ccol + 2), crow + i + 1), startCol = ccol + 1, startRow = crow + i);
      addStyle(wb, sheet, style = styles$pv0, rows = crow + i, cols = ccol + 1, stack = TRUE);
    }
    for (i in c(1,3,5)) {
      writeFormula(wb, sheet, paste0(int2col(2), crow + i - 1, "*", contract$Parameters$ContractData$sumInsured), startCol = ccol + 1, startRow = crow + i);
      addStyle(wb, sheet, style = styles$currency0, rows = crow + i, cols = ccol:(ccol + 1), stack = TRUE, gridExpand = TRUE);
    }


    endrow = (crow + 6 + 1 + nrrow)
    exportBlockID(wb, sheet, id = id, rows = blockid.row, cols = ccol:cl, styles = styles)
    crow = endrow + seprows
  }

  for (b in contract$blocks) {
    res = exportPVTable(
      wb = wb, sheet = sheet, contract = b,
      ccol = ccol, crow = crow, styles = styles, seprows = seprows, freeze = FALSE)
    # exportPVTable returns a list of the next excel table row and a list of all empty columns (to be hidden later on)
    crow = res$crow
    if (is.null(empty.cols)) {
      empty.cols = res$empty.cols
    } else {
      empty.cols = intersect(empty.cols, res$empty.cols)
    }
  }
  list(crow = crow, empty.cols = empty.cols)
}

exportCFTable = function(wb, sheet, contract, ccol = 1, crow = 1, styles = c(), seprows = 5, freeze = TRUE) {
  # Write out only if the contract has unit cash flows (i.e. it is a leaf contract block without children on its own!)
  empty.cols = NULL;
  if (!is.null(contract$Values$cashFlows)) {
    id = contract$Parameters$ContractData$id
    nrrow = contract$Values$int$l
    empty.cols = c()

    blockid.row = crow
    crow = crow + 2
    endrow = (crow + 1 + nrrow)

    if (freeze) {
      freezePane(wb, sheet, firstActiveRow = crow + 2, firstActiveCol = ccol + 2)
    }
    qp = contract$Values$transitionProbabilities[1:nrrow,]; # extract the probabilities once, will be needed in
    cl = ccol

    cl = cl + writeAgeQTable(wb, sheet, probs = qp, crow = crow, ccol = cl, styles = styles);
    empty.cols = c(empty.cols, empty.col.indices(contract$Values$cashFlows)- 1  + cl)
    cl = cl + writeValuesTable(
      wb, sheet, setInsuranceValuesLabels(contract$Values$cashFlows),
      crow = crow, ccol = cl, tableName = tableName("CF_", id), styles = styles,
      caption = "Leistungscashflows", valueStyle = styles$hide0) + 1;
    costCF = costValuesAsDF(setInsuranceValuesLabels(contract$Values$cashFlowsCosts))
    empty.cols = c(empty.cols, empty.col.indices(costCF)- 1  + cl)
    cl = cl + writeValuesTable(
      wb, sheet, costCF,
      crow = crow, ccol = cl, tableName = tableName("CFcosts_", id), styles = styles,
      caption = "Kostencashflows", valueStyle = styles$cost0) + 1;

    exportBlockID(wb, sheet, id = id, rows = blockid.row, cols = ccol:cl, styles = styles)
    crow = endrow + seprows
  }

  for (b in contract$blocks) {
    res = exportCFTable(
      wb = wb, sheet = sheet, contract = b,
      ccol = ccol, crow = crow, styles = styles, seprows = seprows, freeze = FALSE)
    # exportCFTable returns a list of the next excel table row and a list of all empty columns (to be hidden later on)
    crow = res$crow
    if (is.null(empty.cols)) {
      empty.cols = res$empty.cols
    } else {
      empty.cols = intersect(empty.cols, res$empty.cols)
    }
  }
  list(crow = crow, empty.cols = empty.cols)
}





############################################################################### #
#
# The actual export function
#
#    exportInsuranceContract.xlsx(contract, filename)
#
############################################################################### #

#' Export an insurance contract (object of class [InsuranceContract]) to an 'Excel' file
#'
#' @details The function \code{exportInsuranceContract.xlsx} exports an object
#' of class [InsuranceContract] to an 'Excel' file. All basic data, as well as
#' the time series of (absolute and unit) cash flows, reserves, premiums, premium
#' composition and all profit participation scenarios are exported to the file
#' in nicely looking tables.
#'
#' No new calculations are done in this function. It only prints out the values
#' stored in \code{contract$Values}.
#'
#' @param contract The insurance contract to export
#' @param filename Target 'Excel' file name for export
#' @returns None
#'
#' @examples
#' library("MortalityTables")
#' mortalityTables.load("Austria_Annuities_AVOe2005R")
#' # A trivial deferred annuity tariff with no costs:
#' tariff = InsuranceTarif$new(name = "Test Annuity", type = "annuity", tarif = "Annuity 1A",
#'     mortalityTable = AVOe2005R.unisex, i=0.01)
#' contract = InsuranceContract$new(
#'     tariff,
#'     age = 35, YOB = 1981,
#'     policyPeriod = 30, premiumPeriod = 15, deferralPeriod = 15,
#'     sumInsured = 1000,
#'     contractClosing = as.Date("2016-10-01")
#' );
#' \donttest{exportInsuranceContract.xlsx(contract, tempfile("Example_annuity_contract", fileext = ".xlsx"))}
#' @export
exportInsuranceContract.xlsx = function(contract, filename) {
  # TODO: argument checking for contract and filename

  ###
  nrrows = contract$Values$int$l; # Some vectors are longer(e.g. qx), so determine the max nr or rows
  qp = contract$Values$transitionProbabilities[1:nrrows,]; # extract the probabilities once, will be needed in every sheet

  ############################################### #
  # Style information                          ####
  ############################################### #
  styles = list(
    blockID = createStyle(border = "Bottom", borderColour = "#ab6310", fgFill = "#d0d0d0", halign = "left", textDecoration = "bold", fontSize = 14),
    scenarioID = createStyle(halign = "left", textDecoration = "bold", fontSize = 14),
    header = createStyle(border = "TopBottomLeftRight", borderColour = "#DA9694", borderStyle = "medium",
                         fgFill = "#C0504D", fontColour = "#FFFFFF",
                         halign = "center", valign = "center", textDecoration = "bold"),
    tableHeader = createStyle(#border = "To2pLeftRight", borderColour = "#DA9694", borderstyle = "medium",
                              #bgFill = "#2C0504D", fontColour = "#FFFFFF",
                              fgFill = "#E0E0E0",
                              halign = "center", valign = "center", textDecoration = "bold"),
    hide0 = createStyle(numFmt = "General; General; \"\""),
    currency0 = createStyle(numFmt = "[$\u20ac-C07] #,##0.00;[red]-[$\u20ac-C07] #,##0.00;\"\""),
    cost0 = createStyle(numFmt = "0.0##%; 0.0##%; \"\""),
    pv0 = createStyle(numFmt = "0.00000;-0.00000;\"\""),
    qx = createStyle(numFmt = "0.000000"),
    rate = createStyle(numFmt = "0.00####%"),
    digit = createStyle(numFmt = "0;-0;\"\""),
    wrap = createStyle(wrapText = TRUE),
    center = createStyle(halign = "center", valign = "center"),
    separator = createStyle(fgFill = "#000000"),
    unitpremiums = createStyle(numFmt = "0.00000%; -0.00000%;"),
    date = createStyle(numFmt = "DATE")
  );

  ############################################### #
  # General Workbook setup                     ####
  ############################################### #
  wb = openxlsx::createWorkbook();


  ############################################### #
  # Basic parameters                           ####
  ############################################### #

  # Print out general Contract and Tariff information, including results
  sheet = "Tarifinformationen"
  addValuesWorksheet(wb, sheet);
  exportContractDataTable(
    wb, sheet = sheet, contract = contract,
    ccol = 1, crow = 1, styles = styles)


  ################################################# #
  # Print out Basic contract data as time series ####
  ################################################# #

  sheet = "Basisdaten";
  addValuesWorksheet(wb, sheet);
  exportBasicDataTable(
    wb = wb, sheet = sheet, contract = contract,
    ccol = 1, crow = 4, styles = styles)


  ############################################### #
  # Print out Reserves                         ####
  ############################################### #

  sheet = "Reserven";
  addValuesWorksheet(wb, sheet);
  exportReserveTable(
    wb = wb, sheet = sheet, contract = contract,
    ccol = 1, crow = 4, styles = styles)


  ################################################ #
  # Print out Profit Participation              ####
  ################################################ #

  if (!is.null(contract$Values$profitParticipation)) {
    sheet = "Gewinnbeteiligung";
    addValuesWorksheet(wb, sheet);
    exportProfitParticipationTable(
      wb = wb, sheet = sheet, contract = contract,
      ccol = 1, crow = 4, styles = styles)
  }


  ############################################### #
  # Print out premium decomposition            ####
  ############################################### #

  sheet = "Pr\u00e4mienzerlegung";
  addValuesWorksheet(wb, sheet);
  exportPremiumCompositionTable(
    wb = wb, sheet = sheet, contract = contract,
    ccol = 1, crow = 4, styles = styles)


  ################################################ #
  # Print out absolute values of present values ####
  ################################################ #

  sheet = "abs.Barwerte";
  addValuesWorksheet(wb, sheet);
  exportAbsPVTable(
    wb = wb, sheet = sheet, contract = contract,
    ccol = 1, crow = 4,styles = styles)


  ############################################### #
  # Print out absolute values for cash flows   ####
  ############################################### #

  sheet = "abs.Cash-Flows";
  addValuesWorksheet(wb, sheet);
  exportAbsCFTable(
    wb = wb, sheet = sheet, contract = contract,
    ccol = 1, crow = 4, styles = styles)


  ############################################### #
  # Print out present values                   ####
  ############################################### #

  # TODO-blocks
  sheet = "Barwerte";
  addValuesWorksheet(wb, sheet);
  res = exportPVTable(
    wb = wb, sheet = sheet, contract = contract,
    ccol = 1, crow = 4,styles = styles)
  if (length(res$empty.cols) > 0) {
    setColWidths(wb, sheet = sheet, cols = res$empty.cols, hidden = TRUE)
  }


  ############################################## #
  # Print out cash flows                      ####
  ############################################## #

  sheet = "Cash-Flows";
  addValuesWorksheet(wb, sheet);
  res = exportCFTable(
    wb = wb, sheet = sheet, contract = contract,
    ccol = 1, crow = 4, styles = styles)
  if (length(res$empty.cols) > 0) {
    setColWidths(wb, sheet = sheet, cols = res$empty.cols, hidden = TRUE)
  }




  ############################################## #
  # Save file                                 ####
  ############################################## #

  openxlsx::saveWorkbook(wb, filename, overwrite = TRUE)

}
