.onAttach <- function(library, pkg)
{
  # require("stats4")  
  # require("methods")
  # require("mnormt")
  # require("numDeriv")
  if(interactive())
  {
    # pkg <- Package("sn")
    meta <- packageDescription("sn")
    overview <- 'help("overview-sn")'
    packageStartupMessage(
      "Package 'sn', ", meta$Version, " (", meta$Date, "). \n",
      "Type 'help(SN)' and '", overview, "' for basic information.\n",
      "The package redefines function 'sd' but its usual working is unchanged.")
  }
  invisible()
}
