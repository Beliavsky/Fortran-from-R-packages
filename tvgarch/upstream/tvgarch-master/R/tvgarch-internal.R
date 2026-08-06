.onAttach <- function(libname, pkgname)
{
txt <- c("\n",
  paste(sQuote("tvgarch"), "version 2.4.3 \n"),
  "\n",
  paste0("Time Varying GARCH Modelling"),
  "\n",
  paste("CRAN website: https://CRAN.R-project.org/package=tvgarch"),
  paste("Personal webpage: https://sites.google.com/site/susanacamposmartins"),
  "\n")
  
  ##print message:
  if(interactive() || getOption("verbose")){
    packageStartupMessage(paste(strwrap(txt, indent = 2,
      exdent = 4), collapse = "\n"))
  }
} #close .onAttach
