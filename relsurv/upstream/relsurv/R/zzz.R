#.First.lib <- function(lib, pkg) library.dynam("runproba", pkg, lib)

# use .onLoad instead of .First.lib for use with NAMESPACE and R(>= 1.7.0)
.onLoad <- function(lib, pkg) {
  # library.dynam <- function (chname, package, lib.loc, verbose = getOption("verbose"),
  #                            file.ext = .Platform$dynlib.ext, ...)
  # {
  #   dll_list <- .dynLibs()
  #   if (missing(chname) || !nzchar(chname))
  #     return(dll_list)
  #   package
  #   lib.loc
  #   r_arch <- .Platform$r_arch
  #   chname1 <- paste0(chname, file.ext)
  #   # browser()
  #   for (pkg in "C:/Users/dame_/Dropbox (MF Uni LJ)/Damjan Manevski/Research/relsurv/relsurv_2.2-5"){
  #        #find.package('relsurv_2.2-5', lib.loc, verbose = verbose)) {
  #     DLLpath <- if (nzchar(r_arch))
  #       # file.path(pkg, "libs", r_arch)
  #       "C:/Users/dame_/Dropbox (MF Uni LJ)/Damjan Manevski/Research/relsurv/relsurv_2.2-5/src"
  #     else file.path(pkg, "libs")
  #     # file <- file.path(DLLpath, chname1)
  #     file <- "C:/Users/dame_/Dropbox (MF Uni LJ)/Damjan Manevski/Research/relsurv/relsurv_2.2-5/src/relsurv.dll"
  #     # browser()
  #     if (file.exists(file))
  #       break
  #     else file <- ""
  #   }
  #   if (file == "")
  #     if (.Platform$OS.type == "windows")
  #       stop(gettextf("DLL %s not found: maybe not installed for this architecture?",
  #                     sQuote(chname)), domain = NA)
  #   else stop(gettextf("shared object %s not found", sQuote(chname1)),
  #             domain = NA)
  #   # browser()
  #   file <- file.path(normalizePath(DLLpath, "/", TRUE), chname1)
  #   ind <- vapply(dll_list, function(x) x[["path"]] == file,
  #                 NA)
  #   if (length(ind) && any(ind)) {
  #     if (verbose)
  #       if (.Platform$OS.type == "windows")
  #         message(gettextf("DLL %s already loaded", sQuote(chname1)),
  #                 domain = NA)
  #     else message(gettextf("shared object '%s' already loaded",
  #                           sQuote(chname1)), domain = NA)
  #     return(invisible(dll_list[[seq_along(dll_list)[ind]]]))
  #   }
  #   if (.Platform$OS.type == "windows") {
  #     PATH <- Sys.getenv("PATH")
  #     Sys.setenv(PATH = paste(gsub("/", "\\\\", DLLpath),
  #                             PATH, sep = ";"))
  #     on.exit(Sys.setenv(PATH = PATH))
  #   }
  #   if (verbose)
  #     message(gettextf("now dyn.load(\"%s\") ...", file),
  #             domain = NA)
  #   dll <- if ("DLLpath" %in% names(list(...)))
  #     dyn.load(file, ...)
  #   else dyn.load(file, DLLpath = DLLpath, ...)
  #   .dynLibs(c(dll_list, list(dll)))
  #   invisible(dll)
  # }

  library.dynam("relsurv", pkg, lib)

}#end of .onLoad
