# Custom project creation function, copied as a starting point from
#     https://rstudio.github.io/rstudio-extensions/rstudio_project_templates.html
# Copying functions and replacements originally taken from
#     https://blog.devgenius.io/make-your-own-rstudio-project-template-1f77c4888e79
# and heavily modified

create_LIR_project <- function(path, ...) {
    dots <- list(...)
    CompanyName = dots$Company



    ########################################
    #####  1. ensure destination path exists, sanity checks on the company name

    dir.create(path, recursive = TRUE, showWarnings = FALSE)


    ########################################
    #####  2. copy over all files from the template directory

    LIC.src <- function (..., lib.loc = NULL, mustWork = FALSE){
        system.file("rstudio", "templates", "project", ...,
                    package = "LifeInsureR",
                    lib.loc = lib.loc, mustWork = mustWork)
    }
    from <- LIC.src("LifeInsureR")
    fs::dir_copy(path = from, new_path = path, overwrite = TRUE)


    ########################################
    #####  3. Rename all files with XXXCOMPANYXXX -> CompanyName

    patternFiles = list.files(path = path, pattern = "XXXCOMPANYXXX", full.names = TRUE, recursive = TRUE)
    replacementFiles = str_replace(patternFiles, "XXXCOMPANYXXX", CompanyName)
    file.rename(patternFiles, replacementFiles)
    copied_files <- list.files(path = path, recursive = TRUE)


    ########################################
    #####  4. Replace all file contents with XXXCOMPANYXXX -> CompanyName

    LIC.replaceFileContents <- function(file, pattern, replace) {
        suppressWarnings(tx <- readLines(file))
        tx2 <- gsub(pattern = pattern, replacement = replace, x = tx)
        writeLines(tx2, con = file)
    }

    LIC.replaceAllFilesContents <- function(
        copied_files, pattern, replace, path) {
        # Going through copied files to replace package name
        for (f in copied_files) {
            copied_file <- file.path(path, f)

            try({
                LIC.replaceFileContents(
                    file = copied_file,
                    pattern = pattern,
                    replace = replace
                )
            }, silent = TRUE)
        }
    }
    LIC.replaceAllFilesContents(
        copied_files,
        pattern = "XXXCOMPANYXXX",
        replace = CompanyName,
        path
    )


    ########################################
    #####  5. All other small tasks


}
