.onAttach <- function(lib, pkg) {
	ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
        #nick <- "Nickname: \"Stack Smashing Detected\""
        #nick <- "Nickname: \"Morpheus and Euripides\""
        #nick <- "Nickname: \"Dyslexical Scoping\""
        #nick <- "Nickname: \"Mendacious Cosmonaut\""
        #nick <- "Nickname: \"Partial Distinction\""
        #nick <- "Nickname: \"Mephistophelian Transition\""
        nick <- "Nickname: \"Idol Comparison\""
        #nick <- "Nickname: \"Perspicuous Ambivalence\""
        #nick <- "Nickname: \"Cats are not required to wear seatbelts\""
        #nick <- "Nickname: \"There weren't many available options\""
        packageStartupMessage(paste(pkg, ver, "    ",nick))
        msg <- paste("\n     The syntax of deldir() has changed since version",
                     "\n     0.0-10.  Read the help!!!.\n")
        packageStartupMessage(msg)
}
