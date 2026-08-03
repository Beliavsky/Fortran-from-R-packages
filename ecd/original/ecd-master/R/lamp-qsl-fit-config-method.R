#' Read QSLD fit config
#'
#' Read QSLD fit config for plot or custom fit utility.
#' The xtable print utility is also provided to
#' generate high quality latex output for publication purpose.
#'
#' @param key character, the top-level key for config, default to NULL.
#' @param extdata_dir optionally specify user's own extdata folder, default is NULL.
#' @param filename character, optionally specify user's own config file name, default is NULL.
#' @param df the data frame generated from \code{lamp.qsl_fit_config}.
#'
#' @return The data.frame object for the config
#'
#' @keywords sample data
#'
#' @importFrom xtable xtable
#'
#' @export lamp.qsl_fit_config
#' @export lamp.qsl_fit_config_xtable
#'
#' @examples
#' c <- lamp.qsl_fit_config()
#'
### <======================================================================>
"lamp.qsl_fit_config" <- function(key=NULL, extdata_dir=NULL, filename=NULL)
{
    if (is.null(filename)) filename <- "qsl_fit_conf.yml"
    conf_file <- .ecd.locate_file(filename, extdata_dir=extdata_dir)
    # read the conf file
    if (! file.exists(conf_file)) {
        stop(paste("conf_file does not exist:", conf_file))
    }
    conf_all <- yaml.load_file(conf_file)

    stub_key <- "xxxyyyzzz" # a pattern we would never imagine to use
    df <- data.frame(key=stub_key, stringsAsFactors=FALSE) # result template
    df["symbol"] = ""
    df["period"] = ""
    df["ts.on"]  = "days"
    df["lag"]  = 1
    df["drop"]  = 0
    df["t"]  = 0
    df["nu0"]  = 0
    df["theta"]  = 0
    df["convo"]  = 1
    df["beta.a"]  = 0
    df["mu.add"]  = 0
    df["breaks"]  = 200
    df["start.date"]  = ""
    df["end.date"] = ""

    df0 <- df # row template

    for (s in names(conf_all)) {
        cf1 <- conf_all[[s]]
        cf1$t <- eval(parse(text=cf1$t))
        df1 <- df0 # make a copy
        df1[1,"key"] <- s
        for(c in names(cf1)) df1[1,c] <- cf1[[c]]
        df <- rbind(df, df1)
    }
    df <- subset(df, key != stub_key) # remove stub
    if (!is.null(key)) {
        env <- environment()
        df <- subset(df, key==get('key',env))
    }
    rownames(df) <- df$key
    return(df)
}
### <---------------------------------------------------------------------->
#' @rdname lamp.qsl_fit_config
"lamp.qsl_fit_config_xtable" <- function(df)
{
    conf1 <- df
    drops <- c("key", "ts.on", "mu.add", "start.date", "end.date")
    conf <- conf1[ , !(names(conf1) %in% drops)]

    # calculate annualized vol
    m <- as.matrix(conf[ , c("t", "nu0", "theta", "convo", "beta.a")])
    v <- apply(m, 1, function(r) qsl_variance_analytic(
        t=r["t"], nu0=r["nu0"], theta=r["theta"], convo=r["convo"], beta.a=r["beta.a"]))
    v2 <- apply(m, 1, function(r) (r["nu0"]+6*r["theta"])^2+24*r["theta"]^2)
    conf$vol.a <- sqrt(v/m[,"t"])
    conf$vol.a2 <- sqrt(v2)

    conf$t <- ecd.rational(conf$t, as.character=TRUE, pref.denominator=c(12,52,250))
    conf$lag <- as.integer(conf$lag)
    conf$drop <- as.integer(conf$drop)
    conf$convo <- as.integer(conf$convo)
    conf$breaks <- as.integer(conf$breaks)

    tab <- xtable::xtable(conf)
    colnames(tab) <- c("symbol", "period", "lag", "drop",
                       "$t$", "$\\nu_0$", "$\\theta$", "$m$", "$\\beta_a$", "breaks",
                       "$\\sigma_{\\varLambda}$", "$\\sigma_{\\varLambda}^{sym}$")

    ncol <- dim(tab)[1]
    nligne <- dim(tab)[2]
    col <- paste(rep("r",ncol), collapse='')
    col <- paste('{', col ,'}', sep='')
    cat('\\begin{tabular}', col, sep='')
    cat('\\hline')
    print(tab, only.contents=TRUE, include.rownames=FALSE,
          sanitize.text.function=function(x){x})
    cat('\\end{tabular}')

}
### <---------------------------------------------------------------------->
