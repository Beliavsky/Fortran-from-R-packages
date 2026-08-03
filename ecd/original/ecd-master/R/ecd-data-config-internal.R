#' Read sample data config
#' 
#' Read sample data config. This is used in ecd.data() and related tests.
#'
#' @param symbol Optional, if provided, only return config related to it.
#'               default = NULL
#' @param extdata_dir optionally specify user's own extdata folder
#'                    containing user's own data_conf.yml
#' @param filename character, user's own config file name, default is NULL.
#'
#' @return The data.frame object for the config
#'
#' @keywords sample data
#'
#' @examples
#' c <- .ecd.data_config()
#' c <- .ecd.data_config("dji")
#' @noRd
### <======================================================================>
".ecd.data_config" <- function(symbol=NULL, extdata_dir=NULL, filename=NULL)
{
    if (is.null(filename)) filename <- "data_conf.yml"
    conf_file <- .ecd.locate_file(filename, extdata_dir=extdata_dir)
    # read the conf file
    if (! file.exists(conf_file)) {
        stop(paste("conf_file does not exist:", conf_file))
    }
    conf_all <- yaml.load_file(conf_file)
    
    m0 <- function() {
        cols = c("symbol", "cols", "col_in", "date_format", "test_date", "test_val")
        matrix(vector(), 0, 6, dimnames=list(c(), cols)) 
    }
    add <- function(m, c) {
        if (is.null(symbol) || (symbol == c[1])) {
            m2 <- rbind(m, c)
            rownames(m2) <- NULL
            m2
        }else{
            m
        }
    }
    
    mx <- m0()
    for (s in names(conf_all)) {
        cf1 <- conf_all[[s]]
        mx <- add(mx, c(s, cf1$cols, cf1$col_in, cf1$date_format, cf1$test_date, cf1$test_val))
    }
    df <- data.frame(mx, stringsAsFactors=FALSE)
    rownames(df) <- df$symbol
    return(df)
}
### <---------------------------------------------------------------------->
