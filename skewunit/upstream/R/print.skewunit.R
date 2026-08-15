print.skewunit<-function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
    type="est"
    if(!is.null(x$summary)) type="ch"
    name1=switch(x$family1, asin = "ArcSin", triang="triangular",
Uquad="U-quadratic", JSB="Johnson-SB", sbeta="symmetrical beta")
if(!is.null(x$family2)) name2=switch(x$family2, asin = "ArcSin", triang="triangular",
Uquad="U-quadratic", JSB="Johnson-SB", sbeta="symmetrical beta")
coeff=x$coefficients[,1]
names(coeff)=rownames(x$coefficients)
if(!is.null(x$family2))
{
    if (type=="est") {
            cat("Fitting the Family of Skew Distributions\nwith Bounded Support\nwith f ",
                name1," and G ",name2,sep = "")
    }
    if (type=="ch") {
            cat("Selecting a model in the Family of Skew Distributions\nwith Bounded Support.\n",
"The selected combination was f ",name1, " and G ", name2, sep="")
    }
}
if(is.null(x$family2))
{
    if (type=="est") {
            cat("Fitting the Family of Skew Distributions\nwith Bounded Support\nwith f ",
                name1,sep = "")
    }
    if (type=="ch") {
            cat("Selecting a model in the Family of Skew Distributions\nwith Bounded Support.\n",
"The selected combination was f ",name1, sep="")
    }
}
    cat("\n")
    cat("-----------------------------------------------------\n")
    cat("Coefficients:\n")
    print.default(format(coeff, digits = digits), print.gap = 2L,
        quote = FALSE)
    invisible(x)
}
