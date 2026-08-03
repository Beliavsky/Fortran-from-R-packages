## VT::15.09.2013 - this will render the output independent
##  from the version of the package
suppressPackageStartupMessages(library(rrcov))
library(MASS)

## VT::14.01.2020
##  On some platforms minor differences are shown - use
        ## IGNORE_RDIFF_BEGIN
        ## IGNORE_RDIFF_END

dodata <- function(method) {

    options(digits = 5)
    set.seed(101) # <<-- sub-sampling algorithm now based on R's RNG and seed

    tmp <- sys.call()
    cat("\nCall: ", deparse(substitute(tmp)),"\n")
    cat("===================================================\n")

    cat("\nData: ", "hemophilia\n")
    data(hemophilia)
    show(rlda <- Linda(as.factor(gr)~., data=hemophilia, method=method))
    show(predict(rlda))

    cat("\nData: ", "anorexia\n")
    data(anorexia)
    show(rlda <- Linda(Treat~., data=anorexia, method=method))
    show(predict(rlda))

    cat("\nData: ", "Pima\n")
    data(Pima.tr)
    show(rlda <- Linda(type~., data=Pima.tr, method=method))
    show(predict(rlda))

    cat("\nData: ", "Forest soils\n")
    data(soil)
    soil1983 <- soil[soil$D == 0, -2]       # only 1983, remove column D (always 0)

    ## This will not work within the function, of course
    ##  - comment it out
    ## IGNORE_RDIFF_BEGIN
    rlda <- Linda(F~., data=soil1983, method=method)
    ##  show(rlda)
    ## IGNORE_RDIFF_END
    show(predict(rlda))

    cat("\nData: ", "Raven and Miller diabetes data\n")
    data(diabetes)
    show(rlda <- Linda(group~insulin+glucose+sspg, data=diabetes, method=method))
    show(predict(rlda))

    cat("\nData: ", "iris\n")
    data(iris)
    if(method != "mcdA")
    {
        show(rlda <- Linda(Species~., data=iris, method=method, l1med=TRUE))
        show(predict(rlda))
    }

    cat("\nData: ", "crabs\n")
    data(crabs)
    show(rlda <- Linda(sp~., data=crabs, method=method))
    show(predict(rlda))

    cat("\nData: ", "fish\n")
    data(fish)
    fish <- fish[-14,]      # remove observation #14 containing missing value

    # The height and width are calculated as percentages
    #   of the third length variable
    fish[,5] <- fish[,5]*fish[,4]/100
    fish[,6] <- fish[,6]*fish[,4]/100

    ## There is one class with only 6 observations (p=6). Normally
    ##  Linda will fail, therefore use l1med=TRUE.
    ##  This works only for methods mcdB and mcdC

    table(fish$Species)
    if(method != "mcdA")
    {
        ## IGNORE_RDIFF_BEGIN
        rlda <- Linda(Species~., data=fish, method=method, l1med=TRUE)
        ## show(rlda)
        ## IGNORE_RDIFF_END
        show(predict(rlda))
    }

    cat("\nData: ", "pottery\n")
    data(pottery)
    show(rlda <- Linda(origin~., data=pottery, method=method))
    show(predict(rlda))

    cat("\nData: ", "olitos\n")
    data(olitos)
    if(method != "mcdA")
    {
        ## IGNORE_RDIFF_BEGIN
        rlda <- Linda(grp~., data=olitos, method=method, l1med=TRUE)
        ## show(rlda)
        ## IGNORE_RDIFF_END
        show(predict(rlda))
    }

    cat("===================================================\n")
}


## -- now do it:
dodata(method="mcdA")
dodata(method="mcdB")
dodata(method="mcdC")
dodata(method="mrcd")
dodata(method="ogk")
#dodata(method="fsa")
