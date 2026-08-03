getEER <-function(Areas=c("US","JP"),Freq="Monthly",Type="Real",Basket="Broad"){
  Sys.setenv('_R_CHECK_SYSTEM_CLOCK_' = 0)
  cat("Getting", Freq, Type, Basket, "Effective Exchange Rate of",Areas,"\n")

  if (Freq=="Daily" & Type=="Real") stop("Daily frequency does not have REER")

  ID=paste(substr(Freq,1,1),substr(Type,1,1),substr(Basket,1,1),sep=".")

  #Daily, Nominal and Broad
  main=paste0("https://stats.bis.org/api/v2/data/dataflow/BIS/WS_EER/1.0/", ID,".")
  END="?format=csv"
  #areas=c("AR","AU","BR","CA","CN","DE","FR","GB","ID","IN","IT","JP","KR","MX","RU","SA","TR","US","ZA")
  AREA=paste0(main, Areas, END)
  DNB_G20=NULL
  for (i in seq(AREA)) {
    tmp=read.csv(AREA[i])[,c("TIME_PERIOD","OBS_VALUE")]

    if(Freq=="Monthly") {DATE=paste0(tmp[,1],"-01")
    } else {DATE=tmp[,1]}

    dat=as.matrix(tmp[,-1])
    rownames(dat)=DATE
    #dat=timeSeries::as.timeSeries(dat)

    DNB_G20=cbind(DNB_G20, dat)
    cat(Areas[i],"\n")
  }
  colnames(DNB_G20)=Areas

  return(DNB_G20)

}


getFed <- function (var.name="UNRATE", freq="Monthly"){

  mainpath="https://fred.stlouisfed.org/graph/fredgraph.csv?&id="

  url=paste0(mainpath,var.name,"&fq=",freq)
  dat0=read.csv(url)
  dat=dat0[,2,drop=FALSE]
  rownames(dat)=as.character(dat0[,1])
  #dat=timeSeries::as.timeSeries(dat)

  return(dat)
}



getFrench.Factors<-function(filename="F-F_Research_Data_Factors") {

  Url="http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/"

  ##Step 1. Download .zip dataset

#  if (filename %in% c("F-F_Research_Data_Factors_weekly",
#                       "F-F_Research_Data_5_Factors_2x3"))
#    { fileCSV=paste0(filename,".csv") } else {  fileCSV=paste0(filename,".csv") }
    fileCSV=paste0(filename,".csv")
    fileZIP=paste0(filename,"_CSV.zip")
    Link=paste0(Url,fileZIP)
    temp = tempfile(tmpdir=getwd(),fileext = ".zip")
    download.file(Link,temp)



  ##Step 2. Retrieve .csv data


  if (filename %in% c("F-F_Research_Data_Factors_weekly",
             "F-F_Research_Data_Factors_daily",
             "F-F_Research_Data_5_Factors_2x3_daily",
             "F-F_Momentum_Factor_daily"))  {

    tmp0=read.csv(unz(temp, fileCSV),header=FALSE, blank.lines.skip = FALSE)
    S=which(grepl("19", tmp0[1:30,1]))[1];S
    ff.factor0 = read.csv(unz(temp, fileCSV),skip=S-2,header=TRUE, blank.lines.skip = FALSE)

    end.date=grep("Copyright ",ff.factor0[,1])-2
    if (length(end.date)==0) {end.date=nrow(ff.factor0)} else {
      end.date=end.date
      }

    ff.factor=ff.factor0[1:end.date,]
    colnames(ff.factor)=c("Dates", colnames(ff.factor0[-1]))


dat=ff.factor
table.names="Only one table"

  } else if (filename %in% c("F-F_Research_Data_Factors",
                             "F-F_Research_Data_5_Factors_2x3")) {
     #Both datasets have two frequencies: monthly and annual


    tmp0=read.csv(unz(temp, fileCSV), header=FALSE, blank.lines.skip = FALSE)
#    head(tmp0,20)
    S=which(grepl("19", tmp0[1:30,1]))[1]

    ff.factor0 = read.csv(unz(temp, fileCSV), skip=S-2,header=TRUE, blank.lines.skip = FALSE)
#    head(ff.factor0,15)
#    tail(ff.factor0)
    endID=grep("[Aa-zZ] ",ff.factor0[,1])

    if (length(endID)==2) {
    ff.factor.month=ff.factor0[1:(endID[1]-2),]
    ff.factor.annual=ff.factor0[(endID[1]+2):(endID[2]-2),]
      } else if (length(endID)==1) {

    ff.factor.month=ff.factor0[1:(endID[1]-2),]
    ff.factor.annual=ff.factor0[(endID[1]+2):nrow(ff.factor0),]
      }

    table.names=trimws(c("Monthly Factors",ff.factor0[,1][which(grepl("[Aa-zZ]", ff.factor0[,1]))][1]))
    colnames(ff.factor.month)=colnames(ff.factor.annual)=c("Dates",colnames(ff.factor.month)[-1])
#    head(ff.factor.month);tail(ff.factor.month)
#    head(ff.factor.annual);tail(ff.factor.annual)

    dat=list()
    dat[[1]]=ff.factor.month
    dat[[2]]=ff.factor.annual
    names(dat)=table.names


  } else if (filename == "F-F_Momentum_Factor") {

    tmp0=read.csv(unz(temp, fileCSV), header=FALSE, blank.lines.skip = FALSE)
    head(tmp0,20)
    S=which(grepl("19", tmp0[1:30,1]))[1]
    ff.factor0 = read.csv(unz(temp, fileCSV), skip=S-2,header=TRUE, blank.lines.skip = FALSE)
    end.date=grep("Annual Factors",ff.factor0[,1])-2

    ff.factor.month=ff.factor0[1:end.date,]

    colnames(ff.factor.month)=c("Dates",colnames(ff.factor.month)[-1])
    head(ff.factor.month)
    tail(ff.factor.month)

    start.annual=grep("Annual Factors",ff.factor0[,1])+3
    end.annual=grep("Copyright",ff.factor0[,1])
    if (length(end.annual)==0) {end.annual=nrow(ff.factor0)} else {
      end.annual=end.annual-2
    }


    ff.factor.annual=ff.factor0[start.annual:end.annual,]

    colnames(ff.factor.annual)=c("Dates",colnames(ff.factor.annual)[-1])
    head(ff.factor.annual)
    tail(ff.factor.annual)

    table.names=trimws(c("Monthly Factors",ff.factor0[,1][which(grepl("[Aa-zZ]", ff.factor0[,1]))][1]))
    dat=list()
    dat[[1]]=ff.factor.month
    dat[[2]]=ff.factor.annual
    names(dat)=table.names
  }

  file.remove(temp)
  return(list(data=dat,table.names=table.names,file.name=filename))


}




getFrench.Portfolios<-function(filename="Portfolios_Formed_on_ME") {
  url="http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/"

  ##Step 1. Download .zip dataset


  if (grepl("daily",filename)) {
    fileCSV=paste0(filename,".csv")
    fileCSV=gsub(fileCSV,pattern = "daily",replacement = "Daily")

  } else {  fileCSV=paste0(filename,".csv") }


  {
    fileZIP=paste0(filename,"_CSV.zip")
    Link=paste0(url,fileZIP)
    temp = tempfile(tmpdir=getwd(),fileext = ".zip")
    download.file(Link,temp)

  }


  ##Step 2. Retrieve .csv data

  tmp0=read.delim(unz(temp, fileCSV),header=FALSE,blank.lines.skip = FALSE)
  S=which(grepl("19", tmp0[1:50,1]))[1]

  ff.factor0 = read.csv(unz(temp, fileCSV),skip=S-3,header=FALSE, blank.lines.skip = FALSE)
  head(ff.factor0)
  tail(ff.factor0)
  endOfData=grep("For portfolios formed ",ff.factor0[,1])[1]-3
    if (is.na(endOfData)) {
  endOfData=grep("Copyright ",ff.factor0[,1])-2
  }
  if (length(endOfData)==0) {ff.factor0=ff.factor0} else {
    ff.factor0=ff.factor0[1:(endOfData),]
  }
  tail(ff.factor0)

  table.names=trimws(ff.factor0[,1][which(grepl("[Aa-zZ]", ff.factor0[,1]))])
  colNames=c("Dates",as.character(ff.factor0[2,])[-c(1)])

  table.names.where=which(grepl("[Aa-zZ]", ff.factor0[,1]))

  dat=list()
  for (i in 1:length(table.names)) {#i=5
    if (i==length(table.names)) {
      dat0=ff.factor0[(table.names.where[i]+2):nrow(ff.factor0),]
    } else {

      dat0=ff.factor0[(table.names.where[i]+2):(table.names.where[i+1]-3),]

    }
    colnames(dat0)=colNames

    dat[[i]]=dat0

  }

  names(dat) <- trimws(table.names)

  file.remove(temp)
  return(list(data=dat,table.names=table.names,file.name=filename))

}




.seriesPlot <-
  function(x, type = "l", col = "indianred2", MAIN=MAIN,grid = TRUE, box = TRUE, rug = TRUE, ...)

  {
    # Series Plots:
      plot(x, type = type, col = col, ...)
      title(main = MAIN)
      if(grid) grid()
      if(box) box()
      if(rug) rug(as.vector(x), ticksize = 0.01, side = 2, quiet = TRUE)


    # Return Value:
    invisible()
  }


