rate.conv = function(rate, conv=1, type="interest", nom=1){

  all=list(rate,conv,type,nom)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All variables must be of length 1.")
  #NA
  if(any(lapply(all,is.na)==T)) stop("Cannot input any variables as NA.")
  #Numeric
  num2=list(rate,conv,nom)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("rate, conv, and nom must be numeric.")
  #Infinite
  if(any(c(rate,conv,nom)==Inf)){
    stop('"rate", "conv", and "nom" cannot be infinite.')}
  #Positive
  if(rate <= 0 | conv <= 0 | nom <= 0){
    stop('"rate", "conv", and "nom" must be positive')}

  if(type != "interest" & type != "discount" & type != "force"){
    stop('type must be one of "interest", "discount", or "force"')}


  if(type == "interest") eff.i=(1+rate/conv)^(conv)-1
  if(type == "discount") eff.i=(1-rate/conv)^(-conv)-1
  if(type == "force") eff.i=exp(rate)-1
  eff.d=eff.i/(1+eff.i)
  force=log(1+eff.i)
  out = matrix(c(eff.i, eff.i, force), nrow=1, ncol=3,byrow=F,
               dimnames=list("eff",c("i", "d", "delta")))

  if(conv != 1) {
    conv.i=((1+eff.i)^(1/conv)-1)*conv
    conv.d=((1-eff.d)^(1/conv)-1)*-conv
    out = matrix(c(eff.i,eff.d,force,conv.i,conv.d,
                   NA), nrow=2, ncol=3, byrow=T,
                 dimnames=list(c("eff",sprintf("nom(%s)",conv)),
                               c("i", "d", "delta")))
  }

  if(nom != conv & nom != 1) {
    nom.i=((1+eff.i)^(1/nom)-1)*nom
    nom.d=((1-eff.d)^(1/nom)-1)*-nom
    out = rbind(out,matrix(c(nom.i,nom.d,NA),ncol=3,
                           dimnames=list(c(sprintf("nom(%s)",nom)),
                                         c("i", "d", "delta"))))
  }
  return(out)
}

