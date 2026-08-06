polynomial <- function(n, coeficientes=numeric()){
  aux1 <- length(coeficientes)
  u <- (1:n)/n

  eta <- matrix(NA, ncol = aux1, nrow = n)
  for(j in 1:(aux1)){
    eta[,j] <- coeficientes[j]*u^(j-1)
  }
  result <- numeric(n)
  for (i in 1:n) {
    result[i] <- sum(eta[i, ], na.rm = TRUE)
  }
  return(result)
}

NoLineal <- function(n, coeficientes=numeric(), exponentes=numeric()){
  aux1 <- length(coeficientes)
  u <- (1:n)/n

  eta <- matrix(NA, ncol = aux1, nrow = n)
  for(j in 1:(aux1)){
    eta[,j] <- coeficientes[j]*u^(exponentes[j])
  }
  result <- numeric(n)
  for (i in 1:n) {
    result[i] <- sum(eta[i, ], na.rm = TRUE)
  }
  return(result)
}

trigonometric <- function(n, coeficientes = numeric(), trig = c("sin","cos"), arg = "u"){

  aux1 <- length(coeficientes)
  u <- (1:n)/n

  eta <- matrix(NA, ncol = aux1, nrow = n)
  eta[,1] <- coeficientes[1]
  arg <- parse(text = arg)

  for(j in 2:(aux1)){
    if (trig == "cos"){
      eta[,j] <- coeficientes[j]*cos(eval(arg, envir = list(u = u)))
    }else{
      eta[,j] <- coeficientes[j]*sin(eval(arg, envir = list(u = u)))
    }
  }

  result <- numeric(n)
  for (i in 1:n) {
    result[i] <- sum(eta[i, ], na.rm = TRUE)
  }
  return(result)
}

checkInput <- function(type=c("polynomial","NoLineal","trigonometric")) {
  if(length(type) != 3) stop("invalid type length")
  if(type[1] != "polynomial" && type[1] != "NoLineal" && type[1] != "trigonometric") stop("invalid type for c")
  if(type[2] != "polynomial" && type[2] != "NoLineal" && type[2] != "trigonometric") stop("invalid type for a")
  if(type[3] != "polynomial" && type[3] != "NoLineal" && type[3] != "trigonometric") stop("invalid type for b")
}

