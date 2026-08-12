

LE = function(x, y)
{
  i = length(x)
  while(i >= 1L)
  {
    if(x[i] < y[i]) return(T)
    if(x[i] > y[i]) return(F)
    i = i - 1L
  }
  return(T)
}


GE = function(x, y) { return( LE(y, x) ) }


updateLB = function(len, v, target, me, lb, ub)
{
  if(len < 2) stop("len should be no less than 2")
  target = target - me
  # cat("target = ", target, "\n")
  for(b in lb[1]:ub[1])
  {
    s = colSums(v[c(b, ub[2:len]), , drop = F])
    if ( GE(s, target) ) break
    else
    {
      if ( b >= ub[1] ) return ( integer(0) )
    }
  }
  lb[1] = b
  # cat(lb)


  for(i in 2:len)
  {
    b = max(lb[i - 1] + 1L, lb[i])
    if(i < len) ubind = ub[(i + 1):len]
    else ubind = integer(0)
    for(b in b:ub[i])
    {
      btrail = pmin( (b - i + 1L):b, ub[1:i] )
      sumInd = c(btrail, ubind)
      s = colSums(v[sumInd, , drop = F])
      if ( GE(s, target) ) break
      else
      {
        if ( b >= ub[i] ) return ( integer(0) )
      }
    }
    lb[i] = b
  }


  return(lb)
}


updateUB = function(len, v, target, me, lb, ub)
{
  # cat("updateUB\n")
  if(len < 2) stop("len should be no less than 2")
  target = target + me
  for(b in ub[len]:lb[len])
  {
    s = colSums(v[c(b, lb[1:(len - 1L)]), , drop = F])
    if ( LE(s, target) ) break
    else
    {
      if ( b <= lb[len] ) return ( integer(0) )
    }
  }
  ub[len] = b


  for(i in (len - 1L):1)
  {
    b = min(ub[i + 1L] - 1L, ub[i])
    if ( i > 1L ) lbind = lb[ 1:(i - 1L) ]
    else lbind = integer(0)
    for (b in b:lb[i])
    {
      btrail = pmax( b:(b + len - i), lb[i:len] )
      sumInd = c(btrail, lbind)
      s = colSums(v[sumInd, , drop = F])
      if ( LE(s, target) ) break
      else
      {
        if ( b <= lb[i] ) return ( integer(0) )
      }
    }
    ub[i] = b
  }


  return(ub)
}


findBoundsR = function(len, v, target, me = 0L, lb, ub)
{
  firstIter = T
  foundOrNot = 0L
  while(T)
  {
    # cat("LB =", lb - 1L, "\n")
    # cat("UB =", ub - 1L, "\n\n")


    LBnew = updateLB ( len, v, target, me, lb, ub )
    if( length(LBnew) == 0 ) return (
      list( foundOrNot = 0L, lb = integer(0), ub = integer(0) ) )
    if( !firstIter && sum( abs( LBnew - lb ) ) == 0 )
    {
      lb = LBnew
      foundOrNot = 1L
      break
    }
    # cat("LBnew = ", LBnew, "\n")
    firstIter = F
    lb = LBnew


    # cat("LB =", lb - 1L, "\n")
    # cat("UB =", ub - 1L, "\n\n")


    UBnew = updateUB ( len, v, target, me, lb, ub )
    if ( length(UBnew) == 0 ) return (
      list( foundOrNot = 0L, lb = integer(0), ub = integer(0) ) )
    if( sum( abs( UBnew - ub ) ) == 0)
    {
      ub = UBnew
      foundOrNot = 1L
      break
    }
    ub = UBnew
  }
  if(sum(abs(ub - lb)) == 0) foundOrNot = 2L
  list( foundOrNot = foundOrNot, lb = lb, ub = ub)
}




# Tests.
if(F)
{

  # 2^64 = "18446744073709551616"
  #        "18446744073709551615"


  Rcpp::sourceCpp("src/arbitraryDimFLSSStests/findBound.cpp", verbose = T)


  options(scipen = 999)
  set.seed(42)


  while(T){


    d = sample(10, 1) + 1L
    len = sample(15, 1) + 1L
    N = sample(60, 1) + len
    roundN = sample(5, 1) + 1L


    # d = 5L
    # len = 5L
    # N = 20L
    # roundN = 5L


    X = matrix(round(runif(d * N, -1, 1), roundN), nrow = N)
    sol = sort(sample(N, len))
    target = round(colSums(X[sol, , drop = F]), roundN)
    Xc = matrix(as.character(X), nrow = N)
    targetc = as.character(target)


    tmpX = round(X * 10 ^ roundN)
    tmpXmin = apply(tmpX, 2, function(x) min(x))
    tmpX = t(t(tmpX) - tmpXmin)
    tmpTarget = round(target * 10 ^ roundN)
    tmpTarget = round(tmpTarget - tmpXmin * len)
    tmpXrank = apply(tmpX, 2, function(x) rank(x, ties.method = "first"))
    colOdr = order(colSums(t(tmpXrank) %*% tmpXrank))
    Xc = Xc[, colOdr, drop = F]
    targetc = targetc[colOdr]
    tmpX = tmpX[, colOdr, drop = F]
    tmpTarget = as.integer(tmpTarget[colOdr])
    rowOdr = do.call(order, as.data.frame(tmpX[, ncol(tmpX):1, drop = F]))
    tmpX = matrix(as.integer(tmpX[rowOdr, , drop = F]), nrow = nrow(tmpX))
    Xc = Xc[rowOdr, , drop = F]
    tmpsol = sort(match(sol, rowOdr))


    rst = findBound(Xval = Xc, tsum = targetc, len = len, maxCore = 1)
    # rst = findBound(Xval = matrix(as.character(tmpX), ncol = d),
    #                 tsum = as.character(tmpTarget), len = len, maxCore = 1)


    if( (sum(abs(rst$order - 0:(N - 1L))) +
         sum(abs(rst$colOrder - 0:(d-1L)))) != 0 ) next


    lb = 1:len
    ub = (N - len + 1L):N
    correctRst = findBoundsR(len, tmpX, tmpTarget, me = 0L, lb, ub)


    if(rst$foundOrNot == 0 && rst$foundOrNot != 0)
    {
      cat("rst$foundOrNot == 0 && rst$foundOrNot != 0\n")
      break
    }
    if(rst$foundOrNot != 0 && rst$foundOrNot == 0)
    {
      cat("rst$foundOrNot == 0 && rst$foundOrNot != 0\n")
      break
    }
    if(rst$foundOrNot == 0 && rst$foundOrNot == 0) next


    err = sum(abs(unlist(rst[1:3]) - unlist(correctRst)))


    if(err != 0) break


  }





















}




























