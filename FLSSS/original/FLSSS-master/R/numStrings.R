

addPositiveStrings = function(x)
{
  if(length(x) == 1) return(x)
  maxNum = max(nchar(x))
  x = sapply(x, function(u)
  {
    paste0(paste0(rep('0', maxNum - nchar(u)), collapse = ''), u)
  })
  names(x) = NULL
  x = strsplit(x, split = '')
  z = rowSums(as.data.frame(lapply(x, function(u) as.integer(u))))
  i = length(z)
  z = as.list(z)
  while(i >= 2L)
  {
    if(z[[i]] >= 10L)
    {
      q = z[[i]] %/% 10L
      r = z[[i]] %% 10L
      z[[i]] = r
      z[[i - 1]] = z[[i - 1]] + q
    }
    i = i - 1L
  }
  paste0(unlist(z), collapse = '')
}


subtract2positiveStrings = function (x, y)
{
  if (x == y) return('0')
  maxNdigit = max(nchar(x), nchar(y))
  x = paste0(paste0(rep('0', maxNdigit - nchar(x)), collapse = ''), x)
  y = paste0(paste0(rep('0', maxNdigit - nchar(y)), collapse = ''), y)
  asign = 1L
  if(x < y)
  {
    asign = -1L
    tmp = x; x = y; y = tmp
  }
  x = as.integer(strsplit(x, split = '')[[1]])
  y = as.integer(strsplit(y, split = '')[[1]])
  z = as.list(x - y)
  i = length(z)
  while (i >= 2L)
  {
    if(z[[i]] < 0)
    {
      z[[i]] = z[[i]] + 10L
      z[[i - 1]] = z[[i - 1]] - 1L
    }
    i = i - 1L
  }
  i = 1L
  while(i <= length(z))
  {
    if(z[i] != 0L) break
    i = i + 1L
  }
  if(i > length(z)) return('0')
  if(i > 1L) z = z[-(1:(i - 1L))]
  z = paste0(unlist(z), collapse = '')
  if(asign == -1L) z = paste0('-', z)
  z
}


asIntegerString = function(s)
{
  ssplit = strsplit(s, split = '[.]')
  Nd = max(unlist(lapply(ssplit, function(x) # max number of fraction digits.
  {
    if(length(x) == 1) return(0L)
    nchar(x[2])
  })))
  if(Nd != 0)
  {
    s = unlist(lapply(ssplit, function(x)
    {
      if(length(x) == 1) paste0(x, paste0(rep('0', Nd), collapse = ''))
      else
        paste0(x[1], x[2], paste0(rep('0', Nd - nchar(x[2])), collapse = ''))
    }))
  }
  list(NfracDigits = Nd, s = s)
}


addIntStrings = function(s)
{
  tmp = substr(s, 1, 1) != '-'
  addnums = s[tmp]
  if(length(addnums) != 0) Saddnums = addPositiveStrings(addnums)
  else Saddnums = '0'


  subnums = s[!tmp]
  if(length(subnums) != 0)
  {
    subnums = substr(subnums, 2L, max(nchar(subnums)))
    Ssubnums = addPositiveStrings(subnums)
  }
  else Ssubnums = '0'


  subtract2positiveStrings(Saddnums, Ssubnums)
}


addNumStrings = function(s)
{
  tmp = asIntegerString(s)
  s = tmp$s
  Nd = tmp$NfracDigits


  rst = addIntStrings(s)


  if(substr(rst, 1, 1) == '-')
  {
    rst = substr(rst, 2, nchar(rst))
    asign = -1L
  }
  else asign = 1L


  if(Nd != 0)
  {
    if(Nd < nchar(rst)) rst = paste0(
      substr(rst, 1, nchar(rst) - Nd), '.',
      substr(rst, nchar(rst) - Nd + 1L, nchar(rst)))
    else rst = paste0(
      "0.",  paste0(rep('0', Nd - nchar(rst)), collapse = ''), rst)
    tmp = strsplit(rst, split = '[.]')[[1]]
    frac = tmp[2]
    i = nchar(frac)
    while(i >= 1)
    {
      if(substr(frac, i, i) != '0') break
      i = i - 1L
    }
    if(i >= 1)
    {
      frac = substr(frac, 1, i)
      rst = paste0(tmp[1], '.', frac)
    }
    else rst = tmp[1]
  }


  if(asign == -1) rst = paste0('-', rst)
  rst
}























