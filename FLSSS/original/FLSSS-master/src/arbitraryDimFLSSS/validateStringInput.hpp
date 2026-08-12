#pragma once


// #include <Rcpp.h>
// using namespace Rcpp;
// Test kit in ../arbitraryDimFLSSS/validateStringInput.cpp
inline bool validateNumStrings(StringVector x)
{
  std::string a = "";
  for(int i = 0, iend = x.size(); i < iend; ++i)
  {
    const char *s = &x[i][0];
    for(int k = 0, kend = x[i].size(); k < kend; ++k)
    {
      char c = s[k];
      if(c != '.' and c != '-' and (c < 48 or c > 57)) { a = x[i]; break; }
      if(c == '-' and k != 0) { a = x[i]; break; }
      if(c == '.' and ((k == kend - 1 or k == 0) or (
        s[k - 1] < 48 or s[k - 1] > 57)) ) { a = x[i]; break; }
    }
    if(a != "") break;
  }
  if(a != "")
  {
    std::string msg = a + " is not in standard format. \
Consider setting 'options(scipen = 999)' in R before casting numerics to strings.\n";
    Rcpp::stop(msg);
    return false;
  }
  return true;
}




