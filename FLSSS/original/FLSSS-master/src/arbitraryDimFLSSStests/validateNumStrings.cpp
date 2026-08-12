#include <Rcpp.h>
using namespace Rcpp;
#include "../arbitraryDimFLSSS/validateNumStrings.hpp"


// [[Rcpp::export]]
bool validateNumStringsTest(StringVector x)
{
  return validateNumStrings(x);
}





