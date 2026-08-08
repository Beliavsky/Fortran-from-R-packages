
// common_R.hpp -- this file is part of trustOptim, a nonlinear optimization package
// for the R statistical programming platform.
//
// Copyright (C) 2013-2020 Michael Braun
//
// This Source Code Form is subject to the license terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, you can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef COMMON_R
#define COMMON_R

#ifndef TRUST_COUT
#define TRUST_COUT Rcpp::Rcout
#endif

#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstring>


void throw_exception(const std::string reason,
		     const std::string file,
		     const int line)
{
  
  std::ostringstream oss;
  
  oss << "\nException in File " << file;
  oss << "  at Line " << line <<".\n";
  oss << "Reason : " << reason << "\n";
  
  std::string message = oss.str();

  Rcpp::stop(message);

};

inline void check_interrupt() {
  Rcpp::checkUserInterrupt();
}

inline bool my_ret_bool(bool x) {return(x);}

template<typename T>
bool my_finite(const T& x) {
  return( (std::abs(x) <= __DBL_MAX__ ) && ( x == x ) );
}



enum MB_Status {
  
  SUCCESS  = 0, 
  FAILURE  = -1,
  CONTINUE = -2,  /* iteration has not converged */
  ENOPROG  = 1,  /* iteration is not making progress towards solution */
  ETOLF    = 3,  /* cannot reach the specified tolerance in F */
  ETOLX    = 4,  /* cannot reach the specified tolerance in X */
  ETOLG    = 5,  /* cannot reach the specified tolerance in gradient */
  EMAXITER = 6,  /* exceeded max number of iterations */
  EBADLEN  = 7,  /* matrix, vector lengths are not conformant */
  ENOTSQR  = 8,  /* matrix not square */
  ESING    = 9,  /* apparent singularity detected */
  ENOMOVE  = 10,  /* no movement on last iteration */
  FAILEDCG = 11,  /* CG step returned a bad proposed step.*/
  MOVED = 12,     /* Iteration resulted in a successful move */
  CONTRACT = 13,   /* Iteration resulted in contracting the trust region */
  EXPAND = 14,     /* Iteration resulted in a successful move and expanding the trust region */
  UNKNOWN = 15,    /* Unknown/unspecified status */
  ENEGMOVE = 16   /* Negative prediction */
};


const char * MB_strerror (const MB_Status & code) {
  
  
  switch (code)
    {
    case SUCCESS:
      return "Success";
    case FAILURE:  
      return "Unspecified failure";
    case CONTINUE: 
      return "Continuing";
    case ENOPROG:  
      return "Not making any progress";
    case ETOLF:    
      return "Cannot reach tolerance in F";
    case ETOLX:    
      return "Cannot reach tolerance in X";
    case ETOLG:    
      return "Radius of trust region is less than stop.trust.radius";
    case EMAXITER: 
      return "Exceeded max iterations";
    case EBADLEN:  
      return "Matrix, vector lengths not conformant";
    case ENOTSQR:  
      return "Matrix is not square";
    case ESING:    
      return "Matrix is apparently singular";
    case FAILEDCG:
      return "CG failed";
    case MOVED:
      return "Continuing";
    case CONTRACT:
      return "Continuing - TR contract";
    case EXPAND:
      return "Continuing - TR expand";
    case UNKNOWN:
      return "Unspecified status";
    case ENOMOVE:
      return "No proposed movement";
    case ENEGMOVE:
      return "Negative predicted move";
    default:
      return "Unspecified error";
    }
}



#endif
