/*
	Data.c
		initialized data for Cuba
		by Thomas Hahn
		last modified 30 Oct 20 th
*/

/* This #ifndef is to selectively exclude the definition of cubafun_ in stddecl
   so that R compiler flags used don't generate a warning about unusued variable
   for cubafun_.
*/

#ifndef _data_c_source
#define _data_c_source
#endif

#include "stddecl.h"

/* Commented out for gcc 10.0 -fno-common default */
/* coreinit cubafun_; */

int cubaverb_ = uninitialized;

#ifdef HAVE_FORK
corespec cubaworkers_ = {
  uninitialized, uninitialized, 
  uninitialized, uninitialized };
#endif

