/*
    kz.h
    Copyright (C) 2005  Brian D. Close

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#include <Rdefines.h>

SEXP kza(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP kz(SEXP, SEXP, SEXP);
SEXP kz1d(SEXP, SEXP, SEXP);
SEXP kz2d(SEXP, SEXP, SEXP);
void kza1d_test(double *, long *, long *, long *, long *);
SEXP kza1d(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP kza2d(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
void differenced(double *, double *, double *, long, int);
void R_differenced(SEXP, SEXP, SEXP, int);
double R_maximum(SEXP);
SEXP kzsv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_kzsv(SEXP, SEXP, SEXP, SEXP, SEXP);
void copyArray(SEXP, SEXP);
SEXP kzs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
//SEXP kzm(SEXP, SEXP, SEXP, SEXP, SEXP);
 
#ifndef MIN
#define MIN(x,y) ((x)<(y)?(x):(y))
#endif
#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#endif
#ifndef ABS
#define ABS(x)  ((x<0)?(-x):(x))
#endif
