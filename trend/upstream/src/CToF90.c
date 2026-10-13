#include <R.h>
#include <Rmath.h>
/*
    Copyright (C) 2017-2018 Thorsten Pohlert

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// This is to get the state of RNG
void F77_SUB(rndstart)(void) { GetRNGstate(); }

// This is to put the state of RNG
void F77_SUB(rndend)(void) { PutRNGstate(); }

// This is the random generator for standard normal variates 
double F77_SUB(normrand)(void) { return norm_rand(); }
