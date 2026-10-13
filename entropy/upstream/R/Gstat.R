
### Gstat.R  (2021-04-25)
###
###    Computation of the G statistic and the chi-squared statistic
###    and corresponding p-values
###
### Copyright 2021 Korbinian Strimmer
###
###
### This file is part of the `entropy' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


# internal function for computing one-sided p-values
pvt.chisq.pval = function(stat, df, unit=c("log", "log2", "log10"))
{
    unit = match.arg(unit)

    if (unit == "log") pval = 1-pchisq(stat, df) 
    if (unit == "log2") pval = 1-pchisq(stat*log(2), df) 
    if (unit == "log10") pval = 1-pchisq(stat*log(10), df) 

    return( pval )
}

Gstat = function(y, freqs, unit=c("log", "log2", "log10"))
{
   unit = match.arg(unit)

   n = sum(y)
   df = length(y)-1

   stat = 2*n*KL.empirical(y1=y, y2=freqs, unit=unit)
   pval = pvt.chisq.pval(stat, df, unit=unit)

   return( list(stat=stat, df=df, pval=pval) )
}

chi2stat = function(y, freqs, unit=c("log", "log2", "log10"))
{
   unit = match.arg(unit)

   n = sum(y)
   df = length(y)-1
   
   stat = n*chi2.empirical(y1=y, y2=freqs, unit=unit)
   pval = pvt.chisq.pval(stat, df, unit=unit)

   return( list(stat=stat, df=df, pval=pval) )
}


Gstatindep = function(y2d, unit=c("log", "log2", "log10"))
{
   unit = match.arg(unit)

   n = sum(y2d)
   d1 = nrow(y2d)
   d2 = ncol(y2d)
   df = d1*d2-d1-d2+1

   stat = 2*n*mi.empirical(y2d, unit=unit)
   pval = pvt.chisq.pval(stat, df, unit=unit)

   return( list(stat=stat, df=df, pval=pval) )
}


chi2statindep = function(y2d, unit=c("log", "log2", "log10"))
{
   unit = match.arg(unit)

   n = sum(y2d)
   d1 = nrow(y2d)
   d2 = ncol(y2d)
   df = d1*d2-d1-d2+1

   stat = n*chi2indep.empirical(y2d, unit=unit)
   pval = pvt.chisq.pval(stat, df, unit=unit)

   return( list(stat=stat, df=df, pval=pval) )
}

