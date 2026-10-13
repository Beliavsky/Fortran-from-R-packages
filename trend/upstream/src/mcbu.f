      subroutine mcbu (stat, n, m, pval)
C
C     Copyright (C) 2017-2018 Thorsten Pohlert
C
C     This program is free software: you can redistribute it and/or modify
C     it under the terms of the GNU General Public License as published by
C     the Free Software Foundation, either version 3 of the License, or
C     (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see <http://www.gnu.org/licenses/>.
C
C     DESCRIPTION
C     This subroutine estimates the p-value of the Buishand U test
C     with Monte-Carlo simulation 

      implicit none
      double precision, intent(in) :: stat
      integer, intent(in) :: n, m
      double precision, intent(out) :: pval
C     local variables
      integer :: j, i, k
      double precision, parameter :: zero = 0.0d0, one = 1.0d0 
      double precision, parameter :: two=2.0d0
      double precision, dimension(n) :: sk, x
      double precision, dimension(m) :: statm
      double precision :: mu, sigma, tmp
C     Intrinsic functions:
      intrinsic :: sum, sqrt
C     External functions:
      double precision :: normrand, getpval
      external :: normrand, getpval


      call rndstart()
C     Monte Carlo loop
      do j = 1, m
         
C     user interrupt
         call rchkusr()

C     Get normal random deviates
         do i = 1, n
            x(i) = normrand()
         end do
         
C     mean(x)
         mu = sum(x) / dble(n)
         
C     sd(x)
         tmp = zero
         do i = 1, n
            tmp = tmp + (x(i) - mu)**two
         end do
         sigma = sqrt(tmp / dble(n))

C     get Sk
         do k = 1, n
            tmp = zero
            do i = 1, k
               tmp = tmp + x(i) - mu
            end do
            sk(k) = tmp
         end do


         tmp = zero
         do k = 1, n-1
            tmp = tmp + (sk(k) / sigma)**two
         end do
         statm(j) = one / (dble(n) * 
     *        (dble(n) + one)) * tmp  
  
      end do
      call rndend()

      pval = getpval(statm, stat, m)
      end subroutine
