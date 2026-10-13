      subroutine mcsnht(stat, n, m, pval)
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
C     This subroutine calculates the pvalue of the SNH-Test
C     with a Monte-Carlo simulation
C
      implicit none
      double precision, intent(in) :: stat
      integer, intent(in) :: n, m
      double precision, intent(out) :: pval

C     local variables
      integer :: j, i, k
      double precision, parameter :: zero = 0.0d0, one = 1.0d0 
      double precision, parameter :: two=2.0d0
      double precision, dimension(n) :: x
      double precision, dimension(:),allocatable :: Tk
      double precision, dimension(m) :: T
      double precision :: mu, sigma, z1, z2, tmp
C     Intrinsic functions:
      intrinsic :: sum, sqrt
C     External functions:
      double precision :: normrand, getpval
      external :: normrand, getpval

      allocate(Tk(n-1))


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

C     get Tk
         do k = 1, n - 1
            z1 = zero
            z2 = zero
            do i = 1, k
               z1 = z1 + (x(i) - mu) / sigma
            end do
            do i = k + 1, n
               z2 = z2 + (x(i) - mu) / sigma
            end do
            z1 = z1 / dble(k)
            z2 = z2 / dble((n - k))

            
            Tk(k) = dble(k) * z1**two +
     *           (dble(n - k)) * z2**two
         end do

         T(j) = maxval(Tk)
  
      end do
      call rndend()
      
      pval = getpval(T, stat, m)

      deallocate(Tk)
      end subroutine
