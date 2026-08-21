! SPDX-License-Identifier: GPL-2.0-or-later
module quadprog_core_mod
  use quadprog_kinds, only: dp
  implicit none
  private
  public :: qpgen2, qpgen1
contains
!
!  Copyright (C) 1995-2010 Berwin A. Turlach <Berwin.Turlach@gmail.com>
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
!  USA.
!
!  this routine uses the Goldfarb/Idnani algorithm to solve the
!  following minimization problem:
!
!        minimize  -d^T x + 1/2 *  x^T D x
!        where   A1^T x  = b1
!                A2^T x >= b2
!
!  the matrix D is assumed to be positive definite.  Especially,
!  w.l.o.g. D is assumed to be symmetric.
!  
!  Input parameter:
!  dmat   nxn matrix, the matrix D from above (dp)
!         *** WILL BE DESTROYED ON EXIT ***
!         The user has two possibilities:
!         a) Give D (ierr=0), in this case we use routines from LINPACK
!            to decompose D.
!         b) To get the algorithm started we need R^-1, where D=R^TR.
!            So if it is cheaper to calculate R^-1 in another way (D may
!            be a band matrix) then with the general routine, the user
!            may pass R^{-1}.  Indicated by ierr not equal to zero.
!  dvec   nx1 vector, the vector d from above (dp)
!         *** WILL BE DESTROYED ON EXIT ***
!         contains on exit the solution to the initial, i.e.,
!         unconstrained problem
!  fddmat scalar, the leading dimension of the matrix dmat
!  n      the dimension of dmat and dvec (int)
!  amat   nxq matrix, the matrix A from above (dp) [ A=(A1 A2)^T ]
!         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
!             CHANGED SIGNES ON EXIT ***
!  bvec   qx1 vector, the vector of constants b in the constraints (dp)
!         [ b = (b1^T b2^T)^T ]
!         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
!             CHANGED SIGNES ON EXIT ***
!  fdamat the first dimension of amat as declared in the calling program. 
!         fdamat >= n !!
!  q      integer, the number of constraints.
!  meq    integer, the number of equality constraints, 0 <= meq <= q.
!  ierr   integer, code for the status of the matrix D:
!            ierr =  0, we have to decompose D
!            ierr != 0, D is already decomposed into D=R^TR and we were
!                       given R^{-1}.
!
!  Output parameter:
!  sol   nx1 the final solution (x in the notation above)
!  lagr  qx1 the final Lagrange multipliers
!  crval scalar, the value of the criterion at the minimum      
!  iact  qx1 vector, the constraints which are active in the final
!        fit (int)
!  nact  scalar, the number of constraints active in the final fit (int)
!  iter  2x1 vector, first component gives the number of "main" 
!        iterations, the second one says how many constraints were
!        deleted after they became active
!  ierr  integer, error code on exit, if
!           ierr = 0, no problems
!           ierr = 1, the minimization problem has no solution
!           ierr = 2, problems with decomposing D, in this case sol
!                     contains garbage!!
!
!  Working space:
!  work  vector with length at least 2*n+r*(r+5)/2 + 2*q +1
!        where r=min(n,q)
!
subroutine qpgen2(dmat, dvec, fddmat, n, sol, lagr, crval, amat, bvec, fdamat, q, meq, iact, nact, iter, work, ierr)
implicit none
integer n, i, j, l, l1, info, q, iact(*), iter(*), it1, ierr, nact, iwzv, iwrv, iwrm, iwsv, iwuv, nvl, r, fdamat, iwnbv, meq, fddmat
real(dp) :: dmat(fddmat, *), dvec(*), lagr(*), sol(*), bvec(*)
real(dp) :: work(*), temp, sum, t1, tt, gc, gs, crval, nu
real(dp) :: amat(fdamat, *), vsmall, tmpa, tmpb
logical t1inf, t2min
it1 = 0
r = min(n,q)
l = 2*n + (r*(r+5))/2 + 2*q + 1
!
!     code gleaned from Powell's ZQPCVX routine to determine a small
!     number  that can be assumed to be an upper bound on the relative
!     precision of the computer arithmetic.
!
vsmall = 1.0d-60
1 vsmall = vsmall + vsmall
tmpa = 1.0d0 + 0.1d0*vsmall
tmpb = 1.0d0 + 0.2d0*vsmall
if( tmpa .LE. 1.0d0 ) goto 1
if( tmpb .LE. 1.0d0 ) goto 1
! 
! store the initial dvec to calculate below the unconstrained minima of
! the critical value.
!
loop_10: do i=1,n
work(i) = dvec(i)
end do loop_10
loop_11: do i=n+1,l
work(i) = 0.d0
end do loop_11
loop_12: do i=1,q
iact(i) = 0
lagr(i) = 0.d0
end do loop_12
!
! get the initial solution
!
if( ierr .EQ. 0 )then
call dpofa(dmat,fddmat,n,info)
if( info .NE. 0 )then
ierr = 2
goto 999
endif
call dposl(dmat,fddmat,n,dvec)
call dpori(dmat,fddmat,n)
else
!        
! Matrix D is already factorized, so we have to multiply d first with 
! R^-T and then with R^-1.  R^-1 is stored in the upper half of the
! array dmat.
!
loop_20: do j=1,n
sol(j)  = 0.d0
loop_21: do i=1,j
sol(j) = sol(j) + dmat(i,j)*dvec(i)
end do loop_21
end do loop_20
loop_22: do j=1,n
dvec(j) = 0.d0
loop_23: do i=j,n
dvec(j) = dvec(j) + dmat(j,i)*sol(i)
end do loop_23
end do loop_22
endif
!
! set lower triangular of dmat to zero, store dvec in sol and
! calculate value of the criterion at unconstrained minima
!
crval = 0.d0
loop_30: do j=1,n
sol(j)  = dvec(j)
crval   = crval + work(j)*sol(j)
work(j) = 0.d0
loop_32: do i=j+1,n
dmat(i,j) = 0.d0
end do loop_32
end do loop_30
crval = -crval/2.d0
ierr = 0
!
! calculate some constants, i.e., from which index on the different
! quantities are stored in the work matrix
!
iwzv  = n
iwrv  = iwzv + n
iwuv  = iwrv + r
iwrm  = iwuv + r+1
iwsv  = iwrm + (r*(r+1))/2
iwnbv = iwsv + q
!
! calculate the norm of each column of the A matrix
!
loop_51: do i=1,q
sum = 0.d0
loop_52: do j=1,n
sum = sum + amat(j,i)*amat(j,i)
end do loop_52
work(iwnbv+i) = sqrt(sum)
end do loop_51
nact = 0
iter(1) = 0
iter(2) = 0
50 continue
!
! start a new iteration      
!
iter(1) = iter(1)+1
!
! calculate all constraints and check which are still violated
! for the equality constraints we have to check whether the normal
! vector has to be negated (as well as bvec in that case)
!
l = iwsv
loop_60: do i=1,q
l = l+1
sum = -bvec(i)
loop_61: do j = 1,n
sum = sum + amat(j,i)*sol(j)
end do loop_61
if ( abs(sum) .LT. vsmall ) then
sum = 0.0d0
endif
if (i .GT. meq) then
work(l) = sum
else
work(l) = -abs(sum)
if (sum .GT. 0.d0) then
loop_62: do j=1,n
amat(j,i) = -amat(j,i)
end do loop_62
bvec(i) = -bvec(i)
endif
endif
end do loop_60
!
! as safeguard against rounding errors set already active constraints
! explicitly to zero
!
loop_70: do i=1,nact
work(iwsv+iact(i)) = 0.d0
end do loop_70
! 
! we weight each violation by the number of non-zero elements in the
! corresponding row of A. then we choose the violated constraint which
! has maximal absolute value, i.e., the minimum.
! by obvious commenting and uncommenting we can choose the strategy to      
! take always the first constraint which is violated. ;-)
!
nvl = 0
temp = 0.d0
loop_71: do i=1,q
if (work(iwsv+i) .LT. temp*work(iwnbv+i)) then
nvl = i
temp = work(iwsv+i)/work(iwnbv+i)
endif
!         if (work(iwsv+i) .LT. 0.d0) then
!            nvl = i
!            goto 72
!         endif
end do loop_71
if (nvl .EQ. 0) then
loop_73: do i=1,nact
lagr(iact(i))=work(iwuv+i)
end do loop_73
goto 999
endif
!     
! calculate d=J^Tn^+ where n^+ is the normal vector of the violated
! constraint. J is stored in dmat in this implementation!!
! if we drop a constraint, we have to jump back here.
!
55 continue
loop_80: do i=1,n
sum = 0.d0
loop_81: do j=1,n
sum = sum + dmat(j,i)*amat(j,nvl)
end do loop_81
work(i) = sum
end do loop_80
!
! Now calculate z = J_2 d_2
!
l1 = iwzv
loop_90: do i=1,n
work(l1+i) =0.d0
end do loop_90
loop_92: do j=nact+1,n
loop_93: do i=1,n
work(l1+i) = work(l1+i) + dmat(i,j)*work(j)
end do loop_93
end do loop_92
!
! and r = R^{-1} d_1, check also if r has positive elements (among the 
! entries corresponding to inequalities constraints).
!
t1inf = .TRUE.
loop_95: do i=nact,1,-1
sum = work(i)
l  = iwrm+(i*(i+3))/2
l1 = l-i
loop_96: do j=i+1,nact
sum = sum - work(l)*work(iwrv+j)
l   = l+j
end do loop_96
sum = sum / work(l1)
work(iwrv+i) = sum
if (iact(i) .LE. meq) cycle loop_95
if (sum .LE. 0.d0) cycle loop_95
t1inf = .FALSE.
it1 = i
end do loop_95
!
! if r has positive elements, find the partial step length t1, which is
! the maximum step in dual space without violating dual feasibility.
! it1  stores in which component t1, the min of u/r, occurs.
! 
if ( .NOT. t1inf) then
t1   = work(iwuv+it1)/work(iwrv+it1)
loop_100: do i=1,nact
if (iact(i) .LE. meq) cycle loop_100
if (work(iwrv+i) .LE. 0.d0) cycle loop_100
temp = work(iwuv+i)/work(iwrv+i)
if (temp .LT. t1) then
t1   = temp
it1  = i
endif
end do loop_100
endif
!
! test if the z vector is equal to zero
!
sum = 0.d0
loop_110: do i=iwzv+1,iwzv+n
sum = sum + work(i)*work(i)
end do loop_110
if (abs(sum) .LE. vsmall) then
!         
! No step in primal space such that the new constraint becomes
! feasible. Take step in dual space and drop a constant.
!     
if (t1inf) then
!            
! No step in dual space possible either, problem is not solvable
!
ierr = 1
goto 999
else
!
! we take a partial step in dual space and drop constraint it1,
! that is, we drop the it1-th active constraint.
! then we continue at step 2(a) (marked by label 55)
!
loop_111: do i=1,nact
work(iwuv+i) = work(iwuv+i) - t1*work(iwrv+i)
end do loop_111
work(iwuv+nact+1) = work(iwuv+nact+1) + t1
goto 700
endif
else
!
! compute full step length t2, minimum step in primal space such that
! the constraint becomes feasible.
! keep sum (which is z^Tn^+) to update crval below!
!
sum = 0.d0
loop_120: do i = 1,n
sum = sum + work(iwzv+i)*amat(i,nvl)
end do loop_120
tt = -work(iwsv+nvl)/sum
t2min = .TRUE.
if (.NOT. t1inf) then
if (t1 .LT. tt) then
tt    = t1
t2min = .FALSE.
endif
endif
!
! take step in primal and dual space
!
loop_130: do i=1,n
sol(i) = sol(i) + tt*work(iwzv+i)
end do loop_130
crval = crval + tt*sum*(tt/2.d0 + work(iwuv+nact+1))
loop_131: do i=1,nact
work(iwuv+i) = work(iwuv+i) - tt*work(iwrv+i)
end do loop_131
work(iwuv+nact+1) = work(iwuv+nact+1) + tt
!
! if it was a full step, then we check wheter further constraints are
! violated otherwise we can drop the current constraint and iterate once
! more 
if(t2min) then
!
! we took a full step. Thus add constraint nvl to the list of active
! constraints and update J and R
!
nact = nact + 1
iact(nact) = nvl
!
! to update R we have to put the first nact-1 components of the d vector
! into column (nact) of R
!
l = iwrm + ((nact-1)*nact)/2 + 1
loop_150: do i=1,nact-1
work(l) = work(i)
l = l+1
end do loop_150
!
! if now nact=n, then we just have to add the last element to the new
! row of R.
! Otherwise we use Givens transformations to turn the vector d(nact:n)
! into a multiple of the first unit vector. That multiple goes into the
! last element of the new row of R and J is accordingly updated by the
! Givens transformations. 
!
if (nact .EQ. n) then
work(l) = work(n)
else
loop_160: do i=n,nact+1,-1
!
! we have to find the Givens rotation which will reduce the element
! (l1) of d to zero.
! if it is already zero we don't have to do anything, except of
! decreasing l1
!
if (work(i) .EQ. 0.d0) cycle loop_160
gc   = max(abs(work(i-1)),abs(work(i)))
gs   = min(abs(work(i-1)),abs(work(i)))
temp = sign(gc*sqrt(1+(gs/gc)*(gs/gc)), work(i-1))
gc   = work(i-1)/temp
gs   = work(i)/temp
! 
! The Givens rotation is done with the matrix (gc gs, gs -gc).
! If gc is one, then element (i) of d is zero compared with element
! (l1-1). Hence we don't have to do anything. 
! If gc is zero, then we just have to switch column (i) and column (i-1) 
! of J. Since we only switch columns in J, we have to be careful how we
! update d depending on the sign of gs.
! Otherwise we have to apply the Givens rotation to these columns.
! The i-1 element of d has to be updated to temp.                  
!
if (gc .EQ. 1.d0) cycle loop_160
if (gc .EQ. 0.d0) then
work(i-1) = gs * temp
loop_170: do j=1,n
temp        = dmat(j,i-1)
dmat(j,i-1) = dmat(j,i)
dmat(j,i)   = temp
end do loop_170
else
work(i-1) = temp
nu = gs/(1.d0+gc)
loop_180: do j=1,n
temp        = gc*dmat(j,i-1) + gs*dmat(j,i)
dmat(j,i)   = nu*(dmat(j,i-1)+temp) - dmat(j,i)
dmat(j,i-1) = temp
end do loop_180
endif
end do loop_160
!     
! l is still pointing to element (nact,nact) of the matrix R.
! So store d(nact) in R(nact,nact)
work(l) = work(nact)
endif
else
!
! we took a partial step in dual space. Thus drop constraint it1,
! that is, we drop the it1-th active constraint.
! then we continue at step 2(a) (marked by label 55)
! but since the fit changed, we have to recalculate now "how much"
! the fit violates the chosen constraint now.
!
sum = -bvec(nvl)
loop_190: do j = 1,n
sum = sum + sol(j)*amat(j,nvl)
end do loop_190
if( nvl .GT. meq ) then
work(iwsv+nvl) = sum
else
work(iwsv+nvl) = -abs(sum)
if( sum .GT. 0.d0) then
loop_191: do j=1,n
amat(j,nvl) = -amat(j,nvl)
end do loop_191
bvec(nvl) = -bvec(nvl)
endif
endif
goto 700
endif
endif
goto 50
!
! Drop constraint it1
!
700 continue
!
! if it1 = nact it is only necessary to update the vector u and nact
!
if (it1 .EQ. nact) goto 799
!
! After updating one row of R (column of J) we will also come back here
!
797 continue
!
! we have to find the Givens rotation which will reduce the element
! (it1+1,it1+1) of R to zero.
! if it is already zero we don't have to do anything except of updating
! u, iact, and shifting column (it1+1) of R to column (it1)
! l  will point to element (1,it1+1) of R
! l1 will point to element (it1+1,it1+1) of R
!
l  = iwrm + (it1*(it1+1))/2 + 1
l1 = l+it1
if (work(l1) .EQ. 0.d0) goto 798
gc   = max(abs(work(l1-1)),abs(work(l1)))
gs   = min(abs(work(l1-1)),abs(work(l1)))
temp = sign(gc*sqrt(1+(gs/gc)*(gs/gc)), work(l1-1))
gc   = work(l1-1)/temp
gs   = work(l1)/temp
! 
! The Givens rotatin is done with the matrix (gc gs, gs -gc).
! If gc is one, then element (it1+1,it1+1) of R is zero compared with
! element (it1,it1+1). Hence we don't have to do anything.
! if gc is zero, then we just have to switch row (it1) and row (it1+1)
! of R and column (it1) and column (it1+1) of J. Since we swithc rows in
! R and columns in J, we can ignore the sign of gs.
! Otherwise we have to apply the Givens rotation to these rows/columns.
!
if (gc .EQ. 1.d0) goto 798
if (gc .EQ. 0.d0) then
loop_710: do i=it1+1,nact
temp       = work(l1-1)
work(l1-1) = work(l1)
work(l1)   = temp
l1 = l1+i
end do loop_710
loop_711: do i=1,n
temp          = dmat(i,it1)
dmat(i,it1)   = dmat(i,it1+1)
dmat(i,it1+1) = temp
end do loop_711
else
nu = gs/(1.d0+gc)
loop_720: do i=it1+1,nact
temp       = gc*work(l1-1) + gs*work(l1)
work(l1)   = nu*(work(l1-1)+temp) - work(l1)
work(l1-1) = temp
l1 = l1+i
end do loop_720
loop_721: do i=1,n
temp          = gc*dmat(i,it1) + gs*dmat(i,it1+1)
dmat(i,it1+1) = nu*(dmat(i,it1)+temp) - dmat(i,it1+1)
dmat(i,it1)   = temp
end do loop_721
endif
!
! shift column (it1+1) of R to column (it1) (that is, the first it1
! elements). The posit1on of element (1,it1+1) of R was calculated above
! and stored in l.
!
798 continue
l1 = l-it1
loop_730: do i=1,it1
work(l1)=work(l)
l  = l+1
l1 = l1+1
end do loop_730
!      
! update vector u and iact as necessary
! Continue with updating the matrices J and R
!
work(iwuv+it1) = work(iwuv+it1+1)
iact(it1)      = iact(it1+1)
it1 = it1+1
if (it1 .LT. nact) goto 797
799 work(iwuv+nact)   = work(iwuv+nact+1)
work(iwuv+nact+1) = 0.d0
iact(nact)        = 0
nact = nact-1
iter(2) = iter(2)+1
goto 55
999 continue
return
end subroutine qpgen2


  subroutine dpofa(a, lda, n, info)
    integer, intent(in) :: lda, n
    integer, intent(out) :: info
    real(dp), intent(inout) :: a(lda, *)
    integer :: i, j, k
    real(dp) :: s

    info = 0
    do j = 1, n
      do k = 1, j - 1
        s = a(k, j)
        do i = 1, k - 1
          s = s - a(i, k) * a(i, j)
        end do
        s = s / a(k, k)
        a(k, j) = s
      end do
      s = a(j, j)
      do i = 1, j - 1
        s = s - a(i, j) * a(i, j)
      end do
      if (s <= 0.0d0) then
        info = j
        return
      end if
      a(j, j) = sqrt(s)
    end do
  end subroutine dpofa

  subroutine dposl(a, lda, n, b)
    integer, intent(in) :: lda, n
    real(dp), intent(in) :: a(lda, *)
    real(dp), intent(inout) :: b(*)
    integer :: i, j
    real(dp) :: s

    do i = 1, n
      s = b(i)
      do j = 1, i - 1
        s = s - a(j, i) * b(j)
      end do
      b(i) = s / a(i, i)
    end do
    do i = n, 1, -1
      s = b(i)
      do j = i + 1, n
        s = s - a(i, j) * b(j)
      end do
      b(i) = s / a(i, i)
    end do
  end subroutine dposl

  subroutine dpori(a, lda, n)
    integer, intent(in) :: lda, n
    real(dp), intent(inout) :: a(lda, *)
    integer :: i, j
    real(dp), allocatable :: rinv(:, :)

    allocate(rinv(n, n))
    rinv = 0.0d0
    do j = 1, n
      rinv(j, j) = 1.0d0 / a(j, j)
      do i = j - 1, 1, -1
        rinv(i, j) = -sum(a(i, i + 1:j) * rinv(i + 1:j, j)) / a(i, i)
      end do
    end do
    do j = 1, n
      do i = 1, j
        a(i, j) = rinv(i, j)
      end do
    end do
  end subroutine dpori
!
!  Copyright (C) 1995-2010 Berwin A. Turlach <Berwin.Turlach@gmail.com>
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
!  USA.
!
!  this routine uses the Goldfarb/Idnani algorithm to solve the
!  following minimization problem:
!
!        minimize  -d^T x + 1/2 *  x^T D x
!        where   A1^T x  = b1
!                A2^T x >= b2
!
!  the matrix D is assumed to be positive definite.  Especially,
!  w.l.o.g. D is assumed to be symmetric.
!  
!  Input parameter:
!  dmat   nxn matrix, the matrix D from above (dp)
!         *** WILL BE DESTROYED ON EXIT ***
!         The user has two possibilities:
!         a) Give D (ierr=0), in this case we use routines from LINPACK
!            to decompose D.
!         b) To get the algorithm started we need R^-1, where D=R^TR.
!            So if it is cheaper to calculate R^-1 in another way (D may
!            be a band matrix) then with the general routine, the user
!            may pass R^{-1}.  Indicated by ierr not equal to zero.
!  dvec   nx1 vector, the vector d from above (dp)
!         *** WILL BE DESTROYED ON EXIT ***
!         contains on exit the solution to the initial, i.e.,
!         unconstrained problem
!  fddmat scalar, the leading dimension of the matrix dmat
!  n      the dimension of dmat and dvec (int)
!  amat   lxq matrix (dp)
!         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
!             CHANGED SIGNES ON EXIT ***
!  iamat  (l+1)xq matrix (int)
!         these two matrices store the matrix A in compact form. the format
!         is: [ A=(A1 A2)^T ]
!           iamat(1,i) is the number of non-zero elements in column i of A
!           iamat(k,i) for k>=2, is equal to j if the (k-1)-th non-zero 
!                      element in column i of A is A(i,j)
!            amat(k,i) for k>=1, is equal to the k-th non-zero element
!                      in column i of A.
!           
!  bvec   qx1 vector, the vector of constants b in the constraints (dp)
!         [ b = (b1^T b2^T)^T ]
!         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
!             CHANGED SIGNES ON EXIT ***
!  fdamat the first dimension of amat as declared in the calling program. 
!         fdamat >= n (and iamat must have fdamat+1 as first dimension)
!  q      integer, the number of constraints.
!  meq    integer, the number of equality constraints, 0 <= meq <= q.
!  ierr   integer, code for the status of the matrix D:
!            ierr =  0, we have to decompose D
!            ierr != 0, D is already decomposed into D=R^TR and we were
!                       given R^{-1}.
!
!  Output parameter:
!  sol   nx1 the final solution (x in the notation above)
!  lagr  qx1 the final Lagrange multipliers
!  crval scalar, the value of the criterion at the minimum      
!  iact  qx1 vector, the constraints which are active in the final
!        fit (int)
!  nact  scalar, the number of constraints active in the final fit (int)
!  iter  2x1 vector, first component gives the number of "main" 
!        iterations, the second one says how many constraints were
!        deleted after they became active
!  ierr  integer, error code on exit, if
!           ierr = 0, no problems
!           ierr = 1, the minimization problem has no solution
!           ierr = 2, problems with decomposing D, in this case sol
!                     contains garbage!!
!
!  Working space:
!  work  vector with length at least 2*n+r*(r+5)/2 + 2*q +1
!        where r=min(n,q)
!
subroutine qpgen1(dmat, dvec, fddmat, n, sol, lagr, crval, amat, &
  iamat, bvec, fdamat, q, meq, iact, nact, iter, work, ierr)
implicit none
integer :: n, i, j, l, l1, info, q, fdamat, it1, ierr, nact
integer :: iwzv, iwrv, iwrm, iwsv, iwuv, nvl, r, iwnbv, meq, fddmat
integer :: iamat(fdamat + 1, *), iact(*), iter(*)
real(dp) :: dmat(fddmat, *), dvec(*), sol(*), lagr(*), bvec(*)
real(dp) :: work(*), temp, sum, t1, tt, gc, gs, crval, nu
real(dp) :: amat(fdamat, *), vsmall, tmpa, tmpb
logical :: t1inf, t2min
it1 = 0
r = min(n,q)
l = 2*n + (r*(r+5))/2 + 2*q + 1
!
!     code gleaned from Powell's ZQPCVX routine to determine a small
!     number  that can be assumed to be an upper bound on the relative
!     precision of the computer arithmetic.
!
vsmall = 1.0d-60
1 vsmall = vsmall + vsmall
tmpa = 1.0d0 + 0.1d0*vsmall
tmpb = 1.0d0 + 0.2d0*vsmall
if( tmpa .LE. 1.0d0 ) goto 1
if( tmpb .LE. 1.0d0 ) goto 1
! 
! store the initial dvec to calculate below the unconstrained minima of
! the critical value.
!
loop_10: do i=1,n
work(i) = dvec(i)
end do loop_10
loop_11: do i=n+1,l
work(i) = 0.d0
end do loop_11
loop_12: do i=1,q
iact(i) = 0
lagr(i) = 0.d0
end do loop_12
!
! get the initial solution
!
if( ierr .EQ. 0 )then
call dpofa(dmat,fddmat,n,info)
if( info .NE. 0 )then
ierr = 2
goto 999
endif
call dposl(dmat,fddmat,n,dvec)
call dpori(dmat,fddmat,n)
else
!        
! Matrix D is already factorized, so we have to multiply d first with 
! R^-T and then with R^-1.  R^-1 is stored in the upper half of the
! array dmat.
!
loop_20: do j=1,n
sol(j)  = 0.d0
loop_21: do i=1,j
sol(j) = sol(j) + dmat(i,j)*dvec(i)
end do loop_21
end do loop_20
loop_22: do j=1,n
dvec(j) = 0.d0
loop_23: do i=j,n
dvec(j) = dvec(j) + dmat(j,i)*sol(i)
end do loop_23
end do loop_22
endif
!
! set lower triangular of dmat to zero, store dvec in sol and
! calculate value of the criterion at unconstrained minima
!
crval = 0.d0
loop_30: do j=1,n
sol(j)  = dvec(j)
crval   = crval + work(j)*sol(j)
work(j) = 0.d0
loop_32: do i=j+1,n
dmat(i,j) = 0.d0
end do loop_32
end do loop_30
crval = -crval/2.d0
ierr  = 0
!
! calculate some constants, i.e., from which index on the different
! quantities are stored in the work matrix
!
iwzv  = n
iwrv  = iwzv + n
iwuv  = iwrv + r
iwrm  = iwuv + r+1
iwsv  = iwrm + (r*(r+1))/2
iwnbv = iwsv + q
!
! calculate the norm of each column of the A matrix
!
loop_51: do i=1,q
sum = 0.d0
loop_52: do j=1,iamat(1,i)
sum = sum + amat(j,i)*amat(j,i)
end do loop_52
work(iwnbv+i) = sqrt(sum)
end do loop_51
nact = 0
iter(1) = 0
iter(2) = 0
50 continue
!
! start a new iteration      
!
iter(1) = iter(1)+1
!
! calculate all constraints and check which are still violated
! for the equality constraints we have to check whether the normal
! vector has to be negated (as well as bvec in that case)
!
l = iwsv
loop_60: do i=1,q
l = l+1
sum = -bvec(i)
loop_61: do j = 1,iamat(1,i)
sum = sum + amat(j,i)*sol(iamat(j+1,i))
end do loop_61
if ( abs(sum) .LT. vsmall ) then
sum = 0.0d0
endif
if (i .GT. meq) then
work(l) = sum
else
work(l) = -abs(sum)
if (sum .GT. 0.d0) then
loop_62: do j=1,iamat(1,i)
amat(j,i) = -amat(j,i)
end do loop_62
bvec(i) = -bvec(i)
endif
endif
end do loop_60
!
! as safeguard against rounding errors set already active constraints
! explicitly to zero
!
loop_70: do i=1,nact
work(iwsv+iact(i)) = 0.d0
end do loop_70
! 
! we weight each violation by the number of non-zero elements in the
! corresponding row of A. then we choose the violated constraint which
! has maximal absolute value, i.e., the minimum.
! by obvious commenting and uncommenting we can choose the strategy to      
! take always the first constraint which is violated. ;-)
!
nvl = 0
temp = 0.d0
loop_71: do i=1,q
if (work(iwsv+i) .LT. temp*work(iwnbv+i)) then
nvl = i
temp = work(iwsv+i)/work(iwnbv+i)
endif
!         if (work(iwsv+i) .LT. 0.d0) then
!            nvl = i
!            goto 72
!         endif
end do loop_71
if (nvl .EQ. 0) then
loop_73: do i=1,nact
lagr(iact(i))=work(iwuv+i)
end do loop_73
goto 999
endif
!     
! calculate d=J^Tn^+ where n^+ is the normal vector of the violated
! constraint. J is stored in dmat in this implementation!!
! if we drop a constraint, we have to jump back here.
!
55 continue
loop_80: do i=1,n
sum = 0.d0
loop_81: do j=1,iamat(1,nvl)
sum = sum + dmat(iamat(j+1,nvl),i)*amat(j,nvl)
end do loop_81
work(i) = sum
end do loop_80
!
! Now calculate z = J_2 d_2
!
l1 = iwzv
loop_90: do i=1,n
work(l1+i) =0.d0
end do loop_90
loop_92: do j=nact+1,n
loop_93: do i=1,n
work(l1+i) = work(l1+i) + dmat(i,j)*work(j)
end do loop_93
end do loop_92
!
! and r = R^{-1} d_1, check also if r has positive elements (among the 
! entries corresponding to inequalities constraints).
!
t1inf = .TRUE.
loop_95: do i=nact,1,-1
sum = work(i)
l  = iwrm+(i*(i+3))/2
l1 = l-i
loop_96: do j=i+1,nact
sum = sum - work(l)*work(iwrv+j)
l   = l+j
end do loop_96
sum = sum / work(l1)
work(iwrv+i) = sum
if (iact(i) .LE. meq) cycle loop_95
if (sum .LE. 0.d0) cycle loop_95
t1inf = .FALSE.
it1 = i
end do loop_95
!
! if r has positive elements, find the partial step length t1, which is
! the maximum step in dual space without violating dual feasibility.
! it1  stores in which component t1, the min of u/r, occurs.
! 
if ( .NOT. t1inf) then
t1   = work(iwuv+it1)/work(iwrv+it1)
loop_100: do i=1,nact
if (iact(i) .LE. meq) cycle loop_100
if (work(iwrv+i) .LE. 0.d0) cycle loop_100
temp = work(iwuv+i)/work(iwrv+i)
if (temp .LT. t1) then
t1   = temp
it1  = i
endif
end do loop_100
endif
!
! test if the z vector is equal to zero
!
sum = 0.d0
loop_110: do i=iwzv+1,iwzv+n
sum = sum + work(i)*work(i)
end do loop_110
if (abs(sum) .LE. vsmall) then
!         
! No step in primal space such that the new constraint becomes
! feasible. Take step in dual space and drop a constant.
!     
if (t1inf) then
!            
! No step in dual space possible either, problem is not solvable
!
ierr = 1
goto 999
else
!
! we take a partial step in dual space and drop constraint it1,
! that is, we drop the it1-th active constraint.
! then we continue at step 2(a) (marked by label 55)
!
loop_111: do i=1,nact
work(iwuv+i) = work(iwuv+i) - t1*work(iwrv+i)
end do loop_111
work(iwuv+nact+1) = work(iwuv+nact+1) + t1
goto 700
endif
else
!
! compute full step length t2, minimum step in primal space such that
! the constraint becomes feasible.
! keep sum (which is z^Tn^+) to update crval below!
!
sum = 0.d0
loop_120: do i = 1,iamat(1,nvl)
sum = sum + work(iwzv+iamat(i+1,nvl))*amat(i,nvl)
end do loop_120
tt = -work(iwsv+nvl)/sum
t2min = .TRUE.
if (.NOT. t1inf) then
if (t1 .LT. tt) then
tt    = t1
t2min = .FALSE.
endif
endif
!
! take step in primal and dual space
!
loop_130: do i=1,n
sol(i) = sol(i) + tt*work(iwzv+i)
end do loop_130
crval = crval + tt*sum*(tt/2.d0 + work(iwuv+nact+1))
loop_131: do i=1,nact
work(iwuv+i) = work(iwuv+i) - tt*work(iwrv+i)
end do loop_131
work(iwuv+nact+1) = work(iwuv+nact+1) + tt
!
! if it was a full step, then we check wheter further constraints are
! violated otherwise we can drop the current constraint and iterate once
! more 
if(t2min) then
!
! we took a full step. Thus add constraint nvl to the list of active
! constraints and update J and R
!
nact = nact + 1
iact(nact) = nvl
!
! to update R we have to put the first nact-1 components of the d vector
! into column (nact) of R
!
l = iwrm + ((nact-1)*nact)/2 + 1
loop_150: do i=1,nact-1
work(l) = work(i)
l = l+1
end do loop_150
!
! if now nact=n, then we just have to add the last element to the new
! row of R.
! Otherwise we use Givens transformations to turn the vector d(nact:n)
! into a multiple of the first unit vector. That multiple goes into the
! last element of the new row of R and J is accordingly updated by the
! Givens transformations. 
!
if (nact .EQ. n) then
work(l) = work(n)
else
loop_160: do i=n,nact+1,-1
!
! we have to find the Givens rotation which will reduce the element
! (l1) of d to zero.
! if it is already zero we don't have to do anything, except of
! decreasing l1
!
if (work(i) .EQ. 0.d0) cycle loop_160
gc   = max(abs(work(i-1)),abs(work(i)))
gs   = min(abs(work(i-1)),abs(work(i)))
temp = sign(gc*sqrt(1+(gs/gc)*(gs/gc)), work(i-1))
gc   = work(i-1)/temp
gs   = work(i)/temp
! 
! The Givens rotation is done with the matrix (gc gs, gs -gc).
! If gc is one, then element (i) of d is zero compared with element
! (l1-1). Hence we don't have to do anything. 
! If gc is zero, then we just have to switch column (i) and column (i-1) 
! of J. Since we only switch columns in J, we have to be careful how we
! update d depending on the sign of gs.
! Otherwise we have to apply the Givens rotation to these columns.
! The i-1 element of d has to be updated to temp.                  
!
if (gc .EQ. 1.d0) cycle loop_160
if (gc .EQ. 0.d0) then
work(i-1) = gs * temp
loop_170: do j=1,n
temp        = dmat(j,i-1)
dmat(j,i-1) = dmat(j,i)
dmat(j,i)   = temp
end do loop_170
else
work(i-1) = temp
nu = gs/(1.d0+gc)
loop_180: do j=1,n
temp        = gc*dmat(j,i-1) + gs*dmat(j,i)
dmat(j,i)   = nu*(dmat(j,i-1)+temp) - dmat(j,i)
dmat(j,i-1) = temp
end do loop_180
endif
end do loop_160
!     
! l is still pointing to element (nact,nact) of the matrix R.
! So store d(nact) in R(nact,nact)
work(l) = work(nact)
endif
else
!
! we took a partial step in dual space. Thus drop constraint it1,
! that is, we drop the it1-th active constraint.
! then we continue at step 2(a) (marked by label 55)
! but since the fit changed, we have to recalculate now "how much"
! the fit violates the chosen constraint now.
!
sum = -bvec(nvl)
loop_190: do j = 1,iamat(1,nvl)
sum = sum + sol(iamat(j+1,nvl))*amat(j,nvl)
end do loop_190
if( nvl .GT. meq ) then
work(iwsv+nvl) = sum
else
work(iwsv+nvl) = -abs(sum)
if( sum .GT. 0.d0) then
loop_191: do j=1,iamat(1,nvl)
amat(j,nvl) = -amat(j,nvl)
end do loop_191
bvec(nvl) = -bvec(nvl)
endif
endif
goto 700
endif
endif
goto 50
!
! Drop constraint it1
!
700 continue
!
! if it1 = nact it is only necessary to update the vector u and nact
!
if (it1 .EQ. nact) goto 799
!
! After updating one row of R (column of J) we will also come back here
!
797 continue
!
! we have to find the Givens rotation which will reduce the element
! (it1+1,it1+1) of R to zero.
! if it is already zero we don't have to do anything except of updating
! u, iact, and shifting column (it1+1) of R to column (it1)
! l  will point to element (1,it1+1) of R
! l1 will point to element (it1+1,it1+1) of R
!
l  = iwrm + (it1*(it1+1))/2 + 1
l1 = l+it1
if (work(l1) .EQ. 0.d0) goto 798
gc   = max(abs(work(l1-1)),abs(work(l1)))
gs   = min(abs(work(l1-1)),abs(work(l1)))
temp = sign(gc*sqrt(1+(gs/gc)*(gs/gc)), work(l1-1))
gc   = work(l1-1)/temp
gs   = work(l1)/temp
! 
! The Givens rotatin is done with the matrix (gc gs, gs -gc).
! If gc is one, then element (it1+1,it1+1) of R is zero compared with
! element (it1,it1+1). Hence we don't have to do anything.
! if gc is zero, then we just have to switch row (it1) and row (it1+1)
! of R and column (it1) and column (it1+1) of J. Since we swithc rows in
! R and columns in J, we can ignore the sign of gs.
! Otherwise we have to apply the Givens rotation to these rows/columns.
!
if (gc .EQ. 1.d0) goto 798
if (gc .EQ. 0.d0) then
loop_710: do i=it1+1,nact
temp       = work(l1-1)
work(l1-1) = work(l1)
work(l1)   = temp
l1 = l1+i
end do loop_710
loop_711: do i=1,n
temp          = dmat(i,it1)
dmat(i,it1)   = dmat(i,it1+1)
dmat(i,it1+1) = temp
end do loop_711
else
nu = gs/(1.d0+gc)
loop_720: do i=it1+1,nact
temp       = gc*work(l1-1) + gs*work(l1)
work(l1)   = nu*(work(l1-1)+temp) - work(l1)
work(l1-1) = temp
l1 = l1+i
end do loop_720
loop_721: do i=1,n
temp          = gc*dmat(i,it1) + gs*dmat(i,it1+1)
dmat(i,it1+1) = nu*(dmat(i,it1)+temp) - dmat(i,it1+1)
dmat(i,it1)   = temp
end do loop_721
endif
!
! shift column (it1+1) of R to column (it1) (that is, the first it1
! elements). The posit1on of element (1,it1+1) of R was calculated above
! and stored in l.
!
798 continue
l1 = l-it1
loop_730: do i=1,it1
work(l1)=work(l)
l  = l+1
l1 = l1+1
end do loop_730
!      
! update vector u and iact as necessary
! Continue with updating the matrices J and R
!
work(iwuv+it1) = work(iwuv+it1+1)
iact(it1)      = iact(it1+1)
it1 = it1+1
if (it1 .LT. nact) goto 797
799 work(iwuv+nact)   = work(iwuv+nact+1)
work(iwuv+nact+1) = 0.d0
iact(nact)        = 0
nact = nact-1
iter(2) = iter(2)+1
goto 55
999 continue
return
end subroutine qpgen1

end module quadprog_core_mod
