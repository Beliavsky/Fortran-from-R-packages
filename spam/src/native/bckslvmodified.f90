!     Changes Oct 2019: to address gcc-10 issues with Rank missmatch:
!     in `call blksb[lf]`: b(1,j) -> b(:,j); newrhs -> newrhs(:)



subroutine backsolvef(m,nsuper,nrhs,lindx,xlindx,lnz, &
& xlnz,xsuper,b)
! see below...
implicit none

integer m,nsuper,nrhs,lindx(*),xlindx(m+1), &
& xlnz(m+1),xsuper(m+1)
double precision lnz(*),b(m,nrhs)
integer j
do j = 1,nrhs
call blkslb(nsuper,xsuper,xlindx,lindx,xlnz,lnz,b(:,j))
enddo
return
end


subroutine forwardsolvef(m,nsuper,nrhs,lindx,xlindx, &
& lnz,xlnz,xsuper,b)
! INPUT:
!     m -- the number of column in the matrix
!     lindx -- an nsub-vector of interger which contains, in
!           column major oder, the row subscripts of the nonzero
!           entries in L in a compressed storage format
!     xlindx -- an nsuper-vector of integer of pointers for lindx
!     lnz -- First contains the non-zero entries of d; later
!            contains the entries of the Cholesky factor
!     xlnz -- column pointer for L stored in lnz
!     xsuper -- array of length m+1 containing the supernode
!               partitioning
!     b -- the rhs of the equality constraint
! OUTPUT:
!     b -- the solution

implicit none

integer m,nsuper,nrhs,lindx(*),xlindx(m+1), &
& xlnz(m+1),xsuper(m+1)
double precision lnz(*),b(m,nrhs)
integer j
!
do j = 1,nrhs
call blkslf(nsuper,xsuper,xlindx,lindx,xlnz,lnz,b(:,j))
enddo
return
end
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!                   Slight modification by Reinhard Furrer
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!
!***********************************************************************
subroutine pivotforwardsolve(m,nsuper,nrhs,lindx,xlindx,lnz, &
& xlnz,invp,perm,xsuper,newrhs,sol,b)
! Sparse least squares solver via Ng-Peyton's sparse Cholesky
!    factorization for sparse symmetric positive definite
! INPUT:
!     m -- the number of column in the design matrix X
!     nsubmax -- upper bound of the dimension of lindx
!     lindx -- an nsub-vector of interger which contains, in
!           column major oder, the row subscripts of the nonzero
!           entries in L in a compressed storage format
!     xlindx -- an nsuper-vector of integer of pointers for lindx
!     lnz -- First contains the non-zero entries of d; later
!            contains the entries of the Cholesky factor
!     xlnz -- column pointer for L stored in lnz
!     invp -- an m-vector of integer of inverse permutation
!             vector
!     perm -- an m-vector of integer of permutation vector
!     xsuper -- array of length m+1 containing the supernode
!               partitioning
!     newrhs -- extra work vector for right-hand side and
!               solution
!     sol -- the least squares solution
!     b -- an m-vector, usualy the rhs of the equality constraint
!          X'a = (1-tau)X'e in the rq setting
! OUTPUT:
!     y -- an m-vector of least squares solution
! WORK ARRAYS:
!     b -- an m-vector, usually the rhs of the equality constraint
!          X'a = (1-tau)X'e in the rq setting
implicit none

integer m,nsuper,nrhs,lindx(*),xlindx(m+1), &
& invp(m),perm(m),xlnz(m+1), xsuper(m+1)
integer i,j
double precision lnz(*),b(m,nrhs),newrhs(m),sol(m,nrhs)

do j = 1,nrhs
do i = 1,m
newrhs(i) = b(perm(i),j)
enddo
call blkslf(nsuper,xsuper,xlindx,lindx,xlnz,lnz,newrhs(:))
do i = 1,m
sol(i,j) = newrhs(invp(i))
enddo
enddo
return
end
!***********************************************************************
subroutine pivotbacksolve(m,nsuper,nrhs,lindx,xlindx,lnz, &
& xlnz,invp,perm,xsuper,newrhs,sol,b)
!     see above
implicit none

integer m, nsuper,nrhs,lindx(*),xlindx(m+1), &
& invp(m),perm(m),xlnz(m+1), xsuper(m+1)
double precision lnz(*),b(m,nrhs),newrhs(m),sol(m,nrhs)
integer i,j
do j = 1,nrhs
do i = 1,m
newrhs(i) = b(perm(i),j)
enddo
call blkslb(nsuper,xsuper,xlindx,lindx,xlnz,lnz,newrhs)
do i = 1,m
sol(i,j) = newrhs(invp(i))
enddo
enddo
return
end
!***********************************************************************
subroutine backsolves(m,nsuper,nrhs,lindx,xlindx,lnz, &
& xlnz,invp,perm,xsuper,newrhs,sol,b)
! Sparse least squares solver via Ng-Peyton's sparse Cholesky
!    factorization for sparse symmetric positive definite
! INPUT:
!     m -- the number of column in the design matrix X
!     nsubmax -- upper bound of the dimension of lindx
!     lindx -- an nsub-vector of interger which contains, in
!           column major oder, the row subscripts of the nonzero
!           entries in L in a compressed storage format
!     xlindx -- an nsuper-vector of integer of pointers for lindx
!     nnzlmax -- the upper bound of the non-zero entries in
!                L stored in lnz, including the diagonal entries
!     lnz -- First contains the non-zero entries of d; later
!            contains the entries of the Cholesky factor
!     xlnz -- column pointer for L stored in lnz
!     invp -- an m-vector of integer of inverse permutation
!             vector
!     perm -- an m-vector of integer of permutation vector
!     xsuper -- array of length m+1 containing the supernode
!               partitioning
!     newrhs -- extra work vector for right-hand side and
!               solution
!     sol -- the least squares solution
!     b -- an m-vector, usualy the rhs of the equality constraint
!          X'a = (1-tau)X'e in the rq setting
! OUTPUT:
!     y -- an m-vector of least squares solution
! WORK ARRAYS:
!     b -- an m-vector, usually the rhs of the equality constraint
!          X'a = (1-tau)X'e in the rq setting
implicit none

integer m,nsuper,nrhs,lindx(*),xlindx(m+1), &
& invp(m),perm(m),xlnz(m+1), xsuper(m+1)
double precision lnz(*),b(m,nrhs),newrhs(m),sol(m,nrhs)

integer i,j
do j = 1,nrhs
do i = 1,m
newrhs(i) = b(perm(i),j)
enddo
call blkslv(nsuper,xsuper,xlindx,lindx,xlnz,lnz,newrhs(:))
do i = 1,m
sol(i,j) = newrhs(invp(i))
enddo
enddo
return
end
