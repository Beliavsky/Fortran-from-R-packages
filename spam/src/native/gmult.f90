!-----------------------------------------------------------------
subroutine gri(i,rpin,rindex)
implicit none
integer j, i, rpin(*), rindex
!-----------------------------------------------------------------
!     purpose:
!     --------
!     get row index for the i-th entrie of the sparse matrix
! on entry:
!----------
!     i      = number of the element
!     rpin   = input rowpointers vector
! on return:
!-----------
!     rindex  = output row index
!-----------------------------------------------------------------
j=1
do while (rpin(j) <= i)
j=j+1
end do
rindex=j-1
return
end
!--------end-of-gri---------------------------------------------

!-----------------------------------------------------------------
subroutine gfact(i,j,splits,fact,nfact,out)
implicit none
integer j, i, nfact, splits(nfact+1)
double precision fact(nfact,nfact), out
!-----------------------------------------------------------------
!     purpose:
!     --------
!     get fact for coordinates i,j,
! on entry:
!----------
!     i,j    = indices of element
!     splits = splits of the matrix
!     fact   = values to be returned
!     nfact  = ncol(fact)=nrow(fact)
! return:
!-----------
!     out    = the correct value of fact
!----------------------------------------------------------------
integer ii, jj
if (i >= splits(nfact+1) .OR. j >= splits(nfact+1)) then
goto 9000
!         stop 'i,j are out of larger than (nfact+1)'
end if
ii=1
do while (splits(ii+1) <= i)
ii=ii+1
end do
jj=1
do while (splits(jj+1) <= j)
jj=jj+1
end do
out=fact(ii,jj)

9000 continue

end
!--------end-of-gfact---------------------------------------------


!-----------------------------------------------------------------
subroutine gmult_f(a, ia, ja, na, splits, fact, nfact, out)
implicit none
integer ia(*), ja(*), na, nfact, splits(nfact+1)
double precision a(na), fact(nfact,nfact), out(na)
!-----------------------------------------------------------------
!     purpose:
!     --------
!     block multipli the entries of a sparse matrix
! on entry:
!----------
!     a      = entries
!     ia     = colindices
!     ja     = rowpointer
!     splits = splits, length is nfact+1
!     fact   = matrix(nfact, nfact)
!     nfact  = ncol(fact)=nrow(fact)
! return:
!-----------
!     out    = modified entries
!----------------------------------------------------------------
integer ii, ri
double precision f
do ii = 1, na
call gri(ii, ja, ri)
call gfact(ri, ia(ii), splits, fact, nfact, f)
out(ii) = a(ii) * f
end do
return
end
!--------end-of-gmult---------------------------------------------
