
!-----------------------------------------------------------------------
subroutine rowsums(a, ia, nrw, rs)
!-----------------------------------------------------------------------
!     purpose:
!     --------
!
!
!     Reinhard Furrer 2012-04-04
!-----------------------------------------------------------------------
!     parameters:
!     -----------
! on entry:
!----------
!     a, ia = the matrix a in compressed sparse row format (input).
!     nrw = number of rows
!
! on return:
!-----------
!     rs     = rowsums of a
!
! note:
!------
!     no error testing is done. It is assumed that b has enough space
!     allocated.
!-----------------------------------------------------------------------
implicit none

integer ia(*), nrw
double precision a(*), rs(*)
!
!     local variables.
!
integer irw, jja
!
do irw = 1,nrw
do jja = ia(irw),ia(irw+1)-1
rs(irw) = rs(irw)+a(jja)
enddo
!     end irw, we've cycled over all lines
enddo

return
!--------end-of-rowsums------------------------------------------------
!-----------------------------------------------------------------------
end

!-----------------------------------------------------------------------
subroutine rowmeans(a, ia, nrw, ncl, flag, rs)
!-----------------------------------------------------------------------
!     purpose:
!     --------
!       see above
!
!     Reinhard Furrer 2012-04-04
!-----------------------------------------------------------------------
implicit none

integer ia(*), nrw, ncl, flag
double precision a(*), rs(*)
!
!     local variables.
!
integer irw, jja
!
do irw = 1,nrw
do jja = ia(irw),ia(irw+1)-1
rs(irw) = rs(irw)+a(jja)
enddo
if (flag.eq.1) then
if ((ia(irw+1)-ia(irw)).gt.0) then
rs(irw) = rs(irw)/(ia(irw+1)-ia(irw))
endif
else
rs(irw) = rs(irw)/ncl
endif
!     end irw, we've cycled over all lines
enddo

return
!--------end-of-rowmeans------------------------------------------------
!-----------------------------------------------------------------------
end

!-----------------------------------------------------------------------
subroutine colsums(a,ja,ia, nrw, cs)
!-----------------------------------------------------------------------
!     purpose:
!     --------
!        see above
!
!     Reinhard Furrer 2012-04-04
!-----------------------------------------------------------------------
implicit none

integer ia(*),ja(*), nrw
double precision a(*), cs(*)
!
!     local variables.
!
integer ij
!
do ij = 1,ia(nrw+1)-1
cs( ja( ij)) = cs( ja( ij)) + a(ij)
enddo

return
!--------end-of-colsums------------------------------------------------
!-----------------------------------------------------------------------
end

!-----------------------------------------------------------------------
subroutine colmeans(a,ja,ia, nrw, ncl, flag, cs,nnzc)
!-----------------------------------------------------------------------
!     purpose:
!     --------
!        see above
!
!       nnzc needs to be initialized by R!!!!
!     Reinhard Furrer 2012-04-04
!-----------------------------------------------------------------------
implicit none

integer ia(*),ja(*), nrw, ncl, flag, nnzc(ncl)
double precision a(*), cs(*)
!
!     local variables.
!
integer ij
!
do ij = 1,ia(nrw+1)-1
cs( ja( ij)) = cs( ja( ij)) + a(ij)
nnzc( ja( ij)) = nnzc( ja( ij)) + 1
enddo

if (flag.eq.1) then
do ij = 1, ncl
if (nnzc(ij).gt.0) then
cs(ij)=cs(ij)/  nnzc(ij)
endif
enddo
else
do ij = 1, ncl
cs(ij)=cs(ij)/nrw
enddo
endif

return
!--------end-of-colmeans------------------------------------------------
!-----------------------------------------------------------------------
end

