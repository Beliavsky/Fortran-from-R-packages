
c-----------------------------------------------------------------------
      subroutine rowsums(a, ia, nrw, rs)
c-----------------------------------------------------------------------
c     purpose:
c     --------
c
c
c     Reinhard Furrer 2012-04-04
c-----------------------------------------------------------------------
c     parameters:
c     -----------
c on entry:
c----------
c     a, ia = the matrix a in compressed sparse row format (input).
c     nrw = number of rows
c
c on return:
c-----------
c     rs     = rowsums of a
c
c note:
c------
c     no error testing is done. It is assumed that b has enough space
c     allocated.
c-----------------------------------------------------------------------
      implicit none

      integer ia(*), nrw
      double precision a(*), rs(*)
c
c     local variables.
c
      integer irw, jja
c
      do irw = 1,nrw
         do jja = ia(irw),ia(irw+1)-1
            rs(irw) = rs(irw)+a(jja)
         enddo
c     end irw, we've cycled over all lines
      enddo

      return
c--------end-of-rowsums------------------------------------------------
c-----------------------------------------------------------------------
      end

c-----------------------------------------------------------------------
      subroutine rowmeans(a, ia, nrw, ncl, flag, rs)
c-----------------------------------------------------------------------
c     purpose:
c     --------
c       see above
c
c     Reinhard Furrer 2012-04-04
c-----------------------------------------------------------------------
      implicit none

      integer ia(*), nrw, ncl, flag
      double precision a(*), rs(*)
c
c     local variables.
c
      integer irw, jja
c
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
c     end irw, we've cycled over all lines
      enddo

      return
c--------end-of-rowmeans------------------------------------------------
c-----------------------------------------------------------------------
      end

c-----------------------------------------------------------------------
      subroutine colsums(a,ja,ia, nrw, cs)
c-----------------------------------------------------------------------
c     purpose:
c     --------
c        see above
c
c     Reinhard Furrer 2012-04-04
c-----------------------------------------------------------------------
      implicit none

      integer ia(*),ja(*), nrw
      double precision a(*), cs(*)
c
c     local variables.
c
      integer ij
c
      do ij = 1,ia(nrw+1)-1
         cs( ja( ij)) = cs( ja( ij)) + a(ij)
       enddo

      return
c--------end-of-colsums------------------------------------------------
c-----------------------------------------------------------------------
      end

c-----------------------------------------------------------------------
      subroutine colmeans(a,ja,ia, nrw, ncl, flag, cs,nnzc)
c-----------------------------------------------------------------------
c     purpose:
c     --------
c        see above
c
c       nnzc needs to be initialized by R!!!!
c     Reinhard Furrer 2012-04-04
c-----------------------------------------------------------------------
      implicit none

      integer ia(*),ja(*), nrw, ncl, flag, nnzc(ncl)
      double precision a(*), cs(*)
c
c     local variables.
c
      integer ij
c
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
c--------end-of-colmeans------------------------------------------------
c-----------------------------------------------------------------------
      end

