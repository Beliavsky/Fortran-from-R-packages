c-----------------------------------------------------------------
      subroutine gri(i,rpin,rindex)
      implicit none
      integer j, i, rpin(*), rindex
c-----------------------------------------------------------------
c     purpose:
c     --------
c     get row index for the i-th entrie of the sparse matrix
c on entry:
c----------
c     i      = number of the element
c     rpin   = input rowpointers vector
c on return:
c-----------
c     rindex  = output row index
c-----------------------------------------------------------------
      j=1
      do while (rpin(j) <= i)
         j=j+1
      end do
      rindex=j-1
      return
      end
c--------end-of-gri---------------------------------------------

c-----------------------------------------------------------------
      subroutine gfact(i,j,splits,fact,nfact,out)
      implicit none
      integer j, i, nfact, splits(nfact+1)
      double precision fact(nfact,nfact), out
c-----------------------------------------------------------------
c     purpose:
c     --------
c     get fact for coordinates i,j,
c on entry:
c----------
c     i,j    = indices of element
c     splits = splits of the matrix
c     fact   = values to be returned
c     nfact  = ncol(fact)=nrow(fact)
c return:
c-----------
c     out    = the correct value of fact
c----------------------------------------------------------------
      integer ii, jj
      if (i >= splits(nfact+1) .OR. j >= splits(nfact+1)) then
         goto 9000
c         stop 'i,j are out of larger than (nfact+1)'
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
c--------end-of-gfact---------------------------------------------


c-----------------------------------------------------------------
      subroutine gmult_f(a, ia, ja, na, splits, fact, nfact, out)
      implicit none
      integer ia(*), ja(*), na, nfact, splits(nfact+1)
      double precision a(na), fact(nfact,nfact), out(na)
c-----------------------------------------------------------------
c     purpose:
c     --------
c     block multipli the entries of a sparse matrix
c on entry:
c----------
c     a      = entries
c     ia     = colindices
c     ja     = rowpointer
c     splits = splits, length is nfact+1
c     fact   = matrix(nfact, nfact)
c     nfact  = ncol(fact)=nrow(fact)
c return:
c-----------
c     out    = modified entries
c----------------------------------------------------------------
      integer ii, ri
      double precision f
      do ii = 1, na
         call gri(ii, ja, ri)
         call gfact(ri, ia(ii), splits, fact, nfact, f)
         out(ii) = a(ii) * f
      end do
      return
      end
c--------end-of-gmult---------------------------------------------
