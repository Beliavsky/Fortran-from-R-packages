      subroutine triplet3csr(nrow,ncol,nnz,a,ir,jc,ao,jao,iao,eps)

      implicit none
      double precision a(*),ao(*),eps
      integer nrow,ncol,nnz,ir(*),jc(*),jao(*),iao(*)

      integer kk,k,i,j,tmpi, cr(nrow), ig(nrow+1), g(nnz)
      integer st(nrow+1)
      double precision tmpa(ncol)

C     We assume that we have the correct dimensions.

c     in case we need to determine the max and min
c     also clean up the vectors containing the elements

c     provide empty arrays:
      do kk = 1,nnz
         g(kk) = 0
      enddo
      do kk = 1,nrow
         cr(kk) = 0
      enddo


c     row need to be determined
      k=0
      do kk=1, nnz
         if ((jc(kk).le.ncol).and.(ir(kk).le.nrow)) then
            k=k+1
            if (k.lt.kk) then
               jc(k)=jc(kk)
               ir(k)=ir(kk)
               a(k)=a(kk)
            endif
         endif
      enddo
      nnz=k

      do kk = 1,nnz
         cr(ir(kk)) = cr(ir(kk)) + 1
c         jao(ir(kk)) = cr(ir(kk))
      enddo


c      return
      st(1) = 1
      do kk = 1, nrow
         st(kk+1) = st(kk) + cr(kk)
      enddo

      do kk = 1, nrow
         ig(kk) = st(kk)
      enddo
      do k=1,nnz
         kk = ir(k)
         g(ig(kk)) = k
         ig(kk) = ig(kk) +1
      enddo
c      return


      kk = 0
      iao(1)=1
      do i = 1,nrow
         do j=1,ncol
            tmpa(j)=0.0
         enddo
         do j=1,cr(i)
            tmpi = g(st(i) + j - 1)
            tmpa(jc(tmpi)) = tmpa(jc(tmpi)) + a(tmpi)
         enddo
         do j=1,ncol
            if( abs(tmpa(j)).gt.eps) then
               kk = kk + 1
               ao(kk) = tmpa(j)
               jao(kk) = j
            endif
         enddo
         iao(i+1)= kk+1
      enddo

      nnz = kk
      return
      end

c-----------------------------------------------------------------------
      subroutine triplet2csr(nrow,ncol,nnz,a,ir,jc,ao,jao,iao,eps)

      implicit none
      double precision a(*),ao(*),eps
      integer nrow,ncol,nnz,ir(*),jc(*),jao(*),iao(*)

      integer     newnnz,   ipos,   k,  i, j, tmp1, tmp2
      double precision tmp
c-----------------------------------------------------------------------
c  Triplet representation     to   Compressed Sparse Row
c  Similar to coocsr from sparsekit
c-----------------------------------------------------------------------
c converts a matrix that is stored in coordinate format
c  a, ir, jc into a row general sparse ao, jao, iao format.
c
c on entry:
c---------
c nrow  = row dimension of matrix
c nrow  = col dimension of matrix
c nnz   = number of nonzero elements in matrix
c a,
c ir,
c jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
c         nonzero elements of the matrix with a(k) = actual real value of
c         the elements, ir(k) = its row number and jc(k) = its column
c         number. The order of the elements is arbitrary.
c
c on return:
c-----------
c nnz   = number of nonzero elements in matrix
c ao, jao, iao = matrix in general sparse matrix format with ao
c       continung the real values, jao containing the column indices,
c       and iao being the pointer to the beginning of the row,
c       in arrays ao, jao.
c
c------------------------------------------------------------------------

c     cycle over all entries and count the number of elements in each row
c     skip if larger than nrow and ncol. newnnz is actual number within
c     matrix(nrow,ncol).
      newnnz = 0
      do  k=1, nnz
         tmp1 = ir(k)
         if (tmp1 .le. nrow) then
            tmp2 = jc(k)
            if (tmp2 .le. ncol) then
               if (abs(a(k)) .gt. eps) then
                  iao(tmp1) = iao(tmp1)+1
                  newnnz = newnnz + 1
                  if (newnnz.lt.k) then
                    jc(newnnz) = tmp2
                    ir(newnnz) = tmp1
                    a(newnnz) = a(k)
                  endif
               endif
            endif
         endif
      enddo
c Starting position of each row, essentially a cumsum of iao
      k = 1
      do j=1,nrow+1
         tmp1 = iao(j)
         iao(j) = k
         k = k + tmp1
      enddo
c Go through the structure  once more. Fill in output matrix.
c iao is miss used.
      do  k=1, newnnz
         i = ir(k)
         tmp1 = iao(i)
         ao(tmp1) = a(k)
         jao(tmp1) = jc(k)
         iao(i) = tmp1+1
      enddo
c Shift back iao
      do j=nrow,1,-1
         iao(j+1) = iao(j)
      enddo
      iao(1) = 1

c Sort the individual rows
      do i = 1, nrow
         do  ipos = iao(i), iao(i+1)-1
            do  j = iao(i+1)-1, ipos+1, -1
               k = j - 1
               if (jao(k).eq.jao(j)) then
                  ao(k) = ao(k)+ao(j)
                  ao(j) = 0.0
               else
                  if (jao(k).gt.jao(j)) then
                     tmp1 = jao(k)
                     jao(k) = jao(j)
                     jao(j) = tmp1
                     tmp = ao(k)
                     ao(k) = ao(j)
                     ao(j) = tmp
                  endif
               endif
            enddo
         enddo
      enddo


      call cleanspam(nrow,ao,jao,iao,eps)
      nnz = iao(nrow+1)-1
      return
c-----------------------------------------------------------------------
      end



      subroutine cleanspam(nrow,a,ja,ia,eps)

      implicit none
      integer nrow, ia(nrow+1), ja(*)
      double precision  a(*), eps
c
c     this routine removes zero entries. for more complicated cleaning
c     use the sparsekit2 subroutine clncsr.
c
c On entry:
c----------
c     nrow    -- row dimension of the matrix
c     a,ja,ia -- input matrix in CSR format
c
c On return:
c-----------
c     a,ja,ia -- cleaned matrix
c
c Notes:
c-------
c     Reinhard Furrer 2006-09-13
c-----------------------------------------------------------------------
c
c     Local
      integer i,j,k, oldia(nrow+1)

      do  i = 1, nrow+1
         oldia(i) = ia(i)
      enddo

      k = 1
      do i = 1, nrow
         ia(i) = k
         do j=oldia(i),oldia(i+1)-1
            if (.not.(dabs(a(j)) .le. eps)) then

               ja(k) = ja(j)
               a(k) = a(j)
               k = k + 1
            endif

         enddo
      enddo

      ia(nrow+1) = k
      return

c---- end of cleanspam -------------------------------------------------
c-----------------------------------------------------------------------
      end



      subroutine circulant(nrow,len, x,j, a,ja,ia)

      implicit none
      integer nrow, len, ia(nrow+1), ja(*), j(len)
      double precision  a(*), x(len)
c
c
c
c On entry:
c----------
c     nrow    -- row dimension of the matrix
c     len     -- #nnz per line
c     x,j     -- input values and indicies
c
c On return:
c-----------
c     a,ja,ia -- cleaned circulant matrix
c
c Notes:
c-------
c     Reinhard Furrer 2011-08-03
c-----------------------------------------------------------------------
c
c     Local
      integer i,k, kk

      kk = 1
      ia(1) = 1

      do  i = 1, nrow
         ia(i+1) = ia(i)+len

         do k = 1, len
            ja(kk) = mod( j(k) +i-2, nrow)+1
            a(kk)  = x(k)
            kk = kk+1
         enddo
      enddo
      call sortrows(nrow,a, ja, ia)
      return

c---- end of circulant -------------------------------------------------
c-----------------------------------------------------------------------
      end

      subroutine toeplitz(nrow,len, x,j, a,ja,ia,kk)

      implicit none
      integer nrow, len, ia(nrow+1), ja(*), j(len), kk
      double precision  a(*), x(len)
c
c
c
c On entry:
c----------
c     nrow    -- row dimension of the matrix
c     len     -- total #nnz per line and column
c     x,j     -- input values and indicies (indices are shifted!)
c
c On return:
c-----------
c     a,ja,ia -- toeplitz matrix
c     kk     -- nonzero elements
c
c Notes:
c-------
c     Reinhard Furrer 2011-08-03
c-----------------------------------------------------------------------
c
c     Local
      integer i,k, newj

      kk = 1
      ia(1) = 1

      do  i = 1, nrow

         do k = 1, len
            newj = j(k) + i - nrow
            if ((newj.ge.1).and.(newj.le.nrow)) then
               ja(kk) = newj
               a(kk)  = x(k)
               kk = kk+1
            endif

         enddo
         ia(i+1) = kk
      enddo
      kk = kk - 1
      return

c---- end of toeplitz -------------------------------------------------
c-----------------------------------------------------------------------
      end

c-----------------------------------------------------------------------
      subroutine sortrows(nrow,a,ja,ia)

      implicit none
      integer nrow
      integer ia(nrow+1),ja(*)
      double precision  a(*)
c
c     sorts the rows according to column entries
c
c On entry:
c----------
c     nrow    -- row dimension of the matrix
c     a,ja,ia -- input matrix in sparse format
c
c On return:
c-----------
c     a,ja,ia -- cleaned matrix
c
c Notes:
c-------
c     Reinhard Furrer 2006-09-13
c-----------------------------------------------------------------------
c     Local variables
      integer i,j,k,ko,ipos
      double precision  tmp
c
c

c     .. order the entries according to column indices
c     burble-sort is used
c
      do 190 i = 1, nrow
         do 160 ipos = ia(i), ia(i+1)-1
            do 150 j = ia(i+1)-1, ipos+1, -1
               k = j - 1
               if (ja(k).gt.ja(j)) then
                  ko = ja(k)
                  ja(k) = ja(j)
                  ja(j) = ko
                  tmp = a(k)
                  a(k) = a(j)
                  a(j) = tmp
               endif
 150        continue
 160     continue
 190  continue
      return
c---- end of sortrows --------------------------------------------------
c-----------------------------------------------------------------------
      end
