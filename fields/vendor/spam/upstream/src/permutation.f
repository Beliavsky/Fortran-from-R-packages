c     It has been tested that embedding the loop over the right hand
c     side into the backsolve routine is not faster.

c-----------------------------------------------------------------------
      subroutine getbwd(n,  ja,ia,ml,mu)
c-----------------------------------------------------------------------
c gets the bandwidth of lower part and upper part of A.
c does not assume that A is sorted.
c-----------------------------------------------------------------------
c on entry:
c----------
c n	= integer = the row dimension of the matrix
c ja, ia = matrix in compressed sparse row format.
c
c on return:
c-----------
c ml	= integer. The bandwidth of the strict lower part of A
c mu	= integer. The bandwidth of the strict upper part of A
c
c Notes:
c ===== ml and mu are allowed to be negative or return. This may be
c       useful since it will tell us whether a band is confined
c       in the strict  upper/lower triangular part.
c       indeed the definitions of ml and mu are
c
c       ml = max ( (i-j)  s.t. a(i,j) .ne. 0  )
c       mu = max ( (j-i)  s.t. a(i,j) .ne. 0  )
c----------------------------------------------------------------------c
c Y. Saad, Sep. 21 1989 | simplified by R. Furrer Sept. 2016           c
c----------------------------------------------------------------------c
      implicit none
c     double precision a(*)
      integer n,ja(*),ia(n+1),ml,mu

      integer ldist,i,k

      ml = - n
      mu = - n
      do 3 i=1,n
         do 31 k=ia(i),ia(i+1)-1
            ldist = i-ja(k)
            ml = max(ml,ldist)
            mu = max(mu,-ldist)
 31      continue
 3    continue
      return
c---------------end-of-getbwd ------------------------------------------
c-----------------------------------------------------------------------
      end


c functions slightly modified from sparsekit:
c cperm,rperm,dperm: job argument is eliminated
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
      subroutine rperm (nrow,a,ja,ia,ao,jao,iao,perm)
      implicit none
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow)
      double precision a(*),ao(*)
c-----------------------------------------------------------------------
c this subroutine permutes the rows of a matrix in CSR format.
c rperm  computes B = P A  where P is a permutation matrix.
c the permutation P is defined through the array perm: for each j,
c perm(j) represents the destination row number of row number j.
c Youcef Saad -- recoded Jan 28, 1991.
c-----------------------------------------------------------------------
c on entry:
c----------
c n 	= dimension of the matrix
c a, ja, ia = input matrix in csr format
c perm 	= integer array of length nrow containing the permutation arrays
c	  for the rows: perm(i) is the destination of row i in the
c         permuted matrix.
c         ---> a(i,j) in the original matrix becomes a(perm(i),j)
c         in the output  matrix.
c
c
c------------
c on return:
c------------
c ao, jao, iao = input matrix in a, ja, ia format
c----------------------------------------------------------------------c
c           Y. Saad, May  2, 1990                                      c
c----------------------------------------------------------------------c

      integer i,j,k,ko,ii
c     determine pointers for output matix.
c
      do 50 j=1,nrow
         i = perm(j)
         iao(i+1) = ia(j+1) - ia(j)
 50   continue
c
c get pointers from lengths
c
      iao(1) = 1
      do 51 j=1,nrow
         iao(j+1)=iao(j+1)+iao(j)
 51   continue
c
c copying
c
      do 100 ii=1,nrow
c
c old row = ii  -- new row = iperm(ii) -- ko = new pointer
c
         ko = iao(perm(ii))
         do 60 k=ia(ii), ia(ii+1)-1
            jao(ko) = ja(k)
            ao(ko) = a(k)
            ko = ko+1
 60      continue
 100  continue
c
      return
c---------end-of-rperm -------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine cperm (nrow,a,ja,ia,ao,jao,iao,perm)
      implicit none
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(*)
      double precision a(*), ao(*)
c-----------------------------------------------------------------------
c this subroutine permutes the columns of a matrix a, ja, ia.
c the result is written in the output matrix  ao, jao, iao.
c cperm computes B = A P, where  P is a permutation matrix
c that maps column j into column perm(j), i.e., on return
c      a(i,j) becomes a(i,perm(j)) in new matrix
c Y. Saad, May 2, 1990 / modified Jan. 28, 1991.
c-----------------------------------------------------------------------
c on entry:
c----------
c nrow 	= row dimension of the matrix
c
c a, ja, ia = input matrix in csr format.
c
c perm	= integer array of length ncol (number of columns of A
c         containing the permutation array  the columns:
c         a(i,j) in the original matrix becomes a(i,perm(j))
c         in the output matrix.
c
c
c------------
c on return:
c------------
c ao, jao, iao = input matrix in a, ja, ia format
c
c Notes:
c-------
c 1. if job=1 then ao, iao are not used.
c 2. This routine is in place: ja, jao can be the same.
c 3. If the matrix is initially sorted (by increasing column number)
c    then ao,jao,iao  may not be on return, hence a call to csort.
c
c----------------------------------------------------------------------c
c local parameters:
      integer k, i, nnz
c
      nnz = ia(nrow+1)-1
      do 100 k=1,nnz
         jao(k) = perm(ja(k))
 100  continue
c
c     done with ja array.
c
      do 1 i=1, nrow+1
         iao(i) = ia(i)
 1    continue
c
      do 2 k=1, nnz
         ao(k) = a(k)
 2    continue
c
      call sortrows(nrow,ao,jao,iao)
c     call csort (nrow,ao,jao,iao,iwork) _does not work_
      return
c---------end-of-cperm--------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine dperm (nrow,a,ja,ia,ao,jao,iao,pperm,qperm)
      implicit none
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),pperm(nrow),
     +        qperm(*)
      double precision a(*),ao(*)
c-----------------------------------------------------------------------
c This routine permutes the rows and columns of a matrix stored in CSR
c format. i.e., it computes P A Q, where P, Q are permutation matrices.
c P maps row i into row perm(i) and Q maps column j into column qperm(j):
c      a(i,j)    becomes   a(pperm(i),qperm(j)) in new matrix
c note that qperm should be of length ncol (number of columns) but this
c is not checked.
c-----------------------------------------------------------------------
c Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991.
c-----------------------------------------------------------------------
c on entry:
c----------
c n 	= dimension of the matrix
c a, ja,
c    ia = input matrix in a, ja, ia format
c pperm 	= integer array of length n containing the permutation arrays
c	  for the rows: pperm(i) is the destination of row i in the
c         permuted matrix
c
c qperm	= same thing for the columns.
c iwork = working array passed to cperm
c
c on return:
c-----------
c ao, jao, iao = input matrix in a, ja, ia format
c
c Notes:
c-------
c  1) algorithm is in place
c----------------------------------------------------------------------c
c local variables
c
c permute rows first
c
      call rperm (nrow,a,ja,ia,   ao,jao,iao,pperm)
c
c then permute columns
c
c
      call cperm (nrow,ao,jao,iao,ao,jao,iao,qperm)
c
      return
c-------end-of-dperm----------------------------------------------------
c-----------------------------------------------------------------------
      end



c-----------------------------------------------------------------------
      subroutine dvperm (n, x, perm)
      implicit none
      integer n
      integer perm(n)
      double precision x(n)
c-----------------------------------------------------------------------
c this subroutine performs an in-place permutation of a real vector x
c according to the permutation array perm(*), i.e., on return,
c the vector x satisfies,
c
c	x(perm(j)) :== x(j), j=1,2,.., n
c
c-----------------------------------------------------------------------
c on entry:
c---------
c n 	= length of vector x.
c perm 	= integer array of length n containing the permutation  array.
c x	= input vector
c
c on return:
c----------
c x	= vector x permuted according to x(perm(*)) :=  x(*)
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
c local variables
      integer init,next,k, ii, j
      double precision tmp, tmp1
c
      init      = 1
      tmp       = x(init)
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
c
c loop
c
 6    k = k+1
c
c save the chased element --
c
      tmp1      = x(ii)
      x(ii)     = tmp
      next      = perm(ii)
      if (next .lt. 0 ) goto 65
c
c test for end
c
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next
c
c end loop
c
      goto 6
c
c reinitilaize cycle --
c
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp = x(init)
      ii = perm(init)
      perm(init)=-perm(init)
      goto 6
c
 101  continue
      do 200 j=1, n
            perm(j) = -perm(j)
 200  continue
c
      return
c-------------------end-of-dvperm---------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine ivperm (n, ix, perm)
      implicit none
      integer n, perm(n+1), ix(n)
c-----------------------------------------------------------------------
c this subroutine performs an in-place permutation of an integer vector
c ix according to the permutation array perm(*), i.e., on return,
c the vector x satisfies,
c
c	ix(perm(j)) :== ix(j), j=1,2,.., n
c
c-----------------------------------------------------------------------
c on entry:
c---------
c n 	= length of vector x.
c perm 	= integer array of length n containing the permutation  array.
c ix	= input vector
c
c on return:
c----------
c ix	= vector x permuted according to ix(perm(*)) :=  ix(*)
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
c local variables
      integer ii,k,j,next,init,tmp, tmp1
c
      init      = 1
      tmp       = ix(init)
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
c
c loop
c
 6    k = k+1
c
c save the chased element --
c
      tmp1     = ix(ii)
      ix(ii)   = tmp
      next     = perm(ii)
      if (next .lt. 0 ) goto 65
c
c test for end
c
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next
c
c end loop
c
      goto 6
c
c reinitilaize cycle --
c
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp       = ix(init)
      ii        = perm(init)
      perm(init)=-perm(init)
      goto 6
c
 101  continue
      do 200 j=1, n
c         if (perm(j) .lt. 0) then
            perm(j) = -perm(j)
c         endif
 200  continue
c
      return
c-------------------end-of-ivperm---------------------------------------
c-----------------------------------------------------------------------
      end
c
c-----------------------------------------------------------------------
      subroutine aplbdg (nrow,ncol,ja,ia,jb,ib,ndegr,nnz,iw)

      implicit none
      integer nrow, ncol, nnz
      integer ja(*),jb(*),ia(nrow+1),ib(nrow+1),iw(ncol),ndegr(nrow)
c-----------------------------------------------------------------------
c gets the number of nonzero elements in each row of A+B and the total
c number of nonzero elements in A+B.
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c b,
c jb,
c ib	=  Matrix B in compressed sparse row format.
c
c iw,nnz,ngegr = zero content
c
c on return:
c----------
c ndegr	= integer array of length nrow containing the degrees (i.e.,
c         the number of nonzeros in  each row of the matrix A + B.
c
c nnz   = total number of nonzero elements found in A * B
c
c work arrays:
c------------
c iw	= integer work array of length equal to ncol.
c
c-----------------------------------------------------------------------
      integer k,j,ii,jr,last,ldg,jc

      do 7 ii=1,nrow
         ldg = 0
c
c    end-of-linked list
c
         last = -1
c
c     row of A
c
         do 5 j = ia(ii),ia(ii+1)-1
            jr = ja(j)
c
c     add element to the linked list
c
            ldg = ldg + 1
            iw(jr) = last
            last = jr
 5       continue
c
c     row of B
c
         do 6 j=ib(ii),ib(ii+1)-1
            jc = jb(j)
            if (iw(jc) .eq. 0) then
c
c     add one element to the linked list
c
               ldg = ldg + 1
               iw(jc) = last
               last = jc
            endif
 6       continue
c     done with row ii.
         ndegr(ii) = ldg
c
c     reset iw to zero
c
         do 61 k=1,ldg
            j = iw(last)
            iw(last) = 0
            last = j
 61      continue
c-----------------------------------------------------------------------
 7    continue
c

      do 8 ii=1, nrow
         nnz = nnz+ndegr(ii)
 8    continue
      return
c----------------end-of-aplbdg -----------------------------------------
c-----------------------------------------------------------------------
      end
