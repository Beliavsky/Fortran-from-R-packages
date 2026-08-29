c-----------------------------------------------------------------------
      subroutine amask (nrow,ncol,a,ja,ia,jmask,imask,
     *                  c,jc,ic,nzmax,ierr)
c---------------------------------------------------------------------
      implicit none
      real(8) a(*),c(*)
      integer nrow, ncol
      integer ia(nrow+1),ja(*),jc(*),ic(nrow+1),jmask(*)
      integer imask(nrow+1), nzmax, ierr
      logical iw(ncol)
c-----------------------------------------------------------------------
c futher used variables
      integer k, len, ii, j, k1, k2
c-----------------------------------------------------------------------
c This subroutine builds a sparse matrix from an input matrix by
c extracting only elements in positions defined by the mask jmask, imask
c-----------------------------------------------------------------------
c On entry:
c---------
c nrow  = integer. row dimension of input matrix
c ncol  = integer. Column dimension of input matrix.
c
c a,
c ja,
c ia    = matrix in Compressed Sparse Row format
c
c jmask,
c imask = matrix defining mask (pattern only) stored in compressed
c         sparse row format.
c
c nzmax = length of arrays c and jc. see ierr.
c
c On return:
c-----------
c
c a, ja, ia and jmask, imask are unchanged.
c
c c
c jc,
c ic    = the output matrix in Compressed Sparse Row format.
c
c ierr  = integer. serving as error message.c
c         ierr = 1  means normal return
c         ierr .gt. 1 means that amask stopped when processing
c         row number ierr, because there was not enough space in
c         c, jc according to the value of nzmax.
c
c work arrays:
c-------------
c iw    = logical work array of length ncol.
c
c note:
c------ the  algorithm is in place: c, jc, ic can be the same as
c a, ja, ia in which cas the code will overwrite the matrix c
c on a, ja, ia
c
c-----------------------------------------------------------------------
      ierr = 0
      len = 0
      do 1 j=1, ncol
         iw(j) = .false.
 1    continue
c     unpack the mask for row ii in iw
      do 100 ii=1, nrow
c     save pointer in order to be able to do things in place
         do 2 k=imask(ii), imask(ii+1)-1
            iw(jmask(k)) = .true.
 2       continue
c     add umasked elemnts of row ii
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         ic(ii) = len+1
         do 200 k=k1,k2
            j = ja(k)
            if (iw(j)) then
               len = len+1
               if (len .gt. nzmax) then
                  ierr = ii
                  return
               endif
               jc(len) = j
               c(len) = a(k)
            endif
 200     continue
c
         do 3 k=imask(ii), imask(ii+1)-1
            iw(jmask(k)) = .false.
 3       continue
 100  continue
      ic(nrow+1)=len+1
c
      return
c-----end-of-amask -----------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine aplsb1 (nrow,ncol,a,ja,ia,s,b,jb,ib,c,jc,ic,
     *     nzmax,ierr)
      implicit none
      real(8) a(*), b(*), c(*), s
      integer nrow, ncol
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1)
      integer nzmax, ierr
c-----------------------------------------------------------------------
c further used vaeriables
      integer i, j1, j2, ka, kamax, kb, kbmax, kc

c-----------------------------------------------------------------------
c performs the operation C = A+s B for matrices in sorted CSR format.
c the difference with aplsb is that the resulting matrix is such that
c the elements of each row are sorted with increasing column indices in
c each row, provided the original matrices are sorted in the same way.
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow  = integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format with entries sorted
c
c s     = real. scalar factor for B.
c
c b,
c jb,
c ib    =  Matrix B in compressed sparse row format with entries sorted
c        ascendly in each row
c
c nzmax = integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number
c         of elements that exceeds exceeds nzmax. See ierr.
c
c on return:
c----------
c c,
c jc,
c ic    = resulting matrix C in compressed sparse row sparse format
c         with entries sorted ascendly in each row.
c
c ierr  = integer. serving as error message.
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number
c         of elements in C exceeds nzmax.
c
c Notes:
c-------
c     this will not work if any of the two input matrices is not sorted
c-----------------------------------------------------------------------
      ierr = 0
      kc = 1
      ic(1) = kc
c
c     the following loop does a merge of two sparse rows + adds  them.
c
      do 6 i=1, nrow
         ka = ia(i)
         kb = ib(i)
         kamax = ia(i+1)-1
         kbmax = ib(i+1)-1
 5       continue
c
c     this is a while  -- do loop --
c
         if (ka .le. kamax .or. kb .le. kbmax) then
c
            if (ka .le. kamax) then
               j1 = ja(ka)
            else
c     take j1 large enough  that always j2 .lt. j1
               j1 = ncol+1
            endif
            if (kb .le. kbmax) then
               j2 = jb(kb)
            else
c     similarly take j2 large enough  that always j1 .lt. j2
               j2 = ncol+1
            endif
c
c     three cases
c
            if (j1 .eq. j2) then
               c(kc) = a(ka)+s*b(kb)
               jc(kc) = j1
               ka = ka+1
               kb = kb+1
               kc = kc+1
            else if (j1 .lt. j2) then
               jc(kc) = j1
               c(kc) = a(ka)
               ka = ka+1
               kc = kc+1
            else if (j1 .gt. j2) then
               jc(kc) = j2
               c(kc) = s*b(kb)
               kb = kb+1
               kc = kc+1
            endif
            if (kc .gt. nzmax) goto 999
            goto 5
c
c     end while loop
c
         endif
         ic(i+1) = kc
 6    continue
      return
 999  ierr = i
      return
c------------end-of-aplsb1 ---------------------------------------------
c-----------------------------------------------------------------------
      end
c

      subroutine submat (job,i1,i2,j1,j2,a,ja,ia,nr,nc,ao,jao,iao)
      implicit none
      integer job,i1,i2,j1,j2,nr,nc,ia(*),ja(*),jao(*),iao(*)
      real(8) a(*),ao(*)
c-----------------------------------------------------------------------
c further used variables
      integer i, ii, j, k, k1, k2, klen
c-----------------------------------------------------------------------
c extracts the submatrix A(i1:i2,j1:j2) and puts the result in
c matrix ao,iao,jao
c---- In place: ao,jao,iao may be the same as a,ja,ia.
c--------------
c on input
c---------
c n     = row dimension of the matrix
c i1,i2 = two integers with i2 .ge. i1 indicating the range of rows to be
c          extracted.
c j1,j2 = two integers with j2 .ge. j1 indicating the range of columns
c         to be extracted.
c         * There is no checking whether the input values for i1, i2, j1,
c           j2 are between 1 and n.
c a,
c ja,
c ia    = matrix in compressed sparse row format.
c
c job   = job indicator: if job .ne. 1 then the real values in a are NOT
c         extracted, only the column indices (i.e. data structure) are.
c         otherwise values as well as column indices are extracted...
c
c on output
c--------------
c nr    = number of rows of submatrix
c nc    = number of columns of submatrix
c         * if either of nr or nc is nonpositive the code will quit.
c
c ao,
c jao,iao = extracted matrix in general sparse format with jao containing
c       the column indices,and iao being the pointer to the beginning
c       of the row,in arrays a,ja.
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
      nr = i2-i1+1
      nc = j2-j1+1
c
      if ( nr .le. 0 .or. nc .le. 0) return
c
      klen = 0
c
c     simple procedure. proceeds row-wise...
c
      do 100 i = 1,nr
         ii = i1+i-1
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         iao(i) = klen+1
c-----------------------------------------------------------------------
         do 60 k=k1,k2
            j = ja(k)
            if (j .ge. j1 .and. j .le. j2) then
               klen = klen+1
               if (job .eq. 1) ao(klen) = a(k)
               jao(klen) = j - j1+1
            endif
 60      continue
 100  continue
      iao(nr+1) = klen+1
      return
c------------end-of submat----------------------------------------------
c-----------------------------------------------------------------------
      end
c
      subroutine amux (n, x, y, a,ja,ia)
      implicit none
      real(8)  x(*), y(*), a(*)
      integer n, ja(*), ia(*)
c-----------------------------------------------------------------------
c         A times a vector
c-----------------------------------------------------------------------
c multiplies a matrix by a vector using the dot product form
c Matrix A is stored in compressed sparse row storage.
c
c on entry:
c----------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=Ax
c
c-----------------------------------------------------------------------
c local variables
c
      real(8) t
      integer i, k
c-----------------------------------------------------------------------
      do 100 i = 1,n
c
c     compute the inner product of row i with vector x
c
         t = 0.0d0
         do 99 k=ia(i), ia(i+1)-1
            t = t + a(k)*x(ja(k))
 99      continue
c
c     store result in y(i)
c
         y(i) = t
 100  continue
c
      return
c---------end-of-amux---------------------------------------------------
c-----------------------------------------------------------------------
      end
c
      subroutine amubdg (nrow,ncol,ncolb,ja,ia,jb,ib,ndegr,nnz,iw)
      implicit none
      integer nrow, ncol, ncolb
      integer ja(*),jb(*),ia(nrow+1),ib(ncol+1)
      integer ndegr(nrow),iw(ncolb)
      integer nnz
c-----------------------------------------------------------------------
c further used variables
      integer ii, j, jc, jr, k, last, ldg
c-----------------------------------------------------------------------
c gets the number of nonzero elements in each row of A*B and the total
c number of nonzero elements in A*B.
c-----------------------------------------------------------------------
c on entry:
c --------
c
c nrow  = integer.  row dimension of matrix A
c ncol  = integer.  column dimension of matrix A = row dimension of
c                   matrix B.
c ncolb = integer. the colum dimension of the matrix B.
c
c ja, ia= row structure of input matrix A: ja = column indices of
c         the nonzero elements of A stored by rows.
c         ia = pointer to beginning of each row  in ja.
c
c jb, ib= row structure of input matrix B: jb = column indices of
c         the nonzero elements of A stored by rows.
c         ib = pointer to beginning of each row  in jb.
c
c on return:
c ---------
c ndegr = integer array of length nrow containing the degrees (i.e.,
c         the number of nonzeros in  each row of the matrix A * B
c
c nnz   = total number of nonzero elements found in A * B
c
c work arrays:
c-------------
c iw    = integer work array of length ncolb.
c-----------------------------------------------------------------------
      do 1 k=1, ncolb
         iw(k) = 0
 1    continue

      do 2 k=1, nrow
         ndegr(k) = 0
 2    continue
c
c     method used: Transp(A) * A = sum [over i=1, nrow]  a(i)^T a(i)
c     where a(i) = i-th row of  A. We must be careful not to add  the
c     elements already accounted for.
c
c
      do 7 ii=1,nrow
c
c     for each row of A
c
         ldg = 0
c
c    end-of-linked list
c
         last = -1
         do 6 j = ia(ii),ia(ii+1)-1
c
c     row number to be added:
c
            jr = ja(j)
            do 5 k=ib(jr),ib(jr+1)-1
               jc = jb(k)
               if (iw(jc) .eq. 0) then
c
c     add one element to the linked list
c
                  ldg = ldg + 1
                  iw(jc) = last
                  last = jc
               endif
 5          continue
 6       continue
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
      nnz = 0
      do 8 ii=1, nrow
         nnz = nnz+ndegr(ii)
 8    continue
c
      return
c---------------end-of-amubdg ------------------------------------------
c-----------------------------------------------------------------------
      end
c


       subroutine amub (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *                  c,jc,ic,nzmax,iw,ierr)
      implicit none
      real(8) a(*), b(*), c(*)
      integer nrow, ncol, job
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(*),ic(*),iw(ncol)
      integer nzmax, ierr
c-----------------------------------------------------------------------
c other used variables
      integer len, k, ka, kb, jpos, jcol, ii, jj, j
c-----------------------------------------------------------------------
c performs the matrix by matrix product C = A B
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow  = integer. The row dimension of A = row dimension of C
c ncol  = integer. The column dimension of B = column dimension of C
c job   = integer. Job indicator. When job = 0, only the structure
c                  (i.e. the arrays jc, ic) is computed and the
c                  real values are ignored.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c b,
c jb,
c ib    =  Matrix B in compressed sparse row format.
c
c nzmax = integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number
c         of elements that exceeds exceeds nzmax. See ierr.
c
c on return:
c----------
c c,
c jc,
c ic    = resulting matrix C in compressed sparse row sparse format.
c
c ierr  = integer. serving as error message.
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw    = integer work array of length equal to the number of
c         columns in A.
c Note:
c-------
c   The row dimension of B is not needed. However there is no checking
c   on the condition that ncol(A) = nrow(B).
c
c-----------------------------------------------------------------------
      real(8) scal
      logical values
      values = (job .ne. 0)
c  the following is not necessary... keep [-Wmaybe-uninitialized] quite
      scal = 0.0
      len = 0
      ic(1) = 1
      ierr = 0
c     initialize array iw.
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
c
      do 500 ii=1, nrow
c     row i
         do 200 ka=ia(ii), ia(ii+1)-1
            if (values) scal = a(ka)
            jj   = ja(ka)
            do 100 kb=ib(jj),ib(jj+1)-1
               jcol = jb(kb)
               jpos = iw(jcol)
               if (jpos .eq. 0) then
                  len = len+1
                  if (len .gt. nzmax) then
                     ierr = ii
                     return
                  endif
                  jc(len) = jcol
                  iw(jcol)= len
                  if (values) c(len)  = scal*b(kb)
               else
                  if (values) c(jpos) = c(jpos) + scal*b(kb)
               endif
 100        continue
 200     continue
         do 201 k=ic(ii), len
            iw(jc(k)) = 0
 201     continue
         ic(ii+1) = len+1
 500  continue
      return
c-------------end-of-amub-----------------------------------------------
c-----------------------------------------------------------------------
      end
c
c------------------------------------------------------------------------
      subroutine getl (n,a,ja,ia,ao,jao,iao)
      implicit none
      integer n, ia(*), ja(*), iao(*), jao(*)
      real(8) a(*), ao(*)
c------------------------------------------------------------------------
c this subroutine extracts the lower triangular part of a matrix
c and writes the result ao, jao, iao. The routine is in place in
c that ao, jao, iao can be the same as a, ja, ia if desired.
c-----------
c on input:
c
c n     = dimension of the matrix a.
c a, ja,
c    ia = matrix stored in compressed sparse row format.
c On return:
c ao, jao,
c    iao = lower triangular matrix (lower part of a)
c       stored in a, ja, ia, format
c note: the diagonal element is the last element in each row.
c i.e. in  a(ia(i+1)-1 )
c ao, jao, iao may be the same as a, ja, ia on entry -- in which case
c getl will overwrite the result on a, ja, ia.
c
c------------------------------------------------------------------------
c local variables
      real(8) t
      integer ko, kold, kdiag, k, i
c
c inititialize ko (pointer for output matrix)
c
      ko = 0
      do  7 i=1, n
         kold = ko
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .gt. i) goto 71
            ko = ko+1
            ao(ko) = a(k)
            jao(ko) = ja(k)
            if (ja(k)  .eq. i) kdiag = ko
 71      continue
         if (kdiag .eq. 0 .or. kdiag .eq. ko) goto 72
c
c     exchange
c
         t = ao(kdiag)
         ao(kdiag) = ao(ko)
         ao(ko) = t
c
         k = jao(kdiag)
         jao(kdiag) = jao(ko)
         jao(ko) = k
 72      iao(i) = kold+1
 7    continue
c     redefine iao(n+1)
      iao(n+1) = ko+1
      return
c----------end-of-getl -------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine getu (n,a,ja,ia,ao,jao,iao)
      implicit none
      integer n, ia(*), ja(*), iao(*), jao(*)
      real(8) a(*), ao(*)
c------------------------------------------------------------------------
c this subroutine extracts the upper triangular part of a matrix
c and writes the result ao, jao, iao. The routine is in place in
c that ao, jao, iao can be the same as a, ja, ia if desired.
c-----------
c on input:
c
c n     = dimension of the matrix a.
c a, ja,
c    ia = matrix stored in a, ja, ia, format
c On return:
c ao, jao,
c    iao = upper triangular matrix (upper part of a)
c       stored in compressed sparse row format
c note: the diagonal element is the last element in each row.
c i.e. in  a(ia(i+1)-1 )
c ao, jao, iao may be the same as a, ja, ia on entry -- in which case
c getu will overwrite the result on a, ja, ia.
c
c------------------------------------------------------------------------
c local variables
      real(8) t
      integer ko, k, i, kdiag, kfirst
      ko = 0
      do  7 i=1, n
         kfirst = ko+1
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .lt. i) goto 71
            ko = ko+1
            ao(ko) = a(k)
            jao(ko) = ja(k)
            if (ja(k)  .eq. i) kdiag = ko
 71      continue
         if (kdiag .eq. 0 .or. kdiag .eq. kfirst) goto 72
c     exchange
         t = ao(kdiag)
         ao(kdiag) = ao(kfirst)
         ao(kfirst) = t
c
         k = jao(kdiag)
         jao(kdiag) = jao(kfirst)
         jao(kfirst) = k
 72      iao(i) = kfirst
 7    continue
c     redefine iao(n+1)
      iao(n+1) = ko+1
      return
c----------end-of-getu -------------------------------------------------
c-----------------------------------------------------------------------
      end
c-
      subroutine csrmsr (n,a,ja,ia,ao,jao,wk,iwk)
      implicit none
      integer n
      real(8) a(*),ao(*),wk(n)
      integer ia(n+1),ja(*),jao(*),iwk(n+1)
c-----------------------------------------------------------------------
c other used variables
      integer k, j, iptr, ii, icount, i
c-----------------------------------------------------------------------
c Compressed Sparse Row   to      Modified - Sparse Row
c                                 Sparse row with separate main diagonal
c-----------------------------------------------------------------------
c converts a general sparse matrix a, ja, ia into
c a compressed matrix using a separated diagonal (referred to as
c the bell-labs format as it is used by bell labs semi conductor
c group. We refer to it here as the modified sparse row format.
c Note: this has been coded in such a way that one can overwrite
c the output matrix onto the input matrix if desired by a call of
c the form
c
c     call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
c
c In case ao, jao, are different from a, ja, then one can
c use ao, jao as the work arrays in the calling sequence:
c
c     call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
c
c-----------------------------------------------------------------------
c
c on entry :
c---------
c a, ja, ia = matrix in csr format. note that the
c            algorithm is in place: ao, jao can be the same
c            as a, ja, in which case it will be overwritten on it
c            upon return.
c
c on return :
c-----------
c
c ao, jao  = sparse matrix in modified sparse row storage format:
c          +  ao(1:n) contains the diagonal of the matrix.
c          +  ao(n+2:nnz) contains the nondiagonal elements of the
c             matrix, stored rowwise.
c          +  jao(n+2:nnz) : their column indices
c          +  jao(1:n+1) contains the pointer array for the nondiagonal
c             elements in ao(n+1:nnz) and jao(n+2:nnz).
c             i.e., for i .le. n+1 jao(i) points to beginning of row i
c             in arrays ao, jao.
c              here nnz = number of nonzero elements+1
c work arrays:
c------------
c wk    = real work array of length n
c iwk   = integer work array of length n+1
c
c notes:
c-------
c        Algorithm is in place.  i.e. both:
c
c          call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
c          (in which  ao, jao, are different from a, ja)
c           and
c          call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
c          (in which  wk, jwk, are different from a, ja)
c        are OK.
c--------
c coded by Y. Saad Sep. 1989. Rechecked Feb 27, 1990.
c-----------------------------------------------------------------------
      icount = 0
c
c store away diagonal elements and count nonzero diagonal elements.
c
      do 1 i=1,n
         wk(i) = 0.0d0
         iwk(i+1) = ia(i+1)-ia(i)
         do 2 k=ia(i),ia(i+1)-1
            if (ja(k) .eq. i) then
               wk(i) = a(k)
               icount = icount + 1
               iwk(i+1) = iwk(i+1)-1
            endif
 2       continue
 1    continue
c
c compute total length
c
      iptr = n + ia(n+1) - icount
c
c     copy backwards (to avoid collisions)
c
      do 500 ii=n,1,-1
         do 100 k=ia(ii+1)-1,ia(ii),-1
            j = ja(k)
            if (j .ne. ii) then
               ao(iptr) = a(k)
               jao(iptr) = j
               iptr = iptr-1
            endif
 100     continue
 500  continue
c
c compute pointer values and copy wk(*)
c
      jao(1) = n+2
      do 600 i=1,n
         ao(i) = wk(i)
         jao(i+1) = jao(i)+iwk(i+1)
 600  continue
      return
c------------ end of subroutine csrmsr ---------------------------------
c-----------------------------------------------------------------------
      end
c
      subroutine getdia (nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff)
      implicit none
      real(8) diag(*),a(*)
      integer nrow, ncol, job, len, ioff, ia(*), ja(*), idiag(*)
c-----------------------------------------------------------------------
c further used variables
      integer izero
c-----------------------------------------------------------------------
c this subroutine extracts a given diagonal from a matrix stored in csr
c format. the output matrix may be transformed with the diagonal removed
c from it if desired (as indicated by job.)
c-----------------------------------------------------------------------
c our definition of a diagonal of matrix is a vector of length nrow
c (always) which contains the elements in rows 1 to nrow of
c the matrix that are contained in the diagonal offset by ioff
c with respect to the main diagonal. if the diagonal element
c falls outside the matrix then it is defined as a zero entry.
c thus the proper definition of diag(*) with offset ioff is
c
c     diag(i) = a(i,ioff+i) i=1,2,...,nrow
c     with elements falling outside the matrix being defined as zero.
c
c-----------------------------------------------------------------------
c
c on entry:
c----------
c
c nrow  = integer. the row dimension of the matrix a.
c ncol  = integer. the column dimension of the matrix a.
c job   = integer. job indicator.  if job = 0 then
c         the matrix a, ja, ia, is not altered on return.
c         if job.ne.0  then getdia will remove the entries
c         collected in diag from the original matrix.
c         this is done in place.
c
c a,ja,
c    ia = matrix stored in compressed sparse row a,ja,ia,format
c ioff  = integer,containing the offset of the wanted diagonal
c         the diagonal extracted is the one corresponding to the
c         entries a(i,j) with j-i = ioff.
c         thus ioff = 0 means the main diagonal
c
c on return:
c-----------
c len   = number of nonzero elements found in diag.
c         (len .le. min(nrow,ncol-ioff)-max(1,1-ioff) + 1 )
c
c diag  = real(8) array of length nrow containing the wanted diagonal.
c         diag contains the diagonal (a(i,j),j-i = ioff ) as defined
c         above.
c
c idiag = integer array of  length len, containing the poisitions
c         in the original arrays a and ja of the diagonal elements
c         collected in diag. a zero entry in idiag(i) means that
c         there was no entry found in row i belonging to the diagonal.
c
c a, ja,
c    ia = if job .ne. 0 the matrix is unchanged. otherwise the nonzero
c         diagonal entries collected in diag are removed from the
c         matrix and therefore the arrays a, ja, ia will change.
c         (the matrix a, ja, ia will contain len fewer elements)
c
c----------------------------------------------------------------------c
c     Y. Saad, sep. 21 1989 - modified and retested Feb 17, 1996.      c
c----------------------------------------------------------------------c
c     local variables
      integer istart, max, iend, i, kold, k, kdiag, ko
c
      izero = 0
      istart = max(izero,-ioff)
      iend = min(nrow,ncol-ioff)
      len = 0
      do 1 i=1,nrow
         idiag(i) = 0
         diag(i) = 0.0d0
 1    continue
c
c     extract  diagonal elements
c
      do 6 i=istart+1, iend
         do 51 k= ia(i),ia(i+1) -1
            if (ja(k)-i .eq. ioff) then
               diag(i)= a(k)
               idiag(i) = k
               len = len+1
               goto 6
            endif
 51      continue
 6    continue
      if (job .eq. 0 .or. len .eq.0) return
c
c     remove diagonal elements and rewind structure
c
      ko = 0
      do  7 i=1, nrow
         kold = ko
         kdiag = idiag(i)
         do 71 k= ia(i), ia(i+1)-1
            if (k .ne. kdiag) then
               ko = ko+1
               a(ko) = a(k)
               ja(ko) = ja(k)
            endif
 71      continue
         ia(i) = kold+1
 7    continue
c
c     redefine ia(nrow+1)
c
      ia(nrow+1) = ko+1
      return
c------------end-of-getdia----------------------------------------------
c-----------------------------------------------------------------------
      end
c

