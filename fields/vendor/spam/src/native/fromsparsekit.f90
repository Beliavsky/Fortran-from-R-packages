!-----------------------------------------------------------------------
subroutine amask (nrow,ncol,a,ja,ia,jmask,imask, &
& c,jc,ic,nzmax,ierr)
!---------------------------------------------------------------------
implicit none
real(8) a(*),c(*)
integer nrow, ncol
integer ia(nrow+1),ja(*),jc(*),ic(nrow+1),jmask(*)
integer imask(nrow+1), nzmax, ierr
logical iw(ncol)
!-----------------------------------------------------------------------
! futher used variables
integer k, len, ii, j, k1, k2
!-----------------------------------------------------------------------
! This subroutine builds a sparse matrix from an input matrix by
! extracting only elements in positions defined by the mask jmask, imask
!-----------------------------------------------------------------------
! On entry:
!---------
! nrow  = integer. row dimension of input matrix
! ncol  = integer. Column dimension of input matrix.
!
! a,
! ja,
! ia    = matrix in Compressed Sparse Row format
!
! jmask,
! imask = matrix defining mask (pattern only) stored in compressed
!         sparse row format.
!
! nzmax = length of arrays c and jc. see ierr.
!
! On return:
!-----------
!
! a, ja, ia and jmask, imask are unchanged.
!
! c
! jc,
! ic    = the output matrix in Compressed Sparse Row format.
!
! ierr  = integer. serving as error message.c
!         ierr = 1  means normal return
!         ierr .gt. 1 means that amask stopped when processing
!         row number ierr, because there was not enough space in
!         c, jc according to the value of nzmax.
!
! work arrays:
!-------------
! iw    = logical work array of length ncol.
!
! note:
!------ the  algorithm is in place: c, jc, ic can be the same as
! a, ja, ia in which cas the code will overwrite the matrix c
! on a, ja, ia
!
!-----------------------------------------------------------------------
ierr = 0
len = 0
do 1 j=1, ncol
iw(j) = .false.
1 continue
!     unpack the mask for row ii in iw
do 100 ii=1, nrow
!     save pointer in order to be able to do things in place
do 2 k=imask(ii), imask(ii+1)-1
iw(jmask(k)) = .true.
2 continue
!     add umasked elemnts of row ii
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
200 continue
!
do 3 k=imask(ii), imask(ii+1)-1
iw(jmask(k)) = .false.
3 continue
100 continue
ic(nrow+1)=len+1
!
return
!-----end-of-amask -----------------------------------------------------
!-----------------------------------------------------------------------
end
!-----------------------------------------------------------------------
subroutine aplsb1 (nrow,ncol,a,ja,ia,s,b,jb,ib,c,jc,ic, &
& nzmax,ierr)
implicit none
real(8) a(*), b(*), c(*), s
integer nrow, ncol
integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1)
integer nzmax, ierr
!-----------------------------------------------------------------------
! further used vaeriables
integer i, j1, j2, ka, kamax, kb, kbmax, kc

!-----------------------------------------------------------------------
! performs the operation C = A+s B for matrices in sorted CSR format.
! the difference with aplsb is that the resulting matrix is such that
! the elements of each row are sorted with increasing column indices in
! each row, provided the original matrices are sorted in the same way.
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow  = integer. The row dimension of A and B
! ncol  = integer. The column dimension of A and B.
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format with entries sorted
!
! s     = real. scalar factor for B.
!
! b,
! jb,
! ib    =  Matrix B in compressed sparse row format with entries sorted
!        ascendly in each row
!
! nzmax = integer. The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!----------
! c,
! jc,
! ic    = resulting matrix C in compressed sparse row sparse format
!         with entries sorted ascendly in each row.
!
! ierr  = integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr .gt. 0 means that amub stopped while computing the
!         i-th row  of C with i=ierr, because the number
!         of elements in C exceeds nzmax.
!
! Notes:
!-------
!     this will not work if any of the two input matrices is not sorted
!-----------------------------------------------------------------------
ierr = 0
kc = 1
ic(1) = kc
!
!     the following loop does a merge of two sparse rows + adds  them.
!
do 6 i=1, nrow
ka = ia(i)
kb = ib(i)
kamax = ia(i+1)-1
kbmax = ib(i+1)-1
5 continue
!
!     this is a while  -- do loop --
!
if (ka .le. kamax .or. kb .le. kbmax) then
!
if (ka .le. kamax) then
j1 = ja(ka)
else
!     take j1 large enough  that always j2 .lt. j1
j1 = ncol+1
endif
if (kb .le. kbmax) then
j2 = jb(kb)
else
!     similarly take j2 large enough  that always j1 .lt. j2
j2 = ncol+1
endif
!
!     three cases
!
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
!
!     end while loop
!
endif
ic(i+1) = kc
6 continue
return
999 ierr = i
return
!------------end-of-aplsb1 ---------------------------------------------
!-----------------------------------------------------------------------
end
!

subroutine submat (job,i1,i2,j1,j2,a,ja,ia,nr,nc,ao,jao,iao)
implicit none
integer job,i1,i2,j1,j2,nr,nc,ia(*),ja(*),jao(*),iao(*)
real(8) a(*),ao(*)
!-----------------------------------------------------------------------
! further used variables
integer i, ii, j, k, k1, k2, klen
!-----------------------------------------------------------------------
! extracts the submatrix A(i1:i2,j1:j2) and puts the result in
! matrix ao,iao,jao
!---- In place: ao,jao,iao may be the same as a,ja,ia.
!--------------
! on input
!---------
! n     = row dimension of the matrix
! i1,i2 = two integers with i2 .ge. i1 indicating the range of rows to be
!          extracted.
! j1,j2 = two integers with j2 .ge. j1 indicating the range of columns
!         to be extracted.
!         * There is no checking whether the input values for i1, i2, j1,
!           j2 are between 1 and n.
! a,
! ja,
! ia    = matrix in compressed sparse row format.
!
! job   = job indicator: if job .ne. 1 then the real values in a are NOT
!         extracted, only the column indices (i.e. data structure) are.
!         otherwise values as well as column indices are extracted...
!
! on output
!--------------
! nr    = number of rows of submatrix
! nc    = number of columns of submatrix
!         * if either of nr or nc is nonpositive the code will quit.
!
! ao,
! jao,iao = extracted matrix in general sparse format with jao containing
!       the column indices,and iao being the pointer to the beginning
!       of the row,in arrays a,ja.
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
nr = i2-i1+1
nc = j2-j1+1
!
if ( nr .le. 0 .or. nc .le. 0) return
!
klen = 0
!
!     simple procedure. proceeds row-wise...
!
do 100 i = 1,nr
ii = i1+i-1
k1 = ia(ii)
k2 = ia(ii+1)-1
iao(i) = klen+1
!-----------------------------------------------------------------------
do 60 k=k1,k2
j = ja(k)
if (j .ge. j1 .and. j .le. j2) then
klen = klen+1
if (job .eq. 1) ao(klen) = a(k)
jao(klen) = j - j1+1
endif
60 continue
100 continue
iao(nr+1) = klen+1
return
!------------end-of submat----------------------------------------------
!-----------------------------------------------------------------------
end
!
subroutine amux (n, x, y, a,ja,ia)
implicit none
real(8)  x(*), y(*), a(*)
integer n, ja(*), ia(*)
!-----------------------------------------------------------------------
!         A times a vector
!-----------------------------------------------------------------------
! multiplies a matrix by a vector using the dot product form
! Matrix A is stored in compressed sparse row storage.
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=Ax
!
!-----------------------------------------------------------------------
! local variables
!
real(8) t
integer i, k
!-----------------------------------------------------------------------
do 100 i = 1,n
!
!     compute the inner product of row i with vector x
!
t = 0.0d0
do 99 k=ia(i), ia(i+1)-1
t = t + a(k)*x(ja(k))
99 continue
!
!     store result in y(i)
!
y(i) = t
100 continue
!
return
!---------end-of-amux---------------------------------------------------
!-----------------------------------------------------------------------
end
!
subroutine amubdg (nrow,ncol,ncolb,ja,ia,jb,ib,ndegr,nnz,iw)
implicit none
integer nrow, ncol, ncolb
integer ja(*),jb(*),ia(nrow+1),ib(ncol+1)
integer ndegr(nrow),iw(ncolb)
integer nnz
!-----------------------------------------------------------------------
! further used variables
integer ii, j, jc, jr, k, last, ldg
!-----------------------------------------------------------------------
! gets the number of nonzero elements in each row of A*B and the total
! number of nonzero elements in A*B.
!-----------------------------------------------------------------------
! on entry:
! --------
!
! nrow  = integer.  row dimension of matrix A
! ncol  = integer.  column dimension of matrix A = row dimension of
!                   matrix B.
! ncolb = integer. the colum dimension of the matrix B.
!
! ja, ia= row structure of input matrix A: ja = column indices of
!         the nonzero elements of A stored by rows.
!         ia = pointer to beginning of each row  in ja.
!
! jb, ib= row structure of input matrix B: jb = column indices of
!         the nonzero elements of A stored by rows.
!         ib = pointer to beginning of each row  in jb.
!
! on return:
! ---------
! ndegr = integer array of length nrow containing the degrees (i.e.,
!         the number of nonzeros in  each row of the matrix A * B
!
! nnz   = total number of nonzero elements found in A * B
!
! work arrays:
!-------------
! iw    = integer work array of length ncolb.
!-----------------------------------------------------------------------
do 1 k=1, ncolb
iw(k) = 0
1 continue

do 2 k=1, nrow
ndegr(k) = 0
2 continue
!
!     method used: Transp(A) * A = sum [over i=1, nrow]  a(i)^T a(i)
!     where a(i) = i-th row of  A. We must be careful not to add  the
!     elements already accounted for.
!
!
do 7 ii=1,nrow
!
!     for each row of A
!
ldg = 0
!
!    end-of-linked list
!
last = -1
do 6 j = ia(ii),ia(ii+1)-1
!
!     row number to be added:
!
jr = ja(j)
do 5 k=ib(jr),ib(jr+1)-1
jc = jb(k)
if (iw(jc) .eq. 0) then
!
!     add one element to the linked list
!
ldg = ldg + 1
iw(jc) = last
last = jc
endif
5 continue
6 continue
ndegr(ii) = ldg
!
!     reset iw to zero
!
do 61 k=1,ldg
j = iw(last)
iw(last) = 0
last = j
61 continue
!-----------------------------------------------------------------------
7 continue
!
nnz = 0
do 8 ii=1, nrow
nnz = nnz+ndegr(ii)
8 continue
!
return
!---------------end-of-amubdg ------------------------------------------
!-----------------------------------------------------------------------
end
!


subroutine amub (nrow,ncol,job,a,ja,ia,b,jb,ib, &
& c,jc,ic,nzmax,iw,ierr)
implicit none
real(8) a(*), b(*), c(*)
integer nrow, ncol, job
integer ja(*),jb(*),jc(*),ia(nrow+1),ib(*),ic(*),iw(ncol)
integer nzmax, ierr
!-----------------------------------------------------------------------
! other used variables
integer len, k, ka, kb, jpos, jcol, ii, jj, j
!-----------------------------------------------------------------------
! performs the matrix by matrix product C = A B
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow  = integer. The row dimension of A = row dimension of C
! ncol  = integer. The column dimension of B = column dimension of C
! job   = integer. Job indicator. When job = 0, only the structure
!                  (i.e. the arrays jc, ic) is computed and the
!                  real values are ignored.
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!
! b,
! jb,
! ib    =  Matrix B in compressed sparse row format.
!
! nzmax = integer. The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!----------
! c,
! jc,
! ic    = resulting matrix C in compressed sparse row sparse format.
!
! ierr  = integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr .gt. 0 means that amub stopped while computing the
!         i-th row  of C with i=ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!------------
! iw    = integer work array of length equal to the number of
!         columns in A.
! Note:
!-------
!   The row dimension of B is not needed. However there is no checking
!   on the condition that ncol(A) = nrow(B).
!
!-----------------------------------------------------------------------
real(8) scal
logical values
values = (job .ne. 0)
!  the following is not necessary... keep [-Wmaybe-uninitialized] quite
scal = 0.0
len = 0
ic(1) = 1
ierr = 0
!     initialize array iw.
do 1 j=1, ncol
iw(j) = 0
1 continue
!
do 500 ii=1, nrow
!     row i
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
100 continue
200 continue
do 201 k=ic(ii), len
iw(jc(k)) = 0
201 continue
ic(ii+1) = len+1
500 continue
return
!-------------end-of-amub-----------------------------------------------
!-----------------------------------------------------------------------
end
!
!------------------------------------------------------------------------
subroutine getl (n,a,ja,ia,ao,jao,iao)
implicit none
integer n, ia(*), ja(*), iao(*), jao(*)
real(8) a(*), ao(*)
!------------------------------------------------------------------------
! this subroutine extracts the lower triangular part of a matrix
! and writes the result ao, jao, iao. The routine is in place in
! that ao, jao, iao can be the same as a, ja, ia if desired.
!-----------
! on input:
!
! n     = dimension of the matrix a.
! a, ja,
!    ia = matrix stored in compressed sparse row format.
! On return:
! ao, jao,
!    iao = lower triangular matrix (lower part of a)
!       stored in a, ja, ia, format
! note: the diagonal element is the last element in each row.
! i.e. in  a(ia(i+1)-1 )
! ao, jao, iao may be the same as a, ja, ia on entry -- in which case
! getl will overwrite the result on a, ja, ia.
!
!------------------------------------------------------------------------
! local variables
real(8) t
integer ko, kold, kdiag, k, i
!
! inititialize ko (pointer for output matrix)
!
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
71 continue
if (kdiag .eq. 0 .or. kdiag .eq. ko) goto 72
!
!     exchange
!
t = ao(kdiag)
ao(kdiag) = ao(ko)
ao(ko) = t
!
k = jao(kdiag)
jao(kdiag) = jao(ko)
jao(ko) = k
72 iao(i) = kold+1
7 continue
!     redefine iao(n+1)
iao(n+1) = ko+1
return
!----------end-of-getl -------------------------------------------------
!-----------------------------------------------------------------------
end
!-----------------------------------------------------------------------
subroutine getu (n,a,ja,ia,ao,jao,iao)
implicit none
integer n, ia(*), ja(*), iao(*), jao(*)
real(8) a(*), ao(*)
!------------------------------------------------------------------------
! this subroutine extracts the upper triangular part of a matrix
! and writes the result ao, jao, iao. The routine is in place in
! that ao, jao, iao can be the same as a, ja, ia if desired.
!-----------
! on input:
!
! n     = dimension of the matrix a.
! a, ja,
!    ia = matrix stored in a, ja, ia, format
! On return:
! ao, jao,
!    iao = upper triangular matrix (upper part of a)
!       stored in compressed sparse row format
! note: the diagonal element is the last element in each row.
! i.e. in  a(ia(i+1)-1 )
! ao, jao, iao may be the same as a, ja, ia on entry -- in which case
! getu will overwrite the result on a, ja, ia.
!
!------------------------------------------------------------------------
! local variables
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
71 continue
if (kdiag .eq. 0 .or. kdiag .eq. kfirst) goto 72
!     exchange
t = ao(kdiag)
ao(kdiag) = ao(kfirst)
ao(kfirst) = t
!
k = jao(kdiag)
jao(kdiag) = jao(kfirst)
jao(kfirst) = k
72 iao(i) = kfirst
7 continue
!     redefine iao(n+1)
iao(n+1) = ko+1
return
!----------end-of-getu -------------------------------------------------
!-----------------------------------------------------------------------
end
!-
subroutine csrmsr (n,a,ja,ia,ao,jao,wk,iwk)
implicit none
integer n
real(8) a(*),ao(*),wk(n)
integer ia(n+1),ja(*),jao(*),iwk(n+1)
!-----------------------------------------------------------------------
! other used variables
integer k, j, iptr, ii, icount, i
!-----------------------------------------------------------------------
! Compressed Sparse Row   to      Modified - Sparse Row
!                                 Sparse row with separate main diagonal
!-----------------------------------------------------------------------
! converts a general sparse matrix a, ja, ia into
! a compressed matrix using a separated diagonal (referred to as
! the bell-labs format as it is used by bell labs semi conductor
! group. We refer to it here as the modified sparse row format.
! Note: this has been coded in such a way that one can overwrite
! the output matrix onto the input matrix if desired by a call of
! the form
!
!     call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
!
! In case ao, jao, are different from a, ja, then one can
! use ao, jao as the work arrays in the calling sequence:
!
!     call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
!
!-----------------------------------------------------------------------
!
! on entry :
!---------
! a, ja, ia = matrix in csr format. note that the
!            algorithm is in place: ao, jao can be the same
!            as a, ja, in which case it will be overwritten on it
!            upon return.
!
! on return :
!-----------
!
! ao, jao  = sparse matrix in modified sparse row storage format:
!          +  ao(1:n) contains the diagonal of the matrix.
!          +  ao(n+2:nnz) contains the nondiagonal elements of the
!             matrix, stored rowwise.
!          +  jao(n+2:nnz) : their column indices
!          +  jao(1:n+1) contains the pointer array for the nondiagonal
!             elements in ao(n+1:nnz) and jao(n+2:nnz).
!             i.e., for i .le. n+1 jao(i) points to beginning of row i
!             in arrays ao, jao.
!              here nnz = number of nonzero elements+1
! work arrays:
!------------
! wk    = real work array of length n
! iwk   = integer work array of length n+1
!
! notes:
!-------
!        Algorithm is in place.  i.e. both:
!
!          call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
!          (in which  ao, jao, are different from a, ja)
!           and
!          call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
!          (in which  wk, jwk, are different from a, ja)
!        are OK.
!--------
! coded by Y. Saad Sep. 1989. Rechecked Feb 27, 1990.
!-----------------------------------------------------------------------
icount = 0
!
! store away diagonal elements and count nonzero diagonal elements.
!
do 1 i=1,n
wk(i) = 0.0d0
iwk(i+1) = ia(i+1)-ia(i)
do 2 k=ia(i),ia(i+1)-1
if (ja(k) .eq. i) then
wk(i) = a(k)
icount = icount + 1
iwk(i+1) = iwk(i+1)-1
endif
2 continue
1 continue
!
! compute total length
!
iptr = n + ia(n+1) - icount
!
!     copy backwards (to avoid collisions)
!
do 500 ii=n,1,-1
do 100 k=ia(ii+1)-1,ia(ii),-1
j = ja(k)
if (j .ne. ii) then
ao(iptr) = a(k)
jao(iptr) = j
iptr = iptr-1
endif
100 continue
500 continue
!
! compute pointer values and copy wk(*)
!
jao(1) = n+2
do 600 i=1,n
ao(i) = wk(i)
jao(i+1) = jao(i)+iwk(i+1)
600 continue
return
!------------ end of subroutine csrmsr ---------------------------------
!-----------------------------------------------------------------------
end
!
subroutine getdia (nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff)
implicit none
real(8) diag(*),a(*)
integer nrow, ncol, job, len, ioff, ia(*), ja(*), idiag(*)
!-----------------------------------------------------------------------
! further used variables
integer izero
!-----------------------------------------------------------------------
! this subroutine extracts a given diagonal from a matrix stored in csr
! format. the output matrix may be transformed with the diagonal removed
! from it if desired (as indicated by job.)
!-----------------------------------------------------------------------
! our definition of a diagonal of matrix is a vector of length nrow
! (always) which contains the elements in rows 1 to nrow of
! the matrix that are contained in the diagonal offset by ioff
! with respect to the main diagonal. if the diagonal element
! falls outside the matrix then it is defined as a zero entry.
! thus the proper definition of diag(*) with offset ioff is
!
!     diag(i) = a(i,ioff+i) i=1,2,...,nrow
!     with elements falling outside the matrix being defined as zero.
!
!-----------------------------------------------------------------------
!
! on entry:
!----------
!
! nrow  = integer. the row dimension of the matrix a.
! ncol  = integer. the column dimension of the matrix a.
! job   = integer. job indicator.  if job = 0 then
!         the matrix a, ja, ia, is not altered on return.
!         if job.ne.0  then getdia will remove the entries
!         collected in diag from the original matrix.
!         this is done in place.
!
! a,ja,
!    ia = matrix stored in compressed sparse row a,ja,ia,format
! ioff  = integer,containing the offset of the wanted diagonal
!         the diagonal extracted is the one corresponding to the
!         entries a(i,j) with j-i = ioff.
!         thus ioff = 0 means the main diagonal
!
! on return:
!-----------
! len   = number of nonzero elements found in diag.
!         (len .le. min(nrow,ncol-ioff)-max(1,1-ioff) + 1 )
!
! diag  = real(8) array of length nrow containing the wanted diagonal.
!         diag contains the diagonal (a(i,j),j-i = ioff ) as defined
!         above.
!
! idiag = integer array of  length len, containing the poisitions
!         in the original arrays a and ja of the diagonal elements
!         collected in diag. a zero entry in idiag(i) means that
!         there was no entry found in row i belonging to the diagonal.
!
! a, ja,
!    ia = if job .ne. 0 the matrix is unchanged. otherwise the nonzero
!         diagonal entries collected in diag are removed from the
!         matrix and therefore the arrays a, ja, ia will change.
!         (the matrix a, ja, ia will contain len fewer elements)
!
!----------------------------------------------------------------------c
!     Y. Saad, sep. 21 1989 - modified and retested Feb 17, 1996.      c
!----------------------------------------------------------------------c
!     local variables
integer istart, max, iend, i, kold, k, kdiag, ko
!
izero = 0
istart = max(izero,-ioff)
iend = min(nrow,ncol-ioff)
len = 0
do 1 i=1,nrow
idiag(i) = 0
diag(i) = 0.0d0
1 continue
!
!     extract  diagonal elements
!
do 6 i=istart+1, iend
do 51 k= ia(i),ia(i+1) -1
if (ja(k)-i .eq. ioff) then
diag(i)= a(k)
idiag(i) = k
len = len+1
goto 6
endif
51 continue
6 continue
if (job .eq. 0 .or. len .eq.0) return
!
!     remove diagonal elements and rewind structure
!
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
71 continue
ia(i) = kold+1
7 continue
!
!     redefine ia(nrow+1)
!
ia(nrow+1) = ko+1
return
!------------end-of-getdia----------------------------------------------
!-----------------------------------------------------------------------
end
!

