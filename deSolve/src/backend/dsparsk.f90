!  The code in this file is was taken from daspk.tgz from
!    https://www.netlib.org/ode/
!  Authors: Brown, Hindmarsh, Petzold
!    originating from SPARSKIT, version 1 by Yousef Saad
!  Adapted for use in R package deSolve by the deSolve authors.
!

!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!        BASIC LINEAR ALGEBRA FOR SPARSE MATRICES. BLASSM MODULE       c
!----------------------------------------------------------------------c
! aplb   :   computes     C = A+B                                      c
! aplb1  :   computes     C = A+B  [Sorted version: A, B, C sorted]    c
! aplsb  :   computes     C = A + s B                                  c
! diamua :   Computes     C = Diag * A                                 c
! amudia :   Computes     C = A* Diag                                  c
! aplsca :   Computes     A:= A + s I    (s = scalar)                  c
!----------------------------------------------------------------------c
subroutine diamua (nrow,job, a, ja, ia, diag, b, jb, ib)
real(kind=kind(0.0d0)) a(*), b(*), diag(nrow), scal
integer ja(*),jb(*), ia(nrow+1),ib(nrow+1)
!-----------------------------------------------------------------------
! performs the matrix by matrix product B = Diag * A  (in place)
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow  = integer. The row dimension of A
!
! job   = integer. job indicator. Job=0 means get array b only
!         job = 1 means get b, and the integer arrays ib, jb.
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!
! diag = diagonal matrix stored as a vector dig(1:n)
!
! on return:
!----------
!
! b,
! jb,
! ib   = resulting matrix B in compressed sparse row sparse format.
!
! Notes:
!-------
! 1)        The column dimension of A is not needed.
! 2)        algorithm in place (B can take the place of A).
!           in this case use job=0.
!-----------------------------------------------------------------
do 1 ii=1,nrow
!
!     normalize each row
!
   k1 = ia(ii)
   k2 = ia(ii+1)-1
   scal = diag(ii)
   do 2 k=k1, k2
      b(k) = a(k)*scal
2 continue
1 continue
!
if (job .EQ. 0) return
!
do 3 ii=1, nrow+1
   ib(ii) = ia(ii)
3 continue
do 31 k=ia(1), ia(nrow+1) -1
   jb(k) = ja(k)
31 continue
return
!----------end-of-diamua------------------------------------------------
!-----------------------------------------------------------------------
end
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!          BASIC MATRIX-VECTOR OPERATIONS - MATVEC MODULE              c
!----------------------------------------------------------------------c
! amux  : A times a vector. Compressed Sparse Row (CSR) format.        c
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!                        INPUT-OUTPUT MODULE                           c
!----------------------------------------------------------------------c
!  prtmt  : prints matrices in the Boeing/Harwell format.              c
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!                    FORMAT CONVERSION MODULE                          c
!----------------------------------------------------------------------c
! csrdns  : converts a row-stored sparse matrix into the dense format. c
! coocsr  : converts coordinate to  to csr format                      c
! coicsr  : in-place conversion of coordinate to csr format            c
! csrcoo  : converts compressed sparse row to coordinate.              c
! csrcsc  : converts compressed sparse row format to compressed sparse c
!           column format (transposition)                              c
! csrcsc2 : rectangular version of csrcsc                              c
! csrdia  : converts a compressed sparse row format into a diagonal    c
!           format.                                                    c
! csrbnd  : converts a compressed sparse row format into a banded      c
!           format (linpack style).                                    c
!----------------------------------------------------------------------c
subroutine csrcsc (n,job,ipos,a,ja,ia,ao,jao,iao)
integer ia(n+1),iao(n+1),ja(*),jao(*)
real(kind=kind(0.0d0))  a(*),ao(*)
!-----------------------------------------------------------------------
! Compressed Sparse Row     to      Compressed Sparse Column
!
! (transposition operation)   Not in place.
!-----------------------------------------------------------------------
! -- not in place --
! this subroutine transposes a matrix stored in a, ja, ia format.
! ---------------
! on entry:
!----------
! n     = dimension of A.
! job   = integer to indicate whether to fill the values (job.eq.1) of the
!         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use for example
!                call csrcsc (n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2))
!         for any other normal usage, enter ipos=1.
! a     = real array of length nnz (nnz=number of nonzero elements in input
!         matrix) containing the nonzero elements.
! ja    = integer array of length nnz containing the column positions
!         of the corresponding elements in a.
! ia    = integer of size n+1. ia(k) contains the position in a, ja of
!         the beginning of the k-th row.
!
! on return:
! ----------
! output arguments:
! ao    = real array of size nzz containing the "a" part of the transpose
! jao   = integer array of size nnz containing the column indices.
! iao   = integer array of size n+1 containing the "ia" index array of
!         the transpose.
!
!-----------------------------------------------------------------------
call csrcsc2 (n,n,job,ipos,a,ja,ia,ao,jao,iao)
end
subroutine csrcsc2 (n,n2,job,ipos,a,ja,ia,ao,jao,iao)
integer ia(n+1),iao(n2+1),ja(*),jao(*)
real(kind=kind(0.0d0))  a(*),ao(*)
!-----------------------------------------------------------------------
! Compressed Sparse Row     to      Compressed Sparse Column
!
! (transposition operation)   Not in place.
!-----------------------------------------------------------------------
! Rectangular version.  n is number of rows of CSR matrix,
!                       n2 (input) is number of columns of CSC matrix.
!-----------------------------------------------------------------------
! -- not in place --
! this subroutine transposes a matrix stored in a, ja, ia format.
! ---------------
! on entry:
!----------
! n     = number of rows of CSR matrix.
! n2    = number of columns of CSC matrix.
! job   = integer to indicate whether to fill the values (job.eq.1) of the
!         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use for example
!                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2))
!         for any other normal usage, enter ipos=1.
! a     = real array of length nnz (nnz=number of nonzero elements in input
!         matrix) containing the nonzero elements.
! ja    = integer array of length nnz containing the column positions
!         of the corresponding elements in a.
! ia    = integer of size n+1. ia(k) contains the position in a, ja of
!         the beginning of the k-th row.
!
! on return:
! ----------
! output arguments:
! ao    = real array of size nzz containing the "a" part of the transpose
! jao   = integer array of size nnz containing the column indices.
! iao   = integer array of size n+1 containing the "ia" index array of
!         the transpose.
!
!-----------------------------------------------------------------------
!----------------- compute lengths of rows of transp(A) ----------------
do 1 i=1,n2+1
   iao(i) = 0
1 continue
do 3 i=1, n
   do 2 k=ia(i), ia(i+1)-1
      j = ja(k)+1
      iao(j) = iao(j)+1
2 continue
3 continue
!---------- compute pointers from lengths ------------------------------
iao(1) = ipos
do 4 i=1,n2
   iao(i+1) = iao(i) + iao(i+1)
4 continue
!--------------- now do the actual copying -----------------------------
do 6 i=1,n
   do 62 k=ia(i),ia(i+1)-1
      j = ja(k)
      next = iao(j)
      if (job .EQ. 1)  ao(next) = a(k)
      jao(next) = i
      iao(j) = next+1
62 continue
6 continue
!-------------------------- reshift iao and leave ----------------------
do 7 i=n2,1,-1
   iao(i+1) = iao(i)
7 continue
iao(1) = ipos
!--------------- end of csrcsc2 ----------------------------------------
!-----------------------------------------------------------------------
end
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!                     UNARY SUBROUTINES MODULE                         c
!----------------------------------------------------------------------c
! rperm  : permutes the rows of a matrix (B = P A)                     c
! cperm  : permutes the columns of a matrix (B = A Q)                  c
! dperm  : permutes both the rows and columns of a matrix (B = P A Q ) c
! dvperm : permutes a real vector (in-place)                           c
! ivperm : permutes an integer vector (in-place)                       c
! diapos : returns the positions of the diagonal elements in A.        c
! getbwd : returns the bandwidth information on a matrix.              c
! infdia : obtains information on the diagonals of A.                  c
! rnrms  : computes the norms of the rows of A                         c
! roscal : scales the rows of a matrix by their norms.                 c
!----------------------------------------------------------------------c
subroutine rperm (nrow,a,ja,ia,ao,jao,iao,perm,job)
integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),job
real(kind=kind(0.0d0)) a(*),ao(*)
!-----------------------------------------------------------------------
! this subroutine permutes the rows of a matrix in CSR format.
! rperm  computes B = P A  where P is a permutation matrix.
! the permutation P is defined through the array perm: for each j,
! perm(j) represents the destination row number of row number j.
! Youcef Saad -- recoded Jan 28, 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! n      = dimension of the matrix
! a, ja, ia = input matrix in csr format
! perm    = integer array of length nrow containing the permutation arrays
!           for the rows: perm(i) is the destination of row i in the
!           permuted matrix.
!         ---> a(i,j) in the original matrix becomes a(perm(i),j)
!         in the output  matrix.
!
! job     = integer indicating the work to be done:
!       job = 1   permute a, ja, ia into ao, jao, iao
!                       (including the copying of real values ao and
!                       the array iao).
!       job .ne. 1 :  ignore real values.
!                     (in which case arrays a and ao are not needed nor
!                      used).
!
!------------
! on return:
!------------
! ao, jao, iao = input matrix in a, ja, ia format
! note :
!        if (job.ne.1)  then the arrays a and ao are not used.
!----------------------------------------------------------------------c
!           Y. Saad, May  2, 1990                                      c
!----------------------------------------------------------------------c
logical values
values = (job .EQ. 1)
!
!     determine pointers for output matix.
!
do 50 j=1,nrow
   i = perm(j)
   iao(i+1) = ia(j+1) - ia(j)
50 continue
!
! get pointers from lengths
!
iao(1) = 1
do 51 j=1,nrow
   iao(j+1)=iao(j+1)+iao(j)
51 continue
!
! copying
!
do 100 ii=1,nrow
!
! old row = ii  -- new row = iperm(ii) -- ko = new pointer
!
   ko = iao(perm(ii))
   do 60 k=ia(ii), ia(ii+1)-1
      jao(ko) = ja(k)
      if (values) ao(ko) = a(k)
      ko = ko+1
60 continue
100 continue
!
return
!---------end-of-rperm -------------------------------------------------
!-----------------------------------------------------------------------
end
subroutine cperm (nrow,a,ja,ia,ao,jao,iao,perm,job)
integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(*), job
real(kind=kind(0.0d0)) a(*), ao(*)
!-----------------------------------------------------------------------
! this subroutine permutes the columns of a matrix a, ja, ia.
! the result is written in the output matrix  ao, jao, iao.
! cperm computes B = A P, where  P is a permutation matrix
! that maps column j into column perm(j), i.e., on return
!      a(i,j) becomes a(i,perm(j)) in new matrix
! Y. Saad, May 2, 1990 / modified Jan. 28, 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! nrow  = row dimension of the matrix
!
! a, ja, ia = input matrix in csr format.
!
! perm  = integer array of length ncol (number of columns of A
!         containing the permutation array  the columns:
!         a(i,j) in the original matrix becomes a(i,perm(j))
!         in the output matrix.
!
! job   = integer indicating the work to be done:
!       job = 1   permute a, ja, ia into ao, jao, iao
!                       (including the copying of real values ao and
!                       the array iao).
!       job .ne. 1 :  ignore real values ao and ignore iao.
!
!------------
! on return:
!------------
! ao, jao, iao = input matrix in a, ja, ia format (array ao not needed)
!
! Notes:
!-------
! 1. if job=1 then ao, iao are not used.
! 2. This routine is in place: ja, jao can be the same.
! 3. If the matrix is initially sorted (by increasing column number)
!    then ao,jao,iao  may not be on return.
!
!----------------------------------------------------------------------c
! local parameters:
integer k, i, nnz
!
nnz = ia(nrow+1)-1
do 100 k=1,nnz
   jao(k) = perm(ja(k))
100 continue
!
!     done with ja array. return if no need to touch values.
!
if (job .NE. 1) return
!
! else get new pointers -- and copy values too.
!
do 1 i=1, nrow+1
   iao(i) = ia(i)
1 continue
!
do 2 k=1, nnz
   ao(k) = a(k)
2 continue
!
return
!---------end-of-cperm--------------------------------------------------
!-----------------------------------------------------------------------
end
subroutine diapos  (n,ja,ia,idiag)
integer ia(n+1), ja(*), idiag(n)
!-----------------------------------------------------------------------
! this subroutine returns the positions of the diagonal elements of a
! sparse matrix a, ja, ia, in the array idiag.
!-----------------------------------------------------------------------
! on entry:
!----------
!
! n     = integer. row dimension of the matrix a.
! a,ja,
!    ia = matrix stored compressed sparse row format. a array skipped.
!
! on return:
!-----------
! idiag  = integer array of length n. The i-th entry of idiag
!          points to the diagonal element a(i,i) in the arrays
!          a, ja. (i.e., a(idiag(i)) = element A(i,i) of matrix A)
!          if no diagonal element is found the entry is set to 0.
!----------------------------------------------------------------------c
!           Y. Saad, March, 1990
!----------------------------------------------------------------------c
do 1 i=1, n
   idiag(i) = 0
1 continue
!
!     sweep through data structure.
!
do  6 i=1,n
   do 51 k= ia(i),ia(i+1) -1
      if (ja(k) .EQ. i) idiag(i) = k
51 continue
6 continue
!----------- -end-of-diapos---------------------------------------------
!-----------------------------------------------------------------------
return
end
subroutine getbwd(n,a,ja,ia,ml,mu)
!-----------------------------------------------------------------------
! gets the bandwidth of lower part and upper part of A.
! does not assume that A is sorted.
!-----------------------------------------------------------------------
! on entry:
!----------
! n = integer = the row dimension of the matrix
! a, ja,
!    ia = matrix in compressed sparse row format.
!
! on return:
!-----------
! ml   = integer. The bandwidth of the strict lower part of A
! mu   = integer. The bandwidth of the strict upper part of A
!
! Notes:
! ===== ml and mu are allowed to be negative or return. This may be
!       useful since it will tell us whether a band is confined
!       in the strict  upper/lower triangular part.
!       indeed the definitions of ml and mu are
!
!       ml = max ( (i-j)  s.t. a(i,j) .ne. 0  )
!       mu = max ( (j-i)  s.t. a(i,j) .ne. 0  )
!----------------------------------------------------------------------c
! Y. Saad, Sep. 21 1989                                                c
!----------------------------------------------------------------------c
real(kind=kind(0.0d0)) a(*)
integer ja(*),ia(n+1),ml,mu,ldist,i,k
ml = - n
mu = - n
do 3 i=1,n
   do 31 k=ia(i),ia(i+1)-1
      ldist = i-ja(k)
      ml = max(ml,ldist)
      mu = max(mu,-ldist)
31 continue
3 continue
return
!---------------end-of-getbwd ------------------------------------------
!-----------------------------------------------------------------------
end
subroutine infdia (n,ja,ia,ind,idiag)
integer ia(*), ind(*), ja(*)
!-----------------------------------------------------------------------
!     obtains information on the diagonals of A.
!-----------------------------------------------------------------------
! this subroutine finds the lengths of each of the 2*n-1 diagonals of A
! it also outputs the number of nonzero diagonals found.
!-----------------------------------------------------------------------
! on entry:
!----------
! n     = dimension of the matrix a.
!
! a,    ..... not needed here.
! ja,
! ia    = matrix stored in csr format
!
! on return:
!-----------
!
! idiag = integer. number of nonzero diagonals found.
!
! ind   = integer array of length at least 2*n-1. The k-th entry in
!         ind contains the number of nonzero elements in the diagonal
!         number k, the numbering beeing from the lowermost diagonal
!         (bottom-left). In other words ind(k) = length of diagonal
!         whose offset wrt the main diagonal is = - n + k.
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
n2= n+n-1
do 1 i=1,n2
   ind(i) = 0
1 continue
do 3 i=1, n
   do 2 k=ia(i),ia(i+1)-1
      j = ja(k)
      ind(n+j-i) = ind(n+j-i) +1
2 continue
3 continue
!     count the nonzero ones.
idiag = 0
do 41 k=1, n2
   if (ind(k) .NE. 0) idiag = idiag+1
41 continue
return
! done
!------end-of-infdia ---------------------------------------------------
!-----------------------------------------------------------------------
end
subroutine rnrms   (nrow, nrm, a, ja, ia, diag)
real(kind=kind(0.0d0)) a(*), diag(nrow), scal
integer ja(*), ia(nrow+1)
!-----------------------------------------------------------------------
! gets the norms of each row of A. (choice of three norms)
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow  = integer. The row dimension of A
!
! nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
!                  means the 2-nrm, nrm = 0 means max norm
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!
! on return:
!----------
!
! diag = real vector of length nrow containing the norms
!
!-----------------------------------------------------------------
do 1 ii=1,nrow
!
!     compute the norm if each element.
!
   scal = 0.0d0
   k1 = ia(ii)
   k2 = ia(ii+1)-1
   if (nrm .EQ. 0) then
      do 2 k=k1, k2
         scal = max(scal,abs(a(k) ) )
2 continue
   elseif (nrm .EQ. 1) then
      do 3 k=k1, k2
         scal = scal + abs(a(k) )
3 continue
   else
      do 4 k=k1, k2
         scal = scal+a(k)**2
4 continue
   endif
   if (nrm .EQ. 2) scal = sqrt(scal)
   diag(ii) = scal
1 continue
return
!-----------------------------------------------------------------------
!-------------end-of-rnrms----------------------------------------------
end
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!                   ITERATIVE SOLVERS MODULE                           c
!----------------------------------------------------------------------c
! ILUT    : Incomplete LU factorization with dual truncation strategy  c
! ILUTP   : ILUT with column  pivoting                                 c
! LUSOL   : forward followed by backward triangular solve (Precond.)   c
! QSPLIT  : quick split routine used by ilut to sort out the k largest c
!           elements in absolute value                                 c
!----------------------------------------------------------------------c
  subroutine qsplit(a,ind,n,ncut)
  real(kind=kind(0.0d0)) a(n)
  integer ind(n), n, ncut
!-----------------------------------------------------------------------
!     does a quick-sort split of a real array.
!     on input a(1:n). is a real array
!     on output a(1:n) is permuted such that its elements satisfy:
!
!     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
!     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
!
!     ind(1:n) is an integer array which permuted in the same way as a(*).
!-----------------------------------------------------------------------
  real(kind=kind(0.0d0)) tmp, abskey
  integer itmp, first, last
!-----
  first = 1
  last = n
  if (ncut .LT. first .OR. ncut .GT. last) return
!
!     outer loop -- while mid .ne. ncut do
!
1 mid = first
  abskey = abs(a(mid))
  do 2 j=first+1, last
     if (abs(a(j)) .GT. abskey) then
        mid = mid+1
!     interchange
        tmp = a(mid)
        itmp = ind(mid)
        a(mid) = a(j)
        ind(mid) = ind(j)
        a(j)  = tmp
        ind(j) = itmp
     endif
2 continue
!
!     interchange
!
  tmp = a(mid)
  a(mid) = a(first)
  a(first)  = tmp
!
  itmp = ind(mid)
  ind(mid) = ind(first)
  ind(first) = itmp
!
!     test for while loop
!
  if (mid .EQ. ncut) return
  if (mid .GT. ncut) then
     last = mid-1
  else
     first = mid+1
  endif
  goto 1
!----------------end-of-qsplit------------------------------------------
!-----------------------------------------------------------------------
  end
!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!               REORDERING ROUTINES -- LEVEL SET BASED ROUTINES        c
!----------------------------------------------------------------------c
! dblstr   : doubled stripe partitioner
! BFS      : Breadth-First search traversal algorithm
! add_lvst : routine to add a level -- used by BFS
! stripes  : finds the level set structure
! perphn   : finds a pseudo-peripheral node and performs a BFS from it.
! rversp   : routine to reverse a given permutation (e.g., for RCMK)
! maskdeg  : integer function to compute the `masked' of a node
!-----------------------------------------------------------------------
subroutine BFS(n,ja,ia,nfirst,iperm,mask,maskval,riord,levels, &
& nlev)
implicit none
integer n,ja(*),ia(*),nfirst,iperm(n),mask(n),riord(*),levels(*), &
& nlev,maskval
!-----------------------------------------------------------------------
! finds the level-structure (breadth-first-search or CMK) ordering for a
! given sparse matrix. Uses add_lvst. Allows an set of nodes to be
! the initial level (instead of just one node).
!-------------------------parameters------------------------------------
! on entry:
!---------
!     n      = number of nodes in the graph
!     ja, ia = pattern of matrix in CSR format (the ja,ia arrays of csr data
!     structure)
!     nfirst = number of nodes in the first level that is input in riord
!     iperm  = integer array indicating in which order to  traverse the graph
!     in order to generate all connected components.
!     if iperm(1) .eq. 0 on entry then BFS will traverse the nodes
!     in the  order 1,2,...,n.
!
!     riord  = (also an ouput argument). On entry riord contains the labels
!     of the nfirst nodes that constitute the first level.
!
!     mask   = array used to indicate whether or not a node should be
!     condidered in the graph. see maskval.
!     mask is also used as a marker of  visited nodes.
!
!     maskval= consider node i only when:  mask(i) .eq. maskval
!     maskval must be .gt. 0.
!     thus, to consider all nodes, take mask(1:n) = 1.
!     maskval=1 (for example)
!
!     on return
!     ---------
!     mask   = on return mask is restored to its initial state.
!     riord  = `reverse permutation array'. Contains the labels of the nodes
!     constituting all the levels found, from the first level to
!     the last.
!     levels = pointer array for the level structure. If lev is a level
!     number, and k1=levels(lev),k2=levels(lev+1)-1, then
!     all the nodes of level number lev are:
!     riord(k1),riord(k1+1),...,riord(k2)
!     nlev   = number of levels found
!-----------------------------------------------------------------------
!
integer j, ii, nod, istart, iend
logical permut
permut = (iperm(1) .NE. 0)
!
!     start pointer structure to levels
!
nlev   = 0
!
!     previous end
!
istart = 0
ii = 0
!
!     current end
!
iend = nfirst
!
!     intialize masks to zero -- except nodes of first level --
!
do 12 j=1, nfirst
   mask(riord(j)) = 0
12 continue
!-----------------------------------------------------------------------
continue
!
1 nlev = nlev+1
levels(nlev) = istart + 1
call add_lvst (istart,iend,nlev,riord,ja,ia,mask,maskval)
if (istart .LT. iend) goto 1
2 ii = ii+1
if (ii .LE. n) then
   nod = ii
   if (permut) nod = iperm(nod)
   if (mask(nod) .EQ. maskval) then
!
!     start a new level
!
      istart = iend
      iend = iend+1
      riord(iend) = nod
      mask(nod) = 0
      goto 1
   else
      goto 2
   endif
endif
!-----------------------------------------------------------------------
levels(nlev+1) = iend+1
do j=1, iend
   mask(riord(j)) = maskval
enddo
!-----------------------------------------------------------------------
return
end
subroutine add_lvst(istart,iend,nlev,riord,ja,ia,mask,maskval)
integer nlev, nod, riord(*), ja(*), ia(*), mask(*)
!-------------------------------------------------------------
!     adds one level set to the previous sets..
!     span all nodes of previous mask
!-------------------------------------------------------------
nod = iend
do 25 ir = istart+1,iend
   i = riord(ir)
   do 24 k=ia(i),ia(i+1)-1
      j = ja(k)
      if (mask(j) .EQ. maskval) then
         nod = nod+1
         mask(j) = 0
         riord(nod) = j
      endif
24 continue
25 continue
istart = iend
iend   = nod
return
end
subroutine stripes (nlev,riord,levels,ip,map,mapptr,ndom)
implicit none
integer nlev,riord(*),levels(nlev+1),ip,map(*), &
& mapptr(*), ndom
!-----------------------------------------------------------------------
!    this is a post processor to BFS. stripes uses the output of BFS to
!    find a decomposition of the adjacency graph by stripes. It fills
!    the stripes level by level until a number of nodes .gt. ip is
!    is reached.
!---------------------------parameters-----------------------------------
! on entry:
! --------
! nlev   = number of levels as found by BFS
! riord  = reverse permutation array produced by BFS --
! levels = pointer array for the level structure as computed by BFS. If
!          lev is a level number, and k1=levels(lev),k2=levels(lev+1)-1,
!          then all the nodes of level number lev are:
!                      riord(k1),riord(k1+1),...,riord(k2)
!  ip    = number of desired partitions (subdomains) of about equal size.
!
! on return
! ---------
! ndom     = number of subgraphs (subdomains) found
! map      = node per processor list. The nodes are listed contiguously
!            from proc 1 to nproc = mpx*mpy.
! mapptr   = pointer array for array map. list for proc. i starts at
!            mapptr(i) and ends at mapptr(i+1)-1 in array map.
!-----------------------------------------------------------------------
! local variables.
!
integer ib,ktr,ilev,k,nsiz,psiz
ndom = 1
ib = 1
! to add: if (ip .le. 1) then ...
nsiz = levels(nlev+1) - levels(1)
psiz = (nsiz-ib)/max(1,(ip - ndom + 1)) + 1
mapptr(ndom) = ib
ktr = 0
do 10 ilev = 1, nlev
!
!     add all nodes of this level to domain
!
   do 3 k=levels(ilev), levels(ilev+1)-1
      map(ib) = riord(k)
      ib = ib+1
      ktr = ktr + 1
      if (ktr .GE. psiz  .OR. k .GE. nsiz) then
         ndom = ndom + 1
         mapptr(ndom) = ib
         psiz = (nsiz-ib)/max(1,(ip - ndom + 1)) + 1
         ktr = 0
      endif
!
3 continue
10 continue
ndom = ndom-1
return
end
integer function maskdeg (ja,ia,nod,mask,maskval)
implicit none
integer ja(*),ia(*),nod,mask(*),maskval
!-----------------------------------------------------------------------
integer deg, k
deg = 0
do k =ia(nod),ia(nod+1)-1
   if (mask(ja(k)) .EQ. maskval) deg = deg+1
enddo
maskdeg = deg
return
end
subroutine perphn(n,ja,ia,init,mask,maskval,nlev,riord,levels)
implicit none
integer n,ja(*),ia(*),init,mask(*),maskval, &
& nlev,riord(*),levels(*)
!-----------------------------------------------------------------------
!     finds a peripheral node and does a BFS search from it.
!-----------------------------------------------------------------------
!     see routine  dblstr for description of parameters
! input:
!-------
! ja, ia  = list pointer array for the adjacency graph
! mask    = array used for masking nodes -- see maskval
! maskval = value to be checked against for determing whether or
!           not a node is masked. If mask(k) .ne. maskval then
!           node k is not considered.
! init    = init node in the pseudo-peripheral node algorithm.
!
! output:
!-------
! init    = actual pseudo-peripherial node found.
! nlev    = number of levels in the final BFS traversal.
! riord   =
! levels  =
!-----------------------------------------------------------------------
integer j,nlevp,deg,nfirst,mindeg,nod,maskdeg
integer iperm(1)
nlevp = 0
1 continue
riord(1) = init
nfirst = 1
iperm(1) = 0
!
call BFS(n,ja,ia,nfirst,iperm,mask,maskval,riord,levels,nlev)
if (nlev .GT. nlevp) then
   mindeg = n+1
   do j=levels(nlev),levels(nlev+1)-1
      nod = riord(j)
      deg = maskdeg(ja,ia,nod,mask,maskval)
      if (deg .LT. mindeg) then
         init = nod
         mindeg = deg
      endif
   enddo
   nlevp = nlev
   goto 1
endif
return
end
!----------------------------------------------------------------------c
!     Non-SPARSKIT utility routine
!----------------------------------------------------------------------c
