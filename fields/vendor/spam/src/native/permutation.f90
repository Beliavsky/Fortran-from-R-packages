!     It has been tested that embedding the loop over the right hand
!     side into the backsolve routine is not faster.

!-----------------------------------------------------------------------
subroutine getbwd(n,  ja,ia,ml,mu)
!-----------------------------------------------------------------------
! gets the bandwidth of lower part and upper part of A.
! does not assume that A is sorted.
!-----------------------------------------------------------------------
! on entry:
!----------
! n	= integer = the row dimension of the matrix
! ja, ia = matrix in compressed sparse row format.
!
! on return:
!-----------
! ml	= integer. The bandwidth of the strict lower part of A
! mu	= integer. The bandwidth of the strict upper part of A
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
! Y. Saad, Sep. 21 1989 | simplified by R. Furrer Sept. 2016           c
!----------------------------------------------------------------------c
implicit none
!     double precision a(*)
integer n,ja(*),ia(n+1),ml,mu

integer ldist,i,k

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


! functions slightly modified from sparsekit:
! cperm,rperm,dperm: job argument is eliminated
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine rperm (nrow,a,ja,ia,ao,jao,iao,perm)
implicit none
integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow)
double precision a(*),ao(*)
!-----------------------------------------------------------------------
! this subroutine permutes the rows of a matrix in CSR format.
! rperm  computes B = P A  where P is a permutation matrix.
! the permutation P is defined through the array perm: for each j,
! perm(j) represents the destination row number of row number j.
! Youcef Saad -- recoded Jan 28, 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! n 	= dimension of the matrix
! a, ja, ia = input matrix in csr format
! perm 	= integer array of length nrow containing the permutation arrays
!	  for the rows: perm(i) is the destination of row i in the
!         permuted matrix.
!         ---> a(i,j) in the original matrix becomes a(perm(i),j)
!         in the output  matrix.
!
!
!------------
! on return:
!------------
! ao, jao, iao = input matrix in a, ja, ia format
!----------------------------------------------------------------------c
!           Y. Saad, May  2, 1990                                      c
!----------------------------------------------------------------------c

integer i,j,k,ko,ii
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
ao(ko) = a(k)
ko = ko+1
60 continue
100 continue
!
return
!---------end-of-rperm -------------------------------------------------
!-----------------------------------------------------------------------
end
!-----------------------------------------------------------------------
subroutine cperm (nrow,a,ja,ia,ao,jao,iao,perm)
implicit none
integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(*)
double precision a(*), ao(*)
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
! nrow 	= row dimension of the matrix
!
! a, ja, ia = input matrix in csr format.
!
! perm	= integer array of length ncol (number of columns of A
!         containing the permutation array  the columns:
!         a(i,j) in the original matrix becomes a(i,perm(j))
!         in the output matrix.
!
!
!------------
! on return:
!------------
! ao, jao, iao = input matrix in a, ja, ia format
!
! Notes:
!-------
! 1. if job=1 then ao, iao are not used.
! 2. This routine is in place: ja, jao can be the same.
! 3. If the matrix is initially sorted (by increasing column number)
!    then ao,jao,iao  may not be on return, hence a call to csort.
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
!     done with ja array.
!
do 1 i=1, nrow+1
iao(i) = ia(i)
1 continue
!
do 2 k=1, nnz
ao(k) = a(k)
2 continue
!
call sortrows(nrow,ao,jao,iao)
!     call csort (nrow,ao,jao,iao,iwork) _does not work_
return
!---------end-of-cperm--------------------------------------------------
!-----------------------------------------------------------------------
end
!-----------------------------------------------------------------------
subroutine dperm (nrow,a,ja,ia,ao,jao,iao,pperm,qperm)
implicit none
integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),pperm(nrow), &
& qperm(*)
double precision a(*),ao(*)
!-----------------------------------------------------------------------
! This routine permutes the rows and columns of a matrix stored in CSR
! format. i.e., it computes P A Q, where P, Q are permutation matrices.
! P maps row i into row perm(i) and Q maps column j into column qperm(j):
!      a(i,j)    becomes   a(pperm(i),qperm(j)) in new matrix
! note that qperm should be of length ncol (number of columns) but this
! is not checked.
!-----------------------------------------------------------------------
! Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991.
!-----------------------------------------------------------------------
! on entry:
!----------
! n 	= dimension of the matrix
! a, ja,
!    ia = input matrix in a, ja, ia format
! pperm 	= integer array of length n containing the permutation arrays
!	  for the rows: pperm(i) is the destination of row i in the
!         permuted matrix
!
! qperm	= same thing for the columns.
! iwork = working array passed to cperm
!
! on return:
!-----------
! ao, jao, iao = input matrix in a, ja, ia format
!
! Notes:
!-------
!  1) algorithm is in place
!----------------------------------------------------------------------c
! local variables
!
! permute rows first
!
call rperm (nrow,a,ja,ia,   ao,jao,iao,pperm)
!
! then permute columns
!
!
call cperm (nrow,ao,jao,iao,ao,jao,iao,qperm)
!
return
!-------end-of-dperm----------------------------------------------------
!-----------------------------------------------------------------------
end



!-----------------------------------------------------------------------
subroutine dvperm (n, x, perm)
implicit none
integer n
integer perm(n)
double precision x(n)
!-----------------------------------------------------------------------
! this subroutine performs an in-place permutation of a real vector x
! according to the permutation array perm(*), i.e., on return,
! the vector x satisfies,
!
!	x(perm(j)) :== x(j), j=1,2,.., n
!
!-----------------------------------------------------------------------
! on entry:
!---------
! n 	= length of vector x.
! perm 	= integer array of length n containing the permutation  array.
! x	= input vector
!
! on return:
!----------
! x	= vector x permuted according to x(perm(*)) :=  x(*)
!
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
! local variables
integer init,next,k, ii, j
double precision tmp, tmp1
!
init      = 1
tmp       = x(init)
ii        = perm(init)
perm(init)= -perm(init)
k         = 0
!
! loop
!
6 k = k+1
!
! save the chased element --
!
tmp1      = x(ii)
x(ii)     = tmp
next      = perm(ii)
if (next .lt. 0 ) goto 65
!
! test for end
!
if (k .gt. n) goto 101
tmp       = tmp1
perm(ii)  = - perm(ii)
ii        = next
!
! end loop
!
goto 6
!
! reinitilaize cycle --
!
65 init      = init+1
if (init .gt. n) goto 101
if (perm(init) .lt. 0) goto 65
tmp = x(init)
ii = perm(init)
perm(init)=-perm(init)
goto 6
!
101 continue
do 200 j=1, n
perm(j) = -perm(j)
200 continue
!
return
!-------------------end-of-dvperm---------------------------------------
!-----------------------------------------------------------------------
end
!-----------------------------------------------------------------------
subroutine ivperm (n, ix, perm)
implicit none
integer n, perm(n+1), ix(n)
!-----------------------------------------------------------------------
! this subroutine performs an in-place permutation of an integer vector
! ix according to the permutation array perm(*), i.e., on return,
! the vector x satisfies,
!
!	ix(perm(j)) :== ix(j), j=1,2,.., n
!
!-----------------------------------------------------------------------
! on entry:
!---------
! n 	= length of vector x.
! perm 	= integer array of length n containing the permutation  array.
! ix	= input vector
!
! on return:
!----------
! ix	= vector x permuted according to ix(perm(*)) :=  ix(*)
!
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
! local variables
integer ii,k,j,next,init,tmp, tmp1
!
init      = 1
tmp       = ix(init)
ii        = perm(init)
perm(init)= -perm(init)
k         = 0
!
! loop
!
6 k = k+1
!
! save the chased element --
!
tmp1     = ix(ii)
ix(ii)   = tmp
next     = perm(ii)
if (next .lt. 0 ) goto 65
!
! test for end
!
if (k .gt. n) goto 101
tmp       = tmp1
perm(ii)  = - perm(ii)
ii        = next
!
! end loop
!
goto 6
!
! reinitilaize cycle --
!
65 init      = init+1
if (init .gt. n) goto 101
if (perm(init) .lt. 0) goto 65
tmp       = ix(init)
ii        = perm(init)
perm(init)=-perm(init)
goto 6
!
101 continue
do 200 j=1, n
!         if (perm(j) .lt. 0) then
perm(j) = -perm(j)
!         endif
200 continue
!
return
!-------------------end-of-ivperm---------------------------------------
!-----------------------------------------------------------------------
end
!
!-----------------------------------------------------------------------
subroutine aplbdg (nrow,ncol,ja,ia,jb,ib,ndegr,nnz,iw)

implicit none
integer nrow, ncol, nnz
integer ja(*),jb(*),ia(nrow+1),ib(nrow+1),iw(ncol),ndegr(nrow)
!-----------------------------------------------------------------------
! gets the number of nonzero elements in each row of A+B and the total
! number of nonzero elements in A+B.
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow	= integer. The row dimension of A and B
! ncol  = integer. The column dimension of A and B.
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!
! b,
! jb,
! ib	=  Matrix B in compressed sparse row format.
!
! iw,nnz,ngegr = zero content
!
! on return:
!----------
! ndegr	= integer array of length nrow containing the degrees (i.e.,
!         the number of nonzeros in  each row of the matrix A + B.
!
! nnz   = total number of nonzero elements found in A * B
!
! work arrays:
!------------
! iw	= integer work array of length equal to ncol.
!
!-----------------------------------------------------------------------
integer k,j,ii,jr,last,ldg,jc

do 7 ii=1,nrow
ldg = 0
!
!    end-of-linked list
!
last = -1
!
!     row of A
!
do 5 j = ia(ii),ia(ii+1)-1
jr = ja(j)
!
!     add element to the linked list
!
ldg = ldg + 1
iw(jr) = last
last = jr
5 continue
!
!     row of B
!
do 6 j=ib(ii),ib(ii+1)-1
jc = jb(j)
if (iw(jc) .eq. 0) then
!
!     add one element to the linked list
!
ldg = ldg + 1
iw(jc) = last
last = jc
endif
6 continue
!     done with row ii.
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

do 8 ii=1, nrow
nnz = nnz+ndegr(ii)
8 continue
return
!----------------end-of-aplbdg -----------------------------------------
!-----------------------------------------------------------------------
end
