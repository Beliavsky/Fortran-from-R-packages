!  The code in this file is was taken from daspk.tgz from
!    https://www.netlib.org/ode/
!  Author: Cleve Moler, University of New Mexico, Argonne National Lab.
!     
!  Adapted for use in R package deSolve by the deSolve authors.
!

subroutine dgefa(a,lda,n,ipvt,info)
integer lda,n,ipvt(*),info
double precision a(lda,*)
!
!     dgefa factors a double precision matrix by gaussian elimination.
!
!     dgefa is usually called by dgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
!
!     on entry
!
!        a       double precision(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that dgesl or dgedi will divide by zero
!                     if called.  use  rcond  in dgeco for a reliable
!                     indication of singularity.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,dscal,idamax
!
!     internal variables
!
double precision t
integer idamax,j,k,kp1,l,nm1
!
!
!     gaussian elimination with partial pivoting
!
info = 0
nm1 = n - 1
if (nm1 .LT. 1) go to 70
do 60 k = 1, nm1
   kp1 = k + 1
!
!        find l = pivot index
!
   l = idamax(n-k+1,a(k,k),1) + k - 1
   ipvt(k) = l
!
!        zero pivot implies this column already triangularized
!
   if (a(l,k) .EQ. 0.0d0) go to 40
!
!           interchange if necessary
!
      if (l .EQ. k) go to 10
         t = a(l,k)
         a(l,k) = a(k,k)
         a(k,k) = t
10 continue
!
!           compute multipliers
!
      t = -1.0d0/a(k,k)
      call dscal(n-k,t,a(k+1,k),1)
!
!           row elimination with column indexing
!
      do 30 j = kp1, n
         t = a(l,j)
         if (l .EQ. k) go to 20
            a(l,j) = a(k,j)
            a(k,j) = t
20 continue
         call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
30 continue
   go to 50
40 continue
      info = k
50 continue
60 continue
70 continue
ipvt(n) = n
if (a(n,n) .EQ. 0.0d0) info = n
return
end
subroutine dgesl(a,lda,n,ipvt,b,job)
integer lda,n,ipvt(*),job
double precision a(lda,*),b(*)
!
!     dgesl solves the double precision system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by dgeco or dgefa.
!
!     on entry
!
!        a       double precision(lda, n)
!                the output from dgeco or dgefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        ipvt    integer(n)
!                the pivot vector from dgeco or dgefa.
!
!        b       double precision(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b  where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if dgeco has set rcond .gt. 0.0
!        or dgefa has set info .eq. 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call dgeco(a,lda,n,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call dgesl(a,lda,n,ipvt,c(1,j),0)
!        10 continue
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,ddot
!
!     internal variables
!
double precision ddot,t
integer k,kb,l,nm1
!
nm1 = n - 1
if (job .NE. 0) go to 50
!
!        job = 0 , solve  a * x = b
!        first solve  l*y = b
!
   if (nm1 .LT. 1) go to 30
   do 20 k = 1, nm1
      l = ipvt(k)
      t = b(l)
      if (l .EQ. k) go to 10
         b(l) = b(k)
         b(k) = t
10 continue
      call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
20 continue
30 continue
!
!        now solve  u*x = y
!
   do 40 kb = 1, n
      k = n + 1 - kb
      b(k) = b(k)/a(k,k)
      t = -b(k)
      call daxpy(k-1,t,a(1,k),1,b(1),1)
40 continue
go to 100
50 continue
!
!        job = nonzero, solve  trans(a) * x = b
!        first solve  trans(u)*y = b
!
   do 60 k = 1, n
      t = ddot(k-1,a(1,k),1,b(1),1)
      b(k) = (b(k) - t)/a(k,k)
60 continue
!
!        now solve trans(l)*x = y
!
   if (nm1 .LT. 1) go to 90
   do 80 kb = 1, nm1
      k = n - kb
      b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
      l = ipvt(k)
      if (l .EQ. k) go to 70
         t = b(l)
         b(l) = b(k)
         b(k) = t
70 continue
80 continue
90 continue
100 continue
return
end
subroutine dgbfa(abd,lda,n,ml,mu,ipvt,info)
integer lda,n,ml,mu,ipvt(*),info
double precision abd(lda,*)
!
!     dgbfa factors a double precision band matrix by elimination.
!
!     dgbfa is usually called by dgbco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!
!     on entry
!
!        abd     double precision(lda, n)
!                contains the matrix in band storage.  the columns
!                of the matrix are stored in the columns of  abd  and
!                the diagonals of the matrix are stored in rows
!                ml+1 through 2*ml+mu+1 of  abd .
!                see the comments below for details.
!
!        lda     integer
!                the leading dimension of the array  abd .
!                lda must be .ge. 2*ml + mu + 1 .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!                0 .le. ml .lt. n .
!
!        mu      integer
!                number of diagonals above the main diagonal.
!                0 .le. mu .lt. n .
!                more efficient if  ml .le. mu .
!     on return
!
!        abd     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that dgbsl will divide by zero if
!                     called.  use  rcond  in dgbco for a reliable
!                     indication of singularity.
!
!     band storage
!
!           if  a  is a band matrix, the following program segment
!           will set up the input.
!
!                   ml = (band width below the diagonal)
!                   mu = (band width above the diagonal)
!                   m = ml + mu + 1
!                   do 20 j = 1, n
!                      i1 = max0(1, j-mu)
!                      i2 = min0(n, j+ml)
!                      do 10 i = i1, i2
!                         k = i - j + m
!                         abd(k,j) = a(i,j)
!                10    continue
!                20 continue
!
!           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
!           in addition, the first  ml  rows in  abd  are used for
!           elements generated during the triangularization.
!           the total number of rows needed in  abd  is  2*ml+mu+1 .
!           the  ml+mu by ml+mu  upper left triangle and the
!           ml by ml  lower right triangle are not referenced.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,dscal,idamax
!     fortran max0,min0
!
!     internal variables
!
double precision t
integer i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
!
!
m = ml + mu + 1
info = 0
!
!     zero initial fill-in columns
!
j0 = mu + 2
j1 = min0(n,m) - 1
if (j1 .LT. j0) go to 30
do 20 jz = j0, j1
   i0 = m + 1 - jz
   do 10 i = i0, ml
      abd(i,jz) = 0.0d0
10 continue
20 continue
30 continue
jz = j1
ju = 0
!
!     gaussian elimination with partial pivoting
!
nm1 = n - 1
if (nm1 .LT. 1) go to 130
do 120 k = 1, nm1
   kp1 = k + 1
!
!        zero next fill-in column
!
   jz = jz + 1
   if (jz .GT. n) go to 50
   if (ml .LT. 1) go to 50
      do 40 i = 1, ml
         abd(i,jz) = 0.0d0
40 continue
50 continue
!
!        find l = pivot index
!
   lm = min0(ml,n-k)
   l = idamax(lm+1,abd(m,k),1) + m - 1
   ipvt(k) = l + k - m
!
!        zero pivot implies this column already triangularized
!
   if (abd(l,k) .EQ. 0.0d0) go to 100
!
!           interchange if necessary
!
      if (l .EQ. m) go to 60
         t = abd(l,k)
         abd(l,k) = abd(m,k)
         abd(m,k) = t
60 continue
!
!           compute multipliers
!
      t = -1.0d0/abd(m,k)
      call dscal(lm,t,abd(m+1,k),1)
!
!           row elimination with column indexing
!
      ju = min0(max0(ju,mu+ipvt(k)),n)
      mm = m
      if (ju .LT. kp1) go to 90
      do 80 j = kp1, ju
         l = l - 1
         mm = mm - 1
         t = abd(l,j)
         if (l .EQ. mm) go to 70
            abd(l,j) = abd(mm,j)
            abd(mm,j) = t
70 continue
         call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
80 continue
90 continue
   go to 110
100 continue
      info = k
110 continue
120 continue
130 continue
ipvt(n) = n
if (abd(m,n) .EQ. 0.0d0) info = n
return
end
subroutine dgbsl(abd,lda,n,ml,mu,ipvt,b,job)
integer lda,n,ml,mu,ipvt(*),job
double precision abd(lda,*),b(*)
!
!     dgbsl solves the double precision band system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by dgbco or dgbfa.
!
!     on entry
!
!        abd     double precision(lda, n)
!                the output from dgbco or dgbfa.
!
!        lda     integer
!                the leading dimension of the array  abd .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!
!        mu      integer
!                number of diagonals above the main diagonal.
!
!        ipvt    integer(n)
!                the pivot vector from dgbco or dgbfa.
!
!        b       double precision(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b , where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if dgbco has set rcond .gt. 0.0
!        or dgbfa has set info .eq. 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
!        10 continue
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,ddot
!     fortran min0
!
!     internal variables
!
double precision ddot,t
integer k,kb,l,la,lb,lm,m,nm1
!
m = mu + ml + 1
nm1 = n - 1
if (job .NE. 0) go to 50
!
!        job = 0 , solve  a * x = b
!        first solve l*y = b
!
   if (ml .EQ. 0) go to 30
   if (nm1 .LT. 1) go to 30
      do 20 k = 1, nm1
         lm = min0(ml,n-k)
         l = ipvt(k)
         t = b(l)
         if (l .EQ. k) go to 10
            b(l) = b(k)
            b(k) = t
10 continue
         call daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
20 continue
30 continue
!
!        now solve  u*x = y
!
   do 40 kb = 1, n
      k = n + 1 - kb
      b(k) = b(k)/abd(m,k)
      lm = min0(k,m) - 1
      la = m - lm
      lb = k - lm
      t = -b(k)
      call daxpy(lm,t,abd(la,k),1,b(lb),1)
40 continue
go to 100
50 continue
!
!        job = nonzero, solve  trans(a) * x = b
!        first solve  trans(u)*y = b
!
   do 60 k = 1, n
      lm = min0(k,m) - 1
      la = m - lm
      lb = k - lm
      t = ddot(lm,abd(la,k),1,b(lb),1)
      b(k) = (b(k) - t)/abd(m,k)
60 continue
!
!        now solve trans(l)*x = y
!
   if (ml .EQ. 0) go to 90
   if (nm1 .LT. 1) go to 90
      do 80 kb = 1, nm1
         k = n - kb
         lm = min0(ml,n-k)
         b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
         l = ipvt(k)
         if (l .EQ. k) go to 70
            t = b(l)
            b(l) = b(k)
            b(k) = t
70 continue
80 continue
90 continue
100 continue
return
end
