!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dstqrb
!
!\Description:
!  Computes all eigenvalues and the last component of the eigenvectors
!  of a symmetric tridiagonal matrix using the implicit QL or QR method.
!
!  This is mostly a modification of the LAPACK routine dsteqr.
!  See Remarks.
!
!\Usage:
!  call dstqrb
!     ( N, D, E, Z, WORK, INFO )
!
!\Arguments
!  N       Integer.  (INPUT)
!          The number of rows and columns in the matrix.  N >= 0.
!
!  D       Double precision array, dimension (N).  (INPUT/OUTPUT)
!          On entry, D contains the diagonal elements of the
!          tridiagonal matrix.
!          On exit, D contains the eigenvalues, in ascending order.
!          If an error exit is made, the eigenvalues are correct
!          for indices 1,2,...,INFO-1, but they are unordered and
!          may not be the smallest eigenvalues of the matrix.
!
!  E       Double precision array, dimension (N-1).  (INPUT/OUTPUT)
!          On entry, E contains the subdiagonal elements of the
!          tridiagonal matrix in positions 1 through N-1.
!          On exit, E has been destroyed.
!
!  Z       Double precision array, dimension (N).  (OUTPUT)
!          On exit, Z contains the last row of the orthonormal
!          eigenvector matrix of the symmetric tridiagonal matrix.
!          If an error exit is made, Z contains the last row of the
!          eigenvector matrix associated with the stored eigenvalues.
!
!  WORK    Double precision array, dimension (max(1,2*N-2)).  (WORKSPACE)
!          Workspace used in accumulating the transformation for
!          computing the last components of the eigenvectors.
!
!  INFO    Integer.  (OUTPUT)
!          = 0:  normal return.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = +i, the i-th eigenvalue has not converged
!                              after a total of  30*N  iterations.
!
!\Remarks
!  1. None.
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     daxpy   Level 1 BLAS that computes a vector triad.
!     dcopy   Level 1 BLAS that copies one vector to another.
!     dswap   Level 1 BLAS that swaps the contents of two vectors.
!     lsame   LAPACK character comparison routine.
!     dlae2   LAPACK routine that computes the eigenvalues of a 2-by-2
!             symmetric matrix.
!     dlaev2  LAPACK routine that eigendecomposition of a 2-by-2 symmetric
!             matrix.
!     dlamch  LAPACK routine that determines machine constants.
!     dlanst  LAPACK routine that computes the norm of a matrix.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     dlartg  LAPACK Givens rotation construction routine.
!     dlascl  LAPACK routine for careful scaling of a matrix.
!     dlaset  LAPACK matrix initialization routine.
!     dlasr   LAPACK routine that applies an orthogonal transformation to
!             a matrix.
!     dlasrt  LAPACK sorting routine.
!     dsteqr  LAPACK routine that computes eigenvalues and eigenvectors
!             of a symmetric tridiagonal matrix.
!     xerbla  LAPACK error handler routine.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: stqrb.F   SID: 2.5   DATE OF SID: 8/27/96   RELEASE: 2
!
!\Remarks
!     1. Starting with version 2.5, this routine is a modified version
!        of LAPACK version 2.0 subroutine SSTEQR. No lines are deleted,
!        only commented out and new lines inserted.
!        All lines commented out have "c$$$" at the beginning.
!        Note that the LAPACK version 1.0 subroutine SSTEQR contained
!        bugs.
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dstqrb ( n, d, e, z, work, info )
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
integer    info, n
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
Double precision &
  &           d( n ), e( n-1 ), z( n ), work( 2*n-2 )
!
!     .. parameters ..
Double precision &
  &                   zero, one, two, three
parameter          ( zero = 0.0D+0, one = 1.0D+0, &
  &                     two = 2.0D+0, three = 3.0D+0 )
integer            maxit
parameter          ( maxit = 30 )
!     ..
!     .. local scalars ..
integer            i, icompz, ii, iscale, j, jtot, k, l, l1, lend, &
  &                   lendm1, lendp1, lendsv, lm1, lsv, m, mm, mm1, &
  &                   nm1, nmaxit
Double precision &
  &                   anorm, b, c, eps, eps2, f, g, p, r, rt1, rt2, &
  &                   s, safmax, safmin, ssfmax, ssfmin, tst
!     ..
!     .. external functions ..
logical            lsame
Double precision &
  &                   dlamch, dlanst, dlapy2
external           lsame, dlamch, dlanst, dlapy2
!     ..
!     .. external subroutines ..
external           dlae2, dlaev2, dlartg, dlascl, dlaset, dlasr, &
  &                   dlasrt, dswap, xerbla
!     ..
!     .. intrinsic functions ..
intrinsic          abs, max, sign, sqrt
!     ..
!     .. executable statements ..
!
!     test the input parameters.
!
info = 0
!
!$$$      IF( LSAME( COMPZ, 'N' ) ) THEN
!$$$         ICOMPZ = 0
!$$$      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
!$$$         ICOMPZ = 1
!$$$      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
!$$$         ICOMPZ = 2
!$$$      ELSE
!$$$         ICOMPZ = -1
!$$$      END IF
!$$$      IF( ICOMPZ.LT.0 ) THEN
!$$$         INFO = -1
!$$$      ELSE IF( N.LT.0 ) THEN
!$$$         INFO = -2
!$$$      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
!$$$     $         N ) ) ) THEN
!$$$         INFO = -6
!$$$      END IF
!$$$      IF( INFO.NE.0 ) THEN
!$$$         CALL XERBLA( 'SSTEQR', -INFO )
!$$$         RETURN
!$$$      END IF
!
!    *** New starting with version 2.5 ***
!
icompz = 2
!    *************************************
!
!     quick return if possible
!
if( n.eq.0 ) &
  &   return
!
if( n.eq.1 ) then
   if( icompz.eq.2 )  z( 1 ) = one
   return
end if
!
!     determine the unit roundoff and over/underflow thresholds.
!
eps = dlamch( 'e' )
eps2 = eps**2
safmin = dlamch( 's' )
safmax = one / safmin
ssfmax = sqrt( safmax ) / three
ssfmin = sqrt( safmin ) / eps2
!
!     compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
!$$      if( icompz.eq.2 )
!$$$     $   call dlaset( 'full', n, n, zero, one, z, ldz )
!
!     *** New starting with version 2.5 ***
!
if ( icompz .eq.2 ) then
   do 5 j = 1, n-1
      z(j) = zero
5    continue
   z( n ) = one
end if
!     *************************************
!
nmaxit = n*maxit
jtot = 0
!
!     determine where the matrix splits and choose ql or qr iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
l1 = 1
nm1 = n - 1
!
10 continue
if( l1.gt.n ) &
  &   go to 160
if( l1.gt.1 ) &
  &   e( l1-1 ) = zero
if( l1.le.nm1 ) then
   do 20 m = l1, nm1
      tst = abs( e( m ) )
      if( tst.eq.zero ) &
  &         go to 30
      if( tst.le.( sqrt( abs( d( m ) ) )*sqrt( abs( d( m+ &
  &          1 ) ) ) )*eps ) then
         e( m ) = zero
         go to 30
      end if
20    continue
end if
m = n
!
30 continue
l = l1
lsv = l
lend = m
lendsv = lend
l1 = m + 1
if( lend.eq.l ) &
  &   go to 10
!
!     scale submatrix in rows and columns l to lend
!
anorm = dlanst( 'i', lend-l+1, d( l ), e( l ) )
iscale = 0
if( anorm.eq.zero ) &
  &   go to 10
if( anorm.gt.ssfmax ) then
   iscale = 1
   call dlascl( 'g', 0, 0, anorm, ssfmax, lend-l+1, 1, d( l ), n, &
  &                info )
   call dlascl( 'g', 0, 0, anorm, ssfmax, lend-l, 1, e( l ), n, &
  &                info )
else if( anorm.lt.ssfmin ) then
   iscale = 2
   call dlascl( 'g', 0, 0, anorm, ssfmin, lend-l+1, 1, d( l ), n, &
  &                info )
   call dlascl( 'g', 0, 0, anorm, ssfmin, lend-l, 1, e( l ), n, &
  &                info )
end if
!
!     choose between ql and qr iteration
!
if( abs( d( lend ) ).lt.abs( d( l ) ) ) then
   lend = lsv
   l = lendsv
end if
!
if( lend.gt.l ) then
!
!        ql iteration
!
!        look for small subdiagonal element.
!
40    continue
   if( l.ne.lend ) then
      lendm1 = lend - 1
      do 50 m = l, lendm1
         tst = abs( e( m ) )**2
         if( tst.le.( eps2*abs( d( m ) ) )*abs( d( m+1 ) )+ &
  &             safmin )go to 60
50       continue
   end if
!
   m = lend
!
60    continue
   if( m.lt.lend ) &
  &      e( m ) = zero
   p = d( l )
   if( m.eq.l ) &
  &      go to 80
!
!        if remaining matrix is 2-by-2, use dlae2 or dlaev2
!        to compute its eigensystem.
!
   if( m.eq.l+1 ) then
      if( icompz.gt.0 ) then
         call dlaev2( d( l ), e( l ), d( l+1 ), rt1, rt2, c, s )
         work( l ) = c
         work( n-1+l ) = s
!$$$               call dlasr( 'r', 'v', 'b', n, 2, work( l ),
!$$$     $                     work( n-1+l ), z( 1, l ), ldz )
!
!              *** New starting with version 2.5 ***
!
         tst      = z(l+1)
         z(l+1) = c*tst - s*z(l)
         z(l)   = s*tst + c*z(l)
!              *************************************
      else
         call dlae2( d( l ), e( l ), d( l+1 ), rt1, rt2 )
      end if
      d( l ) = rt1
      d( l+1 ) = rt2
      e( l ) = zero
      l = l + 2
      if( l.le.lend ) &
  &         go to 40
      go to 140
   end if
!
   if( jtot.eq.nmaxit ) &
  &      go to 140
   jtot = jtot + 1
!
!        form shift.
!
   g = ( d( l+1 )-p ) / ( two*e( l ) )
   r = dlapy2( g, one )
   g = d( m ) - p + ( e( l ) / ( g+sign( r, g ) ) )
!
   s = one
   c = one
   p = zero
!
!        inner loop
!
   mm1 = m - 1
   do 70 i = mm1, l, -1
      f = s*e( i )
      b = c*e( i )
      call dlartg( g, f, c, s, r )
      if( i.ne.m-1 ) &
  &         e( i+1 ) = r
      g = d( i+1 ) - p
      r = ( d( i )-g )*s + two*c*b
      p = s*r
      d( i+1 ) = g + p
      g = c*r - b
!
!           if eigenvectors are desired, then save rotations.
!
      if( icompz.gt.0 ) then
         work( i ) = c
         work( n-1+i ) = -s
      end if
!
70    continue
!
!        if eigenvectors are desired, then apply saved rotations.
!
   if( icompz.gt.0 ) then
      mm = m - l + 1
!$$$            call dlasr( 'r', 'v', 'b', n, mm, work( l ), work( n-1+l ),
!$$$     $                  z( 1, l ), ldz )
!
!             *** New starting with version 2.5 ***
!
        call dlasr( 'r', 'v', 'b', 1, mm, work( l ), &
  &                    work( n-1+l ), z( l ), 1 )
!             *************************************
   end if
!
   d( l ) = d( l ) - p
   e( l ) = g
   go to 40
!
!        eigenvalue found.
!
80    continue
   d( l ) = p
!
   l = l + 1
   if( l.le.lend ) &
  &      go to 40
   go to 140
!
else
!
!        qr iteration
!
!        look for small superdiagonal element.
!
90    continue
   if( l.ne.lend ) then
      lendp1 = lend + 1
      do 100 m = l, lendp1, -1
         tst = abs( e( m-1 ) )**2
         if( tst.le.( eps2*abs( d( m ) ) )*abs( d( m-1 ) )+ &
  &             safmin )go to 110
100       continue
   end if
!
   m = lend
!
110    continue
   if( m.gt.lend ) &
  &      e( m-1 ) = zero
   p = d( l )
   if( m.eq.l ) &
  &      go to 130
!
!        if remaining matrix is 2-by-2, use dlae2 or dlaev2
!        to compute its eigensystem.
!
   if( m.eq.l-1 ) then
      if( icompz.gt.0 ) then
         call dlaev2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2, c, s )
!$$$               work( m ) = c
!$$$               work( n-1+m ) = s
!$$$               call dlasr( 'r', 'v', 'f', n, 2, work( m ),
!$$$     $                     work( n-1+m ), z( 1, l-1 ), ldz )
!
!               *** New starting with version 2.5 ***
!
          tst      = z(l)
          z(l)   = c*tst - s*z(l-1)
          z(l-1) = s*tst + c*z(l-1)
!               *************************************
      else
         call dlae2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2 )
      end if
      d( l-1 ) = rt1
      d( l ) = rt2
      e( l-1 ) = zero
      l = l - 2
      if( l.ge.lend ) &
  &         go to 90
      go to 140
   end if
!
   if( jtot.eq.nmaxit ) &
  &      go to 140
   jtot = jtot + 1
!
!        form shift.
!
   g = ( d( l-1 )-p ) / ( two*e( l-1 ) )
   r = dlapy2( g, one )
   g = d( m ) - p + ( e( l-1 ) / ( g+sign( r, g ) ) )
!
   s = one
   c = one
   p = zero
!
!        inner loop
!
   lm1 = l - 1
   do 120 i = m, lm1
      f = s*e( i )
      b = c*e( i )
      call dlartg( g, f, c, s, r )
      if( i.ne.m ) &
  &         e( i-1 ) = r
      g = d( i ) - p
      r = ( d( i+1 )-g )*s + two*c*b
      p = s*r
      d( i ) = g + p
      g = c*r - b
!
!           if eigenvectors are desired, then save rotations.
!
      if( icompz.gt.0 ) then
         work( i ) = c
         work( n-1+i ) = s
      end if
!
120    continue
!
!        if eigenvectors are desired, then apply saved rotations.
!
   if( icompz.gt.0 ) then
      mm = l - m + 1
!$$$            call dlasr( 'r', 'v', 'f', n, mm, work( m ), work( n-1+m ),
!$$$     $                  z( 1, m ), ldz )
!
!           *** New starting with version 2.5 ***
!
      call dlasr( 'r', 'v', 'f', 1, mm, work( m ), work( n-1+m ), &
  &                  z( m ), 1 )
!           *************************************
   end if
!
   d( l ) = d( l ) - p
   e( lm1 ) = g
   go to 90
!
!        eigenvalue found.
!
130    continue
   d( l ) = p
!
   l = l - 1
   if( l.ge.lend ) &
  &      go to 90
   go to 140
!
end if
!
!     undo scaling if necessary
!
140 continue
if( iscale.eq.1 ) then
   call dlascl( 'g', 0, 0, ssfmax, anorm, lendsv-lsv+1, 1, &
  &                d( lsv ), n, info )
   call dlascl( 'g', 0, 0, ssfmax, anorm, lendsv-lsv, 1, e( lsv ), &
  &                n, info )
else if( iscale.eq.2 ) then
   call dlascl( 'g', 0, 0, ssfmin, anorm, lendsv-lsv+1, 1, &
  &                d( lsv ), n, info )
   call dlascl( 'g', 0, 0, ssfmin, anorm, lendsv-lsv, 1, e( lsv ), &
  &                n, info )
end if
!
!     check for no convergence to an eigenvalue after a total
!     of n*maxit iterations.
!
if( jtot.lt.nmaxit ) &
  &   go to 10
do 150 i = 1, n - 1
   if( e( i ).ne.zero ) &
  &      info = info + 1
150 continue
go to 190
!
!     order eigenvalues and eigenvectors.
!
160 continue
if( icompz.eq.0 ) then
!
!        use quick sort
!
   call dlasrt( 'i', n, d, info )
!
else
!
!        use selection sort to minimize swaps of eigenvectors
!
   do 180 ii = 2, n
      i = ii - 1
      k = i
      p = d( i )
      do 170 j = ii, n
         if( d( j ).lt.p ) then
            k = j
            p = d( j )
         end if
170       continue
      if( k.ne.i ) then
         d( k ) = d( i )
         d( i ) = p
!$$$               call dswap( n, z( 1, i ), 1, z( 1, k ), 1 )
!           *** New starting with version 2.5 ***
!
         p    = z(k)
         z(k) = z(i)
         z(i) = p
!           *************************************
      end if
180    continue
end if
!
190 continue
return
!
!     %---------------%
!     | End of dstqrb |
!     %---------------%
!
end
