!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dsesrt
!
!\Description:
!  Sort the array X in the order specified by WHICH and optionally
!  apply the permutation to the columns of the matrix A.
!
!\Usage:
!  call dsesrt
!     ( WHICH, APPLY, N, X, NA, A, LDA)
!
!\Arguments
!  WHICH   Character*2.  (Input)
!          'LM' -> X is sorted into increasing order of magnitude.
!          'SM' -> X is sorted into decreasing order of magnitude.
!          'LA' -> X is sorted into increasing order of algebraic.
!          'SA' -> X is sorted into decreasing order of algebraic.
!
!  APPLY   Logical.  (Input)
!          APPLY = .TRUE.  -> apply the sorted order to A.
!          APPLY = .FALSE. -> do not apply the sorted order to A.
!
!  N       Integer.  (INPUT)
!          Dimension of the array X.
!
!  X      Double precision array of length N.  (INPUT/OUTPUT)
!          The array to be sorted.
!
!  NA      Integer.  (INPUT)
!          Number of rows of the matrix A.
!
!  A      Double precision array of length NA by N.  (INPUT/OUTPUT)
!
!  LDA     Integer.  (INPUT)
!          Leading dimension of A.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Routines
!     dswap  Level 1 BLAS that swaps the contents of two vectors.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     12/15/93: Version ' 2.1'.
!               Adapted from the sort routine in LANSO and
!               the ARPACK code dsortr
!
!\SCCS Information: @(#)
! FILE: sesrt.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dsesrt (which, apply, n, x, na, a, lda)
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
character*2 which
logical    apply
integer    lda, n, na
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
Double precision &
  &           x(0:n-1), a(lda, 0:n-1)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
integer    i, igap, j
Double precision &
  &           temp
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
external   dswap
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
igap = n / 2
!
if (which .eq. 'SA') then
!
!        X is sorted into decreasing order of algebraic.
!
10    continue
   if (igap .eq.0) go to 9000
   do 30 i = igap, n-1
      j = i-igap
20       continue
!
      if (j.lt.0) go to 30
!
      if (x(j).lt.x(j+igap)) then
         temp = x(j)
         x(j) = x(j+igap)
         x(j+igap) = temp
         if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
      else
         go to 30
      endif
      j = j-igap
      go to 20
30    continue
   igap = igap / 2
   go to 10
!
else if (which .eq. 'SM') then
!
!        X is sorted into decreasing order of magnitude.
!
40    continue
   if (igap .eq.0) go to 9000
   do 60 i = igap, n-1
      j = i-igap
50       continue
!
      if (j.lt.0) go to 60
!
      if (abs(x(j)).lt.abs(x(j+igap))) then
         temp = x(j)
         x(j) = x(j+igap)
         x(j+igap) = temp
         if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
      else
         go to 60
      endif
      j = j-igap
      go to 50
60    continue
   igap = igap / 2
   go to 40
!
else if (which .eq. 'LA') then
!
!        X is sorted into increasing order of algebraic.
!
70    continue
   if (igap .eq.0) go to 9000
   do 90 i = igap, n-1
      j = i-igap
80       continue
!
      if (j.lt.0) go to 90
!
      if (x(j).gt.x(j+igap)) then
         temp = x(j)
         x(j) = x(j+igap)
         x(j+igap) = temp
         if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
      else
         go to 90
      endif
      j = j-igap
      go to 80
90    continue
   igap = igap / 2
   go to 70
!
else if (which .eq. 'LM') then
!
!        X is sorted into increasing order of magnitude.
!
100    continue
   if (igap .eq.0) go to 9000
   do 120 i = igap, n-1
      j = i-igap
110       continue
!
      if (j.lt.0) go to 120
!
      if (abs(x(j)).gt.abs(x(j+igap))) then
         temp = x(j)
         x(j) = x(j+igap)
         x(j+igap) = temp
         if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
      else
         go to 120
      endif
      j = j-igap
      go to 110
120    continue
   igap = igap / 2
   go to 100
end if
!
9000 continue
return
!
!     %---------------%
!     | End of dsesrt |
!     %---------------%
!
end
