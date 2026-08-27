!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dsortr
!
!\Description:
!  Sort the array X1 in the order specified by WHICH and optionally
!  applies the permutation to the array X2.
!
!\Usage:
!  call dsortr
!     ( WHICH, APPLY, N, X1, X2 )
!
!\Arguments
!  WHICH   Character*2.  (Input)
!          'LM' -> X1 is sorted into increasing order of magnitude.
!          'SM' -> X1 is sorted into decreasing order of magnitude.
!          'LA' -> X1 is sorted into increasing order of algebraic.
!          'SA' -> X1 is sorted into decreasing order of algebraic.
!
!  APPLY   Logical.  (Input)
!          APPLY = .TRUE.  -> apply the sorted order to X2.
!          APPLY = .FALSE. -> do not apply the sorted order to X2.
!
!  N       Integer.  (INPUT)
!          Size of the arrays.
!
!  X1      Double precision array of length N.  (INPUT/OUTPUT)
!          The array to be sorted.
!
!  X2      Double precision array of length N.  (INPUT/OUTPUT)
!          Only referenced if APPLY = .TRUE.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     12/16/93: Version ' 2.1'.
!               Adapted from the sort routine in LANSO.
!
!\SCCS Information: @(#)
! FILE: sortr.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dsortr (which, apply, n, x1, x2)
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
character*2 which
logical    apply
integer    n
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
Double precision &
  &           x1(0:n-1), x2(0:n-1)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
integer    i, igap, j
Double precision &
  &           temp
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
igap = n / 2
!
if (which .eq. 'SA') then
!
!        X1 is sorted into decreasing order of algebraic.
!
10    continue
   if (igap .eq.0) go to 9000
   do 30 i = igap, n-1
      j = i-igap
20       continue
!
      if (j.lt.0) go to 30
!
      if (x1(j).lt.x1(j+igap)) then
         temp = x1(j)
         x1(j) = x1(j+igap)
         x1(j+igap) = temp
         if (apply) then
            temp = x2(j)
            x2(j) = x2(j+igap)
            x2(j+igap) = temp
         end if
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
!        X1 is sorted into decreasing order of magnitude.
!
40    continue
   if (igap .eq.0) go to 9000
   do 60 i = igap, n-1
      j = i-igap
50       continue
!
      if (j.lt.0) go to 60
!
      if (abs(x1(j)).lt.abs(x1(j+igap))) then
         temp = x1(j)
         x1(j) = x1(j+igap)
         x1(j+igap) = temp
         if (apply) then
            temp = x2(j)
            x2(j) = x2(j+igap)
            x2(j+igap) = temp
         end if
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
!        X1 is sorted into increasing order of algebraic.
!
70    continue
   if (igap .eq.0) go to 9000
   do 90 i = igap, n-1
      j = i-igap
80       continue
!
      if (j.lt.0) go to 90
!
      if (x1(j).gt.x1(j+igap)) then
         temp = x1(j)
         x1(j) = x1(j+igap)
         x1(j+igap) = temp
         if (apply) then
            temp = x2(j)
            x2(j) = x2(j+igap)
            x2(j+igap) = temp
         end if
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
!        X1 is sorted into increasing order of magnitude.
!
100    continue
   if (igap .eq.0) go to 9000
   do 120 i = igap, n-1
      j = i-igap
110       continue
!
      if (j.lt.0) go to 120
!
      if (abs(x1(j)).gt.abs(x1(j+igap))) then
         temp = x1(j)
         x1(j) = x1(j+igap)
         x1(j+igap) = temp
         if (apply) then
            temp = x2(j)
            x2(j) = x2(j+igap)
            x2(j+igap) = temp
         end if
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
!     | End of dsortr |
!     %---------------%
!
end
