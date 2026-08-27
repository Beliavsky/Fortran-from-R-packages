!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dsortc
!
!\Description:
!  Sorts the complex array in XREAL and XIMAG into the order
!  specified by WHICH and optionally applies the permutation to the
!  real array Y. It is assumed that if an element of XIMAG is
!  nonzero, then its negative is also an element. In other words,
!  both members of a complex conjugate pair are to be sorted and the
!  pairs are kept adjacent to each other.
!
!\Usage:
!  call dsortc
!     ( WHICH, APPLY, N, XREAL, XIMAG, Y )
!
!\Arguments
!  WHICH   Character*2.  (Input)
!          'LM' -> sort XREAL,XIMAG into increasing order of magnitude.
!          'SM' -> sort XREAL,XIMAG into decreasing order of magnitude.
!          'LR' -> sort XREAL into increasing order of algebraic.
!          'SR' -> sort XREAL into decreasing order of algebraic.
!          'LI' -> sort XIMAG into increasing order of magnitude.
!          'SI' -> sort XIMAG into decreasing order of magnitude.
!          NOTE: If an element of XIMAG is non-zero, then its negative
!                is also an element.
!
!  APPLY   Logical.  (Input)
!          APPLY = .TRUE.  -> apply the sorted order to array Y.
!          APPLY = .FALSE. -> do not apply the sorted order to array Y.
!
!  N       Integer.  (INPUT)
!          Size of the arrays.
!
!  XREAL,  Double precision array of length N.  (INPUT/OUTPUT)
!  XIMAG   Real and imaginary part of the array to be sorted.
!
!  Y       Double precision array of length N.  (INPUT/OUTPUT)
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
!     xx/xx/92: Version ' 2.1'
!               Adapted from the sort routine in LANSO.
!
!\SCCS Information: @(#)
! FILE: sortc.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dsortc (which, apply, n, xreal, ximag, y)
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
  &           xreal(0:n-1), ximag(0:n-1), y(0:n-1)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
integer    i, igap, j
Double precision &
  &           temp, temp1, temp2
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
Double precision &
  &           dlapy2
external   dlapy2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
igap = n / 2
!
if (which .eq. 'LM') then
!
!        %------------------------------------------------------%
!        | Sort XREAL,XIMAG into increasing order of magnitude. |
!        %------------------------------------------------------%
!
10    continue
   if (igap .eq.0) go to 9000
!
   do 30 i = igap, n-1
      j = i-igap
20       continue
!
      if (j.lt.0) go to 30
!
      temp1 = dlapy2(xreal(j),ximag(j))
      temp2 = dlapy2(xreal(j+igap),ximag(j+igap))
!
      if (temp1.gt.temp2) then
          temp = xreal(j)
          xreal(j) = xreal(j+igap)
          xreal(j+igap) = temp
!
          temp = ximag(j)
          ximag(j) = ximag(j+igap)
          ximag(j+igap) = temp
!
          if (apply) then
              temp = y(j)
              y(j) = y(j+igap)
              y(j+igap) = temp
          end if
      else
          go to 30
      end if
      j = j-igap
      go to 20
30    continue
   igap = igap / 2
   go to 10
!
else if (which .eq. 'SM') then
!
!        %------------------------------------------------------%
!        | Sort XREAL,XIMAG into decreasing order of magnitude. |
!        %------------------------------------------------------%
!
40    continue
   if (igap .eq.0) go to 9000
!
   do 60 i = igap, n-1
      j = i-igap
50       continue
!
      if (j .lt.0) go to 60
!
      temp1 = dlapy2(xreal(j),ximag(j))
      temp2 = dlapy2(xreal(j+igap),ximag(j+igap))
!
      if (temp1.lt.temp2) then
         temp = xreal(j)
         xreal(j) = xreal(j+igap)
         xreal(j+igap) = temp
!
         temp = ximag(j)
         ximag(j) = ximag(j+igap)
         ximag(j+igap) = temp
!
         if (apply) then
            temp = y(j)
            y(j) = y(j+igap)
            y(j+igap) = temp
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
else if (which .eq. 'LR') then
!
!        %------------------------------------------------%
!        | Sort XREAL into increasing order of algebraic. |
!        %------------------------------------------------%
!
70    continue
   if (igap .eq.0) go to 9000
!
   do 90 i = igap, n-1
      j = i-igap
80       continue
!
      if (j.lt.0) go to 90
!
      if (xreal(j).gt.xreal(j+igap)) then
         temp = xreal(j)
         xreal(j) = xreal(j+igap)
         xreal(j+igap) = temp
!
         temp = ximag(j)
         ximag(j) = ximag(j+igap)
         ximag(j+igap) = temp
!
         if (apply) then
            temp = y(j)
            y(j) = y(j+igap)
            y(j+igap) = temp
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
else if (which .eq. 'SR') then
!
!        %------------------------------------------------%
!        | Sort XREAL into decreasing order of algebraic. |
!        %------------------------------------------------%
!
100    continue
   if (igap .eq.0) go to 9000
   do 120 i = igap, n-1
      j = i-igap
110       continue
!
      if (j.lt.0) go to 120
!
      if (xreal(j).lt.xreal(j+igap)) then
         temp = xreal(j)
         xreal(j) = xreal(j+igap)
         xreal(j+igap) = temp
!
         temp = ximag(j)
         ximag(j) = ximag(j+igap)
         ximag(j+igap) = temp
!
         if (apply) then
            temp = y(j)
            y(j) = y(j+igap)
            y(j+igap) = temp
         end if
      else
         go to 120
      endif
      j = j-igap
      go to 110
120    continue
   igap = igap / 2
   go to 100
!
else if (which .eq. 'LI') then
!
!        %------------------------------------------------%
!        | Sort XIMAG into increasing order of magnitude. |
!        %------------------------------------------------%
!
130    continue
   if (igap .eq.0) go to 9000
   do 150 i = igap, n-1
      j = i-igap
140       continue
!
      if (j.lt.0) go to 150
!
      if (abs(ximag(j)).gt.abs(ximag(j+igap))) then
         temp = xreal(j)
         xreal(j) = xreal(j+igap)
         xreal(j+igap) = temp
!
         temp = ximag(j)
         ximag(j) = ximag(j+igap)
         ximag(j+igap) = temp
!
         if (apply) then
            temp = y(j)
            y(j) = y(j+igap)
            y(j+igap) = temp
         end if
      else
         go to 150
      endif
      j = j-igap
      go to 140
150    continue
   igap = igap / 2
   go to 130
!
else if (which .eq. 'SI') then
!
!        %------------------------------------------------%
!        | Sort XIMAG into decreasing order of magnitude. |
!        %------------------------------------------------%
!
160    continue
   if (igap .eq.0) go to 9000
   do 180 i = igap, n-1
      j = i-igap
170       continue
!
      if (j.lt.0) go to 180
!
      if (abs(ximag(j)).lt.abs(ximag(j+igap))) then
         temp = xreal(j)
         xreal(j) = xreal(j+igap)
         xreal(j+igap) = temp
!
         temp = ximag(j)
         ximag(j) = ximag(j+igap)
         ximag(j+igap) = temp
!
         if (apply) then
            temp = y(j)
            y(j) = y(j+igap)
            y(j+igap) = temp
         end if
      else
         go to 180
      endif
      j = j-igap
      go to 170
180    continue
   igap = igap / 2
   go to 160
end if
!
9000 continue
return
!
!     %---------------%
!     | End of dsortc |
!     %---------------%
!
end
