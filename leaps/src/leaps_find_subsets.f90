module leaps_find_subsets

! a module of routines for finding and recording best-fitting subsets of
! regression variables

!   version 1.10, 17 february 2004
!   author: alan miller
!           formerly of csiro division of mathematical & information sciences
!   phone: (+61) 3 9592-5085
!   e-mail: amiller @ bigpond.net.au
!   www-page: http://users.bigpond.net.au/amiller/

! 17 feb 2004 correction to subroutine efroym for the case in which all the
!             variables are selected. thanks to david jones.
! 12 nov 1999 made changes to routines exadd1 & add2 to prevent the calculation
!             of negative residual sums of squares which could occur in cases
!             in which the true rss is zero.   routine seq2 changed to avoid
!             cycling.
! 24 may 2000 changed lsq_kind to dp (cosmetic change only)
! 4 june 2000 added routine random_pick which picks a random set of variables.
! 29 aug 2002 set value of size in subroutine efroym when ier /= 0.

use leaps_lsq

implicit none

integer, save                 :: max_size, nbest, lopt_dim1
real (dp), allocatable, save  :: bound(:), ress(:,:)
integer, allocatable, save    :: lopt(:,:)


contains

subroutine init_subsets(nvar_max, fit_const, nvar)

integer, intent(in)           :: nvar_max
logical, intent(in)           :: fit_const
integer, intent(in), optional :: nvar

!     local variables

integer    :: i, ier
real (dp)  :: eps = 1.e-14
logical    :: lindep(ncol)

!     the lsq module has probably already been initialized, but just in case ..

if (.not. initialized) then
  if (present(nvar)) call startup(nvar, fit_const)
end if

if (fit_const) then
  max_size = nvar_max + 1
else
  max_size = nvar_max
end if

lopt_dim1 = max_size * (max_size + 1) / 2
if (allocated(bound)) deallocate(bound, ress, lopt)
allocate (bound(max_size), ress(max_size,nbest), lopt(lopt_dim1,nbest))

bound = huge(eps)
ress  = huge(eps)
lopt  = 0

call tolset(eps)
call sing(lindep, ier)

call ss()
do i = 1, max_size
  call report(i, rss(i))
end do

return
end subroutine init_subsets


subroutine add1(first, last, ss, smax, jmax, ier)

! calculate the reduction in residual sum of squares when one variable,
! selected from those in positions first .. last, is added in position first,
! given that the variables in positions 1 .. first-1 (if any) are already
! included.

integer, intent(in)     :: first, last
integer, intent(out)    :: jmax, ier
real (dp), intent(out)  :: ss(:), smax

!     local variables

integer    :: j, inc, pos, row, col
real (dp)  :: zero = 0.0_dp, diag, dy, ssqx, sxx(ncol), sxy(ncol)

!     check call arguments

jmax = 0
smax = zero
ier = 0
if (first > ncol) ier = 1
if (last < first) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0) return

!     accumulate sums of squares & products from row first

sxx(first:last) = zero
sxy(first:last) = zero
inc = ncol - last
pos = row_ptr(first)
do row = first, last
  diag = d(row)
  dy = diag * rhs(row)
  sxx(row) = sxx(row) + diag
  sxy(row) = sxy(row) + dy
  do col = row+1, last
    sxx(col) = sxx(col) + diag * r(pos)**2
    sxy(col) = sxy(col) + dy * r(pos)
    pos = pos + 1
  end do
  pos = pos + inc
end do

!     incremental sum of squares for a variable = sxy * sxy / sxx.
!     calculate whenever sqrt(sxx) > tol for that variable.

do j = first, last
  ssqx = sxx(j)
  if (sqrt(ssqx) > tol(j)) then
    ss(j) = sxy(j)**2 / sxx(j)
    if (ss(j) > smax) then
      smax = ss(j)
      jmax = j
    end if
  else
    ss(j) = zero
  end if
end do

return
end subroutine add1


subroutine add2(first, last, smax, j1, j2, ier)

!     calculate the maximum reduction in residual sum of squares when 2
!     variables, selected from those in positions first .. last, are
!     added, given that the variables in positions 1 .. first-1 (if
!     any) are already included.    j1, j2 are the positions of the two
!     best variables.   n.b. j2 < j1.

integer, intent(in)     :: first, last
integer, intent(out)    :: j1, j2, ier
real (dp), intent(out)  :: smax

!     local variables

integer    :: start, i1, i2, row, pos1, pos2, inc
real (dp)  :: zero = 0.0_dp, temp, det, two = 2.0, sxx(ncol), sxy(ncol), sx1x2

!     check call arguments

smax = zero
j1 = 0
j2 = 0
ier = 0
if (first > ncol) ier = 1
if (last <= first) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0) return

start = row_ptr(first)

!     cycle through all pairs of variables from those between first & last.

do i1 = first, last
  sxx(i1) = d(i1)
  sxy(i1) = d(i1) * rhs(i1)
  pos1 = start + i1 - first - 1
  do row = first, i1-1
    temp = d(row) * r(pos1)
    sxx(i1) = sxx(i1) + temp*r(pos1)
    sxy(i1) = sxy(i1) + temp*rhs(row)
    pos1 = pos1 + ncol - row - 1
  end do

  do i2 = first, i1-1
    pos1 = start + i1 - first - 1
    pos2 = start + i2 - first - 1
    sx1x2 = zero
    do row = first, i2-1
      sx1x2 = sx1x2 + d(row)*r(pos1)*r(pos2)
      inc = ncol - row - 1
      pos1 = pos1 + inc
      pos2 = pos2 + inc
    end do
    sx1x2 = sx1x2 + d(i2)*r(pos1)

!     calculate reduction in rss for pair i1, i2.
!     the sum of squares & cross-products are in:
!              ( sxx(i1)  sx1x2   )      ( sxy(i1) )
!              ( sx1x2    sxx(i2) )      ( sxy(i2) )

    det = max( (sxx(i1) * sxx(i2) - sx1x2**2), zero)
    temp = sqrt(det)
    if (temp < tol(i1)*sqrt(sxx(i2)) .or.             &
        temp < tol(i2)*sqrt(sxx(i1))) cycle
    temp = ((sxx(i2)*sxy(i1) - two*sx1x2*sxy(i2))*sxy(i1) + sxx(i1)*sxy(i2)**2) &
           / det
    if (temp > smax) then
      smax = temp
      j1 = i1
      j2 = i2
    end if
  end do ! i2 = first, i1-1
end do   ! i1 = first, last

return
end subroutine add2


subroutine bakwrd(first, last, ier)

!     backward elimination from variables in positions first .. last.
!     if first > 1, variables in positions prior to this are forced in.
!     if last < ncol, variables in positions after this are forced out.
!     on exit, the array vorder contains the numbers of the variables
!     in the order in which they were deleted.

integer, intent(in)  :: first, last
integer, intent(out) :: ier

!     local variables

integer    :: pos, jmin, i
real (dp)  :: ss(last), smin

!     check call arguments

ier = 0
if (first >= ncol) ier = 1
if (last <= 1) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0) return

!     for pos = last, ..., first+1 call drop1 to find best variable to
!     find which variable to drop next.

do pos = last, first+1, -1
  call drop1(first, pos, ss, smin, jmin, ier)
  call exdrop1(first, pos, ss, smin, jmin)
  if (jmin > 0 .and. jmin < pos) then
    call vmove(jmin, pos, ier)
    if (nbest > 0) then
      do i = jmin, pos-1
        call report(i, rss(i))
      end do
    end if
  end if
end do

return
end subroutine bakwrd


subroutine drop1(first, last, ss, smin, jmin, ier)

! calculate the increase in the residual sum of squares when the variable in
! position j is dropped from the model (i.e. moved to position last),
! for j = first, ..., last-1.

integer, intent(in)     :: first, last
integer, intent(out)    :: jmin, ier
real (dp), intent(out)  :: ss(:), smin

!     local variables

integer    :: j, pos1, inc, pos, row, col, i
real (dp)  :: large = huge(1.0_dp), zero = 0.0_dp, d1, rhs1, d2, x, wk(last), &
              vsmall = tiny(1.0_dp)

!     check call arguments

jmin = 0
smin = large
ier = 0
if (first > ncol) ier = 1
if (last < first) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0) return

!     pos1 = position of first element of row first in r.

pos1 = row_ptr(first)
inc = ncol - last

!     start of outer cycle for the variable to be dropped.

do j = first, last
  d1 = d(j)
  if (sqrt(d1) < tol(j)) then
    ss(j) = zero
    smin = zero
    jmin = j
    go to 50
  end if
  rhs1 = rhs(j)
  if (j == last) go to 40

!     copy row j of r into wk.

  pos = pos1
  do i = j+1, last
    wk(i) = r(pos)
    pos = pos + 1
  end do
  pos = pos + inc

!     lower the variable past each row.

  do row = j+1, last
    x = wk(row)
    d2 = d(row)
    if (abs(x) * sqrt(d1) < tol(row) .or. d2 < vsmall) then
      pos = pos + ncol - row
      cycle
    end if
    d1 = d1 * d2 / (d2 + d1 * x**2)
    do col = row+1, last
      wk(col) = wk(col) - x * r(pos)
      pos = pos + 1
    end do
    rhs1 = rhs1 - x * rhs(row)
    pos = pos + inc
  end do
  40 ss(j) = rhs1 * d1 * rhs1
  if (ss(j) < smin) then
    jmin = j
    smin = ss(j)
  end if

!     update position of first element in row of r.

  50 if (j < last) pos1 = pos1 + ncol - j
end do

return
end subroutine drop1


subroutine efroym(first, last, fin, fout, size, ier, lout)

!     efroymson's stepwise regression from variables in positions first,
!     ..., last.  if first > 1, variables in positions prior to this are
!     forced in.  if last < ncol, variables in positions after this are
!     forced out.

!     a report is written to unit lout if lout >= 0.

integer, intent(in)    :: first, last, lout
integer, intent(out)   :: size, ier
real (dp), intent(in)  :: fin, fout

!     local variables

integer    :: jmax, jmin, i
real (dp)  :: one = 1.0, eps, zero = 0.0, ss(last), smax, base, var, f, smin

!     check call arguments

ier = 0
if (first >= ncol) ier = 1
if (last <= 1) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (fin < fout .or. fin <= zero) ier = ier + 256
if (nobs <= ncol) ier = ier + 512
if (ier /= 0) then
  size = 0
  return
end if

!     eps approximates the smallest quantity such that the calculated value of
!     (1 + eps) is > 1.   it is used to test for a perfect fit (rss = 0).

eps = epsilon(one)

!     size = number of variables in the current subset

size = first - 1

!     find the best variable to add next

20 call add1(size+1, last, ss, smax, jmax, ier)
if (nbest > 0) call exadd1(size+1, smax, jmax, ss, last)

!     calculate 'f-to-enter' value

if (size > 0) then
  base = rss(size)
else
  base = rss(1) + ss(1)
end if
var = (base - smax) / (nobs - size - 1)
if (var < eps*base) then
  ier = -1
  f = zero
else
  f = smax / var
end if
if (lout >= 0) write(lout, 900) vorder(jmax), f
900 format(' best variable to add:  ', i4, '  f-to-enter = ', f10.2)

!     exit if f < fin or ier < 0 (perfect fit)

if (f < fin .or. ier < 0) return

!     add the variable to the subset (in position first).

if (lout >= 0) write(lout, '(50x, "variable added")')
size = size + 1
if (jmax > first) call vmove(jmax, first, ier)
do i = first, min(jmax-1, max_size)
  call report(i, rss(i))
end do

!     see whether a variable entered earlier can be deleted now.

30 if (size <= first) go to 20
call drop1(first+1, size, ss, smin, jmin, ier)
call exdrop1(first+1, size, ss, smin, jmin)
var = rss(size) / (nobs - size)
f = smin / var
if (lout >= 0) write(lout, 910) vorder(jmin), f
910 format(' best variable to drop: ', i4, '  f-to-drop  = ', f10.2)

if (f < fout) then
  if (lout >= 0) write(lout, '(50x, "variable dropped")')
  call vmove(jmin, size, ier)
  if (nbest > 0) then
    do i = jmin, size-1
      call report(i, rss(i))
    end do
  end if
  size = size - 1
  go to 30
end if

if (size >= last) return
go to 20
end subroutine efroym


subroutine exadd1(ivar, smax, jmax, ss, last)

!     update the nbest subsets of ivar variables found from a call
!     to subroutine add1.

integer, intent(in)    :: ivar, jmax, last
real (dp), intent(in)  :: smax, ss(:)

!     local variables

real (dp)  :: zero = 0.0_dp, ssbase, sm, temp, wk(last)
integer    :: i, j, ltemp, jm

if (jmax == 0) return
if (ivar <= 0) return
if (ivar > max_size) return
ltemp = vorder(ivar)
jm = jmax
sm = smax
if (ivar > 1) ssbase = rss(ivar-1)
if (ivar == 1) ssbase = rss(ivar) + ss(1)
wk(ivar:last) = ss(ivar:last)

do i = 1, nbest
  temp = max(ssbase - sm, zero)
  if (temp >= bound(ivar)) exit
  vorder(ivar) = vorder(jm)
  if (jm == ivar) vorder(ivar) = ltemp
  call report(ivar, temp)
  if (i >= nbest) exit
  wk(jm) = zero
  sm = zero
  jm = 0
  do j = ivar, last
    if (wk(j) <= sm) cycle
    jm = j
    sm = wk(j)
  end do
  if (jm == 0) exit
end do

!     restore vorder(ivar)

vorder(ivar) = ltemp

return
end subroutine exadd1


subroutine exdrop1(first, last, ss, smin, jmin)
! record any new subsets of (last-1) variables found from a call to drop1

integer, intent(in)    :: first, last, jmin
real (dp), intent(in)  :: ss(:), smin

! local variables
integer    :: list(1:last), i
real (dp)  :: rss_last, ssq

if (jmin == 0 .or. last < 1 .or. last-1 > max_size) return

rss_last = rss(last)
if (rss_last + smin > bound(last-1)) return

list = vorder(1:last)
do i = first, last-1
  vorder(i:last-1) = list(i+1:last)
  ssq = rss_last + ss(i)
  call report(last-1, ssq)
  vorder(i) = list(i)
end do

return
end subroutine exdrop1


subroutine forwrd(first, last, ier)

!     forward selection from variables in positions first .. last.
!     if first > 1, variables in positions prior to this are forced in.
!     if last < ncol, variables in positions after this are forced out.
!     on exit, the array vorder contains the numbers of the variables
!     in the order in which they were added.

integer, intent(in)  :: first, last
integer, intent(out) :: ier

!     local variables

integer    :: pos, jmax
real (dp)  :: ss(last), smax

!     check call arguments

ier = 0
if (first >= ncol) ier = 1
if (last <= 1) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0) return

!     for pos = first .. max_size, call add1 to find best variable to put
!     into position pos.

do pos = first, max_size
  call add1(pos, last, ss, smax, jmax, ier)
  if (nbest > 0) call exadd1(pos, smax, jmax, ss, last)

!     move the best variable to position pos.

  if (jmax > pos) call vmove(jmax, pos, ier)
end do

return
end subroutine forwrd


subroutine report(nv, ssq)

!     update record of the best nbest subsets of nv variables, if
!     necessary, using ssq.

integer, intent(in)    :: nv
real (dp), intent(in)  :: ssq

!     local variables

integer    :: rank, pos1, j, list(nv)
real (dp)  :: under1 = 0.99999_dp, above1 = 1.00001_dp

!     if residual sum of squares (ssq) for the new subset > the
!     appropriate bound, return.

if(nv > max_size) return
if(ssq >= bound(nv)) return
pos1 = (nv*(nv-1))/2 + 1

!     find rank of the new subset

do rank = 1, nbest
  if(ssq < ress(nv,rank)*above1) then
    list = vorder(1:nv)
    call shell(list, nv)

!     check list of variables if ssq is almost equal to ress(nv,rank) -
!     to avoid including the same subset twice.

    if (ssq > ress(nv,rank)*under1) then
      if (same_vars(list, lopt(pos1:,rank), nv)) return
    end if

!     record the new subset, and move the others down one place.

    do j = nbest-1, rank, -1
      ress(nv,j+1) = ress(nv,j)
      lopt(pos1:pos1+nv-1, j+1) = lopt(pos1:pos1+nv-1, j)
    end do
    ress(nv,rank) = ssq
    lopt(pos1:pos1+nv-1, rank) = list(1:nv)
    bound(nv) = ress(nv,nbest)
    return
  end if
end do

return
end subroutine report



subroutine shell(l, n)

!      perform a shell-sort on integer array l, sorting into increasing order.

!      latest revision - 5 july 1995

integer, intent(in)     :: n
integer, intent(in out) :: l(:)

!     local variables
integer   :: start, finish, temp, new, i1, i2, incr, it

incr = n
do
  incr = incr/3
  if (incr == 2*(incr/2)) incr = incr + 1
  do start = 1, incr
    finish = n

!      temp contains the element being compared; it holds its current
!      location.   it is compared with the elements in locations
!      it+incr, it+2.incr, ... until a larger element is found.   all
!      smaller elements move incr locations towards the start.   after
!      each time through the sequence, the finish is decreased by incr
!      until finish <= incr.

    20 i1 = start
    temp = l(i1)
    it = i1

!      i2 = location of element new to be compared with temp.
!      test i2 <= finish.

    do
      i2 = i1 + incr
      if (i2 > finish) then
        if (i1 > it) l(i1) = temp
        finish = finish - incr
        exit
      end if
      new = l(i2)

!     if temp > new, move new to lower-numbered position.

      if (temp > new) then
        l(i1) = new
        i1 = i2
        cycle
      end if

!     temp <= new so do not swap.
!     use new as the next temp.

      if (i1 > it) l(i1) = temp
      i1 = i2
      temp = new
      it = i1

!     repeat until finish <= incr.
    end do

    if (finish > incr) go to 20
  end do

!      repeat until incr = 1.

  if (incr <= 1) return
end do

return
end subroutine shell



function same_vars(list1, list2, n) result(same)

logical              :: same
integer, intent(in)  :: n, list1(:), list2(:)

same = all(list1(1:n) == list2(1:n))

return
end function same_vars



subroutine seq2(first, last, ier)

! sequential replacement algorithm applied to the variables in positions
! first, ..., last.   2 variables at a time are added or replaced.
! if first > 1, variables in positions prior to this are forced in.
! if last < np, variables in positions after this are left out.

integer, intent(in)  :: first, last
integer, intent(out) :: ier

!     local variables

integer  :: nv, nsize

!     check call arguments

ier = 0
if (first >= ncol) ier = 1
if (last <= 1) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0 .or. nbest <= 0) return

nv = min(max_size, last-1)

!     outer loop; size = current size of subset being considered.

do nsize = first+1, nv
  call replace2(first, last, nsize)
end do

return
end subroutine seq2



subroutine replace2(first, last, nsize)
! replace 2 variables at a time from those in positions first, ..., nsize
! with 2 from positions nsize, .., last - if they reduce the rss.

integer, intent(in)  :: first, last, nsize

! local variables

integer              :: ier, j1, j2, pos1, pos2, best(2), i, iwk(last)
real (dp)            :: smax, rssnew, rssmin, save_rss
real (dp), parameter :: zero = 0.0_dp

10 best(1) = 0
best(2) = 0
rssmin = rss(nsize)

!     two loops to place all pairs of variables in positions nsize-1 and nsize.
!     pos1 = destination for variable from position nsize.
!     pos2 = destination for variable from position nsize-1.

do pos1 = first, nsize
  do pos2 = pos1, nsize-1
    call add2(nsize-1, last, smax, j1, j2, ier)

    if (j1+j2 > nsize + nsize - 1) then
      rssnew = max(rss(nsize-2) - smax, zero)
      if (rssnew < rssmin) then
        best(1) = vorder(j1)
        best(2) = vorder(j2)
        iwk(1:nsize-2) = vorder(1:nsize-2)
        rssmin = rssnew
      end if
    end if

    call vmove(nsize-1, pos2, ier)
  end do
  call vmove(nsize, pos1, ier)
  do i = pos1, nsize
    call report(i, rss(i))
  end do
end do

!     if any replacement reduces the rss, make the best one.

if (best(1) + best(2) > 0) then
  iwk(nsize-1) = best(2)
  iwk(nsize) = best(1)
  save_rss = rss(nsize)
  call reordr(iwk, nsize, 1, ier)
  do i = first, nsize
    call report(i, rss(i))
  end do

!    the calculated value of rssmin above is only a rough approximation to
!    the real residual sum of squares, thiugh usually good enough.
!    the new value of rss(nsize) is more accurate.   it is used below
!    to avoid cycling when several subsets give the same rss.

  if (rss(nsize) < save_rss) go to 10
end if

return
end subroutine replace2



subroutine seqrep(first, last, ier)

!     sequential replacement algorithm applied to the variables in
!     positions first, ..., last.
!     if first > 1, variables in positions prior to this are forced in.
!     if last < ncol, variables in positions after this are forced out.

integer, intent(in)  :: first, last
integer, intent(out) :: ier

!     local variables

integer    :: nv, size, start, best, from, i, jmax, count, j
real (dp)  :: zero = 0.0_dp, ssred, ss(last), smax

!     check call arguments

ier = 0
if (first >= ncol) ier = 1
if (last <= 1) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0 .or. nbest <= 0) return

nv = min(max_size, last-1)

!     outer loop; size = current size of subset being considered.

do size = first, nv
  count = 0
  start = first
  10 ssred = zero
  best = 0
  from = 0

!     find the best variable from those in positions size+1, ..., last
!     to replace the one in position size.   then rotate variables in
!     positions start, ..., size.

  do i = start, size
    call add1(size, last, ss, smax, jmax, ier)
    if (jmax > size) then
      call exadd1(size, smax, jmax, ss, last)
      if (smax > ssred) then
        ssred = smax
        best = jmax
        if (i < size) then
          from = size + start - i - 1
        else
          from = size
        end if
      end if
    end if
    if (i < size) call vmove(size, start, ier)
    do j = start, size-1
      call report(j, rss(j))
    end do
  end do ! i = start, size

!     if any replacement reduces the rss, make the best one.
!     move variable from position from to size.
!     move variable from position best to first.

  if (best > size) then
    if (from < size) call vmove(from, size, ier)
    call vmove(best, first, ier)
    do j = first, best-1
      call report(j, rss(j))
    end do
    count = 0
    start = first + 1
  else
    count = count + 1
  end if

!     repeat until count = size - start + 1

  if (count <= size - start) go to 10
end do

return
end subroutine seqrep


subroutine xhaust(first, last, ier)

!     exhaustive search algorithm, using leaps and bounds, applied to
!     the variables in positions first, ..., last.
!     if first > 1, variables in positions prior to this are forced in.
!     if last < ncol, variables in positions after this are forced out.

integer, intent(in)  :: first, last
integer, intent(out) :: ier

!     local variables

integer    :: row, i, jmax, ipt, newpos, iwk(max_size)
real (dp)  :: ss(last), smax, temp

!     check call arguments

ier = 0
if (first >= ncol) ier = 1
if (last <= 1) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0 .or. nbest <= 0) return

!     record subsets contained in the initial ordering, including check
!     for variables which are linearly related to earlier variables.
!     this should be redundant if the user has first called sing and
!     init_subsets.

do row = first, max_size
  if (d(row) <= tol(row)) then
    ier = -999
    return
  end if
  call report(row, rss(row))
end do

!     iwk(i) contains the upper limit for the i-th simulated do-loop for
!     i = first, ..., max_size-1.
!     ipt points to the current do loop.

iwk(first:max_size) = last

!     innermost loop.
!     find best possible variable for position max_size from those in
!     positions max_size, .., iwk(max_size).

30 call add1(max_size, iwk(max_size), ss, smax, jmax, ier)
call exadd1(max_size, smax, jmax, ss, iwk(max_size))

!     move to next lower numbered loop which has not been exhausted.

ipt = max_size - 1
40 if (ipt >= iwk(ipt)) then
  ipt = ipt - 1
  if (ipt >= first) go to 40
  return
end if

!     lower variable from position ipt to position iwk(ipt).
!     record any good new subsets found by the move.

newpos = iwk(ipt)
call vmove(ipt, newpos, ier)
do i = ipt, min(max_size, newpos-1)
  call report(i, rss(i))
end do

!     reset all ends of loops for i >= ipt.

iwk(ipt:max_size) = newpos - 1

!     if residual sum of squares for all variables above position newpos
!     is greater than bound(i), no better subsets of size i can be found
!     inside the current loop.

temp = rss(newpos-1)
do i = ipt, max_size
  if (temp > bound(i)) go to 80
end do
if (iwk(max_size) > max_size) go to 30
ipt = max_size - 1
go to 40

80 ipt = i - 1
if (ipt < first) return
go to 40

end subroutine xhaust



subroutine random_pick(first, last, npick)
! pick npick variables at random from those in positions first, ..., last
! and move them to occupy positions starting from first.

integer, intent(in)  :: first, last, npick

! local variables

integer  :: first2, i, ilist(1:last), j, k, navail
real (dp) :: r

navail = last + 1 - first
if (npick >= navail .or. npick <= 0) return
do i = first, last
  ilist(i) = vorder(i)
end do

first2 = first
do i = 1, npick
  call random_number(r)
  k = first2 + int(r * real(navail, dp))
  if (k > first2) then
    j = ilist(first2)
    ilist(first2) = ilist(k)
    ilist(k) = j
  end if
  first2 = first2 + 1
  navail = navail - 1
end do

call reordr(ilist(first:), npick, first, i)

return
end subroutine random_pick

end module leaps_find_subsets
