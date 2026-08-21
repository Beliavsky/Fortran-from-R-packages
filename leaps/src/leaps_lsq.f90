module leaps_lsq

!  module for unconstrained linear least-squares calculations.
!  the algorithm is suitable for updating ls calculations as more
!  data are added.   this is sometimes called recursive estimation.
!  only one dependent variable is allowed.
!  based upon applied statistics algorithm as 274.
!  translation from fortran 77 to fortran 90 by alan miller.
!  a function, varprd, has been added for calculating the variances
!  of predicted values, and this uses a subroutine bksub2.

!  version 1.14, 19 august 2002 - elf90 compatible version
!  author: alan miller
!  e-mail : amiller @ bigpond.net.au
!  www-pages: http://www.ozemail.com.au/~milleraj
!             http://users.bigpond.net.au/amiller/

!  bug fixes:
!  1. in regcf a call to tolset has been added in case the user had
!     not set tolerances.
!  2. in sing, each time a singularity is detected, unless it is in the
!     variables in the last position, includ is called.   includ assumes
!     that a new observation is being added and increments the number of
!     cases, nobs.   the line:  nobs = nobs - 1 has been added.
!  3. row_ptr was left out of the deallocate statement in routine startup
!     in version 1.07.
!  4. in cov, now calls ss if rss_set = .false.  29 august 1997
!  5. in tolset, correction to accomodate negative values of d.  19 august 2002

!  other changes:
!  1. array row_ptr added 18 july 1997.   this points to the first element
!     stored in each row thus saving a small amount of time needed to
!     calculate its position.
!  2. optional parameter, eps, added to routine tolset, so that the user
!     can specify the accuracy of the input data.
!  3. cosmetic change of lsq_kind to dp (`double precision')
!  4. change to routine sing to use row_ptr rather than calculate the position
!     of first elements in each row.

!  the public variables are:
!  dp       = a kind parameter for the floating-point quantities calculated
!             in this module.   see the more detailed explanation below.
!             this kind parameter should be used for all floating-point
!             arguments passed to routines in this module.

!  nobs    = the number of observations processed to date.
!  ncol    = the total number of variables, including one for the constant,
!            if a constant is being fitted.
!  r_dim   = the dimension of array r = ncol*(ncol-1)/2
!  vorder  = an integer vector storing the current order of the variables
!            in the qr-factorization.   the initial order is 0, 1, 2, ...
!            if a constant is being fitted, or 1, 2, ... otherwise.
!  initialized = a logical variable which indicates whether space has
!                been allocated for various arrays.
!  tol_set = a logical variable which is set when subroutine tolset has
!            been called to calculate tolerances for use in testing for
!            singularities.
!  rss_set = a logical variable indicating whether residual sums of squares
!            are available and usable.
!  d()     = array of row multipliers for the cholesky factorization.
!            the factorization is x = q.sqrt(d).r where q is an ortho-
!            normal matrix which is not stored, d is a diagonal matrix
!            whose diagonal elements are stored in array d, and r is an
!            upper-triangular matrix with 1's as its diagonal elements.
!  rhs()   = vector of rhs projections (after scaling by sqrt(d)).
!            thus q'y = sqrt(d).rhs
!  r()     = the upper-triangular matrix r.   the upper triangle only,
!            excluding the implicit 1's on the diagonal, are stored by
!            rows.
!  tol()   = array of tolerances used in testing for singularities.
!  rss()   = array of residual sums of squares.   rss(i) is the residual
!            sum of squares with the first i variables in the model.
!            by changing the order of variables, the residual sums of
!            squares can be found for all possible subsets of the variables.
!            the residual sum of squares with no variables in the model,
!            that is the total sum of squares of the y-values, can be
!            calculated as rss(1) + d(1)*rhs(1)^2.   if the first variable
!            is a constant, then rss(1) is the sum of squares of
!            (y - ybar) where ybar is the average value of y.
!  sserr   = residual sum of squares with all of the variables included.
!  row_ptr() = array of indices of first elements in each row of r.
!
!--------------------------------------------------------------------------

!     general declarations

implicit none

integer, save                :: nobs, ncol, r_dim
integer, allocatable, save   :: vorder(:), row_ptr(:)
logical, save                :: initialized = .false.,                  &
                                tol_set = .false., rss_set = .false.

! note. dp is being set to give at least 12 decimal digit
!       representation of floating point numbers.   this should be adequate
!       for most problems except the fitting of polynomials.   dp is
!       being set so that the same code can be run on pcs and unix systems,
!       which will usually represent floating-point numbers in `double
!       precision', and other systems with larger word lengths which will
!       give similar accuracy in `single precision'.

integer, parameter           :: dp = kind(1.0d0)
real (dp), allocatable, save :: d(:), rhs(:), r(:), tol(:), rss(:)
real (dp), save              :: zero = 0.0_dp, one = 1.0_dp, vsmall
real (dp), save              :: sserr, toly

public                       :: dp, nobs, ncol, r_dim, vorder, row_ptr, &
                                initialized, tol_set, rss_set,          &
                                d, rhs, r, tol, rss, sserr
private                      :: zero, one, vsmall


contains

subroutine startup(nvar, fit_const)

!     allocates dimensions for arrays and initializes to zero
!     the calling program must set nvar = the number of variables, and
!     fit_const = .true. if a constant is to be included in the model,
!     otherwise fit_const = .false.
!
!--------------------------------------------------------------------------

implicit none
integer, intent(in)  :: nvar
logical, intent(in)  :: fit_const

!     local variable
integer   :: i

vsmall = 10. * tiny(zero)

nobs = 0
if (fit_const) then
  ncol = nvar + 1
else
  ncol = nvar
end if

if (initialized) deallocate(d, rhs, r, tol, rss, vorder, row_ptr)
r_dim = ncol * (ncol - 1)/2
allocate( d(ncol), rhs(ncol), r(r_dim), tol(ncol), rss(ncol), vorder(ncol),  &
          row_ptr(ncol) )

d = zero
rhs = zero
r = zero
sserr = zero

if (fit_const) then
  do i = 1, ncol
    vorder(i) = i-1
  end do
else
  do i = 1, ncol
    vorder(i) = i
  end do
end if ! (fit_const)

! row_ptr(i) is the position of element r(i,i+1) in array r().

row_ptr(1) = 1
do i = 2, ncol-1
  row_ptr(i) = row_ptr(i-1) + ncol - i + 1
end do
row_ptr(ncol) = 0

initialized = .true.
tol_set = .false.
rss_set = .false.

return
end subroutine startup




subroutine includ(weight, xrow, yelem)

!     algorithm as75.1  appl. statist. (1974) vol.23, no. 3

!     calling this routine updates d, r, rhs and sserr by the
!     inclusion of xrow, yelem with the specified weight.

!     *** warning  array xrow is overwritten.

!     n.b. as this routine will be called many times in most applications,
!          checks have been eliminated.
!
!--------------------------------------------------------------------------


implicit none
real (dp),intent(in)                    :: weight, yelem
real (dp), dimension(:), intent(in out) :: xrow

!     local variables

integer     :: i, k, nextr
real (dp)   :: w, y, xi, di, wxi, dpi, cbar, sbar, xk

nobs = nobs + 1
w = weight
y = yelem
rss_set = .false.
nextr = 1
do i = 1, ncol

!     skip unnecessary transformations.   test on exact zeroes must be
!     used or stability can be destroyed.

  if (abs(w) < vsmall) return
  xi = xrow(i)
  if (abs(xi) < vsmall) then
    nextr = nextr + ncol - i
  else
    di = d(i)
    wxi = w * xi
    dpi = di + wxi*xi
    cbar = di / dpi
    sbar = wxi / dpi
    w = cbar * w
    d(i) = dpi
    do k = i+1, ncol
      xk = xrow(k)
      xrow(k) = xk - xi * r(nextr)
      r(nextr) = cbar * r(nextr) + sbar * xk
      nextr = nextr + 1
    end do
    xk = y
    y = xk - xi * rhs(i)
    rhs(i) = cbar * rhs(i) + sbar * xk
  end if
end do ! i = 1, ncol

!     y * sqrt(w) is now equal to the brown, durbin & evans recursive
!     residual.

sserr = sserr + w * y * y

return
end subroutine includ



subroutine regcf(beta, nreq, ifault)

!     algorithm as274  appl. statist. (1992) vol.41, no. 2

!     modified version of as75.4 to calculate regression coefficients
!     for the first nreq variables, given an orthogonal reduction from
!     as75.1.
!
!--------------------------------------------------------------------------

implicit none
integer, intent(in)                  :: nreq
integer, intent(out)                 :: ifault
real (dp), dimension(:), intent(out) :: beta

!     local variables

integer   :: i, j, nextr

!     some checks.

ifault = 0
if (nreq < 1 .or. nreq > ncol) ifault = ifault + 4
if (ifault /= 0) return

if (.not. tol_set) call tolset()

do i = nreq, 1, -1
  if (sqrt(d(i)) < tol(i)) then
    beta(i) = zero
    d(i) = zero
    ifault = -i
  else
    beta(i) = rhs(i)
    nextr = row_ptr(i)
    do j = i+1, nreq
      beta(i) = beta(i) - r(nextr) * beta(j)
      nextr = nextr + 1
    end do ! j = i+1, nreq
  end if
end do ! i = nreq, 1, -1

return
end subroutine regcf



subroutine tolset(eps)

!     algorithm as274  appl. statist. (1992) vol.41, no. 2

!     sets up array tol for testing for zeroes in an orthogonal
!     reduction formed using as75.1.

real (dp), intent(in), optional :: eps

!     unless the argument eps is set, it is assumed that the input data are
!     recorded to full machine accuracy.   this is often not the case.
!     if, for instance, the data are recorded to `single precision' of about
!     6-7 significant decimal digits, then singularities will not be detected.
!     it is suggested that in this case eps should be set equal to
!     10.0 * epsilon(1.0)
!     if the data are recorded to say 4 significant decimals, then eps should
!     be set to 1.0e-03
!     the above comments apply to the predictor variables, not to the
!     dependent variable.

!     correction - 19 august 2002
!     when negative weights are used, it is possible for an alement of d
!     to be negative.

!     local variables.
!
!--------------------------------------------------------------------------

!     local variables

integer    :: col, row, pos
real (dp)  :: eps1, ten = 10.0, total, work(ncol)

!     eps is a machine-dependent constant.

if (present(eps)) then
  eps1 = max(abs(eps), ten * epsilon(ten))
else
  eps1 = ten * epsilon(ten)
end if

!     set tol(i) = sum of absolute values in column i of r after
!     scaling each element by the square root of its row multiplier,
!     multiplied by eps1.

work = sqrt(abs(d))
do col = 1, ncol
  pos = col - 1
  total = work(col)
  do row = 1, col-1
    total = total + abs(r(pos)) * work(row)
    pos = pos + ncol - row - 1
  end do
  tol(col) = eps1 * total
end do

tol_set = .true.
return
end subroutine tolset




subroutine sing(lindep, ifault)

!     algorithm as274  appl. statist. (1992) vol.41, no. 2

!     checks for singularities, reports, and adjusts orthogonal
!     reductions produced by as75.1.

!     correction - 19 august 2002
!     when negative weights are used, it is possible for an alement of d
!     to be negative.

!     auxiliary routines called: includ, tolset
!
!--------------------------------------------------------------------------

integer, intent(out)                :: ifault
logical, dimension(:), intent(out)  :: lindep

!     local variables

real (dp)  :: temp, x(ncol), work(ncol), y, weight
integer    :: pos, row, pos2

ifault = 0

work = sqrt(abs(d))
if (.not. tol_set) call tolset()

do row = 1, ncol
  temp = tol(row)
  pos = row_ptr(row)         ! pos = location of first element in row

!     if diagonal element is near zero, set it to zero, set appropriate
!     element of lindep, and use includ to augment the projections in
!     the lower rows of the orthogonalization.

  lindep(row) = .false.
  if (work(row) <= temp) then
    lindep(row) = .true.
    ifault = ifault - 1
    if (row < ncol) then
      pos2 = pos + ncol - row - 1
      x = zero
      x(row+1:ncol) = r(pos:pos2)
      y = rhs(row)
      weight = d(row)
      r(pos:pos2) = zero
      d(row) = zero
      rhs(row) = zero
      call includ(weight, x, y)
                             ! includ automatically increases the number
                             ! of cases each time it is called.
      nobs = nobs - 1
    else
      sserr = sserr + d(row) * rhs(row)**2
    end if ! (row < ncol)
  end if ! (work(row) <= temp)
end do ! row = 1, ncol

return
end subroutine sing



subroutine ss()

!     algorithm as274  appl. statist. (1992) vol.41, no. 2

!     calculates partial residual sums of squares from an orthogonal
!     reduction from as75.1.
!
!--------------------------------------------------------------------------

!     local variables

integer    :: i
real (dp)  :: total

total = sserr
rss(ncol) = sserr
do i = ncol, 2, -1
  total = total + d(i) * rhs(i)**2
  rss(i-1) = total
end do

rss_set = .true.
return
end subroutine ss



subroutine cov(nreq, var, covmat, dimcov, sterr, ifault)

!     algorithm as274  appl. statist. (1992) vol.41, no. 2

!     calculate covariance matrix for regression coefficients for the
!     first nreq variables, from an orthogonal reduction produced from
!     as75.1.

!     auxiliary routine called: inv
!
!--------------------------------------------------------------------------

integer, intent(in)                   :: nreq, dimcov
integer, intent(out)                  :: ifault
real (dp), intent(out)                :: var
real (dp), dimension(:), intent(out)  :: covmat, sterr

!     local variables.

integer                :: dim_rinv, pos, row, start, pos2, col, pos1, k
real (dp)              :: total
real (dp), allocatable :: rinv(:)

!     check that dimension of array covmat is adequate.

if (dimcov < nreq*(nreq+1)/2) then
  ifault = 1
  return
end if

!     check for small or zero multipliers on the diagonal.

ifault = 0
do row = 1, nreq
  if (abs(d(row)) < vsmall) ifault = -row
end do
if (ifault /= 0) return

!     calculate estimate of the residual variance.

if (nobs > nreq) then
  if (.not. rss_set) call ss()
  var = rss(nreq) / (nobs - nreq)
else
  ifault = 2
  return
end if

dim_rinv = nreq*(nreq-1)/2
allocate ( rinv(dim_rinv) )

call inv(nreq, rinv)
pos = 1
start = 1
do row = 1, nreq
  pos2 = start
  do col = row, nreq
    pos1 = start + col - row
    if (row == col) then
      total = one / d(col)
    else
      total = rinv(pos1-1) / d(col)
    end if
    do k = col+1, nreq
      total = total + rinv(pos1) * rinv(pos2) / d(k)
      pos1 = pos1 + 1
      pos2 = pos2 + 1
    end do ! k = col+1, nreq
    covmat(pos) = total * var
    if (row == col) sterr(row) = sqrt(covmat(pos))
    pos = pos + 1
  end do ! col = row, nreq
  start = start + nreq - row
end do ! row = 1, nreq

deallocate(rinv)
return
end subroutine cov



subroutine inv(nreq, rinv)

!     algorithm as274  appl. statist. (1992) vol.41, no. 2

!     invert first nreq rows and columns of cholesky factorization
!     produced by as 75.1.
!
!--------------------------------------------------------------------------

integer, intent(in)                  :: nreq
real (dp), dimension(:), intent(out) :: rinv

!     local variables.

integer    :: pos, row, col, start, k, pos1, pos2
real (dp)  :: total

!     invert r ignoring row multipliers, from the bottom up.

pos = nreq * (nreq-1)/2
do row = nreq-1, 1, -1
  start = row_ptr(row)
  do col = nreq, row+1, -1
    pos1 = start
    pos2 = pos
    total = zero
    do k = row+1, col-1
      pos2 = pos2 + nreq - k
      total = total - r(pos1) * rinv(pos2)
      pos1 = pos1 + 1
    end do ! k = row+1, col-1
    rinv(pos) = total - r(pos1)
    pos = pos - 1
  end do ! col = nreq, row+1, -1
end do ! row = nreq-1, 1, -1

return
end subroutine inv



subroutine partial_corr(in, cormat, dimc, ycorr, ifault)

!     replaces subroutines pcorr and cor of:
!     algorithm as274  appl. statist. (1992) vol.41, no. 2

!     calculate partial correlations after the variables in rows
!     1, 2, ..., in have been forced into the regression.
!     if in = 1, and the first row of r represents a constant in the
!     model, then the usual simple correlations are returned.

!     if in = 0, the value returned in array cormat for the correlation
!     of variables xi & xj is:
!       sum ( xi.xj ) / sqrt ( sum (xi^2) . sum (xj^2) )

!     on return, array cormat contains the upper triangle of the matrix of
!     partial correlations stored by rows, excluding the 1's on the diagonal.
!     e.g. if in = 2, the consecutive elements returned are:
!     (3,4) (3,5) ... (3,ncol), (4,5) (4,6) ... (4,ncol), etc.
!     array ycorr stores the partial correlations with the y-variable
!     starting with ycorr(in+1) = partial correlation with the variable in
!     position (in+1).
!
!--------------------------------------------------------------------------

integer, intent(in)                  :: in, dimc
integer, intent(out)                 :: ifault
real (dp), dimension(:), intent(out) :: cormat, ycorr

!     local variables.

integer    :: base_pos, pos, row, col, col1, col2, pos1, pos2
real (dp)  :: rms(in+1:ncol), sumxx, sumxy, sumyy, work(in+1:ncol)

!     some checks.

ifault = 0
if (in < 0 .or. in > ncol-1) ifault = ifault + 4
if (dimc < (ncol-in)*(ncol-in-1)/2) ifault = ifault + 8
if (ifault /= 0) return

!     base position for calculating positions of elements in row (in+1) of r.

base_pos = in*ncol - (in+1)*(in+2)/2

!     calculate 1/rms of elements in columns from in to (ncol-1).

if (d(in+1) > zero) rms(in+1) = one / sqrt(d(in+1))
do col = in+2, ncol
  pos = base_pos + col
  sumxx = d(col)
  do row = in+1, col-1
    sumxx = sumxx + d(row) * r(pos)**2
    pos = pos + ncol - row - 1
  end do ! row = in+1, col-1
  if (sumxx > zero) then
    rms(col) = one / sqrt(sumxx)
  else
    rms(col) = zero
    ifault = -col
  end if ! (sumxx > zero)
end do ! col = in+1, ncol-1

!     calculate 1/rms for the y-variable

sumyy = sserr
do row = in+1, ncol
  sumyy = sumyy + d(row) * rhs(row)**2
end do ! row = in+1, ncol
if (sumyy > zero) sumyy = one / sqrt(sumyy)

!     calculate sums of cross-products.
!     these are obtained by taking dot products of pairs of columns of r,
!     but with the product for each row multiplied by the row multiplier
!     in array d.

pos = 1
do col1 = in+1, ncol
  sumxy = zero
  work(col1+1:ncol) = zero
  pos1 = base_pos + col1
  do row = in+1, col1-1
    pos2 = pos1 + 1
    do col2 = col1+1, ncol
      work(col2) = work(col2) + d(row) * r(pos1) * r(pos2)
      pos2 = pos2 + 1
    end do ! col2 = col1+1, ncol
    sumxy = sumxy + d(row) * r(pos1) * rhs(row)
    pos1 = pos1 + ncol - row - 1
  end do ! row = in+1, col1-1

!     row col1 has an implicit 1 as its first element (in column col1)

  pos2 = pos1 + 1
  do col2 = col1+1, ncol
    work(col2) = work(col2) + d(col1) * r(pos2)
    pos2 = pos2 + 1
    cormat(pos) = work(col2) * rms(col1) * rms(col2)
    pos = pos + 1
  end do ! col2 = col1+1, ncol
  sumxy = sumxy + d(col1) * rhs(col1)
  ycorr(col1) = sumxy * rms(col1) * sumyy
end do ! col1 = in+1, ncol-1

ycorr(1:in) = zero

return
end subroutine partial_corr




subroutine vmove(from, to, ifault)

!     algorithm as274 appl. statist. (1992) vol.41, no. 2

!     move variable from position from to position to in an
!     orthogonal reduction produced by as75.1.
!
!--------------------------------------------------------------------------

integer, intent(in)    :: from, to
integer, intent(out)   :: ifault

!     local variables

real (dp)  :: d1, d2, x, d1new, d2new, cbar, sbar, y
integer    :: m, first, last, inc, m1, m2, mp1, col, pos, row

!     check input parameters

ifault = 0
if (from < 1 .or. from > ncol) ifault = ifault + 4
if (to < 1 .or. to > ncol) ifault = ifault + 8
if (ifault /= 0) return

if (from == to) return

if (.not. rss_set) call ss()

if (from < to) then
  first = from
  last = to - 1
  inc = 1
else
  first = from - 1
  last = to
  inc = -1
end if

do m = first, last, inc

!     find addresses of first elements of r in rows m and (m+1).

  m1 = row_ptr(m)
  m2 = row_ptr(m+1)
  mp1 = m + 1
  d1 = d(m)
  d2 = d(mp1)

!     special cases.

  if (d1 < vsmall .and. d2 < vsmall) go to 40
  x = r(m1)
  if (abs(x) * sqrt(d1) < tol(mp1)) then
    x = zero
  end if
  if (d1 < vsmall .or. abs(x) < vsmall) then
    d(m) = d2
    d(mp1) = d1
    r(m1) = zero
    do col = m+2, ncol
      m1 = m1 + 1
      x = r(m1)
      r(m1) = r(m2)
      r(m2) = x
      m2 = m2 + 1
    end do ! col = m+2, ncol
    x = rhs(m)
    rhs(m) = rhs(mp1)
    rhs(mp1) = x
    go to 40
  else if (d2 < vsmall) then
    d(m) = d1 * x**2
    r(m1) = one / x
    r(m1+1:m1+ncol-m-1) = r(m1+1:m1+ncol-m-1) / x
    rhs(m) = rhs(m) / x
    go to 40
  end if

!     planar rotation in regular case.

  d1new = d2 + d1*x**2
  cbar = d2 / d1new
  sbar = x * d1 / d1new
  d2new = d1 * cbar
  d(m) = d1new
  d(mp1) = d2new
  r(m1) = sbar
  do col = m+2, ncol
    m1 = m1 + 1
    y = r(m1)
    r(m1) = cbar*r(m2) + sbar*y
    r(m2) = y - x*r(m2)
    m2 = m2 + 1
  end do ! col = m+2, ncol
  y = rhs(m)
  rhs(m) = cbar*rhs(mp1) + sbar*y
  rhs(mp1) = y - x*rhs(mp1)

!     swap columns m and (m+1) down to row (m-1).

  40 pos = m
  do row = 1, m-1
    x = r(pos)
    r(pos) = r(pos-1)
    r(pos-1) = x
    pos = pos + ncol - row - 1
  end do ! row = 1, m-1

!     adjust variable order (vorder), the tolerances (tol) and
!     the vector of residual sums of squares (rss).

  m1 = vorder(m)
  vorder(m) = vorder(mp1)
  vorder(mp1) = m1
  x = tol(m)
  tol(m) = tol(mp1)
  tol(mp1) = x
  rss(m) = rss(mp1) + d(mp1) * rhs(mp1)**2
end do

return
end subroutine vmove



subroutine reordr(list, n, pos1, ifault)

!     algorithm as274  appl. statist. (1992) vol.41, no. 2

!     re-order the variables in an orthogonal reduction produced by
!     as75.1 so that the n variables in list start at position pos1,
!     though will not necessarily be in the same order as in list.
!     any variables in vorder before position pos1 are not moved.

!     auxiliary routine called: vmove
!
!--------------------------------------------------------------------------

integer, intent(in)               :: n, pos1
integer, dimension(:), intent(in) :: list
integer, intent(out)              :: ifault

!     local variables.

integer    :: next, i, l, j

!     check n.

ifault = 0
if (n < 1 .or. n > ncol+1-pos1) ifault = ifault + 4
if (ifault /= 0) return

!     work through vorder finding variables which are in list.

next = pos1
i = pos1
10 l = vorder(i)
do j = 1, n
  if (l == list(j)) go to 40
end do
30 i = i + 1
if (i <= ncol) go to 10

!     if this point is reached, one or more variables in list has not
!     been found.

ifault = 8
return

!     variable l is in list; move it up to position next if it is not
!     already there.

40 if (i > next) call vmove(i, next, ifault)
next = next + 1
if (next < n+pos1) go to 30

return
end subroutine reordr



subroutine hdiag(xrow, nreq, hii, ifault)

!     algorithm as274  appl. statist. (1992) vol.41, no. 2
!
!                         -1           -1
! the hat matrix h = x(x'x) x' = x(r'dr) x' = z'dz

!              -1
! where z = x'r

! here we only calculate the diagonal element hii corresponding to one
! row (xrow).   the variance of the i-th least-squares residual is (1 - hii).
!--------------------------------------------------------------------------

integer, intent(in)                  :: nreq
integer, intent(out)                 :: ifault
real (dp), dimension(:), intent(in)  :: xrow
real (dp), intent(out)               :: hii

!     local variables

integer    :: col, row, pos
real (dp)  :: total, wk(ncol)

!     some checks

ifault = 0
if (nreq > ncol) ifault = ifault + 4
if (ifault /= 0) return

!     the elements of xrow.inv(r).sqrt(d) are calculated and stored in wk.

hii = zero
do col = 1, nreq
  if (sqrt(d(col)) <= tol(col)) then
    wk(col) = zero
  else
    pos = col - 1
    total = xrow(col)
    do row = 1, col-1
      total = total - wk(row)*r(pos)
      pos = pos + ncol - row - 1
    end do ! row = 1, col-1
    wk(col) = total
    hii = hii + total**2 / d(col)
  end if
end do ! col = 1, nreq

return
end subroutine hdiag



function varprd(x, nreq) result(fn_val)

!     calculate the variance of x'b where b consists of the first nreq
!     least-squares regression coefficients.
!
!--------------------------------------------------------------------------

integer, intent(in)                  :: nreq
real (dp), dimension(:), intent(in)  :: x
real (dp)                            :: fn_val

!     local variables

integer    :: ifault, row
real (dp)  :: var, wk(nreq)

!     check input parameter values

fn_val = zero
ifault = 0
if (nreq < 1 .or. nreq > ncol) ifault = ifault + 4
if (nobs <= nreq) ifault = ifault + 8
if (ifault /= 0) then
  write(*, '(1x, a, i4)') 'error in function varprd: ifault =', ifault
  return
end if

!     calculate the residual variance estimate.

var = sserr / (nobs - nreq)

!     variance of x'b = var.x'(inv r)(inv d)(inv r')x
!     first call bksub2 to calculate (inv r')x by back-substitution.

call bksub2(x, wk, nreq)
do row = 1, nreq
  if(d(row) > tol(row)) fn_val = fn_val + wk(row)**2 / d(row)
end do

fn_val = fn_val * var

return
end function varprd



subroutine bksub2(x, b, nreq)

!     solve x = r'b for b given x, using only the first nreq rows and
!     columns of r, and only the first nreq elements of r.
!
!--------------------------------------------------------------------------

integer, intent(in)                  :: nreq
real (dp), dimension(:), intent(in)  :: x
real (dp), dimension(:), intent(out) :: b

!     local variables

integer    :: pos, row, col
real (dp)  :: temp

!     solve by back-substitution, starting from the top.

do row = 1, nreq
  pos = row - 1
  temp = x(row)
  do col = 1, row-1
    temp = temp - r(pos)*b(col)
    pos = pos + ncol - col - 1
  end do
  b(row) = temp
end do

return
end subroutine bksub2


end module leaps_lsq
