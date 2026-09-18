module grpreg_spline
   use grpreg_kinds, only : dp
   use grpreg_types, only : spline_expansion_type
   use grpreg_linalg, only : symmetric_eigen
   implicit none
   private
   public :: expand_spline, predict_spline

contains

   subroutine expand_spline(x, expansion, df, degree, spline_type)
      real(dp), intent(in) :: x(:,:) !! Original feature matrix, observations by variables, to expand into spline groups.
      type(spline_expansion_type), intent(out) :: expansion !! Expanded basis matrix plus knots and boundaries for prediction.
      integer, optional, intent(in) :: df !! Number of basis columns per original variable; default three.
      integer, optional, intent(in) :: degree !! Polynomial degree; natural splines force degree three.
      character(len=*), optional, intent(in) :: spline_type !! Spline type, bs or ns; default ns.
      integer :: d, deg, p, n, nk, maxk, j, k
      character(len=2) :: typ
      real(dp), allocatable :: probs(:), basis(:,:)
      d = 3
      if (present(df)) d = df
      deg = 3
      if (present(degree)) deg = degree
      typ = 'ns'
      if (present(spline_type)) typ = trim(spline_type)
      if (typ == 'ns') deg = 3
      if (d < 1) error stop 'expand_spline: df must be positive'
      if (typ == 'bs' .and. d < deg) error stop 'expand_spline: B-spline df must be at least degree'
      if (typ /= 'bs' .and. typ /= 'ns') error stop 'expand_spline: spline_type must be bs or ns'
      n = size(x,1)
      p = size(x,2)
      nk = merge(max(d-deg,0),max(d-1,0),typ == 'bs')
      maxk = max(1,nk)
      allocate(expansion%x(n,d*p),expansion%group(d*p),expansion%knots(maxk,p),expansion%boundary(2,p))
      expansion%knots = 0.0_dp
      expansion%degree = deg
      expansion%df = d
      expansion%spline_type = typ
      do j = 1, p
         expansion%boundary(1,j) = minval(x(:,j))
         expansion%boundary(2,j) = maxval(x(:,j))
         if (nk > 0) then
            allocate(probs(nk))
            do k = 1, nk
               probs(k) = real(k,dp)/real(nk+1,dp)
               expansion%knots(k,j) = quantile_type7(x(:,j),probs(k))
            end do
            deallocate(probs)
         end if
         if (typ == 'bs') then
            call bspline_matrix(x(:,j),expansion%knots(1:nk,j),expansion%boundary(:,j),deg,.false.,basis)
         else
            call natural_spline_matrix(x(:,j),expansion%knots(1:nk,j),expansion%boundary(:,j),basis)
         end if
         expansion%x(:,(j-1)*d+1:j*d) = basis
         expansion%group((j-1)*d+1:j*d) = j
      end do
   end subroutine expand_spline

   subroutine predict_spline(expansion, x, expanded)
      type(spline_expansion_type), intent(in) :: expansion !! Training expansion containing knots, boundaries, degree, and type.
      real(dp), intent(in) :: x(:,:) !! New original-scale feature matrix with the same number of variables as training data.
      real(dp), allocatable, intent(out) :: expanded(:,:) !! New spline basis matrix compatible with expansion%x grouping.
      real(dp), allocatable :: basis(:,:)
      integer :: p, d, nk, j
      p = size(expansion%boundary,2)
      d = expansion%df
      if (size(x,2) /= p) error stop 'predict_spline: number of variables differs from training expansion'
      allocate(expanded(size(x,1),d*p))
      if (expansion%spline_type == 'bs') then
         nk = max(d-expansion%degree,0)
      else
         nk = max(d-1,0)
      end if
      do j = 1, p
         if (expansion%spline_type == 'bs') then
            call bspline_matrix(x(:,j),expansion%knots(1:nk,j),expansion%boundary(:,j), &
               expansion%degree,.false.,basis)
         else
            call natural_spline_matrix(x(:,j),expansion%knots(1:nk,j),expansion%boundary(:,j),basis)
         end if
         expanded(:,(j-1)*d+1:j*d) = basis
      end do
   end subroutine predict_spline

   subroutine bspline_matrix(x, interior, boundary, degree, intercept, basis)
      real(dp), intent(in) :: x(:) !! Evaluation coordinates for one spline variable.
      real(dp), intent(in) :: interior(:) !! Strictly interior spline knots in ascending order.
      real(dp), intent(in) :: boundary(:) !! Lower and upper boundary knots.
      integer, intent(in) :: degree !! Nonnegative polynomial degree of the B-spline basis.
      logical, intent(in) :: intercept !! True retains the leading basis; false drops it like splines::bs default.
      real(dp), allocatable, intent(out) :: basis(:,:) !! B-spline basis evaluated at x.
      real(dp), allocatable :: knots(:), allb(:)
      integer :: ord, nbas, outcols, i, pos
      ord = degree+1
      allocate(knots(2*ord+size(interior)))
      knots(1:ord) = boundary(1)
      if (size(interior) > 0) knots(ord+1:ord+size(interior)) = interior
      pos = ord+size(interior)+1
      knots(pos:) = boundary(2)
      nbas = size(knots)-ord
      outcols = nbas-merge(0,1,intercept)
      allocate(basis(size(x),outcols))
      do i = 1, size(x)
         call bspline_all(x(i),knots,degree,allb)
         if (intercept) then
            basis(i,:) = allb
         else
            basis(i,:) = allb(2:)
         end if
      end do
   end subroutine bspline_matrix

   subroutine natural_spline_matrix(x, interior, boundary, basis)
      real(dp), intent(in) :: x(:) !! Evaluation coordinates for one natural cubic spline variable.
      real(dp), intent(in) :: interior(:) !! Interior knots chosen from sample quantiles.
      real(dp), intent(in) :: boundary(:) !! Lower and upper natural-spline boundary knots.
      real(dp), allocatable, intent(out) :: basis(:,:) !! Natural cubic spline basis spanning the R ns() space.
      real(dp), allocatable :: knots(:), raw(:,:), row(:), cmat(:,:), gram(:,:), eval(:), evec(:,:), nullspace(:,:)
      real(dp) :: h
      integer :: ord, qfull, q, d, i, info
      ord = 4
      allocate(knots(2*ord+size(interior)))
      knots(1:ord) = boundary(1)
      if (size(interior) > 0) knots(ord+1:ord+size(interior)) = interior
      knots(ord+size(interior)+1:) = boundary(2)
      qfull = size(knots)-ord
      q = qfull-1
      d = q-2
      allocate(raw(size(x),q),cmat(2,q),gram(q,q))
      do i = 1, size(x)
         call bspline_all(x(i),knots,3,row)
         raw(i,:) = row(2:)
      end do
      h = max((boundary(2)-boundary(1))*1.0e-5_dp,1.0e-7_dp)
      call boundary_second_derivative(boundary(1),h,knots,.true.,row)
      cmat(1,:) = row(2:)
      call boundary_second_derivative(boundary(2),h,knots,.false.,row)
      cmat(2,:) = row(2:)
      gram = matmul(transpose(cmat),cmat)
      call symmetric_eigen(gram,eval,evec,info)
      if (info /= 0) error stop 'expand_spline: natural-spline constraint eigensystem failed'
      allocate(nullspace(q,d))
      nullspace = evec(:,q-d+1:q)
      allocate(basis(size(x),d))
      basis = matmul(raw,nullspace)
   end subroutine natural_spline_matrix

   subroutine boundary_second_derivative(x0, h, knots, lower, second)
      real(dp), intent(in) :: x0 !! Boundary coordinate at which the second derivative is required.
      real(dp), intent(in) :: h !! Positive one-sided finite-difference spacing.
      real(dp), intent(in) :: knots(:) !! Augmented cubic B-spline knot vector.
      logical, intent(in) :: lower !! True uses a forward stencil; false uses a backward stencil.
      real(dp), allocatable, intent(out) :: second(:) !! Second derivative of every cubic B-spline basis column.
      real(dp), allocatable :: f0(:), f1(:), f2(:), f3(:)
      if (lower) then
         call bspline_all(x0,knots,3,f0)
         call bspline_all(x0+h,knots,3,f1)
         call bspline_all(x0+2.0_dp*h,knots,3,f2)
         call bspline_all(x0+3.0_dp*h,knots,3,f3)
      else
         call bspline_all(x0,knots,3,f0)
         call bspline_all(x0-h,knots,3,f1)
         call bspline_all(x0-2.0_dp*h,knots,3,f2)
         call bspline_all(x0-3.0_dp*h,knots,3,f3)
      end if
      allocate(second(size(f0)))
      second = (2.0_dp*f0-5.0_dp*f1+4.0_dp*f2-f3)/(h*h)
   end subroutine boundary_second_derivative

   subroutine bspline_all(x, knots, degree, basis)
      real(dp), intent(in) :: x !! Scalar coordinate at which all B-spline basis functions are evaluated.
      real(dp), intent(in) :: knots(:) !! Nondecreasing augmented knot vector.
      integer, intent(in) :: degree !! Polynomial degree of the desired B-spline basis.
      real(dp), allocatable, intent(out) :: basis(:) !! Values of all degree-p basis functions at x.
      real(dp), allocatable :: prev(:), curr(:)
      real(dp) :: left, right
      integer :: i, p, n0, nbas
      n0 = size(knots)-1
      allocate(prev(n0))
      prev = 0.0_dp
      if (x <= knots(1)) then
         prev(1) = 1.0_dp
      else if (x >= knots(size(knots))) then
         prev(n0) = 1.0_dp
      else
         do i = 1, n0
            if (x >= knots(i) .and. x < knots(i+1)) prev(i) = 1.0_dp
         end do
      end if
      do p = 1, degree
         allocate(curr(n0-p))
         curr = 0.0_dp
         do i = 1, n0-p
            left = 0.0_dp
            right = 0.0_dp
            if (knots(i+p)-knots(i) > tiny(1.0_dp)) &
               left = (x-knots(i))/(knots(i+p)-knots(i))*prev(i)
            if (knots(i+p+1)-knots(i+1) > tiny(1.0_dp)) &
               right = (knots(i+p+1)-x)/(knots(i+p+1)-knots(i+1))*prev(i+1)
            curr(i) = left+right
         end do
         call move_alloc(curr,prev)
      end do
      nbas = size(knots)-degree-1
      allocate(basis(nbas))
      basis = prev(1:nbas)
      if (x >= knots(size(knots))) then
         basis = 0.0_dp
         basis(nbas) = 1.0_dp
      end if
   end subroutine bspline_all

   real(dp) function quantile_type7(x, probability) result(value)
      real(dp), intent(in) :: x(:) !! Sample values whose R type-7 quantile is requested.
      real(dp), intent(in) :: probability !! Quantile probability in the closed interval [0,1].
      real(dp), allocatable :: work(:)
      real(dp) :: h, frac
      integer :: lo, hi
      allocate(work(size(x)))
      work = x
      call sort_real(work)
      h = 1.0_dp+(real(size(x)-1,dp))*min(max(probability,0.0_dp),1.0_dp)
      lo = max(1,min(size(x),int(floor(h))))
      hi = max(1,min(size(x),lo+1))
      frac = h-real(lo,dp)
      value = (1.0_dp-frac)*work(lo)+frac*work(hi)
   end function quantile_type7

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:) !! Real vector sorted increasingly in place.
      real(dp) :: key
      integer :: i, j
      do i = 2, size(x)
         key = x(i)
         j = i-1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j)
            j = j-1
         end do
         x(j+1) = key
      end do
   end subroutine sort_real

end module grpreg_spline
