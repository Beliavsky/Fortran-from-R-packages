! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_ph
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use r_compat, only: dp
   use matrixdist_linalg
   use matrixdist_types, only: ph_type
   implicit none
   private
   public :: make_ph, ph_density, ph_cdf, ph_survival, ph_hazard, ph_quantile
   public :: ph_moment, ph_mean, ph_variance, ph_laplace, ph_mgf, ph_loglik
   public :: ph_sum, ph_mixture, ph_minimum, ph_maximum, ph_exit_rates

contains

   function make_ph(alpha, s) result(x)
      real(dp), intent(in) :: alpha(:), s(:,:)
      type(ph_type) :: x
      if (size(s,1) /= size(s,2) .or. size(alpha) /= size(s,1)) &
         error stop "make_ph: incompatible dimensions"
      x%alpha = alpha
      x%s = s
   end function make_ph

   function ph_exit_rates(s) result(t)
      real(dp), intent(in) :: s(:,:)
      real(dp) :: t(size(s,1))
      t = -sum(s,dim=2)
   end function ph_exit_rates

   function ph_density(x, alpha, s) result(f)
      real(dp), intent(in) :: x, alpha(:), s(:,:)
      real(dp) :: f
      real(dp), allocatable :: e(:), t(:), m(:,:)
      integer :: p
      p = size(alpha)
      if (x < 0.0_dp) then
         f = 0.0_dp
         return
      end if
      allocate(e(p))
      e = 1.0_dp
      if (x == 0.0_dp) then
         f = max(0.0_dp, 1.0_dp - sum(alpha))
      else
         t = ph_exit_rates(s)
         m = matrix_exponential(s*x)
         f = dot_product(alpha, matmul(m,t))
         if (f < 0.0_dp .and. abs(f) < 100.0_dp*epsilon(f)) f = 0.0_dp
      end if
   end function ph_density

   function ph_survival(x, alpha, s) result(q)
      real(dp), intent(in) :: x, alpha(:), s(:,:)
      real(dp) :: q
      real(dp), allocatable :: m(:,:), e(:)
      if (x < 0.0_dp) then
         q = 1.0_dp
         return
      end if
      allocate(e(size(alpha)))
      e = 1.0_dp
      if (x == 0.0_dp) then
         q = sum(alpha)
      else
         m = matrix_exponential(s*x)
         q = dot_product(alpha, matmul(m,e))
      end if
      q = min(1.0_dp,max(0.0_dp,q))
   end function ph_survival

   function ph_cdf(x, alpha, s, lower_tail) result(pv)
      real(dp), intent(in) :: x, alpha(:), s(:,:)
      logical, intent(in), optional :: lower_tail
      real(dp) :: pv, q
      logical :: lower
      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      q = ph_survival(x,alpha,s)
      if (lower) then
         pv = 1.0_dp-q
      else
         pv = q
      end if
   end function ph_cdf

   function ph_hazard(x, alpha, s) result(h)
      real(dp), intent(in) :: x, alpha(:), s(:,:)
      real(dp) :: h, q
      q = ph_survival(x,alpha,s)
      if (q <= 0.0_dp) then
         h = huge(1.0_dp)
      else
         h = ph_density(x,alpha,s)/q
      end if
   end function ph_hazard

   function ph_quantile(prob, alpha, s, tol, maxit) result(q)
      real(dp), intent(in) :: prob, alpha(:), s(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      real(dp) :: q, lo, hi, mid, eps
      integer :: it, nmax
      if (prob <= ph_cdf(0.0_dp,alpha,s)) then
         q = 0.0_dp
         return
      end if
      if (prob >= 1.0_dp) then
         q = huge(1.0_dp)
         return
      end if
      if (prob < 0.0_dp) then
         q = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      eps = 1.0e-10_dp
      if (present(tol)) eps = tol
      nmax = 200
      if (present(maxit)) nmax = maxit
      lo = 0.0_dp
      hi = 1.0_dp
      do while (ph_cdf(hi,alpha,s) < prob .and. hi < 0.25_dp*huge(hi))
         hi = 2.0_dp*hi
      end do
      do it = 1, nmax
         mid = 0.5_dp*(lo+hi)
         if (ph_cdf(mid,alpha,s) < prob) then
            lo = mid
         else
            hi = mid
         end if
         if (hi-lo <= eps*max(1.0_dp,mid)) exit
      end do
      q = 0.5_dp*(lo+hi)
   end function ph_quantile

   function ph_moment(k, alpha, s) result(mom)
      integer, intent(in) :: k
      real(dp), intent(in) :: alpha(:), s(:,:)
      real(dp) :: mom
      real(dp), allocatable :: u(:,:), uk(:,:), e(:)
      integer :: j
      real(dp) :: fact
      if (k <= 0) error stop "ph_moment: k must be positive"
      u = matrix_inverse(-s)
      uk = matrix_power(k,u)
      allocate(e(size(alpha)))
      e=1.0_dp
      fact = 1.0_dp
      do j=2,k
         fact=fact*real(j,dp)
      end do
      mom = fact*dot_product(alpha,matmul(uk,e))
   end function ph_moment

   function ph_mean(alpha,s) result(m)
      real(dp), intent(in) :: alpha(:), s(:,:)
      real(dp) :: m
      m = ph_moment(1,alpha,s)
   end function ph_mean

   function ph_variance(alpha,s) result(v)
      real(dp), intent(in) :: alpha(:), s(:,:)
      real(dp) :: v, m
      m=ph_mean(alpha,s)
      v=ph_moment(2,alpha,s)-m*m
   end function ph_variance

   function ph_laplace(r, alpha, s) result(v)
      real(dp), intent(in) :: r, alpha(:), s(:,:)
      real(dp) :: v
      real(dp), allocatable :: a(:,:), t(:), z(:)
      a = r*eye_matrix(size(alpha)) - s
      t = ph_exit_rates(s)
      z = solve_vector(a,t)
      v = dot_product(alpha,z)
   end function ph_laplace

   function ph_mgf(r, alpha, s) result(v)
      real(dp), intent(in) :: r, alpha(:), s(:,:)
      real(dp) :: v
      v = ph_laplace(-r,alpha,s)
   end function ph_mgf

   function ph_loglik(alpha,s,obs,weight,rcens,rcweight) result(ll)
      real(dp), intent(in) :: alpha(:), s(:,:), obs(:)
      real(dp), intent(in), optional :: weight(:), rcens(:), rcweight(:)
      real(dp) :: ll, f, w
      integer :: i
      ll=0.0_dp
      do i=1,size(obs)
         w=1.0_dp
         if (present(weight)) w=weight(i)
         f=ph_density(obs(i),alpha,s)
         if (f <= 0.0_dp) then
            ll=-huge(1.0_dp)
            return
         end if
         ll=ll+w*log(f)
      end do
      if (present(rcens)) then
         do i=1,size(rcens)
            w=1.0_dp
            if (present(rcweight)) w=rcweight(i)
            f=ph_survival(rcens(i),alpha,s)
            if (f <= 0.0_dp) then
               ll=-huge(1.0_dp)
               return
            end if
            ll=ll+w*log(f)
         end do
      end if
   end function ph_loglik

   function ph_sum(x1,x2) result(z)
      type(ph_type), intent(in) :: x1,x2
      type(ph_type) :: z
      integer :: p1,p2
      real(dp), allocatable :: t1(:)
      p1=size(x1%alpha)
      p2=size(x2%alpha)
      allocate(z%alpha(p1+p2), z%s(p1+p2,p1+p2))
      z%s=0.0_dp
      z%alpha(1:p1)=x1%alpha
      z%alpha(p1+1:p1+p2)=0.0_dp
      z%s(1:p1,1:p1)=x1%s
      z%s(p1+1:p1+p2,p1+1:p1+p2)=x2%s
      t1=ph_exit_rates(x1%s)
      z%s(1:p1,p1+1:p1+p2)=spread(t1,2,p2)*spread(x2%alpha,1,p1)
   end function ph_sum

   function ph_mixture(x1,x2,prob) result(z)
      type(ph_type), intent(in) :: x1,x2
      real(dp), intent(in) :: prob
      type(ph_type) :: z
      integer :: p1,p2
      p1=size(x1%alpha)
      p2=size(x2%alpha)
      allocate(z%alpha(p1+p2),z%s(p1+p2,p1+p2))
      z%s=0.0_dp
      z%alpha=[prob*x1%alpha,(1.0_dp-prob)*x2%alpha]
      z%s(1:p1,1:p1)=x1%s
      z%s(p1+1:p1+p2,p1+1:p1+p2)=x2%s
   end function ph_mixture

   function ph_minimum(x1,x2) result(z)
      type(ph_type), intent(in) :: x1,x2
      type(ph_type) :: z
      real(dp), allocatable :: a1(:,:), a2(:,:), aa(:,:)
      integer :: p1,p2
      p1=size(x1%alpha)
      p2=size(x2%alpha)
      allocate(a1(1,p1),a2(1,p2))
      a1(1,:)=x1%alpha
      a2(1,:)=x2%alpha
      aa=kronecker(a1,a2)
      z%alpha=aa(1,:)
      z%s=kronecker_sum(x1%s,x2%s)
   end function ph_minimum

   function ph_maximum(x1,x2) result(z)
      type(ph_type), intent(in) :: x1,x2
      type(ph_type) :: z
      integer :: p1,p2,n12,n
      real(dp), allocatable :: a1(:,:),a2(:,:),aa(:,:),t1(:),t2(:)
      p1=size(x1%alpha)
      p2=size(x2%alpha)
      n12=p1*p2
      n=n12+p1+p2
      allocate(a1(1,p1),a2(1,p2))
      a1(1,:)=x1%alpha
      a2(1,:)=x2%alpha
      aa=kronecker(a1,a2)
      allocate(z%alpha(n),z%s(n,n))
      z%alpha=0.0_dp
      z%s=0.0_dp
      z%alpha(1:n12)=aa(1,:)
      z%s(1:n12,1:n12)=kronecker_sum(x1%s,x2%s)
      t1=ph_exit_rates(x1%s)
      t2=ph_exit_rates(x2%s)
      z%s(1:n12,n12+1:n12+p1)=kronecker(eye_matrix(p1),reshape(t2,[p2,1]))
      z%s(1:n12,n12+p1+1:n)=kronecker(reshape(t1,[p1,1]),eye_matrix(p2))
      z%s(n12+1:n12+p1,n12+1:n12+p1)=x1%s
      z%s(n12+p1+1:n,n12+p1+1:n)=x2%s
   end function ph_maximum

end module matrixdist_ph
