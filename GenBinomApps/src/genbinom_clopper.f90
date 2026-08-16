module genbinom_clopper
   use genbinom_kinds, only : dp, i64
   use genbinom_special, only : beta_cdf, beta_quantile
   use genbinom_distribution, only : gbinom_pmf_table
   implicit none
   private

   type, public :: confidence_interval
      character(len=12) :: kind = "upper"
      real(dp) :: lower = 0.0_dp
      real(dp) :: upper = 1.0_dp
      real(dp) :: alpha = 0.1_dp
      integer :: status = 0
   end type confidence_interval

   public :: clopper_pearson_ci, cm_clopper_pearson_ci
   public :: n_clopper_pearson, cm_n_clopper_pearson

contains

   function clopper_pearson_ci(k,n,alpha,ci_kind) result(ci)
      integer(i64), intent(in) :: k,n
      real(dp), intent(in), optional :: alpha
      character(len=*), intent(in), optional :: ci_kind
      type(confidence_interval) :: ci
      real(dp) :: a
      character(len=12) :: kind

      a = 0.1_dp
      if (present(alpha)) a = alpha
      kind = "upper"
      if (present(ci_kind)) kind = adjustl(ci_kind)
      ci%alpha = a
      ci%kind = kind

      if (k < 0_i64 .or. n < k .or. a < 0.0_dp .or. a > 1.0_dp) then
         ci%status = 1
         return
      end if

      select case(trim(kind))
      case("upper")
         ci%lower = 0.0_dp
         if (k == n) then
            ci%upper = 1.0_dp
         else
            ci%upper = beta_quantile(1.0_dp-a,real(k+1_i64,dp),real(n-k,dp))
         end if
      case("lower")
         ci%upper = 1.0_dp
         if (k == 0_i64) then
            ci%lower = 0.0_dp
         else
            ci%lower = beta_quantile(a,real(k,dp),real(n-k+1_i64,dp))
         end if
      case("two.sided")
         if (k == 0_i64) then
            ci%lower = 0.0_dp
         else
            ci%lower = beta_quantile(a/2.0_dp,real(k,dp),real(n-k+1_i64,dp))
         end if
         if (k == n) then
            ci%upper = 1.0_dp
         else
            ci%upper = beta_quantile(1.0_dp-a/2.0_dp,real(k+1_i64,dp),real(n-k,dp))
         end if
      case default
         ci%status = 2
      end select
   end function clopper_pearson_ci

   subroutine countermeasure_weights(sizev,effect,xi)
      integer(i64), intent(in) :: sizev(:)
      real(dp), intent(in) :: effect(:)
      real(dp), allocatable, intent(out) :: xi(:)
      real(dp), allocatable :: pmf(:)
      integer :: k,j
      call gbinom_pmf_table(sizev,effect,pmf)
      k = ubound(pmf,1)
      allocate(xi(0:k))
      do j = 0, k
         xi(j) = pmf(k-j)
      end do
   end subroutine countermeasure_weights

   real(dp) function cm_upper_equation(p,n,xi,target) result(f)
      real(dp), intent(in) :: p,target
      integer(i64), intent(in) :: n
      real(dp), intent(in) :: xi(0:)
      integer :: j,k
      k = ubound(xi,1)
      f = -target
      do j = 0, k
         f = f+xi(j)*beta_cdf(p,real(j+1,dp),real(n-int(j,i64),dp))
      end do
   end function cm_upper_equation

   real(dp) function cm_lower_equation(p,n,xi,target) result(f)
      real(dp), intent(in) :: p,target
      integer(i64), intent(in) :: n
      real(dp), intent(in) :: xi(0:)
      integer :: j,k
      k = ubound(xi,1)
      f = xi(0)*beta_cdf(p,1.0e-100_dp,real(n+1_i64,dp))-target
      do j = 1, k
         f = f+xi(j)*beta_cdf(p,real(j,dp),real(n-int(j,i64)+1_i64,dp))
      end do
   end function cm_lower_equation

   real(dp) function cm_probability_root(n,xi,target,upper_mode,tol,max_iter,status) result(root)
      integer(i64), intent(in) :: n
      real(dp), intent(in) :: xi(0:),target,tol
      logical, intent(in) :: upper_mode
      integer, intent(in) :: max_iter
      integer, intent(out) :: status
      real(dp) :: lo,hi,mid,flo,fhi,fmid
      integer :: it

      status = 0
      lo = 0.0_dp
      hi = 1.0_dp
      if (upper_mode) then
         flo = cm_upper_equation(lo,n,xi,target)
         fhi = cm_upper_equation(hi,n,xi,target)
      else
         flo = cm_lower_equation(lo,n,xi,target)
         fhi = cm_lower_equation(hi,n,xi,target)
      end if

      if (flo*fhi > 0.0_dp) then
         root = 0.5_dp
         status = 2
         return
      end if

      do it = 1, max_iter
         mid = 0.5_dp*(lo+hi)
         if (upper_mode) then
            fmid = cm_upper_equation(mid,n,xi,target)
         else
            fmid = cm_lower_equation(mid,n,xi,target)
         end if
         if (abs(fmid) <= tol .or. hi-lo <= tol*max(1.0_dp,abs(mid))) exit
         if (flo*fmid <= 0.0_dp) then
            hi = mid
            fhi = fmid
         else
            lo = mid
            flo = fmid
         end if
      end do
      root = 0.5_dp*(lo+hi)
      if (it > max_iter) status = 1
   end function cm_probability_root

   function cm_clopper_pearson_ci(n,sizev,effect,alpha,ci_kind,tol,max_iter) result(ci)
      integer(i64), intent(in) :: n,sizev(:)
      real(dp), intent(in) :: effect(:)
      real(dp), intent(in), optional :: alpha,tol
      character(len=*), intent(in), optional :: ci_kind
      integer, intent(in), optional :: max_iter
      type(confidence_interval) :: ci
      real(dp), allocatable :: xi(:)
      real(dp) :: a,eps
      integer(i64) :: k
      integer :: imax,istat
      character(len=12) :: kind

      a = 0.1_dp
      eps = 1.0e-12_dp
      imax = 200
      kind = "upper"
      if (present(alpha)) a = alpha
      if (present(tol)) eps = tol
      if (present(max_iter)) imax = max_iter
      if (present(ci_kind)) kind = adjustl(ci_kind)
      ci%alpha = a
      ci%kind = kind
      k = sum(sizev)

      if (size(sizev) /= size(effect) .or. any(sizev < 0_i64) .or. n < k .or. &
          any(effect < 0.0_dp) .or. any(effect > 1.0_dp) .or. &
          a < 0.0_dp .or. a > 1.0_dp) then
         ci%status = 1
         return
      end if
      call countermeasure_weights(sizev,effect,xi)

      select case(trim(kind))
      case("upper")
         ci%lower = 0.0_dp
         ci%upper = cm_probability_root(n,xi,1.0_dp-a,.true.,eps,imax,istat)
      case("lower")
         ci%upper = 1.0_dp
         if (xi(0) >= 1.0_dp-a) then
            ci%lower = 0.0_dp
            istat = 0
         else
            ci%lower = cm_probability_root(n,xi,a,.false.,eps,imax,istat)
         end if
      case("two.sided")
         if (xi(0) >= 1.0_dp-a/2.0_dp) then
            ci%lower = 0.0_dp
            istat = 0
         else
            ci%lower = cm_probability_root(n,xi,a/2.0_dp,.false.,eps,imax,istat)
         end if
         ci%upper = cm_probability_root(n,xi,1.0_dp-a/2.0_dp,.true.,eps,imax,istat)
      case default
         ci%status = 2
         return
      end select
      ci%status = istat

   end function cm_clopper_pearson_ci

   integer(i64) function n_clopper_pearson(k,p,alpha,max_n,status) result(nreq)
      integer(i64), intent(in) :: k
      real(dp), intent(in) :: p
      real(dp), intent(in), optional :: alpha
      integer(i64), intent(in), optional :: max_n
      integer, intent(out), optional :: status
      real(dp) :: a,target
      integer(i64) :: lo,hi,mid,limit

      a = 0.1_dp
      if (present(alpha)) a = alpha
      limit = ishft(huge(1_i64),-2)
      if (present(max_n)) limit = max_n
      if (present(status)) status = 0

      if (k < 0_i64 .or. p <= 0.0_dp .or. p >= 1.0_dp .or. &
          a < 0.0_dp .or. a >= 1.0_dp) then
         nreq = -1_i64
         if (present(status)) status = 1
         return
      end if

      target = 1.0_dp-a
      lo = k+1_i64
      if (ordinary_condition(lo)) then
         nreq = lo
         return
      end if
      hi = max(lo+1_i64,2_i64*lo)
      do while (.not. ordinary_condition(hi))
         if (hi >= limit/2_i64) then
            nreq = -1_i64
            if (present(status)) status = 2
            return
         end if
         hi = min(limit,2_i64*hi)
      end do
      do while (lo < hi)
         mid = lo+(hi-lo)/2_i64
         if (ordinary_condition(mid)) then
            hi = mid
         else
            lo = mid+1_i64
         end if
      end do
      nreq = lo

   contains
      logical function ordinary_condition(n) result(ok)
         integer(i64), intent(in) :: n
         ok = beta_cdf(p,real(k+1_i64,dp),real(n-k,dp)) >= target
      end function ordinary_condition
   end function n_clopper_pearson

   integer(i64) function cm_n_clopper_pearson(p,sizev,effect,alpha,max_n,status) result(nreq)
      real(dp), intent(in) :: p
      integer(i64), intent(in) :: sizev(:)
      real(dp), intent(in) :: effect(:)
      real(dp), intent(in), optional :: alpha
      integer(i64), intent(in), optional :: max_n
      integer, intent(out), optional :: status
      real(dp), allocatable :: xi(:)
      real(dp) :: a,target
      integer(i64) :: k,lo,hi,mid,limit

      a = 0.1_dp
      if (present(alpha)) a = alpha
      limit = ishft(huge(1_i64),-2)
      if (present(max_n)) limit = max_n
      if (present(status)) status = 0
      k = sum(sizev)

      if (size(sizev) /= size(effect) .or. any(sizev < 0_i64) .or. &
          any(effect < 0.0_dp) .or. any(effect > 1.0_dp) .or. &
          p < 0.0_dp .or. p > 1.0_dp .or. a < 0.0_dp .or. a >= 1.0_dp) then
         nreq = -1_i64
         if (present(status)) status = 1
         return
      end if

      call countermeasure_weights(sizev,effect,xi)
      target = 1.0_dp-a
      lo = k+1_i64
      if (cm_condition(lo)) then
         nreq = lo
         return
      end if
      hi = max(lo+1_i64,2_i64*lo)
      do while (.not. cm_condition(hi))
         if (hi >= limit/2_i64) then
            nreq = -1_i64
            if (present(status)) status = 2
            return
         end if
         hi = min(limit,2_i64*hi)
      end do
      do while (lo < hi)
         mid = lo+(hi-lo)/2_i64
         if (cm_condition(mid)) then
            hi = mid
         else
            lo = mid+1_i64
         end if
      end do
      nreq = lo

   contains
      logical function cm_condition(n) result(ok)
         integer(i64), intent(in) :: n
         integer :: j
         real(dp) :: s
         s = 0.0_dp
         do j = 0, ubound(xi,1)
            s = s+xi(j)*beta_cdf(p,real(j+1,dp),real(n-int(j,i64),dp))
         end do
         ok = s >= target
      end function cm_condition
   end function cm_n_clopper_pearson

end module genbinom_clopper
