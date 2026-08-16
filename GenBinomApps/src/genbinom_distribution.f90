module genbinom_distribution
   use genbinom_kinds, only : dp, i64
   implicit none
   private

   public :: gbinom_pmf_table, gbinom_cdf_table
   public :: dgbinom, dgbinom_vec, pgbinom, pgbinom_vec
   public :: qgbinom, qgbinom_vec, rgbinom
   public :: generalized_binomial_mean, generalized_binomial_variance
   public :: set_genbinom_seed

contains

   pure elemental real(dp) function safe_probability_product(a,b) result(v)
      real(dp), intent(in) :: a,b
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         v = 0.0_dp
      else if (a < tiny(1.0_dp)/b) then
         v = 0.0_dp
      else
         v = a*b
      end if
   end function safe_probability_product

   subroutine validate_parameters(sizev,prob)
      integer(i64), intent(in) :: sizev(:)
      real(dp), intent(in) :: prob(:)
      if (size(sizev) == 0 .or. size(prob) /= size(sizev)) &
         error stop "GenBinomApps: size and prob must be nonempty and have equal lengths"
      if (any(sizev < 0_i64)) error stop "GenBinomApps: size values must be nonnegative"
      if (sum(sizev) < 1_i64) error stop "GenBinomApps: sum(size) must be positive"
      if (any(prob < 0.0_dp) .or. any(prob > 1.0_dp)) &
         error stop "GenBinomApps: probabilities must lie in [0,1]"
      if (sum(sizev) > int(huge(1)-1,i64)) &
         error stop "GenBinomApps: total number of trials is too large for an array"
   end subroutine validate_parameters

   pure real(dp) function generalized_binomial_mean(sizev,prob) result(mu)
      integer(i64), intent(in) :: sizev(:)
      real(dp), intent(in) :: prob(:)
      integer :: i
      mu = 0.0_dp
      do i = 1, size(sizev)
         mu = mu+real(sizev(i),dp)*prob(i)
      end do
   end function generalized_binomial_mean

   pure real(dp) function generalized_binomial_variance(sizev,prob) result(v)
      integer(i64), intent(in) :: sizev(:)
      real(dp), intent(in) :: prob(:)
      integer :: i
      v = 0.0_dp
      do i = 1, size(sizev)
         v = v+real(sizev(i),dp)*prob(i)*(1.0_dp-prob(i))
      end do
   end function generalized_binomial_variance

   subroutine gbinom_pmf_table(sizev,prob,pmf)
      integer(i64), intent(in) :: sizev(:)
      real(dp), intent(in) :: prob(:)
      real(dp), allocatable, intent(out) :: pmf(:)
      integer :: n,current,g,j,k
      real(dp) :: p,q

      call validate_parameters(sizev,prob)
      n = int(sum(sizev))
      allocate(pmf(0:n))
      pmf = 0.0_dp
      pmf(0) = 1.0_dp
      current = 0

      do g = 1, size(sizev)
         p = prob(g)
         q = 1.0_dp-p
         do j = 1, int(sizev(g))
            current = current+1
            pmf(current) = safe_probability_product(pmf(current-1),p)
            do k = current-1, 1, -1
               pmf(k) = safe_probability_product(pmf(k),q) + &
                        safe_probability_product(pmf(k-1),p)
            end do
            pmf(0) = safe_probability_product(pmf(0),q)
         end do
      end do
   end subroutine gbinom_pmf_table

   subroutine gbinom_cdf_table(sizev,prob,cdf)
      integer(i64), intent(in) :: sizev(:)
      real(dp), intent(in) :: prob(:)
      real(dp), allocatable, intent(out) :: cdf(:)
      real(dp), allocatable :: pmf(:)
      integer :: k,n
      call gbinom_pmf_table(sizev,prob,pmf)
      n = ubound(pmf,1)
      allocate(cdf(0:n))
      cdf(0) = pmf(0)
      do k = 1, n
         cdf(k) = min(1.0_dp,cdf(k-1)+pmf(k))
      end do
      cdf(n) = 1.0_dp
   end subroutine gbinom_cdf_table

   real(dp) function dgbinom(x,sizev,prob,log_p) result(value)
      integer(i64), intent(in) :: x
      integer(i64), intent(in) :: sizev(:)
      real(dp), intent(in) :: prob(:)
      logical, intent(in), optional :: log_p
      real(dp), allocatable :: pmf(:)
      logical :: lp
      integer(i64) :: n

      lp = .false.
      if (present(log_p)) lp = log_p
      call validate_parameters(sizev,prob)
      n = sum(sizev)
      if (x < 0_i64 .or. x > n) then
         if (lp) then
            value = -huge(1.0_dp)
         else
            value = 0.0_dp
         end if
         return
      end if
      call gbinom_pmf_table(sizev,prob,pmf)
      if (lp) then
         if (pmf(int(x)) > 0.0_dp) then
            value = log(pmf(int(x)))
         else
            value = -huge(1.0_dp)
         end if
      else
         value = pmf(int(x))
      end if
   end function dgbinom

   subroutine dgbinom_vec(x,sizev,prob,value,log_p)
      integer(i64), intent(in) :: x(:),sizev(:)
      real(dp), intent(in) :: prob(:)
      real(dp), intent(out) :: value(size(x))
      logical, intent(in), optional :: log_p
      real(dp), allocatable :: pmf(:)
      logical :: lp
      integer :: i,n,k

      lp = .false.
      if (present(log_p)) lp = log_p
      call gbinom_pmf_table(sizev,prob,pmf)
      n = ubound(pmf,1)
      do i = 1, size(x)
         if (x(i) < 0_i64 .or. x(i) > int(n,i64)) then
            value(i) = merge(-huge(1.0_dp),0.0_dp,lp)
         else
            k = int(x(i))
            if (lp) then
               if (pmf(k) > 0.0_dp) then
                  value(i) = log(pmf(k))
               else
                  value(i) = -huge(1.0_dp)
               end if
            else
               value(i) = pmf(k)
            end if
         end if
      end do
   end subroutine dgbinom_vec

   real(dp) function pgbinom(q,sizev,prob,lower_tail,log_p) result(value)
      real(dp), intent(in) :: q
      integer(i64), intent(in) :: sizev(:)
      real(dp), intent(in) :: prob(:)
      logical, intent(in), optional :: lower_tail,log_p
      real(dp), allocatable :: cdf(:)
      logical :: lt,lp
      integer :: k,n
      real(dp) :: p

      lt = .true.
      lp = .false.
      if (present(lower_tail)) lt = lower_tail
      if (present(log_p)) lp = log_p
      call gbinom_cdf_table(sizev,prob,cdf)
      n = ubound(cdf,1)

      if (q < 0.0_dp) then
         p = 0.0_dp
      else if (q >= real(n,dp)) then
         p = 1.0_dp
      else
         k = floor(q)
         p = cdf(k)
      end if
      if (.not. lt) p = max(0.0_dp,1.0_dp-p)
      if (lp) then
         if (p > 0.0_dp) then
            value = log(p)
         else
            value = -huge(1.0_dp)
         end if
      else
         value = p
      end if
   end function pgbinom

   subroutine pgbinom_vec(q,sizev,prob,value,lower_tail,log_p)
      real(dp), intent(in) :: q(:)
      integer(i64), intent(in) :: sizev(:)
      real(dp), intent(in) :: prob(:)
      real(dp), intent(out) :: value(size(q))
      logical, intent(in), optional :: lower_tail,log_p
      real(dp), allocatable :: cdf(:)
      logical :: lt,lp
      integer :: i,k,n
      real(dp) :: p

      lt = .true.
      lp = .false.
      if (present(lower_tail)) lt = lower_tail
      if (present(log_p)) lp = log_p
      call gbinom_cdf_table(sizev,prob,cdf)
      n = ubound(cdf,1)

      do i = 1, size(q)
         if (q(i) < 0.0_dp) then
            p = 0.0_dp
         else if (q(i) >= real(n,dp)) then
            p = 1.0_dp
         else
            k = floor(q(i))
            p = cdf(k)
         end if
         if (.not. lt) p = max(0.0_dp,1.0_dp-p)
         if (lp) then
            if (p > 0.0_dp) then
               value(i) = log(p)
            else
               value(i) = -huge(1.0_dp)
            end if
         else
            value(i) = p
         end if
      end do
   end subroutine pgbinom_vec

   integer(i64) function qgbinom(p,sizev,prob,lower_tail,log_p) result(q)
      real(dp), intent(in) :: p
      integer(i64), intent(in) :: sizev(:)
      real(dp), intent(in) :: prob(:)
      logical, intent(in), optional :: lower_tail,log_p
      real(dp), allocatable :: cdf(:)
      logical :: lt,lp
      real(dp) :: target
      integer :: lo,hi,mid,n

      lt = .true.
      lp = .false.
      if (present(lower_tail)) lt = lower_tail
      if (present(log_p)) lp = log_p
      if (lp) then
         if (p > 0.0_dp) then
            q = -huge(1_i64)
            return
         end if
         target = exp(p)
      else
         if (p < 0.0_dp .or. p > 1.0_dp) then
            q = -huge(1_i64)
            return
         end if
         target = p
      end if
      if (.not. lt) target = 1.0_dp-target

      call gbinom_cdf_table(sizev,prob,cdf)
      n = ubound(cdf,1)
      if (target <= 0.0_dp) then
         q = 0_i64
         return
      end if
      if (target >= 1.0_dp) then
         q = int(n,i64)
         return
      end if

      lo = 0
      hi = n
      do while (lo < hi)
         mid = lo+(hi-lo)/2
         if (cdf(mid) >= target) then
            hi = mid
         else
            lo = mid+1
         end if
      end do
      q = int(lo,i64)
   end function qgbinom

   subroutine qgbinom_vec(p,sizev,prob,value,lower_tail,log_p)
      real(dp), intent(in) :: p(:)
      integer(i64), intent(in) :: sizev(:)
      real(dp), intent(in) :: prob(:)
      integer(i64), intent(out) :: value(size(p))
      logical, intent(in), optional :: lower_tail,log_p
      real(dp), allocatable :: cdf(:)
      logical :: lt,lp
      real(dp) :: target
      integer :: i,lo,hi,mid,n

      lt = .true.
      lp = .false.
      if (present(lower_tail)) lt = lower_tail
      if (present(log_p)) lp = log_p
      call gbinom_cdf_table(sizev,prob,cdf)
      n = ubound(cdf,1)

      do i = 1, size(p)
         if (lp) then
            if (p(i) > 0.0_dp) then
               value(i) = -huge(1_i64)
               cycle
            end if
            target = exp(p(i))
         else
            if (p(i) < 0.0_dp .or. p(i) > 1.0_dp) then
               value(i) = -huge(1_i64)
               cycle
            end if
            target = p(i)
         end if
         if (.not. lt) target = 1.0_dp-target
         if (target <= 0.0_dp) then
            value(i) = 0_i64
            cycle
         else if (target >= 1.0_dp) then
            value(i) = int(n,i64)
            cycle
         end if
         lo = 0
         hi = n
         do while (lo < hi)
            mid = lo+(hi-lo)/2
            if (cdf(mid) >= target) then
               hi = mid
            else
               lo = mid+1
            end if
         end do
         value(i) = int(lo,i64)
      end do
   end subroutine qgbinom_vec

   subroutine rgbinom(x,sizev,prob)
      integer(i64), intent(out) :: x(:)
      integer(i64), intent(in) :: sizev(:)
      real(dp), intent(in) :: prob(:)
      real(dp), allocatable :: cdf(:)
      real(dp) :: u
      integer :: i,lo,hi,mid,n

      call gbinom_cdf_table(sizev,prob,cdf)
      n = ubound(cdf,1)
      do i = 1, size(x)
         call random_number(u)
         lo = 0
         hi = n
         do while (lo < hi)
            mid = lo+(hi-lo)/2
            if (cdf(mid) >= u) then
               hi = mid
            else
               lo = mid+1
            end if
         end do
         x(i) = int(lo,i64)
      end do
   end subroutine rgbinom

   subroutine set_genbinom_seed(seed)
      integer, intent(in) :: seed
      integer :: n,i
      integer, allocatable :: state(:)
      call random_seed(size=n)
      allocate(state(n))
      do i = 1, n
         state(i) = mod(seed+104729*i,2147483646)+1
      end do
      call random_seed(put=state)
   end subroutine set_genbinom_seed

end module genbinom_distribution
