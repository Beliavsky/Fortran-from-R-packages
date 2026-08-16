module polya_aeppli_distribution
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use polya_aeppli_kinds, only : dp
   use polya_aeppli_numerics, only : log1p_stable, logaddexp, log1mexp, &
      rand_poisson, rand_geometric_failures
   implicit none
   private

   public :: d_polya_aeppli, p_polya_aeppli, q_polya_aeppli, r_polya_aeppli
   public :: d_polya_aeppli_vec, p_polya_aeppli_vec, q_polya_aeppli_vec
   public :: r_polya_aeppli_vec
   public :: log_pmf_array, log_cdf_array, log_sf_array, log_tail_pa
   public :: polya_aeppli_mean, polya_aeppli_variance

contains

   pure elemental logical function valid_parameters(lambda,prob) result(ok)
      real(dp), intent(in) :: lambda,prob
      ok = lambda > 0.0_dp .and. prob >= 0.0_dp .and. prob < 1.0_dp .and. &
           ieee_is_finite(lambda) .and. ieee_is_finite(prob)
   end function valid_parameters

   pure elemental real(dp) function polya_aeppli_mean(lambda,prob) result(mu)
      real(dp), intent(in) :: lambda,prob
      mu = lambda/(1.0_dp-prob)
   end function polya_aeppli_mean

   pure elemental real(dp) function polya_aeppli_variance(lambda,prob) result(v)
      real(dp), intent(in) :: lambda,prob
      real(dp) :: mu
      mu = polya_aeppli_mean(lambda,prob)
      v = mu*(1.0_dp+prob)/(1.0_dp-prob)
   end function polya_aeppli_variance

   subroutine log_pmf_array(xmax,lambda,prob,lpmf,status)
      integer, intent(in) :: xmax
      real(dp), intent(in) :: lambda,prob
      real(dp), allocatable, intent(out) :: lpmf(:)
      integer, intent(out), optional :: status
      real(dp) :: qprob,next_ratio,term
      integer :: x

      if (present(status)) status = 0
      if (.not. valid_parameters(lambda,prob) .or. xmax < 0) then
         allocate(lpmf(0))
         if (present(status)) status = 1
         return
      end if

      allocate(lpmf(0:max(1,xmax)))
      qprob = 1.0_dp-prob
      lpmf(0) = -lambda
      lpmf(1) = -lambda+log(lambda*qprob)

      do x = 2, xmax
         term = prob*prob*real(x-2,dp)*exp(lpmf(x-2)-lpmf(x-1))
         next_ratio = (lambda*qprob+2.0_dp*prob*real(x-1,dp)-term)/real(x,dp)
         if (next_ratio <= 0.0_dp .or. .not. ieee_is_finite(next_ratio)) then
            if (present(status)) status = 2
            lpmf(x:) = -huge(1.0_dp)
            return
         end if
         lpmf(x) = lpmf(x-1)+log(next_ratio)
      end do
   end subroutine log_pmf_array

   subroutine log_cdf_array(lpmf,lcdf)
      real(dp), intent(in) :: lpmf(0:)
      real(dp), allocatable, intent(out) :: lcdf(:)
      integer :: x,n
      n = ubound(lpmf,1)
      allocate(lcdf(0:n))
      lcdf(0) = lpmf(0)
      do x = 1, n
         lcdf(x) = logaddexp(lcdf(x-1),lpmf(x))
      end do
   end subroutine log_cdf_array

   real(dp) function log_tail_pa(x,lambda,prob,max_iter,status) result(logtail)
      integer, intent(in) :: x
      real(dp), intent(in) :: lambda,prob
      integer, intent(in), optional :: max_iter
      integer, intent(out), optional :: status
      real(dp), allocatable :: lp(:)
      real(dp) :: lminus2,lminus1,lbase,lcur,next_term,sumexp,qprob,ratio
      integer :: i,imax,istat

      if (present(status)) status = 0
      if (x < 0) then
         logtail = 0.0_dp
         return
      end if
      imax = 10000
      if (present(max_iter)) imax = max_iter

      call log_pmf_array(x+1,lambda,prob,lp,istat)
      if (istat /= 0) then
         logtail = -huge(1.0_dp)
         if (present(status)) status = istat
         return
      end if

      lminus2 = lp(x)
      lminus1 = lp(x+1)
      lbase = lminus1
      sumexp = 1.0_dp
      qprob = 1.0_dp-prob

      do i = x+2, x+1+imax
         ratio = (lambda*qprob+2.0_dp*prob*real(i-1,dp) - &
                  prob*prob*real(i-2,dp)*exp(lminus2-lminus1))/real(i,dp)
         if (ratio <= 0.0_dp) exit
         lcur = lminus1+log(ratio)
         next_term = exp(lcur-lbase)
         sumexp = sumexp+next_term
         lminus2 = lminus1
         lminus1 = lcur
         if (next_term <= 2.0_dp*epsilon(1.0_dp)) exit
      end do

      if (i > x+1+imax .and. present(status)) status = 3
      logtail = lbase+log(sumexp)
   end function log_tail_pa

   subroutine log_sf_array(xmax,lambda,prob,lpmf,lsf,status)
      integer, intent(in) :: xmax
      real(dp), intent(in) :: lambda,prob
      real(dp), intent(in) :: lpmf(0:)
      real(dp), allocatable, intent(out) :: lsf(:)
      integer, intent(out), optional :: status
      real(dp) :: top
      integer :: x,istat

      if (present(status)) status = 0
      if (xmax < 0) then
         allocate(lsf(0))
         return
      end if
      top = log_tail_pa(xmax,lambda,prob,status=istat)
      if (istat /= 0 .and. present(status)) status = istat
      allocate(lsf(0:xmax))
      lsf(xmax) = top
      do x = xmax-1, 0, -1
         lsf(x) = logaddexp(lsf(x+1),lpmf(x+1))
      end do
   end subroutine log_sf_array

   real(dp) function d_polya_aeppli(x,lambda,prob,log_p) result(value)
      real(dp), intent(in) :: x,lambda,prob
      logical, intent(in), optional :: log_p
      logical :: lp
      integer :: k,istat
      real(dp), allocatable :: lpmf(:)

      lp = .false.
      if (present(log_p)) lp = log_p
      if (.not. valid_parameters(lambda,prob)) then
         value = merge(-huge(1.0_dp),0.0_dp,lp)
         return
      end if
      if (.not. ieee_is_finite(x) .or. x < 0.0_dp .or. &
          abs(x-anint(x)) >= sqrt(epsilon(1.0_dp))) then
         value = merge(-huge(1.0_dp),0.0_dp,lp)
         return
      end if

      k = nint(x)
      call log_pmf_array(k,lambda,prob,lpmf,istat)
      if (istat /= 0) then
         value = merge(-huge(1.0_dp),0.0_dp,lp)
      else if (lp) then
         value = lpmf(k)
      else
         value = exp(lpmf(k))
      end if
   end function d_polya_aeppli

   real(dp) function p_polya_aeppli(q,lambda,prob,lower_tail,log_p) result(value)
      real(dp), intent(in) :: q,lambda,prob
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      integer :: k,istat
      real(dp), allocatable :: lpmf(:),lcdf(:),lsf(:)
      real(dp) :: lprob,mu

      lt = .true.
      lp = .false.
      if (present(lower_tail)) lt = lower_tail
      if (present(log_p)) lp = log_p

      if (.not. valid_parameters(lambda,prob)) then
         value = merge(-huge(1.0_dp),0.0_dp,lp)
         return
      end if

      if (.not. ieee_is_finite(q)) then
         if (q < 0.0_dp) then
            if (lt) then
               lprob = -huge(1.0_dp)
            else
               lprob = 0.0_dp
            end if
         else
            if (lt) then
               lprob = 0.0_dp
            else
               lprob = -huge(1.0_dp)
            end if
         end if
      else if (q < 0.0_dp) then
         if (lt) then
            lprob = -huge(1.0_dp)
         else
            lprob = 0.0_dp
         end if
      else
         k = floor(q)
         mu = polya_aeppli_mean(lambda,prob)
         if (lt) then
            call log_pmf_array(k,lambda,prob,lpmf,istat)
            if (istat /= 0) then
               lprob = -huge(1.0_dp)
            else
               call log_cdf_array(lpmf,lcdf)
               lprob = lcdf(k)
            end if
         else if (real(k,dp) > mu) then
            call log_pmf_array(k,lambda,prob,lpmf,istat)
            if (istat /= 0) then
               lprob = -huge(1.0_dp)
            else
               call log_sf_array(k,lambda,prob,lpmf,lsf,istat)
               lprob = lsf(k)
            end if
         else
            call log_pmf_array(k,lambda,prob,lpmf,istat)
            if (istat /= 0) then
               lprob = -huge(1.0_dp)
            else
               call log_cdf_array(lpmf,lcdf)
               if (lcdf(k) >= 0.0_dp) then
                  lprob = -huge(1.0_dp)
               else
                  lprob = log1mexp(lcdf(k))
               end if
            end if
         end if
      end if

      if (lp) then
         value = lprob
      else if (lprob <= -0.5_dp*huge(1.0_dp)) then
         value = 0.0_dp
      else
         value = exp(lprob)
      end if
   end function p_polya_aeppli

   integer function q_polya_aeppli(p,lambda,prob,lower_tail,log_p) result(q)
      real(dp), intent(in) :: p,lambda,prob
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      real(dp) :: logtarget,mu,var,lcur
      integer :: lo,hi,mid

      lt = .true.
      lp = .false.
      if (present(lower_tail)) lt = lower_tail
      if (present(log_p)) lp = log_p

      if (.not. valid_parameters(lambda,prob)) then
         q = -huge(1)
         return
      end if

      if (lp) then
         if (p > 0.0_dp) then
            q = -huge(1)
            return
         end if
         logtarget = p
      else
         if (p < 0.0_dp .or. p > 1.0_dp .or. .not. ieee_is_finite(p)) then
            q = -huge(1)
            return
         end if
         if (p <= 0.0_dp) then
            if (lt) then
               q = 0
            else
               q = huge(1)
            end if
            return
         end if
         if (p >= 1.0_dp) then
            if (lt) then
               q = huge(1)
            else
               q = 0
            end if
            return
         end if
         logtarget = log(p)
      end if

      if (lp) then
         if (logtarget <= -0.5_dp*huge(1.0_dp)) then
            if (lt) then
               q = 0
            else
               q = huge(1)
            end if
            return
         end if
         if (logtarget >= 0.0_dp) then
            if (lt) then
               q = huge(1)
            else
               q = 0
            end if
            return
         end if
      end if

      mu = polya_aeppli_mean(lambda,prob)
      var = polya_aeppli_variance(lambda,prob)
      hi = max(1,int(ceiling(mu+8.0_dp*sqrt(var)+10.0_dp)))
      lo = 0

      if (lt) then
         do
            lcur = p_polya_aeppli(real(hi,dp),lambda,prob,.true.,.true.)
            if (lcur >= logtarget) exit
            if (hi > ishft(huge(1),-1)) then
               q = huge(1)
               return
            end if
            hi = 2*hi
         end do
         do while (lo < hi)
            mid = lo+(hi-lo)/2
            lcur = p_polya_aeppli(real(mid,dp),lambda,prob,.true.,.true.)
            if (lcur >= logtarget) then
               hi = mid
            else
               lo = mid+1
            end if
         end do
      else
         do
            lcur = p_polya_aeppli(real(hi,dp),lambda,prob,.false.,.true.)
            if (lcur <= logtarget) exit
            if (hi > ishft(huge(1),-1)) then
               q = huge(1)
               return
            end if
            hi = 2*hi
         end do
         do while (lo < hi)
            mid = lo+(hi-lo)/2
            lcur = p_polya_aeppli(real(mid,dp),lambda,prob,.false.,.true.)
            if (lcur <= logtarget) then
               hi = mid
            else
               lo = mid+1
            end if
         end do
      end if
      q = lo
   end function q_polya_aeppli

   subroutine r_polya_aeppli(x,lambda,prob)
      integer, intent(out) :: x(:)
      real(dp), intent(in) :: lambda,prob
      integer :: i,j,nclusters,total

      if (.not. valid_parameters(lambda,prob)) then
         error stop "r_polya_aeppli: invalid lambda or prob"
      end if

      do i = 1, size(x)
         nclusters = rand_poisson(lambda)
         total = nclusters
         do j = 1, nclusters
            total = total+rand_geometric_failures(1.0_dp-prob)
         end do
         x(i) = total
      end do
   end subroutine r_polya_aeppli

   integer function recycle_length(n1,n2,n3) result(n)
      integer, intent(in) :: n1,n2,n3
      n = max(n1,max(n2,n3))
   end function recycle_length

   subroutine d_polya_aeppli_vec(x,lambda,prob,value,log_p)
      real(dp), intent(in) :: x(:),lambda(:),prob(:)
      real(dp), intent(out) :: value(:)
      logical, intent(in), optional :: log_p
      integer :: i,n
      n = recycle_length(size(x),size(lambda),size(prob))
      if (size(value) /= n) error stop "d_polya_aeppli_vec: wrong output length"
      do i = 1, n
         value(i) = d_polya_aeppli(x(mod(i-1,size(x))+1), &
                    lambda(mod(i-1,size(lambda))+1),prob(mod(i-1,size(prob))+1),log_p)
      end do
   end subroutine d_polya_aeppli_vec

   subroutine p_polya_aeppli_vec(q,lambda,prob,value,lower_tail,log_p)
      real(dp), intent(in) :: q(:),lambda(:),prob(:)
      real(dp), intent(out) :: value(:)
      logical, intent(in), optional :: lower_tail,log_p
      integer :: i,n
      n = recycle_length(size(q),size(lambda),size(prob))
      if (size(value) /= n) error stop "p_polya_aeppli_vec: wrong output length"
      do i = 1, n
         value(i) = p_polya_aeppli(q(mod(i-1,size(q))+1), &
                    lambda(mod(i-1,size(lambda))+1),prob(mod(i-1,size(prob))+1), &
                    lower_tail,log_p)
      end do
   end subroutine p_polya_aeppli_vec

   subroutine q_polya_aeppli_vec(p,lambda,prob,value,lower_tail,log_p)
      real(dp), intent(in) :: p(:),lambda(:),prob(:)
      integer, intent(out) :: value(:)
      logical, intent(in), optional :: lower_tail,log_p
      integer :: i,n
      n = recycle_length(size(p),size(lambda),size(prob))
      if (size(value) /= n) error stop "q_polya_aeppli_vec: wrong output length"
      do i = 1, n
         value(i) = q_polya_aeppli(p(mod(i-1,size(p))+1), &
                    lambda(mod(i-1,size(lambda))+1),prob(mod(i-1,size(prob))+1), &
                    lower_tail,log_p)
      end do
   end subroutine q_polya_aeppli_vec

   subroutine r_polya_aeppli_vec(x,lambda,prob)
      integer, intent(out) :: x(:)
      real(dp), intent(in) :: lambda(:),prob(:)
      integer :: i,j,nclusters,total
      real(dp) :: li,pi
      if (size(lambda) == 0 .or. size(prob) == 0) then
         error stop "r_polya_aeppli_vec: empty parameter vector"
      end if
      do i = 1, size(x)
         li = lambda(mod(i-1,size(lambda))+1)
         pi = prob(mod(i-1,size(prob))+1)
         if (.not. valid_parameters(li,pi)) error stop "r_polya_aeppli_vec: invalid parameter"
         nclusters = rand_poisson(li)
         total = nclusters
         do j = 1, nclusters
            total = total+rand_geometric_failures(1.0_dp-pi)
         end do
         x(i) = total
      end do
   end subroutine r_polya_aeppli_vec

end module polya_aeppli_distribution
