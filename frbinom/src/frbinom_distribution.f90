module frbinom_distribution
   use frbinom_kinds, only : dp
   use frbinom_core, only : frbinom_pmf_table, frbinom_cdf_table, &
      frbinom2_pmf_table, frbinom2_cdf_table
   implicit none
   private

   public :: dfrbinom, pfrbinom, qfrbinom, rfrbinom
   public :: dfrbinom_vec, pfrbinom_vec, qfrbinom_vec
   public :: dfrbinom2, pfrbinom2, qfrbinom2, rfrbinom2
   public :: dfrbinom2_vec, pfrbinom2_vec, qfrbinom2_vec
   public :: set_frbinom_seed

contains

   real(dp) function dfrbinom(x,ntrials,prob,h,c,start) result(value)
      real(dp), intent(in) :: x,prob,h,c
      integer, intent(in) :: ntrials
      logical, intent(in), optional :: start
      real(dp), allocatable :: pmf(:)
      integer :: k,istat

      call frbinom_pmf_table(ntrials,prob,h,c,start,pmf,istat)
      if (istat /= 0 .or. x < 0.0_dp .or. x > real(ntrials,dp) .or. &
          abs(x-anint(x)) > 10.0_dp*epsilon(1.0_dp)) then
         value = 0.0_dp
      else
         k = nint(x)
         value = pmf(k)
      end if
   end function dfrbinom

   subroutine dfrbinom_vec(x,ntrials,prob,h,c,value,start)
      real(dp), intent(in) :: x(:),prob,h,c
      integer, intent(in) :: ntrials
      real(dp), intent(out) :: value(size(x))
      logical, intent(in), optional :: start
      real(dp), allocatable :: pmf(:)
      integer :: i,k,istat

      call frbinom_pmf_table(ntrials,prob,h,c,start,pmf,istat)
      if (istat /= 0) then
         value = 0.0_dp
         return
      end if
      do i = 1, size(x)
         if (x(i) < 0.0_dp .or. x(i) > real(ntrials,dp) .or. &
             abs(x(i)-anint(x(i))) > 10.0_dp*epsilon(1.0_dp)) then
            value(i) = 0.0_dp
         else
            k = nint(x(i))
            value(i) = pmf(k)
         end if
      end do
   end subroutine dfrbinom_vec

   real(dp) function pfrbinom(x,ntrials,prob,h,c,start,lower_tail) result(value)
      real(dp), intent(in) :: x,prob,h,c
      integer, intent(in) :: ntrials
      logical, intent(in), optional :: start,lower_tail
      real(dp), allocatable :: cdf(:)
      logical :: lt
      integer :: k,istat

      lt = .true.
      if (present(lower_tail)) lt = lower_tail
      call frbinom_cdf_table(ntrials,prob,h,c,start,cdf,istat)
      if (istat /= 0) then
         value = 0.0_dp
         return
      end if
      if (x < 0.0_dp) then
         value = 0.0_dp
      else if (x >= real(ntrials,dp)) then
         value = 1.0_dp
      else
         k = floor(x)
         value = cdf(k)
      end if
      if (.not. lt) value = max(0.0_dp,1.0_dp-value)
   end function pfrbinom

   subroutine pfrbinom_vec(x,ntrials,prob,h,c,value,start,lower_tail)
      real(dp), intent(in) :: x(:),prob,h,c
      integer, intent(in) :: ntrials
      real(dp), intent(out) :: value(size(x))
      logical, intent(in), optional :: start,lower_tail
      real(dp), allocatable :: cdf(:)
      logical :: lt
      integer :: i,k,istat

      lt = .true.
      if (present(lower_tail)) lt = lower_tail
      call frbinom_cdf_table(ntrials,prob,h,c,start,cdf,istat)
      if (istat /= 0) then
         value = 0.0_dp
         return
      end if
      do i = 1, size(x)
         if (x(i) < 0.0_dp) then
            value(i) = 0.0_dp
         else if (x(i) >= real(ntrials,dp)) then
            value(i) = 1.0_dp
         else
            k = floor(x(i))
            value(i) = cdf(k)
         end if
      end do
      if (.not. lt) value = max(0.0_dp,1.0_dp-value)
   end subroutine pfrbinom_vec

   integer function qfrbinom(p,ntrials,prob,h,c,start,lower_tail) result(q)
      real(dp), intent(in) :: p,prob,h,c
      integer, intent(in) :: ntrials
      logical, intent(in), optional :: start,lower_tail
      real(dp), allocatable :: cdf(:)
      logical :: lt
      real(dp) :: target
      integer :: lo,hi,mid,istat

      lt = .true.
      if (present(lower_tail)) lt = lower_tail
      if (p < 0.0_dp .or. p > 1.0_dp) then
         q = -huge(1)
         return
      end if
      target = merge(p,1.0_dp-p,lt)
      call frbinom_cdf_table(ntrials,prob,h,c,start,cdf,istat)
      if (istat /= 0) then
         q = -huge(1)
         return
      end if
      if (target <= 0.0_dp) then
         q = 0
         return
      end if
      if (target >= 1.0_dp) then
         q = ntrials
         return
      end if
      lo = 0
      hi = ntrials
      do while (lo < hi)
         mid = lo+(hi-lo)/2
         if (cdf(mid) >= target) then
            hi = mid
         else
            lo = mid+1
         end if
      end do
      q = lo
   end function qfrbinom

   subroutine qfrbinom_vec(p,ntrials,prob,h,c,q,start,lower_tail)
      real(dp), intent(in) :: p(:),prob,h,c
      integer, intent(in) :: ntrials
      integer, intent(out) :: q(size(p))
      logical, intent(in), optional :: start,lower_tail
      real(dp), allocatable :: cdf(:)
      logical :: lt
      real(dp) :: target
      integer :: i,lo,hi,mid,istat

      lt = .true.
      if (present(lower_tail)) lt = lower_tail
      call frbinom_cdf_table(ntrials,prob,h,c,start,cdf,istat)
      if (istat /= 0) then
         q = -huge(1)
         return
      end if
      do i = 1, size(p)
         if (p(i) < 0.0_dp .or. p(i) > 1.0_dp) then
            q(i) = -huge(1)
            cycle
         end if
         target = merge(p(i),1.0_dp-p(i),lt)
         if (target <= 0.0_dp) then
            q(i) = 0
         else if (target >= 1.0_dp) then
            q(i) = ntrials
         else
            lo = 0
            hi = ntrials
            do while (lo < hi)
               mid = lo+(hi-lo)/2
               if (cdf(mid) >= target) then
                  hi = mid
               else
                  lo = mid+1
               end if
            end do
            q(i) = lo
         end if
      end do
   end subroutine qfrbinom_vec

   subroutine rfrbinom(x,ntrials,prob,h,c,start)
      integer, intent(out) :: x(:)
      integer, intent(in) :: ntrials
      real(dp), intent(in) :: prob,h,c
      logical, intent(in), optional :: start
      real(dp), allocatable :: cdf(:)
      real(dp) :: u
      integer :: i,lo,hi,mid,istat

      call frbinom_cdf_table(ntrials,prob,h,c,start,cdf,istat)
      if (istat /= 0) error stop "rfrbinom: invalid parameters"
      do i = 1, size(x)
         call random_number(u)
         lo = 0
         hi = ntrials
         do while (lo < hi)
            mid = lo+(hi-lo)/2
            if (cdf(mid) >= u) then
               hi = mid
            else
               lo = mid+1
            end if
         end do
         x(i) = lo
      end do
   end subroutine rfrbinom

   real(dp) function dfrbinom2(x,ntrials,h,c,la,start) result(value)
      real(dp), intent(in) :: x,h,c,la
      integer, intent(in) :: ntrials
      logical, intent(in), optional :: start
      real(dp), allocatable :: pmf(:)
      integer :: k,istat

      call frbinom2_pmf_table(ntrials,h,c,la,start,pmf,istat)
      if (istat /= 0 .or. x < 0.0_dp .or. x > real(ntrials,dp) .or. &
          abs(x-anint(x)) > 10.0_dp*epsilon(1.0_dp)) then
         value = 0.0_dp
      else
         k = nint(x)
         value = pmf(k)
      end if
   end function dfrbinom2

   subroutine dfrbinom2_vec(x,ntrials,h,c,la,value,start)
      real(dp), intent(in) :: x(:),h,c,la
      integer, intent(in) :: ntrials
      real(dp), intent(out) :: value(size(x))
      logical, intent(in), optional :: start
      real(dp), allocatable :: pmf(:)
      integer :: i,k,istat

      call frbinom2_pmf_table(ntrials,h,c,la,start,pmf,istat)
      if (istat /= 0) then
         value = 0.0_dp
         return
      end if
      do i = 1, size(x)
         if (x(i) < 0.0_dp .or. x(i) > real(ntrials,dp) .or. &
             abs(x(i)-anint(x(i))) > 10.0_dp*epsilon(1.0_dp)) then
            value(i) = 0.0_dp
         else
            k = nint(x(i))
            value(i) = pmf(k)
         end if
      end do
   end subroutine dfrbinom2_vec

   real(dp) function pfrbinom2(x,ntrials,h,c,la,start,lower_tail) result(value)
      real(dp), intent(in) :: x,h,c,la
      integer, intent(in) :: ntrials
      logical, intent(in), optional :: start,lower_tail
      real(dp), allocatable :: cdf(:)
      logical :: lt
      integer :: k,istat

      lt = .true.
      if (present(lower_tail)) lt = lower_tail
      call frbinom2_cdf_table(ntrials,h,c,la,start,cdf,istat)
      if (istat /= 0) then
         value = 0.0_dp
         return
      end if
      if (x < 0.0_dp) then
         value = 0.0_dp
      else if (x >= real(ntrials,dp)) then
         value = 1.0_dp
      else
         k = floor(x)
         value = cdf(k)
      end if
      if (.not. lt) value = max(0.0_dp,1.0_dp-value)
   end function pfrbinom2

   subroutine pfrbinom2_vec(x,ntrials,h,c,la,value,start,lower_tail)
      real(dp), intent(in) :: x(:),h,c,la
      integer, intent(in) :: ntrials
      real(dp), intent(out) :: value(size(x))
      logical, intent(in), optional :: start,lower_tail
      real(dp), allocatable :: cdf(:)
      logical :: lt
      integer :: i,k,istat

      lt = .true.
      if (present(lower_tail)) lt = lower_tail
      call frbinom2_cdf_table(ntrials,h,c,la,start,cdf,istat)
      if (istat /= 0) then
         value = 0.0_dp
         return
      end if
      do i = 1, size(x)
         if (x(i) < 0.0_dp) then
            value(i) = 0.0_dp
         else if (x(i) >= real(ntrials,dp)) then
            value(i) = 1.0_dp
         else
            k = floor(x(i))
            value(i) = cdf(k)
         end if
      end do
      if (.not. lt) value = max(0.0_dp,1.0_dp-value)
   end subroutine pfrbinom2_vec

   integer function qfrbinom2(p,ntrials,h,c,la,start,lower_tail) result(q)
      real(dp), intent(in) :: p,h,c,la
      integer, intent(in) :: ntrials
      logical, intent(in), optional :: start,lower_tail
      real(dp), allocatable :: cdf(:)
      logical :: lt
      real(dp) :: target
      integer :: lo,hi,mid,istat

      lt = .true.
      if (present(lower_tail)) lt = lower_tail
      if (p < 0.0_dp .or. p > 1.0_dp) then
         q = -huge(1)
         return
      end if
      target = merge(p,1.0_dp-p,lt)
      call frbinom2_cdf_table(ntrials,h,c,la,start,cdf,istat)
      if (istat /= 0) then
         q = -huge(1)
         return
      end if
      if (target <= 0.0_dp) then
         q = 0
         return
      end if
      if (target >= 1.0_dp) then
         q = ntrials
         return
      end if
      lo = 0
      hi = ntrials
      do while (lo < hi)
         mid = lo+(hi-lo)/2
         if (cdf(mid) >= target) then
            hi = mid
         else
            lo = mid+1
         end if
      end do
      q = lo
   end function qfrbinom2

   subroutine qfrbinom2_vec(p,ntrials,h,c,la,q,start,lower_tail)
      real(dp), intent(in) :: p(:),h,c,la
      integer, intent(in) :: ntrials
      integer, intent(out) :: q(size(p))
      logical, intent(in), optional :: start,lower_tail
      real(dp), allocatable :: cdf(:)
      logical :: lt
      real(dp) :: target
      integer :: i,lo,hi,mid,istat

      lt = .true.
      if (present(lower_tail)) lt = lower_tail
      call frbinom2_cdf_table(ntrials,h,c,la,start,cdf,istat)
      if (istat /= 0) then
         q = -huge(1)
         return
      end if
      do i = 1, size(p)
         if (p(i) < 0.0_dp .or. p(i) > 1.0_dp) then
            q(i) = -huge(1)
            cycle
         end if
         target = merge(p(i),1.0_dp-p(i),lt)
         if (target <= 0.0_dp) then
            q(i) = 0
         else if (target >= 1.0_dp) then
            q(i) = ntrials
         else
            lo = 0
            hi = ntrials
            do while (lo < hi)
               mid = lo+(hi-lo)/2
               if (cdf(mid) >= target) then
                  hi = mid
               else
                  lo = mid+1
               end if
            end do
            q(i) = lo
         end if
      end do
   end subroutine qfrbinom2_vec

   subroutine rfrbinom2(x,ntrials,h,c,la,start)
      integer, intent(out) :: x(:)
      integer, intent(in) :: ntrials
      real(dp), intent(in) :: h,c,la
      logical, intent(in), optional :: start
      real(dp), allocatable :: cdf(:)
      real(dp) :: u
      integer :: i,lo,hi,mid,istat

      call frbinom2_cdf_table(ntrials,h,c,la,start,cdf,istat)
      if (istat /= 0) error stop "rfrbinom2: invalid parameters"
      do i = 1, size(x)
         call random_number(u)
         lo = 0
         hi = ntrials
         do while (lo < hi)
            mid = lo+(hi-lo)/2
            if (cdf(mid) >= u) then
               hi = mid
            else
               lo = mid+1
            end if
         end do
         x(i) = lo
      end do
   end subroutine rfrbinom2

   subroutine set_frbinom_seed(seed)
      integer, intent(in) :: seed
      integer :: n,i
      integer, allocatable :: state(:)
      call random_seed(size=n)
      allocate(state(n))
      do i = 1, n
         state(i) = mod(seed+104729*i,2147483646)+1
      end do
      call random_seed(put=state)
   end subroutine set_frbinom_seed

end module frbinom_distribution
