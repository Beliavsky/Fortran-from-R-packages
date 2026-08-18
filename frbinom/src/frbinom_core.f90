module frbinom_core
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use frbinom_kinds, only : dp
   implicit none
   private

   public :: frbinom_pmf_table, frbinom_cdf_table
   public :: frbinom2_pmf_table, frbinom2_cdf_table
   public :: frbinom_cmax, frbinom_parameters_ok, frbinom2_parameters_ok

contains

   pure real(dp) function frbinom_cmax(prob,h) result(cmax)
      real(dp), intent(in) :: prob,h
      real(dp) :: a,disc
      a = 2.0_dp**(2.0_dp*h-2.0_dp)
      disc = 4.0_dp*prob-prob*2.0_dp**(2.0_dp*h)+2.0_dp**(4.0_dp*h-4.0_dp)
      if (disc < 0.0_dp) then
         cmax = -1.0_dp
      else
         cmax = min(0.5_dp*(-2.0_dp*prob+a+sqrt(disc)),1.0_dp-prob)
      end if
   end function frbinom_cmax

   pure logical function frbinom_parameters_ok(size,prob,h,c) result(ok)
      integer, intent(in) :: size
      real(dp), intent(in) :: prob,h,c
      real(dp) :: cmax
      cmax = frbinom_cmax(prob,h)
      ok = size >= 1 .and. prob >= 0.0_dp .and. prob <= 1.0_dp .and. &
           h >= 0.0_dp .and. h <= 1.0_dp .and. c >= 0.0_dp .and. c <= cmax
   end function frbinom_parameters_ok

   pure logical function frbinom2_parameters_ok(size,h,c,la) result(ok)
      integer, intent(in) :: size
      real(dp), intent(in) :: h,c,la
      ok = size >= 1 .and. h >= 0.5_dp .and. h <= 1.0_dp .and. &
           c > 0.0_dp .and. c < 2.0_dp**(2.0_dp*h-2.0_dp) .and. &
           la > 0.0_dp .and. la < c
   end function frbinom2_parameters_ok

   subroutine waiting_frbinom(size,prob,h,c,wait)
      integer, intent(in) :: size
      real(dp), intent(in) :: prob,h,c
      real(dp), intent(out) :: wait(0:size-1)
      integer :: i,j
      real(dp) :: d

      wait(0) = prob+c
      do i = 2, size
         d = 0.0_dp
         do j = 1, i-1
            d = d+(prob+c*real(i-j,dp)**(2.0_dp*h-2.0_dp))*wait(j-1)
         end do
         wait(i-1) = prob+c*real(i,dp)**(2.0_dp*h-2.0_dp)-d
      end do
   end subroutine waiting_frbinom

   subroutine initial_frbinom(size,prob,h,c,first)
      integer, intent(in) :: size
      real(dp), intent(in) :: prob,h,c
      real(dp), intent(out) :: first(0:size-1)
      real(dp), allocatable :: wait(:)
      integer :: i,j
      real(dp) :: d

      allocate(wait(0:size-1))
      call waiting_frbinom(size,prob,h,c,wait)
      first(0) = prob
      do i = 2, size
         d = 0.0_dp
         do j = 1, i-1
            d = d+(prob+c*real(i-j,dp)**(2.0_dp*h-2.0_dp))*first(j-1)
         end do
         first(i-1) = prob-d
      end do
   end subroutine initial_frbinom

   subroutine waiting_frbinom2(size,h,c,wait)
      integer, intent(in) :: size
      real(dp), intent(in) :: h,c
      real(dp), intent(out) :: wait(0:size-1)
      integer :: i,j
      real(dp) :: d

      wait(0) = c
      do i = 2, size
         d = 0.0_dp
         do j = 1, i-1
            d = d+c*real(i-j,dp)**(2.0_dp*h-2.0_dp)*wait(j-1)
         end do
         wait(i-1) = c*real(i,dp)**(2.0_dp*h-2.0_dp)-d
      end do
   end subroutine waiting_frbinom2

   subroutine initial_frbinom2(size,h,c,la,first)
      integer, intent(in) :: size
      real(dp), intent(in) :: h,c,la
      real(dp), intent(out) :: first(0:size-1)
      integer :: i,j
      real(dp) :: d,pn

      pn = la*real(size,dp)**(2.0_dp*h-2.0_dp)
      first(0) = pn
      do i = 2, size
         d = 0.0_dp
         do j = 1, i-1
            d = d+c*real(i-j,dp)**(2.0_dp*h-2.0_dp)*first(j-1)
         end do
         first(i-1) = pn-d
      end do
   end subroutine initial_frbinom2

   subroutine count_distribution(size,first,wait,pmf)
      integer, intent(in) :: size
      real(dp), intent(in) :: first(0:size-1),wait(0:size-1)
      real(dp), intent(out) :: pmf(0:size)
      real(dp), allocatable :: a(:,:),surv(:)
      real(dp) :: s
      integer :: k,r,j

      allocate(a(0:size-1,0:size-1),surv(0:size-1))
      a = 0.0_dp
      a(0,:) = first

      surv(0) = 1.0_dp-wait(0)
      do k = 1, size-1
         surv(k) = surv(k-1)-wait(k)
      end do

      do k = 1, size-1
         do r = 1, k
            s = 0.0_dp
            do j = 0, k-1
               s = s+a(r-1,j)*wait(k-j-1)
            end do
            a(r,k) = s
         end do
      end do

      pmf = 0.0_dp
      pmf(0) = 1.0_dp-sum(first)
      do r = 0, size-1
         pmf(r+1) = a(r,size-1)
         if (r <= size-2) then
            s = 0.0_dp
            do j = 0, size-2
               s = s+a(r,j)*surv(size-2-j)
            end do
            pmf(r+1) = pmf(r+1)+s
         end if
      end do

      ! Roundoff in the renewal recurrences can create tiny signed zeros.
      where (pmf < 0.0_dp .and. abs(pmf) <= 100.0_dp*epsilon(1.0_dp))
         pmf = 0.0_dp
      end where
   end subroutine count_distribution

   subroutine frbinom_pmf_table(size,prob,h,c,start,pmf,status)
      integer, intent(in) :: size
      real(dp), intent(in) :: prob,h,c
      logical, intent(in), optional :: start
      real(dp), allocatable, intent(out) :: pmf(:)
      integer, intent(out), optional :: status
      logical :: st
      real(dp), allocatable :: wait(:),first(:)

      st = .false.
      if (present(start)) st = start
      if (present(status)) status = 0

      if (.not. frbinom_parameters_ok(size,prob,h,c)) then
         allocate(pmf(0))
         if (present(status)) status = 1
         return
      end if

      allocate(pmf(0:size))
      if (size == 1) then
         ! The upstream dfrbinom() explicitly reduces to Binomial(1,prob)
         ! regardless of start/c for this special case.
         pmf(0) = 1.0_dp-prob
         pmf(1) = prob
         return
      end if

      allocate(wait(0:size-1),first(0:size-1))
      call waiting_frbinom(size,prob,h,c,wait)
      if (st) then
         first = wait
      else
         call initial_frbinom(size,prob,h,c,first)
      end if
      call count_distribution(size,first,wait,pmf)
   end subroutine frbinom_pmf_table

   subroutine frbinom_cdf_table(size,prob,h,c,start,cdf,status)
      integer, intent(in) :: size
      real(dp), intent(in) :: prob,h,c
      logical, intent(in), optional :: start
      real(dp), allocatable, intent(out) :: cdf(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: pmf(:)
      integer :: k,istat

      call frbinom_pmf_table(size,prob,h,c,start,pmf,istat)
      if (istat /= 0) then
         allocate(cdf(0))
         if (present(status)) status = istat
         return
      end if
      allocate(cdf(0:size))
      cdf(0) = pmf(0)
      do k = 1, size
         cdf(k) = cdf(k-1)+pmf(k)
      end do
      cdf(size) = 1.0_dp
      if (present(status)) status = 0
   end subroutine frbinom_cdf_table

   subroutine frbinom2_pmf_table(size,h,c,la,start,pmf,status)
      integer, intent(in) :: size
      real(dp), intent(in) :: h,c,la
      logical, intent(in), optional :: start
      real(dp), allocatable, intent(out) :: pmf(:)
      integer, intent(out), optional :: status
      logical :: st
      real(dp), allocatable :: wait(:),first(:)
      real(dp) :: pn

      st = .false.
      if (present(start)) st = start
      if (present(status)) status = 0

      if (.not. frbinom2_parameters_ok(size,h,c,la)) then
         allocate(pmf(0))
         if (present(status)) status = 1
         return
      end if

      allocate(pmf(0:size))
      if (size == 1) then
         if (st) then
            pn = c
         else
            pn = la
         end if
         pmf(0) = 1.0_dp-pn
         pmf(1) = pn
         return
      end if

      allocate(wait(0:size-1),first(0:size-1))
      call waiting_frbinom2(size,h,c,wait)
      if (st) then
         first = wait
      else
         call initial_frbinom2(size,h,c,la,first)
      end if
      call count_distribution(size,first,wait,pmf)
   end subroutine frbinom2_pmf_table

   subroutine frbinom2_cdf_table(size,h,c,la,start,cdf,status)
      integer, intent(in) :: size
      real(dp), intent(in) :: h,c,la
      logical, intent(in), optional :: start
      real(dp), allocatable, intent(out) :: cdf(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: pmf(:)
      integer :: k,istat

      call frbinom2_pmf_table(size,h,c,la,start,pmf,istat)
      if (istat /= 0) then
         allocate(cdf(0))
         if (present(status)) status = istat
         return
      end if
      allocate(cdf(0:size))
      cdf(0) = pmf(0)
      do k = 1, size
         cdf(k) = cdf(k-1)+pmf(k)
      end do
      cdf(size) = 1.0_dp
      if (present(status)) status = 0
   end subroutine frbinom2_cdf_table

end module frbinom_core
