module ycevo_types
   use ycevo_kinds, only : dp
   use ycevo_status, only : ycevo_success, ycevo_err_input
   implicit none
   private

   type, public :: bond_panel_t
      integer, allocatable :: day(:)
      integer, allocatable :: id(:)
      integer, allocatable :: tupq(:)
      real(dp), allocatable :: price(:)
      real(dp), allocatable :: cashflow(:)
      integer :: nday = 0
   contains
      procedure :: size => panel_size
      procedure :: validate => validate_panel
      procedure :: sort => sort_panel
   end type bond_panel_t

   type, public :: yield_curve_t
      real(dp) :: xgrid = 0.0_dp
      real(dp), allocatable :: tau(:)
      real(dp), allocatable :: discount(:)
      real(dp), allocatable :: yield(:)
   end type yield_curve_t

   type, public :: yield_surface_t
      real(dp), allocatable :: xgrid(:)
      real(dp), allocatable :: tau(:)
      real(dp), allocatable :: discount(:, :)
      real(dp), allocatable :: yield(:, :)
   end type yield_surface_t

contains

   integer function panel_size(self) result(n)
      class(bond_panel_t), intent(in) :: self

      if (allocated(self%day)) then
         n = size(self%day)
      else
         n = 0
      end if
   end function panel_size

   subroutine validate_panel(self, status, message)
      class(bond_panel_t), intent(in) :: self
      integer, intent(out) :: status
      character(len=*), intent(out), optional :: message
      integer :: n

      status = ycevo_err_input
      if (present(message)) message = ''
      if (.not. allocated(self%day) .or. .not. allocated(self%id) .or. &
          .not. allocated(self%tupq) .or. .not. allocated(self%price) .or. &
          .not. allocated(self%cashflow)) then
         if (present(message)) message = 'All bond-panel arrays must be allocated.'
         return
      end if
      n = size(self%day)
      if (size(self%id) /= n .or. size(self%tupq) /= n .or. &
          size(self%price) /= n .or. size(self%cashflow) /= n) then
         if (present(message)) message = 'Bond-panel arrays must have equal length.'
         return
      end if
      if (n == 0 .or. self%nday <= 0) then
         if (present(message)) message = 'Bond panel cannot be empty.'
         return
      end if
      if (any(self%day < 1) .or. any(self%day > self%nday) .or. &
          any(self%tupq < 1)) then
         if (present(message)) message = 'Invalid day index or time to payment.'
         return
      end if
      status = ycevo_success
   end subroutine validate_panel

   subroutine sort_panel(self)
      class(bond_panel_t), intent(inout) :: self
      integer, allocatable :: order(:), temp(:)
      integer :: i, n

      n = self%size()
      if (n <= 1) return
      allocate(order(n), temp(n))
      order = [(i, i=1,n)]
      call merge_sort(order, temp, 1, n)
      self%day = self%day(order)
      self%id = self%id(order)
      self%tupq = self%tupq(order)
      self%price = self%price(order)
      self%cashflow = self%cashflow(order)

   contains

      recursive subroutine merge_sort(idx, work, left, right)
         integer, intent(inout) :: idx(:), work(:)
         integer, intent(in) :: left, right
         integer :: mid, a, b, k

         if (left >= right) return
         mid = (left + right) / 2
         call merge_sort(idx, work, left, mid)
         call merge_sort(idx, work, mid + 1, right)
         a = left
         b = mid + 1
         do k = left, right
            if (a > mid) then
               work(k) = idx(b)
               b = b + 1
            else if (b > right) then
               work(k) = idx(a)
               a = a + 1
            else if (comes_first(idx(a), idx(b))) then
               work(k) = idx(a)
               a = a + 1
            else
               work(k) = idx(b)
               b = b + 1
            end if
         end do
         idx(left:right) = work(left:right)
      end subroutine merge_sort

      logical function comes_first(a, b)
         integer, intent(in) :: a, b

         comes_first = self%day(a) < self%day(b) .or. &
            (self%day(a) == self%day(b) .and. self%id(a) < self%id(b)) .or. &
            (self%day(a) == self%day(b) .and. self%id(a) == self%id(b) .and. &
             self%tupq(a) <= self%tupq(b))
      end function comes_first
   end subroutine sort_panel

end module ycevo_types
