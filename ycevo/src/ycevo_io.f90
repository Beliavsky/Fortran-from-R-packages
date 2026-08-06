module ycevo_io
   use ycevo_kinds, only : dp
   use ycevo_status, only : ycevo_success, ycevo_err_input
   use ycevo_types, only : bond_panel_t, yield_curve_t
   implicit none
   private

   public :: read_bond_panel_csv, write_yield_curve_csv

contains

   subroutine read_bond_panel_csv(filename, panel, status, message, has_header)
      character(len=*), intent(in) :: filename
      type(bond_panel_t), intent(out) :: panel
      integer, intent(out) :: status
      character(len=*), intent(out), optional :: message
      logical, intent(in), optional :: has_header
      character(len=1024) :: line
      integer :: unit, ios, n, i, day, id, tupq
      real(dp) :: price, cashflow
      logical :: header

      status = ycevo_err_input
      if (present(message)) message = ''
      header = .true.
      if (present(has_header)) header = has_header
      open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         if (present(message)) message = 'Unable to open bond-panel CSV.'
         return
      end if
      if (header) read(unit, '(a)', iostat=ios) line
      n = 0
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) > 0) n = n + 1
      end do
      if (n == 0) then
         close(unit)
         if (present(message)) message = 'Bond-panel CSV contains no data rows.'
         return
      end if
      rewind(unit)
      if (header) read(unit, '(a)') line
      allocate(panel%day(n), panel%id(n), panel%tupq(n), panel%price(n), panel%cashflow(n))
      do i = 1, n
         read(unit, *, iostat=ios) day, id, price, tupq, cashflow
         if (ios /= 0) then
            close(unit)
            if (present(message)) write(message, '(a,i0)') 'Invalid CSV row ', i
            return
         end if
         panel%day(i) = day
         panel%id(i) = id
         panel%price(i) = price
         panel%tupq(i) = tupq
         panel%cashflow(i) = cashflow
      end do
      close(unit)
      panel%nday = maxval(panel%day)
      call panel%sort()
      status = ycevo_success
   end subroutine read_bond_panel_csv

   subroutine write_yield_curve_csv(filename, curve, status)
      character(len=*), intent(in) :: filename
      type(yield_curve_t), intent(in) :: curve
      integer, intent(out) :: status
      integer :: unit, ios, i

      status = ycevo_err_input
      if (.not. allocated(curve%tau) .or. .not. allocated(curve%discount) .or. &
          .not. allocated(curve%yield)) return
      open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
      if (ios /= 0) return
      write(unit, '(a)') 'xgrid,tau,discount,yield'
      do i = 1, size(curve%tau)
         write(unit, '(es24.16,a,es24.16,a,es24.16,a,es24.16)') &
            curve%xgrid, ',', curve%tau(i), ',', curve%discount(i), ',', curve%yield(i)
      end do
      close(unit)
      status = ycevo_success
   end subroutine write_yield_curve_csv

end module ycevo_io
