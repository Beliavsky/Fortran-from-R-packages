! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! Derived from parma 1.7, Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
module parma_timeseries
   use parma_kinds, only: dp
   implicit none
   private
   public :: lag_vector, lag_matrix, lag_matrix_set
   public :: date_to_jdn, jdn_to_date, weekday_number, is_weekday
   public :: sequence_weekdays, last_month_indices

contains

   function lag_vector(x, nlag, pad) result(y)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: nlag
      real(dp), intent(in), optional :: pad
      real(dp) :: y(size(x))
      real(dp) :: fill
      integer :: n

      fill = 0.0_dp
      if (present(pad)) fill = pad
      n = size(x)
      if (nlag <= 0) then
         y = x
      else if (nlag >= n) then
         y = fill
      else
         y(1:nlag) = fill
         y(nlag+1:n) = x(1:n-nlag)
      end if
   end function lag_vector

   function lag_matrix(x, nlag, pad) result(y)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: nlag
      real(dp), intent(in), optional :: pad
      real(dp) :: y(size(x,1),size(x,2))
      real(dp) :: fill
      integer :: n

      fill = 0.0_dp
      if (present(pad)) fill = pad
      n = size(x,1)
      if (nlag <= 0) then
         y = x
      else if (nlag >= n) then
         y = fill
      else
         y(1:nlag,:) = fill
         y(nlag+1:n,:) = x(1:n-nlag,:)
      end if
   end function lag_matrix

   function lag_matrix_set(x, lags, pad) result(y)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: lags(:)
      real(dp), intent(in), optional :: pad
      real(dp) :: y(size(x,1),size(x,2)*size(lags))
      integer :: i, j1, j2

      do i = 1, size(lags)
         j1 = 1 + (i-1)*size(x,2)
         j2 = i*size(x,2)
         if (present(pad)) then
            y(:,j1:j2) = lag_matrix(x,lags(i),pad)
         else
            y(:,j1:j2) = lag_matrix(x,lags(i))
         end if
      end do
   end function lag_matrix_set

   pure function date_to_jdn(year, month, day) result(jdn)
      integer, intent(in) :: year, month, day
      integer :: jdn
      integer :: a, y, m

      a = (14-month)/12
      y = year + 4800 - a
      m = month + 12*a - 3
      jdn = day + (153*m+2)/5 + 365*y + y/4 - y/100 + y/400 - 32045
   end function date_to_jdn

   pure subroutine jdn_to_date(jdn, year, month, day)
      integer, intent(in) :: jdn
      integer, intent(out) :: year, month, day
      integer :: a, b, c, d, e, m

      a = jdn + 32044
      b = (4*a+3)/146097
      c = a - (146097*b)/4
      d = (4*c+3)/1461
      e = c - (1461*d)/4
      m = (5*e+2)/153
      day = e - (153*m+2)/5 + 1
      month = m + 3 - 12*(m/10)
      year = 100*b + d - 4800 + m/10
   end subroutine jdn_to_date

   pure function weekday_number(year, month, day) result(weekday)
      integer, intent(in) :: year, month, day
      integer :: weekday
      weekday = modulo(date_to_jdn(year,month,day),7) + 1
   end function weekday_number

   pure function is_weekday(year, month, day) result(answer)
      integer, intent(in) :: year, month, day
      logical :: answer
      integer :: w
      w = weekday_number(year,month,day)
      answer = w <= 5
   end function is_weekday

   subroutine sequence_weekdays(start_date, end_date, dates, nout, info)
      character(len=*), intent(in) :: start_date, end_date
      character(len=10), allocatable, intent(out) :: dates(:)
      integer, intent(out) :: nout, info
      integer :: ys, ms, ds, ye, me, de, js, je, j, y, m, d, n
      character(len=10), allocatable :: temp(:)

      info = 0
      call parse_iso_date(start_date,ys,ms,ds,info)
      if (info /= 0) return
      call parse_iso_date(end_date,ye,me,de,info)
      if (info /= 0) return
      js = date_to_jdn(ys,ms,ds)
      je = date_to_jdn(ye,me,de)
      if (je < js) then
         info = 2
         nout = 0
         allocate(dates(0))
         return
      end if
      allocate(temp(je-js+1))
      n = 0
      do j = js, je
         call jdn_to_date(j,y,m,d)
         if (is_weekday(y,m,d)) then
            n = n + 1
            write(temp(n),'(i4.4,"-",i2.2,"-",i2.2)') y,m,d
         end if
      end do
      nout = n
      allocate(dates(n))
      if (n > 0) dates = temp(1:n)
   end subroutine sequence_weekdays

   subroutine last_month_indices(dates, indices, nout, info)
      character(len=*), intent(in) :: dates(:)
      integer, allocatable, intent(out) :: indices(:)
      integer, intent(out) :: nout, info
      integer, allocatable :: temp(:)
      integer :: i, y1, m1, d1, y2, m2, d2

      info = 0
      allocate(temp(size(dates)))
      nout = 0
      if (size(dates) < 2) then
         allocate(indices(0))
         return
      end if
      call parse_iso_date(dates(1),y1,m1,d1,info)
      if (info /= 0) return
      do i = 2, size(dates)
         call parse_iso_date(dates(i),y2,m2,d2,info)
         if (info /= 0) return
         if (m2 /= m1 .or. y2 /= y1) then
            nout = nout + 1
            temp(nout) = i - 1
         end if
         y1 = y2
         m1 = m2
         d1 = d2
      end do
      allocate(indices(nout))
      if (nout > 0) indices = temp(1:nout)
   end subroutine last_month_indices

   subroutine parse_iso_date(text, year, month, day, info)
      character(len=*), intent(in) :: text
      integer, intent(out) :: year, month, day, info
      integer :: ios
      character(len=10) :: s

      s = text(1:min(10,len_trim(text)))
      read(s(1:4),*,iostat=ios) year
      if (ios /= 0) then
         info = 1
         return
      end if
      read(s(6:7),*,iostat=ios) month
      if (ios /= 0) then
         info = 1
         return
      end if
      read(s(9:10),*,iostat=ios) day
      if (ios /= 0) then
         info = 1
         return
      end if
      info = 0
   end subroutine parse_iso_date

end module parma_timeseries
