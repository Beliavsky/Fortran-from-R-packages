! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of obAnalytics.
! Copyright (C) 2015,2016 Philip Stubbings.
module ob_io
   use ob_kinds, only : dp, i8
   use ob_types, only : event_t, action_created, action_changed, action_deleted, side_bid, side_ask
   use ob_utils, only : round_digits
   implicit none
   private
   public :: read_event_csv, write_event_csv

contains

   subroutine read_event_csv(filename, events, price_digits, volume_digits, status, message)
      character(len=*), intent(in) :: filename
      type(event_t), allocatable, intent(out) :: events(:)
      integer, intent(in), optional :: price_digits, volume_digits
      integer, intent(out), optional :: status
      character(len=:), allocatable, intent(out), optional :: message
      type(event_t), allocatable :: raw(:), clean(:)
      character(len=4096) :: line
      character(len=256) :: fields(7)
      integer :: unit, ios, nlines, i, n, pd, vd
      logical :: duplicate

      pd = 2
      vd = 8
      if (present(price_digits)) pd = price_digits
      if (present(volume_digits)) vd = volume_digits
      if (present(status)) status = 0
      if (present(message)) message = ''

      open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         allocate(events(0))
         call report_error(1, 'read_event_csv: cannot open file')
         return
      end if
      read(unit,'(a)',iostat=ios) line
      nlines = 0
      do
         read(unit,'(a)',iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) > 0) nlines = nlines + 1
      end do
      rewind(unit)
      read(unit,'(a)',iostat=ios) line
      allocate(raw(nlines))
      n = 0
      do i = 1, nlines
         read(unit,'(a)',iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         call split_csv_line(line, fields, ios)
         if (ios /= 0) then
            close(unit)
            allocate(events(0))
            call report_error(2, 'read_event_csv: malformed CSV row')
            return
         end if
         n = n + 1
         read(fields(1),*,iostat=ios) raw(n)%id
         if (ios /= 0) exit
         read(fields(2),*,iostat=ios) raw(n)%timestamp_ms
         if (ios /= 0) exit
         read(fields(3),*,iostat=ios) raw(n)%exchange_timestamp_ms
         if (ios /= 0) exit
         read(fields(4),*,iostat=ios) raw(n)%price
         if (ios /= 0) exit
         read(fields(5),*,iostat=ios) raw(n)%volume
         if (ios /= 0) exit
         raw(n)%action = parse_action(fields(6))
         raw(n)%side = parse_side(fields(7))
         raw(n)%price = round_digits(raw(n)%price, pd)
         raw(n)%volume = round_digits(raw(n)%volume, vd)
      end do
      close(unit)
      if (ios > 0 .or. n /= nlines) then
         allocate(events(0))
         call report_error(3, 'read_event_csv: invalid numeric or categorical field')
         return
      end if

      allocate(clean(n))
      nlines = 0
      do i = 1, n
         if (raw(i)%volume < 0.0_dp) cycle
         duplicate = .false.
         if (raw(i)%action /= action_changed) then
            duplicate = duplicate_nonchange(raw, i)
         end if
         if (duplicate) cycle
         nlines = nlines + 1
         clean(nlines) = raw(i)
      end do
      if (nlines < size(clean)) clean = clean(1:nlines)
      call sort_events_lifecycle(clean)
      do i = 1, size(clean)
         clean(i)%event_id = i
         if (i == 1) then
            clean(i)%fill = 0.0_dp
         else if (clean(i)%id /= clean(i-1)%id) then
            clean(i)%fill = 0.0_dp
         else
            clean(i)%fill = round_digits(abs(clean(i)%volume-clean(i-1)%volume), vd)
         end if
      end do
      call move_alloc(clean, events)

   contains
      subroutine report_error(code, text)
         integer, intent(in) :: code
         character(len=*), intent(in) :: text
         if (present(status)) status = code
         if (present(message)) message = text
      end subroutine report_error
   end subroutine read_event_csv

   subroutine write_event_csv(filename, events, status)
      character(len=*), intent(in) :: filename
      type(event_t), intent(in) :: events(:)
      integer, intent(out), optional :: status
      integer :: unit, ios, i
      character(len=16) :: action, side
      if (present(status)) status = 0
      open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
      if (ios /= 0) then
         if (present(status)) status = ios
         return
      end if
      write(unit,'(a)') 'id,timestamp,exchange.timestamp,price,volume,action,direction'
      do i = 1, size(events)
         action = action_string(events(i)%action)
         side = side_string(events(i)%side)
         write(unit,'(i0,a,i0,a,i0,a,es24.16,a,es24.16,a,a,a,a)') events(i)%id, ',', &
            events(i)%timestamp_ms, ',', events(i)%exchange_timestamp_ms, ',', events(i)%price, ',', &
            events(i)%volume, ',', trim(action), ',', trim(side)
      end do
      close(unit)
   end subroutine write_event_csv

   subroutine split_csv_line(line, fields, status)
      character(len=*), intent(in) :: line
      character(len=*), intent(out) :: fields(:)
      integer, intent(out) :: status
      integer :: i, start, field, n
      fields = ''
      start = 1
      field = 1
      n = len_trim(line)
      do i = 1, n+1
         if (i == n+1 .or. line(i:i) == ',') then
            if (field > size(fields)) then
               status = 1
               return
            end if
            if (i > start) fields(field) = adjustl(line(start:i-1))
            field = field + 1
            start = i + 1
         end if
      end do
      if (field-1 /= size(fields)) then
         status = 1
      else
         status = 0
      end if
   end subroutine split_csv_line

   pure integer function parse_action(text)
      character(len=*), intent(in) :: text
      character(len=:), allocatable :: word
      word = lowercase(trim(adjustl(text)))
      select case (word)
      case ('created'); parse_action = action_created
      case ('changed','modified'); parse_action = action_changed
      case ('deleted'); parse_action = action_deleted
      case default; parse_action = 0
      end select
   end function parse_action

   pure integer function parse_side(text)
      character(len=*), intent(in) :: text
      character(len=:), allocatable :: word
      word = lowercase(trim(adjustl(text)))
      select case (word)
      case ('bid'); parse_side = side_bid
      case ('ask'); parse_side = side_ask
      case default; parse_side = 0
      end select
   end function parse_side

   pure function lowercase(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i, code
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            lower(i:i) = achar(code + iachar('a') - iachar('A'))
         else
            lower(i:i) = text(i:i)
         end if
      end do
   end function lowercase

   pure logical function duplicate_nonchange(raw, current)
      type(event_t), intent(in) :: raw(:)
      integer, intent(in) :: current
      integer :: j
      duplicate_nonchange = .false.
      do j = 1, current-1
         if (raw(j)%action == action_changed) cycle
         if (raw(j)%id == raw(current)%id .and. raw(j)%price == raw(current)%price .and. &
             raw(j)%volume == raw(current)%volume .and. raw(j)%action == raw(current)%action) then
            duplicate_nonchange = .true.
            return
         end if
      end do
   end function duplicate_nonchange

   subroutine sort_events_lifecycle(events)
      type(event_t), intent(inout) :: events(:)
      type(event_t) :: key
      integer :: i, j
      do i = 2, size(events)
         key = events(i)
         j = i - 1
         do while (j >= 1)
            if (event_before(events(j), key)) exit
            events(j+1) = events(j)
            j = j - 1
         end do
         events(j+1) = key
      end do
   end subroutine sort_events_lifecycle

   pure logical function event_before(a, b)
      type(event_t), intent(in) :: a, b
      event_before = a%id < b%id .or. &
         (a%id == b%id .and. a%action < b%action) .or. &
         (a%id == b%id .and. a%action == b%action .and. a%timestamp_ms <= b%timestamp_ms)
   end function event_before

   pure function action_string(action) result(text)
      integer, intent(in) :: action
      character(len=16) :: text
      select case (action)
      case (action_created); text = 'created'
      case (action_changed); text = 'changed'
      case (action_deleted); text = 'deleted'
      case default; text = 'unknown'
      end select
   end function action_string

   pure function side_string(side) result(text)
      integer, intent(in) :: side
      character(len=16) :: text
      if (side == side_bid) then
         text = 'bid'
      else
         text = 'ask'
      end if
   end function side_string

end module ob_io
