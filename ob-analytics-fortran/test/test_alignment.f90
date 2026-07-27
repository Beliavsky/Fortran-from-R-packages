! SPDX-License-Identifier: GPL-2.0-or-later
program test_alignment
   use ob_analytics
   implicit none
   real(dp), allocatable :: s(:,:)
   integer, allocatable :: a(:,:)
   type(event_t) :: e(5)

   s = similarity_matrix_equal([2.0_dp,4.0_dp,5.0_dp], &
      [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp])
   a = needleman_wunsch(s,-1.0_dp)
   call assert_true(size(a,1)==3,'alignment size')
   call assert_true(all(a(:,1)==[1,2,3]),'alignment left')
   call assert_true(all(a(:,2)==[2,4,5]),'alignment right')

   call init_event(e(1),1,0_i8,side_bid)
   call init_event(e(2),2,10_i8,side_ask)
   call init_event(e(3),3,10000_i8,side_bid)
   call init_event(e(4),4,10080_i8,side_bid)
   call init_event(e(5),5,10090_i8,side_ask)
   call event_match(e,1000_i8)
   call assert_true(e(1)%matching_event==2 .and. e(2)%matching_event==1,'first match')
   call assert_true(e(3)%matching_event==0,'conflicting unmatched')
   call assert_true(e(4)%matching_event==5 .and. e(5)%matching_event==4,'aligned conflict')
   print '(a)', 'test_alignment: PASS'
contains
   subroutine init_event(x,event_id,timestamp,side)
      type(event_t),intent(out)::x
      integer,intent(in)::event_id,side
      integer(i8),intent(in)::timestamp
      x%event_id=event_id
      x%id=int(event_id,i8)
      x%timestamp_ms=timestamp
      x%exchange_timestamp_ms=timestamp
      x%side=side
      x%fill=1234.0_dp
   end subroutine init_event
   subroutine assert_true(ok,label)
      logical,intent(in)::ok
      character(len=*),intent(in)::label
      if(.not.ok) then
         print '(a)', 'failed: '//label
         error stop 1
      end if
   end subroutine assert_true
end program test_alignment
