program test_reference_data
   use cec2005benchmark
   implicit none
   integer :: fid, u, i, ios
   real(dp) :: x(10,50), ref(10), got(10), err, worst
   character(len=256) :: file_name
   character(len=2) :: sid

   worst = 0.0_dp
   do fid = 1, 25
      write(sid,'(i0)') fid
      file_name = 'data/test_data_func'//trim(sid)//'.txt'
      open(newunit=u, file=trim(file_name), status='old', action='read')
      do i = 1, 10
         read(u,*) x(i,:)
      end do
      do i = 1, 10
         read(u,*) ref(i)
      end do
      close(u)

      call cec2005_eval_batch(fid, x, got, 'data', .false., ios)
      if (ios /= 0) error stop 'failed to initialize reference-data case'
      do i = 1, 10
         err = abs(got(i)-ref(i))/max(1.0_dp,abs(ref(i)))
         worst = max(worst,err)
         if (err > 1.0e-8_dp) error stop 'reference-data mismatch'
      end do
   end do
   print '(a,es12.4)', 'PASS test_reference_data; worst relative error = ', worst
end program test_reference_data
