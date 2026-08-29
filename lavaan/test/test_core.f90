program test_core
   use lavaan
   implicit none
   real(dp), allocatable :: a(:, :), v(:), ar(:, :), d(:, :), k(:, :), sq(:, :)
   integer :: info

   a = reshape([1.0_dp, 2.0_dp, 2.0_dp, 5.0_dp], [2,2])
   v = vech(a)
   call check(maxval(abs(v-[1.0_dp,2.0_dp,5.0_dp])) < 1e-12_dp, 'vech')
   ar = vech_reverse(v,2)
   call check(maxval(abs(ar-a)) < 1e-12_dp, 'vech_reverse')
   d = duplication_matrix(2)
   call check(maxval(abs(matmul(d,v)-vec(a))) < 1e-12_dp, 'duplication identity')
   k = commutation_matrix(2,2)
   call check(maxval(abs(matmul(k,vec(a))-vec(transpose(a)))) < 1e-12_dp, 'commutation identity')
   call symmetric_sqrt(a,sq,info=info)
   call check(info==0, 'symmetric sqrt status')
   call check(maxval(abs(matmul(sq,sq)-a)) < 1e-9_dp, 'symmetric sqrt')
   print '(a)', 'test_core: PASS'
contains
   subroutine check(cond,label)
      logical,intent(in)::cond
      character(len=*),intent(in)::label
      if(.not.cond) then
      write(*,'(a,1x,a)') 'FAIL:',trim(label)
      error stop 1
      end if
   end subroutine check
end program test_core
