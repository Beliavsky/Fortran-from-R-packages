program test_v03_pc
   use rfast
   use zigg, only : ziggurat_rng
   implicit none
   integer, parameter :: n=500
   real(dp)::dat(n,4),e1,e2,e3,e4
   integer::i,failures
   type(ziggurat_rng)::rng
   type(pc_skeleton_result)::pc
   failures=0
   call rng%set_seed(123456789)
   do i=1,n
      e1=rng%rnorm();e2=rng%rnorm();e3=rng%rnorm();e4=rng%rnorm()
      dat(i,1)=e1
      dat(i,2)=0.8_dp*dat(i,1)+0.6_dp*e2
      dat(i,3)=0.8_dp*dat(i,2)+0.6_dp*e3
      dat(i,4)=e4
   end do
   pc=pc_skeleton(dat,PC_PEARSON,0.01_dp,1)
   call check(pc%status==0,'pc status')
   call check(pc%graph(1,2)==1.and.pc%graph(2,3)==1,'chain edges retained')
   call check(pc%graph(1,3)==0,'conditional edge removed')
   call check(all(pc%graph(4,:)==0),'independent node removed')
   call check(pc%sep_size(1,3)==1,'separator size')
   call check(pc%sepset(1,3,1)==2,'separator is middle node')
   if(failures==0)then
      print *,'test_v03_pc: PASS'
   else
      print *,'test_v03_pc: FAIL',failures
      error stop 1
   end if
contains
   subroutine check(ok,name)
      logical,intent(in)::ok
      character(*),intent(in)::name
      if(.not.ok)then;print *,'FAIL: ',trim(name);failures=failures+1;end if
   end subroutine check
end program test_v03_pc
