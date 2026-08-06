program test_transition
   use tvgarch, only : dp, tv_spec, make_tv_spec, tv_transition, tv_component, &
                       combinations_binary, pack_tv_parameters, unpack_tv_parameters
   implicit none
   real(dp), allocatable :: x(:),v(:),g(:),par(:)
   integer, allocatable :: b(:,:)
   type(tv_spec) :: spec,spec2
   integer :: st
   x=[0.0_dp,0.5_dp,1.0_dp]
   call tv_transition(0.0_dp,[0.5_dp],x,v,st,0)
   call check(st==0,'transition status')
   call check(maxval(abs(v-0.5_dp))<1e-14_dp,'zero speed logistic')
   call make_tv_spec(spec,[1],2.0_dp,[3.0_dp],[log(10.0_dp)],[0.5_dp],2)
   call tv_component(spec,x,g,st)
   call check(st==0,'component status')
   call check(abs(g(2)-3.5_dp)<1e-12_dp,'component midpoint')
   call pack_tv_parameters(spec,par)
   call unpack_tv_parameters([1],par,2,spec2,st)
   call check(st==0 .and. maxval(abs(spec2%sizes-spec%sizes))<1e-14_dp,'pack roundtrip')
   call combinations_binary(3,b,st)
   call check(st==0 .and. size(b,1)==7 .and. all(sum(b,dim=2)>0),'combinations')
   print '(a)', 'test_transition: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok;character(len=*),intent(in)::msg
      if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
   end subroutine check
end program test_transition
