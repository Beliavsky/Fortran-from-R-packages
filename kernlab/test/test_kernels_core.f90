! SPDX-License-Identifier: GPL-2.0-only
program test_kernels_core
  use kernlab
  implicit none
  real(dp) :: x(3,2), val
  real(dp), allocatable :: k(:,:), km(:,:), estimates(:), probs(:)
  integer :: st
  type(kernel_spec) :: kr, kp, ks
  type(inchol_result) :: ic
  type(ipop_result) :: qp
  real(dp) :: c(2), h(2,2), a(1,2), b(1), l(2), u(2), r(1)
  character(len=8) :: words(3)

  x = reshape([0.0_dp,0.0_dp, 1.0_dp,0.0_dp, 0.0_dp,1.0_dp],[3,2],order=[2,1])
  kr = rbfdot(0.5_dp)
  call check(abs(kernel_value(kr,x(1,:),x(1,:))-1.0_dp)<1.0e-12_dp,'rbf diagonal')
  call check(abs(kernel_value(kr,x(1,:),x(2,:))-exp(-0.5_dp))<1.0e-12_dp,'rbf formula')
  kp = polydot(2,2.0_dp,1.0_dp)
  call check(abs(kernel_value(kp,[1.0_dp,2.0_dp],[2.0_dp,1.0_dp])-81.0_dp)<1.0e-12_dp,'poly')
  call kernel_matrix(kr,x,k,st); call check(st==KL_SUCCESS.and.size(k,1)==3,'kernel matrix')
  allocate(km(3,1)); km(:,1)=[1.0_dp,2.0_dp,3.0_dp]
  call kernel_pol(kr,x,km(:,1),val,st); call check(st==KL_SUCCESS.and.val>0.0_dp,'kernel polynomial')

  words=['abcd    ','abce    ','zzzz    '];ks=stringdot(2,normalized=.true.)
  call string_kernel_matrix(ks,words,k,st)
  call check(st==KL_SUCCESS.and.abs(k(1,1)-1.0_dp)<1.0e-12_dp.and.k(1,2)>0.0_dp,'string spectrum')

  call sigest(x,estimates,st,scaled=.false.)
  call check(st==KL_SUCCESS.and.all(estimates>0.0_dp),'sigest')
  call inchol(x,kr,ic,tol=1.0e-10_dp)
  call check(ic%status==KL_SUCCESS.and.ic%rank>=2,'incomplete Cholesky')

  h=reshape([2.0_dp,0.0_dp,0.0_dp,2.0_dp],[2,2]); c=[-2.0_dp,-4.0_dp]
  a=reshape([1.0_dp,1.0_dp],[1,2]); b=[1.0_dp]; r=[0.0_dp]; l=[0.0_dp,0.0_dp]; u=[2.0_dp,2.0_dp]
  call ipop(c,h,a,b,l,u,r,qp,maxiter=5000,tol=1.0e-8_dp,penalty=1.0e4_dp)
  call check(qp%primal_infeasibility<2.0e-3_dp.and.abs(sum(qp%primal)-1.0_dp)<2.0e-3_dp,'ipop constraint')

  call couple([0.8_dp,0.7_dp,0.6_dp],probs,st,'minpair')
  call check(st==KL_SUCCESS.and.abs(sum(probs)-1.0_dp)<1.0e-10_dp.and.all(probs>=0.0_dp),'coupling')
  print '(a)', 'test_kernels_core: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
  end subroutine check
end program test_kernels_core
