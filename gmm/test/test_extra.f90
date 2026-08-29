program test_extra
use gmm
implicit none
real(dp) :: x(4,5), theta(3), f(3), qt(3,2), vff(3,3), vtf(6,3)
real(dp), allocatable :: gt(:,:), grad(:,:)
type(kleibergen_result_t) :: kr
real(dp) :: tol, u(80,2), gg(2,1), b1, b2
integer :: i

tol=2.0e-10_dp
x(:,1)=[1.0_dp,2.0_dp,2.5_dp,4.0_dp]
x(:,2)=1.0_dp
x(:,3)=[0.0_dp,1.0_dp,0.0_dp,1.0_dp]
x(:,4)=[1.0_dp,2.0_dp,3.0_dp,4.0_dp]
x(:,5)=[0.5_dp,0.2_dp,0.7_dp,0.9_dp]
theta=[1.2_dp,0.8_dp,0.4_dp]
call ate_moments(theta,x,2,ATE_BAL,ATE_LINEAR,gt)
call check(maxval(abs(gt(:,1)-[-0.2_dp,0.0_dp,1.3_dp,2.0_dp]))<tol,'ATE residual moments')
call check(maxval(abs(gt(:,2)-[0.0_dp,0.0_dp,0.0_dp,2.0_dp]))<tol,'ATE treatment moments')
call check(maxval(abs(gt(:,3)-[-0.4_dp,0.6_dp,-0.4_dp,0.6_dp]))<tol,'ATE balance moments')
call check(maxval(abs(gt(:,4)-[-0.4_dp,1.2_dp,-1.2_dp,2.4_dp]))<tol,'ATE covariate 1 moments')
call check(maxval(abs(gt(:,5)-[-0.2_dp,0.12_dp,-0.28_dp,0.54_dp]))<tol,'ATE covariate 2 moments')
call ate_gradient(theta,x,2,ATE_BAL,ATE_LINEAR,grad)
call check(abs(grad(1,1)+1.0_dp)<tol .and. abs(grad(1,2)+0.5_dp)<tol,'ATE gradient top row')
call check(abs(grad(2,1)+0.5_dp)<tol .and. abs(grad(2,2)+0.5_dp)<tol,'ATE gradient treatment row')
call check(abs(grad(3,3)+1.0_dp)<tol,'ATE balance gradient')
call check(abs(grad(4,3)+2.5_dp)<tol .and. abs(grad(5,3)+0.575_dp)<tol,'ATE covariate gradients')

f=[1.2_dp,-0.5_dp,0.8_dp]
qt=reshape([2.0_dp,-0.3_dp,0.8_dp, 0.4_dp,1.5_dp,-0.2_dp],[3,2])
vff=reshape([2.0_dp,0.2_dp,0.1_dp, 0.2_dp,1.5_dp,0.05_dp, 0.1_dp,0.05_dp,1.2_dp],[3,3])
vtf=reshape([ &
 0.1_dp,-0.02_dp,0.05_dp,-0.03_dp,0.02_dp,0.01_dp, &
 0.0_dp,0.1_dp,0.02_dp,0.04_dp,-0.01_dp,0.03_dp, &
 0.05_dp,0.0_dp,0.1_dp,0.02_dp,0.05_dp,-0.02_dp],[6,3])
call kleibergen_k_from_blocks(f,qt,vff,vtf,100,2,kr)
call check(abs(kr%k_stat-0.01407811443190007_dp)<2e-12_dp,'K statistic')
call check(abs(kr%s_stat-0.01459394453876627_dp)<2e-12_dp,'S statistic')
call check(abs(kr%j_stat-0.00051583010686620_dp)<2e-12_dp,'J statistic')


do i=1,80
   u(i,1)=sin(0.17_dp*real(i,dp))+0.4_dp*sin(0.17_dp*real(i-1,dp))
   u(i,2)=cos(0.11_dp*real(i,dp))+0.2_dp*cos(0.11_dp*real(i-1,dp))
end do
gg(:,1)=[-1.0_dp,-0.5_dp]
b1=bw_andrews(u,'Quadratic Spectral','AR(1)')
b2=bw_wilhelm(u,gg,'Quadratic Spectral','AR(1)')
call check(b1>=0.0_dp .and. b1<huge(1.0_dp),'Andrews bandwidth finite')
call check(b2>=0.0_dp .and. b2<huge(1.0_dp),'Wilhelm bandwidth finite')

print '(a)','test_extra: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok
character(len=*),intent(in)::msg
if(.not.ok)then
  print '(a)', 'FAIL: '//trim(msg)
  error stop 1
end if
end subroutine
end program
