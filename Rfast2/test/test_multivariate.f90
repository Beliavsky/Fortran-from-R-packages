program test_multivariate
   use rfast2
   implicit none
   type(pca_result) :: pc
   type(pcr_result) :: pr
   real(dp) :: x(6,2),y(6),d(2),lev(6)
   integer :: k(2)

   x(:,1)=[-3.0_dp,-2.0_dp,-1.0_dp,1.0_dp,2.0_dp,3.0_dp]
   x(:,2)=[1.0_dp,-1.0_dp,1.0_dp,-1.0_dp,1.0_dp,-1.0_dp]
   pc=pca(x,center=.true.,scale=.false.,k=2,vectors=.true.)
   if (pc%status /= 0 .or. pc%values(1) < pc%values(2)) error stop 1
   if (abs(sum(pc%values)-sum((x-spread(sum(x,dim=1)/6.0_dp,1,6))**2)/5.0_dp) > 1.0e-8_dp) error stop 2
   y=2.0_dp+3.0_dp*x(:,1)-0.5_dp*x(:,2)
   k=[1,2]
   pr=pcr(y,x,k,xnew=x)
   if (pr%status /= 0 .or. size(pr%beta,2) /= 2) error stop 3
   if (maxval(abs(pr%fitted(:,2)-y)) > 1.0e-7_dp) error stop 4
   d=discriminability(reshape([0.0_dp,0.0_dp,1.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp],[4,2]))
   if (any(abs(d) > 1.0_dp)) error stop 5
   lev=leverage(reshape([1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],[6,2]))
   if (any(lev < 0.0_dp)) error stop 6
   print '(a)', 'test_multivariate: PASS'
end program test_multivariate
