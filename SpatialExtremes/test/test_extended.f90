program test_extended
   use SpatialExtremes
   use r_compat, only: set_seed_int
   implicit none
   real(dp)::a1(1,1),a2(1,1),c2(2,2),coord2(1,2),grid2(1,2)
   real(dp)::dc(1,1),tc(1,1),obs(1),xs(12000,1),rho,sm,sv
   real(dp)::allx(6),exc(4)
   real(dp)::ef(4,2),th(1),coord4(4,1),dist6(6),theta6(6)
   integer::i
   type(univ_fit_t)::fg1,fg2
   type(maxstab_fit_t)::fls

   a1=maxlinear_design([0.0_dp],[0.0_dp],4.0_dp,0.5_dp)
   call check(abs(a1(1,1)-0.5_dp/sqrt(8.0_dp*pi))<1e-14_dp,'1-D Smith max-linear design')
   coord2=0.0_dp
   grid2=0.0_dp
   c2=0.0_dp
   c2(1,1)=4.0_dp
   c2(2,2)=9.0_dp
   a2=maxlinear_design(coord2,grid2,c2,0.25_dp)
   call check(abs(a2(1,1)-0.25_dp/(12.0_dp*pi))<1e-14_dp,'2-D Smith max-linear design')

   dc(1,1)=0.0_dp
   tc(1,1)=1.0_dp
   obs=2.0_dp
   call set_seed_int(2468)
   xs=conditional_gaussian_process(12000,tc,dc,obs,COV_POWEREXP,1.0_dp,1.0_dp,1.0_dp)
   rho=exp(-1.0_dp)
   sm=sum(xs(:,1))/real(size(xs,1),dp)
   sv=sum((xs(:,1)-sm)**2)/real(size(xs,1)-1,dp)
   call check(abs(sm-2.0_dp*rho)<0.025_dp,'conditional GP mean')
   call check(abs(sv-(1.0_dp-rho*rho))<0.035_dp,'conditional GP variance')

   allx=[0.1_dp,0.5_dp,1.2_dp,1.4_dp,1.8_dp,2.6_dp]
   exc=allx(3:6)
   fg1=gpdmle(allx,1.0_dp)
   fg2=gpdmle(exc,1.0_dp)
   call check(abs(fg1%scale-fg2%scale)<1e-10_dp .and. abs(fg1%shape-fg2%shape)<1e-10_dp,'GPD threshold filtering')

   ef(:,1)=[0.5_dp,1.5_dp,0.5_dp,1.5_dp]
   ef(:,2)=ef(:,1)
   th=smith_extcoeff(ef)
   call check(abs(th(1)-1.0_dp)<1e-14_dp,'Smith extremal coefficient native kernel')
   ef(:,1)=[1.0_dp,2.0_dp,3.0_dp,4.0_dp]
   ef(:,2)=ef(:,1)
   th=st_extcoeff_frechet(ef)
   call check(abs(th(1)-1.0_dp)<1e-14_dp,'Schlather-Tawn complete dependence')

   coord4(:,1)=[0.0_dp,1.0_dp,2.0_dp,4.0_dp]
   dist6=euclidean_distances(coord4)
   do i=1,size(dist6)
      theta6(i)=extremal_coefficient_smith(brown_resnick_a(dist6(i),2.5_dp,1.2_dp))
   end do
   fls=lsfit_brownresnick(theta6,coord4,start=[2.0_dp,1.0_dp])
   call check(fls%convergence==0 .or. fls%iterations>0,'Brown-Resnick LS ran')
   call check(abs(fls%par(1)-2.5_dp)<2e-3_dp .and. abs(fls%par(2)-1.2_dp)<2e-3_dp,'Brown-Resnick LS recovery')
   call check(abs(tic_value(10.0_dp,reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2]), &
                       reshape([2.0_dp,0.0_dp,0.0_dp,3.0_dp],[2,2]))-20.0_dp)<1e-14_dp,'TIC kernel')
   print *,'test_extended: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
      print *,'FAIL: ',msg
      error stop 1
      end if
   end subroutine
end program test_extended
