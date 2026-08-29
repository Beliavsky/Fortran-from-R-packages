program test_design_inference
   use SpatialExtremes
   implicit none
   integer,parameter::nobs=80,nsite=3
   real(dp)::data(nobs,nsite),coord(nsite,2),xloc(nsite,1),xscale(nsite,1),xshape(nsite,1)
   real(dp)::txloc(nobs,1),txscale(nobs,1),txshape(nobs,1)
   real(dp)::bloc(1),bscale(1),bshape(1),btl(1),bts(1),bth(1),cov(2,2)
   real(dp)::z(nobs,nsite),f(nobs,nsite),jac(nobs,nsite),u,loc,sc,sh
   real(dp),allocatable::theta(:)
   logical,allocatable::active(:)
   type(composite_se_t)::se
   integer::i,j,info

   coord=reshape([0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.3_dp,0.9_dp],shape(coord),order=[2,1])
   xloc=1.0_dp
   xscale=1.0_dp
   xshape=1.0_dp
   bloc=[0.0_dp]
   bscale=[1.5_dp]
   bshape=[0.1_dp]
   btl=[0.10_dp]
   bts=[0.04_dp]
   bth=[0.01_dp]
   do i=1,nobs
      txloc(i,1)=sin(0.17_dp*real(i,dp))
      txscale(i,1)=cos(0.11_dp*real(i,dp))
      txshape(i,1)=sin(0.07_dp*real(i,dp))
      do j=1,nsite
         u=0.08_dp+0.84_dp*real(mod(37*i+19*j,101),dp)/100.0_dp
         z(i,j)=1.0_dp/(-log(u))
         loc=bloc(1)+btl(1)*txloc(i,1)
         sc=bscale(1)+bts(1)*txscale(i,1)
         sh=bshape(1)+bth(1)*txshape(i,1)
         data(i,j)=frechet_to_gev(z(i,j),loc,sc,sh)
      end do
   end do

   theta=[1.2_dp,0.8_dp,bloc,bscale,bshape,btl,bts,bth]
   call gev_design_frechet(theta,2,data,xloc,xscale,xshape,txloc,txscale,txshape,f,jac,info)
   call assert_true(info==0,'design Frechet transform status')
   call assert_true(maxval(abs(f-z))<1.0e-10_dp,'design Frechet transform inversion')

   cov=0.0_dp
   cov(1,1)=1.0_dp
   cov(2,2)=1.0_dp
   allocate(active(7))
   active=.false.
   active(1)=.true.
   se=smith_design_standard_errors(data,coord,cov,xloc,xscale,xshape,bloc,bscale,bshape, &
      txloc,txscale,txshape,btl,bts,bth,active=active,isotropic=.true.)
   call assert_true(se%info==0.and.size(se%stderr)==1,'Smith design standard errors')
   deallocate(active)

   allocate(active(9))
   active=.false.
   active(2)=.true.
   se=schlather_design_standard_errors(data,coord,COV_MATERN,0.05_dp,1.2_dp,1.0_dp, &
      xloc,xscale,xshape,bloc,bscale,bshape,temp_xloc=txloc,temp_xscale=txscale, &
      temp_xshape=txshape,btemp_loc=btl,btemp_scale=bts,btemp_shape=bth,active=active)
   call assert_true(se%info==0.and.size(se%stderr)==1,'Schlather design standard errors')
   deallocate(active)

   allocate(active(4))
   active=.false.
   active(4)=.true.
   se=schlather_frechet_standard_errors(z,coord,COV_CAUGEN,0.05_dp,1.2_dp,1.0_dp,1.4_dp,active=active)
   call assert_true(se%info==0.and.size(se%stderr)==1,'generalized-Cauchy smooth2 standard errors')
   deallocate(active)

   allocate(active(10))
   active=.false.
   active(3)=.true.
   se=schlather_ind_design_standard_errors(data,coord,COV_MATERN,0.15_dp,0.05_dp,1.2_dp,1.0_dp, &
      xloc,xscale,xshape,bloc,bscale,bshape,temp_xloc=txloc,temp_xscale=txscale, &
      temp_xshape=txshape,btemp_loc=btl,btemp_scale=bts,btemp_shape=bth,active=active)
   call assert_true(se%info==0.and.size(se%stderr)==1,'Schlather-ind design standard errors')
   deallocate(active)

   allocate(active(8))
   active=.false.
   active(1)=.true.
   se=brownresnick_design_standard_errors(data,coord,1.2_dp,1.0_dp,xloc,xscale,xshape, &
      bloc,bscale,bshape,txloc,txscale,txshape,btl,bts,bth,active=active)
   call assert_true(se%info==0.and.size(se%stderr)==1,'Brown-Resnick design standard errors')
   deallocate(active)

   allocate(active(10))
   active=.false.
   active(3)=.true.
   se=geomgauss_design_standard_errors(data,coord,COV_MATERN,1.0_dp,0.05_dp,1.2_dp,1.0_dp, &
      xloc,xscale,xshape,bloc,bscale,bshape,temp_xloc=txloc,temp_xscale=txscale, &
      temp_xshape=txshape,btemp_loc=btl,btemp_scale=bts,btemp_shape=bth,active=active)
   call assert_true(se%info==0.and.size(se%stderr)==1,'geometric-Gaussian design standard errors')
   deallocate(active)

   allocate(active(10))
   active=.false.
   active(2)=.true.
   se=extremalt_design_standard_errors(data,coord,COV_MATERN,0.05_dp,1.2_dp,1.0_dp,3.0_dp, &
      xloc,xscale,xshape,bloc,bscale,bshape,temp_xloc=txloc,temp_xscale=txscale, &
      temp_xshape=txshape,btemp_loc=btl,btemp_scale=bts,btemp_shape=bth,active=active)
   call assert_true(se%info==0.and.size(se%stderr)==1,'extremal-t design standard errors')
   deallocate(active)

   allocate(active(6))
   active=.false.
   active(4)=.true.
   se=spatgev_design_standard_errors(data,xloc,xscale,xshape,bloc,bscale,bshape, &
      txloc,txscale,txshape,btl,bts,bth,active=active)
   call assert_true(se%info==0.and.size(se%stderr)==1,'spatial-GEV temporal design standard errors')

   print '(a)','design inference parity tests passed'
contains
   subroutine assert_true(ok,msg)
      logical,intent(in)::ok
      character(len=*),intent(in)::msg
      if(.not.ok)then
         write(*,'(a)')'FAIL: '//trim(msg)
         error stop 1
      end if
   end subroutine assert_true
end program test_design_inference
