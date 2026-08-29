program test_remaining_gaps
   use SpatialExtremes
   implicit none
   integer,parameter::nsite=3,nobs=8
   real(dp)::data(nobs,nsite),coord(nsite,2),dl(nsite,1),ds(nsite,1),dh(nsite,1)
   real(dp)::bl(1),bs(1),bh(1),sills(3),ranges(3),smooths(3),gev(nsite,3)
   real(dp)::hsill(2,3),hrange(2,3),hsmooth(2,3),bm(1),bp(1,1),pg(3),pr(3),psm(3)
   real(dp),allocatable::g(:,:)
   type(latent_mcmc_result_t)::fit
   type(conditional_maxstable_result_t)::cfit
   real(dp)::ccov(2,2),cy(1)
   integer::parts(1,1)
   integer::i,j,info

   coord=reshape([0.0_dp,0.0_dp, 1.0_dp,0.0_dp, 0.0_dp,1.0_dp],[nsite,2],order=[2,1])
   dl=1.0_dp
   ds=1.0_dp
   dh=1.0_dp
   bl=[0.0_dp]
   bs=[0.0_dp]
   bh=[0.0_dp]
   sills=0.5_dp
   ranges=2.0_dp
   smooths=1.2_dp
   gev(:,1)=0.0_dp
   gev(:,2)=1.0_dp
   gev(:,3)=0.05_dp
   do j=1,nsite
      do i=1,nobs
         data(i,j)=0.2_dp*real(i-4,dp)+0.05_dp*real(j-2,dp)
      end do
   end do
   hsill(1,:)=2.0_dp
   hsill(2,:)=0.5_dp
   hrange(1,:)=2.0_dp
   hrange(2,:)=2.0_dp
   hsmooth(1,:)=2.0_dp
   hsmooth(2,:)=1.0_dp
   bm=0.0_dp
   bp(1,1)=1.0_dp
   pg=[0.03_dp,0.03_dp,0.01_dp]
   pr=0.0_dp
   psm=0.0_dp
   call latent_gev_mcmc(data,coord,[COV_POWEREXP,COV_POWEREXP,COV_POWEREXP],dl,ds,dh,bl,bs,bh, &
      sills,ranges,smooths,gev,hsill,hrange,hsmooth,bm,bm,bm,bp,bp,bp,pg,pr,psm,5,fit,burn_in=2)
   if(fit%info/=0)error stop 'latent MCMC failed'
   if(any(shape(fit%chain_loc)/=[5,1+3+nsite]))error stop 'latent chain shape'
   if(any(fit%acceptance<0.0_dp) .or. any(fit%acceptance>1.0_dp))error stop 'latent acceptance'

   g=simulate_gaussian_grid_circulant(300,3,[1.0_dp,1.0_dp],COV_POWEREXP,0.1_dp,0.9_dp,2.0_dp,1.5_dp,info)
   if(info/=0)error stop 'circulant embedding failed'
   if(any(.not.(g==g)))error stop 'circulant NaN'
   if(abs(sum(g(:,1)**2)/300.0_dp-(sum(g(:,1))/300.0_dp)**2-1.0_dp)>0.35_dp) &
      error stop 'circulant variance regression'
   ccov=reshape([1.0_dp,0.5_dp,0.5_dp,1.0_dp],[2,2])
   cy=[1.0_dp]
   parts=1
   call conditional_schlather_given_partition(ccov,cy,parts,cfit)
   if(cfit%info/=0)error stop 'conditional Schlather failed'
   if(abs(cfit%sim(1,1)-cy(1))>1.0e-10_dp)error stop 'conditional value not preserved'
   print '(a)','remaining-gap tests passed' 
end program test_remaining_gaps
