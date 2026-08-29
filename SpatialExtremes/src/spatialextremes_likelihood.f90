module spatialextremes_likelihood
   use spatialextremes_base, only: dp,euclidean_distances,distance_vectors,neg_huge,pair_count
   use spatialextremes_univariate, only: gev_to_frechet
   use spatialextremes_covariance
   use spatialextremes_pairwise
   implicit none
   private
   public :: smith_loglik,schlather_loglik,schlather_ind_loglik,brown_resnick_loglik
   public :: geomgauss_loglik,extremalt_loglik
   public :: smith_loglik_contributions,schlather_loglik_contributions
   public :: schlather_ind_loglik_contributions,brown_resnick_loglik_contributions
   public :: geomgauss_loglik_contributions,extremalt_loglik_contributions
contains
   real(dp) function smith_loglik(data,coord,loc,scale,shape,cov,weights) result(ll)
      real(dp),intent(in)::data(:,:),coord(:,:),loc(:),scale(:),shape(:),cov(:,:)
      real(dp),intent(in),optional::weights(:)
      real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2))
      real(dp),allocatable::dv(:,:),a(:)
      integer::info
      call gev_to_frechet(data,loc,scale,shape,f,j,info)
      if(info/=0)then
      ll=neg_huge
      return
      end if
      dv=distance_vectors(coord)
      if(size(coord,2)==2)then
      a=mahalanobis_distances_2d(dv,cov(1:2,1:2))
      else if(size(coord,2)==3)then
      a=mahalanobis_distances_3d(dv,cov(1:3,1:3))
      else
      ll=neg_huge
      return
      end if
      if(any(a<=0.0_dp))then
      ll=neg_huge
      return
      end if
      ll=lplik_smith(f,a,j,weights)
   end function smith_loglik

   real(dp) function schlather_loglik(data,coord,loc,scale,shape,covmod,nugget,range,smooth,smooth2,weights) result(ll)
      real(dp),intent(in)::data(:,:),coord(:,:),loc(:),scale(:),shape(:),nugget,range,smooth
      integer,intent(in)::covmod
      real(dp),intent(in),optional::smooth2,weights(:)
      real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2)),s2
      real(dp),allocatable::d(:),rho(:)
      integer::i,info
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      if(nugget<0.0_dp.or.nugget>=1.0_dp.or.range<=0.0_dp)then
      ll=neg_huge
      return
      end if
      call gev_to_frechet(data,loc,scale,shape,f,j,info)
      if(info/=0)then
      ll=neg_huge
      return
      end if
      d=euclidean_distances(coord)
      allocate(rho(size(d)))
      do i=1,size(d)
      rho(i)=covariance_values(d(i),covmod,nugget,1.0_dp-nugget,range,smooth,s2,size(coord,2))
      end do
      ll=lplik_schlather(f,rho,j,weights)
   end function schlather_loglik

   real(dp) function schlather_ind_loglik(data,coord,loc,scale,shape,alpha,covmod,nugget,range,smooth,smooth2,weights) result(ll)
      real(dp),intent(in)::data(:,:),coord(:,:),loc(:),scale(:),shape(:),alpha,nugget,range,smooth
      integer,intent(in)::covmod
      real(dp),intent(in),optional::smooth2,weights(:)
      real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2)),s2
      real(dp),allocatable::d(:),rho(:)
      integer::i,info
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      if(alpha<0.0_dp.or.alpha>1.0_dp)then
      ll=neg_huge
      return
      end if
      call gev_to_frechet(data,loc,scale,shape,f,j,info)
      if(info/=0)then
      ll=neg_huge
      return
      end if
      d=euclidean_distances(coord)
      allocate(rho(size(d)))
      do i=1,size(d)
      rho(i)=covariance_values(d(i),covmod,nugget,1.0_dp-nugget,range,smooth,s2,size(coord,2))
      end do
      ll=lplik_schlather_ind(f,alpha,rho,j,weights)
   end function schlather_ind_loglik

   real(dp) function brown_resnick_loglik(data,coord,loc,scale,shape,range,smooth,weights) result(ll)
      real(dp),intent(in)::data(:,:),coord(:,:),loc(:),scale(:),shape(:),range,smooth
      real(dp),intent(in),optional::weights(:)
      real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2))
      real(dp),allocatable::d(:),a(:)
      integer::i,info
      if(range<=0.0_dp.or.smooth<=0.0_dp.or.smooth>2.0_dp)then
      ll=neg_huge
      return
      end if
      call gev_to_frechet(data,loc,scale,shape,f,j,info)
      if(info/=0)then
      ll=neg_huge
      return
      end if
      d=euclidean_distances(coord)
      allocate(a(size(d)))
      do i=1,size(d)
      a(i)=brown_resnick_a(d(i),range,smooth)
      end do
      ll=lplik_smith(f,a,j,weights)
   end function brown_resnick_loglik

   real(dp) function geomgauss_loglik(data,coord,loc,scale,shape,covmod,sigma2,nugget,range,smooth,smooth2,weights) result(ll)
      real(dp),intent(in)::data(:,:),coord(:,:),loc(:),scale(:),shape(:),sigma2,nugget,range,smooth
      integer,intent(in)::covmod
      real(dp),intent(in),optional::smooth2,weights(:)
      real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2)),s2
      real(dp),allocatable::d(:),a(:)
      integer::i,info
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      if(sigma2<=0.0_dp.or.nugget<0.0_dp.or.nugget>=1.0_dp)then
      ll=neg_huge
      return
      end if
      call gev_to_frechet(data,loc,scale,shape,f,j,info)
      if(info/=0)then
      ll=neg_huge
      return
      end if
      d=euclidean_distances(coord)
      allocate(a(size(d)))
      do i=1,size(d)
      a(i)=geom_gauss_a(d(i),covmod,sigma2,nugget,range,smooth,s2,size(coord,2))
      end do
      ll=lplik_smith(f,a,j,weights)
   end function geomgauss_loglik

   real(dp) function extremalt_loglik(data,coord,loc,scale,shape,nu,covmod,nugget,range,smooth,smooth2,weights) result(ll)
      real(dp),intent(in)::data(:,:),coord(:,:),loc(:),scale(:),shape(:),nu,nugget,range,smooth
      integer,intent(in)::covmod
      real(dp),intent(in),optional::smooth2,weights(:)
      real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2)),s2
      real(dp),allocatable::d(:),rho(:)
      integer::i,info
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      if(nu<=0.0_dp.or.nugget<0.0_dp.or.nugget>=1.0_dp)then
      ll=neg_huge
      return
      end if
      call gev_to_frechet(data,loc,scale,shape,f,j,info)
      if(info/=0)then
      ll=neg_huge
      return
      end if
      d=euclidean_distances(coord)
      allocate(rho(size(d)))
      do i=1,size(d)
      rho(i)=covariance_values(d(i),covmod,nugget,1.0_dp-nugget,range,smooth,s2,size(coord,2))
      end do
      ll=lplik_extremalt(f,rho,nu,j,weights)
   end function extremalt_loglik

   function smith_loglik_contributions(data,coord,loc,scale,shape,cov,weights) result(c)
      real(dp),intent(in)::data(:,:),coord(:,:),loc(:),scale(:),shape(:),cov(:,:)
      real(dp),intent(in),optional::weights(:)
      real(dp),allocatable::c(:,:),dv(:,:),a(:)
      real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2))
      integer::info
      allocate(c(size(data,1),pair_count(size(data,2))))
      call gev_to_frechet(data,loc,scale,shape,f,j,info)
      if(info/=0)then
      c=neg_huge
      return
      end if
      dv=distance_vectors(coord)
      if(size(coord,2)==2.and.size(cov,1)>=2.and.size(cov,2)>=2)then
         a=mahalanobis_distances_2d(dv,cov(1:2,1:2))
      else if(size(coord,2)==3.and.size(cov,1)>=3.and.size(cov,2)>=3)then
         a=mahalanobis_distances_3d(dv,cov(1:3,1:3))
      else
         c=neg_huge
         return
      end if
      if(any(a<=0.0_dp).or.any(a/=a))then
      c=neg_huge
      return
      end if
      if(present(weights))then
      c=lplik_smith_contributions(f,a,j,weights)
      else
      c=lplik_smith_contributions(f,a,j)
      end if
   end function smith_loglik_contributions

   function schlather_loglik_contributions(data,coord,loc,scale,shape,covmod,nugget,range,smooth,smooth2,weights) result(c)
      real(dp),intent(in)::data(:,:),coord(:,:),loc(:),scale(:),shape(:),nugget,range,smooth
      integer,intent(in)::covmod
      real(dp),intent(in),optional::smooth2,weights(:)
      real(dp),allocatable::c(:,:),d(:),rho(:)
      real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2)),s2
      integer::i,info
      allocate(c(size(data,1),pair_count(size(data,2))))
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      if(nugget<0.0_dp.or.nugget>=1.0_dp.or.range<=0.0_dp.or.smooth<=0.0_dp)then
      c=neg_huge
      return
      end if
      call gev_to_frechet(data,loc,scale,shape,f,j,info)
      if(info/=0)then
      c=neg_huge
      return
      end if
      d=euclidean_distances(coord)
      allocate(rho(size(d)))
      do i=1,size(d)
      rho(i)=covariance_values(d(i),covmod,nugget,1.0_dp-nugget,range,smooth,s2,size(coord,2))
      end do
      if(present(weights))then
      c=lplik_schlather_contributions(f,rho,j,weights)
      else
      c=lplik_schlather_contributions(f,rho,j)
      end if
   end function schlather_loglik_contributions

   function schlather_ind_loglik_contributions(data,coord,loc,scale,shape,alpha,covmod,nugget,range, &
      smooth,smooth2,weights) result(c)
      real(dp),intent(in)::data(:,:),coord(:,:),loc(:),scale(:),shape(:),alpha,nugget,range,smooth
      integer,intent(in)::covmod
      real(dp),intent(in),optional::smooth2,weights(:)
      real(dp),allocatable::c(:,:),d(:),rho(:)
      real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2)),s2
      integer::i,info
      allocate(c(size(data,1),pair_count(size(data,2))))
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      if(alpha<0.0_dp.or.alpha>1.0_dp.or.nugget<0.0_dp.or.nugget>=1.0_dp.or.range<=0.0_dp.or.smooth<=0.0_dp)then
         c=neg_huge
         return
      end if
      call gev_to_frechet(data,loc,scale,shape,f,j,info)
      if(info/=0)then
      c=neg_huge
      return
      end if
      d=euclidean_distances(coord)
      allocate(rho(size(d)))
      do i=1,size(d)
      rho(i)=covariance_values(d(i),covmod,nugget,1.0_dp-nugget,range,smooth,s2,size(coord,2))
      end do
      if(present(weights))then
      c=lplik_schlather_ind_contributions(f,alpha,rho,j,weights)
      else
      c=lplik_schlather_ind_contributions(f,alpha,rho,j)
      end if
   end function schlather_ind_loglik_contributions

   function brown_resnick_loglik_contributions(data,coord,loc,scale,shape,range,smooth,weights) result(c)
      real(dp),intent(in)::data(:,:),coord(:,:),loc(:),scale(:),shape(:),range,smooth
      real(dp),intent(in),optional::weights(:)
      real(dp),allocatable::c(:,:),d(:),a(:)
      real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2))
      integer::i,info
      allocate(c(size(data,1),pair_count(size(data,2))))
      if(range<=0.0_dp.or.smooth<=0.0_dp.or.smooth>2.0_dp)then
      c=neg_huge
      return
      end if
      call gev_to_frechet(data,loc,scale,shape,f,j,info)
      if(info/=0)then
      c=neg_huge
      return
      end if
      d=euclidean_distances(coord)
      allocate(a(size(d)))
      do i=1,size(d)
      a(i)=brown_resnick_a(d(i),range,smooth)
      end do
      if(present(weights))then
      c=lplik_smith_contributions(f,a,j,weights)
      else
      c=lplik_smith_contributions(f,a,j)
      end if
   end function brown_resnick_loglik_contributions

   function geomgauss_loglik_contributions(data,coord,loc,scale,shape,covmod,sigma2,nugget,range,smooth,smooth2,weights) result(c)
      real(dp),intent(in)::data(:,:),coord(:,:),loc(:),scale(:),shape(:),sigma2,nugget,range,smooth
      integer,intent(in)::covmod
      real(dp),intent(in),optional::smooth2,weights(:)
      real(dp),allocatable::c(:,:),d(:),a(:)
      real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2)),s2
      integer::i,info
      allocate(c(size(data,1),pair_count(size(data,2))))
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      if(sigma2<=0.0_dp.or.nugget<0.0_dp.or.nugget>=1.0_dp.or.range<=0.0_dp.or.smooth<=0.0_dp)then
      c=neg_huge
      return
      end if
      call gev_to_frechet(data,loc,scale,shape,f,j,info)
      if(info/=0)then
      c=neg_huge
      return
      end if
      d=euclidean_distances(coord)
      allocate(a(size(d)))
      do i=1,size(d)
      a(i)=geom_gauss_a(d(i),covmod,sigma2,nugget,range,smooth,s2,size(coord,2))
      end do
      if(present(weights))then
      c=lplik_smith_contributions(f,a,j,weights)
      else
      c=lplik_smith_contributions(f,a,j)
      end if
   end function geomgauss_loglik_contributions

   function extremalt_loglik_contributions(data,coord,loc,scale,shape,nu,covmod,nugget,range,smooth,smooth2,weights) result(c)
      real(dp),intent(in)::data(:,:),coord(:,:),loc(:),scale(:),shape(:),nu,nugget,range,smooth
      integer,intent(in)::covmod
      real(dp),intent(in),optional::smooth2,weights(:)
      real(dp),allocatable::c(:,:),d(:),rho(:)
      real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2)),s2
      integer::i,info
      allocate(c(size(data,1),pair_count(size(data,2))))
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      if(nu<=0.0_dp.or.nugget<0.0_dp.or.nugget>=1.0_dp.or.range<=0.0_dp.or.smooth<=0.0_dp)then
      c=neg_huge
      return
      end if
      call gev_to_frechet(data,loc,scale,shape,f,j,info)
      if(info/=0)then
      c=neg_huge
      return
      end if
      d=euclidean_distances(coord)
      allocate(rho(size(d)))
      do i=1,size(d)
      rho(i)=covariance_values(d(i),covmod,nugget,1.0_dp-nugget,range,smooth,s2,size(coord,2))
      end do
      if(present(weights))then
      c=lplik_extremalt_contributions(f,rho,nu,j,weights)
      else
      c=lplik_extremalt_contributions(f,rho,nu,j)
      end if
   end function extremalt_loglik_contributions

end module spatialextremes_likelihood
