module spatialextremes_covariance
   use spatialextremes_base, only: dp,pi,pair_count,pair_indices,distance_to_origin,nan_dp
   use r_compat, only: besselJ,besselK
   implicit none
   private
   integer, parameter, public :: COV_MATERN=1,COV_CAUCHY=2,COV_POWEREXP=3,COV_BESSEL=4,COV_CAUGEN=5
   public :: whittle_matern,cauchy_cov,generalized_cauchy,powered_exponential,bessel_cov
   public :: covariance_values,covariance_matrix,mahalanobis_distances_2d,mahalanobis_distances_3d
   public :: brown_resnick_a,geom_gauss_a,fbm_covariance
contains
   pure real(dp) function whittle_matern(h,nugget,sill,range,smooth) result(v)
      real(dp),intent(in)::h,nugget,sill,range,smooth
      real(dp)::x,c
      if(range<=0.0_dp .or. sill<=0.0_dp .or. nugget<0.0_dp .or. smooth<=epsilon(1.0_dp)) then
         v=nan_dp()
         return
      end if
      if(h==0.0_dp) then
      v=sill+nugget
      return
      end if
      x=h/range
      c=sill*2.0_dp**(1.0_dp-smooth)/gamma(smooth)
      v=c*x**smooth*besselK(x,smooth,.false.)
   end function whittle_matern

   pure real(dp) function cauchy_cov(h,nugget,sill,range,smooth) result(v)
      real(dp),intent(in)::h,nugget,sill,range,smooth
      if(h==0.0_dp) then
      v=sill+nugget
      else
      v=sill*(1.0_dp+(h/range)**2)**(-smooth)
      end if
   end function cauchy_cov

   pure real(dp) function generalized_cauchy(h,nugget,sill,range,smooth,smooth2) result(v)
      real(dp),intent(in)::h,nugget,sill,range,smooth,smooth2
      if(h==0.0_dp) then
      v=sill+nugget
      else
      v=sill*(1.0_dp+(h/range)**smooth2)**(-smooth/smooth2)
      end if
   end function generalized_cauchy

   pure real(dp) function powered_exponential(h,nugget,sill,range,smooth) result(v)
      real(dp),intent(in)::h,nugget,sill,range,smooth
      if(h==0.0_dp) then
      v=sill+nugget
      else
      v=sill*exp(-(h/range)**smooth)
      end if
   end function powered_exponential

   pure real(dp) function bessel_cov(h,dim,nugget,sill,range,smooth) result(v)
      real(dp),intent(in)::h,nugget,sill,range,smooth
      integer,intent(in)::dim
      real(dp)::x,c
      if(h==0.0_dp) then
      v=sill+nugget
      return
      end if
      x=h/range
      c=sill*2.0_dp**smooth*gamma(smooth+1.0_dp)
      if(x<=1.0e5_dp) then
      v=c*x**(-smooth)*besselJ(x,smooth)
      else
      v=c*x**(-smooth)*sqrt(2.0_dp/pi)*cos(x-smooth*pi/2.0_dp-pi/4.0_dp)
      end if
   end function bessel_cov

   pure real(dp) function covariance_values(h,model,nugget,sill,range,smooth,smooth2,dim) result(v)
      real(dp),intent(in)::h,nugget,sill,range,smooth
      real(dp),intent(in),optional::smooth2
      integer,intent(in)::model
      integer,intent(in),optional::dim
      real(dp)::s2
      integer::d
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      d=2
      if(present(dim))d=dim
      select case(model)
      case(COV_MATERN);v=whittle_matern(h,nugget,sill,range,smooth)
      case(COV_CAUCHY);v=cauchy_cov(h,nugget,sill,range,smooth)
      case(COV_POWEREXP);v=powered_exponential(h,nugget,sill,range,smooth)
      case(COV_BESSEL);v=bessel_cov(h,d,nugget,sill,range,smooth)
      case(COV_CAUGEN);v=generalized_cauchy(h,nugget,sill,range,smooth,s2)
      case default;v=nan_dp()
      end select
   end function covariance_values

   function covariance_matrix(coord,model,nugget,sill,range,smooth,smooth2) result(cov)
      real(dp),intent(in)::coord(:,:),nugget,sill,range,smooth
      integer,intent(in)::model
      real(dp),intent(in),optional::smooth2
      real(dp)::cov(size(coord,1),size(coord,1)),h,s2
      integer::i,j
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      do i=1,size(coord,1)
         do j=i,size(coord,1)
            h=sqrt(sum((coord(i,:)-coord(j,:))**2))
            cov(i,j)=covariance_values(h,model,nugget,sill,range,smooth,s2,size(coord,2))
            cov(j,i)=cov(i,j)
         end do
      end do
   end function covariance_matrix

   function mahalanobis_distances_2d(dv,cov) result(m)
      real(dp),intent(in)::dv(:,:),cov(2,2)
      real(dp)::m(size(dv,1)),det
      integer::i
      det=cov(1,1)*cov(2,2)-cov(1,2)**2
      do i=1,size(dv,1)
         m(i)=sqrt((cov(1,1)*dv(i,2)**2-2.0_dp*cov(1,2)*dv(i,1)*dv(i,2)+cov(2,2)*dv(i,1)**2)/det)
      end do
   end function mahalanobis_distances_2d

   function mahalanobis_distances_3d(dv,cov) result(m)
      use spatialextremes_base, only: inverse_spd
      real(dp),intent(in)::dv(:,:),cov(3,3)
      real(dp)::m(size(dv,1)),ic(3,3),x(3)
      integer::i,info
      call inverse_spd(cov,ic,info)
      if(info/=0) then
      m=nan_dp()
      return
      end if
      do i=1,size(dv,1)
      x=dv(i,:)
      m(i)=sqrt(dot_product(x,matmul(ic,x)))
      end do
   end function mahalanobis_distances_3d

   pure real(dp) function brown_resnick_a(h,range,smooth) result(a)
      real(dp),intent(in)::h,range,smooth
      a=sqrt(2.0_dp)*(h/range)**(0.5_dp*smooth)
   end function brown_resnick_a

   pure real(dp) function geom_gauss_a(h,model,sigma2,nugget,range,smooth,smooth2,dim) result(a)
      real(dp),intent(in)::h,sigma2,nugget,range,smooth
      integer,intent(in)::model
      real(dp),intent(in),optional::smooth2
      integer,intent(in),optional::dim
      real(dp)::rho,s2
      integer::d
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      d=2
      if(present(dim))d=dim
      rho=covariance_values(h,model,nugget,1.0_dp-nugget,range,smooth,s2,d)
      a=sqrt(2.0_dp*sigma2*(1.0_dp-rho))
   end function geom_gauss_a

   function fbm_covariance(coord,range,smooth,sill) result(cov)
      real(dp),intent(in)::coord(:,:),range,smooth,sill
      real(dp)::cov(size(coord,1),size(coord,1)),r(size(coord,1)),h
      integer::i,j
      r=distance_to_origin(coord)
      r=(r/range)**smooth
      do i=1,size(coord,1)
         cov(i,i)=sill*2.0_dp*r(i)
         do j=i+1,size(coord,1)
            h=sqrt(sum((coord(i,:)-coord(j,:))**2))
            cov(i,j)=sill*(r(i)+r(j)-(h/range)**smooth)
            cov(j,i)=cov(i,j)
         end do
      end do
   end function fbm_covariance
end module spatialextremes_covariance
