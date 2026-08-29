module spatialextremes_spatgev
   use spatialextremes_base, only: dp,is_finite,neg_huge
   use spatialextremes_univariate, only: dgev
   ! Penalized spline presentation is kept outside this numerical kernel.
   implicit none
   private
   public :: design_to_gev,spatgev_loglik,spatgev_loglik_trend
contains
   subroutine design_to_gev(xloc,xscale,xshape,bloc,bscale,bshape,loc,scale,shape,log_scale)
      real(dp),intent(in)::xloc(:,:),xscale(:,:),xshape(:,:),bloc(:),bscale(:),bshape(:)
      real(dp),intent(out)::loc(size(xloc,1)),scale(size(xscale,1)),shape(size(xshape,1))
      logical,intent(in),optional::log_scale
      logical::lg
      lg=.false.
      if(present(log_scale))lg=log_scale
      loc=matmul(xloc,bloc)
      scale=matmul(xscale,bscale)
      shape=matmul(xshape,bshape)
      if(lg)scale=exp(scale)
   end subroutine

   real(dp) function spatgev_loglik(data,loc,scale,shape) result(ll)
      real(dp),intent(in)::data(:,:),loc(:),scale(:),shape(:)
      integer::i,j
      ll=0.0_dp
      if(any(scale<=0.0_dp))then
      ll=neg_huge
      return
      end if
      do j=1,size(data,2)
         do i=1,size(data,1)
            if(is_finite(data(i,j)))ll=ll+dgev(data(i,j),loc(j),scale(j),shape(j),.true.)
         end do
      end do
   end function

   real(dp) function spatgev_loglik_trend(data,loc,scale,shape,tloc,tscale,tshape) result(ll)
      real(dp),intent(in)::data(:,:),loc(:),scale(:),shape(:),tloc(:),tscale(:),tshape(:)
      integer::i,j
      real(dp)::sc
      ll=0.0_dp
      do j=1,size(data,2)
         do i=1,size(data,1)
            sc=scale(j)+tscale(i)
            if(sc<=0.0_dp)then
            ll=neg_huge
            return
            end if
            if(is_finite(data(i,j)))ll=ll+dgev(data(i,j),loc(j)+tloc(i),sc,shape(j)+tshape(i),.true.)
         end do
      end do
   end function
end module spatialextremes_spatgev
