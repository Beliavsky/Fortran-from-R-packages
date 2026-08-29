module spatialextremes_univariate
   use spatialextremes_base, only: dp,neg_huge,is_finite,nan_dp,exp_rand
   use r_compat, only: runif1
   implicit none
   private
   public :: dgev,pgev,qgev,rgev,dgpd,pgpd,qgpd,rgpd,gev_loglik,gpd_loglik
   public :: gev_to_frechet, gev_to_frechet_trend, frechet_to_gev, gev_to_uniform
contains
   pure real(dp) function pgev(x,loc,scale,shape,lower_tail) result(p)
      real(dp),intent(in)::x,loc,scale,shape
      logical,intent(in),optional::lower_tail
      real(dp)::z,t
      logical::lt
      lt=.true.
      if(present(lower_tail))lt=lower_tail
      if(scale<=0.0_dp) then
      p=nan_dp()
      return
      end if
      z=(x-loc)/scale
      if(abs(shape)<=epsilon(1.0_dp)) then
         p=exp(-exp(-z))
      else
         t=1.0_dp+shape*z
         if(t<=0.0_dp) then
            if(shape>0.0_dp) then
            p=0.0_dp
            else
            p=1.0_dp
            end if
         else
            p=exp(-t**(-1.0_dp/shape))
         end if
      end if
      if(.not.lt)p=1.0_dp-p
   end function pgev

   pure real(dp) function dgev(x,loc,scale,shape,log_density) result(d)
      real(dp),intent(in)::x,loc,scale,shape
      logical,intent(in),optional::log_density
      real(dp)::z,t,ld
      logical::lg
      lg=.false.
      if(present(log_density))lg=log_density
      if(scale<=0.0_dp) then
      d=nan_dp()
      return
      end if
      z=(x-loc)/scale
      if(abs(shape)<=epsilon(1.0_dp)) then
         ld=-log(scale)-z-exp(-z)
      else
         t=1.0_dp+shape*z
         if(t<=0.0_dp) then
         ld=-huge(1.0_dp)
         else
         ld=-log(scale)-t**(-1.0_dp/shape)-(1.0_dp/shape+1.0_dp)*log(t)
         end if
      end if
      if(lg) then
      d=ld
      else
      d=exp(ld)
      end if
   end function dgev

   pure real(dp) function qgev(p,loc,scale,shape,lower_tail) result(x)
      real(dp),intent(in)::p,loc,scale,shape
      logical,intent(in),optional::lower_tail
      real(dp)::pp
      logical::lt
      lt=.true.
      if(present(lower_tail))lt=lower_tail
      pp=p
      if(.not.lt)pp=1.0_dp-p
      if(pp<=0.0_dp .or. pp>=1.0_dp .or. scale<0.0_dp) then
      x=nan_dp()
      return
      end if
      if(abs(shape)<=epsilon(1.0_dp)) then
         x=loc-scale*log(-log(pp))
      else
         x=loc+scale*((-log(pp))**(-shape)-1.0_dp)/shape
      end if
   end function qgev

   real(dp) function rgev(loc,scale,shape) result(x)
      real(dp),intent(in)::loc,scale,shape
      real(dp)::e
      e=exp_rand()
      if(abs(shape)<=epsilon(1.0_dp)) then
      x=loc-scale*log(e)
      else
      x=loc+scale*(e**(-shape)-1.0_dp)/shape
      end if
   end function rgev

   pure real(dp) function pgpd(x,loc,scale,shape,lower_tail,lambda) result(p)
      real(dp),intent(in)::x,loc,scale,shape
      logical,intent(in),optional::lower_tail
      real(dp),intent(in),optional::lambda
      real(dp)::z,lam,t
      logical::lt
      lt=.true.
      if(present(lower_tail))lt=lower_tail
      lam=0.0_dp
      if(present(lambda))lam=lambda
      z=max(x-loc,0.0_dp)/scale
      if(abs(shape)<=epsilon(1.0_dp)) then
      p=1.0_dp-(1.0_dp-lam)*exp(-z)
      else
         t=max(1.0_dp+shape*z,0.0_dp)
         if(t==0.0_dp) then
         p=1.0_dp
         else
         p=1.0_dp-(1.0_dp-lam)*t**(-1.0_dp/shape)
         end if
      end if
      if(.not.lt)p=1.0_dp-p
   end function pgpd

   pure real(dp) function dgpd(x,loc,scale,shape,log_density) result(d)
      real(dp),intent(in)::x,loc,scale,shape
      logical,intent(in),optional::log_density
      real(dp)::z,t,ld
      logical::lg
      lg=.false.
      if(present(log_density))lg=log_density
      z=(x-loc)/scale
      if(scale<=0.0_dp .or. z<=0.0_dp) then
      ld=-huge(1.0_dp)
      else if(abs(shape)<=epsilon(1.0_dp)) then
      ld=-log(scale)-z
      else
         t=1.0_dp+shape*z
         if(t<=0.0_dp) then
         ld=-huge(1.0_dp)
         else
         ld=-log(scale)-(1.0_dp/shape+1.0_dp)*log(t)
         end if
      end if
      if(lg) then
      d=ld
      else
      d=exp(ld)
      end if
   end function dgpd

   pure real(dp) function qgpd(p,loc,scale,shape,lower_tail,lambda) result(x)
      real(dp),intent(in)::p,loc,scale,shape
      logical,intent(in),optional::lower_tail
      real(dp),intent(in),optional::lambda
      real(dp)::pp,lam
      logical::lt
      lt=.true.
      if(present(lower_tail))lt=lower_tail
      lam=0.0_dp
      if(present(lambda))lam=lambda
      pp=p
      if(lt)pp=1.0_dp-p
      pp=pp/(1.0_dp-lam)
      if(abs(shape)<=epsilon(1.0_dp)) then
      x=loc-scale*log(pp)
      else
      x=loc+scale*(pp**(-shape)-1.0_dp)/shape
      end if
   end function qgpd

   real(dp) function rgpd(loc,scale,shape) result(x)
      real(dp),intent(in)::loc,scale,shape
      real(dp)::u
      u=runif1()
      if(abs(shape)<=epsilon(1.0_dp)) then
      x=loc-scale*log(u)
      else
      x=loc+scale*(u**(-shape)-1.0_dp)/shape
      end if
   end function rgpd

   real(dp) function gev_loglik(data,loc,scale,shape) result(ll)
      real(dp),intent(in)::data(:),loc,scale,shape
      integer::i
      ll=0.0_dp
      if(scale<=0.0_dp .or. shape< -1.0_dp) then
      ll=-1.0e6_dp
      return
      end if
      do i=1,size(data)
         if(is_finite(data(i))) ll=ll+dgev(data(i),loc,scale,shape,.true.)
      end do
   end function gev_loglik

   real(dp) function gpd_loglik(data,threshold,scale,shape) result(ll)
      real(dp),intent(in)::data(:),threshold,scale,shape
      integer::i
      ll=0.0_dp
      if(scale<=0.0_dp .or. shape< -1.0_dp) then
      ll=-1.0e6_dp
      return
      end if
      do i=1,size(data)
         if(data(i)<=threshold) then
         ll=-1.0e6_dp
         return
         end if
         ll=ll+dgpd(data(i),threshold,scale,shape,.true.)
      end do
   end function gpd_loglik

   subroutine gev_to_frechet(data,loc,scale,shape,frech,logjac,info)
      real(dp),intent(in)::data(:,:),loc(:),scale(:),shape(:)
      real(dp),intent(out)::frech(size(data,1),size(data,2)),logjac(size(data,1),size(data,2))
      integer,intent(out)::info
      integer::i,j
      real(dp)::t,z
      info=0
      do j=1,size(data,2)
         if(scale(j)<=0.0_dp) then
         info=1
         return
         end if
         do i=1,size(data,1)
            if(.not.is_finite(data(i,j))) then
            frech(i,j)=data(i,j)
            logjac(i,j)=data(i,j)
            cycle
            end if
            z=(data(i,j)-loc(j))/scale(j)
            if(abs(shape(j))<=epsilon(1.0_dp)) then
               frech(i,j)=exp(z)
               logjac(i,j)=z-log(scale(j))
            else
               t=1.0_dp+shape(j)*z
               if(t<=0.0_dp) then
               info=2
               return
               end if
               logjac(i,j)=(1.0_dp/shape(j)-1.0_dp)*log(t)-log(scale(j))
               frech(i,j)=t**(1.0_dp/shape(j))
            end if
         end do
      end do
   end subroutine gev_to_frechet

   subroutine gev_to_frechet_trend(data,loc,scale,shape,tloc,tscale,tshape,frech,logjac,info)
      ! SpatialExtremes gev2frechTrend numerical kernel.  Spatial GEV
      ! parameters are augmented by observation-specific temporal trends.
      real(dp),intent(in)::data(:,:),loc(:),scale(:),shape(:),tloc(:),tscale(:),tshape(:)
      real(dp),intent(out)::frech(size(data,1),size(data,2)),logjac(size(data,1),size(data,2))
      integer,intent(out)::info
      integer::i,j
      real(dp)::lc,sc,sh,t,z
      info=0
      if(size(tloc)/=size(data,1).or.size(tscale)/=size(data,1).or.size(tshape)/=size(data,1))then
         info=3
         return
      end if
      if(size(loc)/=size(data,2).or.size(scale)/=size(data,2).or.size(shape)/=size(data,2))then
         info=4
         return
      end if
      do j=1,size(data,2)
         do i=1,size(data,1)
            if(.not.is_finite(data(i,j)))then
               frech(i,j)=data(i,j)
               logjac(i,j)=data(i,j)
               cycle
            end if
            lc=loc(j)+tloc(i)
            sc=scale(j)+tscale(i)
            sh=shape(j)+tshape(i)
            if(sc<=0.0_dp)then
            info=1
            return
            end if
            z=(data(i,j)-lc)/sc
            if(abs(sh)<=epsilon(1.0_dp))then
               frech(i,j)=exp(z)
               logjac(i,j)=z-log(sc)
            else
               t=1.0_dp+sh*z
               if(t<=0.0_dp)then
               info=2
               return
               end if
               logjac(i,j)=(1.0_dp/sh-1.0_dp)*log(t)-log(sc)
               frech(i,j)=t**(1.0_dp/sh)
            end if
         end do
      end do
   end subroutine gev_to_frechet_trend

   pure real(dp) function frechet_to_gev(z,loc,scale,shape) result(x)
      real(dp),intent(in)::z,loc,scale,shape
      if(abs(shape)<=epsilon(1.0_dp)) then
      x=loc+scale*log(z)
      else
      x=loc+scale*(z**shape-1.0_dp)/shape
      end if
   end function frechet_to_gev

   pure real(dp) function gev_to_uniform(x,loc,scale,shape) result(u)
      real(dp),intent(in)::x,loc,scale,shape
      u=pgev(x,loc,scale,shape)
   end function gev_to_uniform
end module spatialextremes_univariate
