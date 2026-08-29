module spatialextremes_dependence
   use spatialextremes_base, only: dp,pair_count,pair_indices,is_finite
   implicit none
   private
   public :: madogram,variogram,fmadogram_extcoeff,smith_extcoeff,extremal_coeff_from_madogram
   public :: empirical_concurrence, empirical_boot_concurrence, kendall_concurrence, lambda_madogram, logit
contains
   function madogram(data) result(out)
      real(dp),intent(in)::data(:,:)
      real(dp)::out(pair_count(size(data,2)))
      integer::k,i,j,t,n
      out=0.0_dp
      do k=1,size(out)
         call pair_indices(k,size(data,2),i,j)
         n=0
         do t=1,size(data,1)
            if(is_finite(data(t,i)).and.is_finite(data(t,j))) then
               out(k)=out(k)+abs(data(t,i)-data(t,j))
               n=n+1
            end if
         end do
         if(n>0)out(k)=0.5_dp*out(k)/n
      end do
   end function madogram

   function variogram(data) result(out)
      real(dp),intent(in)::data(:,:)
      real(dp)::out(pair_count(size(data,2)))
      integer::k,i,j,t,n
      real(dp)::z
      out=0.0_dp
      do k=1,size(out)
         call pair_indices(k,size(data,2),i,j)
         n=0
         do t=1,size(data,1)
            if(is_finite(data(t,i)).and.is_finite(data(t,j))) then
               z=data(t,i)-data(t,j)
               out(k)=out(k)+z*z
               n=n+1
            end if
         end do
         if(n>0)out(k)=0.5_dp*out(k)/n
      end do
   end function variogram

   pure function fmadogram_extcoeff(fmado) result(theta)
      real(dp),intent(in)::fmado(:)
      real(dp)::theta(size(fmado))
      theta=(1.0_dp+2.0_dp*fmado)/(1.0_dp-2.0_dp*fmado)
   end function fmadogram_extcoeff

   function smith_extcoeff(frech) result(theta)
      real(dp),intent(in)::frech(:,:)
      real(dp)::theta(pair_count(size(frech,2)))
      integer::k,i,j,t
      theta=0.0_dp
      do k=1,size(theta)
         call pair_indices(k,size(frech,2),i,j)
         do t=1,size(frech,1)
         theta(k)=theta(k)+min(frech(t,i),frech(t,j))
         end do
         theta(k)=real(size(frech,1),dp)/theta(k)
      end do
   end function smith_extcoeff

   pure function extremal_coeff_from_madogram(mado,scale,shape) result(theta)
      real(dp),intent(in)::mado(:),scale,shape
      real(dp)::theta(size(mado))
      if(abs(shape)<=epsilon(1.0_dp)) then
         theta=exp(mado/scale)
      else
         theta=(1.0_dp+shape*mado/(scale*gamma(1.0_dp-shape)))**(1.0_dp/shape)
      end if
   end function extremal_coeff_from_madogram

   function empirical_concurrence(data,block_size) result(p)
      real(dp),intent(in)::data(:,:)
      integer,intent(in)::block_size
      real(dp)::p(pair_count(size(data,2)))
      integer::k,i,j,b,t0,t1,nb,mi,mj
      nb=size(data,1)/block_size
      p=0.0_dp
      do k=1,size(p)
         call pair_indices(k,size(data,2),i,j)
         do b=1,nb
            t0=(b-1)*block_size+1
            t1=b*block_size
            mi=maxloc(data(t0:t1,i),1)
            mj=maxloc(data(t0:t1,j),1)
            if(mi==mj)p(k)=p(k)+1.0_dp
         end do
         if(nb>0)p(k)=p(k)/nb
      end do
   end function empirical_concurrence

   function kendall_concurrence(data) result(p)
      real(dp),intent(in)::data(:,:)
      real(dp)::p(pair_count(size(data,2)))
      integer::k,i,j,a,b,n,conc,tot
      p=0.0_dp
      n=size(data,1)
      do k=1,size(p)
         call pair_indices(k,size(data,2),i,j)
         conc=0
         tot=0
         do a=1,n-1
            if(.not.(is_finite(data(a,i)).and.is_finite(data(a,j))))cycle
            do b=a+1,n
               if(.not.(is_finite(data(b,i)).and.is_finite(data(b,j))))cycle
               tot=tot+1
               if((data(a,i)-data(b,i))*(data(a,j)-data(b,j))>0.0_dp)conc=conc+1
            end do
         end do
         if(tot>0)p(k)=2.0_dp*real(conc,dp)/real(tot,dp)-1.0_dp
      end do
   end function kendall_concurrence


   function empirical_boot_concurrence(data,block_size) result(p)
      real(dp),intent(in)::data(:,:)
      integer,intent(in)::block_size
      real(dp)::p(pair_count(size(data,2)))
      integer::k,i,j,t,l,d,n
      real(dp)::normc
      n=size(data,1)
      normc=log_choose(n,block_size)
      p=0.0_dp
      do k=1,size(p)
         call pair_indices(k,size(data,2),i,j)
         do t=1,n
            d=0
            do l=1,n
               if(data(l,i)<data(t,i) .and. data(l,j)<data(t,j))d=d+1
            end do
            if(d>=block_size-1)p(k)=p(k)+exp(log_choose(d,block_size-1)-normc)
         end do
      end do
   contains
      pure real(dp) function log_choose(nn,kk) result(v)
         integer,intent(in)::nn,kk
         if(kk<0.or.kk>nn)then
         v=-huge(1.0_dp)
         else
         v=log_gamma(real(nn+1,dp))-log_gamma(real(kk+1,dp))-log_gamma(real(nn-kk+1,dp))
         end if
      end function
   end function empirical_boot_concurrence

   function lambda_madogram(data,lambda) result(out)
      real(dp),intent(in)::data(:,:),lambda(:)
      real(dp)::out(pair_count(size(data,2)),size(lambda))
      integer::k,i,j,t,l
      real(dp)::lam
      out=0.0_dp
      do k=1,size(out,1)
         call pair_indices(k,size(data,2),i,j)
         do l=1,size(lambda)
            lam=lambda(l)
            do t=1,size(data,1)
               out(k,l)=out(k,l)+abs(data(t,i)**lam-data(t,j)**(1.0_dp-lam)) &
                    -lam*(1.0_dp-data(t,i)**lam)-(1.0_dp-lam)*(1.0_dp-data(t,j)**(1.0_dp-lam))
            end do
            out(k,l)=0.5_dp*out(k,l)/real(size(data,1),dp) &
                 +0.5_dp*(1.0_dp-lam+lam*lam)/(2.0_dp-lam)/(1.0_dp+lam)
         end do
      end do
   end function lambda_madogram

   pure elemental real(dp) function logit(x,inverse) result(y)
      real(dp),intent(in)::x
      logical,intent(in),optional::inverse
      logical::inv
      inv=.false.
      if(present(inverse))inv=inverse
      if(inv)then
      y=1.0_dp/(1.0_dp+exp(-x))
      else
      y=log(x/(1.0_dp-x))
      end if
   end function logit
end module spatialextremes_dependence
