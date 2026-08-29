! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_multi_fit
   use r_compat, only: dp
   use matrixdist_linalg, only: matrix_exponential, matrix_vanloan, matrix_power
   use matrixdist_ph, only: ph_exit_rates
   use matrixdist_dph, only: dph_exit_probs
   use matrixdist_multivariate, only: mph_density_point, mdph_density_point
   implicit none
   private
   public :: mph_loglik, mdph_loglik, emstep_mph_rc, emstep_mdph
contains
   function mph_loglik(y,delta,alpha,s,weight) result(ll)
      real(dp),intent(in)::y(:,:),alpha(:),s(:,:,:)
      logical,intent(in)::delta(:,:)
      real(dp),intent(in),optional::weight(:)
      real(dp)::ll,w,f
      integer::i
      if(any(shape(y)/=shape(delta))) error stop 'mph_loglik: y/delta mismatch'
      ll=0.0_dp
      do i=1,size(y,1)
         w=1.0_dp
         if(present(weight))w=weight(i)
         f=mph_density_point(y(i,:),alpha,s,delta(i,:))
         if(f<=0.0_dp) then
         ll=-huge(1.0_dp)
         return
         end if
         ll=ll+w*log(f)
      end do
   end function mph_loglik

   function mdph_loglik(obs,alpha,s,weight) result(ll)
      integer,intent(in)::obs(:,:)
      real(dp),intent(in)::alpha(:),s(:,:,:)
      real(dp),intent(in),optional::weight(:)
      real(dp)::ll,w,f
      integer::i
      ll=0.0_dp
      do i=1,size(obs,1)
         w=1.0_dp
         if(present(weight))w=weight(i)
         f=mdph_density_point(obs(i,:),alpha,s)
         if(f<=0.0_dp) then
         ll=-huge(1.0_dp)
         return
         end if
         ll=ll+w*log(f)
      end do
   end function mdph_loglik

   subroutine emstep_mph_rc(alpha,s,y,delta,weight)
      ! Exact multivariate-PH EM step.  Conditional on the common initial
      ! state, margins are independent; this formulation is algebraically
      ! equivalent to the native mPH_EM_UNI implementation but reuses the
      ! common Van-Loan kernel and avoids its large temporary tensors.
      real(dp),intent(inout)::alpha(:),s(:,:,:)
      real(dp),intent(in)::y(:,:)
      logical,intent(in)::delta(:,:)
      real(dp),intent(in),optional::weight(:)
      integer::n,p,d,m,j,h,i,l
      real(dp)::w,sumw,total,post,lik
      real(dp),allocatable::comp(:,:),bmean(:),zmean(:,:),nabs(:,:),ntrans(:,:,:)
      real(dp),allocatable::e(:),exitv(:),unit(:),b(:,:),vl(:,:),ev(:,:),c(:,:)
      n=size(y,1)
      d=size(y,2)
      p=size(alpha)
      if(size(delta,1)/=n .or. size(delta,2)/=d .or. size(s,1)/=p .or. size(s,2)/=p .or. size(s,3)/=d) &
         error stop 'emstep_mph_rc: dimension mismatch'
      allocate(comp(p,d),bmean(p),zmean(p,d),nabs(p,d),ntrans(p,p,d),e(p),unit(p))
      bmean=0.0_dp
      zmean=0.0_dp
      nabs=0.0_dp
      ntrans=0.0_dp
      e=1.0_dp
      sumw=0.0_dp
      do m=1,n
         w=1.0_dp
         if(present(weight))w=weight(m)
         sumw=sumw+w
         do j=1,d
            ev=matrix_exponential(s(:,:,j)*y(m,j))
            exitv=ph_exit_rates(s(:,:,j))
            do h=1,p
               if(delta(m,j)) then
                  comp(h,j)=dot_product(ev(h,:),exitv)
               else
                  comp(h,j)=sum(ev(h,:))
               end if
            end do
         end do
         total=0.0_dp
         do h=1,p
         total=total+alpha(h)*product(comp(h,:))
         end do
         if(total<=tiny(1.0_dp)) cycle
         do h=1,p
            post=alpha(h)*product(comp(h,:))/total
            bmean(h)=bmean(h)+w*post
            unit=0.0_dp
            unit(h)=1.0_dp
            do j=1,d
               lik=comp(h,j)
               if(lik<=tiny(1.0_dp)) cycle
               exitv=ph_exit_rates(s(:,:,j))
               if(delta(m,j)) then
                  allocate(b(p,p))
                  b=0.0_dp
                  do i=1,p
                  do l=1,p
                  b(i,l)=exitv(i)*unit(l)
                  end do
                  end do
               else
                  allocate(b(p,p))
                  b=0.0_dp
                  do i=1,p
                  do l=1,p
                  b(i,l)=unit(l)
                  end do
                  end do
               end if
               vl=matrix_exponential(matrix_vanloan(s(:,:,j),s(:,:,j),b)*y(m,j))
               allocate(c(p,p))
               c=vl(1:p,p+1:2*p)
               do i=1,p
                  zmean(i,j)=zmean(i,j)+w*post*c(i,i)/lik
                  do l=1,p
                     if(l/=i) ntrans(i,l,j)=ntrans(i,l,j)+w*post*s(i,l,j)*c(l,i)/lik
                  end do
               end do
               if(delta(m,j)) then
                  ev=vl(1:p,1:p)
                  do i=1,p
                     nabs(i,j)=nabs(i,j)+w*post*ev(h,i)*exitv(i)/lik
                  end do
               end if
               deallocate(b,c)
            end do
         end do
      end do
      if(sumw>0.0_dp) alpha=max(0.0_dp,bmean/sumw)
      if(sum(alpha)>0.0_dp) alpha=alpha/sum(alpha)
      do j=1,d
         do i=1,p
            if(zmean(i,j)<=tiny(1.0_dp)) cycle
            s(i,i,j)=-max(0.0_dp,nabs(i,j)/zmean(i,j))
            do l=1,p
               if(l==i)cycle
               s(i,l,j)=max(0.0_dp,ntrans(i,l,j)/zmean(i,j))
               s(i,i,j)=s(i,i,j)-s(i,l,j)
            end do
         end do
      end do
   end subroutine emstep_mph_rc

   subroutine emstep_mdph(alpha,s,obs,weight)
      ! Port of upstream EMstep_mdph, using direct matrix powers instead of
      ! caching vector_of_powers.  This trades some speed for much less
      ! storage and a simple auditable Fortran implementation.
      real(dp),intent(inout)::alpha(:),s(:,:,:)
      integer,intent(in)::obs(:,:)
      real(dp),intent(in),optional::weight(:)
      integer::n,p,d,m,i,j,l,h,q
      real(dp)::w,sumw,density,factor
      real(dp),allocatable::den(:,:),bmean(:),nabs(:,:),ntrans(:,:,:),exitv(:)
      real(dp),allocatable::powk(:,:),powm(:,:),tail(:)
      n=size(obs,1)
      d=size(obs,2)
      p=size(alpha)
      if(size(s,1)/=p .or. size(s,2)/=p .or. size(s,3)/=d) error stop 'emstep_mdph: dimension mismatch'
      allocate(den(p,d),bmean(p),nabs(p,d),ntrans(p,p,d))
      bmean=0.0_dp
      nabs=0.0_dp
      ntrans=0.0_dp
      sumw=0.0_dp
      do m=1,n
         w=1.0_dp
         if(present(weight))w=weight(m)
         sumw=sumw+w
         do j=1,d
            exitv=dph_exit_probs(s(:,:,j))
            powk=matrix_power(obs(m,j)-1,s(:,:,j))
            do h=1,p
            den(h,j)=dot_product(powk(h,:),exitv)
            end do
         end do
         density=0.0_dp
         do h=1,p
         density=density+alpha(h)*product(den(h,:))
         end do
         if(density<=tiny(1.0_dp))cycle
         do i=1,p
            bmean(i)=bmean(i)+w*alpha(i)*product(den(i,:))/density
            do j=1,d
               ! Absorption from i: sum over common initial states h.
               factor=0.0_dp
               powk=matrix_power(obs(m,j)-1,s(:,:,j))
               do h=1,p
                  other_product: block
                     real(dp)::op
                     integer::qq
                     op=1.0_dp
                     do qq=1,d
                     if(qq/=j)op=op*den(h,qq)
                     end do
                     factor=factor+alpha(h)*op*powk(h,i)
                  end block other_product
               end do
               exitv=dph_exit_probs(s(:,:,j))
               nabs(i,j)=nabs(i,j)+w*exitv(i)*factor/density
               if(obs(m,j)>1) then
                  do l=1,p
                     factor=0.0_dp
                     do q=0,obs(m,j)-2
                        powm=matrix_power(q,s(:,:,j))
                        tail=matmul(matrix_power(obs(m,j)-q-2,s(:,:,j)),exitv)
                        do h=1,p
                           other_prod_h: block
                              real(dp)::op
                              integer::qq
                              op=1.0_dp
                              do qq=1,d
                              if(qq/=j)op=op*den(h,qq)
                              end do
                              factor=factor+alpha(h)*op*powm(h,i)*tail(l)
                           end block other_prod_h
                        end do
                     end do
                     ntrans(i,l,j)=ntrans(i,l,j)+w*s(i,l,j)*factor/density
                  end do
               end if
            end do
         end do
      end do
      if(sumw>0.0_dp)alpha=max(0.0_dp,bmean/sumw)
      if(sum(alpha)>0.0_dp)alpha=alpha/sum(alpha)
      do j=1,d
         do i=1,p
            factor=nabs(i,j)+sum(ntrans(i,:,j))
            if(factor<=tiny(1.0_dp))cycle
            do l=1,p
            s(i,l,j)=max(0.0_dp,ntrans(i,l,j)/factor)
            end do
         end do
      end do
   end subroutine emstep_mdph
end module matrixdist_multi_fit
