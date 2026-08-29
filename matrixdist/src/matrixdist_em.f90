! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_em
   use r_compat, only: dp
   use matrixdist_linalg
   use matrixdist_ph, only: ph_exit_rates, ph_loglik
   use matrixdist_dph, only: dph_exit_probs, dph_loglik
   use matrixdist_iph, only: iph_density, iph_cdf
   implicit none
   private
   public :: emstep_ph, fit_ph_em, emstep_dph, fit_dph_em
   public :: iph_loglik

contains

   function outer_product(a,b) result(c)
      real(dp),intent(in)::a(:),b(:)
      real(dp)::c(size(a),size(b))
      integer::i,j
      do i=1,size(a)
      do j=1,size(b)
      c(i,j)=a(i)*b(j)
      end do
      end do
   end function outer_product

   subroutine emstep_ph(alpha,s,obs,weight,rcens,rcweight)
      real(dp),intent(inout)::alpha(:),s(:,:)
      real(dp),intent(in)::obs(:)
      real(dp),intent(in),optional::weight(:),rcens(:),rcweight(:)
      real(dp),allocatable::e(:),exitv(:),bmean(:),zmean(:),nmean(:,:)
      real(dp),allocatable::jmat(:,:),ex(:,:),c(:,:),avec(:),bvec(:),bprod(:,:)
      real(dp)::density,w,total,rate
      integer::p,k,i,j
      p=size(alpha)
      allocate(e(p),bmean(p),zmean(p),nmean(p,p+1),ex(p,p),c(p,p),avec(p),bvec(p))
      e=1.0_dp
      exitv=ph_exit_rates(s)
      bmean=0.0_dp
      zmean=0.0_dp
      nmean=0.0_dp
      total=0.0_dp
      bprod=outer_product(exitv,alpha)
      do k=1,size(obs)
         w=1.0_dp
         if(present(weight))w=weight(k)
         total=total+w
         jmat=matrix_exponential(matrix_vanloan(s,s,bprod)*obs(k))
         ex=jmat(1:p,1:p)
         c=jmat(1:p,p+1:2*p)
         avec=matmul(transpose(ex),alpha)
         bvec=matmul(ex,exitv)
         density=dot_product(alpha,bvec)
         if(density<=0.0_dp)cycle
         do i=1,p
            bmean(i)=bmean(i)+alpha(i)*bvec(i)*w/density
            nmean(i,p+1)=nmean(i,p+1)+avec(i)*exitv(i)*w/density
            zmean(i)=zmean(i)+c(i,i)*w/density
            do j=1,p
            nmean(i,j)=nmean(i,j)+s(i,j)*c(j,i)*w/density
            end do
         end do
      end do
      if(present(rcens))then
         bprod=outer_product(e,alpha)
         do k=1,size(rcens)
            w=1.0_dp
            if(present(rcweight))w=rcweight(k)
            total=total+w
            jmat=matrix_exponential(matrix_vanloan(s,s,bprod)*rcens(k))
            ex=jmat(1:p,1:p)
            c=jmat(1:p,p+1:2*p)
            bvec=matmul(ex,e)
            density=dot_product(alpha,bvec)
            if(density<=0.0_dp)cycle
            do i=1,p
               bmean(i)=bmean(i)+alpha(i)*bvec(i)*w/density
               zmean(i)=zmean(i)+c(i,i)*w/density
               do j=1,p
               nmean(i,j)=nmean(i,j)+s(i,j)*c(j,i)*w/density
               end do
            end do
         end do
      end if
      if(total<=0.0_dp)return
      alpha=max(0.0_dp,bmean/total)
      do i=1,p
         if(zmean(i)<=0.0_dp)cycle
         rate=max(0.0_dp,nmean(i,p+1)/zmean(i))
         s(i,:)=0.0_dp
         s(i,i)=-rate
         do j=1,p
            if(j/=i)then
            s(i,j)=max(0.0_dp,nmean(i,j)/zmean(i))
            s(i,i)=s(i,i)-s(i,j)
            end if
         end do
      end do
   end subroutine emstep_ph

   subroutine fit_ph_em(alpha,s,obs,steps,weight,rcens,rcweight,ll_history)
      real(dp),intent(inout)::alpha(:),s(:,:)
      real(dp),intent(in)::obs(:)
      integer,intent(in)::steps
      real(dp),intent(in),optional::weight(:),rcens(:),rcweight(:)
      real(dp),allocatable,intent(out),optional::ll_history(:)
      integer::it
      if(present(ll_history))allocate(ll_history(steps+1))
      if(present(ll_history))then
         if(present(rcens))then
            ll_history(1)=ph_loglik(alpha,s,obs,weight,rcens,rcweight)
         else
            ll_history(1)=ph_loglik(alpha,s,obs,weight)
         end if
      end if
      do it=1,steps
         if(present(rcens))then
            call emstep_ph(alpha,s,obs,weight,rcens,rcweight)
         else
            call emstep_ph(alpha,s,obs,weight)
         end if
         if(present(ll_history))then
            if(present(rcens))then
               ll_history(it+1)=ph_loglik(alpha,s,obs,weight,rcens,rcweight)
            else
               ll_history(it+1)=ph_loglik(alpha,s,obs,weight)
            end if
         end if
      end do
   end subroutine fit_ph_em

   subroutine emstep_dph(alpha,s,obs,weight)
      real(dp),intent(inout)::alpha(:),s(:,:)
      integer,intent(in)::obs(:)
      real(dp),intent(in),optional::weight(:)
      real(dp),allocatable::exitv(:),bmean(:),nmean(:,:),j0(:,:),jp(:,:),aux(:,:),c(:,:),avec(:),bvec(:)
      real(dp)::density,w,total,rowsum
      integer::p,k,i,j
      p=size(alpha)
      allocate(bmean(p),nmean(p,p+1),aux(p,p),c(p,p),avec(p),bvec(p))
      exitv=dph_exit_probs(s)
      bmean=0.0_dp
      nmean=0.0_dp
      total=0.0_dp
      j0=matrix_vanloan(s,s,outer_product(exitv,alpha))
      do k=1,size(obs)
         if(obs(k)<1)cycle
         w=1.0_dp
         if(present(weight))w=weight(k)
         total=total+w
         jp=matrix_power(obs(k)-1,j0)
         aux=jp(1:p,1:p)
         c=jp(1:p,p+1:2*p)
         avec=matmul(transpose(aux),alpha)
         bvec=matmul(aux,exitv)
         density=dot_product(alpha,bvec)
         if(density<=0.0_dp)cycle
         do i=1,p
            bmean(i)=bmean(i)+alpha(i)*bvec(i)*w/density
            nmean(i,p+1)=nmean(i,p+1)+avec(i)*exitv(i)*w/density
            if(obs(k)>1)then
               do j=1,p
               nmean(i,j)=nmean(i,j)+s(i,j)*c(j,i)*w/density
               end do
            end if
         end do
      end do
      if(total<=0.0_dp)return
      alpha=max(0.0_dp,bmean/total)
      do i=1,p
         rowsum=sum(nmean(i,:))
         if(rowsum<=0.0_dp)cycle
         do j=1,p
         s(i,j)=max(0.0_dp,nmean(i,j)/rowsum)
         end do
      end do
   end subroutine emstep_dph

   subroutine fit_dph_em(alpha,s,obs,steps,weight,ll_history)
      real(dp),intent(inout)::alpha(:),s(:,:)
      integer,intent(in)::obs(:),steps
      real(dp),intent(in),optional::weight(:)
      real(dp),allocatable,intent(out),optional::ll_history(:)
      integer::it
      if(present(ll_history))then
      allocate(ll_history(steps+1))
      ll_history(1)=dph_loglik(alpha,s,obs,weight)
      end if
      do it=1,steps
         call emstep_dph(alpha,s,obs,weight)
         if(present(ll_history))ll_history(it+1)=dph_loglik(alpha,s,obs,weight)
      end do
   end subroutine fit_dph_em

   function iph_loglik(alpha,s,kind,beta,obs,weight,rcens,rcweight) result(ll)
      real(dp),intent(in)::alpha(:),s(:,:),beta(:),obs(:)
      character(len=*),intent(in)::kind
      real(dp),intent(in),optional::weight(:),rcens(:),rcweight(:)
      real(dp)::ll,f,w
      integer::i
      ll=0.0_dp
      do i=1,size(obs)
         w=1.0_dp
         if(present(weight))w=weight(i)
         f=iph_density(obs(i),alpha,s,kind,beta)
         if(f<=0.0_dp)then
         ll=-huge(1.0_dp)
         return
         end if
         ll=ll+w*log(f)
      end do
      if(present(rcens))then
         do i=1,size(rcens)
            w=1.0_dp
            if(present(rcweight))w=rcweight(i)
            f=iph_cdf(rcens(i),alpha,s,kind,beta,.false.)
            if(f<=0.0_dp)then
            ll=-huge(1.0_dp)
            return
            end if
            ll=ll+w*log(f)
         end do
      end if
   end function iph_loglik

end module matrixdist_em
