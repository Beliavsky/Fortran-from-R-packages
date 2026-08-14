module nspmix_core
   use nspmix_kinds, only : dp
   use nspmix_types
   use nspmix_utils, only : normalize_prob, logsumexp_vec
   use nspmix_families, only : logd_eval, data_weights, nobs_data
   use lsei, only : lsei_solve, ls_result, LSEI_SUCCESS
   implicit none
   private
   public :: loglik, mixture_logdensity, mixture_density, gradient_values, hcnm
contains
   subroutine mixture_logdensity(data,beta,mix,logdmix)
      type(nsp_data),intent(in)::data; real(dp),intent(in)::beta(:); type(disc_dist),intent(in)::mix
      real(dp),allocatable,intent(out)::logdmix(:)
      real(dp),allocatable::ld(:,:),dt(:,:),db(:,:,:),z(:)
      integer::i,j,n,k
      call logd_eval(data,beta,mix%pt,ld,dt,db); n=size(ld,1); k=size(ld,2)
      allocate(logdmix(n),z(k))
      do i=1,n
         do j=1,k
            if(mix%pr(j)>0.0_dp) then; z(j)=ld(i,j)+log(mix%pr(j)); else; z(j)=-huge(1.0_dp); end if
         end do
         logdmix(i)=logsumexp_vec(z)
      end do
   end subroutine

   subroutine mixture_density(data,beta,mix,d)
      type(nsp_data),intent(in)::data; real(dp),intent(in)::beta(:); type(disc_dist),intent(in)::mix
      real(dp),allocatable,intent(out)::d(:)
      real(dp),allocatable::l(:)
      call mixture_logdensity(data,beta,mix,l); allocate(d(size(l))); d=exp(l)
   end subroutine

   real(dp) function loglik(data,beta,mix)
      type(nsp_data),intent(in)::data; real(dp),intent(in)::beta(:); type(disc_dist),intent(in)::mix
      real(dp),allocatable::l(:),w(:)
      call mixture_logdensity(data,beta,mix,l); call data_weights(data,beta,w)
      loglik=sum(w*l)
   end function

   subroutine gradient_values(pt,data,beta,mix,g,dg)
      real(dp),intent(in)::pt(:); type(nsp_data),intent(in)::data; real(dp),intent(in)::beta(:); type(disc_dist),intent(in)::mix
      real(dp),allocatable,intent(out)::g(:),dg(:)
      real(dp),allocatable::ld(:,:),dt(:,:),db(:,:,:),lm(:),w(:)
      integer::i,j
      call logd_eval(data,beta,pt,ld,dt,db); call mixture_logdensity(data,beta,mix,lm); call data_weights(data,beta,w)
      allocate(g(size(pt)),dg(size(pt))); g=-sum(w); dg=0.0_dp
      do j=1,size(pt)
         do i=1,size(w)
            g(j)=g(j)+w(i)*min(exp(ld(i,j)-lm(i)),1.0e100_dp)
            dg(j)=dg(j)+w(i)*min(exp(ld(i,j)-lm(i)),1.0e100_dp)*dt(i,j)
         end do
      end do
   end subroutine

   subroutine hcnm(d,p0,w,res,maxit,tol)
      real(dp),intent(in)::d(:,:)
      real(dp),intent(in),optional::p0(:),w(:)
      type(hcnm_result),intent(out)::res
      integer,intent(in),optional::maxit
      real(dp),intent(in),optional::tol
      real(dp),allocatable::p(:),ww(:),wr(:),prob(:),s(:,:),g(:),a(:,:),b(:),pnew(:),trial(:)
      real(dp),allocatable::ceq(:,:),deq(:),lower(:)
      real(dp)::ll,llold,rise,alpha,eps,gn,sw
      integer::it,lim,m,n
      type(ls_result)::lr
      n=size(d,1); m=size(d,2); lim=1000; if(present(maxit)) lim=maxit; eps=1.0e-6_dp; if(present(tol)) eps=tol
      allocate(ww(n),wr(n),p(m),prob(n),s(n,m),g(m),a(n,m),b(n),pnew(m),trial(m))
      allocate(ceq(1,m),deq(1),lower(m)); ceq=1.0_dp; deq=1.0_dp; lower=0.0_dp
      ww=1.0_dp; if(present(w)) then; if(size(w)==1) then; ww=w(1); else; ww=w; end if; end if
      wr=sqrt(max(ww,0.0_dp)); sw=sum(ww)
      if(present(p0)) then; p=p0; else; p=sum(d*spread(ww,2,m),dim=1); end if
      call normalize_prob(p); prob=max(matmul(d,p),1.0e-300_dp); ll=sum(ww*log(prob)); res%convergence=1
      do it=1,lim
         llold=ll; s=d/spread(prob,2,m); g=matmul(transpose(s),ww)
         a=s*spread(wr,2,m); b=2.0_dp*wr
         call lsei_solve(a,b,c=ceq,d=deq,res=lr,lower=lower,tol=1.0e-12_dp)
         if(lr%mode==LSEI_SUCCESS .and. allocated(lr%x)) then
            pnew=max(lr%x,0.0_dp); call normalize_prob(pnew)
         else
            pnew=p*(g/max(sw,1.0e-300_dp)); call normalize_prob(pnew)
         end if
         rise=dot_product(g,pnew-p); alpha=1.0_dp
         do
            trial=p+alpha*(pnew-p); trial=max(trial,0.0_dp); call normalize_prob(trial)
            prob=max(matmul(d,trial),1.0e-300_dp); ll=sum(ww*log(prob))
            if(ll>=llold+0.33_dp*alpha*rise .or. alpha<1.0e-10_dp) exit
            alpha=0.5_dp*alpha
         end do
         if(ll<llold) then; ll=llold; prob=max(matmul(d,p),1.0e-300_dp); else; p=trial; end if
         if(it>2 .and. ll<=llold+eps) then; res%convergence=0; exit; end if
      end do
      s=d/spread(prob,2,m); g=matmul(transpose(s),ww); gn=maxval(g)-sw
      res%p=p; res%ll=ll; res%maxgrad=gn; res%iterations=min(it,lim)
   end subroutine
end module nspmix_core
