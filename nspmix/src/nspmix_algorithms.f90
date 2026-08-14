module nspmix_algorithms
   use nspmix_kinds, only : dp
   use nspmix_types
   use nspmix_utils, only : make_disc, normalize_prob, collapse_disc
   use nspmix_families, only : default_beta, default_mix_points, gridpoints, logd_eval, data_weights, valid_parameters
   use nspmix_core, only : loglik, gradient_values, hcnm
   implicit none
   private
   public :: cnm, cnm_proportions, cnmms, cnmpl, cnmap, maxgrad
contains
   subroutine fit_masses(data,beta,mix,tol,maxit,ll,maxg)
      type(nsp_data),intent(in)::data; real(dp),intent(in)::beta(:); type(disc_dist),intent(inout)::mix
      real(dp),intent(in)::tol; integer,intent(in)::maxit; real(dp),intent(out)::ll,maxg
      real(dp),allocatable::ld(:,:),dt(:,:),db(:,:,:),ma(:),d(:,:),w(:)
      type(hcnm_result)::hr
      integer::i
      call logd_eval(data,beta,mix%pt,ld,dt,db); allocate(ma(size(ld,1)),d(size(ld,1),size(ld,2)))
      do i=1,size(ld,1); ma(i)=maxval(ld(i,:)); d(i,:)=exp(ld(i,:)-ma(i)); end do
      call data_weights(data,beta,w); call hcnm(d,mix%pr,w,hr,maxit=maxit,tol=tol)
      mix%pr=hr%p; call normalize_prob(mix%pr); ll=hr%ll+dot_product(w,ma); maxg=hr%maxgrad
   end subroutine

   subroutine cnm_proportions(data,points,res,beta,p0,maxit,tol)
      type(nsp_data),intent(in)::data; real(dp),intent(in)::points(:)
      type(nspmix_result),intent(out)::res
      real(dp),intent(in),optional::beta(:),p0(:); integer,intent(in),optional::maxit; real(dp),intent(in),optional::tol
      real(dp),allocatable::b(:),p(:); integer::lim; real(dp)::eps,g
      if(present(beta)) then; b=beta; else; call default_beta(data,b); end if
      allocate(p(size(points))); p=1.0_dp; if(present(p0)) p=p0; call normalize_prob(p)
      call make_disc(points,p,res%mix); lim=1000; if(present(maxit)) lim=maxit; eps=1e-6_dp; if(present(tol)) eps=tol
      call fit_masses(data,b,res%mix,eps,lim,res%ll,g); res%beta=b; res%max_gradient=g; res%iterations=lim; res%convergence=0
   end subroutine

   subroutine cnm(data,res,beta,init_mix,maxit,tol,ngrid,kmax)
      type(nsp_data),intent(in)::data; type(nspmix_result),intent(out)::res
      real(dp),intent(in),optional::beta(:); type(disc_dist),intent(in),optional::init_mix
      integer,intent(in),optional::maxit,ngrid,kmax; real(dp),intent(in),optional::tol
      type(disc_dist)::mix,oldmix
      real(dp),allocatable::b(:),gp(:),g(:),dg(:),pt(:),pr(:)
      real(dp)::ll,llold,maxg,eps,ctol
      integer::it,lim,ng,km,nadd,imax
      if(present(beta)) then; b=beta; else; call default_beta(data,b); end if
      lim=100; if(present(maxit)) lim=maxit; eps=1.0e-6_dp; if(present(tol)) eps=tol
      ng=100; if(present(ngrid)) ng=ngrid; km=huge(1); if(present(kmax)) km=kmax
      if(present(init_mix)) then; mix=init_mix
      else
         call default_mix_points(data,b,min(10,km),pt); call make_disc(pt,d=mix)
         deallocate(pt)
      end if
      if(.not.valid_parameters(data,b,mix%pt)) error stop "cnm: invalid initial parameters"
      call fit_masses(data,b,mix,eps,1000,ll,maxg); res%convergence=1
      do it=1,lim
         llold=ll; oldmix=mix
         call gridpoints(data,b,ng,gp); call gradient_values(gp,data,b,mix,g,dg); maxg=maxval(g)
         if(size(mix%pt)>=km .or. maxg<=max(1.0e-3_dp,10.0_dp*eps)) then
            if(it>1 .and. ll<=llold+eps) then; res%convergence=0; exit; end if
         end if
         if(size(mix%pt)<km) then
            imax=maxloc(g,dim=1)
            if(g(imax)>0.0_dp) then
               nadd=1; allocate(pt(size(mix%pt)+nadd),pr(size(mix%pr)+nadd))
               pt(:size(mix%pt))=mix%pt; pr(:size(mix%pr))=mix%pr
               pt(size(pt))=gp(imax); pr(size(pr))=0.0_dp
               call make_disc(pt,pr,mix,collapse_tol=0.0_dp); deallocate(pt,pr)
            end if
         end if
         call fit_masses(data,b,mix,eps,1000,ll,maxg)
         ctol=max(1.0e-12_dp,eps*1.0e-2_dp); call prune_mix(mix,ctol)
         ll=loglik(data,b,mix)
         if(ll<=llold+eps) then; res%convergence=0; exit; end if
      end do
      call gradient_values(mix%pt,data,b,mix,g,dg)
      res%mix=mix; res%beta=b; res%ll=loglik(data,b,mix); res%max_gradient=maxg
      res%iterations=min(it,lim)
   end subroutine

   subroutine maxgrad(data,beta,mix,ngrid,pt,g,gmax)
      type(nsp_data),intent(in)::data; real(dp),intent(in)::beta(:); type(disc_dist),intent(in)::mix
      integer,intent(in)::ngrid; real(dp),allocatable,intent(out)::pt(:),g(:); real(dp),intent(out)::gmax
      real(dp),allocatable::grid(:),gg(:),dg(:),p2(:),g2(:)
      integer::i,k
      call gridpoints(data,beta,ngrid,grid); call gradient_values(grid,data,beta,mix,gg,dg)
      allocate(p2(size(grid)),g2(size(grid))); k=0
      do i=1,size(grid)
         if((i==1 .or. gg(i)>=gg(max(1,i-1))) .and. (i==size(grid) .or. gg(i)>=gg(min(size(grid),i+1)))) then
            k=k+1; p2(k)=grid(i); g2(k)=gg(i)
         end if
      end do
      pt=p2(:k); g=g2(:k); gmax=maxval(gg)
   end subroutine

   subroutine cnmms(data,res,beta,init_mix,maxit,tol,ngrid,kmax)
      type(nsp_data),intent(in)::data; type(nspmix_result),intent(out)::res
      real(dp),intent(in),optional::beta(:); type(disc_dist),intent(in),optional::init_mix
      integer,intent(in),optional::maxit,ngrid,kmax; real(dp),intent(in),optional::tol
      real(dp),allocatable::b(:); type(nspmix_result)::r0,r1
      integer::it,lim,ng,km; real(dp)::eps,llold
      if(present(beta)) then; b=beta; else; call default_beta(data,b); end if
      lim=100; if(present(maxit)) lim=maxit; ng=100; if(present(ngrid)) ng=ngrid; km=100; if(present(kmax)) km=kmax
      eps=1e-6_dp; if(present(tol)) eps=tol
      if(present(init_mix)) then; call cnm(data,r0,b,init_mix,maxit=20,tol=eps*0.1_dp,ngrid=ng,kmax=km)
      else; call cnm(data,r0,b,maxit=20,tol=eps*0.1_dp,ngrid=ng,kmax=km); end if
      if(size(b)==0) then; res=r0; return; end if
      do it=1,lim
         llold=r0%ll; call profile_beta_step(data,r0,b,ng,km,eps,r1); r0=r1; b=r0%beta
         if(r0%ll<=llold+eps) exit
      end do
      res=r0; res%iterations=it; res%convergence=merge(0,1,it<lim)
   end subroutine

   subroutine cnmpl(data,res,beta,init_mix,maxit,tol,ngrid)
      type(nsp_data),intent(in)::data; type(nspmix_result),intent(out)::res
      real(dp),intent(in),optional::beta(:); type(disc_dist),intent(in),optional::init_mix
      integer,intent(in),optional::maxit,ngrid; real(dp),intent(in),optional::tol
      if(present(init_mix)) then
         call cnmms(data,res,beta=beta,init_mix=init_mix,maxit=maxit,tol=tol,ngrid=ngrid)
      else
         call cnmms(data,res,beta=beta,maxit=maxit,tol=tol,ngrid=ngrid)
      end if
   end subroutine

   subroutine cnmap(data,res,beta,init_mix,maxit,tol,ngrid)
      type(nsp_data),intent(in)::data; type(nspmix_result),intent(out)::res
      real(dp),intent(in),optional::beta(:); type(disc_dist),intent(in),optional::init_mix
      integer,intent(in),optional::maxit,ngrid; real(dp),intent(in),optional::tol
      if(present(init_mix)) then
         call cnmms(data,res,beta=beta,init_mix=init_mix,maxit=maxit,tol=tol,ngrid=ngrid)
      else
         call cnmms(data,res,beta=beta,maxit=maxit,tol=tol,ngrid=ngrid)
      end if
   end subroutine

   subroutine profile_beta_step(data,r0,beta,ng,km,eps,rbest)
      type(nsp_data),intent(in)::data; type(nspmix_result),intent(in)::r0
      real(dp),intent(in)::beta(:); integer,intent(in)::ng,km; real(dp),intent(in)::eps
      type(nspmix_result),intent(out)::rbest
      real(dp),allocatable::btry(:); type(nspmix_result)::rt
      real(dp)::h; integer::j,sgn
      rbest=r0; allocate(btry(size(beta)))
      do j=1,size(beta)
         h=max(1.0e-4_dp,1.0e-3_dp*max(abs(beta(j)),1.0_dp))
         do sgn=-1,1,2
            btry=beta; btry(j)=beta(j)+real(sgn,dp)*h
            if((data%family==NSP_NORMAL .or. data%family==NSP_CVPS) .and. btry(j)<=1.0e-8_dp) cycle
            call cnm(data,rt,btry,r0%mix,maxit=15,tol=eps*0.1_dp,ngrid=ng,kmax=km)
            if(rt%ll>rbest%ll) rbest=rt
         end do
      end do
   end subroutine

   subroutine prune_mix(mix,tol)
      type(disc_dist),intent(inout)::mix; real(dp),intent(in)::tol
      real(dp),allocatable::pt(:),pr(:); integer::i,k
      allocate(pt(size(mix%pt)),pr(size(mix%pr))); k=0
      do i=1,size(mix%pt)
         if(mix%pr(i)>tol) then; k=k+1; pt(k)=mix%pt(i); pr(k)=mix%pr(i); end if
      end do
      if(k==0) then; k=maxloc(mix%pr,dim=1); call make_disc([mix%pt(k)],[1.0_dp],mix)
      else; call make_disc(pt(:k),pr(:k),mix,collapse_tol=tol); end if
   end subroutine
end module nspmix_algorithms
