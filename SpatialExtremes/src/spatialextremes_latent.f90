module spatialextremes_latent
   use spatialextremes_base, only: dp,solve_spd,logdet_spd,chol_upper,pi
   use spatialextremes_univariate, only: gev_loglik
   use spatialextremes_covariance, only: covariance_matrix
   use r_compat, only: runif1,rnorm1,rgamma
   implicit none
   private

   type,public :: dic_result_t
      real(dp)::dic=0.0_dp
      real(dp)::effective_npar=0.0_dp
      real(dp)::dbar=0.0_dp
   end type dic_result_t

   type,public :: latent_mcmc_result_t
      real(dp),allocatable :: chain_loc(:,:),chain_scale(:,:),chain_shape(:,:)
      real(dp) :: acceptance(9)=0.0_dp
      real(dp) :: extreme_rejection(9)=0.0_dp
      integer :: iterations=0
      integer :: info=0
   end type latent_mcmc_result_t

   public :: latent_dic, latent_gev_loglik, gaussian_field_logdensity, latent_gev_mcmc
contains
   real(dp) function latent_gev_loglik(data,loc,scale,shape) result(ll)
      real(dp),intent(in)::data(:,:),loc(:),scale(:),shape(:)
      integer::j
      ll=0.0_dp
      do j=1,size(data,2)
         ll=ll+gev_loglik(data(:,j),loc(j),scale(j),shape(j))
      end do
   end function latent_gev_loglik

   function latent_dic(data,chain_loc,chain_scale,chain_shape,post_loc,post_scale,post_shape) result(ans)
      real(dp),intent(in)::data(:,:),chain_loc(:,:),chain_scale(:,:),chain_shape(:,:)
      real(dp),intent(in),optional::post_loc(:),post_scale(:),post_shape(:)
      type(dic_result_t)::ans
      real(dp)::pl(size(data,2)),ps(size(data,2)),ph(size(data,2)),tmp
      integer::i,j,nchain
      nchain=size(chain_loc,1)
      if(present(post_loc))then
         pl=post_loc
      else
         pl=sum(chain_loc,dim=1)/real(nchain,dp)
      end if
      if(present(post_scale))then
         ps=post_scale
      else
         ps=sum(chain_scale,dim=1)/real(nchain,dp)
      end if
      if(present(post_shape))then
         ph=post_shape
      else
         ph=sum(chain_shape,dim=1)/real(nchain,dp)
      end if
      tmp=0.0_dp
      do i=1,nchain
         do j=1,size(data,2)
            tmp=tmp+gev_loglik(data(:,j),chain_loc(i,j),chain_scale(i,j),chain_shape(i,j))
         end do
      end do
      ans%dbar=-2.0_dp*tmp/real(nchain,dp)
      tmp=latent_gev_loglik(data,pl,ps,ph)
      ans%effective_npar=ans%dbar+2.0_dp*tmp
      ans%dic=ans%dbar+ans%effective_npar
   end function latent_dic

   real(dp) function gaussian_field_logdensity(x,mean,cov,info) result(lp)
      real(dp),intent(in)::x(:),mean(:),cov(:,:)
      integer,intent(out),optional::info
      real(dp)::rhs(size(x),1),sol(size(x),1),ld
      integer::istat
      rhs(:,1)=x-mean
      call solve_spd(cov,rhs,sol,istat)
      if(istat/=0)then
         lp=-huge(1.0_dp)
         if(present(info))info=istat
         return
      end if
      ld=logdet_spd(cov,istat)
      if(istat/=0)then
         lp=-huge(1.0_dp)
         if(present(info))info=istat
         return
      end if
      lp=-0.5_dp*(real(size(x),dp)*log(2.0_dp*pi)+ld+dot_product(rhs(:,1),sol(:,1)))
      if(present(info))info=0
   end function gaussian_field_logdensity

   subroutine latent_gev_mcmc(data,coord,cov_model,design_loc,design_scale,design_shape, &
      beta_loc,beta_scale,beta_shape,sills,ranges,smooths,gev_start,hyper_sill,hyper_range, &
      hyper_smooth,beta_mean_loc,beta_mean_scale,beta_mean_shape,beta_prec_loc,beta_prec_scale, &
      beta_prec_shape,prop_gev,prop_ranges,prop_smooths,n_keep,result,thin,burn_in,use_log_link)
      ! Computational translation of upstream latentgev().
      ! hyper_* rows are (shape, scale) Gamma/Inv-Gamma prior parameters for each margin.
      real(dp),intent(in)::data(:,:),coord(:,:)
      integer,intent(in)::cov_model(3),n_keep
      real(dp),intent(in)::design_loc(:,:),design_scale(:,:),design_shape(:,:)
      real(dp),intent(inout)::beta_loc(:),beta_scale(:),beta_shape(:)
      real(dp),intent(inout)::sills(3),ranges(3),smooths(3),gev_start(:,:)
      real(dp),intent(in)::hyper_sill(2,3),hyper_range(2,3),hyper_smooth(2,3)
      real(dp),intent(in)::beta_mean_loc(:),beta_mean_scale(:),beta_mean_shape(:)
      real(dp),intent(in)::beta_prec_loc(:,:),beta_prec_scale(:,:),beta_prec_shape(:,:)
      real(dp),intent(in)::prop_gev(3),prop_ranges(3),prop_smooths(3)
      type(latent_mcmc_result_t),intent(out)::result
      integer,intent(in),optional::thin,burn_in
      logical,intent(in),optional::use_log_link
      real(dp),allocatable::gpmean(:,:),covs(:,:,:),proposal_cov(:,:),z(:),mu(:), &
         ywork(:),proposal(:),beta_all(:)
      integer::nsite,nobs,ithin,iburn,iter,stored,j,m,istat,total_beta,offset(3),nb(3)
      logical::log_link
      real(dp)::oldll,newll,oldgp,newgp,logratio,u,oldv,newv,propv,lp_old,lp_new

      nobs=size(data,1)
      nsite=size(data,2)
      nb=[size(beta_loc),size(beta_scale),size(beta_shape)]
      if(size(coord,1)/=nsite .or. size(gev_start,1)/=nsite .or. size(gev_start,2)/=3)then
         result%info=1
         return
      end if
      if(size(design_loc,1)/=nsite .or. size(design_scale,1)/=nsite .or. size(design_shape,1)/=nsite)then
         result%info=2
         return
      end if
      if(size(design_loc,2)/=nb(1) .or. size(design_scale,2)/=nb(2) .or. size(design_shape,2)/=nb(3))then
         result%info=3
         return
      end if
      if(any(sills<=0.0_dp) .or. any(ranges<=0.0_dp) .or. any(smooths<=0.0_dp) .or. any(gev_start(:,2)<=0.0_dp))then
         result%info=4
         return
      end if
      ithin=1
      if(present(thin))ithin=max(1,thin)
      iburn=0
      if(present(burn_in))iburn=max(0,burn_in)
      log_link=.false.
      if(present(use_log_link))log_link=use_log_link
      total_beta=sum(nb)
      offset=[0,nb(1),nb(1)+nb(2)]
      allocate(gpmean(nsite,3),covs(nsite,nsite,3),proposal_cov(nsite,nsite),z(nsite),mu(nsite), &
         ywork(nsite),proposal(3),beta_all(total_beta))
      allocate(result%chain_loc(n_keep,nb(1)+3+nsite),result%chain_scale(n_keep,nb(2)+3+nsite), &
         result%chain_shape(n_keep,nb(3)+3+nsite))
      result%chain_loc=0.0_dp
      result%chain_scale=0.0_dp
      result%chain_shape=0.0_dp
      result%acceptance=0.0_dp
      result%extreme_rejection=0.0_dp
      result%info=0
      beta_all=[beta_loc,beta_scale,beta_shape]
      call update_gpmeans()
      do m=1,3
         covs(:,:,m)=covariance_matrix(coord,cov_model(m),0.0_dp,sills(m),ranges(m),smooths(m))
         lp_old=gaussian_field_logdensity(gp_values(m),gpmean(:,m),covs(:,:,m),istat)
         if(istat/=0)then
         result%info=10+m
         return
         end if
      end do

      iter=0
      stored=0
      do while(stored<n_keep)
         iter=iter+1
         ! Sitewise GEV parameters.
         do j=1,nsite
            do m=1,3
               proposal=gev_start(j,:)
               oldv=proposal(m)
               if(m==2)then
                  propv=exp(log(oldv)+prop_gev(m)*rnorm1())
               else
                  propv=oldv+prop_gev(m)*rnorm1()
               end if
               proposal(m)=propv
               newll=gev_loglik(data(:,j),proposal(1),proposal(2),proposal(3))
               if(newll<=-9.0e5_dp)then
                  result%extreme_rejection(m)=result%extreme_rejection(m)+1.0_dp
                  cycle
               end if
               oldll=gev_loglik(data(:,j),gev_start(j,1),gev_start(j,2),gev_start(j,3))
               z=gp_values(m)
               newv=propv
               if(m==2 .and. log_link)then
                  z(j)=log(newv)
               else
                  z(j)=newv
               end if
               newgp=gaussian_field_logdensity(z,gpmean(:,m),covs(:,:,m),istat)
               z=gp_values(m)
               oldgp=gaussian_field_logdensity(z,gpmean(:,m),covs(:,:,m),istat)
               logratio=newll-oldll+newgp-oldgp
               if(m==2)logratio=logratio+log(propv/oldv)
               if(log_link .and. m==2)logratio=logratio-log(propv)+log(oldv)
               if(log(max(runif1(),tiny(1.0_dp)))<min(0.0_dp,logratio))then
                  gev_start(j,m)=propv
                  result%acceptance(m)=result%acceptance(m)+1.0_dp
               end if
            end do
         end do

         ! Conjugate regression updates.
         do m=1,3
            select case(m)
            case(1)
               ywork=gev_start(:,1)
               call update_beta(design_loc,ywork,covs(:,:,1),beta_mean_loc,beta_prec_loc,beta_loc,istat)
            case(2)
               ywork=gev_start(:,2)
               if(log_link)ywork=log(ywork)
               call update_beta(design_scale,ywork,covs(:,:,2),beta_mean_scale,beta_prec_scale,beta_scale,istat)
            case(3)
               ywork=gev_start(:,3)
               call update_beta(design_shape,ywork,covs(:,:,3),beta_mean_shape,beta_prec_shape,beta_shape,istat)
            end select
            if(istat/=0)then
            result%info=20+m
            return
            end if
         end do
         beta_all=[beta_loc,beta_scale,beta_shape]
         call update_gpmeans()

         ! Conjugate inverse-gamma sill updates.
         do m=1,3
            z=gp_values(m)
            call normalized_quad(z-gpmean(:,m),covs(:,:,m),sills(m),oldv,istat)
            if(istat/=0)then
            result%info=30+m
            return
            end if
            newv=0.5_dp*real(nsite,dp)+hyper_sill(1,m)
            lp_new=hyper_sill(2,m)+0.5_dp*oldv
            proposal=rgamma(1,newv,1.0_dp/lp_new)
            if(proposal(1)<=0.0_dp)then
            result%info=40+m
            return
            end if
            sills(m)=1.0_dp/proposal(1)
            covs(:,:,m)=covariance_matrix(coord,cov_model(m),0.0_dp,sills(m),ranges(m),smooths(m))
            lp_old=gaussian_field_logdensity(z,gpmean(:,m),covs(:,:,m),istat)
            if(istat/=0)then
            result%info=50+m
            return
            end if
         end do

         ! Range and smoothness MH steps.
         do m=1,3
            if(prop_ranges(m)>0.0_dp)then
               propv=exp(log(ranges(m))+prop_ranges(m)*rnorm1())
               proposal_cov=covariance_matrix(coord,cov_model(m),0.0_dp,sills(m),propv,smooths(m))
               z=gp_values(m)
               lp_new=gaussian_field_logdensity(z,gpmean(:,m),proposal_cov,istat)
               if(istat/=0)then
                  result%extreme_rejection(3+m)=result%extreme_rejection(3+m)+1.0_dp
               else
                  lp_old=gaussian_field_logdensity(z,gpmean(:,m),covs(:,:,m),istat)
                  logratio=lp_new-lp_old+(hyper_range(1,m)-1.0_dp)*log(propv/ranges(m)) &
                     +(ranges(m)-propv)/hyper_range(2,m)+log(propv/ranges(m))
                  if(log(max(runif1(),tiny(1.0_dp)))<min(0.0_dp,logratio))then
                     ranges(m)=propv
                     covs(:,:,m)=proposal_cov
                     result%acceptance(3+m)=result%acceptance(3+m)+1.0_dp
                  end if
               end if
            end if
            if(prop_smooths(m)>0.0_dp)then
               propv=exp(log(smooths(m))+prop_smooths(m)*rnorm1())
               proposal_cov=covariance_matrix(coord,cov_model(m),0.0_dp,sills(m),ranges(m),propv)
               z=gp_values(m)
               lp_new=gaussian_field_logdensity(z,gpmean(:,m),proposal_cov,istat)
               if(istat/=0)then
                  result%extreme_rejection(6+m)=result%extreme_rejection(6+m)+1.0_dp
               else
                  lp_old=gaussian_field_logdensity(z,gpmean(:,m),covs(:,:,m),istat)
                  logratio=lp_new-lp_old+(hyper_smooth(1,m)-1.0_dp)*log(propv/smooths(m)) &
                     +(smooths(m)-propv)/hyper_smooth(2,m)+log(propv/smooths(m))
                  if(log(max(runif1(),tiny(1.0_dp)))<min(0.0_dp,logratio))then
                     smooths(m)=propv
                     covs(:,:,m)=proposal_cov
                     result%acceptance(6+m)=result%acceptance(6+m)+1.0_dp
                  end if
               end if
            end if
         end do

         if(iter>iburn .and. mod(iter,ithin)==0)then
            stored=stored+1
            result%chain_loc(stored,:)=[beta_loc,sills(1),ranges(1),smooths(1),gev_start(:,1)]
            result%chain_scale(stored,:)=[beta_scale,sills(2),ranges(2),smooths(2),gev_start(:,2)]
            result%chain_shape(stored,:)=[beta_shape,sills(3),ranges(3),smooths(3),gev_start(:,3)]
         end if
      end do
      result%iterations=iter
      result%acceptance=result%acceptance/real(iter,dp)
      result%extreme_rejection=result%extreme_rejection/real(iter,dp)
      result%acceptance(1:3)=result%acceptance(1:3)/real(nsite,dp)
      result%extreme_rejection(1:3)=result%extreme_rejection(1:3)/real(nsite,dp)
   contains
      function gp_values(margin) result(v)
         integer,intent(in)::margin
         real(dp)::v(nsite)
         v=gev_start(:,margin)
         if(margin==2 .and. log_link)v=log(v)
      end function gp_values

      subroutine update_gpmeans()
         gpmean(:,1)=matmul(design_loc,beta_loc)
         gpmean(:,2)=matmul(design_scale,beta_scale)
         gpmean(:,3)=matmul(design_shape,beta_shape)
      end subroutine update_gpmeans
   end subroutine latent_gev_mcmc

   subroutine update_beta(x,y,cov,prior_mean,prior_prec,beta,info)
      real(dp),intent(in)::x(:,:),y(:),cov(:,:),prior_mean(:),prior_prec(:,:)
      real(dp),intent(out)::beta(:)
      integer,intent(out)::info
      real(dp),allocatable::cinvx(:,:),yrhs(:,:),cinvy(:,:),prec(:,:),rhs(:,:),mean(:,:),postcov(:,:),r(:,:),noise(:)
      integer::p,i
      p=size(x,2)
      info=0
      allocate(cinvx(size(x,1),p),yrhs(size(x,1),1),cinvy(size(x,1),1),prec(p,p),rhs(p,1),mean(p,1),postcov(p,p),r(p,p),noise(p))
      call solve_spd(cov,x,cinvx,info)
      if(info/=0)return
      yrhs(:,1)=y
      call solve_spd(cov,yrhs,cinvy,info)
      if(info/=0)return
      prec=prior_prec+matmul(transpose(x),cinvx)
      rhs(:,1)=matmul(prior_prec,prior_mean)+matmul(transpose(x),cinvy(:,1))
      call solve_spd(prec,rhs,mean,info)
      if(info/=0)return
      call solve_spd(prec,identity_matrix(p),postcov,info)
      if(info/=0)return
      call chol_upper(postcov,r,info)
      if(info/=0)return
      do i=1,p
      noise(i)=rnorm1()
      end do
      beta=mean(:,1)+matmul(transpose(r),noise)
   end subroutine update_beta

   subroutine normalized_quad(res,cov,sill,q,info)
      real(dp),intent(in)::res(:),cov(:,:),sill
      real(dp),intent(out)::q
      integer,intent(out)::info
      real(dp)::rhs(size(res),1),sol(size(res),1)
      rhs(:,1)=res
      call solve_spd(cov,rhs,sol,info)
      if(info/=0)then
      q=huge(1.0_dp)
      return
      end if
      q=sill*dot_product(res,sol(:,1))
   end subroutine normalized_quad

   function identity_matrix(n) result(a)
      integer,intent(in)::n
      real(dp)::a(n,n)
      integer::i
      a=0.0_dp
      do i=1,n
      a(i,i)=1.0_dp
      end do
   end function identity_matrix
end module spatialextremes_latent
