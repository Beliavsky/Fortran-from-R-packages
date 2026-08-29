module spatialextremes_conditional
   use spatialextremes_base, only: dp,chol_upper,inverse_spd,logdet_spd,exp_rand,chisq_rand,pi,eye
   use spatialextremes_mvprob, only: mvnorm_cdf_qmc,mvstudent_cdf_qmc
   use spatialextremes_partitions, only: list_set_partitions,canonicalize_partition,partition_block_count
   use spatialextremes_simulation, only: simulate_brownresnick_exact_hitting
   use r_compat, only: rnorm1,runif1
   implicit none
   private
   type,public :: conditional_maxstable_result_t
      real(dp),allocatable :: sim(:,:),sub_extremal(:,:),extremal(:,:),partition_weights(:)
      integer,allocatable :: partitions(:,:)
      integer :: info=0
   end type conditional_maxstable_result_t
   public :: conditional_schlather_given_partition,conditional_extremalt_given_partition
   public :: conditional_brownresnick_given_partition
   public :: schlather_partition_weights,extremalt_partition_weights,brownresnick_partition_weights
   public :: gibbs_partitions_schlather,gibbs_partitions_extremalt,gibbs_partitions_brownresnick
   public :: sample_conditional_schlather,sample_conditional_extremalt,sample_conditional_brownresnick
   public :: starting_partitions_schlather,starting_partitions_extremalt,starting_partitions_brownresnick
contains
   subroutine conditional_schlather_given_partition(cov,cond_data,partitions,result)
      real(dp),intent(in)::cov(:,:),cond_data(:)
      integer,intent(in)::partitions(:,:)
      type(conditional_maxstable_result_t),intent(out)::result
      integer::n,ncond,nsim,i,b,nb,j,istat,iter
      integer,allocatable::part(:,:)
      real(dp),allocatable::z(:),x(:),prop(:),rchol(:,:),gp(:)
      real(dp)::poisson,ip,thresh,normc
      n=size(cov,1)
      ncond=size(cond_data)
      nsim=size(partitions,1)
      if(size(cov,2)/=n .or. size(partitions,2)/=ncond .or. ncond>n .or. any(cond_data<=0.0_dp))then
         result%info=1
         return
      end if
      call normalize_partition_matrix(partitions,part)
      allocate(result%sim(nsim,n),result%sub_extremal(nsim,n),result%extremal(nsim,n),result%partitions(nsim,ncond))
      result%partitions=part
      allocate(z(n),x(n),prop(n),rchol(n,n),gp(n))
      call chol_upper(cov,rchol,istat)
      if(istat/=0)then
      result%info=2
      return
      end if
      normc=sqrt(2.0_dp*pi)
      do i=1,nsim
         z=-huge(1.0_dp)
         nb=maxval(part(i,:))
         do b=1,nb
            call draw_student_extremal(cov,cond_data,part(i,:),b,1.0_dp,.false.,gp,istat)
            if(istat/=0)then
            result%info=10+istat
            return
            end if
            z=max(z,gp)
         end do
         x=-huge(1.0_dp)
         poisson=0.0_dp
         iter=0
         do
            iter=iter+1
            if(iter>2000000)then
            result%info=30
            return
            end if
            poisson=poisson+exp_rand()
            ip=normc/poisson
            thresh=3.5_dp*2.0_dp*pi*ip
            do j=1,n
            prop(j)=rnorm1()
            end do
            prop=matmul(transpose(rchol),prop)*ip
            if(all(prop(1:ncond)<=cond_data))x=max(x,prop)
            if(all(x>=thresh))exit
         end do
         result%sub_extremal(i,:)=x
         result%extremal(i,:)=z
         result%sim(i,:)=max(x,z)
         result%sim(i,1:ncond)=cond_data
      end do
      result%info=0
   end subroutine conditional_schlather_given_partition

   subroutine conditional_extremalt_given_partition(cov,cond_data,nu,partitions,result)
      real(dp),intent(in)::cov(:,:),cond_data(:),nu
      integer,intent(in)::partitions(:,:)
      type(conditional_maxstable_result_t),intent(out)::result
      integer::n,ncond,nsim,i,b,nb,j,istat,iter
      integer,allocatable::part(:,:)
      real(dp),allocatable::z(:),x(:),prop(:),rchol(:,:),gp(:)
      real(dp)::poisson,ip,thresh,cnu
      n=size(cov,1)
      ncond=size(cond_data)
      nsim=size(partitions,1)
      if(nu<=0.0_dp .or. size(cov,2)/=n .or. size(partitions,2)/=ncond .or. ncond>n .or. any(cond_data<=0.0_dp))then
         result%info=1
         return
      end if
      call normalize_partition_matrix(partitions,part)
      allocate(result%sim(nsim,n),result%sub_extremal(nsim,n),result%extremal(nsim,n),result%partitions(nsim,ncond))
      result%partitions=part
      allocate(z(n),x(n),prop(n),rchol(n,n),gp(n))
      call chol_upper(cov,rchol,istat)
      if(istat/=0)then
      result%info=2
      return
      end if
      cnu=sqrt(pi)*2.0_dp**(1.0_dp-0.5_dp*nu)/gamma(0.5_dp*(nu+1.0_dp))
      do i=1,nsim
         z=-huge(1.0_dp)
         nb=maxval(part(i,:))
         do b=1,nb
            call draw_student_extremal(cov,cond_data,part(i,:),b,nu,.true.,gp,istat)
            if(istat/=0)then
            result%info=10+istat
            return
            end if
            z=max(z,gp)
         end do
         x=-huge(1.0_dp)
         poisson=0.0_dp
         iter=0
         do
            iter=iter+1
            if(iter>2000000)then
            result%info=30
            return
            end if
            poisson=poisson+exp_rand()
            ip=cnu/poisson
            thresh=3.5_dp*cnu*ip
            do j=1,n
            prop(j)=rnorm1()
            end do
            prop=matmul(transpose(rchol),prop)
            prop=ip*max(prop,0.0_dp)**nu
            if(all(prop(1:ncond)<=cond_data))x=max(x,prop)
            if(all(x>=thresh))exit
         end do
         result%sub_extremal(i,:)=x
         result%extremal(i,:)=z
         result%sim(i,:)=max(x,z)
         result%sim(i,1:ncond)=cond_data
      end do
      result%info=0
   end subroutine conditional_extremalt_given_partition

   subroutine conditional_brownresnick_given_partition(coord,cond_data,range,smooth,partitions,result,n_subextremal)
      real(dp),intent(in)::coord(:,:),cond_data(:),range,smooth
      integer,intent(in)::partitions(:,:)
      type(conditional_maxstable_result_t),intent(out)::result
      integer,intent(in),optional::n_subextremal
      integer::n,ncond,nsim,i,b,nb,j,istat,nsub,accepted,iter
      integer,allocatable::part(:,:)
      real(dp),allocatable::c(:,:),sigma2(:),h(:,:),mbar(:),cc(:,:),z(:),x(:),prop(:),rchol(:,:),eps(:),yt(:)
      real(dp)::poisson
      n=size(coord,1)
      ncond=size(cond_data)
      nsim=size(partitions,1)
      nsub=500
      if(present(n_subextremal))nsub=max(1,n_subextremal)
      if(size(partitions,2)/=ncond .or. ncond>n .or. size(coord,2)<1 .or. range<=0.0_dp .or. &
         smooth<=0.0_dp .or. smooth>2.0_dp .or. any(cond_data<=0.0_dp))then
         result%info=1
         return
         end if
      call brown_setup(coord,range,smooth,c,sigma2,h,mbar,cc,istat)
      if(istat/=0)then
      result%info=2
      return
      end if
      call normalize_partition_matrix(partitions,part)
      allocate(result%sim(nsim,n),result%sub_extremal(nsim,n),result%extremal(nsim,n),result%partitions(nsim,ncond))
      result%partitions=part
      allocate(z(n),x(n),prop(n),rchol(n,n),eps(n),yt(n))
      yt=0.0_dp
      yt(1:ncond)=log(cond_data)
      call chol_upper(c,rchol,istat)
      if(istat/=0)then
      result%info=3
      return
      end if
      do i=1,nsim
         z=-huge(1.0_dp)
         nb=maxval(part(i,:))
         do b=1,nb
            call draw_brown_extremal(h,mbar,yt,ncond,part(i,:),b,eps,istat)
            if(istat/=0)then
            result%info=10+istat
            return
            end if
            z=max(z,eps)
         end do
         x=-huge(1.0_dp)
         poisson=0.0_dp
         accepted=0
         iter=0
         do while(accepted<nsub)
            iter=iter+1
            if(iter>10000000)then
            result%info=30
            return
            end if
            poisson=poisson+exp_rand()
            do j=1,n
            prop(j)=rnorm1()
            end do
            prop=matmul(transpose(rchol),prop)-0.5_dp*sigma2-log(poisson)
            if(all(prop(1:ncond)<=yt(1:ncond)))then
            x=max(x,prop)
            accepted=accepted+1
            end if
         end do
         result%sub_extremal(i,:)=exp(x)
         result%extremal(i,:)=exp(z)
         result%sim(i,:)=exp(max(x,z))
         result%sim(i,1:ncond)=cond_data
      end do
      result%info=0
   end subroutine conditional_brownresnick_given_partition

   subroutine schlather_partition_weights(cov,y,partitions,weights,info,n_per_dim)
      real(dp),intent(in)::cov(:,:),y(:)
      integer,allocatable,intent(out)::partitions(:,:)
      real(dp),allocatable,intent(out)::weights(:)
      integer,intent(out)::info
      integer,intent(in),optional::n_per_dim
      integer,allocatable::raw(:,:),bc(:)
      integer::p,npart
      real(dp),allocatable::lw(:)
      if(size(cov,1)/=size(y) .or. size(cov,2)/=size(y) .or. any(y<=0.0_dp))then
      info=1
      return
      end if
      call list_set_partitions(size(y),raw,bc,info)
      if(info/=0)return
      npart=size(raw,2)
      allocate(partitions(npart,size(y)),weights(npart),lw(npart))
      do p=1,npart
         partitions(p,:)=raw(:,p)+1
         lw(p)=partition_logweight_student(cov,y,partitions(p,:),1.0_dp,.false.,n_per_dim,info)
         if(info/=0)return
      end do
      call normalize_logweights(lw,weights)
      info=0
   end subroutine schlather_partition_weights

   subroutine extremalt_partition_weights(cov,y,nu,partitions,weights,info,n_per_dim)
      real(dp),intent(in)::cov(:,:),y(:),nu
      integer,allocatable,intent(out)::partitions(:,:)
      real(dp),allocatable,intent(out)::weights(:)
      integer,intent(out)::info
      integer,intent(in),optional::n_per_dim
      integer,allocatable::raw(:,:),bc(:)
      integer::p,npart
      real(dp),allocatable::lw(:)
      if(nu<=0.0_dp .or. size(cov,1)/=size(y) .or. size(cov,2)/=size(y) .or. any(y<=0.0_dp))then
      info=1
      return
      end if
      call list_set_partitions(size(y),raw,bc,info)
      if(info/=0)return
      npart=size(raw,2)
      allocate(partitions(npart,size(y)),weights(npart),lw(npart))
      do p=1,npart
         partitions(p,:)=raw(:,p)+1
         lw(p)=partition_logweight_student(cov,y,partitions(p,:),nu,.true.,n_per_dim,info)
         if(info/=0)return
      end do
      call normalize_logweights(lw,weights)
      info=0
   end subroutine extremalt_partition_weights

   subroutine brownresnick_partition_weights(coord,y,range,smooth,partitions,weights,info,n_per_dim)
      real(dp),intent(in)::coord(:,:),y(:),range,smooth
      integer,allocatable,intent(out)::partitions(:,:)
      real(dp),allocatable,intent(out)::weights(:)
      integer,intent(out)::info
      integer,intent(in),optional::n_per_dim
      integer,allocatable::raw(:,:),bc(:)
      integer::p,npart
      real(dp),allocatable::c(:,:),sigma2(:),h(:,:),mbar(:),cc(:,:),lw(:),ly(:)
      if(size(coord,1)/=size(y) .or. range<=0.0_dp .or. smooth<=0.0_dp .or. smooth>2.0_dp .or. any(y<=0.0_dp))then
         info=1
         return
      end if
      call brown_setup(coord,range,smooth,c,sigma2,h,mbar,cc,info)
      if(info/=0)return
      call list_set_partitions(size(y),raw,bc,info)
      if(info/=0)return
      npart=size(raw,2)
      allocate(partitions(npart,size(y)),weights(npart),lw(npart),ly(size(y)))
      ly=log(y)
      do p=1,npart
         partitions(p,:)=raw(:,p)+1
         lw(p)=partition_logweight_brown(c,sigma2,h,mbar,ly,partitions(p,:),n_per_dim,info)
         if(info/=0)return
      end do
      call normalize_logweights(lw,weights)
      info=0
   end subroutine brownresnick_partition_weights

   subroutine gibbs_partitions_schlather(cov,y,nchain,thin,burnin,chain,info,start,n_per_dim)
      real(dp),intent(in)::cov(:,:),y(:)
      integer,intent(in)::nchain,thin,burnin
      integer,allocatable,intent(out)::chain(:,:)
      integer,intent(out)::info
      integer,intent(in),optional::start(:),n_per_dim
      integer::n,iter,kept,idx,b,ncand,k,sel
      integer,allocatable::part(:),cand(:,:),cp(:)
      real(dp),allocatable::lw(:),pr(:)
      n=size(y)
      if(n<1 .or. nchain<1 .or. thin<1 .or. size(cov,1)/=n .or. size(cov,2)/=n)then
      info=1
      return
      end if
      allocate(part(n),chain(nchain,n))
      call initial_partition(n,start,part,info)
      if(info/=0)return
      iter=0
      kept=0
      do while(kept<nchain)
         idx=1+int(runif1()*real(n,dp))
         idx=min(n,max(1,idx))
         b=partition_block_count(part-1)
         call build_gibbs_candidates(part,idx,cand)
         ncand=size(cand,1)
         allocate(lw(ncand),pr(ncand))
         do k=1,ncand
            lw(k)=partition_logweight_student(cov,y,cand(k,:),1.0_dp,.false.,n_per_dim,info)
            if(info/=0)return
         end do
         call normalize_logweights(lw,pr)
         sel=sample_discrete(pr)
         part=cand(sel,:)
         iter=iter+1
         if(iter>burnin .and. mod(iter,thin)==0)then
         kept=kept+1
         chain(kept,:)=part
         end if
         deallocate(lw,pr,cand)
      end do
      info=0
   end subroutine gibbs_partitions_schlather

   subroutine gibbs_partitions_extremalt(cov,y,nu,nchain,thin,burnin,chain,info,start,n_per_dim)
      real(dp),intent(in)::cov(:,:),y(:),nu
      integer,intent(in)::nchain,thin,burnin
      integer,allocatable,intent(out)::chain(:,:)
      integer,intent(out)::info
      integer,intent(in),optional::start(:),n_per_dim
      integer::n,iter,kept,idx,ncand,k,sel
      integer,allocatable::part(:),cand(:,:)
      real(dp),allocatable::lw(:),pr(:)
      n=size(y)
      if(nu<=0.0_dp .or. n<1 .or. nchain<1 .or. thin<1 .or. &
         size(cov,1)/=n .or. size(cov,2)/=n)then
         info=1
         return
      end if
      allocate(part(n),chain(nchain,n))
      call initial_partition(n,start,part,info)
      if(info/=0)return
      iter=0
      kept=0
      do while(kept<nchain)
         idx=min(n,max(1,1+int(runif1()*real(n,dp))))
         call build_gibbs_candidates(part,idx,cand)
         ncand=size(cand,1)
         allocate(lw(ncand),pr(ncand))
         do k=1,ncand
            lw(k)=partition_logweight_student(cov,y,cand(k,:),nu,.true.,n_per_dim,info)
            if(info/=0)return
         end do
         call normalize_logweights(lw,pr)
         sel=sample_discrete(pr)
         part=cand(sel,:)
         iter=iter+1
         if(iter>burnin .and. mod(iter,thin)==0)then
         kept=kept+1
         chain(kept,:)=part
         end if
         deallocate(lw,pr,cand)
      end do
      info=0
   end subroutine gibbs_partitions_extremalt

   subroutine gibbs_partitions_brownresnick(coord,y,range,smooth,nchain,thin,burnin,chain,info,start,n_per_dim)
      real(dp),intent(in)::coord(:,:),y(:),range,smooth
      integer,intent(in)::nchain,thin,burnin
      integer,allocatable,intent(out)::chain(:,:)
      integer,intent(out)::info
      integer,intent(in),optional::start(:),n_per_dim
      integer::n,iter,kept,idx,ncand,k,sel
      integer,allocatable::part(:),cand(:,:)
      real(dp),allocatable::c(:,:),sigma2(:),h(:,:),mbar(:),cc(:,:),ly(:),lw(:),pr(:)
      n=size(y)
      if(size(coord,1)/=n .or. nchain<1 .or. thin<1 .or. range<=0.0_dp)then
      info=1
      return
      end if
      call brown_setup(coord,range,smooth,c,sigma2,h,mbar,cc,info)
      if(info/=0)return
      allocate(part(n),chain(nchain,n),ly(n))
      ly=log(y)
      call initial_partition(n,start,part,info)
      if(info/=0)return
      iter=0
      kept=0
      do while(kept<nchain)
         idx=min(n,max(1,1+int(runif1()*real(n,dp))))
         call build_gibbs_candidates(part,idx,cand)
         ncand=size(cand,1)
         allocate(lw(ncand),pr(ncand))
         do k=1,ncand
            lw(k)=partition_logweight_brown(c,sigma2,h,mbar,ly,cand(k,:),n_per_dim,info)
            if(info/=0)return
         end do
         call normalize_logweights(lw,pr)
         sel=sample_discrete(pr)
         part=cand(sel,:)
         iter=iter+1
         if(iter>burnin .and. mod(iter,thin)==0)then
         kept=kept+1
         chain(kept,:)=part
         end if
         deallocate(lw,pr,cand)
      end do
      info=0
   end subroutine gibbs_partitions_brownresnick

   subroutine starting_partitions_schlather(cov,nsim,parts,info)
      real(dp),intent(in)::cov(:,:)
      integer,intent(in)::nsim
      integer,allocatable,intent(out)::parts(:,:)
      integer,intent(out)::info
      real(dp),allocatable::r(:,:),x(:),prop(:)
      real(dp)::poisson,ip,thresh,normc
      integer,allocatable::labels(:)
      integer::n,i,j,storm,iter
      logical::changed
      n=size(cov,1)
      if(n<1 .or. size(cov,2)/=n .or. nsim<1)then
      info=1
      return
      end if
      allocate(parts(nsim,n),r(n,n),x(n),prop(n),labels(n))
      call chol_upper(cov,r,info)
      if(info/=0)return
      normc=sqrt(2.0_dp*pi)
      do i=1,nsim
         x=0.0_dp
         labels=0
         poisson=0.0_dp
         storm=0
         iter=0
         do
            iter=iter+1
            if(iter>2000000)then
            info=2
            return
            end if
            poisson=poisson+exp_rand()
            ip=1.0_dp/poisson
            thresh=3.5_dp*normc*ip
            do j=1,n
            prop(j)=rnorm1()
            end do
            prop=max(0.0_dp,normc*ip*matmul(transpose(r),prop))
            changed=.false.
            do j=1,n
               if(prop(j)>x(j))then
               labels(j)=storm+1
               changed=.true.
               end if
               x(j)=max(x(j),prop(j))
            end do
            if(changed)storm=storm+1
            if(all(x>=thresh))exit
         end do
         parts(i,:)=canonicalize_partition(labels)+1
      end do
      info=0
   end subroutine starting_partitions_schlather

   subroutine starting_partitions_extremalt(cov,nu,nsim,parts,info)
      real(dp),intent(in)::cov(:,:),nu
      integer,intent(in)::nsim
      integer,allocatable,intent(out)::parts(:,:)
      integer,intent(out)::info
      real(dp),allocatable::r(:,:),x(:),prop(:)
      real(dp)::poisson,ip,thresh,cnu
      integer,allocatable::labels(:)
      integer::n,i,j,storm,iter
      logical::changed
      n=size(cov,1)
      if(n<1 .or. size(cov,2)/=n .or. nsim<1 .or. nu<=0.0_dp)then
      info=1
      return
      end if
      allocate(parts(nsim,n),r(n,n),x(n),prop(n),labels(n))
      call chol_upper(cov,r,info)
      if(info/=0)return
      cnu=sqrt(pi)*2.0_dp**(1.0_dp-0.5_dp*nu)/gamma(0.5_dp*(nu+1.0_dp))
      do i=1,nsim
         x=0.0_dp
         labels=0
         poisson=0.0_dp
         storm=0
         iter=0
         do
            iter=iter+1
            if(iter>2000000)then
            info=2
            return
            end if
            poisson=poisson+exp_rand()
            ip=1.0_dp/poisson
            thresh=3.5_dp*cnu*ip
            do j=1,n
            prop(j)=rnorm1()
            end do
            prop=cnu*ip*max(0.0_dp,matmul(transpose(r),prop))**nu
            changed=.false.
            do j=1,n
               if(prop(j)>x(j))then
               labels(j)=storm+1
               changed=.true.
               end if
               x(j)=max(x(j),prop(j))
            end do
            if(changed)storm=storm+1
            if(all(x>=thresh))exit
         end do
         parts(i,:)=canonicalize_partition(labels)+1
      end do
      info=0
   end subroutine starting_partitions_extremalt

   subroutine starting_partitions_brownresnick(coord,range,smooth,nsim,parts,info)
      real(dp),intent(in)::coord(:,:),range,smooth
      integer,intent(in)::nsim
      integer,allocatable,intent(out)::parts(:,:)
      integer,intent(out)::info
      real(dp),allocatable::x(:,:),centered(:,:)
      integer,allocatable::hit(:,:)
      integer::i,k
      allocate(centered(size(coord,1),size(coord,2)))
      do k=1,size(coord,2)
      centered(:,k)=coord(:,k)-sum(coord(:,k))/real(size(coord,1),dp)
      end do
      call simulate_brownresnick_exact_hitting(nsim,centered,range,smooth,x,hit,info)
      if(info/=0)return
      allocate(parts(nsim,size(coord,1)))
      do i=1,nsim
      parts(i,:)=canonicalize_partition(hit(i,:))+1
      end do
   end subroutine starting_partitions_brownresnick

   subroutine sample_conditional_schlather(cov,y,nsim,result,thin,burnin,n_per_dim)
      real(dp),intent(in)::cov(:,:),y(:)
      integer,intent(in)::nsim
      type(conditional_maxstable_result_t),intent(out)::result
      integer,intent(in),optional::thin,burnin,n_per_dim
      integer::info,t,b,i
      integer,allocatable::allp(:,:),parts(:,:),startpool(:,:),startp(:)
      real(dp),allocatable::w(:)
      t=size(y)
      if(present(thin))t=thin
      b=50
      if(present(burnin))b=burnin
      if(size(y)<=7)then
         call schlather_partition_weights(cov(1:size(y),1:size(y)),y,allp,w,info,n_per_dim)
         if(info/=0)then
         result%info=info
         return
         end if
         allocate(parts(nsim,size(y)))
         do i=1,nsim
         parts(i,:)=allp(sample_discrete(w),:)
         end do
         call conditional_schlather_given_partition(cov,y,parts,result)
         result%partition_weights=w
      else
         call starting_partitions_schlather(cov(1:size(y),1:size(y)),250,startpool,info)
         if(info/=0)then
         result%info=info
         return
         end if
         call modal_partition(startpool,startp)
         call gibbs_partitions_schlather(cov(1:size(y),1:size(y)),y,nsim,t,b,parts,info, &
            start=startp,n_per_dim=n_per_dim)
         if(info/=0)then
         result%info=info
         return
         end if
         call conditional_schlather_given_partition(cov,y,parts,result)
      end if
   end subroutine sample_conditional_schlather

   subroutine sample_conditional_extremalt(cov,y,nu,nsim,result,thin,burnin,n_per_dim)
      real(dp),intent(in)::cov(:,:),y(:),nu
      integer,intent(in)::nsim
      type(conditional_maxstable_result_t),intent(out)::result
      integer,intent(in),optional::thin,burnin,n_per_dim
      integer::info,t,b,i
      integer,allocatable::allp(:,:),parts(:,:),startpool(:,:),startp(:)
      real(dp),allocatable::w(:)
      t=size(y)
      if(present(thin))t=thin
      b=50
      if(present(burnin))b=burnin
      if(size(y)<=7)then
         call extremalt_partition_weights(cov(1:size(y),1:size(y)),y,nu,allp,w,info,n_per_dim)
         if(info/=0)then
         result%info=info
         return
         end if
         allocate(parts(nsim,size(y)))
         do i=1,nsim
         parts(i,:)=allp(sample_discrete(w),:)
         end do
         call conditional_extremalt_given_partition(cov,y,nu,parts,result)
         result%partition_weights=w
      else
         call starting_partitions_extremalt(cov(1:size(y),1:size(y)),nu,250,startpool,info)
         if(info/=0)then
         result%info=info
         return
         end if
         call modal_partition(startpool,startp)
         call gibbs_partitions_extremalt(cov(1:size(y),1:size(y)),y,nu,nsim,t,b,parts,info, &
            start=startp,n_per_dim=n_per_dim)
         if(info/=0)then
         result%info=info
         return
         end if
         call conditional_extremalt_given_partition(cov,y,nu,parts,result)
      end if
   end subroutine sample_conditional_extremalt

   subroutine sample_conditional_brownresnick(coord,y,range,smooth,nsim,result,thin,burnin,n_per_dim,n_subextremal)
      real(dp),intent(in)::coord(:,:),y(:),range,smooth
      integer,intent(in)::nsim
      type(conditional_maxstable_result_t),intent(out)::result
      integer,intent(in),optional::thin,burnin,n_per_dim,n_subextremal
      integer::info,t,b,i
      integer,allocatable::allp(:,:),parts(:,:),startpool(:,:),startp(:)
      real(dp),allocatable::w(:)
      t=size(y)
      if(present(thin))t=thin
      b=50
      if(present(burnin))b=burnin
      if(size(y)<=7)then
         call brownresnick_partition_weights(coord(1:size(y),:),y,range,smooth,allp,w,info,n_per_dim)
         if(info/=0)then
         result%info=info
         return
         end if
         allocate(parts(nsim,size(y)))
         do i=1,nsim
         parts(i,:)=allp(sample_discrete(w),:)
         end do
         call conditional_brownresnick_given_partition(coord,y,range,smooth,parts,result,n_subextremal)
         result%partition_weights=w
      else
         call starting_partitions_brownresnick(coord(1:size(y),:),range,smooth,250,startpool,info)
         if(info/=0)then
         result%info=info
         return
         end if
         call modal_partition(startpool,startp)
         call gibbs_partitions_brownresnick(coord(1:size(y),:),y,range,smooth,nsim,t,b,parts,info, &
            start=startp,n_per_dim=n_per_dim)
         if(info/=0)then
         result%info=info
         return
         end if
         call conditional_brownresnick_given_partition(coord,y,range,smooth,parts,result,n_subextremal)
      end if
   end subroutine sample_conditional_brownresnick

   subroutine draw_student_extremal(cov,y,partition,block,nu,is_extt,out,info)
      real(dp),intent(in)::cov(:,:),y(:),nu
      integer,intent(in)::partition(:),block
      logical,intent(in)::is_extt
      real(dp),intent(out)::out(:)
      integer,intent(out)::info
      integer,allocatable::tau(:),bar(:)
      real(dp),allocatable::mu(:),scale(:,:),r(:,:),eps(:)
      real(dp)::dof,factor
      integer::k,nbar,ncond,rn
      call block_indices(partition,block,size(cov,1),tau,bar)
      rn=size(tau)
      nbar=size(bar)
      ncond=size(y)
      call student_parameters(cov,y,tau,bar,nu,is_extt,dof,mu,scale,info)
      if(info/=0)return
      allocate(r(nbar,nbar),eps(nbar))
      call chol_upper(scale+1.0e-13_dp*eye(nbar),r,info)
      if(info/=0)return
      out=-huge(1.0_dp)
      out(tau)=y(tau)
      do
         do k=1,nbar
         eps(k)=rnorm1()
         end do
         factor=sqrt(dof/chisq_rand(dof))
         eps=mu+factor*matmul(transpose(r),eps)
         if(is_extt)eps=max(eps,0.0_dp)**nu
         if(ncond-rn==0 .or. all(eps(1:ncond-rn)<=y(bar(1:ncond-rn))))exit
      end do
      out(bar)=eps
      info=0
   end subroutine draw_student_extremal

   subroutine student_parameters(cov,y,tau,bar,nu,is_extt,dof,mu,scale,info)
      real(dp),intent(in)::cov(:,:),y(:),nu
      integer,intent(in)::tau(:),bar(:)
      logical,intent(in)::is_extt
      real(dp),intent(out)::dof
      real(dp),allocatable,intent(out)::mu(:),scale(:,:)
      integer,intent(out)::info
      real(dp),allocatable::cx(:,:),ci(:,:),cs(:,:),csx(:,:),yx(:),tmp(:,:)
      real(dp)::mahal
      integer::r,nr
      r=size(tau)
      nr=size(bar)
      allocate(cx(r,r),ci(r,r),cs(nr,nr),csx(nr,r),yx(r),mu(nr),scale(nr,nr))
      cx=cov(tau,tau)
      cs=cov(bar,bar)
      csx=cov(bar,tau)
      call inverse_spd(cx,ci,info)
      if(info/=0)return
      if(is_extt)then
      yx=y(tau)**(1.0_dp/nu)
      dof=real(r,dp)+nu
      else
      yx=y(tau)
      dof=real(r+1,dp)
      end if
      mu=matmul(csx,matmul(ci,yx))
      mahal=dot_product(yx,matmul(ci,yx))
      scale=(mahal/dof)*(cs-matmul(csx,matmul(ci,transpose(csx))))
      scale=0.5_dp*(scale+transpose(scale))
      info=0
   end subroutine student_parameters

   function partition_logweight_student(cov,y,part,nu,is_extt,n_per_dim,info) result(lw)
      real(dp),intent(in)::cov(:,:),y(:),nu
      integer,intent(in)::part(:)
      logical,intent(in)::is_extt
      integer,intent(in),optional::n_per_dim
      integer,intent(out)::info
      real(dp)::lw,one
      integer::b,nb
      lw=0.0_dp
      nb=maxval(part)
      do b=1,nb
         one=block_logweight_student(cov,y,part,b,nu,is_extt,n_per_dim,info)
         if(info/=0)return
         lw=lw+one
      end do
   end function partition_logweight_student

   function block_logweight_student(cov,y,part,block,nu,is_extt,n_per_dim,info) result(lw)
      real(dp),intent(in)::cov(:,:),y(:),nu
      integer,intent(in)::part(:),block
      logical,intent(in)::is_extt
      integer,intent(in),optional::n_per_dim
      integer,intent(out)::info
      real(dp)::lw,dof,p,f
      integer,allocatable::tau(:),bar(:)
      real(dp),allocatable::mu(:),scale(:,:),bounds(:)
      call block_indices(part,block,size(y),tau,bar)
      if(size(tau)==0)then
      lw=-huge(1.0_dp)
      info=1
      return
      end if
      if(size(bar)>0)then
         call student_parameters(cov,y,tau,bar,nu,is_extt,dof,mu,scale,info)
         if(info/=0)return
         allocate(bounds(size(bar)))
         if(is_extt)then
         bounds=y(bar)**(1.0_dp/nu)
         else
         bounds=y(bar)
         end if
         p=mvstudent_cdf_qmc(bounds,dof,mu,scale,n_per_dim,info)
         if(info/=0)return
         if(p<=tiny(1.0_dp))then
         lw=log(tiny(1.0_dp))
         else
         lw=log(p)
         end if
      else
         lw=0.0_dp
      end if
      if(is_extt)then
      f=extt_intensity_log(cov,y,tau,nu,info)
      else
      f=schlather_intensity_log(cov,y,tau,info)
      end if
      if(info/=0)return
      lw=lw+f
      info=0
   end function block_logweight_student

   function schlather_intensity_log(cov,y,tau,info) result(f)
      real(dp),intent(in)::cov(:,:),y(:)
      integer,intent(in)::tau(:)
      integer,intent(out)::info
      real(dp)::f,ld,mahal
      real(dp),allocatable::a(:,:),ai(:,:),v(:)
      integer::r
      r=size(tau)
      allocate(a(r,r),ai(r,r),v(r))
      a=cov(tau,tau)
      v=y(tau)
      call inverse_spd(a,ai,info)
      if(info/=0)return
      ld=logdet_spd(a,info)
      if(info/=0)return
      mahal=dot_product(v,matmul(ai,v))
      f=0.5_dp*real(1-r,dp)*log(pi)-0.5_dp*ld-0.5_dp*real(r+1,dp)*log(mahal)+log_gamma(0.5_dp*real(r+1,dp))
   end function schlather_intensity_log

   function extt_intensity_log(cov,y,tau,nu,info) result(f)
      real(dp),intent(in)::cov(:,:),y(:),nu
      integer,intent(in)::tau(:)
      integer,intent(out)::info
      real(dp)::f,ld,mahal,logc
      real(dp),allocatable::a(:,:),ai(:,:),v(:),yy(:)
      integer::r
      r=size(tau)
      allocate(a(r,r),ai(r,r),v(r),yy(r))
      a=cov(tau,tau)
      yy=y(tau)
      v=yy**(1.0_dp/nu)
      call inverse_spd(a,ai,info)
      if(info/=0)return
      ld=logdet_spd(a,info)
      if(info/=0)return
      mahal=dot_product(v,matmul(ai,v))
      logc=(1.0_dp-0.5_dp*nu)*log(2.0_dp)+0.5_dp*log(pi)-log_gamma(0.5_dp*(nu+1.0_dp))
      f=sum((1.0_dp/nu-1.0_dp)*log(yy))+logc+real(1-r,dp)*log(nu)+(2.0_dp-nu)*log(2.0_dp) &
         -0.5_dp*real(r,dp)*log(pi)-0.5_dp*ld-0.5_dp*(real(r,dp)+nu)*log(mahal)+log_gamma(0.5_dp*(real(r,dp)+nu))
   end function extt_intensity_log

   subroutine brown_setup(coord,range,smooth,cov,sigma2,h,mbar,centered,info)
      real(dp),intent(in)::coord(:,:),range,smooth
      real(dp),allocatable,intent(out)::cov(:,:),sigma2(:),h(:,:),mbar(:),centered(:,:)
      integer,intent(out)::info
      real(dp),allocatable::ci(:,:),v(:),w(:),one(:),tmp(:,:)
      real(dp)::den,d
      integer::n,i,j,k
      n=size(coord,1)
      allocate(centered(n,size(coord,2)),cov(n,n),sigma2(n),h(n,n),mbar(n),ci(n,n),v(n),w(n),one(n))
      do k=1,size(coord,2)
      centered(:,k)=coord(:,k)-sum(coord(:,k))/real(n,dp)
      end do
      do i=1,n
      sigma2(i)=2.0_dp*(sqrt(sum(centered(i,:)**2))/range)**smooth
      end do
      do i=1,n
      do j=i,n
         d=sqrt(sum((centered(i,:)-centered(j,:))**2))
         cov(i,j)=0.5_dp*(sigma2(i)+sigma2(j))-(d/range)**smooth
         cov(j,i)=cov(i,j)
      end do
      end do
      call inverse_spd(cov,ci,info)
      if(info/=0)then
         cov=cov+1.0e-11_dp*max(1.0_dp,maxval(abs(cov)))*eye(n)
         call inverse_spd(cov,ci,info)
         if(info/=0)return
      end if
      one=1.0_dp
      v=matmul(ci,one)
      w=matmul(ci,sigma2)
      den=sum(v)
      h=ci-outer_product(v,v)/den
      mbar=-0.5_dp*w+(0.5_dp*sum(w)-1.0_dp)*v/den
      info=0
   end subroutine brown_setup

   function partition_logweight_brown(cov,sigma2,h,mbar,y,part,n_per_dim,info) result(lw)
      real(dp),intent(in)::cov(:,:),sigma2(:),h(:,:),mbar(:),y(:)
      integer,intent(in)::part(:)
      integer,intent(in),optional::n_per_dim
      integer,intent(out)::info
      real(dp)::lw,x
      integer::b
      lw=0.0_dp
      do b=1,maxval(part)
         x=block_logweight_brown(cov,sigma2,h,mbar,y,part,b,n_per_dim,info)
         if(info/=0)return
         lw=lw+x
      end do
   end function partition_logweight_brown

   function block_logweight_brown(cov,sigma2,h,mbar,y,part,block,n_per_dim,info) result(lw)
      real(dp),intent(in)::cov(:,:),sigma2(:),h(:,:),mbar(:),y(:)
      integer,intent(in)::part(:),block
      integer,intent(in),optional::n_per_dim
      integer,intent(out)::info
      real(dp)::lw,p,f
      integer,allocatable::tau(:),bar(:)
      real(dp),allocatable::bmat(:,:),mu(:),rhs(:)
      call block_indices(part,block,size(y),tau,bar)
      if(size(bar)>0)then
         call brown_parameters(h,mbar,y,tau,bar,bmat,mu,info)
         if(info/=0)return
         p=mvnorm_cdf_qmc(y(bar),bmat,mu,n_per_dim,info)
         if(info/=0)return
         lw=log(max(p,tiny(1.0_dp)))
      else
      lw=0.0_dp
      end if
      f=brown_intensity_log(cov,sigma2,y,tau,info)
      if(info/=0)return
      lw=lw+f
   end function block_logweight_brown

   function brown_intensity_log(cov,sigma2,y,tau,info) result(f)
      real(dp),intent(in)::cov(:,:),sigma2(:),y(:)
      integer,intent(in)::tau(:)
      integer,intent(out)::info
      real(dp)::f,ld,m1,mm,m1m,den
      real(dp),allocatable::a(:,:),ai(:,:),one(:),mean(:)
      integer::r
      r=size(tau)
      allocate(a(r,r),ai(r,r),one(r),mean(r))
      a=cov(tau,tau)
      one=1.0_dp
      mean=y(tau)+0.5_dp*sigma2(tau)
      call inverse_spd(a,ai,info)
      if(info/=0)return
      ld=logdet_spd(a,info)
      if(info/=0)return
      m1=dot_product(one,matmul(ai,one))
      mm=dot_product(mean,matmul(ai,mean))
      m1m=dot_product(one,matmul(ai,mean))
      f=0.5_dp*real(1-r,dp)*log(2.0_dp*pi)-0.5_dp*(ld+log(m1)+mm-(m1m-1.0_dp)**2/m1)-sum(y(tau))
   end function brown_intensity_log

   subroutine brown_parameters(h,mbar,y,tau,bar,bmat,mu,info)
      real(dp),intent(in)::h(:,:),mbar(:),y(:)
      integer,intent(in)::tau(:),bar(:)
      real(dp),allocatable,intent(out)::bmat(:,:),mu(:)
      integer,intent(out)::info
      real(dp),allocatable::ib(:,:),rhs(:)
      integer::nr
      nr=size(bar)
      allocate(ib(nr,nr),bmat(nr,nr),mu(nr),rhs(nr))
      ib=h(bar,bar)
      call inverse_spd(ib,bmat,info)
      if(info/=0)return
      rhs=mbar(bar)
      if(size(tau)>0)rhs=rhs-matmul(h(bar,tau),y(tau))
      mu=matmul(bmat,rhs)
      info=0
   end subroutine brown_parameters

   subroutine draw_brown_extremal(h,mbar,y,ncond,part,block,out,info)
      real(dp),intent(in)::h(:,:),mbar(:),y(:)
      integer,intent(in)::ncond,part(:),block
      real(dp),intent(out)::out(:)
      integer,intent(out)::info
      integer,allocatable::tau(:),bar(:)
      real(dp),allocatable::bmat(:,:),mu(:),r(:,:),eps(:)
      integer::k,nr,rn
      call block_indices(part,block,size(y),tau,bar)
      rn=size(tau)
      nr=size(bar)
      call brown_parameters(h,mbar,y,tau,bar,bmat,mu,info)
      if(info/=0)return
      allocate(r(nr,nr),eps(nr))
      call chol_upper(bmat+1.0e-13_dp*eye(nr),r,info)
      if(info/=0)return
      out=-huge(1.0_dp)
      out(tau)=y(tau)
      do
         do k=1,nr
         eps(k)=rnorm1()
         end do
         eps=mu+matmul(transpose(r),eps)
         if(ncond-rn==0 .or. all(eps(1:ncond-rn)<=y(bar(1:ncond-rn))))exit
      end do
      out(bar)=eps
      info=0
   end subroutine draw_brown_extremal

   subroutine block_indices(part,block,n_total,tau,bar)
      integer,intent(in)::part(:),block,n_total
      integer,allocatable,intent(out)::tau(:),bar(:)
      integer::r,nr,i,k
      r=count(part==block)
      nr=n_total-r
      allocate(tau(r),bar(nr))
      k=0
      do i=1,size(part)
      if(part(i)==block)then
      k=k+1
      tau(k)=i
      end if
      end do
      k=0
      do i=1,size(part)
      if(part(i)/=block)then
      k=k+1
      bar(k)=i
      end if
      end do
      do i=size(part)+1,n_total
      k=k+1
      bar(k)=i
      end do
   end subroutine block_indices

   subroutine normalize_partition_matrix(pin,pout)
      integer,intent(in)::pin(:,:)
      integer,allocatable,intent(out)::pout(:,:)
      integer::i
      allocate(pout(size(pin,1),size(pin,2)))
      do i=1,size(pin,1)
      pout(i,:)=canonicalize_partition(pin(i,:))+1
      end do
   end subroutine normalize_partition_matrix

   subroutine initial_partition(n,start,part,info)
      integer,intent(in)::n
      integer,intent(in),optional::start(:)
      integer,intent(out)::part(n),info
      integer::i
      if(present(start))then
         if(size(start)/=n)then
         info=1
         return
         end if
         part=canonicalize_partition(start)+1
      else
         do i=1,n
         part(i)=i
         end do
      end if
      info=0
   end subroutine initial_partition

   subroutine build_gibbs_candidates(part,idx,cand)
      integer,intent(in)::part(:),idx
      integer,allocatable,intent(out)::cand(:,:)
      integer::b,old,r,nc,k
      integer,allocatable::tmp(:),labels(:)
      b=maxval(part)
      old=part(idx)
      r=count(part==old)
      nc=b
      if(r>1)nc=b+1
      allocate(cand(nc,size(part)),tmp(size(part)))
      do k=1,nc
         tmp=part
         tmp(idx)=k
         tmp=canonicalize_partition(tmp)+1
         cand(k,:)=tmp
      end do
      call unique_partition_rows(cand)
   end subroutine build_gibbs_candidates

   subroutine unique_partition_rows(a)
      integer,allocatable,intent(inout)::a(:,:)
      integer,allocatable::b(:,:)
      integer::i,j,n,keep
      logical::dup
      allocate(b(size(a,1),size(a,2)))
      keep=0
      do i=1,size(a,1)
         dup=.false.
         do j=1,keep
         if(all(a(i,:)==b(j,:)))then
         dup=.true.
         exit
         end if
         end do
         if(.not.dup)then
         keep=keep+1
         b(keep,:)=a(i,:)
         end if
      end do
      a=b(1:keep,:)
   end subroutine unique_partition_rows

   subroutine modal_partition(samples,part)
      integer,intent(in)::samples(:,:)
      integer,allocatable,intent(out)::part(:)
      integer::i,j,best,best_count,c
      best=1
      best_count=0
      do i=1,size(samples,1)
         c=0
         do j=1,size(samples,1)
            if(all(samples(j,:)==samples(i,:)))c=c+1
         end do
         if(c>best_count)then
         best=i
         best_count=c
         end if
      end do
      allocate(part(size(samples,2)))
      part=samples(best,:)
   end subroutine modal_partition

   subroutine normalize_logweights(lw,w)
      real(dp),intent(in)::lw(:)
      real(dp),intent(out)::w(:)
      real(dp)::m,s
      m=maxval(lw)
      w=exp(lw-m)
      s=sum(w)
      if(s<=0.0_dp .or. .not.(s==s))then
      w=1.0_dp/real(size(w),dp)
      else
      w=w/s
      end if
   end subroutine normalize_logweights

   integer function sample_discrete(w) result(idx)
      real(dp),intent(in)::w(:)
      real(dp)::u,c
      integer::i
      u=runif1()
      c=0.0_dp
      idx=size(w)
      do i=1,size(w)
      c=c+w(i)
      if(u<c)then
      idx=i
      return
      end if
      end do
   end function sample_discrete

   pure function outer_product(a,b) result(c)
      real(dp),intent(in)::a(:),b(:)
      real(dp)::c(size(a),size(b))
      integer::i,j
      do i=1,size(a)
      do j=1,size(b)
      c(i,j)=a(i)*b(j)
      end do
      end do
   end function outer_product
end module spatialextremes_conditional
