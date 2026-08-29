module spatialextremes_design_inference
   use spatialextremes_base, only: dp, neg_huge, pair_count, distance_vectors, &
      euclidean_distances, is_finite
   use spatialextremes_univariate, only: gev_to_frechet_trend, dgev
   use spatialextremes_covariance, only: covariance_values, brown_resnick_a, geom_gauss_a, &
      mahalanobis_distances_2d, mahalanobis_distances_3d, COV_CAUGEN
   use spatialextremes_pairwise, only: lplik_smith_contributions, lplik_schlather_contributions, &
      lplik_schlather_ind_contributions, lplik_extremalt_contributions
   use spatialextremes_standard_errors, only: composite_se_t, composite_sandwich, composite_sandwich_active
   implicit none
   private

   public :: smith_design_standard_errors
   public :: schlather_design_standard_errors, schlather_ind_design_standard_errors
   public :: brownresnick_design_standard_errors, geomgauss_design_standard_errors
   public :: extremalt_design_standard_errors, spatgev_design_standard_errors
   public :: gev_design_frechet
contains
   subroutine canonical_temporal(nobs,xin,bin,x,b,ok)
      integer,intent(in)::nobs
      real(dp),intent(in),optional::xin(:,:),bin(:)
      real(dp),allocatable,intent(out)::x(:,:),b(:)
      logical,intent(out)::ok
      ok=.true.
      if(present(xin).neqv.present(bin))then
         ok=.false.
         allocate(x(nobs,0),b(0))
         return
      end if
      if(.not.present(xin))then
         allocate(x(nobs,0),b(0))
         return
      end if
      if(size(xin,1)/=nobs .or. size(xin,2)/=size(bin))then
         ok=.false.
         allocate(x(nobs,0),b(0))
         return
      end if
      allocate(x(nobs,size(bin)),b(size(bin)))
      x=xin
      b=bin
   end subroutine canonical_temporal

   subroutine make_theta(dep,bloc,bscale,bshape,btl,bts,bth,theta)
      real(dp),intent(in)::dep(:),bloc(:),bscale(:),bshape(:),btl(:),bts(:),bth(:)
      real(dp),allocatable,intent(out)::theta(:)
      integer::k
      allocate(theta(size(dep)+size(bloc)+size(bscale)+size(bshape)+size(btl)+size(bts)+size(bth)))
      k=0
      theta(k+1:k+size(dep))=dep
      k=k+size(dep)
      theta(k+1:k+size(bloc))=bloc
      k=k+size(bloc)
      theta(k+1:k+size(bscale))=bscale
      k=k+size(bscale)
      theta(k+1:k+size(bshape))=bshape
      k=k+size(bshape)
      if(size(btl)>0)then
      theta(k+1:k+size(btl))=btl
      k=k+size(btl)
      end if
      if(size(bts)>0)then
      theta(k+1:k+size(bts))=bts
      k=k+size(bts)
      end if
      if(size(bth)>0)theta(k+1:k+size(bth))=bth
   end subroutine make_theta

   subroutine bad_result(res,code)
      type(composite_se_t),intent(out)::res
      integer,intent(in)::code
      real(dp)::dummy(2,1),theta(1)
      theta=0.0_dp
      call composite_sandwich(dummy_contrib,theta,res)
      res%info=code
   contains
      subroutine dummy_contrib(x,c)
         real(dp),intent(in)::x(:)
         real(dp),allocatable,intent(out)::c(:,:)
         allocate(c(2,1))
         c=neg_huge
      end subroutine dummy_contrib
   end subroutine bad_result

   subroutine run_se(fn,theta,res,rel_step,active)
      interface
         subroutine fn(theta,contrib)
            import dp
            real(dp),intent(in)::theta(:)
            real(dp),allocatable,intent(out)::contrib(:,:)
         end subroutine fn
      end interface
      real(dp),intent(in)::theta(:)
      type(composite_se_t),intent(out)::res
      real(dp),intent(in),optional::rel_step
      logical,intent(in),optional::active(:)
      if(present(active))then
         if(present(rel_step))then
            call composite_sandwich_active(fn,theta,active,res,rel_step)
         else
            call composite_sandwich_active(fn,theta,active,res)
         end if
      else
         if(present(rel_step))then
            call composite_sandwich(fn,theta,res,rel_step)
         else
            call composite_sandwich(fn,theta,res)
         end if
      end if
   end subroutine run_se

   subroutine gev_design_frechet(theta,ndep,data,xloc,xscale,xshape,txloc,txscale,txshape,frech,jac,info)
      integer,intent(in)::ndep
      real(dp),intent(in)::theta(:),data(:,:),xloc(:,:),xscale(:,:),xshape(:,:)
      real(dp),intent(in)::txloc(:,:),txscale(:,:),txshape(:,:)
      real(dp),intent(out)::frech(size(data,1),size(data,2)),jac(size(data,1),size(data,2))
      integer,intent(out)::info
      real(dp)::loc(size(data,2)),scale(size(data,2)),shape(size(data,2))
      real(dp)::tl(size(data,1)),ts(size(data,1)),th(size(data,1))
      integer::k,p1,p2,p3,q1,q2,q3,nneed
      info=0
      p1=size(xloc,2)
      p2=size(xscale,2)
      p3=size(xshape,2)
      q1=size(txloc,2)
      q2=size(txscale,2)
      q3=size(txshape,2)
      if(size(xloc,1)/=size(data,2).or.size(xscale,1)/=size(data,2).or.size(xshape,1)/=size(data,2))then
         info=10
         return
      end if
      if(size(txloc,1)/=size(data,1).or.size(txscale,1)/=size(data,1).or.size(txshape,1)/=size(data,1))then
         info=11
         return
      end if
      nneed=ndep+p1+p2+p3+q1+q2+q3
      if(size(theta)/=nneed)then
      info=12
      return
      end if
      k=ndep
      loc=matmul(xloc,theta(k+1:k+p1))
      k=k+p1
      scale=matmul(xscale,theta(k+1:k+p2))
      k=k+p2
      shape=matmul(xshape,theta(k+1:k+p3))
      k=k+p3
      if(any(scale<=0.0_dp).or.any(shape<=-1.0_dp))then
      info=1
      return
      end if
      tl=0.0_dp
      ts=0.0_dp
      th=0.0_dp
      if(q1>0)then
      tl=matmul(txloc,theta(k+1:k+q1))
      k=k+q1
      end if
      if(q2>0)then
      ts=matmul(txscale,theta(k+1:k+q2))
      k=k+q2
      end if
      if(q3>0)th=matmul(txshape,theta(k+1:k+q3))
      call gev_to_frechet_trend(data,loc,scale,shape,tl,ts,th,frech,jac,info)
   end subroutine gev_design_frechet

   function smith_design_standard_errors(data,coord,cov,xloc,xscale,xshape,bloc,bscale,bshape, &
      temp_xloc,temp_xscale,temp_xshape,btemp_loc,btemp_scale,btemp_shape,weights,rel_step,active,isotropic) result(res)
      real(dp),intent(in)::data(:,:),coord(:,:),cov(:,:),xloc(:,:),xscale(:,:),xshape(:,:)
      real(dp),intent(in)::bloc(:),bscale(:),bshape(:)
      real(dp),intent(in),optional::temp_xloc(:,:),temp_xscale(:,:),temp_xshape(:,:)
      real(dp),intent(in),optional::btemp_loc(:),btemp_scale(:),btemp_shape(:),weights(:),rel_step
      logical,intent(in),optional::active(:),isotropic
      type(composite_se_t)::res
      real(dp),allocatable::tlx(:,:),tsx(:,:),thx(:,:),btl(:),bts(:),bth(:),theta(:),dv(:,:)
      real(dp),allocatable::dep(:)
      logical::ok1,ok2,ok3,iso
      integer::d,ndep
      call canonical_temporal(size(data,1),temp_xloc,btemp_loc,tlx,btl,ok1)
      call canonical_temporal(size(data,1),temp_xscale,btemp_scale,tsx,bts,ok2)
      call canonical_temporal(size(data,1),temp_xshape,btemp_shape,thx,bth,ok3)
      d=size(coord,2)
      iso=.false.
      if(present(isotropic))iso=isotropic
      if(.not.(ok1.and.ok2.and.ok3).or.(d/=2.and.d/=3))then
      call bad_result(res,30)
      return
      end if
      if(size(cov,1)<d.or.size(cov,2)<d)then
      call bad_result(res,31)
      return
      end if
      if(iso)then
         ndep=1
         allocate(dep(1))
         dep(1)=cov(1,1)
      else if(d==2)then
         ndep=3
         allocate(dep(3))
         dep=[cov(1,1),cov(1,2),cov(2,2)]
      else
         ndep=6
         allocate(dep(6))
         dep=[cov(1,1),cov(1,2),cov(1,3),cov(2,2),cov(2,3),cov(3,3)]
      end if
      dv=distance_vectors(coord)
      call make_theta(dep,bloc,bscale,bshape,btl,bts,bth,theta)
      call run_se(contrib,theta,res,rel_step,active)
   contains
      subroutine contrib(x,c)
         real(dp),intent(in)::x(:)
         real(dp),allocatable,intent(out)::c(:,:)
         real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2)),cmat(d,d)
         real(dp),allocatable::a(:)
         integer::info
         cmat=0.0_dp
         if(iso)then
            if(x(1)<=0.0_dp)then
            allocate(c(size(data,1),pair_count(size(data,2))))
            c=neg_huge
            return
            end if
            cmat=0.0_dp
            cmat(1,1)=x(1)
            cmat(2,2)=x(1)
            if(d==3)cmat(3,3)=x(1)
         else if(d==2)then
            cmat(1,1)=x(1)
            cmat(1,2)=x(2)
            cmat(2,1)=x(2)
            cmat(2,2)=x(3)
         else
            cmat(1,1)=x(1)
            cmat(1,2)=x(2)
            cmat(2,1)=x(2)
            cmat(1,3)=x(3)
            cmat(3,1)=x(3)
            cmat(2,2)=x(4)
            cmat(2,3)=x(5)
            cmat(3,2)=x(5)
            cmat(3,3)=x(6)
         end if
         if(d==2)then
         a=mahalanobis_distances_2d(dv,cmat(1:2,1:2))
         else
         a=mahalanobis_distances_3d(dv,cmat)
         end if
         if(any(a/=a).or.any(a<=0.0_dp))then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         call gev_design_frechet(x,ndep,data,xloc,xscale,xshape,tlx,tsx,thx,f,j,info)
         if(info/=0)then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         if(present(weights))then
         c=lplik_smith_contributions(f,a,j,weights)
         else
         c=lplik_smith_contributions(f,a,j)
         end if
      end subroutine contrib
   end function smith_design_standard_errors

   function schlather_design_standard_errors(data,coord,covmod,nugget,range,smooth,xloc,xscale,xshape, &
      bloc,bscale,bshape,smooth2,temp_xloc,temp_xscale,temp_xshape,btemp_loc,btemp_scale,btemp_shape, &
      weights,rel_step,active) result(res)
      real(dp),intent(in)::data(:,:),coord(:,:),nugget,range,smooth,xloc(:,:),xscale(:,:),xshape(:,:)
      integer,intent(in)::covmod
      real(dp),intent(in)::bloc(:),bscale(:),bshape(:)
      real(dp),intent(in),optional::smooth2,temp_xloc(:,:),temp_xscale(:,:),temp_xshape(:,:)
      real(dp),intent(in),optional::btemp_loc(:),btemp_scale(:),btemp_shape(:),weights(:),rel_step
      logical,intent(in),optional::active(:)
      type(composite_se_t)::res
      real(dp),allocatable::tlx(:,:),tsx(:,:),thx(:,:),btl(:),bts(:),bth(:),theta(:),dists(:),dep(:)
      logical::ok1,ok2,ok3
      integer::ndep
      real(dp)::s20
      call canonical_temporal(size(data,1),temp_xloc,btemp_loc,tlx,btl,ok1)
      call canonical_temporal(size(data,1),temp_xscale,btemp_scale,tsx,bts,ok2)
      call canonical_temporal(size(data,1),temp_xshape,btemp_shape,thx,bth,ok3)
      if(.not.(ok1.and.ok2.and.ok3))then
      call bad_result(res,30)
      return
      end if
      s20=1.0_dp
      if(present(smooth2))s20=smooth2
      if(covmod==COV_CAUGEN)then
      ndep=4
      dep=[nugget,range,smooth,s20]
      else
      ndep=3
      dep=[nugget,range,smooth]
      end if
      dists=euclidean_distances(coord)
      call make_theta(dep,bloc,bscale,bshape,btl,bts,bth,theta)
      call run_se(contrib,theta,res,rel_step,active)
   contains
      subroutine contrib(x,c)
         real(dp),intent(in)::x(:)
         real(dp),allocatable,intent(out)::c(:,:)
         real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2)),s2
         real(dp),allocatable::rho(:)
         integer::i,info
         if(x(1)<0.0_dp.or.x(1)>=1.0_dp.or.x(2)<=0.0_dp.or.x(3)<=0.0_dp)then
            allocate(c(size(data,1),pair_count(size(data,2))))
            c=neg_huge
            return
         end if
         s2=1.0_dp
         if(covmod==COV_CAUGEN)then
         s2=x(4)
         if(s2<=0.0_dp)then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         end if
         allocate(rho(size(dists)))
         do i=1,size(dists)
         rho(i)=covariance_values(dists(i),covmod,x(1),1-x(1),x(2),x(3),s2,size(coord,2))
         end do
         if(any(rho/=rho))then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         call gev_design_frechet(x,ndep,data,xloc,xscale,xshape,tlx,tsx,thx,f,j,info)
         if(info/=0)then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         if(present(weights))then
         c=lplik_schlather_contributions(f,rho,j,weights)
         else
         c=lplik_schlather_contributions(f,rho,j)
         end if
      end subroutine contrib
   end function schlather_design_standard_errors

   function schlather_ind_design_standard_errors(data,coord,covmod,alpha,nugget,range,smooth,xloc,xscale,xshape, &
      bloc,bscale,bshape,smooth2,temp_xloc,temp_xscale,temp_xshape,btemp_loc,btemp_scale,btemp_shape, &
      weights,rel_step,active) result(res)
      real(dp),intent(in)::data(:,:),coord(:,:),alpha,nugget,range,smooth,xloc(:,:),xscale(:,:),xshape(:,:)
      integer,intent(in)::covmod
      real(dp),intent(in)::bloc(:),bscale(:),bshape(:)
      real(dp),intent(in),optional::smooth2,temp_xloc(:,:),temp_xscale(:,:),temp_xshape(:,:)
      real(dp),intent(in),optional::btemp_loc(:),btemp_scale(:),btemp_shape(:),weights(:),rel_step
      logical,intent(in),optional::active(:)
      type(composite_se_t)::res
      real(dp),allocatable::tlx(:,:),tsx(:,:),thx(:,:),btl(:),bts(:),bth(:),theta(:),dists(:),dep(:)
      logical::ok1,ok2,ok3
      integer::ndep
      real(dp)::s20
      call canonical_temporal(size(data,1),temp_xloc,btemp_loc,tlx,btl,ok1)
      call canonical_temporal(size(data,1),temp_xscale,btemp_scale,tsx,bts,ok2)
      call canonical_temporal(size(data,1),temp_xshape,btemp_shape,thx,bth,ok3)
      if(.not.(ok1.and.ok2.and.ok3))then
      call bad_result(res,30)
      return
      end if
      s20=1.0_dp
      if(present(smooth2))s20=smooth2
      if(covmod==COV_CAUGEN)then
      ndep=5
      dep=[alpha,nugget,range,smooth,s20]
      else
      ndep=4
      dep=[alpha,nugget,range,smooth]
      end if
      dists=euclidean_distances(coord)
      call make_theta(dep,bloc,bscale,bshape,btl,bts,bth,theta)
      call run_se(contrib,theta,res,rel_step,active)
   contains
      subroutine contrib(x,c)
         real(dp),intent(in)::x(:)
         real(dp),allocatable,intent(out)::c(:,:)
         real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2)),s2
         real(dp),allocatable::rho(:)
         integer::i,info,ir,inug
         ir=3
         inug=2
         if(x(1)<0.0_dp.or.x(1)>1.0_dp.or.x(inug)<0.0_dp.or.x(inug)>=1.0_dp.or.x(ir)<=0.0_dp.or.x(4)<=0.0_dp)then
            allocate(c(size(data,1),pair_count(size(data,2))))
            c=neg_huge
            return
         end if
         s2=1.0_dp
         if(covmod==COV_CAUGEN)then
         s2=x(5)
         if(s2<=0.0_dp)then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         end if
         allocate(rho(size(dists)))
         do i=1,size(dists)
         rho(i)=covariance_values(dists(i),covmod,x(2),1-x(2),x(3),x(4),s2,size(coord,2))
         end do
         if(any(rho/=rho))then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         call gev_design_frechet(x,ndep,data,xloc,xscale,xshape,tlx,tsx,thx,f,j,info)
         if(info/=0)then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         if(present(weights))then
         c=lplik_schlather_ind_contributions(f,x(1),rho,j,weights)
         else
         c=lplik_schlather_ind_contributions(f,x(1),rho,j)
         end if
      end subroutine contrib
   end function schlather_ind_design_standard_errors

   function brownresnick_design_standard_errors(data,coord,range,smooth,xloc,xscale,xshape,bloc,bscale,bshape, &
      temp_xloc,temp_xscale,temp_xshape,btemp_loc,btemp_scale,btemp_shape,weights,rel_step,active) result(res)
      real(dp),intent(in)::data(:,:),coord(:,:),range,smooth,xloc(:,:),xscale(:,:),xshape(:,:)
      real(dp),intent(in)::bloc(:),bscale(:),bshape(:)
      real(dp),intent(in),optional::temp_xloc(:,:),temp_xscale(:,:),temp_xshape(:,:), &
         btemp_loc(:),btemp_scale(:),btemp_shape(:),weights(:),rel_step
      logical,intent(in),optional::active(:)
      type(composite_se_t)::res
      real(dp),allocatable::tlx(:,:),tsx(:,:),thx(:,:),btl(:),bts(:),bth(:),theta(:),dists(:),dep(:)
      logical::ok1,ok2,ok3
      integer,parameter::ndep=2
      call canonical_temporal(size(data,1),temp_xloc,btemp_loc,tlx,btl,ok1)
      call canonical_temporal(size(data,1),temp_xscale,btemp_scale,tsx,bts,ok2)
      call canonical_temporal(size(data,1),temp_xshape,btemp_shape,thx,bth,ok3)
      if(.not.(ok1.and.ok2.and.ok3))then
      call bad_result(res,30)
      return
      end if
      dep=[range,smooth]
      dists=euclidean_distances(coord)
      call make_theta(dep,bloc,bscale,bshape,btl,bts,bth,theta)
      call run_se(contrib,theta,res,rel_step,active)
   contains
      subroutine contrib(x,c)
         real(dp),intent(in)::x(:)
         real(dp),allocatable,intent(out)::c(:,:)
         real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2))
         real(dp),allocatable::a(:)
         integer::i,info
         if(x(1)<=0.0_dp.or.x(2)<=0.0_dp.or.x(2)>2.0_dp)then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         allocate(a(size(dists)))
         do i=1,size(dists)
         a(i)=brown_resnick_a(dists(i),x(1),x(2))
         end do
         call gev_design_frechet(x,ndep,data,xloc,xscale,xshape,tlx,tsx,thx,f,j,info)
         if(info/=0)then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         if(present(weights))then
         c=lplik_smith_contributions(f,a,j,weights)
         else
         c=lplik_smith_contributions(f,a,j)
         end if
      end subroutine contrib
   end function brownresnick_design_standard_errors

   function geomgauss_design_standard_errors(data,coord,covmod,sigma2,nugget,range,smooth,xloc,xscale,xshape, &
      bloc,bscale,bshape,smooth2,temp_xloc,temp_xscale,temp_xshape,btemp_loc,btemp_scale,btemp_shape, &
      weights,rel_step,active) result(res)
      real(dp),intent(in)::data(:,:),coord(:,:),sigma2,nugget,range,smooth,xloc(:,:),xscale(:,:),xshape(:,:)
      integer,intent(in)::covmod
      real(dp),intent(in)::bloc(:),bscale(:),bshape(:)
      real(dp),intent(in),optional::smooth2,temp_xloc(:,:),temp_xscale(:,:),temp_xshape(:,:), &
         btemp_loc(:),btemp_scale(:),btemp_shape(:),weights(:),rel_step
      logical,intent(in),optional::active(:)
      type(composite_se_t)::res
      real(dp),allocatable::tlx(:,:),tsx(:,:),thx(:,:),btl(:),bts(:),bth(:),theta(:),dists(:),dep(:)
      logical::ok1,ok2,ok3
      integer::ndep
      real(dp)::s20
      call canonical_temporal(size(data,1),temp_xloc,btemp_loc,tlx,btl,ok1)
      call canonical_temporal(size(data,1),temp_xscale,btemp_scale,tsx,bts,ok2)
      call canonical_temporal(size(data,1),temp_xshape,btemp_shape,thx,bth,ok3)
      if(.not.(ok1.and.ok2.and.ok3))then
      call bad_result(res,30)
      return
      end if
      s20=1.0_dp
      if(present(smooth2))s20=smooth2
      if(covmod==COV_CAUGEN)then
      ndep=5
      dep=[sigma2,nugget,range,smooth,s20]
      else
      ndep=4
      dep=[sigma2,nugget,range,smooth]
      end if
      dists=euclidean_distances(coord)
      call make_theta(dep,bloc,bscale,bshape,btl,bts,bth,theta)
      call run_se(contrib,theta,res,rel_step,active)
   contains
      subroutine contrib(x,c)
         real(dp),intent(in)::x(:)
         real(dp),allocatable,intent(out)::c(:,:)
         real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2)),s2
         real(dp),allocatable::a(:)
         integer::i,info
         if(x(1)<=0.0_dp.or.x(2)<0.0_dp.or.x(2)>=1.0_dp.or.x(3)<=0.0_dp.or.x(4)<=0.0_dp)then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         s2=1.0_dp
         if(covmod==COV_CAUGEN)then
         s2=x(5)
         if(s2<=0.0_dp)then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         end if
         allocate(a(size(dists)))
         do i=1,size(dists)
         a(i)=geom_gauss_a(dists(i),covmod,x(1),x(2),x(3),x(4),s2,size(coord,2))
         end do
         if(any(a/=a).or.any(a<=0.0_dp))then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         call gev_design_frechet(x,ndep,data,xloc,xscale,xshape,tlx,tsx,thx,f,j,info)
         if(info/=0)then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         if(present(weights))then
         c=lplik_smith_contributions(f,a,j,weights)
         else
         c=lplik_smith_contributions(f,a,j)
         end if
      end subroutine contrib
   end function geomgauss_design_standard_errors

   function extremalt_design_standard_errors(data,coord,covmod,nugget,range,smooth,nu,xloc,xscale,xshape, &
      bloc,bscale,bshape,smooth2,temp_xloc,temp_xscale,temp_xshape,btemp_loc,btemp_scale,btemp_shape, &
      weights,rel_step,active) result(res)
      real(dp),intent(in)::data(:,:),coord(:,:),nugget,range,smooth,nu,xloc(:,:),xscale(:,:),xshape(:,:)
      integer,intent(in)::covmod
      real(dp),intent(in)::bloc(:),bscale(:),bshape(:)
      real(dp),intent(in),optional::smooth2,temp_xloc(:,:),temp_xscale(:,:),temp_xshape(:,:), &
         btemp_loc(:),btemp_scale(:),btemp_shape(:),weights(:),rel_step
      logical,intent(in),optional::active(:)
      type(composite_se_t)::res
      real(dp),allocatable::tlx(:,:),tsx(:,:),thx(:,:),btl(:),bts(:),bth(:),theta(:),dists(:),dep(:)
      logical::ok1,ok2,ok3
      integer::ndep
      real(dp)::s20
      call canonical_temporal(size(data,1),temp_xloc,btemp_loc,tlx,btl,ok1)
      call canonical_temporal(size(data,1),temp_xscale,btemp_scale,tsx,bts,ok2)
      call canonical_temporal(size(data,1),temp_xshape,btemp_shape,thx,bth,ok3)
      if(.not.(ok1.and.ok2.and.ok3))then
      call bad_result(res,30)
      return
      end if
      s20=1.0_dp
      if(present(smooth2))s20=smooth2
      if(covmod==COV_CAUGEN)then
      ndep=5
      dep=[nugget,range,smooth,s20,nu]
      else
      ndep=4
      dep=[nugget,range,smooth,nu]
      end if
      dists=euclidean_distances(coord)
      call make_theta(dep,bloc,bscale,bshape,btl,bts,bth,theta)
      call run_se(contrib,theta,res,rel_step,active)
   contains
      subroutine contrib(x,c)
         real(dp),intent(in)::x(:)
         real(dp),allocatable,intent(out)::c(:,:)
         real(dp)::f(size(data,1),size(data,2)),j(size(data,1),size(data,2)),s2,dof
         real(dp),allocatable::rho(:)
         integer::i,info,idof
         idof=ndep
         dof=x(idof)
         if(x(1)<0.0_dp.or.x(1)>=1.0_dp.or.x(2)<=0.0_dp.or.x(3)<=0.0_dp.or.dof<=0.0_dp)then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         s2=1.0_dp
         if(covmod==COV_CAUGEN)then
         s2=x(4)
         if(s2<=0.0_dp)then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         end if
         allocate(rho(size(dists)))
         do i=1,size(dists)
         rho(i)=covariance_values(dists(i),covmod,x(1),1-x(1),x(2),x(3),s2,size(coord,2))
         end do
         if(any(rho/=rho))then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         call gev_design_frechet(x,ndep,data,xloc,xscale,xshape,tlx,tsx,thx,f,j,info)
         if(info/=0)then
         allocate(c(size(data,1),pair_count(size(data,2))))
         c=neg_huge
         return
         end if
         if(present(weights))then
         c=lplik_extremalt_contributions(f,rho,dof,j,weights)
         else
         c=lplik_extremalt_contributions(f,rho,dof,j)
         end if
      end subroutine contrib
   end function extremalt_design_standard_errors

   function spatgev_design_standard_errors(data,xloc,xscale,xshape,bloc,bscale,bshape,temp_xloc,temp_xscale,temp_xshape, &
      btemp_loc,btemp_scale,btemp_shape,rel_step,active) result(res)
      real(dp),intent(in)::data(:,:),xloc(:,:),xscale(:,:),xshape(:,:),bloc(:),bscale(:),bshape(:)
      real(dp),intent(in),optional::temp_xloc(:,:),temp_xscale(:,:),temp_xshape(:,:), &
         btemp_loc(:),btemp_scale(:),btemp_shape(:),rel_step
      logical,intent(in),optional::active(:)
      type(composite_se_t)::res
      real(dp),allocatable::tlx(:,:),tsx(:,:),thx(:,:),btl(:),bts(:),bth(:),theta(:),dep(:)
      logical::ok1,ok2,ok3
      allocate(dep(0))
      call canonical_temporal(size(data,1),temp_xloc,btemp_loc,tlx,btl,ok1)
      call canonical_temporal(size(data,1),temp_xscale,btemp_scale,tsx,bts,ok2)
      call canonical_temporal(size(data,1),temp_xshape,btemp_shape,thx,bth,ok3)
      if(.not.(ok1.and.ok2.and.ok3))then
      call bad_result(res,30)
      return
      end if
      call make_theta(dep,bloc,bscale,bshape,btl,bts,bth,theta)
      call run_se(contrib,theta,res,rel_step,active)
   contains
      subroutine contrib(x,c)
         real(dp),intent(in)::x(:)
         real(dp),allocatable,intent(out)::c(:,:)
         real(dp)::loc(size(data,2)),scale(size(data,2)),shape(size(data,2)),tl(size(data,1)),ts(size(data,1)),th(size(data,1))
         integer::k,p1,p2,p3,q1,q2,q3,i,jj
         p1=size(xloc,2)
         p2=size(xscale,2)
         p3=size(xshape,2)
         q1=size(tlx,2)
         q2=size(tsx,2)
         q3=size(thx,2)
         k=0
         loc=matmul(xloc,x(k+1:k+p1))
         k=k+p1
         scale=matmul(xscale,x(k+1:k+p2))
         k=k+p2
         shape=matmul(xshape,x(k+1:k+p3))
         k=k+p3
         if(any(scale<=0.0_dp).or.any(shape<=-1.0_dp))then
         allocate(c(size(data,1),size(data,2)))
         c=neg_huge
         return
         end if
         tl=0.0_dp
         ts=0.0_dp
         th=0.0_dp
         if(q1>0)then
         tl=matmul(tlx,x(k+1:k+q1))
         k=k+q1
         end if
         if(q2>0)then
         ts=matmul(tsx,x(k+1:k+q2))
         k=k+q2
         end if
         if(q3>0)th=matmul(thx,x(k+1:k+q3))
         allocate(c(size(data,1),size(data,2)))
         c=0.0_dp
         do jj=1,size(data,2)
         do i=1,size(data,1)
            if(is_finite(data(i,jj)))then
               if(scale(jj)+ts(i)<=0.0_dp)then
               c=neg_huge
               return
               end if
               c(i,jj)=dgev(data(i,jj),loc(jj)+tl(i),scale(jj)+ts(i),shape(jj)+th(i),.true.)
               if(abs(c(i,jj))>=huge(1.0_dp))then
               c=neg_huge
               return
               end if
            end if
         end do
         end do
      end subroutine contrib
   end function spatgev_design_standard_errors
end module spatialextremes_design_inference
