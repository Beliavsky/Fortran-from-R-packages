module spatialextremes_standard_errors
   use spatialextremes_base, only: dp, neg_huge, nan_dp, inverse_spd, pair_count, &
      euclidean_distances, distance_vectors
   use spatialextremes_covariance, only: covariance_values, brown_resnick_a, geom_gauss_a, &
      mahalanobis_distances_2d, mahalanobis_distances_3d, COV_CAUGEN
   use spatialextremes_pairwise, only: lplik_smith_contributions, lplik_schlather_contributions, &
      lplik_schlather_ind_contributions, lplik_extremalt_contributions
   implicit none
   private

   type, public :: composite_se_t
      real(dp), allocatable :: gradient(:)
      real(dp), allocatable :: var_score(:,:)
      real(dp), allocatable :: hessian(:,:)
      real(dp), allocatable :: ihessian(:,:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: stderr(:)
      real(dp), allocatable :: correlation(:,:)
      integer :: info = 0
   end type composite_se_t

   public :: composite_sandwich, composite_sandwich_active
   public :: smith_frechet_standard_errors, schlather_frechet_standard_errors
   public :: schlather_ind_frechet_standard_errors, brownresnick_frechet_standard_errors
   public :: geomgauss_frechet_standard_errors, extremalt_frechet_standard_errors

   abstract interface
      subroutine contribution_function(theta, contrib)
         import dp
         real(dp), intent(in) :: theta(:)
         real(dp), allocatable, intent(out) :: contrib(:,:)
      end subroutine contribution_function
   end interface
contains
   subroutine composite_sandwich(fn,theta,res,rel_step)
      ! Numerical counterpart of SpatialExtremes' standard-error backend.
      ! The callback returns one composite log-likelihood contribution for
      ! every independent replicate (rows) and spatial pair (columns).
      ! As upstream, J is estimated from the covariance of replicate scores
      ! and H from the covariance of pairwise contribution scores.
      procedure(contribution_function) :: fn
      real(dp), intent(in) :: theta(:)
      type(composite_se_t), intent(out) :: res
      real(dp), intent(in), optional :: rel_step
      real(dp), allocatable :: c0(:,:),cp(:,:),cm(:,:),pair_score(:,:),obs_score(:,:)
      real(dp) :: h,step0
      integer :: nobs,npairs,npar,j,tries,info,i,k,row
      logical :: okp,okm

      call fn(theta,c0)
      nobs=size(c0,1)
      npairs=size(c0,2)
      npar=size(theta)
      call allocate_result(res,npar)
      if(nobs<2 .or. npairs<1 .or. npar<1 .or. .not.valid_contrib(c0))then
         res%info=1
         return
      end if
      allocate(pair_score(nobs*npairs,npar),obs_score(nobs,npar))
      pair_score=0.0_dp
      obs_score=0.0_dp
      step0=1.0e-5_dp
      if(present(rel_step))step0=rel_step

      do j=1,npar
         h=step0*max(1.0_dp,abs(theta(j)))
         okp=.false.
         okm=.false.
         do tries=1,12
            call perturbed_contrib(fn,theta,j,h,cp)
            call perturbed_contrib(fn,theta,j,-h,cm)
            okp=valid_contrib(cp)
            okm=valid_contrib(cm)
            if(okp.and.okm)exit
            h=0.25_dp*h
         end do
         if(.not.(okp.and.okm))then
            res%info=2
            return
         end if
         row=0
         do k=1,npairs
            do i=1,nobs
               row=row+1
               pair_score(row,j)=(cp(i,k)-cm(i,k))/(2.0_dp*h)
            end do
         end do
         obs_score(:,j)=sum((cp-cm)/(2.0_dp*h),dim=2)
      end do

      res%gradient=sum(obs_score,dim=1)
      res%var_score=sample_cov(obs_score)*real(nobs,dp)
      res%hessian=sample_cov(pair_score)*real(nobs*npairs,dp)
      call inverse_spd(res%hessian,res%ihessian,info)
      if(info/=0)then
         res%info=10+info
         return
      end if
      res%covariance=matmul(matmul(res%ihessian,res%var_score),res%ihessian)
      do i=1,npar
         if(res%covariance(i,i)>0.0_dp)then
            res%stderr(i)=sqrt(res%covariance(i,i))
         else
            res%stderr(i)=nan_dp()
         end if
      end do
      res%correlation=0.0_dp
      do i=1,npar
         do j=1,npar
            if(res%stderr(i)>0.0_dp .and. res%stderr(j)>0.0_dp) &
               res%correlation(i,j)=res%covariance(i,j)/(res%stderr(i)*res%stderr(j))
         end do
      end do
      res%info=0
   end subroutine composite_sandwich

   subroutine composite_sandwich_active(fn,theta,active,res,rel_step)
      ! Sandwich calculation for a subset of a full parameter vector.
      ! Inactive entries remain fixed at the supplied theta values.
      procedure(contribution_function) :: fn
      real(dp), intent(in) :: theta(:)
      logical, intent(in) :: active(:)
      type(composite_se_t), intent(out) :: res
      real(dp), intent(in), optional :: rel_step
      real(dp), allocatable :: reduced(:)
      integer, allocatable :: idx(:)
      integer :: i,k,nactive

      if(size(active)/=size(theta))then
         call allocate_result(res,1)
         res%info=20
         return
      end if
      nactive=count(active)
      if(nactive<1)then
         call allocate_result(res,1)
         res%info=21
         return
      end if
      allocate(idx(nactive),reduced(nactive))
      k=0
      do i=1,size(theta)
         if(active(i))then
            k=k+1
            idx(k)=i
            reduced(k)=theta(i)
         end if
      end do
      if(present(rel_step))then
         call composite_sandwich(reduced_contrib,reduced,res,rel_step)
      else
         call composite_sandwich(reduced_contrib,reduced,res)
      end if
   contains
      subroutine reduced_contrib(x,c)
         real(dp),intent(in)::x(:)
         real(dp),allocatable,intent(out)::c(:,:)
         real(dp)::full(size(theta))
         integer::j
         full=theta
         do j=1,size(idx)
         full(idx(j))=x(j)
         end do
         call fn(full,c)
      end subroutine reduced_contrib
   end subroutine composite_sandwich_active

   subroutine sandwich_dispatch(fn,theta,res,rel_step,active)
      procedure(contribution_function)::fn
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
   end subroutine sandwich_dispatch

   subroutine allocate_result(res,n)
      type(composite_se_t),intent(inout)::res
      integer,intent(in)::n
      allocate(res%gradient(n),res%var_score(n,n),res%hessian(n,n),res%ihessian(n,n), &
         res%covariance(n,n),res%stderr(n),res%correlation(n,n))
      res%gradient=nan_dp()
      res%var_score=nan_dp()
      res%hessian=nan_dp()
      res%ihessian=nan_dp()
      res%covariance=nan_dp()
      res%stderr=nan_dp()
      res%correlation=nan_dp()
   end subroutine allocate_result

   subroutine perturbed_contrib(fn,theta,j,delta,c)
      procedure(contribution_function)::fn
      real(dp),intent(in)::theta(:),delta
      integer,intent(in)::j
      real(dp),allocatable,intent(out)::c(:,:)
      real(dp)::x(size(theta))
      x=theta
      x(j)=x(j)+delta
      call fn(x,c)
   end subroutine perturbed_contrib

   logical function valid_contrib(c) result(ok)
      real(dp),intent(in)::c(:,:)
      ok=all(c==c).and.all(abs(c)<huge(1.0_dp)).and.all(c>0.5_dp*neg_huge)
   end function valid_contrib

   function sample_cov(x) result(v)
      real(dp),intent(in)::x(:,:)
      real(dp)::v(size(x,2),size(x,2)),z(size(x,1),size(x,2)),mu(size(x,2))
      integer::n
      n=size(x,1)
      if(n<=1)then
      v=0.0_dp
      return
      end if
      mu=sum(x,dim=1)/real(n,dp)
      z=x-spread(mu,1,n)
      v=matmul(transpose(z),z)/real(n-1,dp)
   end function sample_cov

   function smith_frechet_standard_errors(data,coord,cov,weights,rel_step,active) result(res)
      real(dp),intent(in)::data(:,:),coord(:,:),cov(:,:)
      real(dp),intent(in),optional::weights(:),rel_step
      logical,intent(in),optional::active(:)
      type(composite_se_t)::res
      real(dp),allocatable::dv(:,:),theta(:)
      integer::d
      d=size(coord,2)
      dv=distance_vectors(coord)
      if(d==2)then
         theta=[cov(1,1),cov(1,2),cov(2,2)]
      else if(d==3)then
         theta=[cov(1,1),cov(1,2),cov(1,3),cov(2,2),cov(2,3),cov(3,3)]
      else
         call allocate_result(res,1)
         res%info=3
         return
      end if
      call sandwich_dispatch(contrib,theta,res,rel_step,active)
   contains
      subroutine contrib(x,c)
         real(dp),intent(in)::x(:)
         real(dp),allocatable,intent(out)::c(:,:)
         real(dp)::cmat(d,d),jac(size(data,1),size(data,2))
         real(dp),allocatable::a(:)
         jac=0.0_dp
         cmat=0.0_dp
         if(d==2)then
            cmat(1,1)=x(1)
            cmat(1,2)=x(2)
            cmat(2,1)=x(2)
            cmat(2,2)=x(3)
            a=mahalanobis_distances_2d(dv,cmat)
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
            a=mahalanobis_distances_3d(dv,cmat)
         end if
         if(any(a/=a) .or. any(a<=0.0_dp))then
            allocate(c(size(data,1),pair_count(size(data,2))))
            c=neg_huge
            return
         end if
         if(present(weights))then
         c=lplik_smith_contributions(data,a,jac,weights)
         else
         c=lplik_smith_contributions(data,a,jac)
         end if
      end subroutine contrib
   end function smith_frechet_standard_errors

   function brownresnick_frechet_standard_errors(data,coord,range,smooth,weights,rel_step,active) result(res)
      real(dp),intent(in)::data(:,:),coord(:,:),range,smooth
      real(dp),intent(in),optional::weights(:),rel_step
      logical,intent(in),optional::active(:)
      type(composite_se_t)::res
      real(dp),allocatable::d(:)
      real(dp)::theta(2)
      d=euclidean_distances(coord)
      theta=[range,smooth]
      call sandwich_dispatch(contrib,theta,res,rel_step,active)
   contains
      subroutine contrib(x,c)
         real(dp),intent(in)::x(:)
         real(dp),allocatable,intent(out)::c(:,:)
         real(dp)::jac(size(data,1),size(data,2))
         real(dp),allocatable::a(:)
         integer::i
         if(x(1)<=0.0_dp.or.x(2)<=0.0_dp.or.x(2)>2.0_dp)then
            allocate(c(size(data,1),pair_count(size(data,2))))
            c=neg_huge
            return
         end if
         allocate(a(size(d)))
         do i=1,size(d)
         a(i)=brown_resnick_a(d(i),x(1),x(2))
         end do
         jac=0.0_dp
         if(present(weights))then
         c=lplik_smith_contributions(data,a,jac,weights)
         else
         c=lplik_smith_contributions(data,a,jac)
         end if
      end subroutine contrib
   end function brownresnick_frechet_standard_errors

   function schlather_frechet_standard_errors(data,coord,covmod,nugget,range,smooth,smooth2,weights,rel_step,active) result(res)
      real(dp),intent(in)::data(:,:),coord(:,:),nugget,range,smooth
      integer,intent(in)::covmod
      real(dp),intent(in),optional::smooth2,weights(:),rel_step
      logical,intent(in),optional::active(:)
      type(composite_se_t)::res
      real(dp),allocatable::d(:),theta(:)
      real(dp)::s2
      integer::ndep
      d=euclidean_distances(coord)
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      if(covmod==COV_CAUGEN)then
         ndep=4
         theta=[nugget,range,smooth,s2]
      else
         ndep=3
         theta=[nugget,range,smooth]
      end if
      call sandwich_dispatch(contrib,theta,res,rel_step,active)
   contains
      subroutine contrib(x,c)
         real(dp),intent(in)::x(:)
         real(dp),allocatable,intent(out)::c(:,:)
         real(dp)::jac(size(data,1),size(data,2)),ss
         real(dp),allocatable::rho(:)
         integer::i
         if(x(1)<0.0_dp.or.x(1)>=1.0_dp.or.x(2)<=0.0_dp.or.x(3)<=0.0_dp)then
            allocate(c(size(data,1),pair_count(size(data,2))))
            c=neg_huge
            return
         end if
         ss=1.0_dp
         if(ndep==4)then
            ss=x(4)
            if(ss<=0.0_dp)then
            allocate(c(size(data,1),pair_count(size(data,2))))
            c=neg_huge
            return
            end if
         end if
         allocate(rho(size(d)))
         do i=1,size(d)
            rho(i)=covariance_values(d(i),covmod,x(1),1-x(1),x(2),x(3),ss,size(coord,2))
         end do
         jac=0.0_dp
         if(present(weights))then
         c=lplik_schlather_contributions(data,rho,jac,weights)
         else
         c=lplik_schlather_contributions(data,rho,jac)
         end if
      end subroutine contrib
   end function schlather_frechet_standard_errors

   function schlather_ind_frechet_standard_errors(data,coord,covmod,alpha,nugget,range,smooth, &
      smooth2,weights,rel_step,active) result(res)
      real(dp),intent(in)::data(:,:),coord(:,:),alpha,nugget,range,smooth
      integer,intent(in)::covmod
      real(dp),intent(in),optional::smooth2,weights(:),rel_step
      logical,intent(in),optional::active(:)
      type(composite_se_t)::res
      real(dp),allocatable::d(:),theta(:)
      real(dp)::s2
      integer::ndep
      d=euclidean_distances(coord)
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      if(covmod==COV_CAUGEN)then
         ndep=5
         theta=[alpha,nugget,range,smooth,s2]
      else
         ndep=4
         theta=[alpha,nugget,range,smooth]
      end if
      call sandwich_dispatch(contrib,theta,res,rel_step,active)
   contains
      subroutine contrib(x,c)
         real(dp),intent(in)::x(:)
         real(dp),allocatable,intent(out)::c(:,:)
         real(dp)::jac(size(data,1),size(data,2)),ss
         real(dp),allocatable::rho(:)
         integer::i
         if(x(1)<0.0_dp.or.x(1)>1.0_dp.or.x(2)<0.0_dp.or.x(2)>=1.0_dp.or.x(3)<=0.0_dp.or.x(4)<=0.0_dp)then
            allocate(c(size(data,1),pair_count(size(data,2))))
            c=neg_huge
            return
         end if
         ss=1.0_dp
         if(ndep==5)then
            ss=x(5)
            if(ss<=0.0_dp)then
            allocate(c(size(data,1),pair_count(size(data,2))))
            c=neg_huge
            return
            end if
         end if
         allocate(rho(size(d)))
         do i=1,size(d)
            rho(i)=covariance_values(d(i),covmod,x(2),1-x(2),x(3),x(4),ss,size(coord,2))
         end do
         jac=0.0_dp
         if(present(weights))then
         c=lplik_schlather_ind_contributions(data,x(1),rho,jac,weights)
         else
         c=lplik_schlather_ind_contributions(data,x(1),rho,jac)
         end if
      end subroutine contrib
   end function schlather_ind_frechet_standard_errors

   function geomgauss_frechet_standard_errors(data,coord,covmod,sigma2,nugget,range,smooth, &
      smooth2,weights,rel_step,active) result(res)
      real(dp),intent(in)::data(:,:),coord(:,:),sigma2,nugget,range,smooth
      integer,intent(in)::covmod
      real(dp),intent(in),optional::smooth2,weights(:),rel_step
      logical,intent(in),optional::active(:)
      type(composite_se_t)::res
      real(dp),allocatable::d(:),theta(:)
      real(dp)::s2
      integer::ndep
      d=euclidean_distances(coord)
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      if(covmod==COV_CAUGEN)then
         ndep=5
         theta=[sigma2,nugget,range,smooth,s2]
      else
         ndep=4
         theta=[sigma2,nugget,range,smooth]
      end if
      call sandwich_dispatch(contrib,theta,res,rel_step,active)
   contains
      subroutine contrib(x,c)
         real(dp),intent(in)::x(:)
         real(dp),allocatable,intent(out)::c(:,:)
         real(dp)::jac(size(data,1),size(data,2)),ss
         real(dp),allocatable::a(:)
         integer::i
         if(x(1)<=0.0_dp.or.x(2)<0.0_dp.or.x(2)>=1.0_dp.or.x(3)<=0.0_dp.or.x(4)<=0.0_dp)then
            allocate(c(size(data,1),pair_count(size(data,2))))
            c=neg_huge
            return
         end if
         ss=1.0_dp
         if(ndep==5)then
            ss=x(5)
            if(ss<=0.0_dp)then
            allocate(c(size(data,1),pair_count(size(data,2))))
            c=neg_huge
            return
            end if
         end if
         allocate(a(size(d)))
         do i=1,size(d)
            a(i)=geom_gauss_a(d(i),covmod,x(1),x(2),x(3),x(4),ss,size(coord,2))
         end do
         jac=0.0_dp
         if(present(weights))then
         c=lplik_smith_contributions(data,a,jac,weights)
         else
         c=lplik_smith_contributions(data,a,jac)
         end if
      end subroutine contrib
   end function geomgauss_frechet_standard_errors

   function extremalt_frechet_standard_errors(data,coord,covmod,nugget,range,smooth,nu,smooth2,weights,rel_step,active) result(res)
      real(dp),intent(in)::data(:,:),coord(:,:),nugget,range,smooth,nu
      integer,intent(in)::covmod
      real(dp),intent(in),optional::smooth2,weights(:),rel_step
      logical,intent(in),optional::active(:)
      type(composite_se_t)::res
      real(dp),allocatable::d(:),theta(:)
      real(dp)::s2
      integer::ndep
      d=euclidean_distances(coord)
      s2=1.0_dp
      if(present(smooth2))s2=smooth2
      if(covmod==COV_CAUGEN)then
         ndep=5
         theta=[nugget,range,smooth,s2,nu]
      else
         ndep=4
         theta=[nugget,range,smooth,nu]
      end if
      call sandwich_dispatch(contrib,theta,res,rel_step,active)
   contains
      subroutine contrib(x,c)
         real(dp),intent(in)::x(:)
         real(dp),allocatable,intent(out)::c(:,:)
         real(dp)::jac(size(data,1),size(data,2)),ss,dof
         real(dp),allocatable::rho(:)
         integer::i
         dof=x(ndep)
         if(x(1)<0.0_dp.or.x(1)>=1.0_dp.or.x(2)<=0.0_dp.or.x(3)<=0.0_dp.or.dof<=0.0_dp)then
            allocate(c(size(data,1),pair_count(size(data,2))))
            c=neg_huge
            return
         end if
         ss=1.0_dp
         if(ndep==5)then
            ss=x(4)
            if(ss<=0.0_dp)then
            allocate(c(size(data,1),pair_count(size(data,2))))
            c=neg_huge
            return
            end if
         end if
         allocate(rho(size(d)))
         do i=1,size(d)
            rho(i)=covariance_values(d(i),covmod,x(1),1-x(1),x(2),x(3),ss,size(coord,2))
         end do
         jac=0.0_dp
         if(present(weights))then
         c=lplik_extremalt_contributions(data,rho,dof,jac,weights)
         else
         c=lplik_extremalt_contributions(data,rho,dof,jac)
         end if
      end subroutine contrib
   end function extremalt_frechet_standard_errors
end module spatialextremes_standard_errors
