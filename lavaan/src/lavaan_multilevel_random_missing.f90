module lavaan_multilevel_random_missing
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use lavaan_kinds, only : dp
   use lavaan_linalg, only : inverse_general, chol_lower, logdet_spd, solve_linear
   use lavaan_optimizer, only : bfgs_minimize
   use lavaan_multilevel_random, only : random_coefficient_result, random_effects_result
   use numderiv, only : hessian, nd_success
   implicit none
   private

   public :: random_coefficient_loglik_missing
   public :: fit_random_coefficient_missing_ml
   public :: random_effects_eb_missing

contains

   function random_coefficient_loglik_missing(y, cluster, x, beta, z, random_cov, residual_cov) result(ll)
      real(dp), intent(in) :: y(:, :), x(:, :), beta(:, :), z(:, :), random_cov(:, :), residual_cov(:, :)
      integer, intent(in) :: cluster(:)
      real(dp) :: ll
      integer, allocatable :: ids(:), rows(:), orow(:), ovar(:)
      real(dp), allocatable :: v(:, :), d(:), sol(:), mu(:, :)
      integer :: n, p, kr, qg, g, m, nobs, i, a, b, r, s, ii, jj, info
      real(dp) :: cv, ld, pi2

      n = size(y, 1)
      p = size(y, 2)
      kr = size(z, 2)
      qg = p * kr
      if (size(cluster) /= n .or. size(x, 1) /= n .or. size(beta, 1) /= size(x, 2) .or. &
          size(beta, 2) /= p .or. size(z, 1) /= n .or. any(shape(random_cov) /= [qg, qg]) .or. &
          any(shape(residual_cov) /= [p, p])) then
         ll = -huge(1.0_dp)
         return
      end if
      call unique_ids(cluster, ids)
      mu = matmul(x, beta)
      ll = 0.0_dp
      pi2 = 2.0_dp * acos(-1.0_dp)
      do g = 1, size(ids)
         rows = pack([(i, i=1,n)], cluster == ids(g))
         m = size(rows)
         nobs = 0
         do i = 1, m
            do a = 1, p
               if (ieee_is_finite(y(rows(i), a))) nobs = nobs + 1
            end do
         end do
         if (nobs == 0) cycle
         allocate(orow(nobs), ovar(nobs), d(nobs), v(nobs,nobs))
         v = 0.0_dp
         ii = 0
         do i = 1, m
            do a = 1, p
               if (.not.ieee_is_finite(y(rows(i), a))) cycle
               ii = ii + 1
               orow(ii) = rows(i)
               ovar(ii) = a
               d(ii) = y(rows(i), a) - mu(rows(i), a)
            end do
         end do
         do ii = 1, nobs
            a = ovar(ii)
            do jj = 1, nobs
               b = ovar(jj)
               cv = 0.0_dp
               do r = 1, kr
                  do s = 1, kr
                     cv = cv + z(orow(ii),r) * random_cov((a-1)*kr+r,(b-1)*kr+s) * z(orow(jj),s)
                  end do
               end do
               if (orow(ii) == orow(jj)) cv = cv + residual_cov(a,b)
               v(ii,jj) = cv
            end do
         end do
         ld = logdet_spd(v, info)
         if (info /= 0) then
         ll = -huge(1.0_dp)
         return
         end if
         call solve_linear(v, d, sol, info)
         if (info /= 0) then
         ll = -huge(1.0_dp)
         return
         end if
         ll = ll - 0.5_dp * (real(nobs,dp)*log(pi2) + ld + dot_product(d,sol))
         deallocate(orow, ovar, d, v, sol, rows)
      end do
   end function random_coefficient_loglik_missing

   subroutine random_effects_eb_missing(y, cluster, x, beta, z, random_cov, residual_cov, result)
      real(dp), intent(in) :: y(:, :), x(:, :), beta(:, :), z(:, :), random_cov(:, :), residual_cov(:, :)
      integer, intent(in) :: cluster(:)
      type(random_effects_result), intent(out) :: result
      integer, allocatable :: ids(:), rows(:), orow(:), ovar(:)
      real(dp), allocatable :: dmat(:,:), v(:,:), vinv(:,:), resid(:), mu(:,:), gd(:,:)
      integer :: n, p, kr, qg, g, m, nobs, i, a, b, r, s, ii, jj, info
      real(dp) :: cv
      n = size(y,1)
      p = size(y,2)
      kr = size(z,2)
      qg = p*kr
      if (size(cluster) /= n .or. size(x,1) /= n .or. size(beta,1) /= size(x,2) .or. &
          size(beta,2) /= p .or. size(z,1) /= n .or. any(shape(random_cov) /= [qg,qg]) .or. &
          any(shape(residual_cov) /= [p,p])) then
         result%status = -1
         return
      end if
      call unique_ids(cluster, ids)
      result%cluster_id = ids
      allocate(result%mean(size(ids),qg), result%vcov(qg,qg,size(ids)))
      result%mean = 0.0_dp
      mu = matmul(x,beta)
      do g = 1, size(ids)
         rows = pack([(i,i=1,n)], cluster == ids(g))
         m = size(rows)
         nobs = 0
         do i=1,m
         do a=1,p
         if (ieee_is_finite(y(rows(i),a))) nobs=nobs+1
         end do
         end do
         if (nobs == 0) then
            result%vcov(:,:,g) = random_cov
            deallocate(rows)
            cycle
         end if
         allocate(orow(nobs), ovar(nobs), dmat(nobs,qg), v(nobs,nobs), resid(nobs))
         dmat=0.0_dp
         v=0.0_dp
         ii=0
         do i=1,m
            do a=1,p
               if (.not.ieee_is_finite(y(rows(i),a))) cycle
               ii=ii+1
               orow(ii)=rows(i)
               ovar(ii)=a
               resid(ii)=y(rows(i),a)-mu(rows(i),a)
               do r=1,kr
               dmat(ii,(a-1)*kr+r)=z(rows(i),r)
               end do
            end do
         end do
         do ii=1,nobs
            a=ovar(ii)
            do jj=1,nobs
               b=ovar(jj)
               cv=0.0_dp
               do r=1,kr
               do s=1,kr
                  cv=cv+z(orow(ii),r)*random_cov((a-1)*kr+r,(b-1)*kr+s)*z(orow(jj),s)
               end do
               end do
               if (orow(ii)==orow(jj)) cv=cv+residual_cov(a,b)
               v(ii,jj)=cv
            end do
         end do
         call inverse_general(v,vinv,info)
         if(info/=0) then
         result%status=100+info
         return
         end if
         gd=matmul(random_cov,transpose(dmat))
         result%mean(g,:)=matmul(gd,matmul(vinv,resid))
         result%vcov(:,:,g)=random_cov-matmul(gd,matmul(vinv,transpose(gd)))
         result%vcov(:,:,g)=0.5_dp*(result%vcov(:,:,g)+transpose(result%vcov(:,:,g)))
         deallocate(orow,ovar,dmat,v,resid,vinv,gd,rows)
      end do
      result%status=0
   end subroutine random_effects_eb_missing

   subroutine fit_random_coefficient_missing_ml(y, cluster, x, z, result, compute_se)
      real(dp), intent(in) :: y(:, :), x(:, :), z(:, :)
      integer, intent(in) :: cluster(:)
      type(random_coefficient_result), intent(out) :: result
      logical, intent(in), optional :: compute_se
      real(dp), allocatable :: beta0(:,:), g0(:,:), r0(:,:), par(:), hess(:,:), hinv(:,:)
      real(dp), allocatable :: beta(:,:), gcov(:,:), rcov(:,:), xtx(:,:), xty(:), inv(:,:), res(:)
      integer, allocatable :: ids(:), keep(:)
      integer :: n,p,kf,kr,qg,npar,a,i,info,hs,nk
      real(dp) :: fval, vv
      logical :: dose
      n=size(y,1)
      p=size(y,2)
      kf=size(x,2)
      kr=size(z,2)
      qg=p*kr
      dose=.true.
      if(present(compute_se)) dose=compute_se
      if(size(cluster)/=n .or. size(x,1)/=n .or. size(z,1)/=n .or. n<kf+3) then
         result%status=-1
         return
      end if
      allocate(beta0(kf,p))
      beta0=0.0_dp
      do a=1,p
         keep=pack([(i,i=1,n)], [(ieee_is_finite(y(i,a)),i=1,n)])
         nk=size(keep)
         if(nk>=kf) then
            xtx=matmul(transpose(x(keep,:)),x(keep,:))
            call inverse_general(xtx,inv,info)
            if(info==0) then
               xty=matmul(transpose(x(keep,:)),y(keep,a))
               beta0(:,a)=matmul(inv,xty)
            end if
         end if
         if(allocated(keep)) deallocate(keep)
      end do
      allocate(r0(p,p))
      r0=0.0_dp
      do a=1,p
         keep=pack([(i,i=1,n)], [(ieee_is_finite(y(i,a)),i=1,n)])
         if(size(keep)>1) then
            res=y(keep,a)-matmul(x(keep,:),beta0(:,a))
            vv=sum(res*res)/real(max(1,size(keep)-kf),dp)
         else
         vv=1.0_dp
         end if
         r0(a,a)=max(vv,1.0e-4_dp)
         if(allocated(keep)) deallocate(keep)
      end do
      allocate(g0(qg,qg))
      g0=0.0_dp
      do i=1,qg
      a=1+(i-1)/kr
      g0(i,i)=0.05_dp*r0(a,a)
      end do
      npar=kf*p+qg*(qg+1)/2+p*(p+1)/2
      allocate(par(npar))
      call pack_parameters(beta0,g0,r0,par,info)
      if(info/=0) then
      result%status=100+info
      return
      end if
      call bfgs_minimize(nll,par,fval,result%converged,result%iterations,maxiter=900,tol=2.0e-6_dp)
      call unpack_parameters(par,kf,p,kr,beta,gcov,rcov)
      result%beta=beta
      result%random_cov=gcov
      result%residual_cov=rcov
      result%par=par
      result%loglik=-fval
      result%aic=2.0_dp*fval+2.0_dp*real(npar,dp)
      result%bic=2.0_dp*fval+log(real(n,dp))*real(npar,dp)
      call unique_ids(cluster,ids)
      result%ncluster=size(ids)
      allocate(result%vcov(npar,npar),result%se(npar))
      result%vcov=0.0_dp
      result%se=huge(1.0_dp)
      if(dose) then
         call hessian(nll,par,hess,status=hs)
         if(hs==nd_success) then
            call inverse_general(hess,hinv,info)
            if(info==0) then
               result%vcov=hinv
               do i=1,npar
               if(hinv(i,i)>=0.0_dp) result%se(i)=sqrt(hinv(i,i))
               end do
            end if
         end if
      end if
      result%status=0
   contains
      function nll(v) result(f)
         real(dp),intent(in)::v(:)
         real(dp)::f
         real(dp),allocatable::bb(:,:),gg(:,:),rr(:,:)
         call unpack_parameters(v,kf,p,kr,bb,gg,rr)
         f=-random_coefficient_loglik_missing(y,cluster,x,bb,z,gg,rr)
         if(.not.(f<huge(1.0_dp)/10.0_dp)) f=huge(1.0_dp)/100.0_dp
      end function nll
   end subroutine fit_random_coefficient_missing_ml

   subroutine pack_parameters(beta,gcov,rcov,par,info)
      real(dp),intent(in)::beta(:,:),gcov(:,:),rcov(:,:)
      real(dp),intent(out)::par(:)
      integer,intent(out)::info
      real(dp),allocatable::lg(:,:),lr(:,:)
      integer::i,j,pos
      call chol_lower(gcov,lg,info)
      if(info/=0)return
      call chol_lower(rcov,lr,info)
      if(info/=0)return
      pos=0
      do j=1,size(beta,2)
      do i=1,size(beta,1)
      pos=pos+1
      par(pos)=beta(i,j)
      end do
      end do
      do j=1,size(lg,1)
      do i=j,size(lg,1)
      pos=pos+1
         if(i==j) then
         par(pos)=log(max(lg(i,j),1.0e-8_dp))
         else
         par(pos)=lg(i,j)
         end if
      end do
      end do
      do j=1,size(lr,1)
      do i=j,size(lr,1)
      pos=pos+1
         if(i==j) then
         par(pos)=log(max(lr(i,j),1.0e-8_dp))
         else
         par(pos)=lr(i,j)
         end if
      end do
      end do
      info=0
   end subroutine pack_parameters

   subroutine unpack_parameters(par,kf,p,kr,beta,gcov,rcov)
      real(dp),intent(in)::par(:)
      integer,intent(in)::kf,p,kr
      real(dp),allocatable,intent(out)::beta(:,:),gcov(:,:),rcov(:,:)
      real(dp),allocatable::lg(:,:),lr(:,:)
      integer::qg,i,j,pos
      qg=p*kr
      allocate(beta(kf,p),lg(qg,qg),lr(p,p))
      lg=0.0_dp
      lr=0.0_dp
      pos=0
      do j=1,p
      do i=1,kf
      pos=pos+1
      beta(i,j)=par(pos)
      end do
      end do
      do j=1,qg
      do i=j,qg
      pos=pos+1
         if(i==j) then
         lg(i,j)=exp(par(pos))
         else
         lg(i,j)=par(pos)
         end if
      end do
      end do
      do j=1,p
      do i=j,p
      pos=pos+1
         if(i==j) then
         lr(i,j)=exp(par(pos))
         else
         lr(i,j)=par(pos)
         end if
      end do
      end do
      gcov=matmul(lg,transpose(lg))
      rcov=matmul(lr,transpose(lr))
   end subroutine unpack_parameters

   subroutine unique_ids(cluster,ids)
      integer,intent(in)::cluster(:)
      integer,allocatable,intent(out)::ids(:)
      integer,allocatable::tmp(:)
      integer::i,n
      allocate(tmp(size(cluster)))
      n=0
      do i=1,size(cluster)
         if(n==0 .or. .not.any(tmp(1:n)==cluster(i))) then
         n=n+1
         tmp(n)=cluster(i)
         end if
      end do
      allocate(ids(n))
      if(n>0) ids=tmp(1:n)
   end subroutine unique_ids

end module lavaan_multilevel_random_missing
