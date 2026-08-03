! SPDX-License-Identifier: GPL-3.0-only
module fitheavytail_tail
   use fitheavytail_kinds, only: dp
   use fitheavytail_status, only: ht_success, ht_invalid_argument, ht_singular_matrix
   use fitheavytail_linalg, only: column_mean, standardize_columns, sample_covariance, &
      quadratic_forms, weighted_covariance, inverse_matrix, outer_product
   use fitheavytail_rng, only: random_mvt_identity
   implicit none
   private

   real(dp), parameter, public :: default_nu_min = 2.5_dp
   real(dp), parameter, public :: default_nu_max = 100.0_dp

   public :: nu_opp_estimator, nu_pop_estimator, excess_kurtosis_unbiased
   public :: nu_from_average_marginal_kurtosis, nu_from_cross_cumulants
   public :: nu_from_all_cumulants, nu_hill_estimator, nu_pareto_tail_index
   public :: cap_nu

contains

   elemental function cap_nu(nu, nu_min, nu_max) result(value)
      real(dp), intent(in) :: nu
      real(dp), intent(in), optional :: nu_min, nu_max
      real(dp) :: value, lo, hi
      lo = default_nu_min
      hi = default_nu_max
      if (present(nu_min)) lo = nu_min
      if (present(nu_max)) hi = nu_max
      value = min(hi,max(lo,nu))
   end function cap_nu

   function excess_kurtosis_unbiased(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value, centered(size(x)), m2, m4, excess
      integer :: n
      n = size(x)
      if (n < 4) then
         value = 0.0_dp
         return
      end if
      centered = x - sum(x)/real(n,dp)
      m2 = sum(centered**2)/real(n,dp)
      m4 = sum(centered**4)/real(n,dp)
      if (m2 <= tiny(1.0_dp)) then
         value = 0.0_dp
         return
      end if
      excess = real((n-1)*(n+1),dp)/real((n-2)*(n-3),dp) * &
         (m4/(m2*m2)-3.0_dp*real(n-1,dp)/real(n+1,dp))
      value = real(n-1,dp)/real((n-2)*(n-3),dp) * (real(n+1,dp)*excess+6.0_dp)
   end function excess_kurtosis_unbiased

   function nu_from_average_marginal_kurtosis(x, remove_outliers, nu_min, nu_max) result(nu)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: remove_outliers
      real(dp), intent(in), optional :: nu_min, nu_max
      real(dp) :: nu, kappa(size(x,2)), nu_col, lo, hi, ave
      logical :: remove
      integer :: j, count
      lo = default_nu_min
      hi = default_nu_max
      if (present(nu_min)) lo = nu_min
      if (present(nu_max)) hi = nu_max
      remove = .true.
      if (present(remove_outliers)) remove = remove_outliers
      count = 0
      ave = 0.0_dp
      do j = 1, size(x,2)
         kappa(j) = max(-2.0_dp/real(size(x,2)+2,dp)*0.99_dp, &
            excess_kurtosis_unbiased(x(:,j))/3.0_dp)
         nu_col = 2.0_dp/max(1.0e-10_dp,kappa(j))+4.0_dp
         if (.not.remove .or. (nu_col >= lo .and. nu_col <= hi)) then
            ave = ave + kappa(j)
            count = count + 1
         end if
      end do
      if (count == 0) then
         ave = sum(kappa)/real(size(kappa),dp)
      else
         ave = ave/real(count,dp)
      end if
      nu = cap_nu(2.0_dp/max(1.0e-10_dp,ave)+4.0_dp,lo,hi)
   end function nu_from_average_marginal_kurtosis

   function nu_from_cross_cumulants(x, nu_min, nu_max) result(nu)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: nu_min, nu_max
      real(dp) :: nu, z(size(x,1),size(x,2)), z2(size(x,1),size(x,2))
      real(dp) :: s(size(x,2),size(x,2)), s2(size(x,2),size(x,2))
      real(dp) :: sigma2(size(x,2)), kappa, ave
      integer :: i, j, count, n, t
      n = size(x,2)
      t = size(x,1)
      call standardize_columns(x,z)
      z2 = z*z
      sigma2 = sum(z2,dim=1)/real(t,dp)
      s = matmul(transpose(z),z)/real(t,dp)
      s2 = matmul(transpose(z2),z2)/real(t,dp)
      ave = 0.0_dp
      count = 0
      do i = 1, n-1
         do j = i+1, n
            kappa = s2(i,j)/(sigma2(i)*sigma2(j)+2.0_dp*s(i,j)**2)-1.0_dp
            kappa = max(-2.0_dp/real(n+2,dp)*0.99_dp,kappa)
            ave = ave+kappa
            count = count+1
         end do
      end do
      if (count > 0) ave = ave/real(count,dp)
      nu = cap_nu(2.0_dp/max(1.0e-10_dp,ave)+4.0_dp,nu_min,nu_max)
   end function nu_from_cross_cumulants

   function nu_from_all_cumulants(x, nu_min, nu_max) result(nu)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: nu_min, nu_max
      real(dp) :: nu
      nu = 0.5_dp*(nu_from_average_marginal_kurtosis(x,nu_min=nu_min,nu_max=nu_max) + &
                   nu_from_cross_cumulants(x,nu_min,nu_max))
   end function nu_from_all_cumulants

   function nu_hill_estimator(x, nu_min, nu_max) result(nu)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: nu_min, nu_max
      real(dp) :: nu, norms(size(x,1)), temp, b, denominator
      integer :: i, j, k, kn, t
      t = size(x,1)
      do i = 1, t
         norms(i) = sqrt(sum(x(i,:)**2))
      end do
      do i = 1, t-1
         k = i
         do j = i+1, t
            if (norms(j)>norms(k)) k=j
         end do
         if (k/=i) then
            temp=norms(i)
            norms(i)=norms(k)
            norms(k)=temp
         end if
      end do
      b = 4.0_dp/16.0_dp
      kn = max(1,min(t-1,int(real(t,dp)**b)))
      denominator = sum(log(norms(1:kn)/max(norms(kn+1),tiny(1.0_dp))))/real(kn,dp)
      nu = cap_nu(1.0_dp/max(denominator,tiny(1.0_dp)),nu_min,nu_max)
   end function nu_hill_estimator

   function nu_pareto_tail_index(x, center, method, nu_min, nu_max) result(nu)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: center
      character(len=*), intent(in), optional :: method
      real(dp), intent(in), optional :: nu_min, nu_max
      real(dp) :: nu, work(size(x,1),size(x,2)), mu(size(x,2)), minx, total, denom
      character(len=24) :: mode
      integer :: i,j,t,n
      logical :: do_center
      t=size(x,1)
      n=size(x,2)
      do_center=.false.
      if (present(center)) do_center=center
      work=x
      if (do_center) then
         mu=column_mean(work)
         do i=1,t
         work(i,:)=work(i,:)-mu
         end do
      end if
      work=abs(work)
      mode='WLS'
      if (present(method)) mode=adjustl(method)
      total=0.0_dp
      select case(trim(mode))
      case('MLE','MLE-unbiased','WLS')
         do j=1,n
            minx=max(minval(work(:,j)),tiny(1.0_dp))
            if (trim(mode)=='WLS') then
               denom=sum([(log(real(t,dp)/real(i,dp)),i=1,t)])/real(t,dp)
               total=total+(sum(log(max(work(:,j),tiny(1.0_dp))/minx))/real(t,dp))/denom
            else if (trim(mode)=='MLE-unbiased') then
               total=total+real(t,dp)/real(max(1,t-2),dp)* &
                  sum(log(max(work(:,j),tiny(1.0_dp))/minx))/real(t,dp)
            else
               total=total+sum(log(max(work(:,j),tiny(1.0_dp))/minx))/real(t,dp)
            end if
         end do
         total=total/real(n,dp)
      case('WLS-stacked')
         do j=1,n
            minx=max(minval(work(:,j)),tiny(1.0_dp))
            total=total+sum(log(max(work(:,j),tiny(1.0_dp))/minx))
         end do
         denom=sum([(log(real(n*t,dp)/real(i,dp)),i=1,n*t)])/real(n*t,dp)
         total=(total/real(n*t,dp))/denom
      case default
         total=huge(1.0_dp)
      end select
      nu=cap_nu(1.0_dp/max(total,tiny(1.0_dp)),nu_min,nu_max)
   end function nu_pareto_tail_index

   function nu_opp_estimator(var_x, trace_scatter, r2, method, nu_min, nu_max, status) result(nu)
      real(dp), intent(in) :: var_x(:), trace_scatter
      real(dp), intent(in), optional :: r2(:)
      character(len=*), intent(in), optional :: method
      real(dp), intent(in), optional :: nu_min, nu_max
      integer, intent(out), optional :: status
      real(dp) :: nu, eta
      character(len=24) :: mode
      integer :: n, istat
      mode='OPP'
      if (present(method)) mode=adjustl(method)
      n=size(var_x)
      istat=ht_success
      if (trace_scatter<=0.0_dp) then
         eta=1.0_dp
         istat=ht_invalid_argument
      else if (trim(mode)=='OPP-harmonic') then
         if (.not.present(r2) .or. size(r2)==0) then
            eta=1.0_dp
            istat=ht_invalid_argument
         else
            eta=sum(var_x)/trace_scatter*sum(real(n,dp)/max(r2,tiny(1.0_dp)))/real(size(r2),dp)
         end if
      else
         eta=sum(var_x)/trace_scatter
      end if
      if (abs(eta-1.0_dp)<=epsilon(1.0_dp)) then
         nu=default_nu_max
      else
         nu=cap_nu(2.0_dp*eta/(eta-1.0_dp),nu_min,nu_max)
      end if
      if (present(status)) status=istat
   end function nu_opp_estimator

   function nu_pop_estimator(r2, nu_previous, nvar, method, xc, sigma, alpha, &
      nu_min, nu_max, status) result(nu)
      real(dp), intent(in), optional :: r2(:), xc(:,:), sigma(:,:)
      real(dp), intent(in) :: nu_previous
      integer, intent(in), optional :: nvar
      character(len=*), intent(in), optional :: method
      real(dp), intent(in), optional :: alpha, nu_min, nu_max
      integer, intent(out), optional :: status
      real(dp) :: nu, theta, a, sigma_scale
      real(dp), allocatable :: rr(:), u(:), ri(:), scm(:,:), invs(:,:), sim(:,:)
      character(len=32) :: mode
      integer :: n,t,i,j,k,istat,solve_status

      istat=ht_success
      a=1.0_dp
      if(present(alpha)) a=alpha
      mode='POP-approx-2'
      if(present(method)) mode=adjustl(method)
      if(trim(mode)=='POP') mode='POP-approx-2'
      if(present(r2)) then
         t=size(r2)
         allocate(rr(t))
         rr=r2
      else if(present(xc).and.present(sigma)) then
         t=size(xc,1)
         allocate(rr(t))
         call quadratic_forms(xc,sigma,rr,solve_status)
         if(solve_status/=ht_success) istat=solve_status
      else
         allocate(rr(0))
         t=0
         istat=ht_invalid_argument
      end if
      if(present(nvar)) then
         n=nvar
      else if(present(xc)) then
         n=size(xc,2)
      else if(present(sigma)) then
         n=size(sigma,1)
      else
         n=0
         istat=ht_invalid_argument
      end if
      if(t<=n .or. n<=0) istat=ht_invalid_argument
      if(istat/=ht_success) then
         nu=cap_nu(nu_previous,nu_min,nu_max)
         if(present(status)) status=istat
         return
      end if
      allocate(u(t),ri(t))
      select case(trim(mode))
      case('POP-approx-1')
         u=(real(n,dp)+nu_previous)/(nu_previous+rr)
         ri=rr/max(1.0e-12_dp,1.0_dp-rr*u/real(t,dp))
         theta=(1.0_dp-real(n,dp)/real(t,dp))*sum(ri)/real(t*n,dp)
      case('POP-approx-2')
         u=(real(n,dp)+nu_previous)/(nu_previous+rr*real(t,dp)/real(t-1,dp))
         ri=rr/max(1.0e-12_dp,1.0_dp-rr*u/real(t,dp))
         theta=(1.0_dp-real(n,dp)/real(t,dp))*sum(ri)/real(t*n,dp)
      case('POP-approx-3')
         if(.not.present(xc)) then
            theta=1.0_dp
            istat=ht_invalid_argument
         else
            allocate(scm(n,n),invs(n,n))
            scm=matmul(transpose(xc),xc)/real(t,dp)
            call inverse_matrix(scm,invs,solve_status)
            do i=1,t
            ri(i)=dot_product(xc(i,:),matmul(invs,xc(i,:)))
            end do
            ri=ri/max(1.0e-12_dp,1.0_dp-ri/real(t,dp))
            theta=(1.0_dp-real(n,dp)/real(t,dp))*sum(ri)/real(t*n,dp)
         end if
      case('POP-approx-4')
         ri=rr/max(1.0e-12_dp,1.0_dp-rr/real(t,dp))
         theta=(1.0_dp-real(n,dp)/real(t,dp))*sum(ri)/real(t*n,dp)
      case('POP-exact')
         if(.not.present(xc)) then
            theta=1.0_dp
            istat=ht_invalid_argument
         else
            allocate(scm(n,n),invs(n,n))
            u=(real(n,dp)+nu_previous)/(nu_previous+rr)
            do i=1,t
               do k=1,3
                  scm=0.0_dp
                  do j=1,t
                     if(j/=i) scm=scm+u(j)*outer_product(xc(j,:),xc(j,:))
                  end do
                  scm=scm/(real(t,dp)*a)
                  call inverse_matrix(scm,invs,solve_status)
                  do j=1,t
                  ri(j)=dot_product(xc(j,:),matmul(invs,xc(j,:)))
                  end do
                  u=(real(n,dp)+nu_previous)/(nu_previous+ri)
               end do
               rr(i)=dot_product(xc(i,:),matmul(invs,xc(i,:)))
            end do
            theta=(1.0_dp-real(n,dp)/real(t,dp))*sum(rr)/real(t*n,dp)
         end if
      case('POP-sigma-corrected','POP-sigma-corrected-true')
         u=(real(n,dp)+nu_previous)/(nu_previous+rr)
         ri=rr/max(1.0e-12_dp,1.0_dp-rr*u/real(t,dp))
         theta=(1.0_dp-real(n,dp)/real(t,dp))*sum(ri)/real(t*n,dp)
         if(trim(mode)=='POP-sigma-corrected-true') then
            allocate(sim(10000,n),source=0.0_dp)
            call random_mvt_identity(10000,n,nu_previous,sim,24681357)
            deallocate(rr,u,ri)
            allocate(rr(10000),u(10000),ri(10000))
            do i=1,10000
            rr(i)=sum(sim(i,:)**2)
            end do
            t=10000
         end if
         sigma_scale=find_sigma_scale(rr,nu_previous,n)
         theta=theta*sigma_scale
      case default
         theta=1.0_dp
         istat=ht_invalid_argument
      end select
      if(abs(theta-1.0_dp)<=1.0e-12_dp) then
         nu=default_nu_max
      else
         nu=cap_nu(2.0_dp*theta/(theta-1.0_dp),nu_min,nu_max)
      end if
      if(present(status)) status=istat
   end function nu_pop_estimator

   function find_sigma_scale(r2,nu,n) result(root)
      real(dp),intent(in)::r2(:),nu
      integer,intent(in)::n
      real(dp)::root,lo,hi,mid,flo,fmid
      integer::iter
      lo=0.1_dp
      hi=1000.0_dp
      flo=sigma_equation(lo,r2,nu,n)
      do iter=1,100
         mid=sqrt(lo*hi)
         fmid=sigma_equation(mid,r2,nu,n)
         if(abs(fmid)<1.0e-10_dp) exit
         if(flo*fmid<=0.0_dp) then
            hi=mid
         else
            lo=mid
            flo=fmid
         end if
      end do
      root=mid
   end function find_sigma_scale

   function sigma_equation(scale,r2,nu,n) result(value)
      real(dp),intent(in)::scale,r2(:),nu
      integer,intent(in)::n
      real(dp)::value,t(size(r2))
      t=r2/scale
      value=sum((real(n,dp)+nu)/(nu+t)*t)/real(size(t),dp)-real(n,dp)
   end function sigma_equation

end module fitheavytail_tail
