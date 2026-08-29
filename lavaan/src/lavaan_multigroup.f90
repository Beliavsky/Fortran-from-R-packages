module lavaan_multigroup
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model, ram_free_map, ram_get_free, ram_set_free, ram_sigma, ram_mu
   use lavaan_objectives, only : objective_ml, mvn_loglik_complete, mvn_loglik_missing
   use lavaan_linalg, only : inverse_general, logdet_spd, inverse_spd
   use lavaan_optimizer, only : bfgs_minimize
   use numderiv, only : hessian, nd_success
   implicit none
   private

   type, public :: ram_group_spec
      type(ram_model) :: model
      type(ram_free_map) :: map
      integer, allocatable :: link(:)
   end type ram_group_spec

   type, public :: ram_group_data
      real(dp), allocatable :: x(:, :)
   end type ram_group_data

   type, public :: sem_multigroup_result
      real(dp), allocatable :: par(:), se(:), vcov(:, :)
      real(dp), allocatable :: sigma(:, :, :), mu(:, :), group_chisq(:), group_loglik(:)
      real(dp) :: objective = huge(1.0_dp)
      real(dp) :: chisq = huge(1.0_dp)
      real(dp) :: df = 0.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp), bic = huge(1.0_dp)
      logical :: converged = .false.
      integer :: iterations = 0, status = 0
   end type sem_multigroup_result

   public :: independent_group_links, fit_ram_multigroup_cov, fit_ram_multigroup_fiml

contains

   subroutine independent_group_links(groups)
      type(ram_group_spec), intent(inout) :: groups(:)
      integer :: g, k, next
      next = 0
      do g = 1, size(groups)
         if (allocated(groups(g)%link)) deallocate(groups(g)%link)
         allocate(groups(g)%link(size(groups(g)%map%matrix_id)))
         do k = 1, size(groups(g)%link)
            next = next + 1
            groups(g)%link(k) = next
         end do
      end do
   end subroutine independent_group_links

   subroutine fit_ram_multigroup_cov(groups, data_cov, data_mean, nobs, result)
      type(ram_group_spec), intent(in) :: groups(:)
      real(dp), intent(in) :: data_cov(:, :, :), data_mean(:, :)
      integer, intent(in) :: nobs(:)
      type(sem_multigroup_result), intent(out) :: result
      type(ram_model), allocatable :: work(:)
      real(dp), allocatable :: x(:), sumx(:), hess(:, :), hi(:, :), sg(:, :), mm(:), si(:, :)
      integer, allocatable :: cnt(:)
      real(dp) :: fval, f, ld, q
      integer :: g, k, j, nglobal, p, info, status, ntot, nstat
      logical :: meanstructure

      if (size(groups) /= size(nobs) .or. size(groups) /= size(data_cov, 3) .or. &
          size(groups) /= size(data_mean, 2)) then
         result%status = -1
         return
      end if
      if (size(groups) == 0) then
         result%status = -2
         return
      end if
      p = size(data_cov, 1)
      if (size(data_cov, 2) /= p .or. size(data_mean, 1) /= p) then
         result%status = -3
         return
      end if
      nglobal = 0
      do g = 1, size(groups)
         if (.not. allocated(groups(g)%link)) then
            result%status = -4
            return
         end if
         if (size(groups(g)%link) /= size(groups(g)%map%matrix_id)) then
            result%status = -5
            return
         end if
         if (size(groups(g)%link) > 0) nglobal = max(nglobal, maxval(groups(g)%link))
      end do
      if (nglobal <= 0) then
         result%status = -6
         return
      end if

      allocate(work(size(groups)), x(nglobal), sumx(nglobal), cnt(nglobal))
      sumx = 0.0_dp
      cnt = 0
      do g = 1, size(groups)
         work(g) = groups(g)%model
         block
            real(dp), allocatable :: gx(:)
            gx = ram_get_free(groups(g)%model, groups(g)%map)
            do k = 1, size(gx)
               j = groups(g)%link(k)
               if (j < 1 .or. j > nglobal) then
                  result%status = -7
                  return
               end if
               sumx(j) = sumx(j) + gx(k)
               cnt(j) = cnt(j) + 1
            end do
         end block
      end do
      do j = 1, nglobal
         if (cnt(j) == 0) then
            result%status = -8
            return
         end if
         x(j) = sumx(j) / real(cnt(j), dp)
      end do

      call bfgs_minimize(total_nll, x, fval, result%converged, result%iterations, maxiter=1600, tol=1.0e-7_dp)
      result%par = x
      result%objective = 2.0_dp*fval / real(sum(nobs), dp)
      ntot = sum(nobs)
      result%loglik = -fval
      result%chisq = 0.0_dp
      allocate(result%sigma(p,p,size(groups)), result%mu(p,size(groups)))
      allocate(result%group_chisq(size(groups)), result%group_loglik(size(groups)))
      result%group_chisq = 0.0_dp
      result%group_loglik = 0.0_dp
      meanstructure = .false.
      do g = 1, size(groups)
         call apply_group(g, x, work(g))
         call ram_sigma(work(g), sg, info)
         if (info /= 0) then
         result%status = info
         return
         end if
         call ram_mu(work(g), mm, info)
         result%sigma(:,:,g) = sg
         result%mu(:,g) = mm
         meanstructure = meanstructure .or. allocated(groups(g)%model%m)
         f = objective_ml(sg, mm, data_cov(:,:,g), data_mean(:,g), allocated(groups(g)%model%m), info)
         result%group_chisq(g) = real(nobs(g),dp)*f
         result%chisq = result%chisq + result%group_chisq(g)
         ld = logdet_spd(sg, info)
         call inverse_spd(sg, si, info)
         q = sum(data_cov(:,:,g)*si)
         if (allocated(groups(g)%model%m)) then
            q = q + dot_product(data_mean(:,g)-mm, matmul(si, data_mean(:,g)-mm))
         end if
         result%group_loglik(g) = -0.5_dp*real(nobs(g),dp) * &
            (real(p,dp)*log(2.0_dp*acos(-1.0_dp)) + ld + q)
      end do
      result%loglik = sum(result%group_loglik)
      nstat = 0
      do g = 1, size(groups)
         nstat = nstat + p*(p+1)/2 + merge(p,0,allocated(groups(g)%model%m))
      end do
      result%df = real(nstat-nglobal,dp)
      result%aic = -2.0_dp*result%loglik + 2.0_dp*real(nglobal,dp)
      result%bic = -2.0_dp*result%loglik + log(real(ntot,dp))*real(nglobal,dp)

      call hessian(total_nll, x, hess, status=status)
      allocate(result%vcov(nglobal,nglobal), result%se(nglobal))
      result%vcov = 0.0_dp
      result%se = huge(1.0_dp)
      if (status == nd_success) then
         call inverse_general(hess, hi, info)
         if (info == 0) then
            result%vcov = hi
            do j = 1, nglobal
               if (hi(j,j) >= 0.0_dp) result%se(j) = sqrt(hi(j,j))
            end do
         end if
      end if
      result%status = 0

   contains
      subroutine apply_group(gg, z, model)
         integer, intent(in) :: gg
         real(dp), intent(in) :: z(:)
         type(ram_model), intent(inout) :: model
         real(dp), allocatable :: gx(:)
         integer :: kk
         allocate(gx(size(groups(gg)%link)))
         do kk = 1, size(gx)
            gx(kk) = z(groups(gg)%link(kk))
         end do
         model = groups(gg)%model
         call ram_set_free(model, groups(gg)%map, gx)
      end subroutine apply_group

      function total_nll(z) result(v)
         real(dp), intent(in) :: z(:)
         real(dp) :: v, ff
         real(dp), allocatable :: sgg(:,:), mmm(:)
         integer :: gg, istat
         v = 0.0_dp
         do gg = 1, size(groups)
            call apply_group(gg, z, work(gg))
            call ram_sigma(work(gg), sgg, istat)
            if (istat /= 0 .or. any([(sgg(j,j) <= 0.0_dp, j=1,size(sgg,1))])) then
               v = huge(1.0_dp)/100.0_dp
               return
            end if
            call ram_mu(work(gg), mmm, istat)
            ff = objective_ml(sgg, mmm, data_cov(:,:,gg), data_mean(:,gg), &
               allocated(groups(gg)%model%m), istat)
            if (istat /= 0 .or. ff >= huge(1.0_dp)/1000.0_dp) then
               v = huge(1.0_dp)/100.0_dp
               return
            end if
            v = v + 0.5_dp*real(nobs(gg),dp)*ff
         end do
      end function total_nll
   end subroutine fit_ram_multigroup_cov


   subroutine fit_ram_multigroup_fiml(groups,datasets,result)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
      type(ram_group_spec),intent(in)::groups(:)
      type(ram_group_data),intent(in)::datasets(:)
      type(sem_multigroup_result),intent(out)::result
      type(ram_model),allocatable::work(:)
      real(dp),allocatable::x(:),sumx(:),hess(:,:),hi(:,:),sg(:,:),mm(:)
      integer,allocatable::cnt(:)
      integer::g,k,j,nglobal,p,info,status,ntot,nstat
      real(dp)::fval,ll

      if(size(groups)/=size(datasets) .or. size(groups)==0) then
      result%status=-20
      return
      end if
      p=size(datasets(1)%x,2)
      nglobal=0
      ntot=0
      do g=1,size(groups)
         if(.not.allocated(datasets(g)%x) .or. size(datasets(g)%x,2)/=p) then
         result%status=-21
         return
         end if
         if(.not.allocated(groups(g)%link)) then
         result%status=-22
         return
         end if
         if(size(groups(g)%link)/=size(groups(g)%map%matrix_id)) then
         result%status=-23
         return
         end if
         if(size(groups(g)%link)>0) nglobal=max(nglobal,maxval(groups(g)%link))
         ntot=ntot+size(datasets(g)%x,1)
      end do
      if(nglobal<=0 .or. ntot<=0) then
      result%status=-24
      return
      end if
      allocate(work(size(groups)),x(nglobal),sumx(nglobal),cnt(nglobal))
      sumx=0.0_dp
      cnt=0
      do g=1,size(groups)
         work(g)=groups(g)%model
         block
            real(dp),allocatable::gx(:)
            gx=ram_get_free(groups(g)%model,groups(g)%map)
            do k=1,size(gx)
               j=groups(g)%link(k)
               sumx(j)=sumx(j)+gx(k)
               cnt(j)=cnt(j)+1
            end do
         end block
      end do
      do j=1,nglobal
         if(cnt(j)==0) then
         result%status=-25
         return
         end if
         x(j)=sumx(j)/real(cnt(j),dp)
      end do
      call bfgs_minimize(total_nll,x,fval,result%converged,result%iterations,maxiter=1800,tol=1.0e-7_dp)
      result%par=x
      result%objective=2.0_dp*fval/real(ntot,dp)
      result%loglik=-fval
      allocate(result%sigma(p,p,size(groups)),result%mu(p,size(groups)),result%group_loglik(size(groups)))
      allocate(result%group_chisq(size(groups)))
      result%group_chisq=0.0_dp
      do g=1,size(groups)
         call apply_group_f(g,x,work(g))
         call ram_sigma(work(g),sg,info)
         call ram_mu(work(g),mm,info)
         if(info/=0) then
         result%status=info
         return
         end if
         result%sigma(:,:,g)=sg
         result%mu(:,g)=mm
         if(any(ieee_is_nan(datasets(g)%x))) then
            ll=mvn_loglik_missing(datasets(g)%x,mm,sg,info)
         else
            ll=mvn_loglik_complete(datasets(g)%x,mm,sg,info)
         end if
         result%group_loglik(g)=ll
      end do
      result%loglik=sum(result%group_loglik)
      nstat=0
      do g=1,size(groups)
         nstat=nstat+p*(p+1)/2+merge(p,0,allocated(groups(g)%model%m))
      end do
      result%df=real(nstat-nglobal,dp)
      result%chisq=huge(1.0_dp)
      result%aic=-2.0_dp*result%loglik+2.0_dp*real(nglobal,dp)
      result%bic=-2.0_dp*result%loglik+log(real(ntot,dp))*real(nglobal,dp)
      call hessian(total_nll,x,hess,status=status)
      allocate(result%vcov(nglobal,nglobal),result%se(nglobal))
      result%vcov=0.0_dp
      result%se=huge(1.0_dp)
      if(status==nd_success) then
         call inverse_general(hess,hi,info)
         if(info==0) then
            result%vcov=hi
            do j=1,nglobal
            if(hi(j,j)>=0.0_dp) result%se(j)=sqrt(hi(j,j))
            end do
         end if
      end if
      result%status=0
   contains
      subroutine apply_group_f(gg,z,model)
         integer,intent(in)::gg
         real(dp),intent(in)::z(:)
         type(ram_model),intent(inout)::model
         real(dp),allocatable::gx(:)
         integer::kk
         allocate(gx(size(groups(gg)%link)))
         do kk=1,size(gx)
         gx(kk)=z(groups(gg)%link(kk))
         end do
         model=groups(gg)%model
         call ram_set_free(model,groups(gg)%map,gx)
      end subroutine apply_group_f
      function total_nll(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v,lll
         real(dp),allocatable::ss(:,:),muu(:)
         integer::gg,istat
         v=0.0_dp
         do gg=1,size(groups)
            call apply_group_f(gg,z,work(gg))
            call ram_sigma(work(gg),ss,istat)
            call ram_mu(work(gg),muu,istat)
            if(istat/=0) then
            v=huge(1.0_dp)/100.0_dp
            return
            end if
            if(any(ieee_is_nan(datasets(gg)%x))) then
               lll=mvn_loglik_missing(datasets(gg)%x,muu,ss,istat)
            else
               lll=mvn_loglik_complete(datasets(gg)%x,muu,ss,istat)
            end if
            if(istat/=0) then
            v=huge(1.0_dp)/100.0_dp
            return
            end if
            v=v-lll
         end do
      end function total_nll
   end subroutine fit_ram_multigroup_fiml

end module lavaan_multigroup
