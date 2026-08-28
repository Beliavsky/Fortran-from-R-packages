module lavaan_constraints
   use lavaan_kinds, only : dp
   use lavaan_ram, only : ram_model, ram_free_map, ram_get_free, ram_set_free, ram_sigma, ram_mu
   use lavaan_objectives, only : objective_ml
   use lavaan_optimizer, only : bfgs_minimize
   use lavaan_linalg, only : inverse_general, inverse_spd, logdet_spd
   use lavaan_fit, only : sem_fit_result
   use numderiv, only : hessian, nd_success
   implicit none
   private

   abstract interface
      subroutine constraint_callback(x, equality, inequality)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), allocatable, intent(out) :: equality(:), inequality(:)
      end subroutine constraint_callback
   end interface

   public :: constraint_callback, fit_ram_cov_constrained

contains

   subroutine fit_ram_cov_constrained(template, map, data_cov, data_mean, nobs, constraints, result, lower, upper)
      type(ram_model), intent(in) :: template
      type(ram_free_map), intent(in) :: map
      real(dp), intent(in) :: data_cov(:,:), data_mean(:)
      integer, intent(in) :: nobs
      procedure(constraint_callback) :: constraints
      type(sem_fit_result), intent(out) :: result
      real(dp), intent(in), optional :: lower(:), upper(:)
      type(ram_model) :: work
      real(dp), allocatable :: x(:), ceq(:), cineq(:), hess(:,:), hinv(:,:), jac(:,:), middle(:,:), mi(:,:)
      real(dp), allocatable :: si(:,:), d(:)
      real(dp) :: fval, rho, maxviol, ld, q
      integer :: outer, info, status, i, k, p, neq
      logical :: conv

      x = ram_get_free(template, map)
      k = size(x)
      p = size(data_cov,1)
      if (present(lower)) then
         if (size(lower) /= k) then
         result%status=-20
         return
         end if
      end if
      if (present(upper)) then
         if (size(upper) /= k) then
         result%status=-21
         return
         end if
      end if
      rho = 1.0e2_dp
      conv = .false.
      do outer = 1, 8
         call bfgs_minimize(penalized, x, fval, conv, result%iterations, maxiter=800, tol=1.0e-7_dp)
         call constraints(x, ceq, cineq)
         maxviol = 0.0_dp
         if (size(ceq) > 0) maxviol = max(maxviol, maxval(abs(ceq)))
         if (size(cineq) > 0) maxviol = max(maxviol, maxval(max(cineq,0.0_dp)))
         if (present(lower)) maxviol = max(maxviol, maxval(max(lower-x,0.0_dp)))
         if (present(upper)) maxviol = max(maxviol, maxval(max(x-upper,0.0_dp)))
         if (maxviol < 5.0e-7_dp .and. conv) exit
         rho = rho*10.0_dp
      end do
      result%converged = conv .and. maxviol < 1.0e-5_dp
      result%par = x
      work = template
      call ram_set_free(work,map,x)
      call ram_sigma(work,result%sigma,info)
      if(info/=0) then
      result%status=info
      return
      end if
      call ram_mu(work,result%mu,info)
      result%objective = base_objective(x)
      result%chisq = real(nobs,dp)*result%objective
      call inverse_spd(result%sigma,si,info)
      ld=logdet_spd(result%sigma,info)
      q=sum(data_cov*si)
      if(allocated(template%m)) then
         d=data_mean-result%mu
         q=q+dot_product(d,matmul(si,d))
      end if
      result%loglik=-0.5_dp*real(nobs,dp)*(real(p,dp)*log(2.0_dp*acos(-1.0_dp))+ld+q)
      call constraints(x,ceq,cineq)
      neq=size(ceq)
      result%df=real(p*(p+1)/2 + merge(p,0,allocated(template%m)) - k + neq,dp)
      result%aic=-2.0_dp*result%loglik+2.0_dp*real(k-neq,dp)
      result%bic=-2.0_dp*result%loglik+log(real(nobs,dp))*real(k-neq,dp)

      call hessian(total_nll,x,hess,status=status)
      allocate(result%vcov(k,k),result%se(k))
      result%vcov=0.0_dp
      result%se=huge(1.0_dp)
      if(status==nd_success) then
         call inverse_general(hess,hinv,info)
         if(info==0) then
            if(neq>0) then
               call constraint_jacobian(x,jac)
               middle=matmul(jac,matmul(hinv,transpose(jac)))
               call inverse_general(middle,mi,info)
               if(info==0) then
                  result%vcov=hinv-matmul(hinv,matmul(transpose(jac),matmul(mi,matmul(jac,hinv))))
               else
                  result%vcov=hinv
               end if
            else
               result%vcov=hinv
            end if
            do i=1,k
               if(result%vcov(i,i)>=0.0_dp) result%se(i)=sqrt(result%vcov(i,i))
            end do
         end if
      end if
      result%status=merge(0,1,result%converged)

   contains
      function base_objective(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v
         real(dp),allocatable::s(:,:),m(:)
         integer::istat,jj
         work=template
         call ram_set_free(work,map,z)
         call ram_sigma(work,s,istat)
         if(istat/=0 .or. any([(s(jj,jj)<=0.0_dp,jj=1,size(s,1))])) then
            v=huge(1.0_dp)/100.0_dp
            return
         end if
         call ram_mu(work,m,istat)
         v=objective_ml(s,m,data_cov,data_mean,allocated(template%m),istat)
         if(istat/=0) v=huge(1.0_dp)/100.0_dp
      end function base_objective

      function total_nll(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v
         v=0.5_dp*real(nobs,dp)*base_objective(z)
      end function total_nll

      function penalized(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v,pen
         real(dp),allocatable::eq(:),ineq(:)
         call constraints(z,eq,ineq)
         pen=sum(eq*eq)+sum(max(ineq,0.0_dp)**2)
         if(present(lower)) pen=pen+sum(max(lower-z,0.0_dp)**2)
         if(present(upper)) pen=pen+sum(max(z-upper,0.0_dp)**2)
         v=total_nll(z)+rho*pen
      end function penalized

      subroutine constraint_jacobian(z,j)
         real(dp),intent(in)::z(:)
         real(dp),allocatable,intent(out)::j(:,:)
         real(dp),allocatable::zp(:),zm(:),ep(:),em(:),dummy(:)
         real(dp)::h
         integer::ii
         call constraints(z,ep,dummy)
         allocate(j(size(ep),size(z)),zp(size(z)),zm(size(z)))
         do ii=1,size(z)
            h=1.0e-5_dp*max(1.0_dp,abs(z(ii)))
            zp=z
            zm=z
            zp(ii)=zp(ii)+h
            zm(ii)=zm(ii)-h
            call constraints(zp,ep,dummy)
            call constraints(zm,em,dummy)
            j(:,ii)=(ep-em)/(2.0_dp*h)
         end do
      end subroutine constraint_jacobian
   end subroutine fit_ram_cov_constrained
end module lavaan_constraints
