! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_constraints
   use vgam_kinds, only : dp
   use vgam_links, only : link_value, link_inverse, link_derivative
   use vgam_vglm, only : variance_function, default_link
   use vgam_linalg, only : weighted_least_squares
   implicit none
   private

   type, public :: constrained_vglm_result_t
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: free_coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: linear_predictor(:, :)
      integer, allocatable :: families(:)
      integer, allocatable :: links(:)
      real(dp) :: deviance = huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict => predict_constrained_vglm
   end type constrained_vglm_result_t

   public :: fit_constrained_vglm, identity_constraint, parallel_constraint

contains

   subroutine identity_constraint(p, m, constraint)
      integer, intent(in) :: p, m
      real(dp), allocatable, intent(out) :: constraint(:, :)
      integer :: i
      allocate(constraint(p*m,p*m))
      constraint = 0.0_dp
      do i = 1, p*m
         constraint(i,i) = 1.0_dp
      end do
   end subroutine identity_constraint

   subroutine parallel_constraint(p, m, parallel_columns, constraint)
      integer, intent(in) :: p, m
      logical, intent(in) :: parallel_columns(:)
      real(dp), allocatable, intent(out) :: constraint(:, :)
      integer :: j, r, q, col, idx

      if (size(parallel_columns) /= p) error stop "parallel_constraint: size mismatch"
      q = count(parallel_columns) + count(.not.parallel_columns)*m
      allocate(constraint(p*m,q))
      constraint = 0.0_dp
      col = 0
      do j = 1, p
         if (parallel_columns(j)) then
            col = col + 1
            do r = 1, m
               idx = (r-1)*p + j
               constraint(idx,col) = 1.0_dp
            end do
         else
            do r = 1, m
               col = col + 1
               idx = (r-1)*p + j
               constraint(idx,col) = 1.0_dp
            end do
         end if
      end do
   end subroutine parallel_constraint

   subroutine fit_constrained_vglm(y, x, families, constraint, result, links, weights, offsets, max_iter, tol)
      real(dp), intent(in) :: y(:,:), x(:,:), constraint(:,:)
      integer, intent(in) :: families(:)
      type(constrained_vglm_result_t), intent(out) :: result
      integer, intent(in), optional :: links(:), max_iter
      real(dp), intent(in), optional :: weights(:,:), offsets(:,:), tol
      real(dp), allocatable :: w(:,:), off(:,:), eta(:,:), mu(:,:), z(:,:), ww(:,:)
      real(dp), allocatable :: design(:,:), response(:), stack_w(:), theta(:), cov(:,:), full_beta(:)
      real(dp), allocatable :: old_theta(:)
      integer, allocatable :: link_ids(:)
      real(dp) :: tolerance, dlink, var, chg, dev, devold
      integer :: n,p,m,q,niter,iter,i,j,row,stat

      n=size(y,1); m=size(y,2); p=size(x,2); q=size(constraint,2)
      if(n<=0.or.m<=0.or.p<=0.or.size(x,1)/=n.or.size(families)/=m.or. &
         size(constraint,1)/=p*m.or.q<=0)then
         result%status=1; return
      end if
      allocate(link_ids(m))
      do j=1,m
         link_ids(j)=default_link(families(j))
      end do
      if(present(links))then
         if(size(links)/=m)then
            result%status=2; return
         end if
         link_ids=links
      end if
      if(present(weights))then
         if(any(shape(weights)/=shape(y)).or.any(weights<0.0_dp))then
            result%status=3; return
         end if
         w=weights
      else
         allocate(w(n,m)); w=1.0_dp
      end if
      if(present(offsets))then
         if(any(shape(offsets)/=shape(y)))then
            result%status=4; return
         end if
         off=offsets
      else
         allocate(off(n,m)); off=0.0_dp
      end if
      niter=100; if(present(max_iter))niter=max_iter
      tolerance=1.0e-8_dp; if(present(tol))tolerance=tol
      allocate(eta(n,m),mu(n,m),z(n,m),ww(n,m))
      do j=1,m
         do i=1,n
            mu(i,j)=initial_mean(y(i,j),families(j))
            eta(i,j)=link_value(mu(i,j),link_ids(j))
            z(i,j)=eta(i,j)-off(i,j)
            ww(i,j)=max(w(i,j),tiny(1.0_dp))
         end do
      end do
      call build_stacked_design(x,constraint,design)
      allocate(response(n*m),stack_w(n*m),old_theta(q))
      call stack_columns(z,response)
      call stack_columns(ww,stack_w)
      call weighted_least_squares(design,response,stack_w,theta,cov,stat)
      if(stat/=0)then
         result%status=10+stat; return
      end if
      devold=huge(1.0_dp)
      do iter=1,niter
         old_theta=theta
         full_beta=matmul(constraint,theta)
         eta=matmul(x,reshape(full_beta,[p,m]))+off
         do j=1,m
            do i=1,n
               mu(i,j)=link_inverse(eta(i,j),link_ids(j))
               call clamp_mu(mu(i,j),families(j))
               dlink=link_derivative(mu(i,j),link_ids(j))
               var=variance_function(mu(i,j),families(j))
               ww(i,j)=w(i,j)/max(var*dlink*dlink,tiny(1.0_dp))
               z(i,j)=eta(i,j)+(y(i,j)-mu(i,j))*dlink-off(i,j)
            end do
         end do
         call stack_columns(z,response)
         call stack_columns(ww,stack_w)
         call weighted_least_squares(design,response,stack_w,theta,cov,stat)
         if(stat/=0)then
            result%status=20+stat; return
         end if
         chg=maxval(abs(theta-old_theta))/max(1.0_dp,maxval(abs(old_theta)))
         full_beta=matmul(constraint,theta)
         eta=matmul(x,reshape(full_beta,[p,m]))+off
         do j=1,m
            do i=1,n
               mu(i,j)=link_inverse(eta(i,j),link_ids(j))
               call clamp_mu(mu(i,j),families(j))
            end do
         end do
         dev=independent_deviance(y,mu,w,families)
         if(chg<=tolerance.or.abs(dev-devold)<=tolerance*(1.0_dp+abs(devold)))then
            result%converged=.true.; exit
         end if
         devold=dev
      end do
      result%iterations=min(iter,niter)
      result%free_coefficients=theta
      full_beta=matmul(constraint,theta)
      result%coefficients=reshape(full_beta,[p,m])
      result%covariance=cov
      result%fitted=mu
      result%linear_predictor=eta
      result%deviance=independent_deviance(y,mu,w,families)
      result%families=families
      result%links=link_ids
      if(.not.result%converged.and.result%status==0)result%status=100
   end subroutine fit_constrained_vglm

   subroutine build_stacked_design(x,constraint,design)
      real(dp),intent(in)::x(:,:),constraint(:,:)
      real(dp),allocatable,intent(out)::design(:,:)
      real(dp),allocatable::full(:,:)
      integer::n,p,m,q,i,j,row
      n=size(x,1); p=size(x,2); m=size(constraint,1)/p; q=size(constraint,2)
      allocate(full(n*m,p*m)); full=0.0_dp
      row=0
      do j=1,m
         do i=1,n
            row=row+1
            full(row,(j-1)*p+1:j*p)=x(i,:)
         end do
      end do
      design=matmul(full,constraint)
   end subroutine build_stacked_design

   subroutine stack_columns(a,v)
      real(dp),intent(in)::a(:,:)
      real(dp),intent(out)::v(:)
      integer::i,j,k
      k=0
      do j=1,size(a,2)
         do i=1,size(a,1)
            k=k+1; v(k)=a(i,j)
         end do
      end do
   end subroutine stack_columns

   subroutine predict_constrained_vglm(self,x,fitted,offsets)
      class(constrained_vglm_result_t),intent(in)::self
      real(dp),intent(in)::x(:,:)
      real(dp),allocatable,intent(out)::fitted(:,:)
      real(dp),intent(in),optional::offsets(:,:)
      real(dp),allocatable::eta(:,:)
      integer::i,j
      eta=matmul(x,self%coefficients)
      if(present(offsets))eta=eta+offsets
      allocate(fitted(size(eta,1),size(eta,2)))
      do j=1,size(eta,2)
         do i=1,size(eta,1)
            fitted(i,j)=link_inverse(eta(i,j),self%links(j))
            call clamp_mu(fitted(i,j),self%families(j))
         end do
      end do
   end subroutine predict_constrained_vglm

   elemental real(dp) function initial_mean(y,family) result(mu)
      real(dp),intent(in)::y
      integer,intent(in)::family
      select case(family)
      case(3)
         mu=min(0.999999_dp,max(0.000001_dp,(y+0.5_dp)/2.0_dp))
      case(1)
         mu=y
      case default
         mu=max(y,0.1_dp)
      end select
   end function initial_mean

   subroutine clamp_mu(mu,family)
      real(dp),intent(inout)::mu
      integer,intent(in)::family
      if(family==3)then
         mu=min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,mu))
      else if(family/=1)then
         mu=max(mu,1.0e-12_dp)
      end if
   end subroutine clamp_mu

   real(dp) function independent_deviance(y,mu,w,family) result(dev)
      real(dp),intent(in)::y(:,:),mu(:,:),w(:,:)
      integer,intent(in)::family(:)
      real(dp)::yi,mui
      integer::i,j
      dev=0.0_dp
      do j=1,size(y,2)
         do i=1,size(y,1)
            yi=y(i,j); mui=mu(i,j)
            select case(family(j))
            case(1)
               dev=dev+w(i,j)*(yi-mui)**2
            case(2)
               if(yi>0.0_dp)then
                  dev=dev+2.0_dp*w(i,j)*(yi*log(yi/mui)-(yi-mui))
               else
                  dev=dev+2.0_dp*w(i,j)*mui
               end if
            case(3)
               if(yi>0.0_dp)dev=dev+2.0_dp*w(i,j)*yi*log(yi/mui)
               if(yi<1.0_dp)dev=dev+2.0_dp*w(i,j)*(1.0_dp-yi)*log((1.0_dp-yi)/(1.0_dp-mui))
            case(4)
               if(yi>0.0_dp)dev=dev+2.0_dp*w(i,j)*((yi-mui)/mui-log(yi/mui))
            case(5)
               if(yi>0.0_dp)dev=dev+w(i,j)*(yi-mui)**2/(yi*mui*mui)
            end select
         end do
      end do
   end function independent_deviance

end module vgam_constraints
