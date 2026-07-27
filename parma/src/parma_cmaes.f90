! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! The upstream R package translates Nikolaus Hansen's CMA-ES 3.6.
! This module provides an independent modern Fortran implementation.
module parma_cmaes
   use parma_kinds, only: dp
   use parma_types, only: objective_callback, cmaes_result, parma_options
   use parma_rng, only: seed_rng, random_normals
   use parma_linalg, only: jacobi_eigen, vector_norm
   implicit none
   private
   public :: cmaes_minimize

contains

   subroutine cmaes_minimize(objective, x0, result, options, lb, ub)
      procedure(objective_callback) :: objective
      real(dp), intent(in) :: x0(:)
      type(cmaes_result), intent(out) :: result
      type(parma_options), intent(in), optional :: options
      real(dp), intent(in), optional :: lb(:), ub(:)
      type(parma_options) :: opt
      real(dp), allocatable :: mean(:), old_mean(:), pc(:), ps(:)
      real(dp), allocatable :: cmat(:,:), bmat(:,:), invsqrt(:,:), ymat(:,:)
      real(dp), allocatable :: pop(:,:), z(:), y(:), evals(:), weights(:)
      real(dp), allocatable :: eig(:), dvec(:), bestx(:)
      integer, allocatable :: order(:)
      real(dp) :: sigma, mueff, cc, cs, c1, cmu, damps, chin
      real(dp) :: hsig, psnorm, scale, best, previous_best, rel_improve
      real(dp) :: cfac, floor_eig
      integer :: n, lambda, mu, i, k, iter, info, stagnation
      logical :: bounded

      opt = parma_options()
      if (present(options)) opt = options
      n = size(x0)
      if (n <= 0) then
         result%status = -1
         return
      end if
      lambda = opt%population
      if (lambda <= 0) lambda = 4 + int(3.0_dp*log(real(n,dp)))
      lambda = max(lambda,4)
      mu = lambda/2
      allocate(mean(n),old_mean(n),pc(n),ps(n),cmat(n,n),bmat(n,n), &
         invsqrt(n,n),ymat(n,lambda),pop(n,lambda),z(n),y(n),evals(lambda), &
         weights(mu),eig(n),dvec(n),bestx(n),order(lambda))
      mean = x0
      bounded = present(lb) .and. present(ub)
      if (bounded) then
         if (size(lb) /= n .or. size(ub) /= n) then
            result%status = -2
            return
         end if
         mean = min(max(mean,lb),ub)
      end if
      do i = 1, mu
         weights(i) = log(real(mu,dp)+0.5_dp)-log(real(i,dp))
      end do
      weights = weights/sum(weights)
      mueff = 1.0_dp/sum(weights*weights)
      cc = (4.0_dp+mueff/real(n,dp))/(real(n,dp)+4.0_dp+2.0_dp*mueff/real(n,dp))
      cs = (mueff+2.0_dp)/(real(n,dp)+mueff+5.0_dp)
      c1 = 2.0_dp/((real(n,dp)+1.3_dp)**2+mueff)
      cmu = min(1.0_dp-c1,2.0_dp*(mueff-2.0_dp+1.0_dp/mueff)/ &
         ((real(n,dp)+2.0_dp)**2+mueff))
      damps = 1.0_dp+2.0_dp*max(0.0_dp,sqrt((mueff-1.0_dp)/(real(n,dp)+1.0_dp))-1.0_dp)+cs
      chin = sqrt(real(n,dp))*(1.0_dp-1.0_dp/(4.0_dp*real(n,dp))+1.0_dp/(21.0_dp*real(n*n,dp)))
      sigma = max(opt%sigma0,1.0e-12_dp)
      if (bounded) sigma = sigma*max(sum(ub-lb)/real(n,dp),1.0_dp)
      cmat = 0.0_dp
      bmat = 0.0_dp
      invsqrt = 0.0_dp
      do i = 1, n
         cmat(i,i) = 1.0_dp
         bmat(i,i) = 1.0_dp
         invsqrt(i,i) = 1.0_dp
      end do
      dvec = 1.0_dp
      pc = 0.0_dp
      ps = 0.0_dp
      call seed_rng(opt%seed)
      best = objective(mean)
      result%evaluations = 1
      bestx = mean
      previous_best = best
      stagnation = 0
      floor_eig = 1.0e-20_dp

      do iter = 1, opt%max_iter
         call jacobi_eigen(cmat,eig,bmat,info,tol=1.0e-13_dp,max_sweeps=max(100,50*n*n))
         eig = max(eig,floor_eig)
         dvec = sqrt(eig)
         invsqrt = matmul(bmat,matmul(diagonal(1.0_dp/dvec),transpose(bmat)))
         do k = 1, lambda
            call random_normals(z)
            y = matmul(bmat,dvec*z)
            pop(:,k) = mean+sigma*y
            if (bounded) pop(:,k) = reflect_bounds(pop(:,k),lb,ub)
            ymat(:,k) = (pop(:,k)-mean)/sigma
            evals(k) = objective(pop(:,k))
         end do
         result%evaluations = result%evaluations+lambda
         call sort_indices(evals,order)
         if (evals(order(1)) < best) then
            best = evals(order(1))
            bestx = pop(:,order(1))
         end if
         old_mean = mean
         mean = 0.0_dp
         do i = 1, mu
            mean = mean+weights(i)*pop(:,order(i))
         end do
         y = (mean-old_mean)/sigma
         ps = (1.0_dp-cs)*ps+sqrt(cs*(2.0_dp-cs)*mueff)*matmul(invsqrt,y)
         psnorm = vector_norm(ps)
         scale = sqrt(max(1.0_dp-(1.0_dp-cs)**(2*iter),tiny(1.0_dp)))
         if (psnorm/scale/chin < 1.4_dp+2.0_dp/(real(n,dp)+1.0_dp)) then
            hsig = 1.0_dp
         else
            hsig = 0.0_dp
         end if
         pc = (1.0_dp-cc)*pc+hsig*sqrt(cc*(2.0_dp-cc)*mueff)*y
         cfac = 1.0_dp-c1-cmu+c1*(1.0_dp-hsig)*cc*(2.0_dp-cc)
         cmat = cfac*cmat+c1*outer(pc,pc)
         do i = 1, mu
            cmat = cmat+cmu*weights(i)*outer(ymat(:,order(i)),ymat(:,order(i)))
         end do
         cmat = 0.5_dp*(cmat+transpose(cmat))
         sigma = sigma*exp((cs/damps)*(psnorm/chin-1.0_dp))
         sigma = max(sigma,1.0e-15_dp)

         rel_improve = abs(previous_best-best)/max(1.0_dp,abs(best))
         if (rel_improve < opt%tol) then
            stagnation = stagnation+1
         else
            stagnation = 0
         end if
         previous_best = best
         if (opt%print_level > 0 .and. mod(iter,max(1,opt%print_level)) == 0) then
            write(*,'(a,i0,a,es14.6,a,es12.4)') 'cmaes iter=',iter,' value=',best,' sigma=',sigma
         end if
         if (sigma*maxval(dvec) < opt%tol*(1.0_dp+vector_norm(mean))) exit
         if (stagnation >= 50+10*n) exit
      end do
      allocate(result%x(n))
      result%x = bestx
      result%value = best
      result%iterations = min(iter,opt%max_iter)
      if (best < huge(1.0_dp)/10.0_dp) then
         result%status = 0
      else
         result%status = 1
      end if
   end subroutine cmaes_minimize

   function diagonal(x) result(a)
      real(dp), intent(in) :: x(:)
      real(dp) :: a(size(x),size(x))
      integer :: i
      a = 0.0_dp
      do i = 1, size(x)
         a(i,i) = x(i)
      end do
   end function diagonal

   function outer(x,y) result(a)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: a(size(x),size(y))
      a = spread(x,2,size(y))*spread(y,1,size(x))
   end function outer

   function reflect_bounds(x,lb,ub) result(y)
      real(dp), intent(in) :: x(:),lb(:),ub(:)
      real(dp) :: y(size(x)), width
      integer :: i, count
      do i = 1, size(x)
         width = ub(i)-lb(i)
         if (width <= 0.0_dp) then
            y(i) = lb(i)
         else
            y(i) = x(i)
            count = 0
            do while ((y(i) < lb(i) .or. y(i) > ub(i)) .and. count < 20)
               if (y(i) < lb(i)) y(i) = lb(i)+(lb(i)-y(i))
               if (y(i) > ub(i)) y(i) = ub(i)-(y(i)-ub(i))
               count = count+1
            end do
            y(i) = min(max(y(i),lb(i)),ub(i))
         end if
      end do
   end function reflect_bounds

   subroutine sort_indices(values,indices)
      real(dp), intent(in) :: values(:)
      integer, intent(out) :: indices(:)
      integer :: i,j,key
      do i = 1, size(values)
         indices(i) = i
      end do
      do i = 2, size(values)
         key = indices(i)
         j = i-1
         do while (j >= 1)
            if (values(indices(j)) <= values(key)) exit
            indices(j+1) = indices(j)
            j = j-1
         end do
         indices(j+1) = key
      end do
   end subroutine sort_indices

end module parma_cmaes
