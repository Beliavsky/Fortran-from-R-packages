! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_categorical
   use vgam_kinds, only : dp
   use vgam_linalg, only : solve_linear, invert_matrix
   implicit none
   private
   type, public :: multinomial_result_t
      real(dp), allocatable :: coefficients(:,:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: fitted_probabilities(:,:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 0
      logical :: converged = .false.
      integer :: reference_category = 0
   contains
      procedure :: predict_proba => predict_multinomial
   end type multinomial_result_t
   public :: fit_multinomial, multinomial_probabilities, deviance_multinomial
contains

   subroutine multinomial_probabilities(x,beta,prob)
      real(dp),intent(in)::x(:,:),beta(:,:)
      real(dp),allocatable,intent(out)::prob(:,:)
      real(dp),allocatable::eta(:)
      real(dp)::mx,den
      integer::n,k,i,j
      n=size(x,1)
      k=size(beta,2)+1
      allocate(prob(n,k),eta(k))
      prob=0.0_dp
      do i=1,n
         eta(1:k-1)=matmul(x(i,:),beta)
         eta(k)=0.0_dp
         mx=maxval(eta)
         den=sum(exp(eta-mx))
         do j=1,k
      prob(i,j)=exp(eta(j)-mx)/den
      end do
      end do
   end subroutine multinomial_probabilities

   subroutine fit_multinomial(y,x,ncat,result,weights,max_iter,tol)
      integer,intent(in)::y(:),ncat
      real(dp),intent(in)::x(:,:)
      type(multinomial_result_t),intent(out)::result
      real(dp),intent(in),optional::weights(:),tol
      integer,intent(in),optional::max_iter
      real(dp),allocatable::w(:),beta(:,:),trial(:,:),prob(:,:)
      real(dp),allocatable::score(:),info(:,:),delta(:),cov(:,:)
      real(dp),allocatable::bvec(:),btry(:),beta_try(:,:),ptry(:,:)
      real(dp)::ll,lltry,tolerance,step,chg,pjk
      integer::n,p,k,m,i,j,l,a,b,iter,niter,stat,ia,ib
      n=size(y)
      p=size(x,2)
      k=ncat
      m=p*(k-1)
      if(size(x,1)/=n.or.k<2.or.any(y<1).or.any(y>k))then
         result%status=1
      return
      end if
      if(present(weights))then
         if(size(weights)/=n.or.any(weights<0.0_dp))then
            result%status=2
      return
         end if
         w=weights
      else
         allocate(w(n))
      w=1.0_dp
      end if
      niter=100
      if(present(max_iter))niter=max_iter
      tolerance=1.0e-8_dp
      if(present(tol))tolerance=tol
      allocate(beta(p,k-1),trial(n,k),score(m),info(m,m))
      allocate(bvec(m),btry(m),beta_try(p,k-1))
      beta=0.0_dp
      trial=0.0_dp
      do i=1,n
      trial(i,y(i))=1.0_dp
      end do
      call multinomial_probabilities(x,beta,prob)
      ll=multinom_loglik(trial,prob,w)
      do iter=1,niter
         score=0.0_dp
      info=0.0_dp
         do i=1,n
            do j=1,k-1
               ia=(j-1)*p
               score(ia+1:ia+p)=score(ia+1:ia+p)+ &
                  w(i)*x(i,:)*(trial(i,j)-prob(i,j))
               do l=1,k-1
                  ib=(l-1)*p
                  pjk=prob(i,j)*(merge(1.0_dp,0.0_dp,j==l)-prob(i,l))
                  do a=1,p
                     do b=1,p
                        info(ia+a,ib+b)=info(ia+a,ib+b)+ &
                           w(i)*pjk*x(i,a)*x(i,b)
                     end do
                  end do
               end do
            end do
         end do
         call solve_linear(info,score,delta,stat)
         if(stat/=0)then
      result%status=10+stat
      return
      end if
         bvec=reshape(beta,[m])
         step=1.0_dp
         do
            btry=bvec+step*delta
            beta_try=reshape(btry,[p,k-1])
            call multinomial_probabilities(x,beta_try,ptry)
            lltry=multinom_loglik(trial,ptry,w)
            if(lltry>=ll.or.step<1.0e-7_dp)exit
            step=0.5_dp*step
         end do
         chg=maxval(abs(step*delta))/max(1.0_dp,maxval(abs(bvec)))
         beta=beta_try
      prob=ptry
      ll=lltry
         if(chg<tolerance)then
      result%converged=.true.
      exit
      end if
      end do
      ! Recompute expected information at the solution.
      info=0.0_dp
      do i=1,n
         do j=1,k-1
            ia=(j-1)*p
            do l=1,k-1
               ib=(l-1)*p
               pjk=prob(i,j)*(merge(1.0_dp,0.0_dp,j==l)-prob(i,l))
               do a=1,p
                  do b=1,p
                     info(ia+a,ib+b)=info(ia+a,ib+b)+ &
                        w(i)*pjk*x(i,a)*x(i,b)
                  end do
               end do
            end do
         end do
      end do
      call invert_matrix(info,cov,stat)
      if(stat/=0)then
         allocate(cov(m,m))
      cov=0.0_dp
      result%status=20+stat
      end if
      result%coefficients=beta
      result%covariance=cov
      result%fitted_probabilities=prob
      result%loglik=ll
      result%aic=-2.0_dp*ll+2.0_dp*real(m,dp)
      result%iterations=iter
      result%reference_category=k
      if(.not.result%converged.and.result%status==0)result%status=100
   end subroutine fit_multinomial

   real(dp) function multinom_loglik(y,prob,w) result(ll)
      real(dp),intent(in)::y(:,:),prob(:,:),w(:)
      integer::i,j
      ll=0.0_dp
      do i=1,size(y,1)
         do j=1,size(y,2)
            if(y(i,j)>0.0_dp)ll=ll+w(i)*y(i,j)*log(max(prob(i,j),tiny(1.0_dp)))
         end do
      end do
   end function multinom_loglik

   real(dp) function deviance_multinomial(y,prob,weights) result(dev)
      integer,intent(in)::y(:)
      real(dp),intent(in)::prob(:,:)
      real(dp),intent(in),optional::weights(:)
      integer::i
      dev=0.0_dp
      if(present(weights))then
         do i=1,size(y)
      dev=dev-2.0_dp*weights(i)*log(max(prob(i,y(i)),tiny(1.0_dp)))
      end do
      else
         do i=1,size(y)
      dev=dev-2.0_dp*log(max(prob(i,y(i)),tiny(1.0_dp)))
      end do
      end if
   end function deviance_multinomial

   function predict_multinomial(self,x) result(prob)
      class(multinomial_result_t),intent(in)::self
      real(dp),intent(in)::x(:,:)
      real(dp),allocatable::prob(:,:)
      call multinomial_probabilities(x,self%coefficients,prob)
   end function predict_multinomial

end module vgam_categorical
