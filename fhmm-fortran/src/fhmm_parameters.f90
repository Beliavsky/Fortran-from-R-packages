! SPDX-License-Identifier: GPL-3.0-only
! Numerical translation of fHMM 1.4.3.
module fhmm_parameters
   use fhmm_kinds, only: dp
   use fhmm_types, only: hmm_parameters, hhmm_parameters
   use fhmm_types, only: dist_gamma, dist_poisson, dist_student_t
   implicit none
   private
   public :: transition_to_unconstrained, unconstrained_to_transition
   public :: pack_hmm_parameters, unpack_hmm_parameters, hmm_parameter_count
   public :: pack_hhmm_parameters, unpack_hhmm_parameters, hhmm_parameter_count
   public :: validate_hmm_parameters, validate_hhmm_parameters

contains

   function transition_to_unconstrained(gamma) result(u)
      real(dp), intent(in) :: gamma(:, :)
      real(dp), allocatable :: u(:)
      integer :: n,i,j,k
      real(dp) :: diagp
      n=size(gamma,1); allocate(u(n*(n-1))); k=0
      do j=1,n
         do i=1,n
            if(i==j) cycle
            k=k+1
            diagp=max(gamma(i,i),1.0e-12_dp)
            u(k)=log(max(gamma(i,j),1.0e-12_dp)/diagp)
         end do
      end do
   end function transition_to_unconstrained

   function unconstrained_to_transition(u,n,numerical_safeguard) result(gamma)
      real(dp), intent(in) :: u(:)
      integer, intent(in) :: n
      logical, intent(in), optional :: numerical_safeguard
      real(dp), allocatable :: gamma(:, :)
      real(dp) :: row_sum,val
      integer :: i,j,k
      logical :: safe
      safe=.false.; if(present(numerical_safeguard)) safe=numerical_safeguard
      allocate(gamma(n,n)); gamma=0.0_dp; k=0
      do j=1,n
         do i=1,n
            if(i==j) cycle
            k=k+1
            val=exp(min(u(k),700.0_dp))
            if(safe) val=min(val,100.0_dp)
            gamma(i,j)=val
         end do
      end do
      do i=1,n
         gamma(i,i)=1.0_dp
         row_sum=sum(gamma(i,:))
         gamma(i,:)=gamma(i,:)/row_sum
      end do
   end function unconstrained_to_transition

   integer function hmm_parameter_count(n,family) result(k)
      integer,intent(in)::n,family
      k=n*(n-1)+n
      if(family/=dist_poisson) k=k+n
      if(family==dist_student_t) k=k+n
   end function hmm_parameter_count

   function pack_hmm_parameters(par) result(u)
      type(hmm_parameters), intent(in) :: par
      real(dp), allocatable :: u(:),tmp(:)
      integer :: n,k,m
      n=size(par%mu); m=hmm_parameter_count(n,par%distribution); allocate(u(m)); k=0
      tmp=transition_to_unconstrained(par%gamma); u(1:size(tmp))=tmp; k=size(tmp)
      if(par%distribution==dist_gamma .or. par%distribution==dist_poisson) then
         u(k+1:k+n)=log(max(par%mu,1.0e-12_dp))
      else
         u(k+1:k+n)=par%mu
      end if
      k=k+n
      if(par%distribution/=dist_poisson) then
         u(k+1:k+n)=log(max(par%sigma,1.0e-12_dp)); k=k+n
      end if
      if(par%distribution==dist_student_t) then
         u(k+1:k+n)=log(max(par%df,1.0e-12_dp))
      end if
   end function pack_hmm_parameters

   function unpack_hmm_parameters(u,n,family,numerical_safeguard) result(par)
      real(dp), intent(in) :: u(:)
      integer, intent(in) :: n,family
      logical,intent(in),optional::numerical_safeguard
      type(hmm_parameters) :: par
      integer :: k,ng
      logical::safe
      safe=.false.;if(present(numerical_safeguard))safe=numerical_safeguard
      par%distribution=family; ng=n*(n-1); k=ng
      allocate(par%gamma(n,n),par%mu(n),par%sigma(n),par%df(n))
      par%gamma=unconstrained_to_transition(u(1:ng),n,safe)
      if(family==dist_gamma .or. family==dist_poisson) then
         par%mu=exp(min(u(k+1:k+n),700.0_dp))
      else
         par%mu=u(k+1:k+n)
      end if
      k=k+n
      if(family/=dist_poisson) then
         par%sigma=exp(min(u(k+1:k+n),700.0_dp)); k=k+n
         if(safe) par%sigma=min(par%sigma,100.0_dp)
      else
         par%sigma=1.0_dp
      end if
      if(family==dist_student_t) then
         par%df=exp(min(u(k+1:k+n),700.0_dp))
         if(safe) par%df=min(par%df,100.0_dp)
      else
         par%df=10.0_dp
      end if
   end function unpack_hmm_parameters

   integer function hhmm_parameter_count(m,n,coarse_family,fine_family) result(k)
      integer,intent(in)::m,n,coarse_family,fine_family
      k=hmm_parameter_count(m,coarse_family)+m*hmm_parameter_count(n,fine_family)
   end function hhmm_parameter_count

   function pack_hhmm_parameters(par) result(u)
      type(hhmm_parameters),intent(in)::par
      real(dp),allocatable::u(:),tmp(:)
      integer::m,k,s,total
      m=size(par%fine); total=size(pack_hmm_parameters(par%coarse))
      do s=1,m; total=total+size(pack_hmm_parameters(par%fine(s))); end do
      allocate(u(total)); tmp=pack_hmm_parameters(par%coarse); u(1:size(tmp))=tmp; k=size(tmp)
      do s=1,m
         tmp=pack_hmm_parameters(par%fine(s)); u(k+1:k+size(tmp))=tmp; k=k+size(tmp)
      end do
   end function pack_hhmm_parameters

   function unpack_hhmm_parameters(u,m,n,coarse_family,fine_family,numerical_safeguard) result(par)
      real(dp),intent(in)::u(:)
      integer,intent(in)::m,n,coarse_family,fine_family
      logical,intent(in),optional::numerical_safeguard
      type(hhmm_parameters)::par
      integer::kc,kf,k,s
      logical::safe
      safe=.false.;if(present(numerical_safeguard))safe=numerical_safeguard
      kc=hmm_parameter_count(m,coarse_family); kf=hmm_parameter_count(n,fine_family)
      par%coarse=unpack_hmm_parameters(u(1:kc),m,coarse_family,safe)
      allocate(par%fine(m)); k=kc
      do s=1,m
         par%fine(s)=unpack_hmm_parameters(u(k+1:k+kf),n,fine_family,safe); k=k+kf
      end do
   end function unpack_hhmm_parameters

   logical function validate_hmm_parameters(par) result(ok)
      type(hmm_parameters),intent(in)::par
      integer::n
      ok=.false.;if(.not.allocated(par%gamma).or..not.allocated(par%mu))return
      n=size(par%mu);if(any(shape(par%gamma)/=[n,n]))return
      if(any(par%gamma<0.0_dp).or.maxval(abs(sum(par%gamma,dim=2)-1.0_dp))>1.0e-8_dp)return
      if(par%distribution==dist_gamma.or.par%distribution==dist_poisson)then
         if(any(par%mu<=0.0_dp))return
      end if
      if(par%distribution/=dist_poisson)then
         if(.not.allocated(par%sigma).or.size(par%sigma)/=n.or.any(par%sigma<=0.0_dp))return
      end if
      if(par%distribution==dist_student_t)then
         if(.not.allocated(par%df).or.size(par%df)/=n.or.any(par%df<=0.0_dp))return
      end if
      ok=.true.
   end function validate_hmm_parameters

   logical function validate_hhmm_parameters(par) result(ok)
      type(hhmm_parameters),intent(in)::par
      integer::s
      ok=validate_hmm_parameters(par%coarse)
      if(.not.ok.or..not.allocated(par%fine))return
      do s=1,size(par%fine)
         if(.not.validate_hmm_parameters(par%fine(s)))then;ok=.false.;return;end if
      end do
   end function validate_hhmm_parameters

end module fhmm_parameters
