! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_simulation
   use r_compat, only: dp, runif1
   use matrixdist_ph, only: ph_exit_rates
   use matrixdist_dph, only: dph_exit_probs
   use matrixdist_iph, only: iph_inverse_transform
   use matrixdist_transformations, only: merge_matrices
   implicit none
   private
   public :: embedded_mc, cumulate_matrix, cumulate_vector, initial_state, new_state
   public :: rphasetype, riph, rmatrixgev, rdphasetype, rmphstar, rmdphstar
   public :: rmph, rmdph, rmiph, rbivph, rbivdph, rbiviph

contains

   function embedded_mc(s) result(q)
      real(dp),intent(in)::s(:,:)
      real(dp),allocatable::q(:,:)
      real(dp),allocatable::exitv(:)
      integer::p,i,j
      p=size(s,1)
      allocate(q(p+1,p+1))
      q=0.0_dp
      exitv=ph_exit_rates(s)
      do i=1,p
         do j=1,p
            if(j/=i)q(i,j)=-s(i,j)/s(i,i)
         end do
         q(i,p+1)=-exitv(i)/s(i,i)
      end do
      q(p+1,p+1)=1.0_dp
   end function embedded_mc

   function dph_transition_with_absorption(s) result(q)
      real(dp),intent(in)::s(:,:)
      real(dp),allocatable::q(:,:)
      real(dp),allocatable::exitv(:)
      integer::p
      p=size(s,1)
      allocate(q(p+1,p+1))
      q=0.0_dp
      q(1:p,1:p)=s
      exitv=dph_exit_probs(s)
      q(1:p,p+1)=exitv
      q(p+1,p+1)=1.0_dp
   end function dph_transition_with_absorption

   function cumulate_matrix(a) result(c)
      real(dp),intent(in)::a(:,:)
      real(dp)::c(size(a,1),size(a,2))
      integer::i,j
      do i=1,size(a,1)
         c(i,1)=a(i,1)
         do j=2,size(a,2)
         c(i,j)=c(i,j-1)+a(i,j)
         end do
      end do
   end function cumulate_matrix

   function cumulate_vector(a) result(c)
      real(dp),intent(in)::a(:)
      real(dp)::c(size(a))
      integer::i
      if(size(a)==0)return
      c(1)=a(1)
      do i=2,size(a)
      c(i)=c(i-1)+a(i)
      end do
   end function cumulate_vector

   integer function initial_state(cum_alpha,u) result(state)
      real(dp),intent(in)::cum_alpha(:),u
      integer::i
      state=1
      do i=1,size(cum_alpha)
         if(u<=cum_alpha(i))then
         state=i
         return
         end if
      end do
      state=size(cum_alpha)
   end function initial_state

   integer function new_state(prev,cum_q,u) result(state)
      integer,intent(in)::prev
      real(dp),intent(in)::cum_q(:,:),u
      integer::j
      state=size(cum_q,2)
      do j=1,size(cum_q,2)
         if(u<=cum_q(prev,j))then
         state=j
         return
         end if
      end do
   end function new_state

   function rphasetype(n,alpha,s) result(sample)
      integer,intent(in)::n
      real(dp),intent(in)::alpha(:),s(:,:)
      real(dp),allocatable::sample(:)
      real(dp),allocatable::cq(:,:),ca(:)
      real(dp)::time,u
      integer::i,state,p
      p=size(alpha)
      allocate(sample(n))
      cq=cumulate_matrix(embedded_mc(s))
      ca=cumulate_vector(alpha)
      do i=1,n
         time=0.0_dp
         state=initial_state(ca,runif1())
         do while(state/=p+1)
            u=max(tiny(1.0_dp),1.0_dp-runif1())
            time=time+log(u)/s(state,state)
            state=new_state(state,cq,runif1())
         end do
         sample(i)=time
      end do
   end function rphasetype

   function riph(n,kind,alpha,s,beta) result(sample)
      integer,intent(in)::n
      character(len=*),intent(in)::kind
      real(dp),intent(in)::alpha(:),s(:,:),beta(:)
      real(dp),allocatable::sample(:),base(:)
      integer::i
      base=rphasetype(n,alpha,s)
      allocate(sample(n))
      do i=1,n
      sample(i)=iph_inverse_transform(base(i),kind,beta)
      end do
   end function riph

   function rmatrixgev(n,alpha,s,mu,sigma,xi) result(sample)
      integer,intent(in)::n
      real(dp),intent(in)::alpha(:),s(:,:),mu,sigma
      real(dp),intent(in),optional::xi
      real(dp),allocatable::sample(:)
      real(dp)::shape
      shape=0.0_dp
      if(present(xi))shape=xi
      sample=riph(n,'gev',alpha,s,[mu,sigma,shape])
   end function rmatrixgev

   function rdphasetype(n,alpha,s) result(sample)
      integer,intent(in)::n
      real(dp),intent(in)::alpha(:),s(:,:)
      integer,allocatable::sample(:)
      real(dp),allocatable::cq(:,:),ca(:)
      integer::i,state,p,time
      p=size(alpha)
      allocate(sample(n))
      cq=cumulate_matrix(dph_transition_with_absorption(s))
      ca=cumulate_vector(alpha)
      do i=1,n
         time=0
         state=initial_state(ca,runif1())
         do while(state/=p+1)
            time=time+1
            state=new_state(state,cq,runif1())
         end do
         sample(i)=time
      end do
   end function rdphasetype

   function rmphstar(n,alpha,s,reward) result(sample)
      integer,intent(in)::n
      real(dp),intent(in)::alpha(:),s(:,:),reward(:,:)
      real(dp),allocatable::sample(:,:)
      real(dp),allocatable::cq(:,:),ca(:)
      real(dp)::time,u
      integer::i,j,state,p,d
      p=size(alpha)
      d=size(reward,2)
      allocate(sample(n,d))
      sample=0.0_dp
      cq=cumulate_matrix(embedded_mc(s))
      ca=cumulate_vector(alpha)
      do i=1,n
         state=initial_state(ca,runif1())
         do while(state/=p+1)
            u=max(tiny(1.0_dp),1.0_dp-runif1())
            time=log(u)/s(state,state)
            do j=1,d
            sample(i,j)=sample(i,j)+reward(state,j)*time
            end do
            state=new_state(state,cq,runif1())
         end do
      end do
   end function rmphstar

   function rmdphstar(n,alpha,s,reward) result(sample)
      integer,intent(in)::n
      real(dp),intent(in)::alpha(:),s(:,:),reward(:,:)
      real(dp),allocatable::sample(:,:)
      real(dp),allocatable::cq(:,:),ca(:)
      integer::i,j,state,p,d
      p=size(alpha)
      d=size(reward,2)
      allocate(sample(n,d))
      sample=0.0_dp
      cq=cumulate_matrix(dph_transition_with_absorption(s))
      ca=cumulate_vector(alpha)
      do i=1,n
         state=initial_state(ca,runif1())
         do while(state/=p+1)
            do j=1,d
            sample(i,j)=sample(i,j)+reward(state,j)
            end do
            state=new_state(state,cq,runif1())
         end do
      end do
   end function rmdphstar

   function rmph(n,alpha,s) result(sample)
      integer,intent(in)::n
      real(dp),intent(in)::alpha(:),s(:,:,:)
      real(dp),allocatable::sample(:,:),one(:),ca(:)
      real(dp),allocatable::unit(:)
      integer::i,j,h,p,d
      p=size(alpha)
      d=size(s,3)
      allocate(sample(n,d),unit(p))
      ca=cumulate_vector(alpha)
      do i=1,n
         h=initial_state(ca,runif1())
         unit=0.0_dp
         unit(h)=1.0_dp
         do j=1,d
         one=rphasetype(1,unit,s(:,:,j))
         sample(i,j)=one(1)
         end do
      end do
   end function rmph

   function rmdph(n,alpha,s) result(sample)
      integer,intent(in)::n
      real(dp),intent(in)::alpha(:),s(:,:,:)
      integer,allocatable::sample(:,:),one(:)
      real(dp),allocatable::ca(:),unit(:)
      integer::i,j,h,p,d
      p=size(alpha)
      d=size(s,3)
      allocate(sample(n,d),unit(p))
      ca=cumulate_vector(alpha)
      do i=1,n
         h=initial_state(ca,runif1())
         unit=0.0_dp
         unit(h)=1.0_dp
         do j=1,d
         one=rdphasetype(1,unit,s(:,:,j))
         sample(i,j)=one(1)
         end do
      end do
   end function rmdph

   function rmiph(n,alpha,s,kinds,beta) result(sample)
      integer,intent(in)::n
      real(dp),intent(in)::alpha(:),s(:,:,:),beta(:,:)
      character(len=*),intent(in)::kinds(:)
      real(dp),allocatable::sample(:,:),base(:,:)
      integer::i,j,d
      d=size(s,3)
      if(size(kinds)/=d .or. size(beta,2)/=d) error stop 'rmiph: dimension mismatch'
      base=rmph(n,alpha,s)
      allocate(sample(n,d))
      do j=1,d
      do i=1,n
      sample(i,j)=iph_inverse_transform(base(i,j),kinds(j),beta(:,j))
      end do
      end do
   end function rmiph

   function rbivph(n,alpha,s11,s12,s22) result(sample)
      integer,intent(in)::n
      real(dp),intent(in)::alpha(:),s11(:,:),s12(:,:),s22(:,:)
      real(dp),allocatable::sample(:,:),st(:,:),reward(:,:),a(:)
      integer::p1,p2
      p1=size(s11,1)
      p2=size(s22,1)
      st=merge_matrices(s11,s12,s22)
      allocate(a(p1+p2),reward(p1+p2,2))
      a=0.0_dp
      a(1:p1)=alpha
      reward=0.0_dp
      reward(1:p1,1)=1.0_dp
      reward(p1+1:,2)=1.0_dp
      sample=rmphstar(n,a,st,reward)
   end function rbivph

   function rbivdph(n,alpha,s11,s12,s22) result(sample)
      integer,intent(in)::n
      real(dp),intent(in)::alpha(:),s11(:,:),s12(:,:),s22(:,:)
      real(dp),allocatable::st(:,:),reward(:,:),a(:),tmp(:,:)
      integer,allocatable::sample(:,:)
      integer::p1,p2
      p1=size(s11,1)
      p2=size(s22,1)
      st=merge_matrices(s11,s12,s22)
      allocate(a(p1+p2),reward(p1+p2,2))
      a=0.0_dp
      a(1:p1)=alpha
      reward=0.0_dp
      reward(1:p1,1)=1.0_dp
      reward(p1+1:,2)=1.0_dp
      tmp=rmdphstar(n,a,st,reward)
      allocate(sample(n,2))
      sample=nint(tmp)
   end function rbivdph

   function rbiviph(n,alpha,s11,s12,s22,kinds,beta) result(sample)
      integer,intent(in)::n
      real(dp),intent(in)::alpha(:),s11(:,:),s12(:,:),s22(:,:),beta(:,:)
      character(len=*),intent(in)::kinds(2)
      real(dp),allocatable::sample(:,:),base(:,:)
      integer::i,j
      base=rbivph(n,alpha,s11,s12,s22)
      allocate(sample(n,2))
      do j=1,2
      do i=1,n
      sample(i,j)=iph_inverse_transform(base(i,j),kinds(j),beta(:,j))
      end do
      end do
   end function rbiviph

end module matrixdist_simulation
