module forecast_bats
   use forecast_kinds, only : dp, pi
   implicit none
   private
   public :: bats_make_w, bats_make_g, bats_make_f, tbats_make_w, tbats_make_g, tbats_make_f
   public :: state_space_filter, state_space_forecast, make_ci_matrix, make_si_matrix, make_ai_matrix
contains
   function bats_make_w(has_trend,phi,periods,ar,ma) result(w)
      logical,intent(in)::has_trend
      real(dp),intent(in)::phi,ar(:),ma(:)
      integer,intent(in)::periods(:)
      real(dp),allocatable::w(:)
      integer::n,tau,pos,i
      tau=sum(periods)
      n=1+merge(1,0,has_trend)+tau+size(ar)+size(ma)
      allocate(w(n))
      w=0.0_dp
      w(1)=1.0_dp
      pos=1
      if(has_trend)then
      pos=2
      w(2)=phi
      end if
      do i=1,size(periods)
      pos=pos+periods(i)
      w(pos)=1.0_dp
      end do
      if(size(ar)>0)w(2+merge(1,0,has_trend)+tau:1+merge(1,0,has_trend)+tau+size(ar))=ar
      pos=1+merge(1,0,has_trend)+tau+size(ar)
      if(size(ma)>0)w(pos+1:pos+size(ma))=ma
   end function
   function tbats_make_w(has_trend,phi,kvec,ar,ma) result(w)
      logical,intent(in)::has_trend
      real(dp),intent(in)::phi,ar(:),ma(:)
      integer,intent(in)::kvec(:)
      real(dp),allocatable::w(:)
      integer::n,tau,pos,i,j,adj
      tau=2*sum(kvec)
      adj=merge(1,0,has_trend)
      n=1+adj+tau+size(ar)+size(ma)
      allocate(w(n))
      w=0.0_dp
      w(1)=1.0_dp
      if(has_trend)w(2)=phi
      pos=1+adj
      do i=1,size(kvec)
      do j=1,kvec(i)
      w(pos+j)=1.0_dp
      end do
      pos=pos+2*kvec(i)
      end do
      if(size(ar)>0)w(2+adj+tau:1+adj+tau+size(ar))=ar
      pos=1+adj+tau+size(ar)
      if(size(ma)>0)w(pos+1:pos+size(ma))=ma
   end function
   function bats_make_g(alpha,has_trend,beta,gamma,periods,p,q) result(g)
      real(dp),intent(in)::alpha,beta,gamma(:)
      logical,intent(in)::has_trend
      integer,intent(in)::periods(:),p,q
      real(dp),allocatable::g(:)
      integer::adj,tau,n,pos,i
      adj=merge(1,0,has_trend)
      tau=sum(periods)
      n=1+adj+tau+p+q
      allocate(g(n))
      g=0.0_dp
      g(1)=alpha
      if(has_trend)g(2)=beta
      pos=1+adj
      do i=1,size(periods)
      if(i<=size(gamma))g(pos+1)=gamma(i)
      pos=pos+periods(i)
      end do
      pos=1+adj+tau
      if(p>0)g(pos+1)=1.0_dp
      if(q>0)g(pos+p+1)=1.0_dp
   end function
   function tbats_make_g(alpha,has_trend,beta,kvec,gamma1,gamma2,p,q) result(g)
      real(dp),intent(in)::alpha,beta,gamma1(:),gamma2(:)
      logical,intent(in)::has_trend
      integer,intent(in)::kvec(:),p,q
      real(dp),allocatable::g(:)
      integer::adj,tau,n,pos,i,k
      adj=merge(1,0,has_trend)
      tau=2*sum(kvec)
      n=1+adj+tau+p+q
      allocate(g(n))
      g=0.0_dp
      g(1)=alpha
      if(has_trend)g(2)=beta
      pos=1+adj
      do i=1,size(kvec)
         do k=1,kvec(i)
         if(i<=size(gamma1))g(pos+k)=gamma1(i)
         if(i<=size(gamma2))g(pos+kvec(i)+k)=gamma2(i)
         end do
         pos=pos+2*kvec(i)
      end do
      pos=1+adj+tau
      if(p>0)g(pos+1)=1.0_dp
      if(q>0)g(pos+p+1)=1.0_dp
   end function
   function seasonal_shift_matrix(periods) result(A)
      integer,intent(in)::periods(:)
      real(dp),allocatable::A(:,:)
      integer::tau,pos,m,j
      tau=sum(periods)
      allocate(A(tau,tau))
      A=0.0_dp
      pos=0
      do m=1,size(periods)
         A(pos+1,pos+periods(m))=1.0_dp
         do j=2,periods(m)
         A(pos+j,pos+j-1)=1.0_dp
         end do
         pos=pos+periods(m)
      end do
   end function
   function make_ci_matrix(k,m) result(C)
      integer,intent(in)::k
      real(dp),intent(in)::m
      real(dp),allocatable::C(:,:)
      integer::j
      allocate(C(k,k))
      C=0.0_dp
      do j=1,k
      C(j,j)=cos(2.0_dp*pi*real(j,dp)/m)
      end do
   end function
   function make_si_matrix(k,m) result(S)
      integer,intent(in)::k
      real(dp),intent(in)::m
      real(dp),allocatable::S(:,:)
      integer::j
      allocate(S(k,k))
      S=0.0_dp
      do j=1,k
      S(j,j)=sin(2.0_dp*pi*real(j,dp)/m)
      end do
   end function
   function make_ai_matrix(C,S) result(A)
      real(dp),intent(in)::C(:,:),S(:,:)
      real(dp),allocatable::A(:,:)
      integer::k
      k=size(C,1)
      allocate(A(2*k,2*k))
      A=0.0_dp
      A(1:k,1:k)=C
      A(1:k,k+1:)=S
      A(k+1:,1:k)=-S
      A(k+1:,k+1:)=C
   end function
   function trig_seasonal_matrix(periods,kvec) result(A)
      real(dp),intent(in)::periods(:)
      integer,intent(in)::kvec(:)
      real(dp),allocatable::A(:,:),C(:,:),S(:,:),Ai(:,:)
      integer::tau,pos,i,k
      tau=2*sum(kvec)
      allocate(A(tau,tau))
      A=0.0_dp
      pos=0
      do i=1,size(kvec)
      k=kvec(i)
      C=make_ci_matrix(k,periods(i))
      if(nint(periods(i))==2)C=0.0_dp
      S=make_si_matrix(k,periods(i))
      Ai=make_ai_matrix(C,S)
      A(pos+1:pos+2*k,pos+1:pos+2*k)=Ai
      pos=pos+2*k
      end do
   end function
   function bats_make_f(alpha,has_trend,beta,phi,periods,gamma,ar,ma) result(F)
      real(dp),intent(in)::alpha,beta,phi,gamma(:),ar(:),ma(:)
      logical,intent(in)::has_trend
      integer,intent(in)::periods(:)
      real(dp),allocatable::F(:,:),A(:,:),gb(:),outer(:,:)
      integer::adj,tau,p,q,n,pos,i,r
      adj=merge(1,0,has_trend)
      tau=sum(periods)
      p=size(ar)
      q=size(ma)
      n=1+adj+tau+p+q
      allocate(F(n,n))
      F=0.0_dp
      F(1,1)=1.0_dp
      if(has_trend)F(1,2)=phi
      if(p>0)F(1,2+adj+tau:1+adj+tau+p)=alpha*ar
      if(q>0)F(1,2+adj+tau+p:1+adj+tau+p+q)=alpha*ma
      if(has_trend)then
      F(2,2)=phi
      if(p>0)F(2,2+adj+tau:1+adj+tau+p)=beta*ar
      if(q>0)F(2,2+adj+tau+p:)=beta*ma
      end if
      if(tau>0)then
         A=seasonal_shift_matrix(periods)
         F(2+adj:1+adj+tau,2+adj:1+adj+tau)=A
         allocate(gb(tau))
         gb=0.0_dp
         pos=0
         do i=1,size(periods)
         if(i<=size(gamma))gb(pos+1)=gamma(i)
         pos=pos+periods(i)
         end do
         if(p>0)then
         do r=1,tau
         F(1+adj+r,2+adj+tau:1+adj+tau+p)=gb(r)*ar
         end do
         end if
         if(q>0)then
         do r=1,tau
         F(1+adj+r,2+adj+tau+p:)=gb(r)*ma
         end do
         end if
      end if
      call fill_arma_rows(F,1+adj+tau,ar,ma)
   end function
   function tbats_make_f(alpha,has_trend,beta,phi,periods,kvec,gamma1,gamma2,ar,ma) result(F)
      real(dp),intent(in)::alpha,beta,phi,periods(:),gamma1(:),gamma2(:),ar(:),ma(:)
      logical,intent(in)::has_trend
      integer,intent(in)::kvec(:)
      real(dp),allocatable::F(:,:),A(:,:),gb(:)
      integer::adj,tau,p,q,n,pos,i,k,r
      adj=merge(1,0,has_trend)
      tau=2*sum(kvec)
      p=size(ar)
      q=size(ma)
      n=1+adj+tau+p+q
      allocate(F(n,n))
      F=0.0_dp
      F(1,1)=1.0_dp
      if(has_trend)F(1,2)=phi
      if(p>0)F(1,2+adj+tau:1+adj+tau+p)=alpha*ar
      if(q>0)F(1,2+adj+tau+p:)=alpha*ma
      if(has_trend)then
      F(2,2)=phi
      if(p>0)F(2,2+adj+tau:1+adj+tau+p)=beta*ar
      if(q>0)F(2,2+adj+tau+p:)=beta*ma
      end if
      if(tau>0)then
         A=trig_seasonal_matrix(periods,kvec)
         F(2+adj:1+adj+tau,2+adj:1+adj+tau)=A
         allocate(gb(tau))
         gb=0.0_dp
         pos=0
         do i=1,size(kvec)
         k=kvec(i)
         if(i<=size(gamma1))gb(pos+1:pos+k)=gamma1(i)
         if(i<=size(gamma2))gb(pos+k+1:pos+2*k)=gamma2(i)
         pos=pos+2*k
         end do
         if(p>0)then
         do r=1,tau
         F(1+adj+r,2+adj+tau:1+adj+tau+p)=gb(r)*ar
         end do
         end if
         if(q>0)then
         do r=1,tau
         F(1+adj+r,2+adj+tau+p:)=gb(r)*ma
         end do
         end if
      end if
      call fill_arma_rows(F,1+adj+tau,ar,ma)
   end function
   subroutine fill_arma_rows(F,offset,ar,ma)
      real(dp),intent(inout)::F(:,:)
      integer,intent(in)::offset
      real(dp),intent(in)::ar(:),ma(:)
      integer::p,q,i
      p=size(ar)
      q=size(ma)
      if(p>0)then
         F(offset+1,offset+1:offset+p)=ar
         if(q>0)F(offset+1,offset+p+1:offset+p+q)=ma
         do i=2,p
         F(offset+i,offset+i-1)=1.0_dp
         end do
      end if
      if(q>1)then
      do i=2,q
      F(offset+p+i,offset+p+i-1)=1.0_dp
      end do
      end if
   end subroutine
   subroutine state_space_filter(y,w,F,g,x0,fitted,residuals,states,loglik)
      real(dp),intent(in)::y(:),w(:),F(:,:),g(:),x0(:)
      real(dp),allocatable,intent(out)::fitted(:),residuals(:),states(:,:)
      real(dp),intent(out),optional::loglik
      integer::t,n
      real(dp)::rss
      n=size(w)
      if(size(F,1)/=n.or.size(F,2)/=n.or.size(g)/=n.or.size(x0)/=n)error stop 'state_space_filter: dimensions'
      allocate(fitted(size(y)),residuals(size(y)),states(n,size(y)+1))
      states(:,1)=x0
      rss=0.0_dp
      do t=1,size(y)
      fitted(t)=dot_product(w,states(:,t))
      residuals(t)=y(t)-fitted(t)
      states(:,t+1)=matmul(F,states(:,t))+g*residuals(t)
      rss=rss+residuals(t)**2
      end do
      if(present(loglik))loglik=-0.5_dp*real(size(y),dp)*(log(2.0_dp*pi*max(rss/real(size(y),dp),1.0e-300_dp))+1.0_dp)
   end subroutine
   function state_space_forecast(w,F,state,h) result(fc)
      real(dp),intent(in)::w(:),F(:,:),state(:)
      integer,intent(in)::h
      real(dp),allocatable::fc(:),x(:)
      integer::i
      allocate(fc(h))
      x=state
      do i=1,h
      fc(i)=dot_product(w,x)
      x=matmul(F,x)
      end do
   end function
end module forecast_bats
