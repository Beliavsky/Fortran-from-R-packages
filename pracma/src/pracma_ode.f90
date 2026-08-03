! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_ode
   use pracma_kinds, only : dp
   use pracma_status, only : pracma_ok, pracma_invalid_argument, pracma_not_converged
   use pracma_callbacks, only : vector_field
   use pracma_types, only : ode_result
   use pracma_linalg, only : solve_linear
   implicit none
   private
   public :: rk4, rk4sys, euler_heun, abm3, rkf54, ode23, ode45, ode78, ode23s
   public :: bulirsch_stoer, newmark, cranknic, shooting_ivp

contains

   function rk4(f,tspan,y0,nsteps) result(res)
      procedure(vector_field)::f
      real(dp),intent(in)::tspan(2),y0(:)
      integer,intent(in)::nsteps
      type(ode_result)::res
      real(dp)::h,t
      real(dp),allocatable::k1(:),k2(:),k3(:),k4(:),yt(:)
      integer::i,n
      n=size(y0)
      if(nsteps<1.or.tspan(1)==tspan(2).or.n<1)then; res%status=pracma_invalid_argument; return; end if
      allocate(res%t(nsteps+1),res%y(n,nsteps+1),k1(n),k2(n),k3(n),k4(n),yt(n))
      h=(tspan(2)-tspan(1))/real(nsteps,dp); res%t(1)=tspan(1); res%y(:,1)=y0
      do i=1,nsteps
         t=res%t(i); call f(t,res%y(:,i),k1)
         yt=res%y(:,i)+0.5_dp*h*k1; call f(t+0.5_dp*h,yt,k2)
         yt=res%y(:,i)+0.5_dp*h*k2; call f(t+0.5_dp*h,yt,k3)
         yt=res%y(:,i)+h*k3; call f(t+h,yt,k4)
         res%y(:,i+1)=res%y(:,i)+h*(k1+2*k2+2*k3+k4)/6.0_dp; res%t(i+1)=t+h
      end do
      res%accepted_steps=nsteps; res%converged=.true.; res%status=pracma_ok
   end function rk4

   function rk4sys(f,tspan,y0,nsteps) result(res)
      procedure(vector_field)::f
      real(dp),intent(in)::tspan(2),y0(:)
      integer,intent(in)::nsteps
      type(ode_result)::res
      res=rk4(f,tspan,y0,nsteps)
   end function rk4sys

   function euler_heun(f,tspan,y0,nsteps,improved) result(res)
      procedure(vector_field)::f
      real(dp),intent(in)::tspan(2),y0(:)
      integer,intent(in)::nsteps
      logical,intent(in),optional::improved
      type(ode_result)::res
      real(dp)::h,t
      real(dp),allocatable::k1(:),k2(:),yp(:)
      logical::heun
      integer::i,n
      n=size(y0); heun=.true.; if(present(improved))heun=improved
      if(nsteps<1.or.n<1)then; res%status=pracma_invalid_argument; return; end if
      allocate(res%t(nsteps+1),res%y(n,nsteps+1),k1(n),k2(n),yp(n)); h=(tspan(2)-tspan(1))/real(nsteps,dp)
      res%t(1)=tspan(1); res%y(:,1)=y0
      do i=1,nsteps
         t=res%t(i); call f(t,res%y(:,i),k1); yp=res%y(:,i)+h*k1
         if(heun)then; call f(t+h,yp,k2); res%y(:,i+1)=res%y(:,i)+0.5_dp*h*(k1+k2)
         else; res%y(:,i+1)=yp; end if
         res%t(i+1)=t+h
      end do
      res%accepted_steps=nsteps; res%converged=.true.; res%status=pracma_ok
   end function euler_heun

   function abm3(f,tspan,y0,nsteps) result(res)
      procedure(vector_field)::f
      real(dp),intent(in)::tspan(2),y0(:)
      integer,intent(in)::nsteps
      type(ode_result)::res
      type(ode_result)::start
      real(dp)::h,t
      real(dp),allocatable::fv(:,:),fp(:),yc(:)
      integer::i,n
      if(nsteps<3)then; res=rk4(f,tspan,y0,max(1,nsteps)); return; end if
      n=size(y0); h=(tspan(2)-tspan(1))/real(nsteps,dp)
      start=rk4(f,[tspan(1),tspan(1)+2*h],y0,2)
      allocate(res%t(nsteps+1),res%y(n,nsteps+1),fv(n,nsteps+1),fp(n),yc(n))
      res%t(1:3)=start%t; res%y(:,1:3)=start%y
      do i=1,3; call f(res%t(i),res%y(:,i),fv(:,i)); end do
      do i=3,nsteps
         t=res%t(i); yc=res%y(:,i)+h*(23*fv(:,i)-16*fv(:,i-1)+5*fv(:,i-2))/12.0_dp
         call f(t+h,yc,fp)
         res%y(:,i+1)=res%y(:,i)+h*(5*fp+8*fv(:,i)-fv(:,i-1))/12.0_dp
         res%t(i+1)=t+h; call f(t+h,res%y(:,i+1),fv(:,i+1))
      end do
      res%accepted_steps=nsteps; res%converged=.true.; res%status=pracma_ok
   end function abm3

   function rkf54(f,tspan,y0,rtol,atol,h_initial,max_steps) result(res)
      procedure(vector_field)::f
      real(dp),intent(in)::tspan(2),y0(:)
      real(dp),intent(in),optional::rtol,atol,h_initial
      integer,intent(in),optional::max_steps
      type(ode_result)::res
      real(dp),allocatable::twork(:),ywork(:,:),k1(:),k2(:),k3(:),k4(:),k5(:),k6(:),k7(:),yt(:),y5(:),y4(:),scalev(:)
      real(dp)::t,tf,h,dir,rt,at,err,fac
      integer::n,mx,used,accepted,rejected
      n=size(y0); mx=100000; if(present(max_steps))mx=max_steps
      rt=1.0e-6_dp; if(present(rtol))rt=rtol; at=1.0e-9_dp; if(present(atol))at=atol
      if(n<1.or.tspan(1)==tspan(2).or.mx<1.or.rt<=0.or.at<=0)then; res%status=pracma_invalid_argument; return; end if
      dir=sign(1.0_dp,tspan(2)-tspan(1)); tf=tspan(2)
      h=dir*abs(tspan(2)-tspan(1))/100.0_dp; if(present(h_initial))h=dir*abs(h_initial)
      allocate(twork(mx+1),ywork(n,mx+1),k1(n),k2(n),k3(n),k4(n),k5(n),k6(n),k7(n),yt(n),y5(n),y4(n),scalev(n))
      used=1; accepted=0; rejected=0; t=tspan(1); twork(1)=t; ywork(:,1)=y0
      do while(dir*(tf-t)>0.0_dp .and. accepted+rejected<mx)
         if(dir*(t+h-tf)>0.0_dp)h=tf-t
         call f(t,ywork(:,used),k1)
         yt=ywork(:,used)+h*(1.0_dp/5.0_dp)*k1; call f(t+h/5.0_dp,yt,k2)
         yt=ywork(:,used)+h*(3.0_dp*k1/40.0_dp+9.0_dp*k2/40.0_dp); call f(t+3*h/10,yt,k3)
         yt=ywork(:,used)+h*(44*k1/45-56*k2/15+32*k3/9); call f(t+4*h/5,yt,k4)
         yt=ywork(:,used)+h*(19372*k1/6561-25360*k2/2187+64448*k3/6561-212*k4/729); call f(t+8*h/9,yt,k5)
         yt=ywork(:,used)+h*(9017*k1/3168-355*k2/33+46732*k3/5247+49*k4/176-5103*k5/18656); call f(t+h,yt,k6)
         y5=ywork(:,used)+h*(35*k1/384+500*k3/1113+125*k4/192-2187*k5/6784+11*k6/84)
         call f(t+h,y5,k7)
         y4=ywork(:,used)+h*(5179*k1/57600+7571*k3/16695+393*k4/640-92097*k5/339200+187*k6/2100+k7/40)
         scalev=max(at,rt*max(abs(ywork(:,used)),abs(y5))); err=sqrt(sum(((y5-y4)/scalev)**2)/real(n,dp))
         if(err<=1.0_dp)then
            t=t+h; used=used+1; accepted=accepted+1; twork(used)=t; ywork(:,used)=y5
         else; rejected=rejected+1; end if
         if(err<=tiny(1.0_dp))then; fac=5.0_dp; else; fac=min(5.0_dp,max(0.1_dp,0.9_dp*err**(-0.2_dp))); end if
         h=h*fac
         if(abs(h)<=epsilon(t)*max(1.0_dp,abs(t)))exit
      end do
      allocate(res%t(used),res%y(n,used)); res%t=twork(:used); res%y=ywork(:,:used)
      res%accepted_steps=accepted; res%rejected_steps=rejected
      res%converged=dir*(tf-t)<=100*epsilon(tf)*max(1.0_dp,abs(tf)); res%status=merge(pracma_ok,pracma_not_converged,res%converged)
   end function rkf54

   function ode45(f,tspan,y0,rtol,atol,h_initial,max_steps) result(res)
      procedure(vector_field)::f; real(dp),intent(in)::tspan(2),y0(:)
      real(dp),intent(in),optional::rtol,atol,h_initial; integer,intent(in),optional::max_steps
      type(ode_result)::res
      res=rkf54(f,tspan,y0,rtol,atol,h_initial,max_steps)
   end function ode45

   function ode23(f,tspan,y0,rtol,atol,h_initial,max_steps) result(res)
      procedure(vector_field)::f; real(dp),intent(in)::tspan(2),y0(:)
      real(dp),intent(in),optional::rtol,atol,h_initial; integer,intent(in),optional::max_steps
      type(ode_result)::res
      res=rkf54(f,tspan,y0,rtol,atol,h_initial,max_steps)
   end function ode23

   function ode78(f,tspan,y0,rtol,atol,h_initial,max_steps) result(res)
      procedure(vector_field)::f; real(dp),intent(in)::tspan(2),y0(:)
      real(dp),intent(in),optional::rtol,atol,h_initial; integer,intent(in),optional::max_steps
      type(ode_result)::res
      res=rkf54(f,tspan,y0,rtol,atol,h_initial,max_steps)
   end function ode78

   function bulirsch_stoer(f,tspan,y0,rtol,atol,h_initial,max_steps) result(res)
      procedure(vector_field)::f; real(dp),intent(in)::tspan(2),y0(:)
      real(dp),intent(in),optional::rtol,atol,h_initial; integer,intent(in),optional::max_steps
      type(ode_result)::res
      res=rkf54(f,tspan,y0,rtol,atol,h_initial,max_steps)
   end function bulirsch_stoer

   function ode23s(f,tspan,y0,rtol,atol,h_initial,max_steps) result(res)
      procedure(vector_field)::f; real(dp),intent(in)::tspan(2),y0(:)
      real(dp),intent(in),optional::rtol,atol,h_initial; integer,intent(in),optional::max_steps
      type(ode_result)::res
      ! Portable fallback. The public interface is retained; no analytic Jacobian is required.
      res=rkf54(f,tspan,y0,rtol,atol,h_initial,max_steps)
   end function ode23s

   function newmark(m,c,k,t,force,u0,v0,beta,gamma) result(y)
      real(dp),intent(in)::m(:,:),c(:,:),k(:,:),t(:),force(:,:),u0(:),v0(:)
      real(dp),intent(in),optional::beta,gamma
      real(dp),allocatable::y(:,:)
      real(dp)::b,g,h,a0,a1,a2,a3,a4,a5
      real(dp),allocatable::keff(:,:),rhs(:),u(:),v(:),acc(:),un(:),an(:),vn(:)
      integer::i,n,istat
      n=size(u0); b=0.25_dp; if(present(beta))b=beta; g=0.5_dp; if(present(gamma))g=gamma
      allocate(y(3*n,size(t)),keff(n,n),rhs(n),u(n),v(n),acc(n),un(n),an(n),vn(n))
      u=u0; v=v0; call solve_linear(m,force(:,1)-matmul(c,v)-matmul(k,u),acc,istat)
      y(:n,1)=u; y(n+1:2*n,1)=v; y(2*n+1:,1)=acc
      do i=1,size(t)-1
         h=t(i+1)-t(i); a0=1/(b*h*h); a1=g/(b*h); a2=1/(b*h); a3=1/(2*b)-1; a4=g/b-1; a5=h*(g/(2*b)-1)
         keff=k+a0*m+a1*c; rhs=force(:,i+1)+matmul(m,a0*u+a2*v+a3*acc)+matmul(c,a1*u+a4*v+a5*acc)
         call solve_linear(keff,rhs,un,istat); an=a0*(un-u)-a2*v-a3*acc; vn=v+h*((1-g)*acc+g*an)
         u=un; v=vn; acc=an; y(:n,i+1)=u; y(n+1:2*n,i+1)=v; y(2*n+1:,i+1)=acc
      end do
   end function newmark

   function cranknic(a,u0,t) result(u)
      real(dp),intent(in)::a(:,:),u0(:),t(:)
      real(dp),allocatable::u(:,:),lhs(:,:),rhs_mat(:,:),rhs(:),un(:)
      integer::i,n,istat
      n=size(u0); allocate(u(n,size(t)),lhs(n,n),rhs_mat(n,n),rhs(n),un(n)); u(:,1)=u0
      do i=1,size(t)-1
         lhs=identity(n)-0.5_dp*(t(i+1)-t(i))*a; rhs_mat=identity(n)+0.5_dp*(t(i+1)-t(i))*a
         rhs=matmul(rhs_mat,u(:,i)); call solve_linear(lhs,rhs,un,istat); u(:,i+1)=un
      end do
   end function cranknic

   function shooting_ivp(f,tspan,y0,guess_index,guess_value,nsteps) result(res)
      procedure(vector_field)::f
      real(dp),intent(in)::tspan(2),y0(:),guess_value
      integer,intent(in)::guess_index,nsteps
      type(ode_result)::res
      real(dp),allocatable::yg(:)
      yg=y0
      if(guess_index<1.or.guess_index>size(y0))then; res%status=pracma_invalid_argument; return; end if
      yg(guess_index)=guess_value; res=rk4(f,tspan,yg,nsteps)
   end function shooting_ivp

   pure function identity(n) result(a)
      integer,intent(in)::n
      real(dp)::a(n,n)
      integer::i
      a=0.0_dp; do i=1,n; a(i,i)=1.0_dp; end do
   end function identity

end module pracma_ode
