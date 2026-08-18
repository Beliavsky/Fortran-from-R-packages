! SPDX-License-Identifier: GPL-2.0-or-later
module desolve_utilities
  use desolve_kinds, only : dp
  use desolve_types, only : ode_result
  implicit none
  private

  abstract interface
    function scalar_fun(x) result(fx)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: fx
    end function scalar_fun
    subroutine map_fun(t,x,xnext)
      import dp
      real(dp),intent(in)::t,x(:)
      real(dp),intent(out)::xnext(:)
    end subroutine map_fun
  end interface

  type, public :: forcing_table
    real(dp), allocatable :: time(:)
    real(dp), allocatable :: value(:,:)
    integer :: method = 1 ! 1 linear, 2 constant
  contains
    procedure :: eval => forcing_eval
  end type forcing_table

  type, public :: event_table
    real(dp), allocatable :: time(:)
    integer, allocatable :: state(:)
    real(dp), allocatable :: value(:)
    integer, allocatable :: method(:) ! 1 replace, 2 add, 3 multiply
  contains
    procedure :: apply => event_apply
  end type event_table

  type, public :: history_buffer
    integer :: capacity=0, count=0
    real(dp), allocatable :: time(:), y(:,:), dy(:,:)
  contains
    procedure :: init => history_init
    procedure :: append => history_append
    procedure :: lag_value => history_lag_value
    procedure :: lag_deriv => history_lag_deriv
  end type history_buffer

  type, public :: sparsity_pattern
    integer, allocatable :: row_ptr(:)
    integer, allocatable :: col_ind(:)
  end type sparsity_pattern

  public :: brent_root, iterate_map, sparsity_2d, sparsity_3d
  public :: hermite_value, hermite_deriv, nearest_event, clean_event_times

contains

  function brent_root(f,ax,bx,tol,maxit,status) result(root)
    procedure(scalar_fun)::f
    real(dp),intent(in)::ax,bx
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::maxit
    integer,intent(out),optional::status
    real(dp)::root
    real(dp)::a,b,c,fa,fb,fc,prev_step,tol_act,p,q,new_step,t1,cb,t2,eps,tol0,tmp
    integer::it,nit
    a=ax;b=bx;fa=f(a);fb=f(b);c=a;fc=fa
    tol0=1e-10_dp;if(present(tol))tol0=tol
    nit=100;if(present(maxit))nit=maxit
    eps=epsilon(1.0_dp)
    if(abs(fa)<=tiny(1.0_dp))then;root=a;if(present(status))status=0;return;end if
    if(abs(fb)<=tiny(1.0_dp))then;root=b;if(present(status))status=0;return;end if
    if(fa*fb>0.0_dp)then;root=b;if(present(status))status=-1;return;end if
    do it=1,nit
      prev_step=b-a
      if(abs(fc)<abs(fb))then
        tmp=a;a=b;b=c;c=tmp
        tmp=fa;fa=fb;fb=fc;fc=tmp
      end if
      tol_act=2.0_dp*eps*abs(b)+tol0/2.0_dp;new_step=(c-b)/2.0_dp
      if(abs(new_step)<=tol_act.or.abs(fb)<=tiny(1.0_dp))then;root=b;if(present(status))status=0;return;end if
      if(abs(prev_step)>=tol_act.and.abs(fa)>abs(fb))then
        cb=c-b
        if(abs(a-c)<=eps*max(1.0_dp,abs(c)))then
          t1=fb/fa;p=cb*t1;q=1.0_dp-t1
        else
          q=fa/fc;t1=fb/fc;t2=fb/fa
          p=t2*(cb*q*(q-t1)-(b-a)*(t1-1.0_dp));q=(q-1.0_dp)*(t1-1.0_dp)*(t2-1.0_dp)
        end if
        if(p>0.0_dp)then;q=-q;else;p=-p;end if
        if(p<0.75_dp*cb*q-abs(tol_act*q)/2.0_dp.and.p<abs(prev_step*q/2.0_dp))new_step=p/q
      end if
      if(abs(new_step)<tol_act)new_step=sign(tol_act,new_step)
      a=b;fa=fb;b=b+new_step;fb=f(b)
      if((fb>0.0_dp.and.fc>0.0_dp).or.(fb<0.0_dp.and.fc<0.0_dp))then;c=a;fc=fa;end if
    end do
    root=b;if(present(status))status=1
  end function brent_root

  function forcing_eval(self,t) result(v)
    class(forcing_table),intent(in)::self
    real(dp),intent(in)::t
    real(dp),allocatable::v(:)
    integer::i,n
    real(dp)::w
    if(.not.allocated(self%time).or..not.allocated(self%value))error stop 'forcing table not initialized'
    if(size(self%value,2)/=size(self%time))error stop 'forcing table dimensions inconsistent'
    allocate(v(size(self%value,1)));n=size(self%time)
    if(t<=self%time(1))then;v=self%value(:,1);return;end if
    if(t>=self%time(n))then;v=self%value(:,n);return;end if
    i=locate_interval(self%time,t)
    if(self%method==1)then
      w=(t-self%time(i))/(self%time(i+1)-self%time(i));v=(1.0_dp-w)*self%value(:,i)+w*self%value(:,i+1)
    else;v=self%value(:,i);end if
  end function forcing_eval

  subroutine event_apply(self,t,y,n_applied)
    class(event_table),intent(in)::self
    real(dp),intent(in)::t
    real(dp),intent(inout)::y(:)
    integer,intent(out),optional::n_applied
    integer::i,n,s,na
    na=0;if(.not.allocated(self%time))then;if(present(n_applied))n_applied=0;return;end if
    n=size(self%time)
    do i=1,n
      if(abs(self%time(i)-t)<=100.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(t)))then
        s=self%state(i);if(s<1.or.s>size(y))cycle
        select case(self%method(i));case(1);y(s)=self%value(i);case(2);y(s)=y(s)+self%value(i);case(3);y(s)=y(s)*self%value(i)
        end select;na=na+1
      end if
    end do
    if(present(n_applied))n_applied=na
  end subroutine event_apply

  subroutine history_init(self,nvar,capacity)
    class(history_buffer),intent(inout)::self
    integer,intent(in)::nvar,capacity
    self%capacity=max(2,capacity);self%count=0
    allocate(self%time(self%capacity),self%y(nvar,self%capacity),self%dy(nvar,self%capacity))
    self%time=0.0_dp;self%y=0.0_dp;self%dy=0.0_dp
  end subroutine history_init

  subroutine history_append(self,t,y,dy)
    class(history_buffer),intent(inout)::self
    real(dp),intent(in)::t,y(:),dy(:)
    if(.not.allocated(self%time))call self%init(size(y),1000)
    if(size(y)/=size(self%y,1).or.size(dy)/=size(y))error stop 'history append size mismatch'
    if(self%count<self%capacity)then
      self%count=self%count+1
    else
      self%time(1:self%capacity-1)=self%time(2:self%capacity)
      self%y(:,1:self%capacity-1)=self%y(:,2:self%capacity)
      self%dy(:,1:self%capacity-1)=self%dy(:,2:self%capacity)
    end if
    self%time(self%count)=t;self%y(:,self%count)=y;self%dy(:,self%count)=dy
  end subroutine history_append

  function history_lag_value(self,t) result(v)
    class(history_buffer),intent(in)::self
    real(dp),intent(in)::t
    real(dp),allocatable::v(:)
    integer::i,n,j
    allocate(v(size(self%y,1)));n=self%count
    if(n<1)error stop 'history is empty'
    if(t<=self%time(1))then;v=self%y(:,1);return;end if
    if(t>=self%time(n))then;v=self%y(:,n)+self%dy(:,n)*(t-self%time(n));return;end if
    i=locate_interval(self%time(:n),t)
    do j=1,size(v);v(j)=hermite_value(self%time(i),self%time(i+1),self%y(j,i),self%y(j,i+1), &
      self%dy(j,i),self%dy(j,i+1),t);end do
  end function history_lag_value

  function history_lag_deriv(self,t) result(v)
    class(history_buffer),intent(in)::self
    real(dp),intent(in)::t
    real(dp),allocatable::v(:)
    integer::i,n,j
    allocate(v(size(self%y,1)));n=self%count
    if(n<1)error stop 'history is empty'
    if(t<=self%time(1))then;v=self%dy(:,1);return;end if
    if(t>=self%time(n))then;v=self%dy(:,n);return;end if
    i=locate_interval(self%time(:n),t)
    do j=1,size(v);v(j)=hermite_deriv(self%time(i),self%time(i+1),self%y(j,i),self%y(j,i+1), &
      self%dy(j,i),self%dy(j,i+1),t);end do
  end function history_lag_deriv

  pure real(dp) function hermite_value(t0,t1,y0,y1,dy0,dy1,t) result(res)
    real(dp),intent(in)::t0,t1,y0,y1,dy0,dy1,t
    real(dp)::tt0,tt1,tt12,tt02,hh
    tt0=t-t0;tt1=t-t1;tt12=tt1*tt1;tt02=tt0*tt0;hh=t1-t0
    if(abs(hh)>tiny(1.0_dp))then
      res=(dy0*tt0*tt12+dy1*tt1*tt02+(y0*(2.0_dp*tt0+hh)*tt12-y1*(2.0_dp*tt1-hh)*tt02)/hh)/(hh*hh)
    else;res=y0;end if
  end function hermite_value

  pure real(dp) function hermite_deriv(t0,t1,y0,y1,dy0,dy1,t) result(res)
    real(dp),intent(in)::t0,t1,y0,y1,dy0,dy1,t
    real(dp)::tt0,tt1,tt12,tt02,hh
    tt0=t-t0;tt1=t-t1;tt12=tt1*tt1;tt02=tt0*tt0;hh=t1-t0
    if(abs(hh)>tiny(1.0_dp))then
      res=(dy0*(tt12+2.0_dp*tt0*tt1)+dy1*(tt02+2.0_dp*tt0*tt1)+ &
        (y0*2.0_dp*tt1*(2.0_dp*tt0+hh+tt1)-y1*2.0_dp*tt0*(2.0_dp*tt1-hh+tt0))/hh)/(hh*hh)
    else;res=dy0;end if
  end function hermite_deriv

  function iterate_map(fun,x0,times,nsteps) result(sol)
    procedure(map_fun)::fun
    real(dp),intent(in)::x0(:),times(:)
    integer,intent(in),optional::nsteps
    type(ode_result)::sol
    integer::ns,i,j,n;real(dp)::t,dt
    real(dp),allocatable::x(:),xn(:)
    ns=1;if(present(nsteps))ns=max(1,nsteps);n=size(x0)
    allocate(sol%t(size(times)),sol%y(n,size(times)),x(n),xn(n));sol%t=times;sol%y=0.0_dp;sol%y(:,1)=x0;x=x0
    t=times(1)
    do i=2,size(times)
      dt=(times(i)-t)/real(ns,dp)
      do j=1,ns;call fun(t,x,xn);x=xn;t=t+dt;end do
      sol%y(:,i)=x
    end do
    sol%status=0;sol%message='success';sol%stats%n_steps=(size(times)-1)*ns;sol%stats%n_rhs=sol%stats%n_steps
  end function iterate_map

  function sparsity_2d(nspec,nx,ny,cyclic_x,cyclic_y,mask) result(sp)
    integer,intent(in)::nspec,nx,ny
    logical,intent(in)::cyclic_x,cyclic_y
    logical,intent(in),optional::mask(:)
    type(sparsity_pattern)::sp
    logical,allocatable::present_node(:)
    integer,allocatable::map(:),cols(:)
    integer::nt,tot,orig,newn,s,j,k,l,row,nc,m,nb
    nt=nx*ny;tot=nspec*nt;allocate(present_node(tot));present_node=.true.;if(present(mask))present_node=mask
    allocate(map(tot));map=0;newn=0
    do orig=1,tot;if(present_node(orig))then;newn=newn+1;map(orig)=newn;end if;end do
    allocate(sp%row_ptr(newn+1),cols(max(1,newn*(nspec+5))));sp%row_ptr(1)=1;nc=0;row=0
    do s=0,nspec-1;do j=0,nx-1;do k=0,ny-1
      m=s*nt+j*ny+k+1;if(map(m)==0)cycle;row=row+1
      call add_col(m);if(k<ny-1)call add_col(m+1);if(j<nx-1)call add_col(m+ny);if(j>0)call add_col(m-ny);if(k>0)call add_col(m-1)
      if(cyclic_x)then;if(j==0)call add_col(s*nt+(nx-1)*ny+k+1);if(j==nx-1)call add_col(s*nt+k+1);end if
      if(cyclic_y)then;if(k==0)call add_col(s*nt+(j+1)*ny);if(k==ny-1)call add_col(s*nt+j*ny+1);end if
      do l=0,nspec-1;if(l/=s)call add_col(l*nt+j*ny+k+1);end do
      sp%row_ptr(row+1)=nc+1
    end do;end do;end do
    allocate(sp%col_ind(nc));if(nc>0)sp%col_ind=cols(:nc)
  contains
    subroutine add_col(idx)
      integer,intent(in)::idx;integer::mapped
      if(idx<1.or.idx>tot)return;mapped=map(idx);if(mapped==0)return
      do nb=sp%row_ptr(row),nc;if(cols(nb)==mapped)return;end do
      nc=nc+1;if(nc>size(cols))error stop 'internal sparsity_2d capacity';cols(nc)=mapped
    end subroutine add_col
  end function sparsity_2d

  function sparsity_3d(nspec,nx,ny,nz,cyclic_x,cyclic_y,cyclic_z,mask) result(sp)
    integer,intent(in)::nspec,nx,ny,nz
    logical,intent(in)::cyclic_x,cyclic_y,cyclic_z
    logical,intent(in),optional::mask(:)
    type(sparsity_pattern)::sp
    logical,allocatable::present_node(:);integer,allocatable::map(:),cols(:)
    integer::nt,tot,orig,newn,s,j,k,l,q,row,nc,m,nb
    nt=nx*ny*nz;tot=nspec*nt;allocate(present_node(tot));present_node=.true.;if(present(mask))present_node=mask
    allocate(map(tot));map=0;newn=0;do orig=1,tot;if(present_node(orig))then;newn=newn+1;map(orig)=newn;end if;end do
    allocate(sp%row_ptr(newn+1),cols(max(1,newn*(nspec+7))));sp%row_ptr(1)=1;nc=0;row=0
    do s=0,nspec-1;do j=0,nx-1;do k=0,ny-1;do l=0,nz-1
      m=s*nt+j*ny*nz+k*nz+l+1;if(map(m)==0)cycle;row=row+1;call add_col(m)
      if(l<nz-1)call add_col(m+1);if(l>0)call add_col(m-1)
      if(k<ny-1)call add_col(m+nz);if(k>0)call add_col(m-nz)
      if(j<nx-1)call add_col(m+ny*nz);if(j>0)call add_col(m-ny*nz)
      if(cyclic_z)then;if(l==0)call add_col(s*nt+j*ny*nz+k*nz+nz);if(l==nz-1)call add_col(s*nt+j*ny*nz+k*nz+1);end if
      if(cyclic_y)then;if(k==0)call add_col(s*nt+j*ny*nz+(ny-1)*nz+l+1);if(k==ny-1)call add_col(s*nt+j*ny*nz+l+1);end if
      if(cyclic_x)then;if(j==0)call add_col(s*nt+(nx-1)*ny*nz+k*nz+l+1);if(j==nx-1)call add_col(s*nt+k*nz+l+1);end if
      do q=0,nspec-1;if(q/=s)call add_col(q*nt+j*ny*nz+k*nz+l+1);end do;sp%row_ptr(row+1)=nc+1
    end do;end do;end do;end do
    allocate(sp%col_ind(nc));if(nc>0)sp%col_ind=cols(:nc)
  contains
    subroutine add_col(idx)
      integer,intent(in)::idx;integer::mapped
      if(idx<1.or.idx>tot)return;mapped=map(idx);if(mapped==0)return
      do nb=sp%row_ptr(row),nc;if(cols(nb)==mapped)return;end do
      nc=nc+1;if(nc>size(cols))error stop 'internal sparsity_3d capacity';cols(nc)=mapped
    end subroutine add_col
  end function sparsity_3d

  pure real(dp) function nearest_event(t,event_times) result(te)
    real(dp),intent(in)::t,event_times(:);integer::i,k;real(dp)::d,best
    if(size(event_times)==0)then;te=huge(1.0_dp);return;end if
    k=1;best=abs(event_times(1)-t)
    do i=2,size(event_times);d=abs(event_times(i)-t);if(d<best)then;best=d;k=i;end if;end do;te=event_times(k)
  end function nearest_event

  function clean_event_times(times,event_times,tol) result(out)
    real(dp),intent(in)::times(:),event_times(:);real(dp),intent(in),optional::tol
    real(dp),allocatable::out(:),tmp(:);real(dp)::eps0;integer::i,j,n
    eps0=100.0_dp*epsilon(1.0_dp);if(present(tol))eps0=tol
    allocate(tmp(size(times)+size(event_times)));n=0
    do i=1,size(times);n=n+1;tmp(n)=times(i);end do
    do i=1,size(event_times)
      if(all(abs(times-event_times(i))>eps0*max(1.0_dp,abs(event_times(i)))))then;n=n+1;tmp(n)=event_times(i);end if
    end do
    call sort_real(tmp(:n));j=1
    do i=2,n;if(abs(tmp(i)-tmp(j))>eps0*max(1.0_dp,abs(tmp(i))))then;j=j+1;tmp(j)=tmp(i);end if;end do
    allocate(out(j));out=tmp(:j)
  end function clean_event_times

  integer function locate_interval(x,t) result(k)
    real(dp),intent(in)::x(:),t;integer::lo,hi,mid
    lo=1;hi=size(x)
    do while(hi-lo>1);mid=(lo+hi)/2;if(t>=x(mid))then;lo=mid;else;hi=mid;end if;end do;k=lo
  end function locate_interval

  subroutine sort_real(x)
    real(dp),intent(inout)::x(:);integer::i,j;real(dp)::v
    do i=2,size(x);v=x(i);j=i-1;do while(j>=1);if(x(j)<=v)exit;x(j+1)=x(j);j=j-1;end do;x(j+1)=v;end do
  end subroutine sort_real
end module desolve_utilities
