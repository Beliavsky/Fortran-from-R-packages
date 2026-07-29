! SPDX-License-Identifier: MIT
module bayesianou_geometry
  use bayesianou_kinds, only : dp, status_ok, status_bad_input, status_singular, pi
  use bayesianou_math, only : rng_state, rng_seed, rng_uniform, rng_normal, &
                              numerical_hessian, symmetric_eigen, invert_spd, finite_all
  implicit none
  private

  abstract interface
    function geom_log_prob(theta) result(lp)
      import dp
      real(dp), intent(in) :: theta(:)
      real(dp) :: lp
    end function geom_log_prob
    subroutine geom_gradient(theta,g)
      import dp
      real(dp), intent(in) :: theta(:)
      real(dp), intent(out) :: g(size(theta))
    end subroutine geom_gradient
    subroutine geom_hessian(theta,h)
      import dp
      real(dp), intent(in) :: theta(:)
      real(dp), intent(out) :: h(size(theta),size(theta))
    end subroutine geom_hessian
  end interface

  type, public :: ou_geom_target_type
    integer :: dim = 0
    procedure(geom_log_prob), pointer, nopass :: log_prob => null()
    procedure(geom_gradient), pointer, nopass :: gradient => null()
    procedure(geom_hessian), pointer, nopass :: hessian => null()
  end type ou_geom_target_type

  type, public :: ou_geom_metric_type
    integer :: dim = 0
    logical :: position_dependent = .false.
    integer :: curvature = 0       ! 0 Euclidean; 1 SoftAbs
    real(dp), allocatable :: mass0(:,:)
    real(dp) :: alpha = 1.0e6_dp
    real(dp) :: eigen_floor = 1.0e-8_dp
    real(dp) :: fd_step = 1.0e-4_dp
  end type ou_geom_metric_type

  type, public :: ou_geom_hmc_result
    real(dp), allocatable :: draws(:,:), energy(:)
    real(dp) :: accept_rate = 0.0_dp
    real(dp) :: ebfmi = 0.0_dp
    integer :: n_divergent = 0
    integer :: status = status_ok
  end type ou_geom_hmc_result

  public :: ou_geom_target, ou_geom_metric_euclidean, ou_geom_metric_riemannian
  public :: ou_geom_hmc, ou_geom_mass, ou_geom_ebfmi

contains

  function ou_geom_target(dim,log_prob,gradient,hessian) result(target)
    integer, intent(in) :: dim
    procedure(geom_log_prob), pointer, intent(in) :: log_prob
    procedure(geom_gradient), pointer, intent(in) :: gradient
    procedure(geom_hessian), pointer, intent(in), optional :: hessian
    type(ou_geom_target_type) :: target
    target%dim=dim;target%log_prob=>log_prob;target%gradient=>gradient
    if(present(hessian)) target%hessian=>hessian
  end function ou_geom_target

  function ou_geom_metric_euclidean(dim,mass) result(metric)
    integer, intent(in) :: dim
    real(dp), intent(in), optional :: mass(:,:)
    type(ou_geom_metric_type) :: metric
    integer :: i
    metric%dim=dim;metric%position_dependent=.false.;metric%curvature=0
    allocate(metric%mass0(dim,dim));metric%mass0=0.0_dp
    if(present(mass)) then
      metric%mass0=mass
    else
      do i=1,dim;metric%mass0(i,i)=1.0_dp;end do
    end if
  end function ou_geom_metric_euclidean

  function ou_geom_metric_riemannian(target,alpha,eigen_floor,fd_step) result(metric)
    type(ou_geom_target_type), intent(in) :: target
    real(dp), intent(in), optional :: alpha,eigen_floor,fd_step
    type(ou_geom_metric_type) :: metric
    metric%dim=target%dim;metric%position_dependent=.true.;metric%curvature=1
    if(present(alpha))metric%alpha=alpha
    if(present(eigen_floor))metric%eigen_floor=eigen_floor
    if(present(fd_step))metric%fd_step=fd_step
  end function ou_geom_metric_riemannian

  subroutine ou_geom_mass(target,metric,theta,mass,inv_mass,logdet,status)
    type(ou_geom_target_type), intent(in) :: target
    type(ou_geom_metric_type), intent(in) :: metric
    real(dp), intent(in) :: theta(:)
    real(dp), intent(out) :: mass(size(theta),size(theta)),inv_mass(size(theta),size(theta)),logdet
    integer, intent(out) :: status
    real(dp) :: h(size(theta),size(theta)),vals(size(theta)),vecs(size(theta),size(theta)),soft(size(theta))
    integer :: i
    if(.not.metric%position_dependent) then
      mass=metric%mass0
    else
      call target_hessian(target,theta,h)
      h=-0.5_dp*(h+transpose(h))
      call symmetric_eigen(h,vals,vecs,status)
      if(status/=status_ok) then;mass=0;inv_mass=0;logdet=0;return;end if
      do i=1,size(vals)
        if(abs(metric%alpha*vals(i))<1.0e-5_dp) then
          soft(i)=1.0_dp/metric%alpha + metric%alpha*vals(i)**2/3.0_dp
        else if(abs(metric%alpha*vals(i))>40.0_dp) then
          soft(i)=abs(vals(i))
        else
          soft(i)=vals(i)/tanh(metric%alpha*vals(i))
        end if
        soft(i)=max(soft(i),metric%eigen_floor)
      end do
      mass=matmul(vecs,matmul(diag_matrix(soft),transpose(vecs)))
    end if
    call invert_spd(mass,inv_mass,status,logdet)
  end subroutine ou_geom_mass

  subroutine target_hessian(target,theta,h)
    type(ou_geom_target_type), intent(in) :: target
    real(dp), intent(in) :: theta(:)
    real(dp), intent(out) :: h(size(theta),size(theta))
    real(dp) :: gp(size(theta)),gm(size(theta)),xp(size(theta)),xm(size(theta)),delta
    integer :: j
    if(associated(target%hessian)) then
      call target%hessian(theta,h)
    else
      do j=1,size(theta)
        delta=1.0e-4_dp*(1.0_dp+abs(theta(j)))
        xp=theta;xm=theta;xp(j)=xp(j)+delta;xm(j)=xm(j)-delta
        call target%gradient(xp,gp);call target%gradient(xm,gm)
        h(:,j)=(gp-gm)/(2.0_dp*delta)
      end do
      h=0.5_dp*(h+transpose(h))
    end if
  end subroutine target_hessian

  function diag_matrix(x) result(a)
    real(dp), intent(in) :: x(:)
    real(dp) :: a(size(x),size(x))
    integer :: i
    a=0.0_dp;do i=1,size(x);a(i,i)=x(i);end do
  end function diag_matrix

  subroutine draw_momentum(rng,mass,p,status)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: mass(:,:)
    real(dp), intent(out) :: p(size(mass,1))
    integer, intent(out) :: status
    real(dp) :: l(size(mass,1),size(mass,2)),z(size(mass,1))
    integer :: n,info,i
    interface
      subroutine dpotrf(uplo,n,a,lda,info)
        import dp
        character(len=1),intent(in)::uplo
        integer,intent(in)::n,lda
        integer,intent(out)::info
        real(dp),intent(inout)::a(lda,*)
      end subroutine dpotrf
    end interface
    n=size(p);l=mass;call dpotrf('L',n,l,n,info)
    if(info/=0)then;status=status_singular;p=0;return;end if
    do i=1,n;z(i)=rng_normal(rng);end do
    do i=1,n
      if(i<n)l(i,i+1:n)=0.0_dp
    end do
    p=matmul(l,z);status=status_ok
  end subroutine draw_momentum

  function kinetic_value(p,inv_mass,logdet) result(k)
    real(dp), intent(in) :: p(:),inv_mass(:,:),logdet
    real(dp)::k
    k=0.5_dp*dot_product(p,matmul(inv_mass,p))+0.5_dp*logdet
  end function kinetic_value

  subroutine metric_derivatives(target,metric,theta,dm,status)
    type(ou_geom_target_type),intent(in)::target
    type(ou_geom_metric_type),intent(in)::metric
    real(dp),intent(in)::theta(:)
    real(dp),intent(out)::dm(size(theta),size(theta),size(theta))
    integer,intent(out)::status
    real(dp)::xp(size(theta)),xm(size(theta)),mp(size(theta),size(theta)),mm(size(theta),size(theta))
    real(dp)::ip(size(theta),size(theta)),im(size(theta),size(theta)),ld,d
    integer::j,st1,st2
    if(.not.metric%position_dependent)then;dm=0;status=status_ok;return;end if
    do j=1,size(theta)
      d=metric%fd_step*(1+abs(theta(j)));xp=theta;xm=theta;xp(j)=xp(j)+d;xm(j)=xm(j)-d
      call ou_geom_mass(target,metric,xp,mp,ip,ld,st1)
      call ou_geom_mass(target,metric,xm,mm,im,ld,st2)
      if(st1/=status_ok.or.st2/=status_ok)then;status=status_singular;dm=0;return;end if
      dm(:,:,j)=(mp-mm)/(2*d)
    end do
    status=status_ok
  end subroutine metric_derivatives

  subroutine hamiltonian_gradient_theta(target,metric,theta,p,g,status)
    type(ou_geom_target_type),intent(in)::target
    type(ou_geom_metric_type),intent(in)::metric
    real(dp),intent(in)::theta(:),p(:)
    real(dp),intent(out)::g(size(theta))
    integer,intent(out)::status
    real(dp)::mass(size(theta),size(theta)),inv(size(theta),size(theta)),ld
    real(dp)::dm(size(theta),size(theta),size(theta)),glp(size(theta)),v(size(theta))
    integer::j
    call ou_geom_mass(target,metric,theta,mass,inv,ld,status);if(status/=status_ok)return
    call target%gradient(theta,glp)
    if(.not.metric%position_dependent)then;g=-glp;return;end if
    call metric_derivatives(target,metric,theta,dm,status);if(status/=status_ok)return
    v=matmul(inv,p)
    do j=1,size(theta)
      g(j)=-glp(j)+0.5_dp*sum(inv*transpose(dm(:,:,j)))-0.5_dp*dot_product(v,matmul(dm(:,:,j),v))
    end do
  end subroutine hamiltonian_gradient_theta

  subroutine leapfrog(target,metric,theta,p,epsilon,n_steps,converged)
    type(ou_geom_target_type),intent(in)::target
    type(ou_geom_metric_type),intent(in)::metric
    real(dp),intent(inout)::theta(:),p(:)
    real(dp),intent(in)::epsilon
    integer,intent(in)::n_steps
    logical,intent(out)::converged
    real(dp)::g(size(theta)),mass(size(theta),size(theta)),inv(size(theta),size(theta)),ld
    real(dp)::p0(size(theta)),ph(size(theta)),pn(size(theta)),th(size(theta)),tn(size(theta)),dr0(size(theta))
    integer::i,it,status
    real(dp)::delta
    converged=.true.
    do i=1,n_steps
      if(.not.metric%position_dependent)then
        call target%gradient(theta,g);p=p+0.5_dp*epsilon*g
        call ou_geom_mass(target,metric,theta,mass,inv,ld,status);if(status/=status_ok)then;converged=.false.;return;end if
        theta=theta+epsilon*matmul(inv,p)
        call target%gradient(theta,g);p=p+0.5_dp*epsilon*g
      else
        p0=p;ph=p
        do it=1,100
          call hamiltonian_gradient_theta(target,metric,theta,ph,g,status)
          if(status/=status_ok)then;converged=.false.;return;end if
          pn=p0-0.5_dp*epsilon*g;delta=maxval(abs(pn-ph));ph=pn
          if(delta<1e-9_dp)exit
        end do
        if(it>=100)converged=.false.
        call ou_geom_mass(target,metric,theta,mass,inv,ld,status);if(status/=status_ok)then;converged=.false.;return;end if
        dr0=matmul(inv,ph);th=theta
        do it=1,100
          call ou_geom_mass(target,metric,th,mass,inv,ld,status);if(status/=status_ok)then;converged=.false.;return;end if
          tn=theta+0.5_dp*epsilon*(dr0+matmul(inv,ph));delta=maxval(abs(tn-th));th=tn
          if(delta<1e-9_dp)exit
        end do
        if(it>=100)converged=.false.
        theta=th
        call hamiltonian_gradient_theta(target,metric,theta,ph,g,status)
        if(status/=status_ok)then;converged=.false.;return;end if
        p=ph-0.5_dp*epsilon*g
      end if
      if(.not.finite_all(theta).or..not.finite_all(p))then;converged=.false.;return;end if
    end do
  end subroutine leapfrog

  subroutine ou_geom_hmc(target,metric,epsilon,n_steps,n_iter,n_warmup,init,seed,result)
    type(ou_geom_target_type),intent(in)::target
    type(ou_geom_metric_type),intent(in)::metric
    real(dp),intent(in)::epsilon
    integer,intent(in)::n_steps,n_iter,n_warmup,seed
    real(dp),intent(in)::init(:)
    type(ou_geom_hmc_result),intent(out)::result
    type(rng_state)::rng
    real(dp)::theta(size(init)),proposal(size(init)),p0(size(init)),p1(size(init))
    real(dp)::mass(size(init),size(init)),inv(size(init),size(init)),ld,h0,h1,dh,u
    integer::it,total,status,kept,accepted
    logical::conv
    if(size(init)/=target%dim.or.metric%dim/=target%dim)then;result%status=status_bad_input;return;end if
    call rng_seed(rng,seed);theta=init;total=n_iter+n_warmup;kept=0;accepted=0
    allocate(result%draws(n_iter,target%dim),result%energy(n_iter))
    do it=1,total
      call ou_geom_mass(target,metric,theta,mass,inv,ld,status)
      if(status/=status_ok)then;result%status=status;return;end if
      call draw_momentum(rng,mass,p0,status);proposal=theta;p1=p0
      h0=-target%log_prob(theta)+kinetic_value(p0,inv,ld)
      call leapfrog(target,metric,proposal,p1,epsilon,n_steps,conv)
      if(conv)then
        call ou_geom_mass(target,metric,proposal,mass,inv,ld,status)
        if(status==status_ok)then;h1=-target%log_prob(proposal)+kinetic_value(p1,inv,ld);else;h1=huge(1.0_dp);end if
      else;h1=huge(1.0_dp);end if
      dh=h1-h0;u=log(rng_uniform(rng))
      if(conv.and.abs(dh)<=1000.0_dp.and.u< -dh)then;theta=proposal;accepted=accepted+1;end if
      if(.not.conv.or.abs(dh)>1000.0_dp.or..not.finite_all([dh]))result%n_divergent=result%n_divergent+1
      if(it>n_warmup)then;kept=kept+1;result%draws(kept,:)=theta;result%energy(kept)=h0;end if
    end do
    result%accept_rate=real(accepted,dp)/real(total,dp);result%ebfmi=ou_geom_ebfmi(result%energy);result%status=status_ok
  end subroutine ou_geom_hmc

  function ou_geom_ebfmi(energy) result(v)
    real(dp),intent(in)::energy(:)
    real(dp)::v,m,den
    if(size(energy)<2)then;v=0;return;end if
    m=sum(energy)/real(size(energy),dp);den=sum((energy-m)**2)
    if(den<=0)then;v=0;else;v=sum((energy(2:)-energy(:size(energy)-1))**2)/den;end if
  end function ou_geom_ebfmi

end module bayesianou_geometry
