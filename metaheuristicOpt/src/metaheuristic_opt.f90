! Modern Fortran translation of the computational code in R package metaheuristicOpt 2.0.0.
! SPDX-License-Identifier: GPL-2.0-or-later
module metaheuristic_opt
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use mh_support, only : dp, mh_pi, mh_rng, shuffle_int, sample_distinct, weighted_index
   implicit none
   private
   abstract interface
      function objective_function(x) result(f)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: f
      end function objective_function
   end interface
   type, public :: mh_control
      integer :: num_population = 40
      integer :: max_iter = 500
      integer(int64) :: seed = 1234567_int64
      logical :: maximize = .false.
      logical :: legacy_quirks = .true.
      ! PSO
      real(dp) :: vmax = 2.0_dp
      real(dp) :: ci = 1.49445_dp
      real(dp) :: cg = 1.49445_dp
      real(dp) :: w = 0.729_dp
      ! FFA
      real(dp) :: b0 = 1.0_dp
      real(dp) :: gamma = 1.0_dp
      real(dp) :: alpha_ffa = 0.2_dp
      ! GA
      real(dp) :: pm = 0.1_dp
      real(dp) :: pc = 0.8_dp
      ! HS
      real(dp) :: par = 0.3_dp
      real(dp) :: hmcr = 0.95_dp
      real(dp) :: bandwidth = 0.05_dp
      ! CLONALG
      integer :: selection_size = 0
      real(dp) :: multiplication_factor = 0.5_dp
      real(dp) :: hypermutation_rate = 0.1_dp
      ! DE
      real(dp) :: scaling_vector = 0.8_dp
      real(dp) :: crossover_rate = 0.5_dp
      character(len=24) :: de_strategy = 'best 1'
      ! SFL
      integer :: num_memeplex = 0
      integer :: frog_leaping_iteration = 10
      ! CSO
      real(dp) :: mixture_ratio = 0.5_dp
      real(dp) :: tracing_constant = 0.1_dp
      real(dp) :: maximum_velocity = 1.0_dp
      integer :: smp = 20
      real(dp) :: srd = 20.0_dp
      integer :: cdc = 0
      logical :: spc = .true.
      ! ABC
      integer :: cycle_limit = 0
      ! KH
      real(dp) :: max_motion_induced = 0.01_dp
      real(dp) :: inertia_motion = 0.01_dp
      real(dp) :: epsilon = 1.0e-5_dp
      real(dp) :: foraging_speed = 0.02_dp
      real(dp) :: inertia_foraging = 0.01_dp
      real(dp) :: max_diffusion_speed = 0.01_dp
      real(dp) :: constant_space = 1.0_dp
      real(dp) :: mu = 0.1_dp
      ! CS
      real(dp) :: abandoned_fraction = 0.5_dp
      ! BA
      real(dp) :: max_frequency = 0.1_dp
      real(dp) :: min_frequency = -0.1_dp
      real(dp) :: gama = 1.0_dp
      real(dp) :: alpha_ba = 0.1_dp
      ! GBS
      real(dp) :: gravitational_const = -1.0_dp
      real(dp) :: kbest = 0.1_dp
   end type mh_control
   type, public :: mh_result
      real(dp), allocatable :: par(:)
      real(dp) :: value = huge(1.0_dp)
      real(dp), allocatable :: history(:)
      integer :: iterations = 0
      integer :: evaluations = 0
      character(len=16) :: algorithm = ''
   end type mh_result
   public :: dp, objective_function, metaopt, metaopt_many
   public :: pso, alo, gwo, da, ffa, ga, goa, hs, mfo, sca, woa
   public :: clonalg, de, sfl, cso, abc, kh, cs, ba, gbs, bho
   public :: sphere, schwefel, rastrigin, cantilever_beam, centilever_beam, mae
contains
subroutine metaopt(algorithm, fun, lower, upper, result, control)
   character(len=*), intent(in) :: algorithm
   procedure(objective_function) :: fun
   real(dp), intent(in) :: lower(:), upper(:)
   type(mh_result), intent(out) :: result
   type(mh_control), intent(in), optional :: control
   type(mh_control) :: c
   character(len=:), allocatable :: a
   c = mh_control()
   if (present(control)) c = control
   a = upper_string(trim(adjustl(algorithm)))
   select case(a)
   case('PSO')
    call pso(fun,lower,upper,result,c)
   case('ALO')
    call alo(fun,lower,upper,result,c)
   case('GWO')
    call gwo(fun,lower,upper,result,c)
   case('DA')
    call da(fun,lower,upper,result,c)
   case('FFA')
    call ffa(fun,lower,upper,result,c)
   case('GA')
    call ga(fun,lower,upper,result,c)
   case('GOA')
    call goa(fun,lower,upper,result,c)
   case('HS','IHS')
    call hs(fun,lower,upper,result,c)
   case('MFO')
    call mfo(fun,lower,upper,result,c)
   case('SCA')
    call sca(fun,lower,upper,result,c)
   case('WOA')
    call woa(fun,lower,upper,result,c)
   case('CLONALG')
    call clonalg(fun,lower,upper,result,c)
   case('DE')
    call de(fun,lower,upper,result,c)
   case('SFL')
    call sfl(fun,lower,upper,result,c)
   case('CSO')
    call cso(fun,lower,upper,result,c)
   case('ABC')
    call abc(fun,lower,upper,result,c)
   case('KH')
    call kh(fun,lower,upper,result,c)
   case('CS')
    call cs(fun,lower,upper,result,c)
   case('BA')
    call ba(fun,lower,upper,result,c)
   case('GBS')
    call gbs(fun,lower,upper,result,c)
   case('BHO')
    call bho(fun,lower,upper,result,c)
   case default
      error stop 'metaopt: unknown algorithm'
   end select
end subroutine metaopt

subroutine metaopt_many(algorithms, fun, lower, upper, results, control)
   character(len=*), intent(in) :: algorithms(:)
   procedure(objective_function) :: fun
   real(dp), intent(in) :: lower(:), upper(:)
   type(mh_result), intent(out) :: results(size(algorithms))
   type(mh_control), intent(in), optional :: control
   type(mh_control) :: c
   integer :: i

   c = mh_control()
   if (present(control)) c = control
   do i = 1, size(algorithms)
      ! The R metaOpt wrapper resets set.seed(seed) for each requested method.
      call metaopt(trim(algorithms(i)), fun, lower, upper, results(i), c)
   end do
end subroutine metaopt_many
function upper_string(s) result(out)
   character(len=*), intent(in) :: s
   character(len=len(s)) :: out
   integer :: i,k
   out=s
   do i=1,len(s)
      k=iachar(out(i:i))
      if (k>=iachar('a') .and. k<=iachar('z')) out(i:i)=achar(k-32)
   end do
end function upper_string
real(dp) function internal_value(fun, x, c) result(v)
   procedure(objective_function) :: fun
   real(dp), intent(in) :: x(:)
   type(mh_control), intent(in) :: c
   if (c%maximize) then
      v = -fun(x)
   else
      v = fun(x)
   end if
end function internal_value
real(dp) function sign_value(c) result(s)
   type(mh_control), intent(in) :: c
   if (c%maximize) then
      s=-1.0_dp
   else
      s=1.0_dp
   end if
end function sign_value
subroutine validate_bounds(lower,upper)
   real(dp),intent(in)::lower(:),upper(:)
   if(size(lower)/=size(upper)) error stop 'lower/upper sizes differ'
   if(size(lower)<1) error stop 'zero-dimensional problem'
   if(any(upper<lower)) error stop 'upper bound below lower bound'
end subroutine validate_bounds
subroutine init_result(result,n,c,name)
   type(mh_result),intent(out)::result
   integer,intent(in)::n
   type(mh_control),intent(in)::c
   character(len=*),intent(in)::name
   allocate(result%par(n), result%history(max(0,c%max_iter)))
   result%par=0.0_dp
    result%history=0.0_dp
   result%value=huge(1.0_dp)
    result%iterations=0
    result%evaluations=0
   result%algorithm=name
end subroutine init_result
subroutine random_population(rng,x,lower,upper)
   type(mh_rng),intent(inout)::rng
   real(dp),intent(out)::x(:,:)
   real(dp),intent(in)::lower(:),upper(:)
   integer::i,j
   do i=1,size(x,1)
      do j=1,size(x,2)
         x(i,j)=lower(j)+rng%uniform()*(upper(j)-lower(j))
      end do
   end do
end subroutine random_population
subroutine clamp_vector(x,lower,upper)
   real(dp),intent(inout)::x(:)
   real(dp),intent(in)::lower(:),upper(:)
   x=max(lower,min(upper,x))
end subroutine clamp_vector
subroutine clamp_population(x,lower,upper)
   real(dp),intent(inout)::x(:,:)
   real(dp),intent(in)::lower(:),upper(:)
   integer::i
   do i=1,size(x,1)
    call clamp_vector(x(i,:),lower,upper)
    end do
end subroutine clamp_population
subroutine eval_population(fun,c,x,fit,evals)
   procedure(objective_function)::fun
   type(mh_control),intent(in)::c
   real(dp),intent(in)::x(:,:)
   real(dp),intent(out)::fit(:)
   integer,intent(inout)::evals
   integer::i
   do i=1,size(x,1)
      fit(i)=internal_value(fun,x(i,:),c)
       evals=evals+1
   end do
end subroutine eval_population
integer function argmin(v) result(k)
   real(dp),intent(in)::v(:)
   integer::i
   k=1
   do i=2,size(v)
      if(v(i)<v(k)) k=i
   end do
end function argmin
integer function argmax(v) result(k)
   real(dp),intent(in)::v(:)
   integer::i
   k=1
   do i=2,size(v)
      if(v(i)>v(k)) k=i
   end do
end function argmax
subroutine sort_population(x,fit)
   real(dp),intent(inout)::x(:,:),fit(:)
   real(dp),allocatable::row(:)
   real(dp)::fv
   integer::i,j
   allocate(row(size(x,2)))
   do i=2,size(fit)
      fv=fit(i)
       row=x(i,:)
       j=i-1
      do while(j>=1)
         if(fit(j)<=fv) exit
         fit(j+1)=fit(j)
          x(j+1,:)=x(j,:)
          j=j-1
      end do
      fit(j+1)=fv
       x(j+1,:)=row
   end do
end subroutine sort_population
subroutine update_best(fun,c,x,fit,best,fbest,evals)
   procedure(objective_function)::fun
   type(mh_control),intent(in)::c
   real(dp),intent(in)::x(:,:)
   real(dp),intent(in),optional::fit(:)
   real(dp),intent(inout)::best(:),fbest
   integer,intent(inout)::evals
   real(dp)::f
   integer::i
   if(present(fit)) then
      i=argmin(fit)
      if(fit(i)<fbest) then
       fbest=fit(i)
      best=x(i,:)
      end if
   else
      do i=1,size(x,1)
         f=internal_value(fun,x(i,:),c)
         evals=evals+1
         if(f<fbest) then
         fbest=f
         best=x(i,:)
         end if
      end do
   end if
end subroutine update_best
subroutine finish_result(result,c,best,fbest,it,evals)
   type(mh_result),intent(inout)::result
   type(mh_control),intent(in)::c
   real(dp),intent(in)::best(:),fbest
   integer,intent(in)::it,evals
   result%par=best
   result%value=sign_value(c)*fbest
   result%iterations=it
   result%evaluations=evals
end subroutine finish_result
real(dp) function sphere(x) result(f)
   real(dp),intent(in)::x(:)
    f=sum(x*x)
end function sphere
real(dp) function schwefel(x) result(f)
   real(dp),intent(in)::x(:)
    f=sum(-x*sin(sqrt(abs(x))))
end function schwefel
real(dp) function rastrigin(x) result(f)
   real(dp),intent(in)::x(:)
    f=sum(x*x-10.0_dp*cos(2.0_dp*mh_pi*x)+10.0_dp)
end function rastrigin
real(dp) function cantilever_beam(x) result(f)
   real(dp),intent(in)::x(:)
   real(dp),parameter::u5(5)=[61.0_dp,37.0_dp,19.0_dp,7.0_dp,1.0_dp]
   real(dp)::pen
   if(size(x)/=5) error stop 'cantilever_beam requires 5 variables'
   pen=max(sum(u5/(x*x*x))-1.0_dp,0.0_dp)
   f=0.6224_dp*sum(x)+1.0e18_dp*pen
end function cantilever_beam

real(dp) function centilever_beam(x) result(f)
   real(dp),intent(in)::x(:)
   f=cantilever_beam(x)
end function centilever_beam

real(dp) function mae(data,x) result(f)
   real(dp),intent(in)::data(:,:),x(:)
   real(dp)::pred
   integer::i
   if(size(data,2)<4) error stop 'mae: data must have at least four columns'
   if(size(x)/=4) error stop 'mae: coefficient vector must have length four'
   f=0.0_dp
   do i=1,size(data,1)
      pred=data(i,1)*x(1)+data(i,2)*x(2)+data(i,3)*x(3)+x(4)
      f=f+abs(pred-data(i,4))
   end do
   if(size(data,1)>0) f=f/real(size(data,1),dp)
end function mae
subroutine pso(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),v(:,:),pbest(:,:),pfit(:)
   real(dp),allocatable::gbest(:)
   real(dp)::fbest,fnew,ri,rg,newv
   integer::np,n,i,j,t,evals
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   if(np<1) error stop 'PSO: num_population must be positive'
   call init_result(result,n,c,'PSO')
    call rng%seed(c%seed)
   allocate(x(np,n),v(np,n),pbest(np,n),pfit(np),gbest(n))
   call random_population(rng,x,lower,upper)
   do i=1,np
   do j=1,n
   v(i,j)=-c%vmax+2.0_dp*c%vmax*rng%uniform()
   end do
   end do
   pbest=x
   evals=0
   call eval_population(fun,c,x,pfit,evals)
   if (c%legacy_quirks) then
      ! The R wrapper initializes Gbest with calcBest(FUN, optimType, ...),
      ! whose which.max convention selects the worst internal fitness.
      i=argmax(pfit)
   else
      i=argmin(pfit)
   end if
   gbest=x(i,:)
   fbest=pfit(i)
   do t=1,c%max_iter
      do i=1,np
         do j=1,n
            ri=rng%uniform()
            rg=rng%uniform()
            newv=c%w*v(i,j)+c%ci*ri*(pbest(i,j)-x(i,j))+c%cg*rg*(gbest(j)-x(i,j))
            v(i,j)=max(-c%vmax,min(c%vmax,newv))
            x(i,j)=max(lower(j),min(upper(j),x(i,j)+v(i,j)))
            ! R source evaluates the complete position after every coordinate update.
            fnew=internal_value(fun,x(i,:),c)
            evals=evals+1
            if(fnew<pfit(i)) then
               pbest(i,:)=x(i,:)
               pfit(i)=fnew
               if(fnew<fbest) then
               fbest=fnew
               gbest=x(i,:)
               end if
            end if
         end do
      end do
      result%history(t)=sign_value(c)*fbest
   end do
   call finish_result(result,c,gbest,fbest,c%max_iter,evals)
end subroutine pso
subroutine gwo(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
    type(mh_rng)::rng
   real(dp),allocatable::x(:,:),fit(:),alpha(:),beta(:),delta(:)
   real(dp)::fa,fb,fd,a,r1,r2,a1,c1,a2,c2,a3,c3,x1,x2,x3,f
   integer::n,np,i,j,t,evals
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   if(np<3) error stop 'GWO: population must be at least 3'
   call init_result(result,n,c,'GWO')
   call rng%seed(c%seed)
   allocate(x(np,n),fit(np),alpha(n),beta(n),delta(n))
   call random_population(rng,x,lower,upper)
   evals=0
   call eval_population(fun,c,x,fit,evals)
   call sort_population(x,fit)
   alpha=x(1,:)
   fa=fit(1)
   beta=x(2,:)
   fb=fit(2)
   delta=x(3,:)
   fd=fit(3)
   do t=1,c%max_iter
      a=2.0_dp-2.0_dp*real(t,dp)/max(1.0_dp,real(c%max_iter,dp))
      do i=1,np
         do j=1,n
            r1=rng%uniform()
            r2=rng%uniform()
            a1=2*a*r1-a
            c1=2*r2
            x1=alpha(j)-a1*abs(c1*alpha(j)-x(i,j))
            r1=rng%uniform()
            r2=rng%uniform()
            a2=2*a*r1-a
            c2=2*r2
            x2=beta(j)-a2*abs(c2*beta(j)-x(i,j))
            r1=rng%uniform()
            r2=rng%uniform()
            a3=2*a*r1-a
            c3=2*r2
            x3=delta(j)-a3*abs(c3*delta(j)-x(i,j))
            x(i,j)=(x1+x2+x3)/3.0_dp
         end do
         call clamp_vector(x(i,:),lower,upper)
         f=internal_value(fun,x(i,:),c)
         evals=evals+1
         if(f<fa) then
            fa=f
            alpha=x(i,:)
         else if(f>fa .and. f<fb) then
            fb=f
            beta=x(i,:)
         else if(f>fa .and. f>fb .and. f<fd) then
            fd=f
            delta=x(i,:)
         end if
      end do
      result%history(t)=sign_value(c)*fa
   end do
   call finish_result(result,c,alpha,fa,c%max_iter,evals)
end subroutine gwo
subroutine sca(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),fit(:),best(:)
   real(dp)::fbest,r1,r2,r3,r4,f
   integer::n,np,i,j,t,evals
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   call init_result(result,n,c,'SCA')
   call rng%seed(c%seed)
   allocate(x(np,n),fit(np),best(n))
   call random_population(rng,x,lower,upper)
   evals=0
   call eval_population(fun,c,x,fit,evals)
   i=argmin(fit)
   best=x(i,:)
   fbest=fit(i)
   do t=1,c%max_iter
      r1=2.0_dp-2.0_dp*real(t,dp)/max(1.0_dp,real(c%max_iter,dp))
      do i=1,np
         do j=1,n
            r2=2.0_dp*mh_pi*rng%uniform()
            r3=2.0_dp*rng%uniform()
            r4=rng%uniform()
            if(r4<0.5_dp) then
               x(i,j)=x(i,j)+r1*sin(r2)*abs(r3*best(j)-x(i,j))
            else
               x(i,j)=x(i,j)+r1*cos(r2)*abs(r3*best(j)-x(i,j))
            end if
         end do
         call clamp_vector(x(i,:),lower,upper)
         f=internal_value(fun,x(i,:),c)
         evals=evals+1
         if(f<fbest) then
         fbest=f
         best=x(i,:)
         end if
      end do
      result%history(t)=sign_value(c)*fbest
   end do
   call finish_result(result,c,best,fbest,c%max_iter,evals)
end subroutine sca
subroutine woa(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),fit(:),best(:),xrand(:)
   real(dp)::fbest,a,a2,r1,r2,aa,cc,b,l,p,dist,f
   integer::n,np,i,j,t,k,evals
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   call init_result(result,n,c,'WOA')
   call rng%seed(c%seed)
   allocate(x(np,n),fit(np),best(n),xrand(n))
   call random_population(rng,x,lower,upper)
   evals=0
   call eval_population(fun,c,x,fit,evals)
   i=argmin(fit)
   best=x(i,:)
   fbest=fit(i)
   do t=1,c%max_iter
      a=2.0_dp-2.0_dp*real(t,dp)/real(max(1,c%max_iter),dp)
      a2=-1.0_dp-real(t,dp)/real(max(1,c%max_iter),dp)
      do i=1,np
         r1=rng%uniform()
         r2=rng%uniform()
         aa=2*a*r1-a
         cc=2*r2
         b=1.0_dp
         l=(a2-1.0_dp)*rng%uniform()+1.0_dp
         p=rng%uniform()
         if(p<0.5_dp .and. abs(aa)>=1.0_dp) then
            k=rng%randint(1,np)
            xrand=x(k,:)
            do j=1,n
            x(i,j)=xrand(j)-aa*abs(cc*xrand(j)-x(i,j))
            end do
         else if(p<0.5_dp) then
            do j=1,n
            x(i,j)=best(j)-aa*abs(cc*best(j)-x(i,j))
            end do
         else
            do j=1,n
            dist=abs(best(j)-x(i,j))
            x(i,j)=dist*exp(b*l)*cos(l*2*mh_pi)+best(j)
            end do
         end if
         call clamp_vector(x(i,:),lower,upper)
         f=internal_value(fun,x(i,:),c)
         evals=evals+1
         if(f<fbest) then
         fbest=f
         best=x(i,:)
         end if
      end do
      result%history(t)=sign_value(c)*fbest
   end do
   call finish_result(result,c,best,fbest,c%max_iter,evals)
end subroutine woa
subroutine mfo(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::moth(:,:),flames(:,:),fu(:,:),mf(:),ff(:),allfit(:),best(:)
   real(dp)::fbest,a,r,dist
   integer::n,np,i,j,t,nfl,evals
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   call init_result(result,n,c,'MFO')
   call rng%seed(c%seed)
   allocate(moth(np,n),flames(np,n),fu(2*np,n),mf(np),ff(np),allfit(2*np),best(n))
   call random_population(rng,moth,lower,upper)
   evals=0
   call eval_population(fun,c,moth,mf,evals)
   flames=moth
   ff=mf
   call sort_population(flames,ff)
   best=flames(1,:)
   fbest=ff(1)
   do t=1,c%max_iter
      nfl=nint(real(np,dp)-real(t,dp)*real(np-1,dp)/real(max(1,c%max_iter),dp))
      nfl=max(1,min(np,nfl))
      a=-1.0_dp-real(t,dp)/real(max(1,c%max_iter),dp)
      do i=1,np
         do j=1,n
            if(i<=nfl) then
            dist=abs(flames(i,j)-moth(i,j))
            else
            dist=abs(flames(nfl,j)-moth(i,j))
            end if
            r=(a-1.0_dp)*rng%uniform()+1.0_dp
            if(i<=nfl) then
               moth(i,j)=dist*exp(r)*cos(2*mh_pi*r)+flames(i,j)
            else
               moth(i,j)=dist*exp(r)*cos(2*mh_pi*r)+flames(nfl,j)
            end if
         end do
         call clamp_vector(moth(i,:),lower,upper)
      end do
      call eval_population(fun,c,moth,mf,evals)
      fu(1:np,:)=moth
      fu(np+1:2*np,:)=flames
      allfit(1:np)=mf
      allfit(np+1:2*np)=ff
      call sort_population(fu,allfit)
      flames=fu(1:np,:)
      ff=allfit(1:np)
      best=flames(1,:)
      fbest=ff(1)
      result%history(t)=sign_value(c)*fbest
   end do
   call finish_result(result,c,best,fbest,c%max_iter,evals)
end subroutine mfo
subroutine hs(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::hm(:,:),newh(:,:),allx(:,:),hf(:),nf(:),af(:),best(:)
   real(dp)::fbest,r,temp
   integer::n,np,i,j,k,t,evals
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   call init_result(result,n,c,'HS')
   call rng%seed(c%seed)
   allocate(hm(np,n),newh(np,n),allx(2*np,n),hf(np),nf(np),af(2*np),best(n))
   call random_population(rng,hm,lower,upper)
   evals=0
   call eval_population(fun,c,hm,hf,evals)
   i=argmin(hf)
   best=hm(i,:)
   fbest=hf(i)
   do t=1,c%max_iter
      do i=1,np
         do j=1,n
            if(rng%uniform()<c%hmcr) then
               k=rng%randint(1,np)
               newh(i,j)=hm(k,j)
               if(rng%uniform()<c%par) then
                  r=rng%uniform()
                  if(r<0.5_dp) then
                  temp=newh(i,j)-r*c%bandwidth
                  else
                  temp=newh(i,j)+r*c%bandwidth
                  end if
                  if(temp>=lower(j).and.temp<=upper(j)) newh(i,j)=temp
               end if
            else
               newh(i,j)=lower(j)+rng%uniform()*(upper(j)-lower(j))
            end if
         end do
      end do
      call eval_population(fun,c,newh,nf,evals)
      allx(1:np,:)=hm
      allx(np+1:2*np,:)=newh
      af(1:np)=hf
      af(np+1:2*np)=nf
      call sort_population(allx,af)
      hm=allx(1:np,:)
      hf=af(1:np)
      if(hf(1)<fbest) then
      fbest=hf(1)
      best=hm(1,:)
      end if
      result%history(t)=sign_value(c)*fbest
   end do
   call finish_result(result,c,best,fbest,c%max_iter,evals)
end subroutine hs
subroutine ffa(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),light(:),best(:)
   real(dp)::r,beta,randomization,fbest
   integer::n,np,i,j,k,t,evals,bi
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   call init_result(result,n,c,'FFA')
   call rng%seed(c%seed)
   allocate(x(np,n),light(np),best(n))
   call random_population(rng,x,lower,upper)
   evals=0
   call eval_population(fun,c,x,light,evals)
   call sort_population(x,light)
   best=x(1,:)
   fbest=light(1)
   do t=1,c%max_iter
      do i=1,np
         do j=1,np
            if(light(j)<light(i)) then
               r=sqrt(sum((x(i,:)-x(j,:))**2))
               beta=(1.0_dp-c%b0)*exp(-c%gamma*r*r)+c%b0
               do k=1,n
                  randomization=c%alpha_ffa*(rng%uniform()-0.5_dp)
                  x(i,k)=x(i,k)*(1.0_dp-beta)+x(j,k)*beta+randomization
               end do
               call clamp_vector(x(i,:),lower,upper)
            end if
         end do
      end do
      call eval_population(fun,c,x,light,evals)
      if(c%legacy_quirks) then
         bi=argmax(light) ! exact R source uses which.max here
      else
         bi=argmin(light)
      end if
      best=x(bi,:)
      fbest=light(bi)
      result%history(t)=sign_value(c)*fbest
   end do
   call finish_result(result,c,best,fbest,c%max_iter,evals)
end subroutine ffa
subroutine ga(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),fit(:),offs(:,:),newp(:,:),nf(:),best(:),weights(:)
   real(dp)::fbest
   integer::n,np,nsel,i,j,p1,p2,t,evals,k,idx,mcount
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   call init_result(result,n,c,'GA')
   call rng%seed(c%seed)
   allocate(x(np,n),fit(np),best(n),weights(np))
   call random_population(rng,x,lower,upper)
   evals=0
   call eval_population(fun,c,x,fit,evals)
   i=argmin(fit)
   best=x(i,:)
   fbest=fit(i)
   do t=1,c%max_iter
      nsel=nint(c%pc*real(np,dp))
      if(mod(nsel,2)/=0)nsel=nsel+1
      nsel=max(2,min(np,nsel))
      if(mod(nsel,2)/=0)nsel=nsel-1
      allocate(offs(nsel,n),newp(np+nsel,n),nf(np+nsel))
      ! R rouletteWhell is called on raw fitness values
      ! preserve that weighting.
      weights=fit
      do i=1,nsel
      idx=weighted_index(rng,weights)
      offs(i,:)=x(idx,:)
      end do
      do i=1,nsel,2
         p1=rng%randint(1,n)
         p2=rng%randint(1,n)
         do j=min(p1,p2),max(p1,p2)
            call swap_real(offs(i,j),offs(i+1,j))
         end do
      end do
      newp(1:np,:)=x
      newp(np+1:np+nsel,:)=offs
      mcount=max(1,ceiling(0.1_dp*real(n,dp)))
      do i=1,np+nsel
         if(rng%uniform()<c%pm) then
            do k=1,mcount
               j=rng%randint(1,n)
               newp(i,j)=lower(j)+rng%uniform()*(upper(j)-lower(j))
            end do
         end if
      end do
      call eval_population(fun,c,newp,nf,evals)
      call sort_population(newp,nf)
      x=newp(1:np,:)
      fit=nf(1:np)
      best=x(1,:)
      fbest=fit(1)
      result%history(t)=sign_value(c)*fbest
      deallocate(offs,newp,nf)
   end do
   call finish_result(result,c,best,fbest,c%max_iter,evals)
contains
   subroutine swap_real(a,b)
   real(dp),intent(inout)::a,b
   real(dp)::tmp
   tmp=a
   a=b
   b=tmp
   end subroutine
end subroutine ga
subroutine de(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),fit(:),mut(:,:),mfit(:),best(:)
   integer,allocatable::idx1(:),idx2(:),idx3(:),idx4(:),idx5(:)
   integer::n,np,i,j,t,evals,need
   character(len=:),allocatable::strategy

   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   strategy=trim(adjustl(c%de_strategy))
   select case(strategy)
   case('best 2','rand 2')
      need=5
   case('rand 2 dir')
      need=3
   case default
      need=2
   end select
   if(np<need) error stop 'DE: population too small for selected strategy'
   call init_result(result,n,c,'DE')
   call rng%seed(c%seed)
   allocate(x(np,n),fit(np),mut(np,n),mfit(np),best(n))
   allocate(idx1(np),idx2(np),idx3(np),idx4(np),idx5(np))
   call random_population(rng,x,lower,upper)
   evals=0
   call eval_population(fun,c,x,fit,evals)
   best=x(argmin(fit),:)

   do t=1,c%max_iter
      best=x(argmin(fit),:)
      call make_permutation(rng,idx1)
      call make_exclusive_permutation(rng,idx2,idx1)
      if(need>=3) call make_exclusive_permutation(rng,idx3,idx1,idx2)
      if(need>=4) call make_exclusive_permutation(rng,idx4,idx1,idx2,idx3)
      if(need>=5) call make_exclusive_permutation(rng,idx5,idx1,idx2,idx3,idx4)

      do i=1,np
         select case(strategy)
         case('clasical','classical')
            mut(i,:)=x(i,:)+c%scaling_vector*(x(idx1(i),:)-x(idx2(i),:))
         case('best 1')
            mut(i,:)=best+c%scaling_vector*(x(idx1(i),:)-x(idx2(i),:))
         case('target to best')
            mut(i,:)=x(i,:)+c%scaling_vector*(best-x(i,:))+ &
               c%scaling_vector*(x(idx1(i),:)-x(idx2(i),:))
         case('best 2')
            mut(i,:)=best+c%scaling_vector*(x(idx1(i),:)-x(idx2(i),:))+ &
               c%scaling_vector*(x(idx3(i),:)-x(idx4(i),:))
         case('rand 2')
            mut(i,:)=x(idx1(i),:)+c%scaling_vector*(x(idx2(i),:)-x(idx3(i),:))+ &
               c%scaling_vector*(x(idx4(i),:)-x(idx5(i),:))
         case('rand 2 dir')
            mut(i,:)=x(idx1(i),:)+0.5_dp*c%scaling_vector* &
               (x(idx1(i),:)-x(idx2(i),:)-x(idx3(i),:))
         case default
            error stop 'DE: unknown strategy'
         end select
         do j=1,n
            if(rng%uniform()>c%crossover_rate) mut(i,j)=x(i,j)
         end do
      end do

      call eval_population(fun,c,mut,mfit,evals)
      do i=1,np
         if(mfit(i)<fit(i)) then
            x(i,:)=mut(i,:)
            fit(i)=mfit(i)
         end if
      end do
      call clamp_population(x,lower,upper)
      ! The R source checks bounds after selection without refreshing cached fitness.
      if(.not.c%legacy_quirks) call eval_population(fun,c,x,fit,evals)
      best=x(argmin(fit),:)
      result%history(t)=sign_value(c)*minval(fit)
   end do
   i=argmin(fit)
   call finish_result(result,c,x(i,:),fit(i),c%max_iter,evals)
contains
   subroutine make_permutation(rng,p)
      type(mh_rng),intent(inout)::rng
      integer,intent(out)::p(:)
      integer::q
      do q=1,size(p)
         p(q)=q
      end do
      call shuffle_int(rng,p)
   end subroutine make_permutation

   subroutine make_exclusive_permutation(rng,p,a,b,d,e)
      type(mh_rng),intent(inout)::rng
      integer,intent(out)::p(:)
      integer,intent(in)::a(:)
      integer,intent(in),optional::b(:),d(:),e(:)
      logical::ok
      integer::tries
      tries=0
      do
         call make_permutation(rng,p)
         ok=all(p/=a)
         if(present(b))ok=ok.and.all(p/=b)
         if(present(d))ok=ok.and.all(p/=d)
         if(present(e))ok=ok.and.all(p/=e)
         if(ok)exit
         tries=tries+1
         if(tries>100000) error stop 'DE: could not construct exclusive permutations'
      end do
   end subroutine make_exclusive_permutation
end subroutine de
subroutine bho(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),fit(:),bh(:)
   real(dp)::fbh,eventh
   logical,allocatable::cross(:)
   integer::n,np,i,j,t,evals,bi
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   call init_result(result,n,c,'BHO')
   call rng%seed(c%seed)
   allocate(x(np,n),fit(np),bh(n),cross(np))
   call random_population(rng,x,lower,upper)
   evals=0
   do t=1,c%max_iter
      call eval_population(fun,c,x,fit,evals)
      bi=argmin(fit)
      bh=x(bi,:)
      do i=1,np
      do j=1,n
      x(i,j)=x(i,j)+rng%uniform()*(bh(j)-x(i,j))
      end do
      end do
      call eval_population(fun,c,x,fit,evals)
      bi=argmin(fit)
      bh=x(bi,:)
      fbh=fit(bi)
      if(abs(sum(fit))<=tiny(1.0_dp)) then
      eventh=0.0_dp
      else
      eventh=fbh/sum(fit)
      if(ieee_is_nan(eventh))eventh=0.0_dp
      end if
      cross=abs(fit-fbh)<eventh
      cross(bi)=.false.
      do i=1,np
         if(cross(i)) then
            do j=1,n
            x(i,j)=lower(j)+rng%uniform()*(upper(j)-lower(j))
            end do
         end if
      end do
      result%history(t)=sign_value(c)*fbh
   end do
   call eval_population(fun,c,x,fit,evals)
   bi=argmin(fit)
   call finish_result(result,c,x(bi,:),fit(bi),c%max_iter,evals)
end subroutine bho
subroutine gbs(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),v(:,:),fit(:),mass(:),force(:,:),best(:),worst(:),gbest(:)
   real(dp)::fbest,fworst,totalgm,gconst,dist,eps
   integer::n,np,i,j,k,t,evals,nk
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   call init_result(result,n,c,'GBS')
   call rng%seed(c%seed)
   allocate(x(np,n),v(np,n),fit(np),mass(np),force(np,n),best(n),worst(n),gbest(n))
   v=0
   call random_population(rng,x,lower,upper)
   evals=0
   call eval_population(fun,c,x,fit,evals)
   i=argmin(fit)
   gbest=x(i,:)
   fbest=fit(i)
   eps=8.854e-12_dp
   gconst=c%gravitational_const
   if(gconst<0)gconst=maxval(upper)
   do t=1,c%max_iter
      call eval_population(fun,c,x,fit,evals)
      call sort_population(x,fit)
      best=x(1,:)
      worst=x(np,:)
      fbest=fit(1)
      fworst=fit(np)
      if(abs(fbest-fworst)<=tiny(1.0_dp)) exit
      mass=(fit-fworst)/(fbest-fworst)
      totalgm=sum(mass)
      if(abs(totalgm)<=tiny(1.0_dp))exit
      mass=mass/totalgm
      gconst=gconst ! base control value is used through t formula below
      nk=max(1,min(np,nint(real(np,dp)*c%kbest)))
      force=0.0_dp
      do i=1,np
         do k=1,nk
            dist=sqrt(sum((x(i,:)-x(k,:))**2))
            if (dist > eps) then
               force(i,:) = force(i,:) + rng%uniform() * &
                  (merge(c%gravitational_const, maxval(upper), c%gravitational_const >= 0.0_dp) / &
                  exp(0.01_dp*t)) * mass(k) * (x(k,:) - x(i,:)) / dist
            end if
         end do
      end do
      do i=1,np
      do j=1,n
      v(i,j)=rng%uniform()*v(i,j)+force(i,j)
      x(i,j)=x(i,j)+v(i,j)
      end do
      end do
      call update_best(fun,c,x,best=gbest,fbest=fbest,evals=evals)
      result%history(t)=sign_value(c)*fbest
   end do
   call clamp_vector(gbest,lower,upper)
   fbest=internal_value(fun,gbest,c)
   evals=evals+1
   call finish_result(result,c,gbest,fbest,min(t,c%max_iter),evals)
end subroutine gbs
subroutine ba(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),v(:,:),freq(:,:),amp(:),pulse(:),fit(:),fly(:,:),ff(:),best(:)
   real(dp)::fbest
   integer::n,np,i,j,t,evals
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   call init_result(result,n,c,'BA')
   call rng%seed(c%seed)
   allocate(x(np,n),v(np,n),freq(np,n),amp(np),pulse(np),fit(np),fly(np,n),ff(np),best(n))
   v=0
   call random_population(rng,x,lower,upper)
   do i=1,np
   amp(i)=rng%uniform()
   pulse(i)=rng%uniform()
   end do
   evals=0
   call eval_population(fun,c,x,fit,evals)
   i=argmin(fit)
   best=x(i,:)
   fbest=fit(i)
   do t=1,c%max_iter
      do i=1,np
      do j=1,n
      freq(i,j)=c%min_frequency+rng%uniform()*(c%max_frequency-c%min_frequency)
      v(i,j)=v(i,j)+(x(i,j)-best(j))*freq(i,j)
      x(i,j)=x(i,j)+v(i,j)
      end do
      end do
      call clamp_population(x,lower,upper)
      call eval_population(fun,c,x,fit,evals)
      i=argmin(fit)
      if(fit(i)<fbest)then
      best=x(i,:)
      fbest=fit(i)
      end if
      do i=1,np
         if(rng%uniform()>pulse(i)) x(i,:)=amp(i)*best
      end do
      call random_population(rng,fly,lower,upper)
      call eval_population(fun,c,fly,ff,evals)
      do i=1,np
         if (ff(i) < fbest) then
            if (rng%uniform() < amp(i)) then
               x(i,:)=fly(i,:)
               pulse(i)=pulse(i)*(1.0_dp-exp(-c%gama))
               amp(i)=c%alpha_ba*amp(i)
               fbest=ff(i)
               best=fly(i,:)
            end if
         end if
      end do
      result%history(t)=sign_value(c)*fbest
   end do
   call clamp_vector(best,lower,upper)
   fbest=internal_value(fun,best,c)
   evals=evals+1
   call finish_result(result,c,best,fbest,c%max_iter,evals)
end subroutine ba
subroutine cs(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),fit(:),best(:),egg(:)
   real(dp)::fbest,fegg,sigma,beta,u,vv,step
   integer::n,np,i,j,t,evals,chosen,other,nremove
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   call init_result(result,n,c,'CS')
   call rng%seed(c%seed)
   allocate(x(np,n),fit(np),best(n),egg(n))
   call random_population(rng,x,lower,upper)
   evals=0
   call eval_population(fun,c,x,fit,evals)
   i=argmin(fit)
   best=x(i,:)
   fbest=fit(i)
   beta=1.5_dp
   sigma = (gamma(1.0_dp + beta) * sin(mh_pi*beta/2.0_dp) / &
      (gamma((1.0_dp + beta)/2.0_dp) * beta * &
      2.0_dp**((beta - 1.0_dp)/2.0_dp)))**(1.0_dp/beta)
   do t=1,c%max_iter
      i=argmin(fit)
      if(fit(i)<fbest)then
      fbest=fit(i)
      best=x(i,:)
      end if
      chosen=rng%randint(1,np)
      egg=x(chosen,:)
      do j=1,n
         u=rng%normal()*sigma
         vv=rng%normal()
         step=u/(abs(vv)+tiny(1.0_dp))**(1.0_dp/beta)
         egg(j)=egg(j)+0.01_dp*step*(x(chosen,j)-best(j))*rng%normal()
      end do
      other=rng%randint(1,np)
      fegg=internal_value(fun,egg,c)
      evals=evals+1
      if(fegg<fit(other)) then
      x(other,:)=egg
      fit(other)=fegg
      end if
      call sort_population(x,fit)
      nremove=nint(real(np,dp)*c%abandoned_fraction)
      nremove=max(0,min(np,nremove))
      if(nremove>0) then
         do i=np-nremove+1,np
            do j=1,n
            x(i,j)=lower(j)+rng%uniform()*(upper(j)-lower(j))
            end do
            fit(i)=internal_value(fun,x(i,:),c)
            evals=evals+1
         end do
      end if
      call clamp_population(x,lower,upper)
      i=argmin(fit)
      if(fit(i)<fbest)then
      fbest=fit(i)
      best=x(i,:)
      end if
      result%history(t)=sign_value(c)*fbest
   end do
   call finish_result(result,c,best,fbest,c%max_iter,evals)
end subroutine cs
subroutine abc(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),fit(:),newx(:),best(:),prob(:)
   integer,allocatable::lim(:)
   real(dp)::fnew,fbest,s
   integer::n,np,i,j,k,t,evals,cycle,bi
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   cycle=c%cycle_limit
   if(cycle<=0)cycle=n*np
   call init_result(result,n,c,'ABC')
   call rng%seed(c%seed)
   allocate(x(np,n),fit(np),newx(n),best(n),prob(np),lim(np))
   lim=cycle
   ! R generateRandomABC uses min(lower), max(upper) for all coordinates.
   do i=1,np
   do j=1,n
   x(i,j)=minval(lower)+rng%uniform()*(maxval(upper)-minval(lower))
   end do
   end do
   evals=0
   call eval_population(fun,c,x,fit,evals)
   bi=argmin(fit)
   best=x(bi,:)
   fbest=fit(bi)
   do t=1,c%max_iter
      ! Employed bees: R constructs a full random matrix but only one coordinate difference is nonzero.
      do i=1,np
         j=rng%randint(1,n)
         do
         k=rng%randint(1,np)
         if(k/=i)exit
         end do
         newx=x(i,:)
         newx(j)=x(i,j)+(2.0_dp*rng%uniform()-1.0_dp)*(x(i,j)-x(k,j))
         fnew=internal_value(fun,newx,c)
         evals=evals+1
         if(fnew<=fit(i)) then
         x(i,:)=newx
         fit(i)=fnew
         else
         lim(i)=lim(i)-1
         end if
      end do
      s=sum(fit)
      if(abs(s)<=tiny(1.0_dp))then
      prob=0.0_dp
      else
      prob=fit/s
      end if
      do i=1,np
         if(rng%uniform()<prob(i)) then
            j=rng%randint(1,n)
            do
            k=rng%randint(1,np)
            if(k/=i)exit
            end do
            newx=x(i,:)
            newx(j)=x(i,j)+(2.0_dp*rng%uniform()-1.0_dp)*(x(i,j)-x(k,j))
            fnew=internal_value(fun,newx,c)
            evals=evals+1
            if(fnew<=fit(i)) then
            x(i,:)=newx
            fit(i)=fnew
            else
            lim(i)=lim(i)-1
            end if
         end if
      end do
      bi=argmin(fit)
      lim(bi)=cycle
      do i=1,np
         if(lim(i)<=0) then
            do j=1,n
            x(i,j)=minval(lower)+rng%uniform()*(maxval(upper)-minval(lower))
            end do
            fit(i)=internal_value(fun,x(i,:),c)
            evals=evals+1
            lim(i)=cycle
         end if
      end do
      call clamp_population(x,lower,upper)
      if(.not.c%legacy_quirks) call eval_population(fun,c,x,fit,evals)
      bi=argmin(fit)
      if(fit(bi)<fbest)then
      fbest=fit(bi)
      best=x(bi,:)
      end if
      result%history(t)=sign_value(c)*fbest
   end do
   call finish_result(result,c,best,fbest,c%max_iter,evals)
end subroutine abc
subroutine alo(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::antlion(:,:),ant(:,:),fit(:),afit(:),best(:),ra(:),re(:),weights(:),allx(:,:),allf(:)
   real(dp)::fbest
   integer::n,np,i,t,k,evals
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   call init_result(result,n,c,'ALO')
   call rng%seed(c%seed)
   allocate(antlion(np,n),ant(np,n),fit(np),afit(np),best(n),ra(n),re(n),weights(np),allx(2*np,n),allf(2*np))
   call random_population(rng,antlion,lower,upper)
   call random_population(rng,ant,lower,upper)
   evals=0
   call eval_population(fun,c,antlion,fit,evals)
   call sort_population(antlion,fit)
   best=antlion(1,:)
   fbest=fit(1)
   do t=1,c%max_iter
      do i=1,np
         where(abs(fit)>tiny(1.0_dp))
         weights=1.0_dp/fit
         elsewhere
         weights=huge(1.0_dp)/real(np,dp)
         end where
         k=weighted_index(rng,weights)
         call alo_walk_at(rng,c%max_iter,lower,upper,antlion(k,:),t,ra)
         call alo_walk_at(rng,c%max_iter,lower,upper,best,t,re)
         ant(i,:)=(ra+re)/2.0_dp
      end do
      call clamp_population(ant,lower,upper)
      call eval_population(fun,c,ant,afit,evals)
      allx(1:np,:)=antlion
      allx(np+1:2*np,:)=ant
      allf(1:np)=fit
      allf(np+1:2*np)=afit
      call sort_population(allx,allf)
      antlion=allx(1:np,:)
      fit=allf(1:np)
      if(fit(1)<fbest) then
      fbest=fit(1)
      best=antlion(1,:)
      end if
      result%history(t)=sign_value(c)*fbest
   end do
   call finish_result(result,c,best,fbest,c%max_iter,evals)
end subroutine alo
subroutine alo_walk_at(rng,maxiter,lower,upper,pos,t,out)
   type(mh_rng),intent(inout)::rng
   integer,intent(in)::maxiter,t
   real(dp),intent(in)::lower(:),upper(:),pos(:)
   real(dp),intent(out)::out(:)
   real(dp),allocatable::lb(:),ub(:)
   real(dp)::ratio,cur,xmin,xmax,xt,step
   integer::j,s
   allocate(lb(size(pos)),ub(size(pos)))
   ratio=1.0_dp
   if(t>maxiter*0.1_dp)ratio=1.0_dp+100.0_dp*real(t,dp)/real(maxiter,dp)
   if(t>maxiter*0.5_dp)ratio=1.0_dp+1000.0_dp*real(t,dp)/real(maxiter,dp)
   if(t>maxiter*0.75_dp)ratio=1.0_dp+10000.0_dp*real(t,dp)/real(maxiter,dp)
   if(t>maxiter*0.9_dp)ratio=1.0_dp+100000.0_dp*real(t,dp)/real(maxiter,dp)
   if(t>maxiter*0.95_dp)ratio=1.0_dp+1000000.0_dp*real(t,dp)/real(maxiter,dp)
   lb=lower/ratio
   ub=upper/ratio
   if(rng%uniform()<0.5_dp)then
   lb=lb+pos
   else
   lb=-lb+pos
   end if
   if(rng%uniform()<0.5_dp)then
   ub=ub+pos
   else
   ub=-ub+pos
   end if
   do j=1,size(pos)
      cur=0
      xmin=0
      xmax=0
      xt=0
      do s=1,maxiter
         step=merge(1.0_dp,-1.0_dp,rng%uniform()>0.5_dp)
         cur=cur+step
         xmin=min(xmin,cur)
         xmax=max(xmax,cur)
         if(s==t-1)xt=cur
      end do
      if(t==1)xt=0.0_dp
      if(xmax>xmin)then
      out(j)=((xt-xmin)*(ub(j)-lb(j)))/(xmax-xmin)+lb(j)
      else
      out(j)=lb(j)
      end if
   end do
end subroutine alo_walk_at
subroutine da(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),delta(:,:),dmax(:),fit(:),food(:),enemy(:),rad(:)
   real(dp),allocatable::sep(:),align(:),cohes(:),fav(:),ene(:),best(:)
   real(dp)::ffood,fenemy,w,myc,sw,aw,cw,fw,ew,fval,beta,sigma,u,vv,step
   real(dp)::dist
   integer::n,np,i,j,k,t,evals,nn
   logical::neighbor,foodnear,enemynear
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   call init_result(result,n,c,'DA')
   call rng%seed(c%seed)
   allocate(x(np,n),delta(np,n),dmax(n),fit(np),food(n),enemy(n),rad(n),sep(n),align(n),cohes(n),fav(n),ene(n),best(n))
   call random_population(rng,x,lower,upper)
   dmax=(upper-lower)/20.0_dp
   do i=1,np
   do j=1,n
   delta(i,j)=-dmax(j)+2*dmax(j)*rng%uniform()
   end do
   end do
   evals=0
   call eval_population(fun,c,x,fit,evals)
   call sort_population(x,fit)
   food=x(1,:)
   ffood=fit(1)
   enemy=x(np,:)
   fenemy=fit(np)
   best=food
   beta=1.5_dp
   sigma = (gamma(1.0_dp + beta) * sin(mh_pi*beta/2.0_dp) / &
      (gamma((1.0_dp + beta)/2.0_dp) * beta * &
      2.0_dp**((beta - 1.0_dp)/2.0_dp)))**(1.0_dp/beta)
   do t=1,c%max_iter
      rad=(upper-lower)/4.0_dp+(upper-lower)*(real(t,dp)/real(max(1,c%max_iter),dp))*2.0_dp
      w=0.9_dp-real(t,dp)*(0.5_dp/real(max(1,c%max_iter),dp))
      myc=max(0.0_dp,0.1_dp-real(t,dp)*(0.1_dp/(real(max(1,c%max_iter),dp)/2.0_dp)))
      sw=2*rng%uniform()*myc
      aw=2*rng%uniform()*myc
      cw=2*rng%uniform()*myc
      fw=2*rng%uniform()
      ew=myc
      do i=1,np
         sep=0
         align=delta(i,:)
         cohes=0
         nn=0
         do k=1,np
            if(k==i)cycle
            neighbor=.true.
            do j=1,n
               dist=abs(x(i,j)-x(k,j))
               if(dist>rad(j).or.dist<=tiny(1.0_dp))neighbor=.false.
            end do
            if(neighbor)then
               nn=nn+1
               sep=sep+(x(k,:)-x(i,:))
               align=align+delta(k,:)
               cohes=cohes+x(k,:)
            end if
         end do
         if(nn>1)then
         sep=-sep
         align=(align-delta(i,:))/real(nn,dp)
         cohes=cohes/real(nn,dp)-x(i,:)
         else
         sep=0
         align=delta(i,:)
         cohes=0
         end if
         foodnear=all(abs(x(i,:)-food)<=rad)
         enemynear=all(abs(x(i,:)-enemy)<=rad)
         fav=0
         if(foodnear)fav=food-x(i,:)
         ene=0
         if(enemynear)ene=enemy+x(i,:)
         if(.not.foodnear)then
            if(nn>1)then
               do j=1,n
               delta(i,j)=w*delta(i,j)+rng%uniform()*align(j)+rng%uniform()*cohes(j)+rng%uniform()*sep(j)
               delta(i,j)=max(-dmax(j),min(dmax(j),delta(i,j)))
               x(i,j)=x(i,j)+delta(i,j)
               end do
            else
               do j=1,n
                  u=rng%uniform()*sigma
                  vv=max(rng%uniform(),tiny(1.0_dp))
                  step=u/abs(vv)**(1.0_dp/beta)
                  x(i,j)=x(i,j)+0.01_dp*step*x(i,j)
                  delta(i,j)=0
               end do
            end if
         else
            do j=1,n
            delta(i,j)=aw*align(j)+cw*cohes(j)+sw*sep(j)+fw*fav(j)+ew*ene(j)+w*delta(i,j)
            delta(i,j)=max(-dmax(j),min(dmax(j),delta(i,j)))
            x(i,j)=x(i,j)+delta(i,j)
            end do
         end if
         call clamp_vector(x(i,:),lower,upper)
      end do
      do i=1,np
         fval=internal_value(fun,x(i,:),c)
         evals=evals+1
         if(fval<ffood)then
         ffood=fval
         food=x(i,:)
         end if
         if(fval>fenemy)then
         fenemy=fval
         enemy=x(i,:)
         end if
      end do
      best=food
      result%history(t)=sign_value(c)*ffood
   end do
   call finish_result(result,c,best,ffood,c%max_iter,evals)
end subroutine da
subroutine goa(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::lo(:),hi(:),x(:,:),newx(:,:),fit(:),best(:),si(:)
   real(dp)::fbest,cc,r,xx,social,dx1,dx2,normr
   integer::n,n2,np,i,j,k,t,evals,bi
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   n2=n
   if(mod(n,2)/=0)n2=n+1
   np=c%num_population
   call init_result(result,n,c,'GOA')
   call rng%seed(c%seed)
   allocate(lo(n2),hi(n2),x(np,n2),newx(np,n2),fit(np),best(n2),si(n2))
   lo(1:n)=lower
   hi(1:n)=upper
   if(n2>n)then
   lo(n2)=sum(lower)/real(n,dp)
   hi(n2)=sum(upper)/real(n,dp)
   end if
   call random_population(rng,x,lo,hi) ! R reinitializes if odd dimension, same net effect.
   evals=0
   do i=1,np
   fit(i)=internal_value(fun,x(i,1:n),c)
   evals=evals+1
   end do
   bi=argmin(fit)
   best=x(bi,:)
   fbest=fit(bi)
   do t=1,c%max_iter
      cc=1.0_dp-real(t,dp)*(1.0_dp-0.00004_dp)/real(max(1,c%max_iter),dp)
      newx=x
      do i=1,np
         si=0.0_dp
         do j=1,n2,2
            do k=1,np
               if(k==i)cycle
               dx1=x(k,j)-x(i,j)
               dx2=x(k,j+1)-x(i,j+1)
               r=sqrt(dx1*dx1+dx2*dx2)
               normr=r+2.2204e-16_dp
               xx=2.0_dp+modulo(r,2.0_dp)
               social=0.5_dp*exp(-xx/1.5_dp)-exp(-xx)
               si(j)=si(j)+((hi(j)-lo(j))*cc/2.0_dp)*social*dx1/normr
               si(j+1)=si(j+1)+((hi(j+1)-lo(j+1))*cc/2.0_dp)*social*dx2/normr
            end do
         end do
         newx(i,:)=cc*si+best
      end do
      x=newx
      call clamp_population(x,lo,hi)
      do i=1,np
         fit(i)=internal_value(fun,x(i,1:n),c)
         evals=evals+1
         if(fit(i)<fbest)then
         fbest=fit(i)
         best=x(i,:)
         end if
      end do
      result%history(t)=sign_value(c)*fbest
   end do
   call finish_result(result,c,best(1:n),fbest,c%max_iter,evals)
end subroutine goa
subroutine clonalg(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),fit(:),clone(:,:),cf(:),allx(:,:),allf(:),best(:)
   integer::n,np,sel,i,j,t,evals,nc,total,pos,bi
   real(dp)::fbest
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   sel=c%selection_size
   if(sel<=0)sel=np/4
   sel=max(1,min(np,sel))
   call init_result(result,n,c,'CLONALG')
   call rng%seed(c%seed)
   allocate(x(np,n),fit(np),best(n))
   call random_population(rng,x,lower,upper)
   evals=0
   call eval_population(fun,c,x,fit,evals)
   bi=argmin(fit)
   best=x(bi,:)
   fbest=fit(bi)
   do t=1,c%max_iter
      call sort_population(x,fit)
      total=0
      do i=1,sel
      total=total+nint(c%multiplication_factor*real(np,dp)/real(i,dp))
      end do
      total=max(1,total)
      allocate(clone(total,n),cf(total))
      pos=0
      do i=1,sel
         nc=nint(c%multiplication_factor*real(np,dp)/real(i,dp))
         do j=1,nc
         pos=pos+1
         clone(pos,:)=x(i,:)
         end do
      end do
      do i=1,total
      do j=1,n
      if(rng%uniform()<=c%hypermutation_rate)clone(i,j)=lower(j)+rng%uniform()*(upper(j)-lower(j))
      end do
      end do
      call eval_population(fun,c,clone,cf,evals)
      allocate(allx(np+total,n),allf(np+total))
      allx(1:np,:)=x
      allx(np+1:,:)=clone
      allf(1:np)=fit
      allf(np+1:)=cf
      call sort_population(allx,allf)
      x=allx(1:np,:)
      fit=allf(1:np)
      deallocate(clone,cf,allx,allf)
      bi=argmin(fit)
      if(fit(bi)<fbest)then
      fbest=fit(bi)
      best=x(bi,:)
      end if
      result%history(t)=sign_value(c)*fbest
   end do
   if(c%legacy_quirks)then
   bi=argmax(fit)
   call finish_result(result,c,x(bi,:),fit(bi),c%max_iter,evals)
   else
   call finish_result(result,c,best,fbest,c%max_iter,evals)
   end if
end subroutine clonalg
subroutine sfl(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),fit(:),best(:),trial(:)
   integer,allocatable::mid(:),bestid(:),worstid(:)
   real(dp)::fbest,oldf,newf
   integer::n,np,nm,j,k,t,fl,evals,idx,ib,iw
   logical::moved
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   nm=c%num_memeplex
   if(nm<=0)nm=np/3
   nm=max(1,min(np,nm))
   call init_result(result,n,c,'SFL')
   call rng%seed(c%seed)
   allocate(x(np,n),fit(np),best(n),trial(n),mid(np),bestid(nm),worstid(nm))
   if (np < 1) error stop 'SFL: num_population must be positive'
   fit=huge(1.0_dp)
   call random_population(rng,x,lower,upper)
   evals=0
   call eval_population(fun,c,x,fit,evals)
   idx=argmin(fit)
   best=x(idx,:)
   fbest=fit(idx)
   call setup_memeplex_ids(np,nm,mid,bestid,worstid)
   do t=1,c%max_iter
      call reorder_sfl(fun,c,x,fit,mid,evals)
      do fl=1,c%frog_leaping_iteration
         do k=1,nm
            ib=bestid(k)
            iw=worstid(k)
            oldf=internal_value(fun,x(iw,:),c)
            evals=evals+1
            trial=x(iw,:)+random_vec(rng,n)*(x(ib,:)-x(iw,:))
            newf=internal_value(fun,trial,c)
            evals=evals+1
            moved=.false.
            if(newf<=oldf)then
            x(iw,:)=trial
            moved=.true.
            end if
            if(.not.moved)then
               trial=x(iw,:)+random_vec(rng,n)*(best-x(iw,:))
               newf=internal_value(fun,trial,c)
               evals=evals+1
               if(newf<=oldf)then
               x(iw,:)=trial
               moved=.true.
               end if
            end if
            if(.not.moved)then
            do j=1,n
            x(iw,j)=lower(j)+rng%uniform()*(upper(j)-lower(j))
            end do
            end if
         end do
         call reorder_sfl(fun,c,x,fit,mid,evals)
         idx=argmin(fit)
         if(fit(idx)<fbest)then
         fbest=fit(idx)
         best=x(idx,:)
         end if
      end do
      result%history(t)=sign_value(c)*fbest
   end do
   call finish_result(result,c,best,fbest,c%max_iter,evals)
contains
   function random_vec(rng,n) result(v)
      type(mh_rng),intent(inout)::rng
      integer,intent(in)::n
      real(dp)::v(n)
      integer::q
      do q=1,n
      v(q)=rng%uniform()
      end do
   end function
   subroutine setup_memeplex_ids(np,nm,mid,bestid,worstid)
      integer,intent(in)::np,nm
      integer,intent(out)::mid(np),bestid(nm),worstid(nm)
      integer::q,k
      do q=1,np
      mid(q)=mod(q-1,nm)+1
      end do
      call sort_int(mid)
      do k=1,nm
      bestid(k)=0
      worstid(k)=0
      do q=1,np
      if(mid(q)==k)then
      if(bestid(k)==0)bestid(k)=q
      worstid(k)=q
      end if
      end do
      end do
   end subroutine
   subroutine sort_int(a)
   integer,intent(inout)::a(:)
   integer::q,r,tmp
   do q=2,size(a)
   tmp=a(q)
   r=q-1
   do while(r>=1)
   if(a(r)<=tmp)exit
   a(r+1)=a(r)
   r=r-1
   end do
   a(r+1)=tmp
   end do
   end subroutine
   subroutine reorder_sfl(fun,c,x,fit,mid,evals)
      procedure(objective_function)::fun
      type(mh_control),intent(in)::c
      real(dp),intent(inout)::x(:,:),fit(:)
      integer,intent(inout)::mid(:)
      integer,intent(inout)::evals
      integer::q,np,nm
      real(dp),allocatable::key(:),tmpx(:,:)
      integer,allocatable::tmpm(:)
      np=size(x,1)
      nm=maxval(mid)
      allocate(key(np))
      do q=1,np
         if(c%legacy_quirks)then
         key(q)=fun(x(q,:))
         else
         key(q)=internal_value(fun,x(q,:),c)
         end if
         evals=evals+1
      end do
      call sort_population(x,key)
      do q=1,np
      mid(q)=mod(q-1,nm)+1
      end do
      allocate(tmpx(np,size(x,2)),tmpm(np))
      call group_by_mid(x,mid,tmpx,tmpm)
      x=tmpx
      mid=tmpm
      call eval_population(fun,c,x,fit,evals)
   end subroutine
   subroutine group_by_mid(x,mid,y,my)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::mid(:)
      real(dp),intent(out)::y(:,:)
      integer,intent(out)::my(:)
      integer::k,q,p
      p=0
      do k=1,maxval(mid)
      do q=1,size(mid)
      if(mid(q)==k)then
      p=p+1
      y(p,:)=x(q,:)
      my(p)=k
      end if
      end do
      end do
   end subroutine
end subroutine sfl
subroutine cso(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),vel(:,:),fit(:),best(:),worst(:),copies(:,:),pfit(:),prob(:)
   integer,allocatable::flag(:),pick(:)
   real(dp)::fbest,fworst,den,r
   integer::n,np,i,j,k,t,evals,nseek,ntrace,ncopy,chosen,bi,wi,nmod,q
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   if(c%cdc<=0)c%cdc=n
   call init_result(result,n,c,'CSO')
   call rng%seed(c%seed)
   allocate(x(np,n),vel(np,n),fit(np),best(n),worst(n),flag(np))
   call random_population(rng,x,lower,upper)
   do i=1,np
   do j=1,n
   vel(i,j)=rng%uniform()*c%maximum_velocity
   end do
   end do
   evals=0
   call eval_population(fun,c,x,fit,evals)
   bi=argmin(fit)
   best=x(bi,:)
   fbest=fit(bi)
   do t=1,c%max_iter
      ! R flagingCSO creates ceil(mixture*np) seeking and ceil((1-mixture)*np) tracing,
      ! shuffles, then truncates back to np.
      nseek=ceiling(c%mixture_ratio*real(np,dp))
      ntrace=ceiling((1.0_dp-c%mixture_ratio)*real(np,dp))
      call build_flags(rng,np,nseek,ntrace,flag)
      bi=argmin(fit)
      wi=argmax(fit)
      if(fit(bi)<fbest)then
      fbest=fit(bi)
      best=x(bi,:)
      end if
      worst=x(wi,:)
      fworst=fit(wi)
      ! tracing mode
      do i=1,np
         if(flag(i)==2)then
            do j=1,n
               vel(i,j)=vel(i,j)+rng%uniform()*c%tracing_constant*(best(j)-x(i,j))
               if(vel(i,j)>c%maximum_velocity)vel(i,j)=c%maximum_velocity
               x(i,j)=x(i,j)+vel(i,j)
            end do
         end if
      end do
      ! seeking mode, independently for each cat
      do i=1,np
         if(flag(i)/=1)cycle
         if(c%spc)then
         ncopy=c%smp+1
         else
         ncopy=max(1,c%smp)
         end if
         allocate(copies(ncopy,n),pfit(ncopy),prob(ncopy))
         copies(1,:)=x(i,:)
         do k=2,ncopy
         copies(k,:)=x(i,:)
         end do
         do k=2,ncopy
            nmod=max(0,min(n,c%cdc))
            if(nmod>0)then
               allocate(pick(nmod))
               call sample_distinct(rng,n,nmod,pick)
               do q=1,nmod
                  j=pick(q)
                  if(rng%uniform()<0.5_dp)then
                  r=1.0_dp
                  else
                  r=-1.0_dp
                  end if
                  copies(k,j)=copies(k,j)*r*c%srd/100.0_dp
               end do
               deallocate(pick)
            end if
         end do
         call eval_population(fun,c,copies,pfit,evals)
         den=fworst-fbest
         if(abs(den)<=tiny(1.0_dp))then
         prob=0.0_dp
         else
         prob=abs(pfit-fbest)/den
         end if
         if(all(abs(prob)<=tiny(1.0_dp)))then
         chosen=rng%randint(1,ncopy)
         else
         chosen=weighted_index(rng,prob)
         end if
         x(i,:)=copies(chosen,:)
         deallocate(copies,pfit,prob)
      end do
      call eval_population(fun,c,x,fit,evals)
      call clamp_population(x,lower,upper)
      call eval_population(fun,c,x,fit,evals)
      bi=argmin(fit)
      if(.not.c%legacy_quirks .and. fit(bi)<fbest)then
      fbest=fit(bi)
      best=x(bi,:)
      end if
      result%history(t)=sign_value(c)*fbest
   end do
   call finish_result(result,c,best,fbest,c%max_iter,evals)
contains
   subroutine build_flags(rng,np,ns,nt,flag)
      type(mh_rng),intent(inout)::rng
      integer,intent(in)::np,ns,nt
      integer,intent(out)::flag(np)
      integer,allocatable::tmp(:)
      integer::m,q
      m=max(np,ns+nt)
      allocate(tmp(m))
      tmp=2
      do q=1,min(ns,m)
      tmp(q)=1
      end do
      do q=ns+1,min(ns+nt,m)
      tmp(q)=2
      end do
      call shuffle_int(rng,tmp)
      flag=tmp(1:np)
   end subroutine
end subroutine cso
subroutine kh(fun,lower,upper,result,control)
   procedure(objective_function)::fun
   real(dp),intent(in)::lower(:),upper(:)
   type(mh_result),intent(out)::result
   type(mh_control),intent(in),optional::control
   type(mh_control)::c
   type(mh_rng)::rng
   real(dp),allocatable::x(:,:),fit(:),nmotion(:,:),forage(:,:),alpha(:,:),beta_v(:,:)
   real(dp),allocatable::best(:),worst(:),gbest(:),xfood(:),dir(:),distmat(:,:),sense(:)
   real(dp)::fbest,fworst,fgbest,foodfit,klocal,kbestv,c_best,c_food,diff,delta_t
   real(dp)::cr,mu_prob,wsum
   integer::n,np,i,j,k,t,evals,p,q,bi,wi
   c=mh_control()
   if(present(control))c=control
   call validate_bounds(lower,upper)
   n=size(lower)
   np=c%num_population
   call init_result(result,n,c,'KH')
   call rng%seed(c%seed)
   allocate(x(np,n),fit(np),nmotion(np,n),forage(np,n),alpha(np,n),beta_v(np,n))
   allocate(best(n),worst(n),gbest(n),xfood(n),dir(n),distmat(np,np),sense(np))
   nmotion=0.0_dp
   forage=0.0_dp
   call random_population(rng,x,lower,upper)
   evals=0
   call eval_population(fun,c,x,fit,evals)
   bi=argmin(fit)
   gbest=x(bi,:)
   fgbest=fit(bi)
   do t=1,c%max_iter
      call eval_population(fun,c,x,fit,evals)
      bi=argmin(fit)
      wi=argmax(fit)
      best=x(bi,:)
      worst=x(wi,:)
      fbest=fit(bi)
      fworst=fit(wi)
      if(fbest<fgbest)then
      fgbest=fbest
      gbest=best
      end if
      do i=1,np
      do k=1,np
      distmat(i,k)=sqrt(sum((x(i,:)-x(k,:))**2))
      end do
      sense(i)=real(n,dp)*sum(distmat(i,:))/5.0_dp
      end do
      alpha=0.0_dp
      do i=1,np
         do k=1,np
            if(distmat(i,k)<sense(i))then
               dir=(x(k,:)-x(i,:))/(distmat(i,k)+c%epsilon)
               if(abs(fworst-fbest)<=tiny(1.0_dp))then
               klocal=0.0_dp
               else
               klocal=(fit(i)-fit(k))/(fworst-fbest)
               end if
               alpha(i,:)=alpha(i,:)+klocal*dir
            end if
         end do
         c_best=2.0_dp*(rng%uniform()+real(t,dp)/real(max(1,c%max_iter),dp))
         dir=(best-x(i,:))/(sqrt(sum((best-x(i,:))**2))+c%epsilon)
         if(abs(fworst-fbest)<=tiny(1.0_dp))then
         kbestv=0.0_dp
         else
         kbestv=(fit(i)-fbest)/(fworst-fbest)
         end if
         alpha(i,:)=alpha(i,:)+c_best*kbestv*dir
      end do
      nmotion=c%max_motion_induced*alpha+c%inertia_motion*nmotion
      ! food location weighted by reciprocal fitness, as in the R source
      wsum=0.0_dp
      xfood=0.0_dp
      do i=1,np
         if(abs(fit(i))>tiny(1.0_dp))then
         wsum=wsum+1.0_dp/fit(i)
         xfood=xfood+x(i,:)/fit(i)
         end if
      end do
      if(abs(wsum)>tiny(1.0_dp))then
      xfood=xfood/wsum
      else
      xfood=0.0_dp
      end if
      foodfit=internal_value(fun,xfood,c)
      evals=evals+1
      c_food=2.0_dp*(1.0_dp-real(t,dp)/real(max(1,c%max_iter),dp))
      do i=1,np
         if(abs(fworst-fbest)<=tiny(1.0_dp))then
            beta_v(i,:)=0.0_dp
         else
            klocal=(fit(i)-foodfit)/(fworst-fbest)
            dir=(xfood-x(i,:))/(sqrt(sum((xfood-x(i,:))**2))+c%epsilon)
            beta_v(i,:)=c_food*klocal*dir
            kbestv=(fit(i)-fbest)/(fworst-fbest)
            dir=(best-x(i,:))/(sqrt(sum((best-x(i,:))**2))+c%epsilon)
            beta_v(i,:)=beta_v(i,:)+kbestv*dir
         end if
      end do
      forage=c%foraging_speed*beta_v+c%inertia_foraging*forage
      diff=c%max_diffusion_speed*(1.0_dp-real(t,dp)/real(max(1,c%max_iter),dp))*(2.0_dp*rng%uniform()-1.0_dp)
      delta_t=c%constant_space*sum(upper-lower)
      x=x+delta_t*(nmotion+forage+diff)
      ! Genetic crossover
      call eval_population(fun,c,x,fit,evals)
      bi=argmin(fit)
      wi=argmax(fit)
      fbest=fit(bi)
      fworst=fit(wi)
      if(fbest<fgbest)then
      fgbest=fbest
      gbest=x(bi,:)
      end if
      do i=1,np
         if(abs(fworst-fbest)<=tiny(1.0_dp))then
         kbestv=0.0_dp
         else
         kbestv=(fit(i)-fbest)/(fworst-fbest)
         end if
         cr=0.2_dp*kbestv
         do j=1,n
            if(rng%uniform()<cr)then
            k=rng%randint(1,np)
            x(i,j)=x(k,j)
            end if
         end do
      end do
      ! Genetic mutation
      call eval_population(fun,c,x,fit,evals)
      bi=argmin(fit)
      wi=argmax(fit)
      fbest=fit(bi)
      fworst=fit(wi)
      if(fbest<fgbest)then
      fgbest=fbest
      gbest=x(bi,:)
      end if
      do i=1,np
         if(abs(fworst-fbest)<=tiny(1.0_dp))then
         kbestv=0.0_dp
         else
         kbestv=(fit(i)-fbest)/(fworst-fbest)
         end if
         mu_prob=0.05_dp*kbestv
         do j=1,n
            if(rng%uniform()<mu_prob)then
            p=rng%randint(1,np)
            q=rng%randint(1,np)
            x(i,j)=gbest(j)+c%mu*(x(p,j)-x(q,j))
            end if
         end do
      end do
      result%history(t)=sign_value(c)*fgbest
   end do
   call eval_population(fun,c,x,fit,evals)
   bi=argmin(fit)
   if(fit(bi)<fgbest)then
   fgbest=fit(bi)
   gbest=x(bi,:)
   end if
   call clamp_vector(gbest,lower,upper)
   fgbest=internal_value(fun,gbest,c)
   evals=evals+1
   call finish_result(result,c,gbest,fgbest,c%max_iter,evals)
end subroutine kh
end module metaheuristic_opt
