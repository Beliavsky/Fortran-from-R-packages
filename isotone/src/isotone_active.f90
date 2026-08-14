module isotone_active
   use isotone_kinds, only : dp
   use isotone_utils, only : weighted_mean_value, weighted_median_value, &
      weighted_fractile_value, weighted_midrange_value, constraint_values, &
      component_labels, lagrange_multipliers, transpose_constraint_product, kkt_values
   use isotone_linalg, only : solve_linear
   implicit none
   private

   integer, parameter, public :: ISO_LS = 1
   integer, parameter, public :: ISO_L1 = 2
   integer, parameter, public :: ISO_QUANTILE = 3
   integer, parameter, public :: ISO_GLS = 4
   integer, parameter, public :: ISO_POISSON = 5
   integer, parameter, public :: ISO_LP = 6
   integer, parameter, public :: ISO_ASYLS = 7
   integer, parameter, public :: ISO_L1EPS = 8
   integer, parameter, public :: ISO_HUBER = 9
   integer, parameter, public :: ISO_SILF = 10
   integer, parameter, public :: ISO_CHEBYSHEV = 11
   integer, parameter, public :: ISO_CUSTOM = 12

   integer, parameter, public :: ISO_SUCCESS = 0
   integer, parameter, public :: ISO_MAXITER = 1
   integer, parameter, public :: ISO_BAD_INPUT = 2
   integer, parameter, public :: ISO_SOLVER_FAILURE = 3

   type, public :: active_set_options
      integer :: solver = ISO_LS
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: weight_matrix(:,:)
      real(dp) :: aw = 1.0_dp
      real(dp) :: bw = 1.0_dp
      real(dp) :: eps = 1.0e-4_dp
      real(dp) :: beta = 0.5_dp
      real(dp) :: p = 2.0_dp
      integer :: bfgs_maxiter = 300
      real(dp) :: bfgs_tol = 1.0e-10_dp
   end type active_set_options


   type, public :: isotone_solver_result
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: lambda(:)
      real(dp), allocatable :: gradient(:)
      real(dp) :: f = 0.0_dp
      integer :: status = ISO_BAD_INPUT
   end type isotone_solver_result

   type, public :: active_set_result
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: lambda(:)
      real(dp), allocatable :: constr_val(:)
      real(dp), allocatable :: alambda(:)
      real(dp), allocatable :: gradient(:)
      real(dp) :: fval = 0.0_dp
      real(dp) :: kkt(4) = 0.0_dp
      integer :: niter = 0
      integer :: status = ISO_BAD_INPUT
   end type active_set_result

   abstract interface
      subroutine loss_callback(x, f, g)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: f
         real(dp), intent(out) :: g(:)
      end subroutine loss_callback
   end interface

   public :: active_set, active_set_custom
   public :: ls_solver, d_solver, p_solver, lf_solver, s_solver, o_solver
   public :: a_solver, e_solver, h_solver, i_solver, m_solver, f_solver
contains
   subroutine ls_solver(z,a,y,weights,result)
      real(dp),intent(in)::z(:),y(:),weights(:)
      integer,intent(in)::a(:,:)
      type(isotone_solver_result),intent(out)::result
      type(active_set_options)::opt
      allocate(opt%weights(size(weights)));opt%weights=weights;opt%solver=ISO_LS
      call standalone_solver(z,a,y,opt,result)
   end subroutine ls_solver

   subroutine d_solver(z,a,y,weights,result)
      real(dp),intent(in)::z(:),y(:),weights(:)
      integer,intent(in)::a(:,:)
      type(isotone_solver_result),intent(out)::result
      type(active_set_options)::opt
      allocate(opt%weights(size(weights)));opt%weights=weights;opt%solver=ISO_L1
      call standalone_solver(z,a,y,opt,result)
   end subroutine d_solver

   subroutine p_solver(z,a,y,weights,aw,bw,result)
      real(dp),intent(in)::z(:),y(:),weights(:),aw,bw
      integer,intent(in)::a(:,:)
      type(isotone_solver_result),intent(out)::result
      type(active_set_options)::opt
      allocate(opt%weights(size(weights)));opt%weights=weights;opt%solver=ISO_QUANTILE
      opt%aw=aw;opt%bw=bw;call standalone_solver(z,a,y,opt,result)
   end subroutine p_solver

   subroutine lf_solver(z,a,y,weight_matrix,result)
      real(dp),intent(in)::z(:),y(:),weight_matrix(:,:)
      integer,intent(in)::a(:,:)
      type(isotone_solver_result),intent(out)::result
      type(active_set_options)::opt
      allocate(opt%weight_matrix(size(weight_matrix,1),size(weight_matrix,2)))
      opt%weight_matrix=weight_matrix;opt%solver=ISO_GLS
      call standalone_solver(z,a,y,opt,result)
   end subroutine lf_solver

   subroutine s_solver(z,a,y,result)
      real(dp),intent(in)::z(:),y(:)
      integer,intent(in)::a(:,:)
      type(isotone_solver_result),intent(out)::result
      type(active_set_options)::opt
      opt%solver=ISO_POISSON;call standalone_solver(z,a,y,opt,result)
   end subroutine s_solver

   subroutine o_solver(z,a,y,weights,pow,result)
      real(dp),intent(in)::z(:),y(:),weights(:),pow
      integer,intent(in)::a(:,:)
      type(isotone_solver_result),intent(out)::result
      type(active_set_options)::opt
      allocate(opt%weights(size(weights)));opt%weights=weights;opt%solver=ISO_LP;opt%p=pow
      call standalone_solver(z,a,y,opt,result)
   end subroutine o_solver

   subroutine a_solver(z,a,y,weights,aw,bw,result)
      real(dp),intent(in)::z(:),y(:),weights(:),aw,bw
      integer,intent(in)::a(:,:)
      type(isotone_solver_result),intent(out)::result
      type(active_set_options)::opt
      allocate(opt%weights(size(weights)));opt%weights=weights;opt%solver=ISO_ASYLS
      opt%aw=aw;opt%bw=bw;call standalone_solver(z,a,y,opt,result)
   end subroutine a_solver

   subroutine e_solver(z,a,y,weights,eps,result)
      real(dp),intent(in)::z(:),y(:),weights(:),eps
      integer,intent(in)::a(:,:)
      type(isotone_solver_result),intent(out)::result
      type(active_set_options)::opt
      allocate(opt%weights(size(weights)));opt%weights=weights;opt%solver=ISO_L1EPS;opt%eps=eps
      call standalone_solver(z,a,y,opt,result)
   end subroutine e_solver

   subroutine h_solver(z,a,y,weights,eps,result)
      real(dp),intent(in)::z(:),y(:),weights(:),eps
      integer,intent(in)::a(:,:)
      type(isotone_solver_result),intent(out)::result
      type(active_set_options)::opt
      allocate(opt%weights(size(weights)));opt%weights=weights;opt%solver=ISO_HUBER;opt%eps=eps
      call standalone_solver(z,a,y,opt,result)
   end subroutine h_solver

   subroutine i_solver(z,a,y,weights,eps,beta,result)
      real(dp),intent(in)::z(:),y(:),weights(:),eps,beta
      integer,intent(in)::a(:,:)
      type(isotone_solver_result),intent(out)::result
      type(active_set_options)::opt
      allocate(opt%weights(size(weights)));opt%weights=weights;opt%solver=ISO_SILF
      opt%eps=eps;opt%beta=beta;call standalone_solver(z,a,y,opt,result)
   end subroutine i_solver

   subroutine m_solver(z,a,y,weights,result)
      real(dp),intent(in)::z(:),y(:),weights(:)
      integer,intent(in)::a(:,:)
      type(isotone_solver_result),intent(out)::result
      type(active_set_options)::opt
      allocate(opt%weights(size(weights)));opt%weights=weights;opt%solver=ISO_CHEBYSHEV
      call standalone_solver(z,a,y,opt,result)
   end subroutine m_solver

   subroutine f_solver(z,a,y,callback,result)
      real(dp),intent(in)::z(:),y(:)
      integer,intent(in)::a(:,:)
      procedure(loss_callback)::callback
      type(isotone_solver_result),intent(out)::result
      type(active_set_options)::opt
      opt%solver=ISO_CUSTOM;call standalone_solver(z,a,y,opt,result,callback)
   end subroutine f_solver

   subroutine standalone_solver(z,a,y,opt,result,callback)
      real(dp),intent(in)::z(:),y(:)
      integer,intent(in)::a(:,:)
      type(active_set_options),intent(in)::opt
      type(isotone_solver_result),intent(out)::result
      procedure(loss_callback),optional::callback
      integer,allocatable::active(:)
      real(dp),allocatable::x(:),lbd(:),g(:)
      real(dp)::f
      logical::ok
      integer::i
      if(size(a,2)/=2 .or. size(z)/=size(y)) then;result%status=ISO_BAD_INPUT;return;end if
      allocate(active(size(a,1)),x(size(y)),g(size(y)))
      active=[(i,i=1,size(a,1))]
      if(present(callback)) then
         call solve_active(y,a,active,z,opt,x,lbd,f,g,ok,callback)
      else
         call solve_active(y,a,active,z,opt,x,lbd,f,g,ok)
      end if
      allocate(result%x(size(x)),result%gradient(size(g)),result%lambda(size(lbd)))
      result%x=x;result%gradient=g;result%lambda=lbd;result%f=f
      if(ok) then;result%status=ISO_SUCCESS;else;result%status=ISO_SOLVER_FAILURE;end if
   end subroutine standalone_solver

   subroutine active_set(isomat, yobs, result, options, x0, ups, maxiter, check)
      integer, intent(in) :: isomat(:,:)
      real(dp), intent(in) :: yobs(:)
      type(active_set_result), intent(out) :: result
      type(active_set_options), intent(in), optional :: options
      real(dp), intent(in), optional :: x0(:), ups
      integer, intent(in), optional :: maxiter
      logical, intent(in), optional :: check
      type(active_set_options) :: opt
      real(dp) :: tol
      integer :: itmax
      logical :: do_check

      opt = active_set_options()
      if (present(options)) opt = options
      tol = 1.0e-12_dp; if (present(ups)) tol = ups
      itmax = 100; if (present(maxiter)) itmax = maxiter
      do_check = .true.; if (present(check)) do_check = check
      call active_set_driver(isomat,yobs,result,opt,x0,tol,itmax,do_check)
   end subroutine active_set

   subroutine active_set_custom(isomat, yobs, callback, result, x0, ups, maxiter, check)
      integer, intent(in) :: isomat(:,:)
      real(dp), intent(in) :: yobs(:)
      procedure(loss_callback) :: callback
      type(active_set_result), intent(out) :: result
      real(dp), intent(in), optional :: x0(:), ups
      integer, intent(in), optional :: maxiter
      logical, intent(in), optional :: check
      type(active_set_options) :: opt
      real(dp) :: tol
      integer :: itmax
      logical :: do_check
      opt%solver = ISO_CUSTOM
      tol=1.0e-12_dp; if(present(ups)) tol=ups
      itmax=100; if(present(maxiter)) itmax=maxiter
      do_check=.true.; if(present(check)) do_check=check
      call active_set_driver(isomat,yobs,result,opt,x0,tol,itmax,do_check,callback)
   end subroutine active_set_custom

   subroutine active_set_driver(isomat,yobs,result,opt,x0,ups,maxiter,do_check,callback)
      integer,intent(in)::isomat(:,:),maxiter
      real(dp),intent(in)::yobs(:),ups
      type(active_set_result),intent(out)::result
      type(active_set_options),intent(in)::opt
      real(dp),intent(in),optional::x0(:)
      logical,intent(in)::do_check
      procedure(loss_callback),optional::callback
      real(dp),allocatable::xold(:),xnew(:),candidate(:),ax(:),ay(:),grad(:), &
         lbd(:),full_lbd(:),alambda(:),rat(:)
      integer,allocatable::active(:),cand_idx(:)
      integer::n,m,iter,iy,il,j,kidx,ncand
      real(dp)::my,ml,alpha,denom
      logical::ok,feasible,dual_ok

      n=size(yobs);m=size(isomat,1)
      result%status=ISO_BAD_INPUT
      if(size(isomat,2)/=2 .or. n<=0) return
      if(any(isomat<1) .or. any(isomat>n)) return
      allocate(xold(n),xnew(n),candidate(n),ax(m),ay(m),grad(n))
      if(present(x0)) then
         if(size(x0)/=n) return
         xold=x0
      else
         xold=0.0_dp
      end if
      xnew=xold
      call constraint_values(isomat,xold,ax)
      if(any(ax < -ups)) return
      call active_indices(ax,ups,active)

      do iter=1,maxiter
         if(present(callback)) then
            call solve_active(yobs,isomat,active,xold,opt,candidate,lbd, &
               result%fval,grad,ok,callback)
         else
            call solve_active(yobs,isomat,active,xold,opt,candidate,lbd, &
               result%fval,grad,ok)
         end if
         if(.not.ok) then
            result%status=ISO_SOLVER_FAILURE
            exit
         end if
         call constraint_values(isomat,candidate,ay)
         iy=minloc(ay,dim=1);my=ay(iy)
         if(size(lbd)==0) then
            ml=huge(1.0_dp);il=0
         else
            il=minloc(lbd,dim=1);ml=lbd(il)
         end if
         feasible = my > -ups
         dual_ok = ml > -ups
         if(feasible .and. dual_ok) then
            xnew=candidate
            result%status=ISO_SUCCESS
            exit
         else if(feasible) then
            xnew=candidate;ax=ay
            if(il>0) call remove_position(active,il)
         else
            allocate(cand_idx(m),rat(m));ncand=0
            do j=1,m
               if(ax(j)>-ups .and. ay(j)<ups) then
                  denom=ax(j)-ay(j)
                  if(denom>tiny(1.0_dp)) then
                     ncand=ncand+1;cand_idx(ncand)=j;rat(ncand)=-ay(j)/denom
                  end if
               end if
            end do
            if(ncand==0) then
               result%status=ISO_SOLVER_FAILURE
               deallocate(cand_idx,rat)
               exit
            end if
            kidx=maxloc(rat(1:ncand),dim=1);alpha=rat(kidx)
            xnew=candidate+alpha*(xold-candidate)
            call constraint_values(isomat,xnew,ax)
            call add_unique_sorted(active,cand_idx(kidx))
            deallocate(cand_idx,rat)
         end if
         xold=xnew
         if(iter==maxiter) result%status=ISO_MAXITER
      end do

      result%niter=min(iter,maxiter)
      allocate(result%x(n),result%y(n),result%lambda(m),result%constr_val(m), &
         result%alambda(n),result%gradient(n))
      result%x=xnew;result%y=yobs
      call constraint_values(isomat,result%x,result%constr_val)
      ! Recompute final solver information at the final active set.
      if(present(callback)) then
         call solve_active(yobs,isomat,active,result%x,opt,candidate,lbd, &
            result%fval,grad,ok,callback)
      else
         call solve_active(yobs,isomat,active,result%x,opt,candidate,lbd, &
            result%fval,grad,ok)
      end if
      result%gradient=grad
      allocate(full_lbd(m),alambda(n));full_lbd=0.0_dp
      if(allocated(lbd)) then
         do j=1,min(size(active),size(lbd))
            full_lbd(active(j))=lbd(j)
         end do
      end if
      result%lambda=full_lbd
      call transpose_constraint_product(isomat,full_lbd,alambda)
      result%alambda=alambda
      if(do_check) then
         call kkt_values(grad,result%constr_val,full_lbd,alambda,result%kkt)
      else
         result%kkt=0.0_dp
      end if
   end subroutine active_set_driver

   subroutine solve_active(yobs,isomat,active,xstart,opt,x,lbd,f,grad,ok,callback)
      real(dp),intent(in)::yobs(:),xstart(:)
      integer,intent(in)::isomat(:,:),active(:)
      type(active_set_options),intent(in)::opt
      real(dp),intent(out)::x(:),f,grad(:)
      real(dp),allocatable,intent(out)::lbd(:)
      logical,intent(out)::ok
      procedure(loss_callback),optional::callback
      integer,allocatable::labels(:),idx(:)
      integer::n,ncomp,c,i,nc
      real(dp),allocatable::w(:),z(:),wm(:,:),b(:),h(:,:),r(:)
      logical::lin_ok

      n=size(yobs);allocate(labels(n))
      call component_labels(n,isomat,active,labels,ncomp)
      allocate(w(n));
      if(allocated(opt%weights)) then
         if(size(opt%weights)/=n) then;ok=.false.;return;end if
         w=opt%weights
      else
         w=1.0_dp
      end if
      select case(opt%solver)
      case(ISO_LS,ISO_L1,ISO_QUANTILE,ISO_CHEBYSHEV,ISO_LP)
         do c=1,ncomp
            nc=count(labels==c);allocate(z(nc),idx(nc));idx=pack([(i,i=1,n)],labels==c)
            z=yobs(idx)
            select case(opt%solver)
            case(ISO_LS);b=[weighted_mean_value(z,w(idx))]
            case(ISO_L1);b=[weighted_median_value(z,w(idx))]
            case(ISO_QUANTILE)
               b=[weighted_fractile_value(z,w(idx),opt%bw/max(tiny(1.0_dp),opt%aw+opt%bw))]
            case(ISO_CHEBYSHEV);b=[weighted_midrange_value(z,w(idx))]
            case(ISO_LP);b=[weighted_lp_value(z,w,opt%p)]
            end select
            x(idx)=b(1)
            deallocate(z,idx,b)
         end do
         if(opt%solver==ISO_CHEBYSHEV) then
            call chebyshev_full_kkt(x,yobs,w,isomat,active,f,grad,lbd,ok)
         else
            call evaluate_builtin(x,yobs,opt,f,grad)
            call lagrange_multipliers(isomat,active,grad,lbd,ok)
         end if
      case(ISO_GLS)
         if(.not.allocated(opt%weight_matrix)) then;ok=.false.;return;end if
         if(size(opt%weight_matrix,1)/=n .or. size(opt%weight_matrix,2)/=n) then
            ok=.false.;return
         end if
         allocate(wm(n,ncomp));wm=0.0_dp
         do i=1,n;wm(i,labels(i))=1.0_dp;end do
         allocate(h(ncomp,ncomp),r(ncomp),b(ncomp))
         h=matmul(transpose(wm),matmul(opt%weight_matrix,wm))
         r=matmul(transpose(wm),matmul(opt%weight_matrix,yobs))
         call solve_linear(h,r,b,lin_ok)
         if(.not.lin_ok) then;ok=.false.;return;end if
         x=matmul(wm,b)
         grad=2.0_dp*matmul(opt%weight_matrix,x-yobs)
         f=dot_product(x-yobs,matmul(opt%weight_matrix,x-yobs))
         call lagrange_multipliers(isomat,active,grad,lbd,ok)
      case default
         if(present(callback)) then
            call reduced_bfgs(labels,ncomp,xstart,yobs,opt,x,f,grad,ok,callback)
         else
            call reduced_bfgs(labels,ncomp,xstart,yobs,opt,x,f,grad,ok)
         end if
         if(ok) call lagrange_multipliers(isomat,active,grad,lbd,ok)
      end select
   end subroutine solve_active

   subroutine evaluate_builtin(x,y,opt,f,g)
      real(dp),intent(in)::x(:),y(:)
      type(active_set_options),intent(in)::opt
      real(dp),intent(out)::f,g(:)
      real(dp),allocatable::w(:)
      real(dp)::d,ad,fac
      integer::i,n
      n=size(x);allocate(w(n))
      if(allocated(opt%weights)) then;w=opt%weights;else;w=1.0_dp;end if
      f=0.0_dp;g=0.0_dp
      select case(opt%solver)
      case(ISO_LS)
         g=2.0_dp*w*(x-y);f=sum(w*(x-y)**2)
      case(ISO_L1)
         do i=1,n
            d=x(i)-y(i);f=f+w(i)*abs(d)
            if(d>0.0_dp) g(i)=w(i)
            if(d<0.0_dp) g(i)=-w(i)
         end do
      case(ISO_QUANTILE)
         do i=1,n
            d=x(i)-y(i)
            if(d<=0.0_dp) then;f=f+w(i)*opt%aw*(-d);g(i)=-w(i)*opt%aw
            else;f=f+w(i)*opt%bw*d;g(i)=w(i)*opt%bw;end if
         end do
      case(ISO_CHEBYSHEV)
         call chebyshev_fg(x,y,w,f,g)
      case(ISO_POISSON)
         do i=1,n
            if(x(i)<=tiny(1.0_dp)) then;f=huge(1.0_dp);g(i)=-huge(1.0_dp);return;end if
            f=f+x(i)-y(i)*log(x(i));g(i)=1.0_dp-y(i)/x(i)
         end do
      case(ISO_LP)
         do i=1,n
            d=x(i)-y(i);ad=abs(d);f=f+w(i)*ad**opt%p
            if(ad>tiny(1.0_dp)) g(i)=opt%p*w(i)*sign(ad,d)*ad**(opt%p-1.0_dp)
         end do
      case(ISO_ASYLS)
         do i=1,n
            d=x(i)-y(i);if(d<0.0_dp) then;fac=opt%aw;else;fac=opt%bw;end if
            f=f+w(i)*d*d*fac;g(i)=2.0_dp*w(i)*d*fac
         end do
      case(ISO_L1EPS)
         do i=1,n
            d=x(i)-y(i);fac=sqrt(d*d+opt%eps);f=f+w(i)*fac;g(i)=w(i)*d/fac
         end do
      case(ISO_HUBER)
         do i=1,n
            d=x(i)-y(i);ad=abs(d)
            if(ad<2.0_dp*opt%eps) then
               f=f+w(i)*d*d/(4.0_dp*opt%eps);g(i)=w(i)*d/(2.0_dp*opt%eps)
            else
               f=f+w(i)*(ad-opt%eps);g(i)=w(i)*sign(1.0_dp,d)
            end if
         end do
      case(ISO_SILF)
         do i=1,n
            d=x(i)-y(i);ad=abs(d)
            if(ad<(1.0_dp-opt%beta)*opt%eps) then
               cycle
            else if(ad>(1.0_dp+opt%beta)*opt%eps) then
               f=f+w(i)*(ad-opt%eps);g(i)=w(i)*sign(1.0_dp,d);cycle
            else
               fac=(ad-(1.0_dp-opt%beta)*opt%eps)/(2.0_dp*opt%beta*opt%eps)
            end if
            f=f+w(i)*(ad-(1.0_dp-opt%beta)*opt%eps)**2/(4.0_dp*opt%beta*opt%eps)
            if(ad>tiny(1.0_dp)) g(i)=w(i)*fac*sign(1.0_dp,d)
         end do
      case default
         g=0.0_dp;f=0.0_dp
      end select
   end subroutine evaluate_builtin

   subroutine chebyshev_fg(x,y,w,f,g)
      real(dp),intent(in)::x(:),y(:),w(:)
      real(dp),intent(out)::f,g(:)
      real(dp),allocatable::dv(:)
      integer::i1,i2
      allocate(dv(size(x)));dv=w*(x-y)
      i1=maxloc(dv,dim=1);i2=minloc(dv,dim=1);f=maxval(abs(dv));g=0.0_dp
      if(abs(dv(i1))>=abs(dv(i2))) then;g(i1)=w(i1);else;g(i2)=-w(i2);end if
   end subroutine chebyshev_fg

   real(dp) function weighted_lp_value(y,w,pow) result(v)
      real(dp),intent(in)::y(:),w(:),pow
      real(dp)::lo,hi,mid,gmid,d
      integer::i,it
      if(pow<=1.0_dp+10.0_dp*epsilon(1.0_dp)) then
         v=weighted_median_value(y,w);return
      end if
      lo=minval(y);hi=maxval(y)
      if(hi-lo<=100.0_dp*epsilon(1.0_dp)*max(1.0_dp,max(abs(lo),abs(hi)))) then
         v=0.5_dp*(lo+hi);return
      end if
      do it=1,120
         mid=0.5_dp*(lo+hi);gmid=0.0_dp
         do i=1,size(y)
            d=mid-y(i)
            if(abs(d)>tiny(1.0_dp)) gmid=gmid+w(i)*sign(1.0_dp,d)*abs(d)**(pow-1.0_dp)
         end do
         if(gmid>0.0_dp) then;hi=mid;else;lo=mid;end if
         if(hi-lo<=1.0e-13_dp*(1.0_dp+abs(mid))) exit
      end do
      v=0.5_dp*(lo+hi)
   end function weighted_lp_value

   subroutine chebyshev_full_kkt(x,y,w,isomat,active,f,g,lbd,ok)
      real(dp),intent(in)::x(:),y(:),w(:)
      integer,intent(in)::isomat(:,:),active(:)
      real(dp),intent(out)::f,g(:)
      real(dp),allocatable,intent(out)::lbd(:)
      logical,intent(out)::ok
      real(dp),allocatable::dv(:),g1(:),g2(:),l1(:),l2(:)
      real(dp)::den
      integer::i1,i2
      allocate(dv(size(x)),g1(size(x)),g2(size(x)))
      dv=w*(x-y);i1=maxloc(dv,dim=1);i2=minloc(dv,dim=1)
      f=maxval(abs(dv));g1=0.0_dp;g2=0.0_dp
      g1(i1)=w(i1);g2(i2)=-w(i2)
      call lagrange_multipliers(isomat,active,g1,l1,ok)
      if(.not.ok) return
      call lagrange_multipliers(isomat,active,g2,l2,ok)
      if(.not.ok) return
      den=w(i1)+w(i2)
      if(abs(den)<=tiny(1.0_dp)) then
         g=0.5_dp*(g1+g2);lbd=0.5_dp*(l1+l2)
      else
         g=(w(i2)*g1+w(i1)*g2)/den
         lbd=(w(i2)*l1+w(i1)*l2)/den
      end if
   end subroutine chebyshev_full_kkt

   subroutine reduced_bfgs(labels,ncomp,xstart,yobs,opt,x,f,g,ok,callback)
      integer,intent(in)::labels(:),ncomp
      real(dp),intent(in)::xstart(:),yobs(:)
      type(active_set_options),intent(in)::opt
      real(dp),intent(out)::x(:),f,g(:)
      logical,intent(out)::ok
      procedure(loss_callback),optional::callback
      real(dp),allocatable::b(:),gb(:),h(:,:),dir(:),bn(:),gn(:),s(:),yv(:), &
         xx(:),gx(:),gxn(:),iunit(:,:)
      real(dp)::fn,step,gd,ys,rho,tol
      integer::i,j,it,n
      n=size(x);allocate(b(ncomp),gb(ncomp),h(ncomp,ncomp),dir(ncomp),bn(ncomp), &
         gn(ncomp),s(ncomp),yv(ncomp),xx(n),gx(n),gxn(n),iunit(ncomp,ncomp))
      do j=1,ncomp
         b(j)=sum(pack(xstart,labels==j))/real(count(labels==j),dp)
      end do
      call expand_components(labels,b,xx)
      if(present(callback)) then;call callback(xx,f,gx);else;call evaluate_builtin(xx,yobs,opt,f,gx);end if
      call reduce_gradient(labels,ncomp,gx,gb)
      h=0.0_dp;iunit=0.0_dp
      do i=1,ncomp;h(i,i)=1.0_dp;iunit(i,i)=1.0_dp;end do
      tol=opt%bfgs_tol;ok=.false.
      do it=1,opt%bfgs_maxiter
         if(maxval(abs(gb))<=tol*(1.0_dp+abs(f))) then;ok=.true.;exit;end if
         dir=-matmul(h,gb);gd=dot_product(gb,dir)
         if(gd>=0.0_dp) then;dir=-gb;h=iunit;gd=-dot_product(gb,gb);end if
         step=1.0_dp
         do
            bn=b+step*dir;call expand_components(labels,bn,xx)
            if(present(callback)) then;call callback(xx,fn,gxn);else;call evaluate_builtin(xx,yobs,opt,fn,gxn);end if
            if(fn < huge(1.0_dp)/10.0_dp .and. fn<=f+1.0e-4_dp*step*gd) exit
            step=0.5_dp*step
            if(step<1.0e-14_dp) exit
         end do
         if(step<1.0e-14_dp) then;ok=maxval(abs(gb))<1.0e-6_dp;exit;end if
         call reduce_gradient(labels,ncomp,gxn,gn)
         s=bn-b;yv=gn-gb;ys=dot_product(yv,s)
         if(ys>sqrt(epsilon(1.0_dp))*dot_product(s,s)) then
            rho=1.0_dp/ys
            h=matmul(iunit-rho*outer(s,yv),matmul(h,iunit-rho*outer(yv,s)))+rho*outer(s,s)
         else
            h=iunit
         end if
         b=bn;gb=gn;gx=gxn;f=fn
      end do
      call expand_components(labels,b,x)
      if(present(callback)) then;call callback(x,f,g);else;call evaluate_builtin(x,yobs,opt,f,g);end if
      if(.not.ok) ok=maxval(abs(gb))<1.0e-6_dp
   contains
      function outer(a,bv) result(c)
         real(dp),intent(in)::a(:),bv(:)
         real(dp)::c(size(a),size(bv))
         integer::ii,jj
         do jj=1,size(bv);do ii=1,size(a);c(ii,jj)=a(ii)*bv(jj);end do;end do
      end function outer
   end subroutine reduced_bfgs

   subroutine expand_components(labels,b,x)
      integer,intent(in)::labels(:)
      real(dp),intent(in)::b(:)
      real(dp),intent(out)::x(:)
      integer::i
      do i=1,size(x);x(i)=b(labels(i));end do
   end subroutine expand_components

   subroutine reduce_gradient(labels,ncomp,g,gb)
      integer,intent(in)::labels(:),ncomp
      real(dp),intent(in)::g(:)
      real(dp),intent(out)::gb(:)
      integer::i
      if(size(gb)/=ncomp) error stop "reduce_gradient: bad size"
      gb=0.0_dp
      do i=1,size(g);gb(labels(i))=gb(labels(i))+g(i);end do
   end subroutine reduce_gradient

   subroutine active_indices(ax,ups,active)
      real(dp),intent(in)::ax(:),ups
      integer,allocatable,intent(out)::active(:)
      integer::i,k
      allocate(active(count(abs(ax)<ups)));k=0
      do i=1,size(ax);if(abs(ax(i))<ups) then;k=k+1;active(k)=i;end if;end do
   end subroutine active_indices

   subroutine remove_position(a,pos)
      integer,allocatable,intent(inout)::a(:)
      integer,intent(in)::pos
      integer,allocatable::b(:)
      integer::n
      n=size(a);allocate(b(max(0,n-1)))
      if(pos>1)b(1:pos-1)=a(1:pos-1)
      if(pos<n)b(pos:n-1)=a(pos+1:n)
      call move_alloc(b,a)
   end subroutine remove_position

   subroutine add_unique_sorted(a,value)
      integer,allocatable,intent(inout)::a(:)
      integer,intent(in)::value
      integer,allocatable::b(:)
      integer::i,j,n
      if(any(a==value)) return
      n=size(a);allocate(b(n+1));b(1:n)=a;b(n+1)=value
      do i=2,n+1
         value_insert: block
            integer::key
            key=b(i);j=i-1
            do while(j>=1)
               if(b(j)<=key) exit
               b(j+1)=b(j);j=j-1
            end do
            b(j+1)=key
         end block value_insert
      end do
      call move_alloc(b,a)
   end subroutine add_unique_sorted
end module isotone_active
