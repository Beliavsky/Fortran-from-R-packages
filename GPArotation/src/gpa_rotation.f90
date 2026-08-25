module gpa_rotation
  use gpa_kinds, only: dp
  use gpa_linalg, only: eye, inverse_matrix, polar_orthogonal, normalize_columns
  use gpa_linalg, only: frobenius_norm, random_orthogonal
  use gpa_criteria, only: criterion_options, criterion_result, evaluate_criterion
  implicit none
  private

  type, public :: rotation_options
    real(dp) :: eps = 1.0e-5_dp
    integer :: maxit = 2000
    character(len=12) :: algorithm = 'bb'
    integer :: fwindow = 10
    integer :: normalize = 0       ! 0 none, 1 Kaiser, 2 Cureton-Mulaik
    integer :: random_starts = 0
  end type rotation_options

  type, public :: rotation_result
    real(dp), allocatable :: loadings(:,:)
    real(dp), allocatable :: th(:,:)
    real(dp), allocatable :: phi(:,:)
    real(dp), allocatable :: gq(:,:)
    real(dp), allocatable :: table(:,:)  ! iter, f, log10(s), alpha
    real(dp), allocatable :: qvalues(:)
    real(dp) :: objective = huge(1.0_dp)
    integer :: iterations = 0
    integer :: info = 0
    logical :: converged = .false.
    logical :: orthogonal = .true.
    character(len=64) :: method = ''
  end type rotation_result

  public :: gpforth, gpfoblq, rotate_random_starts, lp_rotate
  public :: normalizing_weights, apply_normalization

contains

  subroutine gpforth(a, method, result, crit_opts, rot_opts, tinit)
    real(dp), intent(in) :: a(:,:)
    character(len=*), intent(in) :: method
    type(rotation_result), intent(out) :: result
    type(criterion_options), intent(in), optional :: crit_opts
    type(rotation_options), intent(in), optional :: rot_opts
    real(dp), intent(in), optional :: tinit(:,:)
    type(rotation_options) :: o
    type(criterion_options) :: co
    type(criterion_result) :: cr, crt
    real(dp), allocatable :: b(:,:), t(:,:), ttrial(:,:), l(:,:), g(:,:), gp(:,:)
    real(dp), allocatable :: tprev(:,:), gpprev(:,:), hist(:,:), wrow(:)
    real(dp) :: alpha, s, targetf, denom, bbnum, oldf
    real(dp) :: m(size(a,2),size(a,2)), sym(size(a,2),size(a,2))
    real(dp) :: ww(size(a,2),size(a,2)), invx(size(a,2),size(a,2))
    integer :: k, iter, ii, j0, info, nused
    logical :: haveprev

    if(present(rot_opts)) o=rot_opts
    if(present(crit_opts)) co=crit_opts
    k=size(a,2)
    call init_rotation_result(result,a,k,.true.)
    if(k<2 .or. size(a,1)<1) then; result%info=-1; return; end if
    allocate(b(size(a,1),k),t(k,k),ttrial(k,k),l(size(a,1),k),g(k,k),gp(k,k))
    allocate(tprev(k,k),gpprev(k,k),hist(o%maxit+1,4),wrow(size(a,1)))
    call apply_normalization(a,o%normalize,b,wrow)
    if(present(tinit)) then
      if(any(shape(tinit)/=[k,k])) then; result%info=-2; return; end if
      t=tinit
    else
      t=eye(k)
    end if
    l=matmul(b,t)
    call evaluate_criterion(method,l,cr,co)
    if(cr%info/=0) then; result%info=cr%info; return; end if
    g=matmul(transpose(b),cr%gq)
    alpha=1.0_dp; haveprev=.false.; hist=0.0_dp; oldf=cr%f

    do iter=0,o%maxit
      m=matmul(transpose(t),g); sym=0.5_dp*(m+transpose(m))
      gp=g-matmul(t,sym); s=frobenius_norm(gp)
      hist(iter+1,:)=[real(iter,dp),oldf,log10(max(s,tiny(1.0_dp))),alpha]
      if(s<o%eps) exit
      if((trim(lower(o%algorithm))=='bb' .or. trim(lower(o%algorithm))=='cayley') .and. haveprev) then
        denom=sum((gp-gpprev)**2); bbnum=sum((t-tprev)**2)
        if(denom>tiny(1.0_dp)) then
          alpha=bbnum/max(abs(sum((t-tprev)*(gp-gpprev))),tiny(1.0_dp))
          alpha=max(1.0e-10_dp,min(alpha,20.0_dp))
        end if
      else
        alpha=2.0_dp*alpha
      end if
      j0=max(1,iter+2-max(1,o%fwindow)); targetf=maxval(hist(j0:iter+1,2))
      do ii=0,10
        if(trim(lower(o%algorithm))=='cayley') then
          ww=matmul(gp,transpose(t))-matmul(t,transpose(gp))
          call inverse_matrix(eye(k)+(alpha/2.0_dp)*ww,invx,info)
          if(info/=0) then; alpha=alpha/2.0_dp; cycle; end if
          ttrial=matmul(matmul(invx,eye(k)-(alpha/2.0_dp)*ww),t)
        else
          call polar_orthogonal(t-alpha*gp,ttrial,info)
          if(info/=0) then; alpha=alpha/2.0_dp; cycle; end if
        end if
        l=matmul(b,ttrial)
        call evaluate_criterion(method,l,crt,co)
        if(crt%info/=0) then; result%info=crt%info; return; end if
        if((targetf-crt%f)>0.5_dp*s*s*alpha) exit
        alpha=alpha/2.0_dp
      end do
      tprev=t; gpprev=gp; haveprev=.true.; t=ttrial; cr=crt; oldf=cr%f
      g=matmul(transpose(b),cr%gq)
    end do

    nused=min(iter,o%maxit)+1
    result%iterations=min(iter,o%maxit); result%converged=(s<o%eps)
    result%th=t; result%phi=eye(k); result%objective=cr%f; result%gq=cr%gq
    result%method=cr%method; result%table=hist(1:nused,:); result%info=merge(0,1,result%converged)
    l=matmul(b,t)
    do ii=1,size(a,1); result%loadings(ii,:)=l(ii,:)*wrow(ii); end do
  end subroutine gpforth

  subroutine gpfoblq(a, method, result, crit_opts, rot_opts, tinit)
    real(dp), intent(in) :: a(:,:)
    character(len=*), intent(in) :: method
    type(rotation_result), intent(out) :: result
    type(criterion_options), intent(in), optional :: crit_opts
    type(rotation_options), intent(in), optional :: rot_opts
    real(dp), intent(in), optional :: tinit(:,:)
    type(rotation_options) :: o
    type(criterion_options) :: co
    type(criterion_result) :: cr, crt
    real(dp), allocatable :: b(:,:),t(:,:),tt(:,:),ti(:,:),tit(:,:),l(:,:),g(:,:),gp(:,:)
    real(dp), allocatable :: tp(:,:),gpp(:,:),hist(:,:),wrow(:)
    real(dp) :: alpha,s,targetf,denom,bbnum,oldf,nrm
    real(dp) :: d(size(a,2))
    integer :: k,iter,ii,j,j0,info,nused
    logical :: haveprev
    if(present(rot_opts)) o=rot_opts
    if(present(crit_opts)) co=crit_opts
    if(trim(lower(o%algorithm))=='cayley') then
      call init_rotation_result(result,a,size(a,2),.false.); result%info=-3; return
    end if
    k=size(a,2); call init_rotation_result(result,a,k,.false.)
    if(k<2) then; result%info=-1; return; end if
    allocate(b(size(a,1),k),t(k,k),tt(k,k),ti(k,k),tit(k,k),l(size(a,1),k),g(k,k),gp(k,k))
    allocate(tp(k,k),gpp(k,k),hist(o%maxit+1,4),wrow(size(a,1)))
    call apply_normalization(a,o%normalize,b,wrow)
    if(present(tinit)) then
      if(any(shape(tinit)/=[k,k])) then; result%info=-2; return; end if
      t=tinit
    else
      t=eye(k)
    end if
    call inverse_matrix(t,ti,info); if(info/=0) then; result%info=info; return; end if
    l=matmul(b,transpose(ti)); call evaluate_criterion(method,l,cr,co)
    if(cr%info/=0) then; result%info=cr%info; return; end if
    g=-transpose(matmul(matmul(transpose(l),cr%gq),ti))
    alpha=1.0_dp; haveprev=.false.; hist=0.0_dp; oldf=cr%f
    do iter=0,o%maxit
      d=sum(t*g,dim=1)
      gp=g
      do j=1,k; gp(:,j)=gp(:,j)-t(:,j)*d(j); end do
      s=frobenius_norm(gp)
      hist(iter+1,:)=[real(iter,dp),oldf,log10(max(s,tiny(1.0_dp))),alpha]
      if(s<o%eps) exit
      if(trim(lower(o%algorithm))=='bb' .and. haveprev) then
        denom=sum((gp-gpp)**2); bbnum=sum((t-tp)**2)
        if(denom>tiny(1.0_dp)) then
          alpha=bbnum/max(abs(sum((t-tp)*(gp-gpp))),tiny(1.0_dp))
          alpha=max(1.0e-10_dp,min(alpha,20.0_dp))
        end if
      else
        alpha=2.0_dp*alpha
      end if
      j0=max(1,iter+2-max(1,o%fwindow)); targetf=maxval(hist(j0:iter+1,2))
      do ii=0,10
        tt=t-alpha*gp
        do j=1,k
          nrm=sqrt(sum(tt(:,j)*tt(:,j)))
          if(nrm<=tiny(1.0_dp)) exit
          tt(:,j)=tt(:,j)/nrm
        end do
        call inverse_matrix(tt,tit,info)
        if(info/=0) then; alpha=alpha/2.0_dp; cycle; end if
        l=matmul(b,transpose(tit)); call evaluate_criterion(method,l,crt,co)
        if(crt%info/=0) then; result%info=crt%info; return; end if
        if((targetf-crt%f)>0.5_dp*s*s*alpha) exit
        alpha=alpha/2.0_dp
      end do
      tp=t; gpp=gp; haveprev=.true.; t=tt; ti=tit; cr=crt; oldf=cr%f
      g=-transpose(matmul(matmul(transpose(l),cr%gq),ti))
    end do
    nused=min(iter,o%maxit)+1
    result%iterations=min(iter,o%maxit); result%converged=(s<o%eps)
    result%th=t; result%phi=matmul(transpose(t),t); result%objective=cr%f
    result%gq=cr%gq; result%method=cr%method; result%table=hist(1:nused,:)
    result%info=merge(0,1,result%converged)
    call inverse_matrix(t,ti,info); l=matmul(b,transpose(ti))
    do ii=1,size(a,1); result%loadings(ii,:)=l(ii,:)*wrow(ii); end do
  end subroutine gpfoblq

  subroutine rotate_random_starts(a,method,orthogonal,nstarts,result,crit_opts,rot_opts)
    real(dp),intent(in)::a(:,:)
    character(len=*),intent(in)::method
    logical,intent(in)::orthogonal
    integer,intent(in)::nstarts
    type(rotation_result),intent(out)::result
    type(criterion_options),intent(in),optional::crit_opts
    type(rotation_options),intent(in),optional::rot_opts
    type(rotation_result)::cur
    type(rotation_options)::o
    real(dp),allocatable::t(:,:),qv(:)
    integer::s,info,k,ns
    if(present(rot_opts)) o=rot_opts
    k=size(a,2); ns=max(1,nstarts); allocate(t(k,k),qv(ns))
    result%objective=huge(1.0_dp)
    do s=1,ns
      if(s==1 .and. nstarts<=0) then
        t=eye(k)
      else
        call random_orthogonal(k,t,info); if(info/=0) cycle
      end if
      if(orthogonal) then
        call gpforth(a,method,cur,crit_opts,o,t)
      else
        call gpfoblq(a,method,cur,crit_opts,o,t)
      end if
      qv(s)=cur%objective
      if(cur%objective<result%objective) result=cur
    end do
    result%qvalues=qv
  end subroutine rotate_random_starts

  subroutine lp_rotate(a,p,orthogonal,result,rot_opts,tinit)
    real(dp),intent(in)::a(:,:),p
    logical,intent(in)::orthogonal
    type(rotation_result),intent(out)::result
    type(rotation_options),intent(in),optional::rot_opts
    real(dp),intent(in),optional::tinit(:,:)
    type(rotation_options)::o,inner_o
    type(criterion_options)::co
    type(rotation_result)::rr
    real(dp),allocatable::t(:,:),lprev(:,:),w(:,:)
    real(dp)::epslp
    integer::it,k
    if(present(rot_opts)) o=rot_opts
    k=size(a,2); allocate(t(k,k),lprev(size(a,1),k),w(size(a,1),k))
    if(present(tinit)) then; t=tinit; else; t=eye(k); end if
    epslp=max(o%eps,1.0e-10_dp); lprev=huge(1.0_dp)
    inner_o=o; inner_o%maxit=5; inner_o%random_starts=0
    do it=1,o%maxit
      if(it==1) then
        if(orthogonal) then
          w=(matmul(a,t)**2+epslp)**(p/2.0_dp-1.0_dp)
        else
          call oblique_loadings(a,t,lprev)
          w=(lprev*lprev+epslp)**(p/2.0_dp-1.0_dp)
        end if
      end if
      co%weight=w
      if(orthogonal) then
        call gpforth(a,'lp.wls',rr,co,inner_o,t)
      else
        call gpfoblq(a,'lp.wls',rr,co,inner_o,t)
      end if
      if(maxval(abs(rr%loadings-lprev))<o%eps) then
        result=rr; result%converged=.true.; exit
      end if
      lprev=rr%loadings; w=(lprev*lprev+epslp)**(p/2.0_dp-1.0_dp); t=rr%th; result=rr
    end do
    result%objective=sum(abs(result%loadings)**p)/real(size(a,1),dp)
    result%iterations=it
  contains
    subroutine oblique_loadings(x,tt,ll)
      real(dp),intent(in)::x(:,:),tt(:,:)
      real(dp),intent(out)::ll(size(x,1),size(x,2))
      real(dp)::invtt(size(tt,1),size(tt,2)); integer::ierr
      call inverse_matrix(tt,invtt,ierr)
      if(ierr==0) then; ll=matmul(x,transpose(invtt)); else; ll=0.0_dp; end if
    end subroutine oblique_loadings
  end subroutine lp_rotate

  subroutine normalizing_weights(a,mode,w)
    real(dp),intent(in)::a(:,:)
    integer,intent(in)::mode
    real(dp),intent(out)::w(size(a,1))
    real(dp)::h,fpls,acosi,alpha,num,dem,wt
    integer::i,m
    if(mode==0) then; w=1.0_dp; return; end if
    m=size(a,2); acosi=acos(real(m,dp)**(-0.5_dp))
    do i=1,size(a,1)
      h=max(sqrt(sum(a(i,:)*a(i,:))),epsilon(1.0_dp))
      if(mode==1) then
        w(i)=h
      else
        fpls=max(-1.0_dp,min(1.0_dp,a(i,1)/h)); alpha=acos(abs(fpls))
        if(abs(fpls)<real(m,dp)**(-0.5_dp)) then; dem=acosi-0.5_dp*acos(-1.0_dp); else; dem=acosi; end if
        num=acosi-alpha
        if(abs(dem)<=epsilon(1.0_dp)) then; wt=0.001_dp; else; wt=cos((num/dem)*0.5_dp*acos(-1.0_dp))**2+0.001_dp; end if
        w(i)=h/wt
      end if
    end do
  end subroutine normalizing_weights

  subroutine apply_normalization(a,mode,b,w)
    real(dp),intent(in)::a(:,:)
    integer,intent(in)::mode
    real(dp),intent(out)::b(size(a,1),size(a,2)),w(size(a,1))
    integer::i
    call normalizing_weights(a,mode,w)
    do i=1,size(a,1); b(i,:)=a(i,:)/w(i); end do
  end subroutine apply_normalization

  subroutine init_rotation_result(r,a,k,orth)
    type(rotation_result),intent(out)::r
    real(dp),intent(in)::a(:,:)
    integer,intent(in)::k
    logical,intent(in)::orth
    allocate(r%loadings(size(a,1),k),r%th(k,k),r%phi(k,k),r%gq(size(a,1),k))
    r%loadings=0.0_dp; r%th=eye(k); r%phi=eye(k); r%gq=0.0_dp
    r%orthogonal=orth; r%objective=huge(1.0_dp); r%info=0
  end subroutine init_rotation_result

  pure function lower(s) result(t)
    character(len=*),intent(in)::s
    character(len=len(s))::t
    integer::i,c
    do i=1,len(s); c=iachar(s(i:i)); if(c>=65 .and. c<=90) then; t(i:i)=achar(c+32); else; t(i:i)=s(i:i); end if; end do
  end function lower

end module gpa_rotation
