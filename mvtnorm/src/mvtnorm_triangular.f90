! SPDX-License-Identifier: GPL-2.0-only
module mvtnorm_triangular
  use mvtnorm_kinds, only : dp
  use mvtnorm_linalg, only : cholesky_lower, inverse_lower, inverse_spd, &
    covariance_to_correlation, symmetrize, identity_matrix
  implicit none
  private
  public :: pack_lower, unpack_lower, lower_tri_size
  public :: lower_multiply, lower_solve, lower_logdet
  public :: tcrossprod_lower, crossprod_lower
  public :: chol2cov, cov2chol, chol2invchol, invchol2chol
  public :: invchol2cov, cov2invchol, invchol2pre, chol2pre
  public :: chol2cor, invchol2cor, chol2pc, invchol2pc
  public :: dchol, invchold, vectrick, permute_cholesky, deperma_score, destandardize_score

contains

  integer function lower_tri_size(n,include_diagonal) result(k)
    integer,intent(in)::n
    logical,intent(in),optional::include_diagonal
    logical::diag
    diag=.true.; if(present(include_diagonal)) diag=include_diagonal
    if(diag) then; k=n*(n+1)/2; else; k=n*(n-1)/2; end if
  end function lower_tri_size

  function pack_lower(a,include_diagonal,byrow) result(v)
    real(dp),intent(in)::a(:,:)
    logical,intent(in),optional::include_diagonal,byrow
    real(dp),allocatable::v(:)
    logical::diag,row
    integer::n,i,j,k
    n=size(a,1); diag=.true.; row=.false.
    if(present(include_diagonal)) diag=include_diagonal
    if(present(byrow)) row=byrow
    allocate(v(lower_tri_size(n,diag))); k=0
    if(row) then
      do i=1,n
        do j=1,i
          if(.not.diag .and. i==j) cycle
          k=k+1; v(k)=a(i,j)
        end do
      end do
    else
      do j=1,n
        do i=j,n
          if(.not.diag .and. i==j) cycle
          k=k+1; v(k)=a(i,j)
        end do
      end do
    end if
  end function pack_lower

  function unpack_lower(v,n,include_diagonal,byrow,unit_diagonal) result(a)
    real(dp),intent(in)::v(:)
    integer,intent(in)::n
    logical,intent(in),optional::include_diagonal,byrow,unit_diagonal
    real(dp),allocatable::a(:,:)
    logical::diag,row,unitd
    integer::i,j,k
    diag=.true.; row=.false.; unitd=.false.
    if(present(include_diagonal)) diag=include_diagonal
    if(present(byrow)) row=byrow
    if(present(unit_diagonal)) unitd=unit_diagonal
    allocate(a(n,n)); a=0.0_dp; k=0
    if(row) then
      do i=1,n
        do j=1,i
          if(.not.diag .and. i==j) cycle
          k=k+1; if(k<=size(v)) a(i,j)=v(k)
        end do
      end do
    else
      do j=1,n
        do i=j,n
          if(.not.diag .and. i==j) cycle
          k=k+1; if(k<=size(v)) a(i,j)=v(k)
        end do
      end do
    end if
    if(unitd .or. .not.diag) then
      do i=1,n; a(i,i)=1.0_dp; end do
    end if
  end function unpack_lower

  function lower_multiply(l,y,transpose_l) result(z)
    real(dp),intent(in)::l(:,:),y(:,:)
    logical,intent(in),optional::transpose_l
    real(dp),allocatable::z(:,:)
    logical::tr
    tr=.false.; if(present(transpose_l)) tr=transpose_l
    if(tr) then; z=matmul(transpose(l),y); else; z=matmul(l,y); end if
  end function lower_multiply

  function lower_solve(l,b,transpose_l,ok) result(x)
    real(dp),intent(in)::l(:,:),b(:,:)
    logical,intent(in),optional::transpose_l
    logical,intent(out),optional::ok
    real(dp),allocatable::x(:,:)
    logical::lok,tr
    integer::i,j,n
    n=size(l,1); tr=.false.; if(present(transpose_l)) tr=transpose_l
    allocate(x(n,size(b,2))); x=b; lok=.true.
    if(.not.tr) then
      do i=1,n
        if(i>1) x(i,:)=x(i,:)-matmul(l(i,1:i-1),x(1:i-1,:))
        if(abs(l(i,i))<=tiny(1.0_dp)) then; lok=.false.; exit; end if
        x(i,:)=x(i,:)/l(i,i)
      end do
    else
      do i=n,1,-1
        if(i<n) then
          do j=1,size(b,2)
            x(i,j)=x(i,j)-dot_product(l(i+1:n,i),x(i+1:n,j))
          end do
        end if
        if(abs(l(i,i))<=tiny(1.0_dp)) then; lok=.false.; exit; end if
        x(i,:)=x(i,:)/l(i,i)
      end do
    end if
    if(present(ok)) ok=lok
  end function lower_solve

  real(dp) function lower_logdet(l) result(v)
    real(dp),intent(in)::l(:,:)
    integer::i
    v=0.0_dp
    do i=1,min(size(l,1),size(l,2)); v=v+log(abs(l(i,i))); end do
  end function lower_logdet

  function tcrossprod_lower(l,diag_only) result(a)
    real(dp),intent(in)::l(:,:)
    logical,intent(in),optional::diag_only
    real(dp),allocatable::a(:,:)
    logical::d
    integer::i,n
    d=.false.; if(present(diag_only)) d=diag_only
    n=size(l,1)
    if(d) then
      allocate(a(n,1)); do i=1,n; a(i,1)=dot_product(l(i,:),l(i,:)); end do
    else
      a=matmul(l,transpose(l)); call symmetrize(a)
    end if
  end function tcrossprod_lower

  function crossprod_lower(l,diag_only) result(a)
    real(dp),intent(in)::l(:,:)
    logical,intent(in),optional::diag_only
    real(dp),allocatable::a(:,:)
    logical::d
    integer::i,n
    d=.false.; if(present(diag_only)) d=diag_only
    n=size(l,2)
    if(d) then
      allocate(a(n,1)); do i=1,n; a(i,1)=dot_product(l(:,i),l(:,i)); end do
    else
      a=matmul(transpose(l),l); call symmetrize(a)
    end if
  end function crossprod_lower

  function chol2cov(l) result(cov)
    real(dp),intent(in)::l(:,:)
    real(dp),allocatable::cov(:,:)
    cov=matmul(l,transpose(l)); call symmetrize(cov)
  end function chol2cov

  function cov2chol(cov,ok,message) result(l)
    real(dp),intent(in)::cov(:,:)
    logical,intent(out),optional::ok
    character(len=*),intent(out),optional::message
    real(dp),allocatable::l(:,:)
    logical::lok; character(len=256)::msg
    call cholesky_lower(cov,l,lok,msg)
    if(present(ok)) ok=lok; if(present(message)) message=msg
  end function cov2chol

  function chol2invchol(l,ok) result(linv)
    real(dp),intent(in)::l(:,:)
    logical,intent(out),optional::ok
    real(dp),allocatable::linv(:,:)
    logical::lok
    call inverse_lower(l,linv,lok); if(present(ok)) ok=lok
  end function chol2invchol

  function invchol2chol(linv,ok) result(l)
    real(dp),intent(in)::linv(:,:)
    logical,intent(out),optional::ok
    real(dp),allocatable::l(:,:)
    logical::lok
    call inverse_lower(linv,l,lok); if(present(ok)) ok=lok
  end function invchol2chol

  function invchol2cov(linv) result(cov)
    real(dp),intent(in)::linv(:,:)
    real(dp),allocatable::cov(:,:),l(:,:)
    logical::ok
    call inverse_lower(linv,l,ok); cov=matmul(l,transpose(l)); call symmetrize(cov)
  end function invchol2cov

  function cov2invchol(cov,ok,message) result(linv)
    real(dp),intent(in)::cov(:,:)
    logical,intent(out),optional::ok
    character(len=*),intent(out),optional::message
    real(dp),allocatable::linv(:,:),l(:,:)
    logical::lok; character(len=256)::msg
    call cholesky_lower(cov,l,lok,msg)
    if(lok) call inverse_lower(l,linv,lok)
    if(.not.allocated(linv)) then; allocate(linv(size(cov,1),size(cov,2))); linv=0.0_dp; end if
    if(present(ok)) ok=lok; if(present(message)) message=msg
  end function cov2invchol

  function invchol2pre(linv) result(pre)
    real(dp),intent(in)::linv(:,:)
    real(dp),allocatable::pre(:,:)
    pre=matmul(transpose(linv),linv); call symmetrize(pre)
  end function invchol2pre

  function chol2pre(l) result(pre)
    real(dp),intent(in)::l(:,:)
    real(dp),allocatable::pre(:,:),linv(:,:)
    logical::ok
    call inverse_lower(l,linv,ok); pre=matmul(transpose(linv),linv); call symmetrize(pre)
  end function chol2pre

  function chol2cor(l) result(cor)
    real(dp),intent(in)::l(:,:)
    real(dp),allocatable::cor(:,:),cov(:,:),sd(:)
    logical::ok; character(len=256)::msg
    cov=chol2cov(l); call covariance_to_correlation(cov,cor,sd,ok,msg)
  end function chol2cor

  function invchol2cor(linv) result(cor)
    real(dp),intent(in)::linv(:,:)
    real(dp),allocatable::cor(:,:),cov(:,:),sd(:)
    logical::ok; character(len=256)::msg
    cov=invchol2cov(linv); call covariance_to_correlation(cov,cor,sd,ok,msg)
  end function invchol2cor

  function invchol2pc(linv) result(pc)
    real(dp),intent(in)::linv(:,:)
    real(dp),allocatable::pc(:,:),pre(:,:)
    integer::i,j,n
    pre=invchol2pre(linv); n=size(pre,1); allocate(pc(n,n)); pc=0.0_dp
    do i=1,n
      do j=1,n
        if(i/=j) pc(i,j)=-pre(i,j)/sqrt(pre(i,i)*pre(j,j))
      end do
    end do
    call symmetrize(pc)
  end function invchol2pc

  function chol2pc(l) result(pc)
    real(dp),intent(in)::l(:,:)
    real(dp),allocatable::pc(:,:),linv(:,:)
    logical::ok
    call inverse_lower(l,linv,ok); pc=invchol2pc(linv)
  end function chol2pc

  function dchol(l,d) result(out)
    real(dp),intent(in)::l(:,:)
    real(dp),intent(in),optional::d(:)
    real(dp),allocatable::out(:,:)
    real(dp),allocatable::scale(:)
    integer::i,n
    n=size(l,1); allocate(scale(n)); out=l
    if(present(d)) then; scale=d; else; do i=1,n; scale(i)=1.0_dp/sqrt(dot_product(l(i,:),l(i,:))); end do; end if
    do i=1,n; out(i,:)=scale(i)*out(i,:); end do
  end function dchol

  function invchold(linv,d) result(out)
    real(dp),intent(in)::linv(:,:)
    real(dp),intent(in),optional::d(:)
    real(dp),allocatable::out(:,:),l(:,:),scale(:)
    logical::ok
    integer::j,n
    n=size(linv,1); allocate(scale(n)); out=linv
    if(present(d)) then
      scale=d
    else
      call inverse_lower(linv,l,ok)
      do j=1,n; scale(j)=sqrt(dot_product(l(:,j),l(:,j))); end do
    end if
    do j=1,n; out(:,j)=out(:,j)*scale(j); end do
  end function invchold

  function vectrick(c,s,a,transpose_c,transpose_a) result(v)
    real(dp),intent(in)::c(:,:),s(:,:),a(:,:)
    logical,intent(in),optional::transpose_c,transpose_a
    real(dp),allocatable::v(:,:)
    real(dp),allocatable::cc(:,:),aa(:,:)
    logical::tc,ta
    tc=.true.; ta=.true.
    if(present(transpose_c)) tc=transpose_c
    if(present(transpose_a)) ta=transpose_a
    allocate(cc(size(c,1),size(c,2)),aa(size(a,1),size(a,2)))
    if(tc) then
      cc=transpose(c)
    else
      cc=c
    end if
    if(ta) then
      aa=transpose(a)
    else
      aa=a
    end if
    v=matmul(cc,matmul(s,aa))
  end function vectrick

  function permute_cholesky(l,perm,ok,message) result(lp)
    real(dp),intent(in)::l(:,:)
    integer,intent(in)::perm(:)
    logical,intent(out),optional::ok
    character(len=*),intent(out),optional::message
    real(dp),allocatable::lp(:,:),cov(:,:),pcov(:,:)
    logical::lok
    character(len=256)::msg
    integer::i,j,n
    n=size(l,1); cov=chol2cov(l); allocate(pcov(n,n))
    do i=1,n
      do j=1,n
        pcov(i,j)=cov(perm(i),perm(j))
      end do
    end do
    call cholesky_lower(pcov,lp,lok,msg)
    if(present(ok)) ok=lok
    if(present(message)) message=msg
  end function permute_cholesky

  function deperma_score(chol,perm,score_permuted,step) result(score)
    real(dp),intent(in)::chol(:,:),score_permuted(:)
    integer,intent(in)::perm(:)
    real(dp),intent(in),optional::step
    real(dp),allocatable::score(:)
    real(dp),allocatable::theta(:),tp(:),tm(:),lp(:,:),lm(:,:),vp(:),vm(:)
    real(dp)::h
    logical::ok
    character(len=256)::message
    integer::j,n,k
    n=size(chol,1); theta=pack_lower(chol,.true.,.false.); k=size(theta)
    allocate(score(k)); score=0.0_dp
    do j=1,k
      tp=theta; tm=theta
      h=sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(theta(j)))
      if(present(step)) h=step
      tp(j)=tp(j)+h; tm(j)=tm(j)-h
      lp=permute_cholesky(unpack_lower(tp,n,.true.,.false.),perm,ok,message)
      if(.not.ok) cycle
      lm=permute_cholesky(unpack_lower(tm,n,.true.,.false.),perm,ok,message)
      if(.not.ok) cycle
      vp=pack_lower(lp,.true.,.false.); vm=pack_lower(lm,.true.,.false.)
      score(j)=dot_product(score_permuted,(vp-vm)/(2.0_dp*h))
    end do
  end function deperma_score

  function destandardize_score(chol,score_standardized,step) result(score)
    real(dp),intent(in)::chol(:,:),score_standardized(:)
    real(dp),intent(in),optional::step
    real(dp),allocatable::score(:)
    real(dp),allocatable::theta(:),tp(:),tm(:),sp(:,:),sm(:,:),vp(:),vm(:)
    real(dp)::h
    integer::j,n,k
    n=size(chol,1); theta=pack_lower(chol,.true.,.false.); k=size(theta)
    allocate(score(k)); score=0.0_dp
    do j=1,k
      tp=theta; tm=theta
      h=sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(theta(j)))
      if(present(step)) h=step
      tp(j)=tp(j)+h; tm(j)=tm(j)-h
      sp=dchol(unpack_lower(tp,n,.true.,.false.))
      sm=dchol(unpack_lower(tm,n,.true.,.false.))
      vp=pack_lower(sp,.true.,.false.); vm=pack_lower(sm,.true.,.false.)
      score(j)=dot_product(score_standardized,(vp-vm)/(2.0_dp*h))
    end do
  end function destandardize_score

end module mvtnorm_triangular
