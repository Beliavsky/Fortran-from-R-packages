! SPDX-License-Identifier: GPL-2.0-or-later
module compositions_geometry
  use compositions_kinds, only: dp
  use compositions_linalg, only: pseudoinverse
  implicit none
  private
  public :: closure, closure_rows, perturb, power_comp
  public :: clr, clr_rows, clr_inv, clr_inv_rows
  public :: ilr_base, ilr, ilr_rows, ilr_inv, ilr_inv_rows
  public :: alr, alr_rows, alr_inv, alr_inv_rows
  public :: apt, apt_rows, apt_inv, apt_inv_rows
  public :: cpt, cpt_rows, cpt_inv, cpt_inv_rows, ipt, ipt_rows, ipt_inv, ipt_inv_rows
  public :: ilt, ilt_rows, ilt_inv, ilt_inv_rows
  public :: clrvar_to_ilr, ilrvar_to_clr, clrvar_to_variation, variation_to_clrvar
  public :: geometric_mean, geometric_mean_rows, geometric_mean_cols
  public :: acomp_mean, variation_matrix, pairwise_logratios
  public :: build_ilr_base, balance_coordinate, balance_basis_balanced
  public :: missing_projector_acomp, sum_missing_projector_acomp

contains
  function closure(x,total) result(y)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: total
    real(dp) :: y(size(x)), s,t
    if(any(x<0.0_dp)) error stop 'closure: negative values are not valid'
    s=sum(x)
    if(s<=0.0_dp) error stop 'closure: sum must be positive'
    t=1.0_dp; if(present(total)) t=total
    y=x*t/s
  end function closure

  function closure_rows(x,total) result(y)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(in), optional :: total
    real(dp) :: y(size(x,1),size(x,2))
    integer :: i
    do i=1,size(x,1)
      if(present(total)) then
        y(i,:)=closure(x(i,:),total)
      else
        y(i,:)=closure(x(i,:))
      end if
    end do
  end function closure_rows

  function perturb(x,y) result(z)
    real(dp), intent(in) :: x(:),y(:)
    real(dp) :: z(size(x))
    if(size(y)/=size(x)) error stop 'perturb: dimension mismatch'
    z=closure(x*y)
  end function perturb

  function power_comp(x,a) result(z)
    real(dp), intent(in) :: x(:),a
    real(dp) :: z(size(x))
    if(any(x<=0.0_dp)) error stop 'power_comp: all parts must be positive'
    z=closure(x**a)
  end function power_comp

  function clr(x) result(z)
    real(dp), intent(in) :: x(:)
    real(dp) :: z(size(x)), lx(size(x)),m
    if(any(x<=0.0_dp)) error stop 'clr: all parts must be positive'
    lx=log(x); m=sum(lx)/real(size(x),dp); z=lx-m
  end function clr

  function clr_rows(x) result(z)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: z(size(x,1),size(x,2))
    integer :: i
    do i=1,size(x,1); z(i,:)=clr(x(i,:)); end do
  end function clr_rows

  function clr_inv(z) result(x)
    real(dp), intent(in) :: z(:)
    real(dp) :: x(size(z)),m
    m=maxval(z); x=closure(exp(z-m))
  end function clr_inv

  function clr_inv_rows(z) result(x)
    real(dp), intent(in) :: z(:,:)
    real(dp) :: x(size(z,1),size(z,2))
    integer :: i
    do i=1,size(z,1); x(i,:)=clr_inv(z(i,:)); end do
  end function clr_inv_rows

  function ilr_base(d) result(v)
    integer, intent(in) :: d
    real(dp), allocatable :: v(:,:)
    integer :: j
    if(d<1) error stop 'ilr_base: d must be positive'
    if(d==1) then
      allocate(v(1,0)); return
    end if
    allocate(v(d,d-1)); v=0.0_dp
    ! Normalized Helmert contrasts, equivalent to compositions basic ilr base.
    do j=1,d-1
      v(1:j,j)=1.0_dp/sqrt(real(j*(j+1),dp))
      v(j+1,j)=-real(j,dp)/sqrt(real(j*(j+1),dp))
    end do
  end function ilr_base

  function ilr(x,v) result(z)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: v(:,:)
    real(dp), allocatable :: z(:),vb(:,:)
    if(present(v)) then
      if(size(v,1)/=size(x)) error stop 'ilr: basis row mismatch'
      allocate(z(size(v,2))); z=matmul(clr(x),v)
    else
      vb=ilr_base(size(x)); allocate(z(size(x)-1)); z=matmul(clr(x),vb)
    end if
  end function ilr

  function ilr_rows(x,v) result(z)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(in), optional :: v(:,:)
    real(dp), allocatable :: z(:,:),vb(:,:)
    if(present(v)) then
      allocate(z(size(x,1),size(v,2))); z=matmul(clr_rows(x),v)
    else
      vb=ilr_base(size(x,2)); allocate(z(size(x,1),size(x,2)-1)); z=matmul(clr_rows(x),vb)
    end if
  end function ilr_rows

  function ilr_inv(z,v) result(x)
    real(dp), intent(in) :: z(:)
    real(dp), intent(in), optional :: v(:,:)
    real(dp), allocatable :: x(:),vb(:,:),c(:)
    if(present(v)) then
      allocate(c(size(v,1))); c=matmul(v,z)
    else
      vb=ilr_base(size(z)+1); allocate(c(size(z)+1)); c=matmul(vb,z)
    end if
    x=clr_inv(c)
  end function ilr_inv

  function ilr_inv_rows(z,v) result(x)
    real(dp), intent(in) :: z(:,:)
    real(dp), intent(in), optional :: v(:,:)
    real(dp), allocatable :: x(:,:),vb(:,:),c(:,:)
    integer :: i
    if(present(v)) then
      allocate(c(size(z,1),size(v,1))); c=matmul(z,transpose(v))
    else
      vb=ilr_base(size(z,2)+1); allocate(c(size(z,1),size(z,2)+1)); c=matmul(z,transpose(vb))
    end if
    allocate(x(size(c,1),size(c,2)))
    do i=1,size(c,1); x(i,:)=clr_inv(c(i,:)); end do
  end function ilr_inv_rows

  function alr(x,ivar) result(z)
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: ivar
    real(dp), allocatable :: z(:)
    integer :: d,j,k,iv
    d=size(x); iv=d; if(present(ivar)) iv=ivar
    if(iv<1 .or. iv>d) error stop 'alr: invalid denominator index'
    if(any(x<=0.0_dp)) error stop 'alr: all parts must be positive'
    allocate(z(d-1)); k=0
    do j=1,d
      if(j/=iv) then; k=k+1; z(k)=log(x(j)/x(iv)); end if
    end do
  end function alr

  function alr_rows(x,ivar) result(z)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in), optional :: ivar
    real(dp), allocatable :: z(:,:)
    integer :: i,iv
    iv=size(x,2); if(present(ivar)) iv=ivar
    allocate(z(size(x,1),size(x,2)-1))
    do i=1,size(x,1); z(i,:)=alr(x(i,:),iv); end do
  end function alr_rows

  function alr_inv(z,ivar,d) result(x)
    real(dp), intent(in) :: z(:)
    integer, intent(in), optional :: ivar,d
    real(dp), allocatable :: x(:)
    integer :: dd,iv,j,k
    dd=size(z)+1; if(present(d)) dd=d
    if(dd/=size(z)+1) error stop 'alr_inv: d must equal size(z)+1'
    iv=dd; if(present(ivar)) iv=ivar
    allocate(x(dd)); x=1.0_dp; k=0
    do j=1,dd
      if(j/=iv) then; k=k+1; x(j)=exp(z(k)); end if
    end do
    x=closure(x)
  end function alr_inv

  function alr_inv_rows(z,ivar,d) result(x)
    real(dp), intent(in) :: z(:,:)
    integer, intent(in), optional :: ivar,d
    real(dp), allocatable :: x(:,:)
    integer :: i,dd,iv
    dd=size(z,2)+1; if(present(d)) dd=d
    iv=dd; if(present(ivar)) iv=ivar
    allocate(x(size(z,1),dd))
    do i=1,size(z,1); x(i,:)=alr_inv(z(i,:),iv,dd); end do
  end function alr_inv_rows

  function apt(x) result(z)
    real(dp), intent(in) :: x(:)
    real(dp) :: z(size(x)-1),c(size(x))
    if(size(x)<2) error stop 'apt: at least two parts required'
    c=closure(x); z=c(1:size(x)-1)
  end function apt

  function apt_rows(x) result(z)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: z(size(x,1),size(x,2)-1)
    integer :: i
    do i=1,size(x,1); z(i,:)=apt(x(i,:)); end do
  end function apt_rows

  function apt_inv(z) result(x)
    real(dp), intent(in) :: z(:)
    real(dp) :: x(size(z)+1)
    x(1:size(z))=z; x(size(x))=1.0_dp-sum(z)
  end function apt_inv

  function apt_inv_rows(z) result(x)
    real(dp), intent(in) :: z(:,:)
    real(dp) :: x(size(z,1),size(z,2)+1)
    integer :: i
    do i=1,size(z,1); x(i,:)=apt_inv(z(i,:)); end do
  end function apt_inv_rows

  function cpt(x) result(z)
    real(dp), intent(in) :: x(:)
    real(dp) :: z(size(x))
    z=closure(x)-1.0_dp/real(size(x),dp)
  end function cpt

  function cpt_rows(x) result(z)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: z(size(x,1),size(x,2))
    integer :: i
    do i=1,size(x,1); z(i,:)=cpt(x(i,:)); end do
  end function cpt_rows

  function cpt_inv(z) result(x)
    real(dp), intent(in) :: z(:)
    real(dp) :: x(size(z))
    x=z+1.0_dp/real(size(z),dp)
    if(any(x<0.0_dp)) error stop 'cpt_inv: point outside simplex'
    x=closure(x)
  end function cpt_inv

  function cpt_inv_rows(z) result(x)
    real(dp), intent(in) :: z(:,:)
    real(dp) :: x(size(z,1),size(z,2))
    integer :: i
    do i=1,size(z,1); x(i,:)=cpt_inv(z(i,:)); end do
  end function cpt_inv_rows

  function ipt(x,v) result(z)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: v(:,:)
    real(dp), allocatable :: z(:),vb(:,:)
    if(present(v)) then; allocate(z(size(v,2))); z=matmul(cpt(x),v)
    else; vb=ilr_base(size(x)); allocate(z(size(x)-1)); z=matmul(cpt(x),vb); end if
  end function ipt

  function ipt_rows(x,v) result(z)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(in), optional :: v(:,:)
    real(dp), allocatable :: z(:,:),vb(:,:)
    if(present(v)) then; allocate(z(size(x,1),size(v,2))); z=matmul(cpt_rows(x),v)
    else; vb=ilr_base(size(x,2)); allocate(z(size(x,1),size(x,2)-1)); z=matmul(cpt_rows(x),vb); end if
  end function ipt_rows

  function ipt_inv(z,v) result(x)
    real(dp), intent(in) :: z(:)
    real(dp), intent(in), optional :: v(:,:)
    real(dp), allocatable :: x(:),vb(:,:),c(:)
    if(present(v)) then; allocate(c(size(v,1))); c=matmul(v,z)
    else; vb=ilr_base(size(z)+1); allocate(c(size(z)+1)); c=matmul(vb,z); end if
    x=cpt_inv(c)
  end function ipt_inv

  function ipt_inv_rows(z,v) result(x)
    real(dp), intent(in) :: z(:,:)
    real(dp), intent(in), optional :: v(:,:)
    real(dp), allocatable :: x(:,:),vb(:,:),c(:,:)
    integer :: i
    if(present(v)) then; allocate(c(size(z,1),size(v,1))); c=matmul(z,transpose(v))
    else; vb=ilr_base(size(z,2)+1); allocate(c(size(z,1),size(z,2)+1)); c=matmul(z,transpose(vb)); end if
    allocate(x(size(c,1),size(c,2)))
    do i=1,size(c,1); x(i,:)=cpt_inv(c(i,:)); end do
  end function ipt_inv_rows

  function ilt(x) result(z)
    real(dp), intent(in) :: x(:)
    real(dp) :: z(size(x))
    if(any(x<=0.0_dp)) error stop 'ilt: all components must be positive'
    z=log(x)
  end function ilt

  function ilt_rows(x) result(z)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: z(size(x,1),size(x,2))
    if(any(x<=0.0_dp)) error stop 'ilt_rows: all components must be positive'
    z=log(x)
  end function ilt_rows

  function ilt_inv(z) result(x)
    real(dp), intent(in) :: z(:)
    real(dp) :: x(size(z))
    x=exp(z)
  end function ilt_inv

  function ilt_inv_rows(z) result(x)
    real(dp), intent(in) :: z(:,:)
    real(dp) :: x(size(z,1),size(z,2))
    x=exp(z)
  end function ilt_inv_rows

  function clrvar_to_ilr(sigma,v) result(s)
    real(dp), intent(in) :: sigma(:,:)
    real(dp), intent(in), optional :: v(:,:)
    real(dp), allocatable :: s(:,:),vb(:,:)
    if(present(v)) then; s=matmul(transpose(v),matmul(sigma,v))
    else; vb=ilr_base(size(sigma,1)); s=matmul(transpose(vb),matmul(sigma,vb)); end if
  end function clrvar_to_ilr

  function ilrvar_to_clr(sigma,v) result(s)
    real(dp), intent(in) :: sigma(:,:)
    real(dp), intent(in), optional :: v(:,:)
    real(dp), allocatable :: s(:,:),vb(:,:)
    if(present(v)) then; s=matmul(v,matmul(sigma,transpose(v)))
    else; vb=ilr_base(size(sigma,1)+1); s=matmul(vb,matmul(sigma,transpose(vb))); end if
  end function ilrvar_to_clr

  function clrvar_to_variation(sigma) result(t)
    real(dp), intent(in) :: sigma(:,:)
    real(dp) :: t(size(sigma,1),size(sigma,2))
    integer :: i,j,d
    d=size(sigma,1)
    do j=1,d; do i=1,d
      t(i,j)=sigma(i,i)+sigma(j,j)-2.0_dp*sigma(i,j)
    end do; end do
  end function clrvar_to_variation

  function variation_to_clrvar(t) result(sigma)
    real(dp), intent(in) :: t(:,:)
    real(dp) :: sigma(size(t,1),size(t,2)),h(size(t,1),size(t,2))
    integer :: d,i
    d=size(t,1); h=-1.0_dp/real(d,dp)
    do i=1,d; h(i,i)=h(i,i)+1.0_dp; end do
    sigma=-0.5_dp*matmul(h,matmul(t,h))
  end function variation_to_clrvar

  real(dp) function geometric_mean(x) result(g)
    real(dp), intent(in) :: x(:)
    if(any(x<=0.0_dp)) error stop 'geometric_mean: values must be positive'
    g=exp(sum(log(x))/real(size(x),dp))
  end function geometric_mean

  function geometric_mean_rows(x) result(g)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: g(size(x,1))
    integer :: i
    do i=1,size(x,1); g(i)=geometric_mean(x(i,:)); end do
  end function geometric_mean_rows

  function geometric_mean_cols(x) result(g)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: g(size(x,2))
    integer :: j
    do j=1,size(x,2); g(j)=geometric_mean(x(:,j)); end do
  end function geometric_mean_cols

  function acomp_mean(x) result(mu)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: mu(size(x,2)),logs(size(x,2))
    if(any(x<=0.0_dp)) error stop 'acomp_mean: values must be positive'
    logs=sum(log(x),dim=1)/real(size(x,1),dp)
    mu=closure(exp(logs-maxval(logs)))
  end function acomp_mean

  function variation_matrix(x) result(t)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: t(size(x,2),size(x,2)),lr(size(x,1)),m
    integer :: i,j,n,d
    n=size(x,1); d=size(x,2)
    if(any(x<=0.0_dp)) error stop 'variation_matrix: values must be positive'
    t=0.0_dp
    do j=1,d; do i=1,d
      lr=log(x(:,i)/x(:,j)); m=sum(lr)/real(n,dp)
      if(n>1) t(i,j)=sum((lr-m)**2)/real(n-1,dp)
    end do; end do
  end function variation_matrix

  function pairwise_logratios(x) result(y)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: y(:,:)
    integer :: d,p,i,j,k
    d=size(x,2); p=d*(d-1)/2
    allocate(y(size(x,1),p)); k=0
    do i=1,d-1; do j=i+1,d
      k=k+1; y(:,k)=log(x(:,j)/x(:,i))
    end do; end do
  end function pairwise_logratios

  function build_ilr_base(signary) result(v)
    integer, intent(in) :: signary(:,:)
    real(dp), allocatable :: v(:,:)
    integer :: d,k,j,np,nn
    real(dp) :: normv
    d=size(signary,1); k=size(signary,2)
    allocate(v(d,k)); v=0.0_dp
    do j=1,k
      np=count(signary(:,j)>0); nn=count(signary(:,j)<0)
      if(np==0 .or. nn==0) cycle
      where(signary(:,j)>0) v(:,j)=real(nn,dp)
      where(signary(:,j)<0) v(:,j)=-real(np,dp)
      normv=sqrt(sum(v(:,j)**2)); if(normv>0.0_dp) v(:,j)=v(:,j)/normv
    end do
  end function build_ilr_base

  real(dp) function balance_coordinate(x,signs) result(b)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: signs(:)
    integer :: r,s
    r=count(signs>0); s=count(signs<0)
    if(r==0 .or. s==0) error stop 'balance_coordinate: both groups required'
    b=sqrt(real(r*s,dp)/real(r+s,dp))* &
      (sum(log(x),mask=signs>0)/real(r,dp)-sum(log(x),mask=signs<0)/real(s,dp))
  end function balance_coordinate

  function balance_basis_balanced(d) result(v)
    integer, intent(in) :: d
    real(dp), allocatable :: v(:,:)
    integer, allocatable :: sgn(:,:),left(:),right(:),queue_lo(:),queue_hi(:)
    integer :: ncol,head,tail,lo,hi,mid,j
    if(d<2) then; allocate(v(d,0)); return; end if
    allocate(sgn(d,d-1)); sgn=0
    allocate(queue_lo(d),queue_hi(d)); head=1; tail=1; queue_lo(1)=1; queue_hi(1)=d; ncol=0
    do while(head<=tail .and. ncol<d-1)
      lo=queue_lo(head); hi=queue_hi(head); head=head+1
      if(lo>=hi) cycle
      mid=(lo+hi)/2; ncol=ncol+1
      sgn(lo:mid,ncol)=1; sgn(mid+1:hi,ncol)=-1
      if(mid>lo) then; tail=tail+1; queue_lo(tail)=lo; queue_hi(tail)=mid; end if
      if(hi>mid+1) then; tail=tail+1; queue_lo(tail)=mid+1; queue_hi(tail)=hi; end if
    end do
    v=build_ilr_base(sgn)
  end function balance_basis_balanced

  function missing_projector_acomp(has) result(p)
    logical, intent(in) :: has(:)
    real(dp) :: p(size(has),size(has))
    integer :: i,j,nobs
    nobs=count(has); p=0.0_dp
    if(nobs==0) return
    do i=1,size(has)
      if(has(i)) p(i,i)=1.0_dp
    end do
    do j=1,size(has); do i=1,size(has)
      if(has(i).and.has(j)) p(i,j)=p(i,j)-1.0_dp/real(nobs,dp)
    end do; end do
  end function missing_projector_acomp

  function sum_missing_projector_acomp(has) result(p)
    logical, intent(in) :: has(:,:)
    real(dp) :: p(size(has,2),size(has,2))
    integer :: i,j,k,nobs,d
    d=size(has,2); p=0.0_dp
    do k=1,size(has,1)
      nobs=count(has(k,:)); if(nobs==0) cycle
      do i=1,d; if(has(k,i)) p(i,i)=p(i,i)+1.0_dp; end do
      do j=1,d; do i=1,d
        if(has(k,i).and.has(k,j)) p(i,j)=p(i,j)-1.0_dp/real(nobs,dp)
      end do; end do
    end do
  end function sum_missing_projector_acomp
end module compositions_geometry
