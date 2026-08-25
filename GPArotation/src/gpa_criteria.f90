module gpa_criteria
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use gpa_kinds, only: dp
  use gpa_linalg, only: inverse_matrix, logabsdet
  implicit none
  private

  type, public :: criterion_options
    real(dp) :: gam = 0.0_dp
    real(dp) :: kappa = 0.0_dp
    real(dp) :: delta = 0.01_dp
    integer :: simplimax_k = -1
    real(dp), allocatable :: target(:,:)
    real(dp), allocatable :: weight(:,:)
  end type criterion_options

  type, public :: criterion_result
    real(dp), allocatable :: gq(:,:)
    real(dp) :: f = 0.0_dp
    character(len=64) :: method = ''
    integer :: info = 0
  end type criterion_result

  public :: evaluate_criterion
  public :: vgq_oblimin, vgq_quartimin, vgq_cf, vgq_target, vgq_pst
  public :: vgq_entropy, vgq_infomax, vgq_mccammon, vgq_geomin
  public :: vgq_simplimax, vgq_bifactor, vgq_bigeomin
  public :: vgq_tandem1, vgq_tandem2, vgq_oblimax, vgq_bentler
  public :: vgq_quartimax, vgq_varimax, vgq_binormamin, vgq_varimin, vgq_lp_wls

contains

  subroutine init_result(l, r)
    real(dp), intent(in) :: l(:,:)
    type(criterion_result), intent(out) :: r
    allocate(r%gq(size(l,1),size(l,2)))
    r%gq=0.0_dp; r%f=0.0_dp; r%method=''; r%info=0
  end subroutine init_result

  subroutine evaluate_criterion(name,l,r,opts)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: l(:,:)
    type(criterion_result), intent(out) :: r
    type(criterion_options), intent(in), optional :: opts
    type(criterion_options) :: o
    if(present(opts)) o=opts
    select case(trim(adjustl(lower(name))))
    case('oblimin'); call vgq_oblimin(l,o%gam,r)
    case('quartimin'); call vgq_quartimin(l,r)
    case('cf'); call vgq_cf(l,o%kappa,r)
    case('target')
      if(.not.allocated(o%target)) then; call init_result(l,r); r%info=-2; return; end if
      call vgq_target(l,o%target,r)
    case('pst')
      if(.not.allocated(o%target) .or. .not.allocated(o%weight)) then
        call init_result(l,r); r%info=-2; return
      end if
      call vgq_pst(l,o%weight,o%target,r)
    case('entropy'); call vgq_entropy(l,r)
    case('infomax'); call vgq_infomax(l,r)
    case('mccammon'); call vgq_mccammon(l,r)
    case('geomin'); call vgq_geomin(l,o%delta,r)
    case('simplimax'); call vgq_simplimax(l,o%simplimax_k,r)
    case('bifactor'); call vgq_bifactor(l,r)
    case('bigeomin'); call vgq_bigeomin(l,o%delta,r)
    case('tandemi','tandem1'); call vgq_tandem1(l,r)
    case('tandemii','tandem2'); call vgq_tandem2(l,r)
    case('oblimax'); call vgq_oblimax(l,r)
    case('bentler'); call vgq_bentler(l,r)
    case('quartimax'); call vgq_quartimax(l,r)
    case('varimax'); call vgq_varimax(l,r)
    case('binormamin'); call vgq_binormamin(l,r)
    case('varimin'); call vgq_varimin(l,r)
    case('lp.wls','lp_wls')
      if(.not.allocated(o%weight)) then; call init_result(l,r); r%info=-2; return; end if
      call vgq_lp_wls(l,o%weight,r)
    case default
      call init_result(l,r); r%info=-1; r%method='unknown criterion'
    end select
  end subroutine evaluate_criterion

  subroutine vgq_oblimin(l,gam,r)
    real(dp),intent(in)::l(:,:),gam
    type(criterion_result),intent(out)::r
    real(dp)::x(size(l,1),size(l,2)),cs(size(l,2))
    integer::i,j
    call init_result(l,r)
    do i=1,size(l,1)
      x(i,:)=sum(l(i,:)**2)-l(i,:)**2
    end do
    if(abs(gam)>tiny(1.0_dp)) then
      cs=sum(x,dim=1)
      do j=1,size(l,2); x(:,j)=x(:,j)-(gam/real(size(l,1),dp))*cs(j); end do
    end if
    r%gq=l*x; r%f=sum(l*l*x)/4.0_dp; r%method='Oblimin'
  end subroutine vgq_oblimin

  subroutine vgq_quartimin(l,r)
    real(dp),intent(in)::l(:,:)
    type(criterion_result),intent(out)::r
    real(dp)::l2(size(l,1),size(l,2)),x(size(l,1),size(l,2))
    integer::i
    call init_result(l,r); l2=l*l
    do i=1,size(l,1); x(i,:)=sum(l2(i,:))-l2(i,:); end do
    r%gq=l*x; r%f=sum(l2*x)/4.0_dp; r%method='Quartimin'
  end subroutine vgq_quartimin

  subroutine vgq_cf(l,kappa,r)
    real(dp),intent(in)::l(:,:),kappa
    type(criterion_result),intent(out)::r
    real(dp)::l2(size(l,1),size(l,2)),a(size(l,1),size(l,2)),b(size(l,1),size(l,2))
    integer::i,j
    call init_result(l,r); l2=l*l
    do i=1,size(l,1); a(i,:)=sum(l2(i,:))-l2(i,:); end do
    do j=1,size(l,2); b(:,j)=sum(l2(:,j))-l2(:,j); end do
    r%gq=(1.0_dp-kappa)*l*a+kappa*l*b
    r%f=((1.0_dp-kappa)*sum(l2*a)+kappa*sum(l2*b))/4.0_dp
    r%method='Crawford-Ferguson'
  end subroutine vgq_cf

  subroutine vgq_target(l,target,r)
    real(dp),intent(in)::l(:,:),target(:,:)
    type(criterion_result),intent(out)::r
    integer::i,j
    real(dp)::d
    call init_result(l,r); r%method='Target rotation'
    if(any(shape(l)/=shape(target))) then; r%info=-2; return; end if
    do j=1,size(l,2); do i=1,size(l,1)
      if(ieee_is_nan(target(i,j))) then
        r%gq(i,j)=0.0_dp
      else
        d=l(i,j)-target(i,j); r%gq(i,j)=2.0_dp*d; r%f=r%f+d*d
      end if
    end do; end do
  end subroutine vgq_target

  subroutine vgq_pst(l,w,target,r)
    real(dp),intent(in)::l(:,:),w(:,:),target(:,:)
    type(criterion_result),intent(out)::r
    real(dp)::d(size(l,1),size(l,2))
    call init_result(l,r); r%method='Partially specified target'
    if(any(shape(l)/=shape(w)) .or. any(shape(l)/=shape(target))) then; r%info=-2; return; end if
    d=w*l-w*target; r%gq=2.0_dp*d; r%f=sum(d*d)
  end subroutine vgq_pst

  subroutine vgq_entropy(l,r)
    real(dp),intent(in)::l(:,:)
    type(criterion_result),intent(out)::r
    real(dp)::l2(size(l,1),size(l,2)),lg(size(l,1),size(l,2))
    call init_result(l,r); l2=l*l; lg=log(l2+epsilon(1.0_dp))
    r%gq=-(l*lg+l); r%f=-sum(l2*lg)/2.0_dp; r%method='Minimum entropy'
  end subroutine vgq_entropy

  subroutine vgq_infomax(l,r)
    real(dp),intent(in)::l(:,:)
    type(criterion_result),intent(out)::r
    real(dp)::s2(size(l,1),size(l,2)),e(size(l,1),size(l,2)),h(size(l,1),size(l,2))
    real(dp)::g0(size(l,1),size(l,2)),g1(size(l,1),size(l,2)),g2(size(l,1),size(l,2))
    real(dp)::e1(size(l,1)),e2(size(l,2)),h1(size(l,1)),h2(size(l,2)),ss,q0,q1,q2,a0,a1,a2
    integer::i,j
    call init_result(l,r); s2=l*l; ss=sum(s2)
    if(ss<=0.0_dp) then; r%info=1; return; end if
    e=s2/ss; e1=sum(s2,dim=2)/ss; e2=sum(s2,dim=1)/ss
    q0=-sum(e*log(e+epsilon(1.0_dp)))
    q1=-sum(e1*log(e1+epsilon(1.0_dp)))
    q2=-sum(e2*log(e2+epsilon(1.0_dp)))
    h=-(log(e+epsilon(1.0_dp))+1.0_dp)
    h1=-(log(e1+epsilon(1.0_dp))+1.0_dp)
    h2=-(log(e2+epsilon(1.0_dp))+1.0_dp)
    a0=sum(s2*h)/ss; a1=sum(sum(s2,dim=2)*h1)/ss; a2=sum(sum(s2,dim=1)*h2)/ss
    g0=(h-a0)/ss
    do j=1,size(l,2); do i=1,size(l,1)
      g1(i,j)=(h1(i)-a1)/ss; g2(i,j)=(h2(j)-a2)/ss
    end do; end do
    r%gq=2.0_dp*l*(g0-g1-g2); r%f=log(real(size(l,2),dp))+q0-q1-q2; r%method='Infomax'
  end subroutine vgq_infomax

  subroutine vgq_mccammon(l,r)
    real(dp),intent(in)::l(:,:)
    type(criterion_result),intent(out)::r
    real(dp)::s(size(l,1),size(l,2)),p(size(l,1),size(l,2)),logp(size(l,1),size(l,2))
    real(dp)::h1(size(l,1),size(l,2)),g1(size(l,1),size(l,2)),g2(size(l,1),size(l,2))
    real(dp)::cs(size(l,2)),p2(size(l,2)),lp2(size(l,2)),hh2(size(l,2)),al1(size(l,2))
    real(dp)::ts,q1,q2,al2
    integer::i,j
    call init_result(l,r); s=l*l; cs=sum(s,dim=1); ts=sum(s)
    if(minval(cs)<=0.0_dp .or. ts<=0.0_dp) then; r%info=1; return; end if
    do j=1,size(l,2); p(:,j)=s(:,j)/cs(j); end do
    where(p>0.0_dp); logp=log(p); elsewhere; logp=0.0_dp; end where
    q1=-sum(p*logp); p2=cs/ts
    where(p2>0.0_dp); lp2=log(p2); elsewhere; lp2=0.0_dp; end where
    q2=-sum(p2*lp2)
    if(q1<=0.0_dp .or. q2<=0.0_dp) then; r%info=2; return; end if
    h1=-(logp+1.0_dp); al1=sum(p*h1,dim=1)
    do j=1,size(l,2); g1(:,j)=(h1(:,j)-al1(j))/cs(j); end do
    hh2=-(lp2+1.0_dp); al2=sum(p2*hh2)
    do j=1,size(l,2); do i=1,size(l,1); g2(i,j)=(hh2(j)-al2)/ts; end do; end do
    r%gq=2.0_dp*l*(g1/q1-g2/q2); r%f=log(q1)-log(q2); r%method='McCammon entropy'
  end subroutine vgq_mccammon

  subroutine vgq_geomin(l,delta,r)
    real(dp),intent(in)::l(:,:),delta
    type(criterion_result),intent(out)::r
    real(dp)::l2(size(l,1),size(l,2)),pro(size(l,1))
    integer::i,j,k
    call init_result(l,r); k=size(l,2); l2=l*l+delta
    do i=1,size(l,1); pro(i)=exp(sum(log(l2(i,:)))/real(k,dp)); end do
    do j=1,k; r%gq(:,j)=(2.0_dp/real(k,dp))*(l(:,j)/l2(:,j))*pro; end do
    r%f=sum(pro); r%method='Geomin'
  end subroutine vgq_geomin

  subroutine vgq_simplimax(l,kin,r)
    real(dp),intent(in)::l(:,:)
    integer,intent(in)::kin
    type(criterion_result),intent(out)::r
    real(dp),allocatable::vals(:)
    logical,allocatable::pick(:)
    integer::n,k,j,idx
    call init_result(l,r); n=size(l); k=kin
    if(k<0) k=size(l,1)*(size(l,2)-1); k=max(0,min(k,n))
    allocate(vals(n),pick(n)); vals=reshape(l*l,[n]); pick=.false.
    do j=1,k
      idx=minloc(vals,dim=1,mask=.not.pick); pick(idx)=.true.
    end do
    r%gq=reshape(merge(2.0_dp*reshape(l,[n]),0.0_dp,pick),shape(l))
    r%f=sum(vals,mask=pick); r%method='Simplimax'
  end subroutine vgq_simplimax

  subroutine vgq_bifactor(l,r)
    real(dp),intent(in)::l(:,:)
    type(criterion_result),intent(out)::r
    real(dp)::lt2(size(l,1),max(1,size(l,2)-1)),rows(size(l,1))
    integer::j
    call init_result(l,r); r%method='Bifactor Biquartimin'
    if(size(l,2)<2) then; r%info=-1; return; end if
    lt2=l(:,2:)**2; rows=sum(lt2(:,1:size(l,2)-1),dim=2)
    do j=2,size(l,2); r%gq(:,j)=4.0_dp*l(:,j)*(rows-l(:,j)**2); end do
    r%f=sum(rows*rows)-sum(lt2(:,1:size(l,2)-1)**2)
  end subroutine vgq_bifactor

  subroutine vgq_bigeomin(l,delta,r)
    real(dp),intent(in)::l(:,:),delta
    type(criterion_result),intent(out)::r
    type(criterion_result)::sub
    call init_result(l,r); r%method='Bi-Geomin'
    if(size(l,2)<2) then; r%info=-1; return; end if
    call vgq_geomin(l(:,2:),delta,sub); r%gq(:,2:)=sub%gq; r%f=sub%f; r%info=sub%info
  end subroutine vgq_bigeomin

  subroutine vgq_tandem1(l,r)
    real(dp),intent(in)::l(:,:)
    type(criterion_result),intent(out)::r
    real(dp)::l2(size(l,1),size(l,2)),ll(size(l,1),size(l,1)),ll2(size(l,1),size(l,1))
    real(dp)::mp(size(l,1),size(l,2)),tc(size(l,1),size(l,1))
    call init_result(l,r); l2=l*l; ll=matmul(l,transpose(l)); ll2=ll*ll
    mp=matmul(ll2,l2); tc=matmul(l2,transpose(l2))
    r%f=-sum(l2*mp); r%gq=-(4.0_dp*l*mp+4.0_dp*matmul(ll*tc,l)); r%method='Tandem I'
  end subroutine vgq_tandem1

  subroutine vgq_tandem2(l,r)
    real(dp),intent(in)::l(:,:)
    type(criterion_result),intent(out)::r
    real(dp)::l2(size(l,1),size(l,2)),ll(size(l,1),size(l,1)),ill2(size(l,1),size(l,1))
    real(dp)::mp(size(l,1),size(l,2)),tc(size(l,1),size(l,1))
    call init_result(l,r); l2=l*l; ll=matmul(l,transpose(l)); ill2=1.0_dp-ll*ll
    mp=matmul(ill2,l2); tc=matmul(l2,transpose(l2))
    r%f=sum(l2*mp); r%gq=4.0_dp*l*mp-4.0_dp*matmul(ll*tc,l); r%method='Tandem II'
  end subroutine vgq_tandem2

  subroutine vgq_oblimax(l,r)
    real(dp),intent(in)::l(:,:)
    type(criterion_result),intent(out)::r
    real(dp)::s2,s4
    call init_result(l,r); s2=sum(l*l); s4=sum(l**4)
    if(s2<=0.0_dp .or. s4<=0.0_dp) then; r%info=1; return; end if
    r%gq=-(4.0_dp*l**3/s4-4.0_dp*l/s2); r%f=-(log(s4)-2.0_dp*log(s2)); r%method='Oblimax'
  end subroutine vgq_oblimax

  subroutine vgq_bentler(l,r)
    real(dp),intent(in)::l(:,:)
    type(criterion_result),intent(out)::r
    real(dp)::l2(size(l,1),size(l,2)),m(size(l,2),size(l,2)),mi(size(l,2),size(l,2))
    real(dp)::d(size(l,2)),dinv(size(l,2),size(l,2)),ldm,ldd
    integer::j,info
    call init_result(l,r); l2=l*l; m=matmul(transpose(l2),l2); d=0.0_dp
    do j=1,size(l,2); d(j)=m(j,j); end do
    if(minval(d)<=0.0_dp) then; r%info=1; return; end if
    ldm=logabsdet(m,info=info); if(info/=0) then; r%info=info; return; end if
    ldd=sum(log(d)); call inverse_matrix(m,mi,info); if(info/=0) then; r%info=info; return; end if
    dinv=0.0_dp; do j=1,size(l,2); dinv(j,j)=1.0_dp/d(j); end do
    r%f=-(ldm-ldd)/4.0_dp; r%gq=-l*matmul(l2,mi-dinv); r%method='Bentler'
  end subroutine vgq_bentler

  subroutine vgq_quartimax(l,r)
    real(dp),intent(in)::l(:,:)
    type(criterion_result),intent(out)::r
    call init_result(l,r); r%gq=-l**3; r%f=-sum(l**4)/4.0_dp; r%method='Quartimax'
  end subroutine vgq_quartimax

  subroutine vgq_varimax(l,r)
    real(dp),intent(in)::l(:,:)
    type(criterion_result),intent(out)::r
    real(dp)::ql(size(l,1),size(l,2)),cm
    integer::j
    call init_result(l,r)
    do j=1,size(l,2); cm=sum(l(:,j)**2)/real(size(l,1),dp); ql(:,j)=l(:,j)**2-cm; end do
    r%gq=-l*ql; r%f=-sum(ql*ql)/4.0_dp; r%method='Varimax'
  end subroutine vgq_varimax

  subroutine vgq_binormamin(l,r)
    real(dp),intent(in)::l(:,:)
    type(criterion_result),intent(out)::r
    integer::p,m,i,j
    real(dp)::l2(size(l,1),size(l,2)),u(size(l,2)),invu(size(l,2))
    real(dp)::nmat(size(l,2),size(l,2)),c(size(l,2),size(l,2))
    real(dp)::l2u(size(l,1),size(l,2)),sumrow(size(l,1)),t1(size(l,1),size(l,2))
    real(dp)::noveru(size(l,2),size(l,2)),sumk(size(l,2)),t2(size(l,2))
    call init_result(l,r); p=size(l,1); m=size(l,2); l2=l*l; u=sum(l2,dim=1)
    if(minval(u)<=0.0_dp) then; r%info=1; return; end if
    invu=1.0_dp/u; nmat=matmul(transpose(l2),l2)
    do j=1,m; do i=1,m; c(i,j)=nmat(i,j)/(u(i)*u(j)); end do; end do
    r%f=sum(c)-sum([(c(i,i),i=1,m)])
    do j=1,m; l2u(:,j)=l2(:,j)*invu(j); end do; sumrow=sum(l2u,dim=2)
    do j=1,m; t1(:,j)=(sumrow-l2u(:,j))*invu(j); noveru(:,j)=nmat(:,j)*invu(j); end do
    do i=1,m; sumk(i)=sum(noveru(i,:))-noveru(i,i); t2(i)=sumk(i)*invu(i)**2; end do
    do j=1,m; r%gq(:,j)=4.0_dp*l(:,j)*(t1(:,j)-t2(j)); end do
    r%method='Binormamin'
  end subroutine vgq_binormamin

  subroutine vgq_varimin(l,r)
    real(dp),intent(in)::l(:,:)
    type(criterion_result),intent(out)::r
    real(dp)::ql(size(l,1),size(l,2)),cm
    integer::j
    call init_result(l,r)
    do j=1,size(l,2); cm=sum(l(:,j)**2)/real(size(l,1),dp); ql(:,j)=l(:,j)**2-cm; end do
    r%gq=l*ql; r%f=sum(ql*ql)/4.0_dp; r%method='Varimin'
  end subroutine vgq_varimin

  subroutine vgq_lp_wls(l,w,r)
    real(dp),intent(in)::l(:,:),w(:,:)
    type(criterion_result),intent(out)::r
    call init_result(l,r); r%method='Weighted least squares for Lp rotation'
    if(any(shape(l)/=shape(w))) then; r%info=-2; return; end if
    r%gq=2.0_dp*w*l/real(size(l,1),dp); r%f=sum(w*l*l)/real(size(l,1),dp)
  end subroutine vgq_lp_wls

  pure function lower(s) result(t)
    character(len=*),intent(in)::s
    character(len=len(s))::t
    integer::i,c
    do i=1,len(s)
      c=iachar(s(i:i)); if(c>=65 .and. c<=90) then; t(i:i)=achar(c+32); else; t(i:i)=s(i:i); end if
    end do
  end function lower

end module gpa_criteria
