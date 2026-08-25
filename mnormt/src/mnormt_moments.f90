module mnormt_moments
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use mnormt_special, only: dp, normal_pdf, normal_cdf
  use mnormt_linalg, only: pd_solve
  use mnormt_core, only: sadmvn_prob, probability_result
  implicit none
  private

  type :: boundary_table
    real(dp), allocatable :: values(:)
    integer, allocatable :: shape(:)
  end type boundary_table

  type, public :: trunc_moment_result
    real(dp) :: probability = 0.0_dp
    integer, allocatable :: shape(:)
    real(dp), allocatable :: raw(:)
    real(dp), allocatable :: mean(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: cum3(:,:,:)
    real(dp), allocatable :: cum4(:,:,:,:)
    real(dp), allocatable :: marginal_skewness(:)
    real(dp), allocatable :: marginal_excess_kurtosis(:)
    real(dp) :: gamma1_mardia = 0.0_dp
    real(dp) :: gamma2_mardia = 0.0_dp
    real(dp) :: beta2_mardia = 0.0_dp
    integer :: status = 0
  end type trunc_moment_result

  type, public :: mardia_result
    real(dp) :: b1 = 0.0_dp, b2 = 0.0_dp
    real(dp) :: g1 = 0.0_dp, g2 = 0.0_dp
    real(dp) :: p_b1 = 0.0_dp, p_b2 = 0.0_dp
    integer :: n = 0
    integer :: status = 0
  end type mardia_result

  public :: recintab, mom_mtruncnorm, mom2cum, raw_moment_at, sample_mardia_measures
contains

  pure integer function flat_index(alpha,shape) result(idx)
    integer, intent(in) :: alpha(:),shape(:)
    integer :: i,stride
    idx=1; stride=1
    do i=1,size(alpha)
      idx=idx+alpha(i)*stride
      stride=stride*shape(i)
    end do
  end function flat_index

  pure subroutine index_to_alpha(idx,shape,alpha)
    integer, intent(in) :: idx,shape(:)
    integer, intent(out) :: alpha(size(shape))
    integer :: i,q
    q=idx-1
    do i=1,size(shape)
      alpha(i)=mod(q,shape(i)); q=q/shape(i)
    end do
  end subroutine index_to_alpha

  pure real(dp) function pow_boundary(x,k) result(v)
    real(dp), intent(in) :: x
    integer, intent(in) :: k
    if(k==0) then
      v=1.0_dp
    else
      v=x**k
    end if
  end function pow_boundary

  recursive subroutine recintab(kappa,a,b,mu,s,m,status)
    integer, intent(in) :: kappa(:)
    real(dp), intent(in) :: a(:),b(:),mu(:),s(:,:)
    real(dp), allocatable, intent(out) :: m(:)
    integer, intent(out), optional :: status
    integer :: n,pk,idx,i,j,k,ii,istat
    integer, allocatable :: shape(:),alpha(:),beta(:),subk(:),subalpha(:),map(:)
    real(dp), allocatable :: suba(:),subb(:),submu(:),subs(:,:),cvec(:)
    type(boundary_table), allocatable :: g(:),h(:)
    type(probability_result) :: pr
    real(dp) :: sd,aa,bb,pdfa,pdfb,term

    n=size(kappa)
    if(present(status)) status=0
    if(size(a)/=n .or. size(b)/=n .or. size(mu)/=n .or. size(s,1)/=n .or. size(s,2)/=n) then
      allocate(m(1)); m=ieee_value(1.0_dp,ieee_quiet_nan); if(present(status)) status=-1; return
    end if
    if(any(kappa<0) .or. any(a>=b)) then
      allocate(m(1)); m=ieee_value(1.0_dp,ieee_quiet_nan); if(present(status)) status=-2; return
    end if
    shape=kappa+1; pk=product(shape); allocate(m(pk)); m=0.0_dp
    if(n==1) then
      sd=sqrt(s(1,1)); aa=(a(1)-mu(1))/sd; bb=(b(1)-mu(1))/sd
      m(1)=normal_cdf(bb)-normal_cdf(aa)
      if(kappa(1)>0) then
        pdfa=merge(sd*normal_pdf(aa),0.0_dp,ieee_is_finite(a(1)))
        pdfb=merge(sd*normal_pdf(bb),0.0_dp,ieee_is_finite(b(1)))
        m(2)=mu(1)*m(1)+pdfa-pdfb
        do i=2,kappa(1)
          if(ieee_is_finite(a(1))) pdfa=pdfa*a(1)
          if(ieee_is_finite(b(1))) pdfb=pdfb*b(1)
          m(i+1)=mu(1)*m(i)+real(i-1,dp)*s(1,1)*m(i-1)+pdfa-pdfb
        end do
      end if
      return
    end if

    pr=sadmvn_prob(a,b,mu,s)
    m(1)=pr%value
    if(pr%status>1) then; if(present(status)) status=pr%status; return; end if

    allocate(g(n),h(n))
    do j=1,n
      allocate(subk(n-1),suba(n-1),subb(n-1),submu(n-1),subs(n-1,n-1),cvec(n-1),map(n-1))
      k=0
      do i=1,n
        if(i/=j) then; k=k+1; map(k)=i; subk(k)=kappa(i); suba(k)=a(i); subb(k)=b(i); cvec(k)=s(i,j); end if
      end do
      do i=1,n-1
        submu(i)=mu(map(i))
        do k=1,n-1
          subs(i,k)=s(map(i),map(k))-cvec(i)*cvec(k)/s(j,j)
        end do
      end do
      g(j)%shape=subk+1; h(j)%shape=subk+1
      if(ieee_is_finite(a(j))) then
        submu=[(mu(map(i))+cvec(i)/s(j,j)*(a(j)-mu(j)),i=1,n-1)]
        call recintab(subk,suba,subb,submu,subs,g(j)%values,istat)
        g(j)%values=g(j)%values*normal_pdf((a(j)-mu(j))/sqrt(s(j,j)))/sqrt(s(j,j))
      else
        allocate(g(j)%values(product(subk+1))); g(j)%values=0.0_dp
      end if
      if(ieee_is_finite(b(j))) then
        submu=[(mu(map(i))+cvec(i)/s(j,j)*(b(j)-mu(j)),i=1,n-1)]
        call recintab(subk,suba,subb,submu,subs,h(j)%values,istat)
        h(j)%values=h(j)%values*normal_pdf((b(j)-mu(j))/sqrt(s(j,j)))/sqrt(s(j,j))
      else
        allocate(h(j)%values(product(subk+1))); h(j)%values=0.0_dp
      end if
      deallocate(subk,suba,subb,submu,subs,cvec,map)
    end do

    allocate(alpha(n),beta(n))
    do idx=2,pk
      call index_to_alpha(idx,shape,alpha)
      i=1; do while(i<=n .and. alpha(i)==0); i=i+1; end do
      if(i>n) cycle
      beta=alpha; beta(i)=beta(i)-1
      m(idx)=mu(i)*m(flat_index(beta,shape))
      do j=1,n
        term=0.0_dp
        if(beta(j)>0) then
          beta(j)=beta(j)-1
          term=term+real(beta(j)+1,dp)*m(flat_index(beta,shape))
          beta(j)=beta(j)+1
        end if
        allocate(subalpha(n-1)); k=0
        do ii=1,n
          if(ii/=j) then; k=k+1; subalpha(k)=beta(ii); end if
        end do
        term=term+pow_boundary(a(j),beta(j))*g(j)%values(flat_index(subalpha,g(j)%shape)) &
                 -pow_boundary(b(j),beta(j))*h(j)%values(flat_index(subalpha,h(j)%shape))
        m(idx)=m(idx)+s(i,j)*term
        deallocate(subalpha)
      end do
    end do
  end subroutine recintab

  pure real(dp) function raw_moment_at(raw,shape,alpha) result(v)
    real(dp), intent(in) :: raw(:)
    integer, intent(in) :: shape(:),alpha(:)
    if(any(alpha<0) .or. any(alpha>=shape)) then
      v=ieee_value(1.0_dp,ieee_quiet_nan)
    else
      v=raw(flat_index(alpha,shape))
    end if
  end function raw_moment_at

  function mom_mtruncnorm(kappa,mean,varcov,lower,upper) result(out)
    integer, intent(in) :: kappa(:)
    real(dp), intent(in) :: mean(:),varcov(:,:),lower(:),upper(:)
    type(trunc_moment_result) :: out
    real(dp), allocatable :: m(:),inv(:,:)
    integer, allocatable :: alpha(:)
    real(dp) :: m2ij,m2kl,m3ijk,m3jkl,m3ikl,m3ijl,m3ijk2
    real(dp) :: m1i,m1j,m1k,m1l
    integer :: d,i,j,k,l,info
    call recintab(kappa,lower,upper,mean,varcov,m,out%status)
    out%shape=kappa+1; out%probability=m(1)
    if(m(1)<=0.0_dp) then; out%raw=m; return; end if
    out%raw=m/m(1); d=size(kappa); allocate(alpha(d),out%mean(d)); alpha=0
    do i=1,d; alpha=0; alpha(i)=1; out%mean(i)=raw_moment_at(out%raw,out%shape,alpha); end do
    if(all(kappa>=2)) then
      allocate(out%covariance(d,d)); out%covariance=0.0_dp
      do i=1,d; do j=1,d
        alpha=0; alpha(i)=alpha(i)+1; alpha(j)=alpha(j)+1
        out%covariance(i,j)=raw_moment_at(out%raw,out%shape,alpha)-out%mean(i)*out%mean(j)
      end do; end do
    end if
    if(all(kappa>=3)) then
      allocate(out%cum3(d,d,d),out%marginal_skewness(d)); out%cum3=0.0_dp
      do i=1,d; do j=1,d; do k=1,d
        alpha=0; alpha(i)=alpha(i)+1; alpha(j)=alpha(j)+1; alpha(k)=alpha(k)+1
        m3ijk=raw_moment_at(out%raw,out%shape,alpha)
        alpha=0; alpha(j)=alpha(j)+1; alpha(k)=alpha(k)+1; m2ij=raw_moment_at(out%raw,out%shape,alpha)
        alpha=0; alpha(i)=alpha(i)+1; alpha(k)=alpha(k)+1; m2kl=raw_moment_at(out%raw,out%shape,alpha)
        alpha=0; alpha(i)=alpha(i)+1; alpha(j)=alpha(j)+1; m3ijk2=raw_moment_at(out%raw,out%shape,alpha)
        out%cum3(i,j,k)=m3ijk-out%mean(i)*m2ij-out%mean(j)*m2kl-out%mean(k)*m3ijk2 &
                         +2.0_dp*out%mean(i)*out%mean(j)*out%mean(k)
      end do; end do; end do
      do i=1,d
        out%marginal_skewness(i)=out%cum3(i,i,i)/out%covariance(i,i)**1.5_dp
      end do
      allocate(inv(d,d)); call pd_solve(out%covariance,inv,info=info)
      if(info==0) then
        out%gamma1_mardia = contract_cum3(out%cum3,inv)
      end if
      deallocate(inv)
    end if
    if(all(kappa>=4)) then
      allocate(out%cum4(d,d,d,d),out%marginal_excess_kurtosis(d)); out%cum4=0.0_dp
      do i=1,d; do j=1,d; do k=1,d; do l=1,d
        alpha=0; alpha(i)=alpha(i)+1; alpha(j)=alpha(j)+1; alpha(k)=alpha(k)+1; alpha(l)=alpha(l)+1
        m3ijk2=raw_moment_at(out%raw,out%shape,alpha)
        m1i=out%mean(i); m1j=out%mean(j); m1k=out%mean(k); m1l=out%mean(l)
        m3jkl=raw3(out%raw,out%shape,j,k,l); m3ikl=raw3(out%raw,out%shape,i,k,l)
        m3ijl=raw3(out%raw,out%shape,i,j,l); m3ijk=raw3(out%raw,out%shape,i,j,k)
        m2ij=raw2(out%raw,out%shape,i,j); m2kl=raw2(out%raw,out%shape,k,l)
        out%cum4(i,j,k,l)=m3ijk2-(m1i*m3jkl+m1j*m3ikl+m1k*m3ijl+m1l*m3ijk) &
          -(m2ij*m2kl+raw2(out%raw,out%shape,i,k)*raw2(out%raw,out%shape,j,l) &
          +raw2(out%raw,out%shape,i,l)*raw2(out%raw,out%shape,j,k)) &
          +2.0_dp*(m1i*m1j*m2kl+m1i*m1k*raw2(out%raw,out%shape,j,l) &
          +m1i*m1l*raw2(out%raw,out%shape,j,k)+m1j*m1k*raw2(out%raw,out%shape,i,l) &
          +m1j*m1l*raw2(out%raw,out%shape,i,k)+m1k*m1l*m2ij)-6.0_dp*m1i*m1j*m1k*m1l
      end do; end do; end do; end do
      do i=1,d; out%marginal_excess_kurtosis(i)=out%cum4(i,i,i,i)/out%covariance(i,i)**2; end do
      allocate(inv(d,d)); call pd_solve(out%covariance,inv,info=info)
      if(info==0) then
        out%gamma2_mardia=0.0_dp
        do i=1,d; do j=1,d; do k=1,d; do l=1,d
          out%gamma2_mardia=out%gamma2_mardia+out%cum4(i,j,k,l)*inv(i,j)*inv(k,l)
        end do; end do; end do; end do
        out%beta2_mardia=out%gamma2_mardia+real(d*(d+2),dp)
      end if
      deallocate(inv)
    end if
  contains
    pure real(dp) function raw2(raw,shape,i,j) result(v)
      real(dp),intent(in)::raw(:); integer,intent(in)::shape(:),i,j; integer::a(size(shape))
      a=0;a(i)=a(i)+1;a(j)=a(j)+1;v=raw_moment_at(raw,shape,a)
    end function raw2
    pure real(dp) function raw3(raw,shape,i,j,k) result(v)
      real(dp),intent(in)::raw(:); integer,intent(in)::shape(:),i,j,k; integer::a(size(shape))
      a=0;a(i)=a(i)+1;a(j)=a(j)+1;a(k)=a(k)+1;v=raw_moment_at(raw,shape,a)
    end function raw3
    pure real(dp) function contract_cum3(c,ainv) result(v)
      real(dp),intent(in)::c(:,:,:),ainv(:,:)
      integer::a,b,c1,d1,e,f,n
      v=0.0_dp;n=size(ainv,1)
      do a=1,n;do b=1,n;do c1=1,n;do d1=1,n;do e=1,n;do f=1,n
        v=v+c(a,b,c1)*c(d1,e,f)*ainv(a,d1)*ainv(b,e)*ainv(c1,f)
      end do;end do;end do;end do;end do;end do
    end function contract_cum3
  end function mom_mtruncnorm

  function mom2cum(raw,shape) result(out)
    real(dp), intent(in) :: raw(:)
    integer, intent(in) :: shape(:)
    type(trunc_moment_result) :: out
    integer, allocatable :: alpha(:)
    real(dp), allocatable :: inv(:,:)
    real(dp) :: m2ij,m2kl,m3ijk,m3jkl,m3ikl,m3ijl,m4ijkl
    real(dp) :: m1i,m1j,m1k,m1l
    integer :: d,i,j,k,l,info
    allocate(out%shape(size(shape)),out%raw(size(raw)))
    out%shape=shape; out%raw=raw; out%probability=1.0_dp
    if(size(raw)/=product(shape) .or. abs(raw(1)-1.0_dp)>1.0e-10_dp) then
      out%status=-1; return
    end if
    d=size(shape); allocate(alpha(d),out%mean(d))
    if(any(shape<2)) then; out%status=-2; return; end if
    do i=1,d
      alpha=0; alpha(i)=1; out%mean(i)=raw_moment_at(raw,shape,alpha)
    end do
    if(all(shape>=3)) then
      allocate(out%covariance(d,d)); out%covariance=0.0_dp
      do i=1,d; do j=1,d
        out%covariance(i,j)=raw2(raw,shape,i,j)-out%mean(i)*out%mean(j)
      end do; end do
    end if
    if(all(shape>=4)) then
      allocate(out%cum3(d,d,d),out%marginal_skewness(d)); out%cum3=0.0_dp
      do i=1,d; do j=1,d; do k=1,d
        m3ijk=raw3(raw,shape,i,j,k)
        out%cum3(i,j,k)=m3ijk-out%mean(i)*raw2(raw,shape,j,k) &
          -out%mean(j)*raw2(raw,shape,i,k)-out%mean(k)*raw2(raw,shape,i,j) &
          +2.0_dp*out%mean(i)*out%mean(j)*out%mean(k)
      end do; end do; end do
      do i=1,d
        out%marginal_skewness(i)=out%cum3(i,i,i)/out%covariance(i,i)**1.5_dp
      end do
      allocate(inv(d,d)); call pd_solve(out%covariance,inv,info=info)
      if(info==0) out%gamma1_mardia=contract_cum3_local(out%cum3,inv)
      deallocate(inv)
    end if
    if(all(shape>=5)) then
      allocate(out%cum4(d,d,d,d),out%marginal_excess_kurtosis(d)); out%cum4=0.0_dp
      do i=1,d; do j=1,d; do k=1,d; do l=1,d
        m4ijkl=raw4(raw,shape,i,j,k,l)
        m1i=out%mean(i); m1j=out%mean(j); m1k=out%mean(k); m1l=out%mean(l)
        m3jkl=raw3(raw,shape,j,k,l); m3ikl=raw3(raw,shape,i,k,l)
        m3ijl=raw3(raw,shape,i,j,l); m3ijk=raw3(raw,shape,i,j,k)
        m2ij=raw2(raw,shape,i,j); m2kl=raw2(raw,shape,k,l)
        out%cum4(i,j,k,l)=m4ijkl-(m1i*m3jkl+m1j*m3ikl+m1k*m3ijl+m1l*m3ijk) &
          -(m2ij*m2kl+raw2(raw,shape,i,k)*raw2(raw,shape,j,l) &
          +raw2(raw,shape,i,l)*raw2(raw,shape,j,k)) &
          +2.0_dp*(m1i*m1j*m2kl+m1i*m1k*raw2(raw,shape,j,l) &
          +m1i*m1l*raw2(raw,shape,j,k)+m1j*m1k*raw2(raw,shape,i,l) &
          +m1j*m1l*raw2(raw,shape,i,k)+m1k*m1l*m2ij)-6.0_dp*m1i*m1j*m1k*m1l
      end do; end do; end do; end do
      do i=1,d
        out%marginal_excess_kurtosis(i)=out%cum4(i,i,i,i)/out%covariance(i,i)**2
      end do
      allocate(inv(d,d)); call pd_solve(out%covariance,inv,info=info)
      if(info==0) then
        do i=1,d; do j=1,d; do k=1,d; do l=1,d
          out%gamma2_mardia=out%gamma2_mardia+out%cum4(i,j,k,l)*inv(i,j)*inv(k,l)
        end do; end do; end do; end do
        out%beta2_mardia=out%gamma2_mardia+real(d*(d+2),dp)
      end if
      deallocate(inv)
    end if
  contains
    pure real(dp) function raw2(r,sh,i,j) result(v)
      real(dp),intent(in)::r(:); integer,intent(in)::sh(:),i,j; integer::a(size(sh))
      a=0;a(i)=a(i)+1;a(j)=a(j)+1;v=raw_moment_at(r,sh,a)
    end function raw2
    pure real(dp) function raw3(r,sh,i,j,k) result(v)
      real(dp),intent(in)::r(:); integer,intent(in)::sh(:),i,j,k; integer::a(size(sh))
      a=0;a(i)=a(i)+1;a(j)=a(j)+1;a(k)=a(k)+1;v=raw_moment_at(r,sh,a)
    end function raw3
    pure real(dp) function raw4(r,sh,i,j,k,l) result(v)
      real(dp),intent(in)::r(:); integer,intent(in)::sh(:),i,j,k,l; integer::a(size(sh))
      a=0;a(i)=a(i)+1;a(j)=a(j)+1;a(k)=a(k)+1;a(l)=a(l)+1;v=raw_moment_at(r,sh,a)
    end function raw4
    pure real(dp) function contract_cum3_local(c,ainv) result(v)
      real(dp),intent(in)::c(:,:,:),ainv(:,:)
      integer::a,b,c1,d1,e,f,n
      v=0.0_dp;n=size(ainv,1)
      do a=1,n;do b=1,n;do c1=1,n;do d1=1,n;do e=1,n;do f=1,n
        v=v+c(a,b,c1)*c(d1,e,f)*ainv(a,d1)*ainv(b,e)*ainv(c1,f)
      end do;end do;end do;end do;end do;end do
    end function contract_cum3_local
  end function mom2cum

  function sample_mardia_measures(data,correct) result(out)
    real(dp), intent(in) :: data(:,:)
    logical, intent(in), optional :: correct
    type(mardia_result) :: out
    real(dp), allocatable :: mean(:),y(:,:),s(:,:),inv(:,:),qmat(:,:)
    real(dp) :: f,kfac,r,std_b2,chisq_x,dfchi
    integer :: n,d,i,info
    logical :: corr
    n=size(data,1); d=size(data,2); out%n=n; corr=.false.; if(present(correct)) corr=correct
    if(n<4) then; out%status=-1; return; end if
    allocate(mean(d),y(n,d),s(d,d),inv(d,d),qmat(n,n)); mean=sum(data,dim=1)/real(n,dp)
    do i=1,n; y(i,:)=data(i,:)-mean; end do
    s=matmul(transpose(y),y)/real(n-1,dp)
    f=merge(1.0_dp,real(n-1,dp)/real(n,dp),corr); s=s*f
    call pd_solve(s,inv,info=info); if(info/=0) then; out%status=info; return; end if
    qmat=matmul(matmul(y,inv),transpose(y)); out%b1=sum(qmat**3)/real(n*n,dp)
    out%b2=0.0_dp
    do i=1,n; out%b2=out%b2+dot_product(y(i,:),matmul(inv,y(i,:)))**2; end do
    out%b2=out%b2/real(n,dp)
    out%g1=real(n*(n-1),dp)*out%b1/real((n-2)**2,dp)
    out%g2=real(n-1,dp)*(real(n+1,dp)*out%b2-real((n-1)*d*(d+2),dp))/real((n-2)*(n-3),dp)
    kfac=real((d+1)*(n+1)*(n+3),dp)/real(n*((n+1)*(d+1)-6),dp)
    dfchi=real(d*(d+1)*(d+2),dp)/6.0_dp; chisq_x=out%b1*real(n,dp)*kfac/6.0_dp
    out%p_b1=1.0_dp-chi_square_cdf(chisq_x,dfchi)
    if(n-d-1>0) then
      r=sqrt(real((n+3)*(n+5),dp)/real(n-d-1,dp))
      std_b2=(real(n+1,dp)*out%b2-real(d*(d+2)*(n-1),dp))*r / &
             real(8*d*(d+2)*(n-3)*(n-d-1),dp)
      out%p_b2=erfc(abs(std_b2)/sqrt(2.0_dp))
    else
      out%p_b2=ieee_value(1.0_dp,ieee_quiet_nan)
    end if
  contains
    pure real(dp) function chi_square_cdf(x,df) result(p)
      real(dp),intent(in)::x,df
      p=regularized_gamma_p(0.5_dp*df,0.5_dp*x)
    end function chi_square_cdf
    pure real(dp) function regularized_gamma_p(a,x) result(p)
      real(dp),intent(in)::a,x
      integer::i
      real(dp)::ap,del,sumv,b,c,d,h,an
      if(x<=0.0_dp) then;p=0.0_dp;return;end if
      if(x<a+1.0_dp) then
        ap=a;del=1.0_dp/a;sumv=del
        do i=1,1000;ap=ap+1.0_dp;del=del*x/ap;sumv=sumv+del;if(abs(del)<abs(sumv)*1e-15_dp)exit;end do
        p=sumv*exp(-x+a*log(x)-log_gamma(a))
      else
        b=x+1.0_dp-a;c=1.0e300_dp;d=1.0_dp/b;h=d
        do i=1,1000
          an=-real(i,dp)*(real(i,dp)-a);b=b+2.0_dp;d=an*d+b;if(abs(d)<1e-300_dp)d=1e-300_dp
          c=b+an/c;if(abs(c)<1e-300_dp)c=1e-300_dp;d=1.0_dp/d;del=d*c;h=h*del;if(abs(del-1.0_dp)<1e-15_dp)exit
        end do
        p=1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
      end if
      p=max(0.0_dp,min(1.0_dp,p))
    end function regularized_gamma_p
  end function sample_mardia_measures
end module mnormt_moments
