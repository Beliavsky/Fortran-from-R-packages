! SPDX-License-Identifier: GPL-2.0-only
module ks_tests
  use ks_kinds, only: dp, pi
  use ks_linalg, only: determinant_spd, covariance_matrix
  use ks_normal, only: mvn_pdf, normal_cdf
  use ks_kde, only: kde_model, fit_kde, kdde_eval
  implicit none
  private
  public :: qr0, kde_two_sample_test, kde_test_result
  type :: kde_test_result
    real(dp) :: tstat=0.0_dp,zstat=0.0_dp,pvalue=1.0_dp,mean=0.0_dp,variance=0.0_dp
    real(dp) :: psi1=0.0_dp,psi12=0.0_dp,psi21=0.0_dp,psi2=0.0_dp
  end type
contains
  function qr0(x,y,sigma) result(q)
    real(dp),intent(in)::x(:,:),y(:,:),sigma(:,:)
    real(dp)::q,zero(size(x,2)),delta(size(x,2))
    integer::i,j,nx,ny,d
    nx=size(x,1);ny=size(y,1);d=size(x,2)
    if(size(y,2)/=d.or.size(sigma,1)/=d.or.size(sigma,2)/=d) error stop 'qr0: shape'
    zero=0.0_dp;q=0.0_dp
    do i=1,nx
      do j=1,ny
        delta=x(i,:)-y(j,:)
        q=q+mvn_pdf(delta,zero,sigma)
      end do
    end do
    q=q/real(nx*ny,dp)
  end function

  function kde_two_sample_test(x1,x2,H1,H2) result(res)
    real(dp),intent(in)::x1(:,:),x2(:,:),H1(:,:),H2(:,:)
    type(kde_test_result)::res
    integer::n1,n2,d,j
    real(dp)::k0,det1,det2,var1,var2
    real(dp),allocatable::s1(:,:),s2(:,:),mean1(:,:),mean2(:,:),g1(:,:),g2(:,:)
    type(kde_model)::m1,m2
    n1=size(x1,1);n2=size(x2,1);d=size(x1,2)
    if(size(x2,2)/=d.or.n1<2.or.n2<2) error stop 'kde_two_sample_test: shape'
    res%psi1=qr0(x1,x1,H1);res%psi2=qr0(x2,x2,H2)
    res%psi12=qr0(x1,x2,H1);res%psi21=qr0(x2,x1,H2)
    res%tstat=res%psi1+res%psi2-res%psi12-res%psi21
    det1=determinant_spd(H1);det2=determinant_spd(H2)
    k0=(2.0_dp*pi)**(-0.5_dp*real(d,dp))
    res%mean=(1.0_dp/real(n1,dp)/sqrt(det1)+1.0_dp/real(n2,dp)/sqrt(det2))*k0
    allocate(s1(d,d),s2(d,d),mean1(1,d),mean2(1,d))
    call covariance_matrix(x1,s1);call covariance_matrix(x2,s2)
    do j=1,d
      mean1(1,j)=sum(x1(:,j))/real(n1,dp);mean2(1,j)=sum(x2(:,j))/real(n2,dp)
    end do
    call fit_kde(x1,m1,H=H1);call fit_kde(x2,m2,H=H2)
    call kdde_eval(m1,mean1,1,g1);call kdde_eval(m2,mean2,1,g2)
    var1=dot_product(g1(1,:),matmul(s1,g1(1,:)))
    var2=dot_product(g2(1,:),matmul(s2,g2(1,:)))
    res%variance=3.0_dp*(real(n1,dp)*var1+real(n2,dp)*var2)/real(n1+n2,dp)*(1.0_dp/real(n1,dp)+1.0_dp/real(n2,dp))
    if(res%variance>0.0_dp) then
      res%zstat=(res%tstat-res%mean)/sqrt(res%variance)
      res%pvalue=1.0_dp-normal_cdf(res%zstat,0.0_dp,1.0_dp)
      if(res%pvalue<=0.0_dp) res%pvalue=normal_cdf(-abs(res%zstat),0.0_dp,1.0_dp)
    end if
  end function
end module ks_tests
