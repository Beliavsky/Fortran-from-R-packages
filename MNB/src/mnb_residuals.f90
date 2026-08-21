! SPDX-License-Identifier: GPL-2.0-or-later
module mnb_residuals
  use mnb_kinds, only : dp
  use mnb_types, only : mnb_fit_result,mnb_residual_result
  use mnb_core, only : fit_mnb,mnb_cluster_sums
  use mnb_math, only : invert_matrix,normal_quantile
  implicit none
  private
  public :: residuals_mnb, randomized_quantile_residuals, nb_total_pmf
contains
  pure real(dp) function nb_total_pmf(y,phi,mu) result(p)
    integer,intent(in)::y
    real(dp),intent(in)::phi,mu
    real(dp)::q
    if(y<0 .or. phi<=0.0_dp .or. mu<0.0_dp)then;p=0.0_dp;return;end if
    q=phi/(phi+mu)
    p=exp(log_gamma(real(y,dp)+phi)-log_gamma(real(y+1,dp))-log_gamma(phi)+phi*log(q)+real(y,dp)*log(max(1.0_dp-q,tiny(1.0_dp))))
  end function nb_total_pmf

  subroutine randomized_quantile_residuals(par,y,x,n,mi,rq,offset)
    real(dp),intent(in)::par(:),y(:),x(:,:)
    integer,intent(in)::n,mi
    real(dp),intent(out)::rq(n)
    real(dp),intent(in),optional::offset(:)
    real(dp),allocatable::eta(:),mu(:),ys(:),ms(:)
    real(dp)::flo,fhi,u
    integer::i,k,yi
    allocate(eta(n*mi),mu(n*mi),ys(n),ms(n));eta=matmul(x,par(2:));if(present(offset))eta=eta+offset
    mu=exp(eta);call mnb_cluster_sums(y,mu,n,mi,ys,ms)
    do i=1,n
      yi=nint(ys(i));flo=0.0_dp
      do k=0,yi-1;flo=flo+nb_total_pmf(k,par(1),ms(i));end do
      fhi=flo+nb_total_pmf(yi,par(1),ms(i));call random_number(u);u=flo+u*(fhi-flo)
      u=min(1.0_dp-1.0e-15_dp,max(1.0e-15_dp,u));rq(i)=normal_quantile(u)
    end do
  end subroutine randomized_quantile_residuals

  function residuals_mnb(start,y,x,n,mi,offset) result(r)
    real(dp),intent(in)::start(:),y(:),x(:,:)
    integer,intent(in)::n,mi
    real(dp),intent(in),optional::offset(:)
    type(mnb_residual_result)::r
    type(mnb_fit_result)::fit
    real(dp),allocatable::eta(:),mu(:),ys(:),ms(:),w(:,:),xtwx(:,:),inv(:,:),a(:,:),block(:,:)
    real(dp),allocatable::ai(:),hmean(:)
    real(dp)::phi,aux1,aux2,term
    integer::nn,p,i,j,l,u,k
    logical::ok
    fit=fit_mnb(start,y,x,n,mi,offset);nn=n*mi;p=size(x,2);phi=fit%par(1)
    allocate(eta(nn),mu(nn),ys(n),ms(n),w(nn,nn),xtwx(p,p),inv(p,p),ai(nn),hmean(n))
    allocate(r%weighted(nn),r%standardized_weighted(nn),r%pearson(nn),r%standardized_pearson(nn),r%deviance(n),r%leverage(nn))
    eta=matmul(x,fit%par(2:));if(present(offset))eta=eta+offset;mu=exp(eta);call mnb_cluster_sums(y,mu,n,mi,ys,ms)
    w=0.0_dp
    do i=1,nn;w(i,i)=mu(i);end do
    ! Source-compatible R formula: global outer product, with cluster denominator by row.
    do i=1,nn
      k=(i-1)/mi+1
      do j=1,nn;w(i,j)=w(i,j)+mu(i)*mu(j)/(phi+ms(k));end do
    end do
    xtwx=matmul(transpose(x),matmul(w,x));call invert_matrix(xtwx,inv,ok)
    if(.not.ok)error stop 'residuals_mnb: singular X''WX'
    do i=1,n
      l=(i-1)*mi+1;u=i*mi;allocate(block(mi,mi),a(mi,p));block=sqrt(max(w(l:u,l:u),0.0_dp))
      a=matmul(block,x(l:u,:))
      do j=1,mi;r%leverage(l+j-1)=dot_product(a(j,:),matmul(inv,a(j,:)));end do
      deallocate(block,a)
    end do
    do i=1,n
      l=(i-1)*mi+1;u=i*mi;ai(l:u)=(phi+ys(i))/(phi+ms(i))
    end do
    r%weighted=(y-ai*mu)/sqrt(max([(w(i,i),i=1,nn)],tiny(1.0_dp)))
    r%standardized_weighted=r%weighted/sqrt(max(1.0_dp-r%leverage,tiny(1.0_dp)))
    r%pearson=sqrt(phi)*(y-mu)/sqrt(max(mu*(phi+mu),tiny(1.0_dp)))
    r%standardized_pearson=r%pearson/sqrt(max(1.0_dp-r%leverage,tiny(1.0_dp)))
    do i=1,n
      l=(i-1)*mi+1;u=i*mi;aux1=0.0_dp
      do j=l,u
        if(y(j)>0.0_dp)then
          term=y(j)*(phi+ms(i))/(mu(j)*(phi+ys(i)));aux1=aux1+y(j)*log(term)
        end if
      end do
      aux2=phi*log((phi+ms(i))/(phi+ys(i)));term=max(0.0_dp,2.0_dp*(aux2+aux1))
      hmean(i)=sum(r%leverage(l:u))/real(mi,dp)
      r%deviance(i)=sign(1.0_dp,ys(i)-ms(i))*sqrt(term)/sqrt(max(1.0_dp-hmean(i),tiny(1.0_dp)))
    end do
  end function residuals_mnb
end module mnb_residuals
