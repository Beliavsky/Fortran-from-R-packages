! SPDX-License-Identifier: GPL-3.0-only
module pa_factor_models
  use pa_kinds, only: dp
  use pa_types, only: factor_model_result
  use pa_linalg, only: jacobi_eigen_sym
  use pa_statistics, only: sample_moments, sample_coskewness, sample_cokurtosis
  implicit none
  private
  public :: fit_statistical_factor_model, factor_model_covariance
  public :: covariance_sf, coskewness_sf, cokurtosis_sf
  public :: covariance_mf, coskewness_mf, cokurtosis_mf

contains

  subroutine fit_statistical_factor_model(returns,nfactors,model,info)
    real(dp),intent(in)::returns(:,:)
    integer,intent(in)::nfactors
    type(factor_model_result),intent(out)::model
    integer,intent(out),optional::info
    real(dp),allocatable::mu(:),sigma(:,:),eval(:),evec(:,:),xc(:,:)
    integer::nobs,nassets,i,istat
    nobs=size(returns,1)
    nassets=size(returns,2)
    if (present(info)) info=0
    if (nfactors<1 .or. nfactors>nassets .or. nobs<2) then
      if (present(info)) info=-1
      return
    end if
    allocate(mu(nassets),sigma(nassets,nassets),eval(nassets),evec(nassets,nassets),xc(nobs,nassets))
    call sample_moments(returns,mu,sigma)
    call jacobi_eigen_sym(sigma,eval,evec,istat)
    if (istat/=0 .and. present(info)) info=istat
    do i=1,nobs
      xc(i,:)=returns(i,:)-mu
    end do
    allocate(model%means(nassets),model%loadings(nassets,nfactors), &
             model%factors(nobs,nfactors),model%residuals(nobs,nassets), &
             model%eigenvalues(nfactors))
    model%means=mu
    model%loadings=evec(:,1:nfactors)
    model%eigenvalues=eval(1:nfactors)
    model%factors=matmul(xc,model%loadings)
    model%residuals=xc-matmul(model%factors,transpose(model%loadings))
    model%nobs=nobs
    model%nassets=nassets
    model%nfactors=nfactors
  end subroutine fit_statistical_factor_model

  subroutine factor_model_covariance(model,sigma)
    type(factor_model_result),intent(in)::model
    real(dp),intent(out)::sigma(:,:)
    real(dp),allocatable::fm(:),fs(:,:),rm(:),rs(:,:),rv(:)
    integer::n,k,i
    n=model%nassets
    k=model%nfactors
    allocate(fm(k),fs(k,k),rm(n),rs(n,n),rv(n))
    call sample_moments(model%factors,fm,fs)
    call sample_moments(model%residuals,rm,rs)
    do i=1,n
      rv(i)=rs(i,i)
    end do
    call covariance_mf(model%loadings,rv,fs,sigma)
  end subroutine factor_model_covariance

  subroutine covariance_sf(beta,stock_m2,factor_m2,sigma)
    real(dp),intent(in)::beta(:),stock_m2(:),factor_m2
    real(dp),intent(out)::sigma(:,:)
    integer::i,n
    n=size(beta)
    do i=1,n
      sigma(i,:)=beta(i)*beta*factor_m2
      sigma(i,i)=sigma(i,i)+stock_m2(i)
    end do
  end subroutine covariance_sf

  subroutine coskewness_sf(beta,stock_m3,factor_m3,m3)
    real(dp),intent(in)::beta(:),stock_m3(:),factor_m3
    real(dp),intent(out)::m3(:,:)
    integer::n,i,j,k,col
    n=size(beta)
    do i=1,n
      do j=1,n
        do k=1,n
          col=(j-1)*n+k
          m3(i,col)=beta(i)*beta(j)*beta(k)*factor_m3
          if (i==j .and. j==k) m3(i,col)=m3(i,col)+stock_m3(i)
        end do
      end do
    end do
  end subroutine coskewness_sf

  subroutine cokurtosis_sf(beta,stock_m2,stock_m4,factor_m2,factor_m4,m4)
    real(dp),intent(in)::beta(:),stock_m2(:),stock_m4(:),factor_m2,factor_m4
    real(dp),intent(out)::m4(:,:)
    real(dp),allocatable::b(:,:),fm2(:,:),fm4(:,:)
    integer::n
    n=size(beta)
    allocate(b(n,1),fm2(1,1),fm4(1,1))
    b(:,1)=beta
    fm2(1,1)=factor_m2
    fm4(1,1)=factor_m4
    call cokurtosis_mf(b,stock_m2,stock_m4,fm2,fm4,m4)
  end subroutine cokurtosis_sf

  subroutine covariance_mf(beta,stock_m2,factor_m2,sigma)
    real(dp),intent(in)::beta(:,:),stock_m2(:),factor_m2(:,:)
    real(dp),intent(out)::sigma(:,:)
    integer::i
    sigma=matmul(matmul(beta,factor_m2),transpose(beta))
    do i=1,size(beta,1)
      sigma(i,i)=sigma(i,i)+stock_m2(i)
    end do
    sigma=0.5_dp*(sigma+transpose(sigma))
  end subroutine covariance_mf

  subroutine coskewness_mf(beta,stock_m3,factor_m3,m3)
    real(dp),intent(in)::beta(:,:),stock_m3(:),factor_m3(:,:)
    real(dp),intent(out)::m3(:,:)
    integer::n,kf,i,j,k,a,b,c,col,fcol
    real(dp)::value
    n=size(beta,1)
    kf=size(beta,2)
    m3=0.0_dp
    do i=1,n
      do j=1,n
        do k=1,n
          value=0.0_dp
          do a=1,kf
            do b=1,kf
              do c=1,kf
                fcol=(b-1)*kf+c
                value=value+beta(i,a)*beta(j,b)*beta(k,c)*factor_m3(a,fcol)
              end do
            end do
          end do
          if (i==j .and. j==k) value=value+stock_m3(i)
          col=(j-1)*n+k
          m3(i,col)=value
        end do
      end do
    end do
  end subroutine coskewness_mf

  subroutine cokurtosis_mf(beta,stock_m2,stock_m4,factor_m2,factor_m4,m4)
    real(dp),intent(in)::beta(:,:),stock_m2(:),stock_m4(:),factor_m2(:,:),factor_m4(:,:)
    real(dp),intent(out)::m4(:,:)
    real(dp),allocatable::syscov(:,:)
    real(dp)::value
    integer::n,kf,i,j,k,l,a,b,c,d,col,fcol
    n=size(beta,1)
    kf=size(beta,2)
    allocate(syscov(n,n))
    syscov=matmul(matmul(beta,factor_m2),transpose(beta))
    m4=0.0_dp
    do i=1,n
      do j=1,n
        do k=1,n
          do l=1,n
            value=0.0_dp
            do a=1,kf
              do b=1,kf
                do c=1,kf
                  do d=1,kf
                    fcol=((b-1)*kf+(c-1))*kf+d
                    value=value+beta(i,a)*beta(j,b)*beta(k,c)*beta(l,d)*factor_m4(a,fcol)
                  end do
                end do
              end do
            end do
            if (i==j) value=value+stock_m2(i)*syscov(k,l)
            if (i==k) value=value+stock_m2(i)*syscov(j,l)
            if (i==l) value=value+stock_m2(i)*syscov(j,k)
            if (j==k) value=value+stock_m2(j)*syscov(i,l)
            if (j==l) value=value+stock_m2(j)*syscov(i,k)
            if (k==l) value=value+stock_m2(k)*syscov(i,j)
            if (i==j .and. j==k .and. k==l) then
              value=value+stock_m4(i)
            else
              if (i==j .and. k==l) value=value+stock_m2(i)*stock_m2(k)
              if (i==k .and. j==l) value=value+stock_m2(i)*stock_m2(j)
              if (i==l .and. j==k) value=value+stock_m2(i)*stock_m2(j)
            end if
            col=((j-1)*n+(k-1))*n+l
            m4(i,col)=value
          end do
        end do
      end do
    end do
  end subroutine cokurtosis_mf

end module pa_factor_models
