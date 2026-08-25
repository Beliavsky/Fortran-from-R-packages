module greybox_stats
    use greybox_kinds, only: dp
    use greybox_special, only: gamma_p, beta_inc
    use greybox_linalg, only: least_squares, invert_matrix, covariance_matrix, pearson_cor
    implicit none
    private
    public :: cramer_v, mcor, pcor_matrix, determination, association_numeric
    public :: correlation_result

    type :: correlation_result
        real(dp) :: value=0.0_dp
        real(dp) :: statistic=0.0_dp
        real(dp) :: p_value=1.0_dp
        integer :: df=0
        integer :: df_residual=0
    end type correlation_result

contains

    function cramer_v(x,y) result(res)
        integer,intent(in)::x(:),y(:)
        type(correlation_result)::res
        integer::xmin,xmax,ymin,ymax,nx,ny,i,j,n
        real(dp),allocatable::tab(:,:),rs(:),cs(:)
        real(dp)::expected,chi2
        n=size(x);xmin=minval(x);xmax=maxval(x);ymin=minval(y);ymax=maxval(y)
        nx=xmax-xmin+1;ny=ymax-ymin+1;allocate(tab(nx,ny),rs(nx),cs(ny));tab=0.0_dp
        do i=1,n;tab(x(i)-xmin+1,y(i)-ymin+1)=tab(x(i)-xmin+1,y(i)-ymin+1)+1.0_dp;end do
        rs=sum(tab,dim=2);cs=sum(tab,dim=1);chi2=0.0_dp
        do i=1,nx;do j=1,ny
            expected=rs(i)*cs(j)/real(n,dp)
            if(expected>0.0_dp)chi2=chi2+(tab(i,j)-expected)**2/expected
        end do;end do
        res%statistic=chi2;res%df=(nx-1)*(ny-1);res%df_residual=n-res%df
        if(min(nx-1,ny-1)>0)res%value=sqrt(chi2/(real(n,dp)*real(min(nx-1,ny-1),dp)))
        res%p_value=1.0_dp-gamma_p(0.5_dp*real(res%df,dp),0.5_dp*chi2)
    end function cramer_v

    function mcor(x,y) result(res)
        real(dp),intent(in)::x(:,:),y(:)
        type(correlation_result)::res
        real(dp),allocatable::xx(:,:),beta(:),fit(:)
        real(dp)::sst,sse,r2,f
        integer::n,p,info
        n=size(x,1);p=size(x,2);allocate(xx(n,p+1),beta(p+1),fit(n));xx(:,1)=1.0_dp;xx(:,2:)=x
        call least_squares(xx,y,beta,info)
        fit=matmul(xx,beta);sst=sum((y-sum(y)/real(n,dp))**2);sse=sum((y-fit)**2)
        if(sst>0.0_dp)then;r2=max(0.0_dp,min(1.0_dp,1.0_dp-sse/sst));else;r2=0.0_dp;end if
        res%value=sqrt(r2);res%df=p;res%df_residual=n-p-1
        if(res%df>0.and.res%df_residual>0.and.r2<1.0_dp)then
            f=(r2/real(res%df,dp))/((1.0_dp-r2)/real(res%df_residual,dp));res%statistic=f
            res%p_value=1.0_dp-beta_inc(0.5_dp*real(res%df,dp),0.5_dp*real(res%df_residual,dp), &
                real(res%df,dp)*f/(real(res%df,dp)*f+real(res%df_residual,dp)))
        end if
    end function mcor

    function pcor_matrix(x) result(r)
        real(dp),intent(in)::x(:,:)
        real(dp)::r(size(x,2),size(x,2))
        real(dp),allocatable::c(:,:),d(:,:),p(:,:)
        integer::nvar,i,j,info
        nvar=size(x,2);allocate(c(nvar,nvar),d(nvar,nvar),p(nvar,nvar));c=covariance_matrix(x)
        do i=1,nvar;do j=1,nvar
            if(c(i,i)>0.0_dp.and.c(j,j)>0.0_dp)then;d(i,j)=c(i,j)/sqrt(c(i,i)*c(j,j));else;d(i,j)=0.0_dp;end if
        end do;end do
        call invert_matrix(d,p,info);r=0.0_dp
        if(info/=0)return
        do i=1,nvar;r(i,i)=1.0_dp;do j=i+1,nvar
            if(p(i,i)>0.0_dp.and.p(j,j)>0.0_dp)then;r(i,j)=-p(i,j)/sqrt(p(i,i)*p(j,j));r(j,i)=r(i,j);end if
        end do;end do
    end function pcor_matrix

    function determination(x) result(r2)
        real(dp),intent(in)::x(:,:)
        real(dp)::r2(size(x,2))
        real(dp),allocatable::others(:,:)
        type(correlation_result)::res
        integer::p,j,k,c
        p=size(x,2)
        if(p==1)then;r2=0.0_dp;return;end if
        allocate(others(size(x,1),p-1))
        do j=1,p;c=0;do k=1,p;if(k==j)cycle;c=c+1;others(:,c)=x(:,k);end do;res=mcor(others,x(:,j));r2(j)=res%value**2;end do
    end function determination

    pure function association_numeric(x) result(r)
        real(dp),intent(in)::x(:,:)
        real(dp)::r(size(x,2),size(x,2))
        integer::i,j,p
        p=size(x,2);r=0.0_dp
        do i=1,p;r(i,i)=1.0_dp;do j=i+1,p;r(i,j)=pearson_cor(x(:,i),x(:,j));r(j,i)=r(i,j);end do;end do
    end function association_numeric

end module greybox_stats
