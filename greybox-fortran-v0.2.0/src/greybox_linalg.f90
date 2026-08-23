module greybox_linalg
    use greybox_kinds, only: dp
    implicit none
    private
    public :: solve_linear, least_squares, invert_matrix, covariance_matrix, pearson_cor

contains

    subroutine solve_linear(a,b,x,info)
        real(dp),intent(in)::a(:,:),b(:)
        real(dp),intent(out)::x(:)
        integer,intent(out),optional::info
        real(dp),allocatable::aa(:,:),bb(:),row(:)
        real(dp)::fac,tmp
        integer::n,i,k,piv,istat
        n=size(b);istat=0
        allocate(aa(n,n),bb(n),row(n));aa=a;bb=b
        do k=1,n-1
            piv=k
            do i=k+1,n;if(abs(aa(i,k))>abs(aa(piv,k)))piv=i;end do
            if(abs(aa(piv,k))<=epsilon(1.0_dp)*max(1.0_dp,maxval(abs(aa))))then;istat=1;exit;end if
            if(piv/=k)then;row=aa(k,:);aa(k,:)=aa(piv,:);aa(piv,:)=row;tmp=bb(k);bb(k)=bb(piv);bb(piv)=tmp;end if
            do i=k+1,n
                fac=aa(i,k)/aa(k,k);aa(i,k:n)=aa(i,k:n)-fac*aa(k,k:n);bb(i)=bb(i)-fac*bb(k)
            end do
        end do
        if(istat==0 .and. abs(aa(n,n))<=epsilon(1.0_dp)*max(1.0_dp,maxval(abs(aa))))istat=1
        if(istat/=0)then;x=0.0_dp;if(present(info))info=istat;return;end if
        x(n)=bb(n)/aa(n,n)
        do i=n-1,1,-1;x(i)=(bb(i)-dot_product(aa(i,i+1:n),x(i+1:n)))/aa(i,i);end do
        if(present(info))info=0
    end subroutine solve_linear

    subroutine least_squares(xmat,y,beta,info,ridge)
        real(dp),intent(in)::xmat(:,:),y(:)
        real(dp),intent(out)::beta(:)
        integer,intent(out),optional::info
        real(dp),intent(in),optional::ridge
        real(dp),allocatable::a(:,:),b(:)
        real(dp)::lam
        integer::p,i,istat
        p=size(xmat,2);lam=0.0_dp;if(present(ridge))lam=ridge
        allocate(a(p,p),b(p));a=matmul(transpose(xmat),xmat);b=matmul(transpose(xmat),y)
        if(lam>0.0_dp)then;do i=1,p;a(i,i)=a(i,i)+lam;end do;end if
        call solve_linear(a,b,beta,istat)
        if(present(info))info=istat
    end subroutine least_squares

    subroutine invert_matrix(a,ainv,info)
        real(dp),intent(in)::a(:,:)
        real(dp),intent(out)::ainv(:,:)
        integer,intent(out),optional::info
        real(dp),allocatable::e(:),x(:)
        integer::n,j,istat
        n=size(a,1);allocate(e(n),x(n));ainv=0.0_dp;istat=0
        do j=1,n
            e=0.0_dp;e(j)=1.0_dp;call solve_linear(a,e,x,istat)
            if(istat/=0)exit
            ainv(:,j)=x
        end do
        if(present(info))info=istat
    end subroutine invert_matrix

    pure function covariance_matrix(x) result(c)
        real(dp),intent(in)::x(:,:)
        real(dp)::c(size(x,2),size(x,2)),m(size(x,2))
        integer::n,i
        n=size(x,1);m=sum(x,dim=1)/real(n,dp);c=0.0_dp
        do i=1,n;c=c+spread(x(i,:)-m,2,size(m))*spread(x(i,:)-m,1,size(m));end do
        c=c/real(max(1,n-1),dp)
    end function covariance_matrix

    pure real(dp) function pearson_cor(x,y) result(r)
        real(dp),intent(in)::x(:),y(:)
        real(dp)::xm,ym,sx,sy
        xm=sum(x)/real(size(x),dp);ym=sum(y)/real(size(y),dp)
        sx=sqrt(sum((x-xm)**2));sy=sqrt(sum((y-ym)**2))
        if(sx<=tiny(1.0_dp).or.sy<=tiny(1.0_dp))then;r=0.0_dp;else;r=sum((x-xm)*(y-ym))/(sx*sy);end if
    end function pearson_cor

end module greybox_linalg
