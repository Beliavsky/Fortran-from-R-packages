module greybox_utils
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use greybox_kinds, only: dp
    use greybox_special, only: normal_quantile
    implicit none
    private
    public :: polyprod, backshift, xreg_expander, xreg_multiplier, xreg_transformer
    public :: dyn_mult_calc, outlier_dummy, temporal_dummy

contains

    pure function polyprod(x,y) result(z)
        real(dp),intent(in)::x(:),y(:)
        real(dp)::z(size(x)+size(y)-1)
        integer::i,j
        z=0.0_dp
        do i=1,size(x);do j=1,size(y);z(i+j-1)=z(i+j-1)+x(i)*y(j);end do;end do
    end function polyprod

    pure function backshift(x,k,gap_mode) result(y)
        real(dp),intent(in)::x(:)
        integer,intent(in)::k
        integer,intent(in),optional::gap_mode ! 0 zero, 1 naive, 2 NaN
        real(dp)::y(size(x)),fill_left,fill_right
        integer::i,gm,j
        gm=1;if(present(gap_mode))gm=gap_mode
        select case(gm)
        case(0);fill_left=0.0_dp;fill_right=0.0_dp
        case(2);fill_left=ieee_value(0.0_dp,ieee_quiet_nan);fill_right=fill_left
        case default;fill_left=x(1);fill_right=x(size(x))
        end select
        if(k==0)then;y=x;return;end if
        do i=1,size(x)
            j=i-k
            if(j<1)then;y(i)=fill_left;else if(j>size(x))then;y(i)=fill_right;else;y(i)=x(j);end if
        end do
    end function backshift

    function xreg_expander(x,lags,gap_mode) result(out)
        real(dp),intent(in)::x(:,:)
        integer,intent(in)::lags(:)
        integer,intent(in),optional::gap_mode
        real(dp),allocatable::out(:,:)
        integer::n,p,nl,i,j,col
        n=size(x,1);p=size(x,2);nl=count(lags/=0)
        allocate(out(n,p*(nl+1)));col=0
        do j=1,p
            col=col+1;out(:,col)=x(:,j)
            do i=1,size(lags)
                if(lags(i)==0)cycle
                col=col+1
                ! R xregExpander: negative = lag, positive = lead. backshift positive = lag.
                if(present(gap_mode))then
                    out(:,col)=backshift(x(:,j),-lags(i),gap_mode)
                else
                    out(:,col)=backshift(x(:,j),-lags(i))
                end if
            end do
        end do
    end function xreg_expander

    function xreg_multiplier(x) result(out)
        real(dp),intent(in)::x(:,:)
        real(dp),allocatable::out(:,:)
        integer::n,p,nc,i,j,k
        n=size(x,1);p=size(x,2);nc=p*(p-1)/2
        allocate(out(n,p+nc));out(:,1:p)=x;k=p
        do i=1,p-1;do j=i+1,p;k=k+1;out(:,k)=x(:,i)*x(:,j);end do;end do
    end function xreg_multiplier

    function xreg_transformer(x,codes) result(out)
        real(dp),intent(in)::x(:,:)
        integer,intent(in)::codes(:) ! 1 log, 2 exp, 3 inv, 4 sqrt, 5 square
        real(dp),allocatable::out(:,:)
        integer::n,p,nf,i,j,k
        n=size(x,1);p=size(x,2);nf=size(codes);allocate(out(n,p*(nf+1)));k=0
        do j=1,p
            k=k+1;out(:,k)=x(:,j)
            do i=1,nf
                k=k+1
                select case(codes(i))
                case(1);out(:,k)=log(x(:,j))
                case(2);out(:,k)=exp(x(:,j))
                case(3);out(:,k)=1.0_dp/x(:,j)
                case(4);out(:,k)=sqrt(x(:,j))
                case(5);out(:,k)=x(:,j)**2
                case default;out(:,k)=x(:,j)
                end select
            end do
        end do
    end function xreg_transformer

    pure function dyn_mult_calc(phi,beta,h) result(c)
        real(dp),intent(in)::phi(:),beta(:)
        integer,intent(in)::h
        real(dp)::c(max(0,h))
        integer::s,i,p
        c=0.0_dp;if(h<1)return;p=size(phi);c(1)=beta(1)
        do s=1,h-1
            if(s+1<=size(beta))c(s+1)=beta(s+1)
            do i=1,min(p,s);c(s+1)=c(s+1)+phi(i)*c(s-i+1);end do
        end do
    end function dyn_mult_calc

    pure function outlier_dummy(residuals,level) result(dummy)
        real(dp),intent(in)::residuals(:),level
        integer::dummy(size(residuals))
        real(dp)::m,s,z
        integer::n
        n=size(residuals);m=sum(residuals)/real(n,dp)
        s=sqrt(sum((residuals-m)**2)/real(max(1,n-1),dp))
        z=normal_quantile(0.5_dp+0.5_dp*level,0.0_dp,1.0_dp)
        dummy=merge(1,0,abs(residuals-m)>z*s)
    end function outlier_dummy

    pure function temporal_dummy(n,frequency,drop_first) result(d)
        integer,intent(in)::n,frequency
        logical,intent(in),optional::drop_first
        real(dp),allocatable::d(:,:)
        logical::drop
        integer::p,i,k,start
        drop=.false.;if(present(drop_first))drop=drop_first
        p=frequency;if(drop)p=p-1
        allocate(d(n,max(0,p)));d=0.0_dp
        start=merge(2,1,drop)
        do i=1,n
            k=mod(i-1,frequency)+1
            if(k>=start .and. k-start+1<=p)d(i,k-start+1)=1.0_dp
        end do
    end function temporal_dummy

end module greybox_utils
