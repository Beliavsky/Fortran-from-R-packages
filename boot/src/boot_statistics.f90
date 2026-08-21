module boot_statistics
    use boot_kinds, only : dp
    use boot_special, only : logistic_cdf, logit_fn
    implicit none
    private
    public :: weighted_corr, var_linear, k3_linear, cum3, logit, inv_logit
    public :: mean_dp, variance_dp, covariance_dp, sort_real, sample_quantile

contains
    real(dp) function mean_dp(x) result(m)
        real(dp),intent(in)::x(:)
        if(size(x)==0) error stop "mean_dp: empty vector"
        m=sum(x)/real(size(x),dp)
    end function

    real(dp) function variance_dp(x) result(v)
        real(dp),intent(in)::x(:)
        real(dp)::m
        if(size(x)<2) then
        v=0.0_dp
        return
        end if
        m=mean_dp(x)
        v=sum((x-m)**2)/real(size(x)-1,dp)
    end function

    real(dp) function covariance_dp(x,y) result(v)
        real(dp),intent(in)::x(:),y(:)
        real(dp)::mx,my
        if(size(x)/=size(y) .or. size(x)<2) error stop "covariance_dp: invalid sizes"
        mx=mean_dp(x)
        my=mean_dp(y)
        v=sum((x-mx)*(y-my))/real(size(x)-1,dp)
    end function

    real(dp) function weighted_corr(d,w) result(rho)
        real(dp),intent(in)::d(:,:),w(:)
        real(dp)::s,m1,m2,v1,v2,c12
        if(size(d,2)/=2 .or. size(d,1)/=size(w)) error stop "weighted_corr: size mismatch"
        s=sum(w)
        if(s<=0.0_dp) error stop "weighted_corr: nonpositive weight sum"
        m1=sum(d(:,1)*w)/s
        m2=sum(d(:,2)*w)/s
        v1=sum(d(:,1)**2*w)/s-m1*m1
        v2=sum(d(:,2)**2*w)/s-m2*m2
        c12=sum(d(:,1)*d(:,2)*w)/s-m1*m2
        rho=c12/sqrt(v1*v2)
    end function

    real(dp) function var_linear(l,strata) result(v)
        real(dp),intent(in)::l(:)
        integer,intent(in),optional::strata(:)
        integer::i,j,n,ns
        logical::first
        v=0.0_dp
        n=size(l)
        if(.not.present(strata)) then
            v=sum(l*l)/real(n*n,dp)
            return
        end if
        if(size(strata)/=n) error stop "var_linear: strata mismatch"
        do i=1,n
            first=.true.
            do j=1,i-1
                if(strata(j)==strata(i)) first=.false.
            end do
            if(.not.first) cycle
            ns=count(strata==strata(i))
            do j=1,n
                if(strata(j)==strata(i)) v=v+l(j)*l(j)/real(ns*ns,dp)
            end do
        end do
    end function

    real(dp) function k3_linear(l,strata) result(k3)
        real(dp),intent(in)::l(:)
        integer,intent(in),optional::strata(:)
        integer::i,j,n,ns
        logical::first
        k3=0.0_dp
        n=size(l)
        if(.not.present(strata)) then
            k3=sum(l**3)/real(n**3,dp)
            return
        end if
        do i=1,n
            first=.true.
            do j=1,i-1
                if(strata(j)==strata(i)) first=.false.
            end do
            if(.not.first) cycle
            ns=count(strata==strata(i))
            do j=1,n
                if(strata(j)==strata(i)) k3=k3+l(j)**3/real(ns**3,dp)
            end do
        end do
    end function

    real(dp) function cum3(a,b,c,unbiased) result(k)
        real(dp),intent(in)::a(:),b(:),c(:)
        logical,intent(in),optional::unbiased
        logical::ub
        integer::n
        real(dp)::mult
        n=size(a)
        if(size(b)/=n .or. size(c)/=n) error stop "cum3: size mismatch"
        ub=.true.
        if(present(unbiased))ub=unbiased
        if(ub) then
            if(n<3) error stop "cum3: n < 3"
            mult=real(n,dp)/real((n-1)*(n-2),dp)
        else
            mult=1.0_dp/real(n,dp)
        end if
        k=mult*sum((a-mean_dp(a))*(b-mean_dp(b))*(c-mean_dp(c)))
    end function

    elemental real(dp) function logit(p) result(x)
        real(dp),intent(in)::p
        x=logit_fn(p)
    end function

    elemental real(dp) function inv_logit(x) result(p)
        real(dp),intent(in)::x
        p=logistic_cdf(x)
    end function

    subroutine sort_real(x)
        real(dp),intent(inout)::x(:)
        integer::i,j
        real(dp)::tmp
        do i=2,size(x)
            tmp=x(i)
            j=i-1
            do while(j>=1)
                if(x(j)<=tmp)exit
                x(j+1)=x(j)
                j=j-1
            end do
            x(j+1)=tmp
        end do
    end subroutine

    real(dp) function sample_quantile(x,p) result(q)
        real(dp),intent(in)::x(:),p
        real(dp),allocatable::z(:)
        real(dp)::h,f
        integer::i,n
        n=size(x)
        allocate(z(n))
        z=x
        call sort_real(z)
        if(p<=0.0_dp)then
        q=z(1)
        else if(p>=1.0_dp)then
        q=z(n)
        else
            h=1.0_dp+real(n-1,dp)*p
            i=floor(h)
            f=h-real(i,dp)
            if(i>=n)then
            q=z(n)
            else
            q=z(i)+f*(z(i+1)-z(i))
            end if
        end if
    end function
end module boot_statistics
