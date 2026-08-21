module boot_influence
    use boot_kinds, only : dp
    use boot_linalg, only : least_squares
    implicit none
    private
    public :: infinitesimal_jackknife, delete1_jackknife, positive_jackknife
    public :: regression_influence, linear_approximation
    abstract interface
        function weighted_stat(data,weights) result(value)
            import dp
            real(dp),intent(in)::data(:,:),weights(:)
            real(dp)::value
        end function weighted_stat
    end interface
contains
    subroutine infinitesimal_jackknife(data,statistic,strata,l,eps)
        real(dp),intent(in)::data(:,:)
        procedure(weighted_stat)::statistic
        integer,intent(in)::strata(:)
        real(dp),intent(out)::l(size(data,1))
        real(dp),intent(in),optional::eps
        integer::n,i,j,ns
        real(dp)::e,tobs
        real(dp),allocatable::w0(:),w(:)
        n=size(data,1)
        allocate(w0(n),w(n))
        e=0.001_dp/real(n,dp)
        if(present(eps))e=eps/real(n,dp)
        do i=1,n
        ns=count(strata==strata(i))
        w0(i)=1.0_dp/real(ns,dp)
        end do
        tobs=statistic(data,w0)
        do i=1,n
            w=w0
            do j=1,n
            if(strata(j)==strata(i))w(j)=(1.0_dp-e)*w(j)
            end do
            w(i)=w(i)+e
            l(i)=(statistic(data,w)-tobs)/e
        end do
    end subroutine infinitesimal_jackknife

    subroutine delete1_jackknife(data,statistic,strata,l)
        real(dp),intent(in)::data(:,:)
        procedure(weighted_stat)::statistic
        integer,intent(in)::strata(:)
        real(dp),intent(out)::l(size(data,1))
        integer::n,i,j,ns
        real(dp)::tobs,s
        real(dp),allocatable::w0(:),w(:)
        n=size(data,1)
        allocate(w0(n),w(n))
        do i=1,n
        ns=count(strata==strata(i))
        w0(i)=1.0_dp/real(ns,dp)
        end do
        tobs=statistic(data,w0)
        do i=1,n
            w=w0
            w(i)=0.0_dp
            s=sum(w,mask=strata==strata(i))
            do j=1,n
            if(strata(j)==strata(i))w(j)=w(j)/s
            end do
            ns=count(strata==strata(i))
            l(i)=real(ns-1,dp)*(tobs-statistic(data,w))
        end do
    end subroutine delete1_jackknife

    subroutine positive_jackknife(data,statistic,strata,l)
        real(dp),intent(in)::data(:,:)
        procedure(weighted_stat)::statistic
        integer,intent(in)::strata(:)
        real(dp),intent(out)::l(size(data,1))
        integer::n,i,j,ns
        real(dp)::tobs,s
        real(dp),allocatable::w0(:),w(:)
        n=size(data,1)
        allocate(w0(n),w(n))
        do i=1,n
        ns=count(strata==strata(i))
        w0(i)=1.0_dp/real(ns,dp)
        end do
        tobs=statistic(data,w0)
        do i=1,n
            ns=count(strata==strata(i))
            w=w0
            do j=1,n
            if(strata(j)==strata(i))w(j)=1.0_dp/real(ns+1,dp)
            end do
            w(i)=2.0_dp/real(ns+1,dp)
            s=sum(w,mask=strata==strata(i))
            do j=1,n
            if(strata(j)==strata(i))w(j)=w(j)/s
            end do
            l(i)=real(ns+1,dp)*(statistic(data,w)-tobs)
        end do
    end subroutine positive_jackknife

    subroutine regression_influence(freq,t,strata,l)
        integer,intent(in)::freq(:,:)
        real(dp),intent(in)::t(:)
        integer,intent(in)::strata(:)
        real(dp),intent(out)::l(size(freq,2))
        integer::r,n,i,j,p,k,info
        logical::first
        integer,allocatable::inc(:)
        real(dp),allocatable::x(:,:),beta(:),y(:)
        r=size(freq,1)
        n=size(freq,2)
        if(size(t)/=r .or. size(strata)/=n)error stop "regression_influence: size mismatch"
        allocate(inc(n))
        p=0
        do i=1,n
            first=.true.
            do j=1,i-1
            if(strata(j)==strata(i))first=.false.
            end do
            if(.not.first)then
            p=p+1
            inc(p)=i
            end if
        end do
        allocate(x(r,p+1),beta(p+1),y(r))
        x(:,1)=1.0_dp
        y=t
        do k=1,p
            i=inc(k)
            x(:,k+1)=real(freq(:,i),dp)/real(count(strata==strata(i)),dp)
        end do
        call least_squares(x,y,beta,info,1.0e-12_dp)
        if(info/=0)error stop "regression_influence: singular regression"
        l=0.0_dp
        do k=1,p
        l(inc(k))=beta(k+1)
        end do
        do i=1,n
        l(i)=l(i)-sum(l,mask=strata==strata(i))/real(count(strata==strata(i)),dp)
        end do
    end subroutine regression_influence

    subroutine linear_approximation(freq,l,strata,t0,tlin)
        integer,intent(in)::freq(:,:)
        real(dp),intent(in)::l(:),t0
        integer,intent(in)::strata(:)
        real(dp),intent(out)::tlin(size(freq,1))
        integer::r,i,n
        n=size(freq,2)
        if(size(l)/=n .or. size(strata)/=n)error stop "linear_approximation: size mismatch"
        do r=1,size(freq,1)
            tlin(r)=t0
            do i=1,n
                tlin(r)=tlin(r)+real(freq(r,i),dp)*l(i)/real(count(strata==strata(i)),dp)
            end do
        end do
    end subroutine linear_approximation
end module boot_influence
