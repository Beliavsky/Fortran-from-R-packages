module boot_ci
    use boot_kinds, only : dp
    use boot_special, only : normal_cdf, normal_quantile
    use boot_statistics, only : mean_dp, variance_dp, sort_real
    implicit none
    private
    public :: norm_inter, normal_ci, basic_ci, percentile_ci, studentized_ci, bca_ci, abc_ci

    abstract interface
        function weighted_scalar_stat(data,weights) result(value)
            import dp
            real(dp),intent(in)::data(:,:),weights(:)
            real(dp)::value
        end function weighted_scalar_stat
    end interface
contains
    subroutine norm_inter(t,alpha,quantile,ranks)
        real(dp),intent(in)::t(:),alpha(:)
        real(dp),intent(out)::quantile(size(alpha))
        real(dp),intent(out),optional::ranks(size(alpha))
        real(dp),allocatable::z(:)
        real(dp)::rk,zp,z0,z1,frac
        integer::r,k,i
        r=size(t)
        allocate(z(r))
        z=t
        call sort_real(z)
        do i=1,size(alpha)
            rk=real(r+1,dp)*alpha(i)
            if(present(ranks))ranks(i)=rk
            k=int(rk)
            if(k<=0)then
                quantile(i)=z(1)
            else if(k>=r)then
                quantile(i)=z(r)
            else if(abs(rk-real(k,dp))<1.0e-14_dp)then
                quantile(i)=z(k)
            else
                zp=normal_quantile(alpha(i))
                z0=normal_quantile(real(k,dp)/real(r+1,dp))
                z1=normal_quantile(real(k+1,dp)/real(r+1,dp))
                frac=(zp-z0)/(z1-z0)
                quantile(i)=z(k)+frac*(z(k+1)-z(k))
            end if
        end do
    end subroutine norm_inter

    subroutine normal_ci(t0,t,conf,lo,hi,var_t0)
        real(dp),intent(in)::t0,t(:),conf(:)
        real(dp),intent(out)::lo(size(conf)),hi(size(conf))
        real(dp),intent(in),optional::var_t0
        real(dp)::v,bias,merr
        integer::i
        if(present(var_t0))then
        v=var_t0
        else
        v=variance_dp(t)
        end if
        bias=mean_dp(t)-t0
        do i=1,size(conf)
            merr=sqrt(v)*normal_quantile((1.0_dp+conf(i))/2.0_dp)
            lo(i)=t0-bias-merr
            hi(i)=t0-bias+merr
        end do
    end subroutine normal_ci

    subroutine basic_ci(t0,t,conf,lo,hi)
        real(dp),intent(in)::t0,t(:),conf(:)
        real(dp),intent(out)::lo(size(conf)),hi(size(conf))
        real(dp)::alpha(2*size(conf)),q(2*size(conf))
        integer::i
        do i=1,size(conf)
            alpha(2*i-1)=(1.0_dp+conf(i))/2.0_dp
            alpha(2*i)=(1.0_dp-conf(i))/2.0_dp
        end do
        call norm_inter(t,alpha,q)
        do i=1,size(conf)
        lo(i)=2*t0-q(2*i-1)
        hi(i)=2*t0-q(2*i)
        end do
    end subroutine basic_ci

    subroutine percentile_ci(t,conf,lo,hi)
        real(dp),intent(in)::t(:),conf(:)
        real(dp),intent(out)::lo(size(conf)),hi(size(conf))
        real(dp)::alpha(2*size(conf)),q(2*size(conf))
        integer::i
        do i=1,size(conf)
            alpha(2*i-1)=(1.0_dp-conf(i))/2.0_dp
            alpha(2*i)=(1.0_dp+conf(i))/2.0_dp
        end do
        call norm_inter(t,alpha,q)
        do i=1,size(conf)
        lo(i)=q(2*i-1)
        hi(i)=q(2*i)
        end do
    end subroutine percentile_ci

    subroutine studentized_ci(t0,var0,t,var_t,conf,lo,hi)
        real(dp),intent(in)::t0,var0,t(:),var_t(:),conf(:)
        real(dp),intent(out)::lo(size(conf)),hi(size(conf))
        real(dp),allocatable::z(:),alpha(:),q(:)
        integer::i,m
        if(size(t)/=size(var_t))error stop "studentized_ci: size mismatch"
        m=size(conf)
        allocate(z(size(t)),alpha(2*m),q(2*m))
        z=(t-t0)/sqrt(var_t)
        do i=1,m
            alpha(2*i-1)=(1.0_dp+conf(i))/2.0_dp
            alpha(2*i)=(1.0_dp-conf(i))/2.0_dp
        end do
        call norm_inter(z,alpha,q)
        do i=1,m
        lo(i)=t0-sqrt(var0)*q(2*i-1)
        hi(i)=t0-sqrt(var0)*q(2*i)
        end do
    end subroutine studentized_ci

    subroutine bca_ci(t0,t,influence,conf,lo,hi,adjusted_alpha)
        real(dp),intent(in)::t0,t(:),influence(:),conf(:)
        real(dp),intent(out)::lo(size(conf)),hi(size(conf))
        real(dp),intent(out),optional::adjusted_alpha(2*size(conf))
        real(dp)::w,a,zalpha,aa(2*size(conf)),q(2*size(conf)),p_less
        integer::i
        p_less=real(count(t<t0),dp)/real(size(t),dp)
        if(p_less<=0.0_dp .or. p_less>=1.0_dp) error stop "bca_ci: infinite bias adjustment"
        w=normal_quantile(p_less)
        a=sum(influence**3)/(6.0_dp*sum(influence**2)**1.5_dp)
        do i=1,size(conf)
            zalpha=normal_quantile((1.0_dp-conf(i))/2.0_dp)
            aa(2*i-1)=normal_cdf(w+(w+zalpha)/(1.0_dp-a*(w+zalpha)))
            zalpha=normal_quantile((1.0_dp+conf(i))/2.0_dp)
            aa(2*i)=normal_cdf(w+(w+zalpha)/(1.0_dp-a*(w+zalpha)))
        end do
        call norm_inter(t,aa,q)
        do i=1,size(conf)
        lo(i)=q(2*i-1)
        hi(i)=q(2*i)
        end do
        if(present(adjusted_alpha))adjusted_alpha=aa
    end subroutine bca_ci

    subroutine abc_ci(data,statistic,strata,conf,lo,hi,eps)
        real(dp),intent(in)::data(:,:)
        procedure(weighted_scalar_stat)::statistic
        integer,intent(in)::strata(:)
        real(dp),intent(in)::conf(:)
        real(dp),intent(out)::lo(size(conf)),hi(size(conf))
        real(dp),intent(in),optional::eps
        integer::n,i,j,ns
        real(dp)::e,t0,t1,t2,temp1,sigmahat,ahat,bhat,chat,bprime,z,lalpha
        real(dp),allocatable::w0(:),w1(:),w2(:),l(:),l2(:),dhat(:),w3(:),w4(:),wf(:)
        n=size(data,1)
        if(size(strata)/=n)error stop "abc_ci: strata mismatch"
        e=0.001_dp/real(n,dp)
        if(present(eps))e=eps
        allocate(w0(n),w1(n),w2(n),l(n),l2(n),dhat(n),w3(n),w4(n),wf(n))
        do i=1,n
        ns=count(strata==strata(i))
        w0(i)=1.0_dp/real(ns,dp)
        end do
        t0=statistic(data,w0)
        do i=1,n
            w1=w0
            w2=w0
            do j=1,n
                if(strata(j)==strata(i))then
                    w1(j)=(1.0_dp-e)*w0(j)
                    w2(j)=(1.0_dp+e)*w0(j)
                end if
            end do
            w1(i)=w1(i)+e
            w2(i)=w2(i)-e
            t1=statistic(data,w1)
            t2=statistic(data,w2)
            l(i)=(t1-t2)/(2.0_dp*e)
            l2(i)=(t1-2.0_dp*t0+t2)/(e*e)
        end do
        temp1=sum(l*l)
        sigmahat=sqrt(temp1)/real(n,dp)
        ahat=sum(l**3)/(6.0_dp*temp1**1.5_dp)
        bhat=sum(l2)/(2.0_dp*real(n*n,dp))
        dhat=l/(real(n*n,dp)*sigmahat)
        w3=w0+e*dhat
        w4=w0-e*dhat
        call renorm_by_strata(w3,strata)
        call renorm_by_strata(w4,strata)
        chat=(statistic(data,w3)-2.0_dp*t0+statistic(data,w4))/(2.0_dp*e*e*sigmahat)
        bprime=ahat-(bhat/sigmahat-chat)
        do i=1,size(conf)
            z=normal_quantile((1.0_dp-conf(i))/2.0_dp)
            lalpha=(bprime+z)/(1.0_dp-ahat*(bprime+z))**2
            wf=w0+lalpha*dhat
            call renorm_by_strata(wf,strata)
            lo(i)=statistic(data,wf)
            z=normal_quantile((1.0_dp+conf(i))/2.0_dp)
            lalpha=(bprime+z)/(1.0_dp-ahat*(bprime+z))**2
            wf=w0+lalpha*dhat
            call renorm_by_strata(wf,strata)
            hi(i)=statistic(data,wf)
        end do
    end subroutine abc_ci

    subroutine renorm_by_strata(w,strata)
        real(dp),intent(inout)::w(:)
        integer,intent(in)::strata(:)
        integer::i
        real(dp)::s
        do i=1,size(w)
            s=sum(w,mask=strata==strata(i))
            w(i)=w(i)/s
        end do
    end subroutine renorm_by_strata
end module boot_ci
