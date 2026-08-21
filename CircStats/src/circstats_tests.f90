module circstats_tests
    use circstats_kinds, only: dp, pi, twopi
    use circstats_types, only: test_result, rao_homogeneity_result
    use circstats_utils, only: wrap_2pi, sort_real, sample_variance, sample_covariance
    use circstats_core, only: circ_mean, est_kappa, deg
    use circstats_distributions, only: pvm
    use circstats_special, only: normal_pdf, normal_cdf, chi_square_cdf, chi_square_quantile
    use circstats_rao_table, only: rao_critical_table
    implicit none
    private
    public :: rayleigh_test, v0_test, kuiper_test, watson_uniform_test
    public :: watson_vm_test, watson_two_test, rao_spacing_test, rao_homogeneity_test

contains

    function rayleigh_test(x,degree) result(res)
        real(dp), intent(in) :: x(:)
        logical, intent(in), optional :: degree
        type(test_result) :: res
        real(dp), allocatable :: w(:)
        real(dp) :: ss,cc,rbar,z,temp
        integer :: n
        logical :: use_degree
        n=size(x)
        allocate(w(n))
        w=x
        use_degree=.false.
        if (present(degree)) use_degree=degree
        if (use_degree) w=w*pi/180.0_dp
        ss=sum(sin(w))
        cc=sum(cos(w))
        rbar=hypot(ss,cc)/real(n,dp)
        z=real(n,dp)*rbar*rbar
        temp=1.0_dp
        if (n < 50) then
            temp=1.0_dp+(2.0_dp*z-z*z)/(4.0_dp*real(n,dp)) - &
                (24.0_dp*z-132.0_dp*z*z+76.0_dp*z**3-9.0_dp*z**4)/(288.0_dp*real(n*n,dp))
        end if
        res%statistic=rbar
        res%p_value=max(0.0_dp,min(1.0_dp,exp(-z)*temp))
    end function rayleigh_test

    function v0_test(x,mu0,degree) result(res)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in), optional :: mu0
        logical, intent(in), optional :: degree
        type(test_result) :: res
        real(dp), allocatable :: w(:)
        real(dp) :: m,r0,z,pz,fz
        integer :: n
        logical :: use_degree
        n=size(x)
        allocate(w(n))
        w=x
        m=0.0_dp
        if (present(mu0)) m=mu0
        use_degree=.false.
        if (present(degree)) use_degree=degree
        if (use_degree) then
            w=w*pi/180.0_dp
            m=m*pi/180.0_dp
        end if
        r0=sum(cos(w-m))/real(n,dp)
        z=sqrt(2.0_dp*real(n,dp))*r0
        pz=normal_cdf(z)
        fz=normal_pdf(z)
        res%statistic=r0
        res%p_value=1.0_dp-pz+fz*((3.0_dp*z-z**3)/(16.0_dp*real(n,dp)) + &
            (15.0_dp*z+305.0_dp*z**3-125.0_dp*z**5+9.0_dp*z**7)/(4608.0_dp*real(n*n,dp)))
        res%p_value=max(0.0_dp,min(1.0_dp,res%p_value))
    end function v0_test

    function kuiper_test(x,alpha) result(res)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in), optional :: alpha
        type(test_result) :: res
        real(dp), parameter :: levels(5)=[0.15_dp,0.10_dp,0.05_dp,0.025_dp,0.01_dp]
        real(dp), parameter :: crits(5)=[1.537_dp,1.620_dp,1.747_dp,1.862_dp,2.001_dp]
        real(dp), allocatable :: w(:)
        real(dp) :: dpv,dmv,a
        integer :: n,i,j
        n=size(x)
        allocate(w(n))
        w=wrap_2pi(x)/twopi
        call sort_real(w)
        dpv=0.0_dp
        dmv=0.0_dp
        do i=1,n
            dpv=max(dpv,real(i,dp)/real(n,dp)-w(i))
            dmv=max(dmv,w(i)-real(i-1,dp)/real(n,dp))
        end do
        res%statistic=(dpv+dmv)*(sqrt(real(n,dp))+0.155_dp+0.24_dp/sqrt(real(n,dp)))
        if (present(alpha)) then
            a=alpha
            if (a > 0.0_dp) then
                j=0
                do i=1,5
                    if (abs(a-levels(i)) < 1.0e-12_dp) j=i
                end do
                if (j == 0) error stop "kuiper_test: invalid alpha"
                res%critical=crits(j)
                res%reject=res%statistic>res%critical
                return
            end if
        end if
        call bracket_from_crit(res%statistic,crits,levels,res%p_lower,res%p_upper)
    end function kuiper_test

    function watson_uniform_test(x,alpha) result(res)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in), optional :: alpha
        type(test_result) :: res
        real(dp), parameter :: levels(4)=[0.01_dp,0.025_dp,0.05_dp,0.10_dp]
        real(dp), parameter :: crits(4)=[0.267_dp,0.221_dp,0.187_dp,0.152_dp]
        real(dp), allocatable :: u(:)
        real(dp) :: ubar,a
        integer :: n,i,j
        n=size(x)
        allocate(u(n))
        u=wrap_2pi(x)/twopi
        call sort_real(u)
        ubar=sum(u)/real(n,dp)
        res%statistic=0.0_dp
        do i=1,n
            res%statistic=res%statistic+(u(i)-ubar-real(2*i-1,dp)/(2.0_dp*real(n,dp))+0.5_dp)**2
        end do
        res%statistic=res%statistic+1.0_dp/(12.0_dp*real(n,dp))
        res%statistic=(res%statistic-0.1_dp/real(n,dp)+0.1_dp/real(n*n,dp))*(1.0_dp+0.8_dp/real(n,dp))
        if (present(alpha)) then
            a=alpha
            if (a > 0.0_dp) then
                j=0
                do i=1,4
                    if (abs(a-levels(i)) < 1.0e-12_dp) j=i
                end do
                if (j == 0) error stop "watson_uniform_test: invalid alpha"
                res%critical=crits(j)
                res%reject=res%statistic>res%critical
                return
            end if
        end if
        call bracket_from_crit(res%statistic,crits,levels,res%p_lower,res%p_upper)
    end function watson_uniform_test

    function watson_vm_test(x,alpha) result(res)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in), optional :: alpha
        type(test_result) :: res
        real(dp), parameter :: levels(3)=[0.10_dp,0.05_dp,0.01_dp]
        real(dp), parameter :: crits(7,3)=reshape([ &
            0.052_dp,0.056_dp,0.066_dp,0.077_dp,0.084_dp,0.093_dp,0.096_dp, &
            0.061_dp,0.066_dp,0.079_dp,0.092_dp,0.101_dp,0.113_dp,0.117_dp, &
            0.081_dp,0.090_dp,0.110_dp,0.128_dp,0.142_dp,0.158_dp,0.164_dp], [7,3])
        real(dp), allocatable :: z(:)
        real(dp) :: mu,kappa,zbar,a
        integer :: n,i,row,j
        n=size(x)
        mu=circ_mean(x)
        kappa=est_kappa(x)
        allocate(z(n))
        do i=1,n
            z(i)=pvm(wrap_2pi(x(i)-mu),0.0_dp,kappa)
        end do
        call sort_real(z)
        zbar=sum(z)/real(n,dp)
        res%statistic=0.0_dp
        do i=1,n
            res%statistic=res%statistic+(z(i)-real(2*i-1,dp)/(2.0_dp*real(n,dp)))**2
        end do
        res%statistic=res%statistic-real(n,dp)*(zbar-0.5_dp)**2+1.0_dp/(12.0_dp*real(n,dp))
        if (kappa < 0.25_dp) then
            row=1
        else if (kappa < 0.75_dp) then
            row=2
        else if (kappa < 1.25_dp) then
            row=3
        else if (kappa < 1.75_dp) then
            row=4
        else if (kappa < 3.0_dp) then
            row=5
        else if (kappa < 5.0_dp) then
            row=6
        else
            row=7
        end if
        if (present(alpha)) then
            a=alpha
            if (a > 0.0_dp) then
                j=0
                do i=1,3
                    if (abs(a-levels(i)) < 1.0e-12_dp) j=i
                end do
                if (j == 0) error stop "watson_vm_test: invalid alpha"
                res%critical=crits(row,j)
                res%reject=res%statistic>res%critical
                return
            end if
        end if
        call bracket_from_crit(res%statistic,crits(row,:),levels,res%p_lower,res%p_upper)
    end function watson_vm_test

    function watson_two_test(x,y,alpha) result(res)
        real(dp), intent(in) :: x(:),y(:)
        real(dp), intent(in), optional :: alpha
        type(test_result) :: res
        real(dp), parameter :: levels(4)=[0.001_dp,0.01_dp,0.05_dp,0.10_dp]
        real(dp), parameter :: crits(4)=[0.385_dp,0.268_dp,0.187_dp,0.152_dp]
        real(dp), allocatable :: values(:),d(:)
        integer, allocatable :: group(:),idx(:)
        real(dp) :: dbar,a,key
        integer :: n1,n2,n,i,j,itmp,ca,cb
        n1=size(x)
        n2=size(y)
        n=n1+n2
        allocate(values(n),group(n),idx(n),d(n))
        values(1:n1)=wrap_2pi(x)
        values(n1+1:n)=wrap_2pi(y)
        group(1:n1)=1
        group(n1+1:n)=2
        do i=1,n
            idx(i)=i
        end do
        do i=2,n
            key=values(idx(i))
            itmp=idx(i)
            j=i-1
            do while (j >= 1)
                if (values(idx(j)) <= key) exit
                idx(j+1)=idx(j)
                j=j-1
            end do
            idx(j+1)=itmp
        end do
        ca=0
        cb=0
        do i=1,n
            if (group(idx(i)) == 1) ca=ca+1
            if (group(idx(i)) == 2) cb=cb+1
            d(i)=real(cb,dp)/real(n2,dp)-real(ca,dp)/real(n1,dp)
        end do
        dbar=sum(d)/real(n,dp)
        res%statistic=real(n1*n2,dp)/real(n*n,dp)*sum((d-dbar)**2)
        if (present(alpha)) then
            a=alpha
            if (a > 0.0_dp) then
                j=0
                do i=1,4
                    if (abs(a-levels(i)) < 1.0e-12_dp) j=i
                end do
                if (j == 0) error stop "watson_two_test: invalid alpha"
                res%critical=crits(j)
                res%reject=res%statistic>res%critical
                return
            end if
        end if
        call bracket_from_crit(res%statistic,crits,levels,res%p_lower,res%p_upper)
    end function watson_two_test

    function rao_spacing_test(x,alpha,radians) result(res)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in), optional :: alpha
        logical, intent(in), optional :: radians
        type(test_result) :: res
        real(dp), parameter :: levels(4)=[0.001_dp,0.01_dp,0.05_dp,0.10_dp]
        real(dp), allocatable :: w(:)
        real(dp) :: expected,a
        integer :: n,i,row,col
        logical :: radian_input
        n=size(x)
        if (n < 4) error stop "rao_spacing_test: sample size too small"
        allocate(w(n))
        radian_input=.true.
        if (present(radians)) radian_input=radians
        if (radian_input) then
            w=modulo(deg(x),360.0_dp)
        else
            w=modulo(x,360.0_dp)
        end if
        call sort_real(w)
        expected=360.0_dp/real(n,dp)
        res%statistic=0.5_dp*abs(w(1)-w(n)+360.0_dp-expected)
        do i=1,n-1
            res%statistic=res%statistic+0.5_dp*abs(w(i+1)-w(i)-expected)
        end do
        row=rao_row(n)
        if (present(alpha)) then
            a=alpha
            if (a > 0.0_dp) then
                col=0
                do i=1,4
                    if (abs(a-levels(i)) < 1.0e-12_dp) col=i
                end do
                if (col == 0) error stop "rao_spacing_test: invalid alpha"
                res%critical=rao_critical_table(row,col)
                res%reject=res%statistic>res%critical
                return
            end if
        end if
        call bracket_from_crit(res%statistic,rao_critical_table(row,:),levels,res%p_lower,res%p_upper)
    end function rao_spacing_test

    function rao_homogeneity_test(data,group,alpha) result(res)
        real(dp), intent(in) :: data(:)
        integer, intent(in) :: group(:)
        real(dp), intent(in), optional :: alpha
        type(rao_homogeneity_result) :: res
        real(dp), allocatable :: cm(:),sm(:),sco(:),sss(:),scs(:),sp(:),sd(:),tanv(:),u(:)
        real(dp), allocatable :: cx(:),sx(:)
        real(dp) :: a
        integer :: k,g,nj,i,j,maxn,df
        k=maxval(group)
        allocate(cm(k),sm(k),sco(k),sss(k),scs(k),sp(k),sd(k),tanv(k),u(k))
        cm=0.0_dp
        sm=0.0_dp
        sco=0.0_dp
        sss=0.0_dp
        scs=0.0_dp
        sp=0.0_dp
        sd=0.0_dp
        maxn=size(data)
        allocate(cx(maxn),sx(maxn))
        do g=1,k
            j=0
            do i=1,size(data)
                if (group(i) == g) then
                    j=j+1
                    cx(j)=cos(data(i))
                    sx(j)=sin(data(i))
                end if
            end do
            nj=j
            if (nj < 2) error stop "rao_homogeneity_test: each group needs at least two observations"
            cm(g)=sum(cx(1:nj))/real(nj,dp)
            sm(g)=sum(sx(1:nj))/real(nj,dp)
            sco(g)=sample_variance(cx(1:nj))
            sss(g)=sample_variance(sx(1:nj))
            scs(g)=sample_covariance(cx(1:nj),sx(1:nj))
            sp(g)=(sss(g)/(cm(g)**2)+(sm(g)**2*sco(g))/(cm(g)**4)- &
                2.0_dp*sm(g)*scs(g)/(cm(g)**3))/real(nj,dp)
            u(g)=cm(g)**2+sm(g)**2
            sd(g)=4.0_dp*(cm(g)**2*sco(g)+sm(g)**2*sss(g)+2.0_dp*cm(g)*sm(g)*scs(g))/real(nj,dp)
            tanv(g)=sm(g)/cm(g)
        end do
        res%polar%statistic=sum(tanv*tanv/sp)-(sum(tanv/sp)**2)/sum(1.0_dp/sp)
        res%dispersion%statistic=sum(u*u/sd)-(sum(u/sd)**2)/sum(1.0_dp/sd)
        df=k-1
        res%polar%df=df
        res%dispersion%df=df
        res%polar%p_value=1.0_dp-chi_square_cdf(res%polar%statistic,df)
        res%dispersion%p_value=1.0_dp-chi_square_cdf(res%dispersion%statistic,df)
        if (present(alpha)) then
            a=alpha
            if (a > 0.0_dp) then
                res%polar%critical=chi_square_quantile(1.0_dp-a,df)
                res%dispersion%critical=res%polar%critical
                res%polar%reject=res%polar%statistic>res%polar%critical
                res%dispersion%reject=res%dispersion%statistic>res%dispersion%critical
            end if
        end if
    end function rao_homogeneity_test

    pure integer function rao_row(n) result(row)
        integer, intent(in) :: n
        if (n <= 30) then
            row=n-3
        else if (n <= 32) then
            row=27
        else if (n <= 37) then
            row=28
        else if (n <= 42) then
            row=29
        else if (n <= 47) then
            row=30
        else if (n <= 62) then
            row=31
        else if (n <= 87) then
            row=32
        else if (n <= 125) then
            row=33
        else if (n <= 175) then
            row=34
        else if (n <= 250) then
            row=35
        else if (n <= 350) then
            row=36
        else if (n <= 450) then
            row=37
        else if (n <= 550) then
            row=38
        else if (n <= 650) then
            row=39
        else if (n <= 750) then
            row=40
        else if (n <= 850) then
            row=41
        else if (n <= 950) then
            row=42
        else
            row=43
        end if
    end function rao_row

    subroutine bracket_from_crit(stat,crits,levels,plo,phi)
        real(dp), intent(in) :: stat,crits(:),levels(:)
        real(dp), intent(out) :: plo,phi
        integer :: i,n
        n=size(crits)
        plo=-1.0_dp
        phi=-1.0_dp
        if (stat > crits(1)) then
            plo=0.0_dp
            phi=levels(1)
            return
        end if
        do i=1,n-1
            if (stat > crits(i+1) .and. stat <= crits(i)) then
                plo=levels(i)
                phi=levels(i+1)
                return
            end if
        end do
        plo=levels(n)
        phi=1.0_dp
    end subroutine bracket_from_crit
end module circstats_tests
