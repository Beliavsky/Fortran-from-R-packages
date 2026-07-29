! SPDX-License-Identifier: GPL-2.0-or-later
module evir_eda
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use evir_kinds, only : dp
    use evir_types, only : xy_result, band_result, records_result, matrix_result, &
        tail_curve_result, gpd_fit_result, evir_ok, evir_invalid_input, &
        evir_domain_error
    use evir_math, only : sort_ascending, sort_descending, ppoints, normal_quantile, &
        harmonic_numbers, cumulative_max, mean_value, safe_nan
    use evir_data, only : find_threshold, block_maxima
    use evir_distributions, only : qgpd, pgpd
    use evir_fitting, only : fit_gpd, gpd_quantile_estimate, gpd_quantile_wald
    implicit none
    private

    public :: empirical_tail, mean_excess, quantile_plot_data, records_analysis
    public :: hill_estimates, extremal_index, shape_stability, quantile_stability
    public :: tail_curve

contains

    function empirical_tail(data) result(out)
        real(dp), intent(in) :: data(:)
        type(xy_result) :: out
        integer :: i, n
        n = size(data)
        if (n == 0) then
            out%status = evir_invalid_input
            allocate(out%x(0), out%y(0))
            return
        end if
        allocate(out%x(n), out%y(n))
        out%x = data
        call sort_ascending(out%x)
        do i = 1, n
            out%y(i) = 1.0_dp - ppoints(i, n)
        end do
    end function empirical_tail

    function mean_excess(data, omit) result(out)
        real(dp), intent(in) :: data(:)
        integer, intent(in), optional :: omit
        type(xy_result) :: out
        real(dp), allocatable :: sorted(:), unique_values(:)
        integer :: n, i, nu, nkeep, om, n_exc

        n = size(data)
        om = 3
        if (present(omit)) om = omit
        if (n < 2 .or. om < 0) then
            out%status = evir_invalid_input
            allocate(out%x(0), out%y(0))
            return
        end if
        sorted = data
        call sort_ascending(sorted)
        allocate(unique_values(n))
        nu = 0
        do i = 1, n
            if (i == 1) then
                nu = nu+1
                unique_values(nu) = sorted(i)
            else if (abs(sorted(i)-sorted(i-1)) > 10.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(sorted(i)))) then
                nu = nu+1
                unique_values(nu) = sorted(i)
            end if
        end do
        nkeep = nu-1-om
        if (nkeep <= 0) then
            out%status = evir_invalid_input
            allocate(out%x(0), out%y(0))
            return
        end if
        allocate(out%x(nkeep), out%y(nkeep))
        do i = 1, nkeep
            out%x(i) = unique_values(i)
            n_exc = count(data > out%x(i))
            if (n_exc > 0) then
                out%y(i) = sum(pack(data-out%x(i),data>out%x(i)))/real(n_exc,dp)
            else
                out%y(i) = safe_nan()
            end if
        end do
    end function mean_excess

    function quantile_plot_data(data, xi, trim_value, threshold) result(out)
        real(dp), intent(in) :: data(:)
        real(dp), intent(in), optional :: xi, trim_value, threshold
        type(xy_result) :: out
        real(dp), allocatable :: work(:)
        real(dp) :: shape
        logical, allocatable :: mask(:)
        integer :: i, n

        shape = 0.0_dp
        if (present(xi)) shape = xi
        allocate(mask(size(data)))
        mask = .true.
        if (present(threshold)) mask = mask .and. data >= threshold
        if (present(trim_value)) mask = mask .and. data < trim_value
        work = pack(data,mask)
        n = size(work)
        if (n == 0) then
            out%status = evir_invalid_input
            allocate(out%x(0),out%y(0))
            return
        end if
        call sort_ascending(work)
        allocate(out%x(n), out%y(n))
        out%x = work
        do i = 1,n
            if (abs(shape) <= 1.0e-12_dp) then
                out%y(i) = -log(1.0_dp-ppoints(i,n))
            else
                out%y(i) = qgpd(ppoints(i,n),shape)
            end if
        end do
    end function quantile_plot_data

    function records_analysis(data) result(out)
        real(dp), intent(in) :: data(:)
        type(records_result) :: out
        real(dp), allocatable :: rec_all(:), h(:), h2(:)
        integer, allocatable :: trial_tmp(:)
        integer :: n, nr, i
        n = size(data)
        if (n == 0) then
            out%status = evir_invalid_input
            allocate(out%number(0),out%trial(0),out%record(0),out%expected(0),out%se(0))
            return
        end if
        allocate(rec_all(n),h(n),h2(n),trial_tmp(n))
        call cumulative_max(data,rec_all)
        call harmonic_numbers(n,h,h2)
        nr = 1
        trial_tmp(1) = 1
        do i = 2,n
            if (rec_all(i) > rec_all(i-1)) then
                nr = nr+1
                trial_tmp(nr) = i
            end if
        end do
        allocate(out%number(nr),out%trial(nr),out%record(nr),out%expected(nr),out%se(nr))
        do i = 1,nr
            out%number(i) = i
            out%trial(i) = trial_tmp(i)
            out%record(i) = rec_all(trial_tmp(i))
            out%expected(i) = h(trial_tmp(i))
            out%se(i) = sqrt(max(h(trial_tmp(i))-h2(trial_tmp(i)),0.0_dp))
        end do
    end function records_analysis

    function hill_estimates(data, option, start_k, end_k, p, ci) result(out)
        real(dp), intent(in) :: data(:)
        character(len=*), intent(in), optional :: option
        integer, intent(in), optional :: start_k, end_k
        real(dp), intent(in), optional :: p, ci
        type(band_result) :: out
        real(dp), allocatable :: ordered(:), positive(:), logx(:), ave(:), xi_hat(:), alpha_hat(:), yall(:)
        character(len=12) :: opt
        integer :: n, startv, endv, m, i, k
        real(dp) :: prob, civ, z, se

        opt = 'alpha'
        if (present(option)) opt = lowercase(option)
        ordered = data
        call sort_descending(ordered)
        positive = pack(ordered,ordered>0.0_dp)
        n = size(positive)
        startv = 15
        if (present(start_k)) startv = start_k
        endv = n
        if (present(end_k)) endv = min(end_k,n)
        if (n < 2 .or. startv < 2 .or. endv < startv) then
            out%status = evir_invalid_input
            allocate(out%index(0),out%estimate(0),out%lower(0),out%upper(0),out%threshold(0))
            return
        end if
        prob = 0.99_dp
        if (present(p)) prob = p
        civ = 0.95_dp
        if (present(ci)) civ = ci
        allocate(logx(n),ave(n),xi_hat(n),alpha_hat(n),yall(n))
        logx = log(positive)
        ave(1) = logx(1)
        do i = 2,n
            ave(i) = (ave(i-1)*real(i-1,dp)+logx(i))/real(i,dp)
        end do
        xi_hat(1) = safe_nan()
        alpha_hat(1) = safe_nan()
        do i = 2,n
            xi_hat(i) = ave(i)-logx(i)
            if (abs(xi_hat(i)) > epsilon(1.0_dp)) then
                alpha_hat(i) = 1.0_dp/xi_hat(i)
            else
                alpha_hat(i) = huge(1.0_dp)
            end if
        end do
        select case(trim(opt))
        case('alpha')
            yall = alpha_hat
        case('xi')
            yall = xi_hat
        case('quantile')
            do i = 2,n
                yall(i) = positive(i)*(real(n,dp)*(1.0_dp-prob)/real(i,dp))**(-1.0_dp/alpha_hat(i))
            end do
            yall(1) = safe_nan()
        case default
            out%status = evir_invalid_input
            allocate(out%index(0),out%estimate(0),out%lower(0),out%upper(0),out%threshold(0))
            return
        end select
        m = endv-startv+1
        allocate(out%index(m),out%estimate(m),out%lower(m),out%upper(m),out%threshold(m))
        z = normal_quantile(1.0_dp-(1.0_dp-civ)/2.0_dp)
        do i = 1,m
            k = endv-i+1
            out%index(i) = real(k,dp)
            out%estimate(i) = yall(k)
            out%threshold(i) = find_threshold(data,k)
            if (trim(opt) == 'quantile') then
                out%lower(i) = out%estimate(i)
                out%upper(i) = out%estimate(i)
            else
                se = out%estimate(i)/sqrt(real(k,dp))
                out%lower(i) = out%estimate(i)-z*se
                out%upper(i) = out%estimate(i)+z*se
            end if
        end do
    end function hill_estimates

    function extremal_index(data, block_size, start_k, end_k) result(out)
        real(dp), intent(in) :: data(:)
        integer, intent(in) :: block_size
        integer, intent(in), optional :: start_k, end_k
        type(matrix_result) :: out
        real(dp), allocatable :: sorted(:), bmax(:), unique_b(:)
        integer :: n, kblocks, startv, endv, nb, nu, i, j, kcount, ncount, m
        real(dp) :: r

        n = size(data)
        if (n < 2 .or. block_size <= 0) then
            out%status = evir_invalid_input
            allocate(out%values(0,0))
            return
        end if
        sorted = data
        call sort_descending(sorted)
        bmax = block_maxima(data,block_size)
        call sort_descending(bmax)
        nb = size(bmax)
        allocate(unique_b(nb))
        nu = 0
        do i = 1,nb
            if (i == 1) then
                nu = nu+1
                unique_b(nu) = bmax(i)
            else if (abs(bmax(i)-bmax(i-1)) > 10.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(bmax(i)))) then
                nu = nu+1
                unique_b(nu) = bmax(i)
            end if
        end do
        if (nu <= 1) then
            out%status = evir_domain_error
            allocate(out%values(0,0))
            return
        end if
        startv = 5
        if (present(start_k)) startv = start_k
        kblocks = nint(real(n,dp)/real(block_size,dp))
        endv = kblocks
        if (present(end_k)) endv = end_k
        r = real(block_size,dp)
        m = 0
        do i = 2,nu
            kcount = first_position(bmax,unique_b(i))-1
            if (kcount < endv .and. kcount >= startv) m = m+1
        end do
        allocate(out%values(m,5))
        j = 0
        do i = 2,nu
            kcount = first_position(bmax,unique_b(i))-1
            ncount = first_position(sorted,unique_b(i))-1
            if (kcount < endv .and. kcount >= startv .and. ncount > 0) then
                j = j+1
                out%values(j,1) = real(ncount,dp)
                out%values(j,2) = real(kcount,dp)
                out%values(j,3) = unique_b(i)
                out%values(j,4) = real(kcount,dp)/real(ncount,dp)
                out%values(j,5) = log(1.0_dp-real(kcount,dp)/real(kblocks,dp))/ &
                    (r*log(1.0_dp-real(ncount,dp)/real(n,dp)))
            end if
        end do
    end function extremal_index

    function shape_stability(data, models, start_n, end_n, ci) result(out)
        real(dp), intent(in) :: data(:)
        integer, intent(in), optional :: models,start_n,end_n
        real(dp), intent(in), optional :: ci
        type(band_result) :: out
        integer :: nm,s,e,i,nex
        real(dp) :: civ,z
        type(gpd_fit_result) :: fit
        nm=30; if(present(models)) nm=models
        s=15; if(present(start_n)) s=start_n
        e=min(500,size(data)); if(present(end_n)) e=min(end_n,size(data))
        civ=0.95_dp; if(present(ci)) civ=ci
        if(nm<1.or.s<2.or.e<s) then
            out%status=evir_invalid_input
            allocate(out%index(0),out%estimate(0),out%lower(0),out%upper(0),out%threshold(0))
            return
        end if
        allocate(out%index(nm),out%estimate(nm),out%lower(nm),out%upper(nm),out%threshold(nm))
        z=normal_quantile(1.0_dp-(1.0_dp-civ)/2.0_dp)
        do i=1,nm
            if(nm==1) then
                nex=e
            else
                nex=int(real(e,dp)+(real(s-e,dp)*real(i-1,dp)/real(nm-1,dp)))
            end if
            fit=fit_gpd(data,nextremes=nex,information='expected')
            out%index(i)=real(nex,dp)
            out%threshold(i)=fit%threshold
            out%estimate(i)=fit%xi
            out%lower(i)=fit%xi-z*fit%se(1)
            out%upper(i)=fit%xi+z*fit%se(1)
            if(fit%status/=evir_ok) out%status=fit%status
        end do
    end function shape_stability

    function quantile_stability(data,p,models,start_n,end_n,ci) result(out)
        real(dp), intent(in) :: data(:),p
        integer, intent(in), optional :: models,start_n,end_n
        real(dp), intent(in), optional :: ci
        type(band_result) :: out
        integer :: nm,s,e,i,nex,status
        real(dp) :: civ,se
        type(gpd_fit_result) :: fit
        nm=30; if(present(models)) nm=models
        s=15; if(present(start_n)) s=start_n
        e=min(500,size(data)); if(present(end_n)) e=min(end_n,size(data))
        civ=0.95_dp; if(present(ci)) civ=ci
        if(nm<1.or.s<2.or.e<s.or.p<=0.0_dp.or.p>=1.0_dp) then
            out%status=evir_invalid_input
            allocate(out%index(0),out%estimate(0),out%lower(0),out%upper(0),out%threshold(0))
            return
        end if
        allocate(out%index(nm),out%estimate(nm),out%lower(nm),out%upper(nm),out%threshold(nm))
        do i=1,nm
            if(nm==1) then
                nex=e
            else
                nex=int(real(e,dp)+(real(s-e,dp)*real(i-1,dp)/real(nm-1,dp)))
            end if
            fit=fit_gpd(data,nextremes=nex,information='expected')
            out%index(i)=real(nex,dp)
            out%threshold(i)=fit%threshold
            call gpd_quantile_wald(fit,p,out%lower(i),out%estimate(i),se,out%upper(i),civ,.true.,status)
            if(status/=evir_ok) out%status=status
        end do
    end function quantile_stability

    function tail_curve(fit,n_points,extend) result(out)
        type(gpd_fit_result), intent(in) :: fit
        integer, intent(in), optional :: n_points
        real(dp), intent(in), optional :: extend
        type(tail_curve_result) :: out
        real(dp), allocatable :: sorted(:)
        integer :: n,i,np
        real(dp) :: ext,plotmax,prob,p
        n=size(fit%exceedances)
        np=1000; if(present(n_points)) np=n_points
        ext=1.5_dp; if(present(extend)) ext=extend
        if(n==0.or.np<2.or.ext<=1.0_dp) then
            out%status=evir_invalid_input
            allocate(out%empirical_x(0),out%empirical_y(0),out%model_x(0),out%model_y(0))
            return
        end if
        sorted=fit%exceedances
        call sort_ascending(sorted)
        allocate(out%empirical_x(n),out%empirical_y(n),out%model_x(np),out%model_y(np))
        out%empirical_x=sorted
        prob=1.0_dp-fit%p_less_threshold
        do i=1,n
            out%empirical_y(i)=prob*(1.0_dp-ppoints(i,n))
        end do
        plotmax=maxval(sorted)*ext
        do i=1,np
            p=real(i-1,dp)/real(np-1,dp)
            out%model_x(i)=min(max(qgpd(p,fit%xi,fit%threshold,fit%beta),fit%threshold),plotmax)
            out%model_y(i)=prob*(1.0_dp-pgpd(out%model_x(i),fit%xi,fit%threshold,fit%beta))
        end do
        out%threshold=fit%threshold
        out%shape=fit%xi
        out%scale=fit%beta*prob**fit%xi
        if(abs(fit%xi)<=1.0e-10_dp) then
            out%location=fit%threshold+out%scale*log(prob)
        else
            out%location=fit%threshold-out%scale*(prob**(-fit%xi)-1.0_dp)/fit%xi
        end if
    end function tail_curve

    integer function first_position(x,value) result(pos)
        real(dp),intent(in)::x(:),value
        integer::i
        pos=size(x)+1
        do i=1,size(x)
            if(abs(x(i)-value)<=10.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(value))) then
                pos=i; return
            end if
        end do
    end function first_position

    pure character(len=12) function lowercase(text) result(out)
        character(len=*),intent(in)::text
        integer::i,c,n
        out=' '
        n=min(len_trim(text),len(out))
        do i=1,n
            c=iachar(text(i:i))
            if(c>=iachar('A').and.c<=iachar('Z')) c=c+32
            out(i:i)=achar(c)
        end do
    end function lowercase

end module evir_eda
