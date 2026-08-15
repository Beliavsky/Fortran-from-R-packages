module actuar_aggregate_v02
    use actuar_kinds, only: dp
    use actuar_special, only: normal_cdf, normal_quantile
    use actuar_risk, only: aggregate_dist_t, convolve_pmf
    implicit none
    private
    abstract interface
        integer function frequency_rng_iface()
        end function frequency_rng_iface
        function severity_rng_iface() result(x)
            import dp
            real(dp) :: x
        end function severity_rng_iface
    end interface
    public :: aggregate_exact, aggregate_empirical_grid
    public :: aggregate_normal_cdf, aggregate_normal_quantile
    public :: aggregate_npower_cdf, aggregate_npower_quantile
    public :: compound_sums, rcompound_callbacks, mixture_sample
contains
    function aggregate_exact(severity_pmf, frequency_pmf, x_scale) result(dist)
        real(dp), intent(in) :: severity_pmf(:), frequency_pmf(:)
        real(dp), intent(in), optional :: x_scale
        type(aggregate_dist_t) :: dist
        real(dp), allocatable :: fs(:), fxc(:), next(:)
        integer :: nmax, i, r, j
        if (size(frequency_pmf) < 1 .or. size(severity_pmf) < 1) then
            allocate(dist%pmf(0),dist%cdf(0))
            return
        end if
        nmax = size(frequency_pmf)-1
        r = nmax*(size(severity_pmf)-1)+1
        allocate(fs(r)); fs=0.0_dp
        fs(1)=frequency_pmf(1)
        allocate(fxc(1)); fxc=1.0_dp
        do i=1,nmax
            next=convolve_pmf(severity_pmf,fxc)
            call move_alloc(next,fxc)
            fs(1:size(fxc))=fs(1:size(fxc))+frequency_pmf(i+1)*fxc
        end do
        allocate(dist%pmf(r),dist%cdf(r));dist%pmf=fs
        if (present(x_scale)) dist%x_scale=x_scale
        if(r>0) then
            dist%cdf(1)=dist%pmf(1)
            do j=2,r
                dist%cdf(j)=min(1.0_dp,dist%cdf(j-1)+dist%pmf(j))
            end do
        end if
    end function aggregate_exact

    function aggregate_empirical_grid(samples, step, max_value) result(dist)
        real(dp), intent(in) :: samples(:), step
        real(dp), intent(in), optional :: max_value
        type(aggregate_dist_t) :: dist
        real(dp) :: xmax
        integer :: nbin, i, k
        if (size(samples)==0 .or. step<=0.0_dp) then
            allocate(dist%pmf(0),dist%cdf(0));return
        end if
        xmax=maxval(samples);if(present(max_value)) xmax=max_value
        nbin=max(1,nint(xmax/step)+1)
        allocate(dist%pmf(nbin),dist%cdf(nbin));dist%pmf=0.0_dp
        do i=1,size(samples)
            k=nint(max(0.0_dp,samples(i))/step)+1
            if(k<=nbin) dist%pmf(k)=dist%pmf(k)+1.0_dp/real(size(samples),dp)
        end do
        dist%cdf(1)=dist%pmf(1)
        do i=2,nbin;dist%cdf(i)=min(1.0_dp,dist%cdf(i-1)+dist%pmf(i));end do
        dist%x_scale=step
    end function aggregate_empirical_grid

    pure real(dp) function aggregate_normal_cdf(x, mean, variance) result(p)
        real(dp), intent(in) :: x, mean, variance
        if (variance<=0.0_dp) then
            p=merge(0.0_dp,1.0_dp,x<mean)
        else
            p=normal_cdf((x-mean)/sqrt(variance))
        end if
    end function aggregate_normal_cdf

    pure real(dp) function aggregate_normal_quantile(p, mean, variance) result(x)
        real(dp), intent(in) :: p, mean, variance
        if(variance<=0.0_dp) then
            x=mean
        else
            x=mean+sqrt(variance)*normal_quantile(p)
        end if
    end function aggregate_normal_quantile

    pure real(dp) function aggregate_npower_cdf(x, mean, variance, skewness) result(p)
        real(dp), intent(in) :: x, mean, variance, skewness
        real(dp) :: zarg
        if (x<=mean .or. variance<=0.0_dp .or. skewness==0.0_dp) then
            p=0.0_dp
            return
        end if
        zarg=1.0_dp+9.0_dp/skewness**2+6.0_dp*(x-mean)/(sqrt(variance)*skewness)
        if(zarg<=0.0_dp) then
            p=0.0_dp
        else
            p=normal_cdf(sqrt(zarg)-3.0_dp/skewness)
        end if
    end function aggregate_npower_cdf

    pure real(dp) function aggregate_npower_quantile(p, mean, variance, skewness) result(x)
        real(dp), intent(in) :: p, mean, variance, skewness
        real(dp) :: z
        if(p<=0.0_dp) then
            x=mean
        else if(p>=1.0_dp) then
            x=huge(1.0_dp)
        else if(skewness==0.0_dp) then
            x=aggregate_normal_quantile(p,mean,variance)
        else
            z=normal_quantile(p)
            x=mean+sqrt(variance)*skewness*((z+3.0_dp/skewness)**2-1.0_dp-9.0_dp/skewness**2)/6.0_dp
            x=max(mean,x)
        end if
    end function aggregate_npower_quantile

    pure function compound_sums(counts, severity_draws) result(sums)
        integer, intent(in) :: counts(:)
        real(dp), intent(in) :: severity_draws(:,:)
        real(dp), allocatable :: sums(:)
        integer :: i,n
        allocate(sums(size(counts)));sums=0.0_dp
        do i=1,size(counts)
            n=max(0,min(counts(i),size(severity_draws,2)))
            if(n>0) sums(i)=sum(severity_draws(i,1:n))
        end do
    end function compound_sums


    function rcompound_callbacks(nsim,frequency_rng,severity_rng) result(sums)
        integer,intent(in)::nsim
        procedure(frequency_rng_iface)::frequency_rng
        procedure(severity_rng_iface)::severity_rng
        real(dp),allocatable::sums(:)
        integer::i,j,n
        allocate(sums(max(0,nsim)));sums=0.0_dp
        do i=1,nsim
            n=max(0,frequency_rng())
            do j=1,n
                sums(i)=sums(i)+severity_rng()
            end do
        end do
    end function rcompound_callbacks

    function mixture_sample(weights, component_samples) result(x)
        real(dp), intent(in) :: weights(:), component_samples(:,:)
        real(dp), allocatable :: x(:)
        real(dp) :: sw,u,c
        integer :: i,j,k
        allocate(x(size(component_samples,1)));x=0.0_dp
        sw=sum(max(weights,0.0_dp))
        if(sw<=0.0_dp .or. size(weights)/=size(component_samples,2)) return
        do i=1,size(x)
            call random_number(u);u=u*sw;c=0.0_dp;k=size(weights)
            do j=1,size(weights)
                c=c+max(0.0_dp,weights(j))
                if(u<=c) then;k=j;exit;end if
            end do
            x(i)=component_samples(i,k)
        end do
    end function mixture_sample
end module actuar_aggregate_v02
