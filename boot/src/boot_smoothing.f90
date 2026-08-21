module boot_smoothing
    use boot_kinds, only : dp, pi
    use boot_statistics, only : variance_dp
    implicit none
    private
    public :: smooth_frequencies
contains
    subroutine smooth_frequencies(theta,t,freq,strata,width,out)
        real(dp),intent(in)::theta(:),t(:),width
        integer,intent(in)::freq(:,:),strata(:)
        real(dp),intent(out)::out(size(theta),size(freq,2))
        real(dp)::eps,wgt,s
        real(dp),allocatable::kernel(:),f(:)
        integer::m,r,n,i,j,k
        r=size(freq,1)
        n=size(freq,2)
        m=size(theta)
        if(size(t)/=r .or. size(strata)/=n)error stop "smooth_frequencies: size mismatch"
        eps=width*sqrt(variance_dp(t))
        if(eps<=0.0_dp)error stop "smooth_frequencies: zero bandwidth"
        allocate(kernel(r),f(n))
        out=0.0_dp
        do k=1,m
            do i=1,r
                wgt=(theta(k)-t(i))/eps
                kernel(i)=exp(-0.5_dp*wgt*wgt)/(sqrt(2.0_dp*pi)*eps)
            end do
            do j=1,n
            f(j)=sum(real(freq(:,j),dp)*kernel)
            end do
            do j=1,n
                s=sum(f,mask=strata==strata(j))
                out(k,j)=f(j)/s
            end do
        end do
    end subroutine smooth_frequencies
end module boot_smoothing
