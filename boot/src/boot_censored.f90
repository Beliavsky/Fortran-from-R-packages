module boot_censored
    use boot_kinds, only : dp
    use boot_resampling, only : ordinary_array
    implicit none
    private
    public :: censored_case_indices, sample_product_limit, sample_conditional_censoring
contains
    subroutine censored_case_indices(n,r,strata,idx)
        integer,intent(in)::n,r,strata(n)
        integer,intent(out)::idx(r,n)
        call ordinary_array(n,r,strata,idx)
    end subroutine censored_case_indices

    subroutine sample_product_limit(time,survival,r,n,draws)
        real(dp),intent(in)::time(:),survival(:)
        integer,intent(in)::r,n
        real(dp),intent(out)::draws(r,n)
        real(dp),allocatable::tt(:),prob(:)
        real(dp)::prev,u,s
        integer::m,i,j,k
        m=size(time)
        if(size(survival)/=m)error stop "sample_product_limit: mismatch"
        if(survival(m)>0.0_dp)then
            allocate(tt(m+1),prob(m+1))
            tt(1:m)=time
            tt(m+1)=huge(1.0_dp)
        else
            allocate(tt(m),prob(m))
            tt=time
        end if
        prev=1.0_dp
        do i=1,m
        prob(i)=prev-survival(i)
        prev=survival(i)
        end do
        if(size(prob)>m)prob(m+1)=prev
        s=sum(prob)
        if(s<=0.0_dp)error stop "sample_product_limit: invalid survival curve"
        prob=prob/s
        do i=1,r
        do j=1,n
            call random_number(u)
            s=0.0_dp
            k=size(prob)
            do m=1,size(prob)
            s=s+prob(m)
            if(u<=s)then
            k=m
            exit
            end if
            end do
            draws(i,j)=tt(k)
        end do
        end do
    end subroutine sample_product_limit

    subroutine sample_conditional_censoring(obs_time,event,time,survival,r,draws)
        real(dp),intent(in)::obs_time(:),time(:),survival(:)
        integer,intent(in)::event(:),r
        real(dp),intent(out)::draws(r,size(obs_time))
        real(dp),allocatable::tt(:),prob(:)
        real(dp)::prev,s,u
        integer::m,i,j,k
        m=size(time)
        if(size(survival)/=m .or. size(event)/=size(obs_time))error stop "sample_conditional_censoring: mismatch"
        if(survival(m)>0.0_dp)then
        allocate(tt(m+1),prob(m+1))
        tt(1:m)=time
        tt(m+1)=huge(1.0_dp)
        else
        allocate(tt(m),prob(m))
        tt=time
        end if
        prev=1.0_dp
        do i=1,m
        prob(i)=prev-survival(i)
        prev=survival(i)
        end do
        if(size(prob)>m)prob(m+1)=prev
        do j=1,size(obs_time)
            if(event(j)==0)then
            draws(:,j)=obs_time(j)
            else
                s=sum(prob,mask=tt>obs_time(j))
                if(s<=0.0_dp)then
                draws(:,j)=huge(1.0_dp)
                cycle
                end if
                do i=1,r
                    call random_number(u)
                    u=u*s
                    s=0.0_dp
                    k=size(tt)
                    do m=1,size(tt)
                        if(tt(m)>obs_time(j))then
                        s=s+prob(m)
                        if(u<=s)then
                        k=m
                        exit
                        end if
                        end if
                    end do
                    draws(i,j)=tt(k)
                end do
            end if
        end do
    end subroutine sample_conditional_censoring
end module boot_censored
