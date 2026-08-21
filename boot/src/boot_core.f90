module boot_core
    use boot_kinds, only : dp
    use boot_resampling, only : ordinary_array, balanced_array, antithetic_array, importance_array, &
                               balanced_importance_array, frequency_array
    implicit none
    private
    public :: bootstrap_result, bootstrap_weighted
    type :: bootstrap_result
        real(dp) :: t0 = 0.0_dp
        real(dp), allocatable :: t(:)
        integer, allocatable :: indices(:,:), frequencies(:,:)
    end type bootstrap_result
    abstract interface
        function weighted_statistic(data,weights) result(value)
            import dp
            real(dp),intent(in)::data(:,:),weights(:)
            real(dp)::value
        end function weighted_statistic
    end interface
contains
    subroutine bootstrap_weighted(data,statistic,r,sim,result,strata,influence,weights)
        real(dp),intent(in)::data(:,:)
        procedure(weighted_statistic)::statistic
        integer,intent(in)::r
        character(len=*),intent(in)::sim
        type(bootstrap_result),intent(out)::result
        integer,intent(in),optional::strata(:)
        real(dp),intent(in),optional::influence(:),weights(:)
        integer::n,i,j
        integer,allocatable::st(:)
        real(dp),allocatable::w(:)
        n=size(data,1)
        allocate(st(n),w(n),result%indices(r,n),result%frequencies(r,n),result%t(r))
        if(present(strata))then
        st=strata
        else
        st=1
        end if
        do i=1,n
        w(i)=1.0_dp/real(count(st==st(i)),dp)
        end do
        result%t0=statistic(data,w)
        select case(trim(sim))
        case('ordinary')
            if(present(weights))then
            call importance_array(n,r,weights,st,result%indices)
            else
            call ordinary_array(n,r,st,result%indices)
            end if
        case('balanced')
            if(present(weights))then
            call balanced_importance_array(n,r,weights,st,result%indices)
            else
            call balanced_array(n,r,st,result%indices)
            end if
        case('antithetic')
            if(.not.present(influence))error stop "bootstrap_weighted: influence required for antithetic"
            call antithetic_array(n,r,influence,st,result%indices)
        case default
            error stop "bootstrap_weighted: unsupported sim"
        end select
        call frequency_array(result%indices,result%frequencies)
        do i=1,r
            do j=1,n
            w(j)=real(result%frequencies(i,j),dp)/real(count(st==st(j)),dp)
            end do
            result%t(i)=statistic(data,w)
        end do
    end subroutine bootstrap_weighted
end module boot_core
