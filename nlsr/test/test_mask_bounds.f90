program test_mask_bounds
    use nlsr
    implicit none
    real(dp),parameter :: y(12)=[5.308_dp,7.24_dp,9.638_dp,12.866_dp,17.069_dp,23.192_dp, &
        31.443_dp,38.558_dp,50.156_dp,62.948_dp,75.995_dp,91.972_dp]
    real(dp)::start(3),lo(3),up(3)
    type(nlsr_result)::fit
    type(nlsr_control)::ctl
    start=[100.0_dp,10.0_dp,0.25_dp]; lo=[0.0_dp,0.0_dp,0.25_dp]; up=[500.0_dp,100.0_dp,0.25_dp]
    ctl=nlsr_control(); ctl%stepredn=0.5_dp; ctl%nbtlim=12
    call nlfb(start,12,res,fit,ctl,jac,lo,up)
    if (abs(fit%coefficients(3)-0.25_dp)>1.0e-15_dp) error stop 'masked parameter moved'
    if (.not.fit%masked(3)) error stop 'mask not reported'
    if (fit%ssquares>50.0_dp) error stop 'bounded fit poor'
    print *, 'PASS test_mask_bounds',fit%coefficients,fit%ssquares
contains
    subroutine res(p,r,ierr)
        real(dp),intent(in)::p(:); real(dp),intent(out)::r(:); integer,intent(out)::ierr
        integer::i; real(dp)::e,d
        do i=1,12; e=exp(-p(3)*real(i,dp)); d=1+p(2)*e; r(i)=y(i)-p(1)/d; end do; ierr=0
    end subroutine
    subroutine jac(p,j,ierr)
        real(dp),intent(in)::p(:); real(dp),intent(out)::j(:,:); integer,intent(out)::ierr
        integer::i; real(dp)::e,d
        do i=1,12
            e=exp(-p(3)*real(i,dp)); d=1+p(2)*e
            j(i,1)=-1/d; j(i,2)=p(1)*e/(d*d); j(i,3)=-p(1)*p(2)*real(i,dp)*e/(d*d)
        end do; ierr=0
    end subroutine
end program
