! SPDX-License-Identifier: GPL-2.0-only
module nlsr_stats
    use nlsr_kinds, only : dp
    use nlsr_types, only : nlsr_result
    use nlsr_linalg, only : invert_symmetric_positive
    implicit none
    private
    public :: nlsr_standard_errors
contains
    subroutine nlsr_standard_errors(result, covariance, se, sigma, ok)
        type(nlsr_result), intent(in) :: result
        real(dp), intent(out) :: covariance(:,:), se(:), sigma
        logical, intent(out) :: ok
        integer :: p, n, nf, i, j, ii, jj
        integer, allocatable :: idx(:)
        real(dp), allocatable :: jw(:,:), xtx(:,:), inv(:,:)
        logical :: invok
        p=size(result%coefficients); n=size(result%residuals)
        covariance=0.0_dp; se=0.0_dp; sigma=huge(1.0_dp); ok=.false.
        if (size(covariance,1)/=p .or. size(covariance,2)/=p .or. size(se)/=p) return
        nf=count(.not.result%masked)
        if (nf<1 .or. n<=nf) return
        allocate(idx(nf)); ii=0
        do i=1,p
            if (.not.result%masked(i)) then; ii=ii+1; idx(ii)=i; end if
        end do
        allocate(jw(n,nf),xtx(nf,nf),inv(nf,nf))
        do j=1,nf; jw(:,j)=result%jacobian(:,idx(j)); end do
        xtx=matmul(transpose(jw),jw)
        call invert_symmetric_positive(xtx,inv,invok)
        if (.not.invok) return
        sigma=sqrt(result%ssquares/real(n-nf,dp))
        do ii=1,nf
            i=idx(ii)
            do jj=1,nf
                j=idx(jj)
                covariance(i,j)=sigma*sigma*inv(ii,jj)
            end do
            se(i)=sqrt(max(0.0_dp,covariance(i,i)))
        end do
        ok=.true.
    end subroutine nlsr_standard_errors
end module nlsr_stats
