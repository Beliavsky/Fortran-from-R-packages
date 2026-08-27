! Spectral matrix-function methods corresponding to expm's R_Eigen/Eigen paths.
! GPL-3.0-or-later; see LICENSE and LICENSES.md.
module expm_eigen_methods
    use expm_kinds, only : dp
    use expm_linalg, only : eye_complex, solve_complex, normf_complex
    use expm_matrix_functions, only : expm_ward77
    use r_linalg, only : general_complex_eigen
    implicit none
    private
    public :: expm_eigen, logm_eigen, expm_hybrid_eigen_ward
contains
    function expm_eigen(a,info) result(x)
        real(dp), intent(in) :: a(:,:)
        integer, intent(out), optional :: info
        complex(dp), allocatable :: x(:,:)
        integer :: istat
        call spectral_function(a,1,x,istat)
        if(present(info)) info=istat
        if(istat/=0 .and. .not.present(info)) error stop "expm_eigen: eigendecomposition failed"
    end function expm_eigen

    function logm_eigen(a,info) result(x)
        real(dp), intent(in) :: a(:,:)
        integer, intent(out), optional :: info
        complex(dp), allocatable :: x(:,:)
        integer :: istat
        call spectral_function(a,2,x,istat)
        if(present(info)) info=istat
        if(istat/=0 .and. .not.present(info)) error stop "logm_eigen: eigendecomposition failed"
    end function logm_eigen

    function expm_hybrid_eigen_ward(a,tol) result(x)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(in), optional :: tol
        real(dp), allocatable :: x(:,:)
        complex(dp), allocatable :: z(:,:)
        real(dp) :: threshold
        integer :: info
        threshold=epsilon(1.0_dp); if(present(tol)) threshold=tol
        z=expm_eigen(a,info)
        if(info==0 .and. maxval(abs(aimag(z)))<=max(threshold,1.0e-13_dp)*(1.0_dp+maxval(abs(real(z,dp))))) then
            x=real(z,dp)
        else
            x=expm_ward77(a)
        end if
    end function expm_hybrid_eigen_ward

    subroutine spectral_function(a,which,x,info)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in) :: which
        complex(dp), allocatable, intent(out) :: x(:,:)
        integer, intent(out) :: info
        complex(dp), allocatable :: aa(:,:),w(:),vr(:,:),vinv(:,:),fv(:,:)
        real(dp) :: resid
        integer :: n,i
        n=size(a,1); if(size(a,2)/=n) error stop "spectral_function: matrix must be square"
        allocate(aa(n,n)); aa=cmplx(a,0.0_dp,dp)
        call general_complex_eigen(aa,w,vr,info)
        if(info/=0) then; allocate(x(n,n)); x=cmplx(0.0_dp,0.0_dp,dp); return; end if
        vinv=solve_complex(vr,eye_complex(n),info)
        if(info/=0) then; allocate(x(n,n)); x=cmplx(0.0_dp,0.0_dp,dp); return; end if
        fv=vr
        select case(which)
        case(1)
            do i=1,n; fv(:,i)=fv(:,i)*exp(w(i)); end do
        case(2)
            do i=1,n; fv(:,i)=fv(:,i)*log(w(i)); end do
        end select
        x=matmul(fv,vinv)
        resid=normf_complex(matmul(vr,matmul(diagonal_matrix(w),vinv))-cmplx(a,0.0_dp,dp))/(1.0_dp+sqrt(sum(a*a)))
        if(resid>1.0e-7_dp) info=1
    end subroutine spectral_function

    function diagonal_matrix(d) result(a)
        complex(dp), intent(in) :: d(:)
        complex(dp), allocatable :: a(:,:)
        integer :: i,n
        n=size(d); allocate(a(n,n)); a=cmplx(0.0_dp,0.0_dp,dp)
        do i=1,n; a(i,i)=d(i); end do
    end function diagonal_matrix
end module expm_eigen_methods
