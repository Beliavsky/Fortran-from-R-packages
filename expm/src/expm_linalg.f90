! Linear algebra helpers for expm-fortran.
! GPL-3.0-or-later; see LICENSE and LICENSES.md.
module expm_linalg
    use expm_kinds, only : dp
    use r_linalg, only : balance_matrix
    use r_linalg, only : shared_singular_values => singular_values
    use r_linalg, only : shared_solve_system => solve_system
    implicit none
    private

    type, public :: balance_real_result
        real(dp), allocatable :: z(:,:)
        real(dp), allocatable :: scale(:)
        integer :: ilo = 1
        integer :: ihi = 0
        integer :: info = 0
    end type balance_real_result

    type, public :: balance_complex_result
        complex(dp), allocatable :: z(:,:)
        real(dp), allocatable :: scale(:)
        integer :: ilo = 1
        integer :: ihi = 0
        integer :: info = 0
    end type balance_complex_result

    public :: eye_real, eye_complex, norm1_real, norm1_complex
    public :: norminf_real, normf_real, normf_complex
    public :: solve_real, solve_complex, inverse_complex
    public :: balance_real, balance_complex, reverse_balance_real, reverse_balance_complex
    public :: spectral_norm_real

contains
    function eye_real(n) result(a)
        integer, intent(in) :: n
        real(dp), allocatable :: a(:,:)
        integer :: i
        allocate(a(n,n)); a = 0.0_dp
        do i=1,n; a(i,i)=1.0_dp; end do
    end function eye_real

    function eye_complex(n) result(a)
        integer, intent(in) :: n
        complex(dp), allocatable :: a(:,:)
        integer :: i
        allocate(a(n,n)); a = cmplx(0.0_dp,0.0_dp,dp)
        do i=1,n; a(i,i)=cmplx(1.0_dp,0.0_dp,dp); end do
    end function eye_complex

    pure real(dp) function norm1_real(a)
        real(dp), intent(in) :: a(:,:)
        integer :: j
        norm1_real = 0.0_dp
        do j=1,size(a,2); norm1_real=max(norm1_real,sum(abs(a(:,j)))); end do
    end function norm1_real

    pure real(dp) function norm1_complex(a)
        complex(dp), intent(in) :: a(:,:)
        integer :: j
        norm1_complex = 0.0_dp
        do j=1,size(a,2); norm1_complex=max(norm1_complex,sum(abs(a(:,j)))); end do
    end function norm1_complex

    pure real(dp) function norminf_real(a)
        real(dp), intent(in) :: a(:,:)
        integer :: i
        norminf_real=0.0_dp
        do i=1,size(a,1); norminf_real=max(norminf_real,sum(abs(a(i,:)))); end do
    end function norminf_real

    pure real(dp) function normf_real(a)
        real(dp), intent(in) :: a(:,:)
        normf_real = sqrt(sum(a*a))
    end function normf_real

    pure real(dp) function normf_complex(a)
        complex(dp), intent(in) :: a(:,:)
        normf_complex = sqrt(sum(abs(a)**2))
    end function normf_complex

    function solve_real(a,b,info) result(x)
        real(dp), intent(in) :: a(:,:), b(:,:)
        integer, intent(out), optional :: info
        real(dp), allocatable :: x(:,:)
        integer :: n,nrhs,istat
        n=size(a,1); nrhs=size(b,2)
        if(size(a,2)/=n .or. size(b,1)/=n) error stop "solve_real: nonconformable arrays"
        allocate(x(n,nrhs))
        call shared_solve_system(a,b,x,istat)
        if(present(info)) info=istat
        if(istat/=0 .and. .not.present(info)) error stop "solve_real: dgesv failed"
    end function solve_real

    function solve_complex(a,b,info) result(x)
        complex(dp), intent(in) :: a(:,:), b(:,:)
        integer, intent(out), optional :: info
        complex(dp), allocatable :: x(:,:)
        integer :: n,nrhs,istat
        n=size(a,1); nrhs=size(b,2)
        if(size(a,2)/=n .or. size(b,1)/=n) error stop "solve_complex: nonconformable arrays"
        allocate(x(n,nrhs))
        call shared_solve_system(a,b,x,istat)
        if(present(info)) info=istat
        if(istat/=0 .and. .not.present(info)) error stop "solve_complex: zgesv failed"
    end function solve_complex

    function inverse_complex(a,info) result(ai)
        complex(dp), intent(in) :: a(:,:)
        integer, intent(out), optional :: info
        complex(dp), allocatable :: ai(:,:)
        ai = solve_complex(a,eye_complex(size(a,1)),info)
    end function inverse_complex

    function balance_real(a,job) result(r)
        real(dp), intent(in) :: a(:,:)
        character(len=1), intent(in), optional :: job
        type(balance_real_result) :: r
        character(len=1) :: jb
        integer :: n
        n=size(a,1); if(size(a,2)/=n) error stop "balance_real: matrix must be square"
        jb='B'; if(present(job)) jb=job
        call balance_matrix(a,r%z,r%scale,r%ilo,r%ihi,r%info,job=jb)
    end function balance_real

    function balance_complex(a,job) result(r)
        complex(dp), intent(in) :: a(:,:)
        character(len=1), intent(in), optional :: job
        type(balance_complex_result) :: r
        character(len=1) :: jb
        integer :: n
        n=size(a,1); if(size(a,2)/=n) error stop "balance_complex: matrix must be square"
        jb='B'; if(present(job)) jb=job
        call balance_matrix(a,r%z,r%scale,r%ilo,r%ihi,r%info,job=jb)
    end function balance_complex

    subroutine reverse_balance_real(x,perm,scal)
        real(dp), intent(inout) :: x(:,:)
        type(balance_real_result), intent(in) :: perm,scal
        integer :: n,i,p
        real(dp), allocatable :: tmp(:)
        n=size(x,1); allocate(tmp(n))
        do i=1,n
            x(i,:) = x(i,:) * scal%scale(i)
            x(:,i) = x(:,i) / scal%scale(i)
        end do
        if(perm%ilo>1) then
            do i=perm%ilo-1,1,-1
                p=nint(perm%scale(i)); if(p/=i) then
                    tmp=x(:,i); x(:,i)=x(:,p); x(:,p)=tmp
                    tmp=x(i,:); x(i,:)=x(p,:); x(p,:)=tmp
                end if
            end do
        end if
        if(perm%ihi<n) then
            do i=perm%ihi+1,n
                p=nint(perm%scale(i)); if(p/=i) then
                    tmp=x(:,i); x(:,i)=x(:,p); x(:,p)=tmp
                    tmp=x(i,:); x(i,:)=x(p,:); x(p,:)=tmp
                end if
            end do
        end if
    end subroutine reverse_balance_real

    subroutine reverse_balance_complex(x,perm,scal)
        complex(dp), intent(inout) :: x(:,:)
        type(balance_complex_result), intent(in) :: perm,scal
        integer :: n,i,p
        complex(dp), allocatable :: tmp(:)
        n=size(x,1); allocate(tmp(n))
        do i=1,n
            x(i,:) = x(i,:) * scal%scale(i)
            x(:,i) = x(:,i) / scal%scale(i)
        end do
        if(perm%ilo>1) then
            do i=perm%ilo-1,1,-1
                p=nint(perm%scale(i)); if(p/=i) then
                    tmp=x(:,i); x(:,i)=x(:,p); x(:,p)=tmp
                    tmp=x(i,:); x(i,:)=x(p,:); x(p,:)=tmp
                end if
            end do
        end if
        if(perm%ihi<n) then
            do i=perm%ihi+1,n
                p=nint(perm%scale(i)); if(p/=i) then
                    tmp=x(:,i); x(:,i)=x(:,p); x(:,p)=tmp
                    tmp=x(i,:); x(i,:)=x(p,:); x(p,:)=tmp
                end if
            end do
        end if
    end subroutine reverse_balance_complex

    real(dp) function spectral_norm_real(a) result(sn)
        real(dp), intent(in) :: a(:,:)
        real(dp), allocatable :: s(:)
        integer :: info
        call shared_singular_values(a,s,info)
        if(info/=0) error stop "spectral_norm_real: dgesvd failed"
        sn=maxval(s)
    end function spectral_norm_real
end module expm_linalg
