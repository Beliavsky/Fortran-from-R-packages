! Matrix square root and logarithm support for expm-fortran.
! sqrtm uses a complex Schur decomposition and triangular recurrence, mirroring
! the Schur-based algorithm used by the R package while avoiding R/Matrix glue.
! logm uses inverse scaling-and-squaring with repeated Schur square roots and
! the atanh series for log(I+X).
! GPL-3.0-or-later; see LICENSE and LICENSES.md.
module expm_logsqrt
    use expm_kinds, only : dp
    use expm_linalg, only : eye_complex, norm1_complex, normf_complex, solve_complex
    implicit none
    private
    public :: sqrtm, logm
    interface sqrtm
        module procedure sqrtm_real
        module procedure sqrtm_complex
    end interface sqrtm
    interface logm
        module procedure logm_real
        module procedure logm_complex
    end interface logm

    abstract interface
        logical function zselect_iface(z)
            import dp
            complex(dp), intent(in) :: z
        end function zselect_iface
    end interface

    interface
        subroutine zgees(jobvs,sort,select,n,a,lda,sdim,w,vs,ldvs,work,lwork,rwork,bwork,info)
            import dp,zselect_iface
            character(len=1), intent(in) :: jobvs,sort
            procedure(zselect_iface) :: select
            integer, intent(in) :: n,lda,ldvs,lwork
            complex(dp), intent(inout) :: a(lda,*)
            integer, intent(out) :: sdim,info
            complex(dp), intent(out) :: w(*),vs(ldvs,*),work(*)
            real(dp), intent(out) :: rwork(*)
            logical, intent(out) :: bwork(*)
        end subroutine zgees
    end interface
contains
    function sqrtm_real(a) result(x)
        real(dp), intent(in) :: a(:,:)
        complex(dp), allocatable :: x(:,:),ac(:,:)
        ac=cmplx(a,0.0_dp,dp); x=sqrtm_complex(ac)
    end function sqrtm_real

    function sqrtm_complex(a) result(x)
        complex(dp), intent(in) :: a(:,:)
        complex(dp), allocatable :: x(:,:),t(:,:),q(:,:),w(:),work(:),r(:,:)
        complex(dp) :: wk(1),denom,sm
        real(dp), allocatable :: rwork(:)
        logical, allocatable :: bwork(:)
        integer :: n,lwork,info,sdim,i,j,k,p
        n=size(a,1); if(size(a,2)/=n) error stop "sqrtm: matrix must be square"
        if(n==0) then; allocate(x(0,0)); return; end if
        allocate(t(n,n),q(n,n),w(n),rwork(max(1,n)),bwork(max(1,n))); t=a
        lwork=-1
        call zgees('V','N',select_none,n,t,n,sdim,w,q,n,wk,lwork,rwork,bwork,info)
        if(info/=0) error stop "sqrtm: zgees workspace query failed"
        lwork=max(2*n,int(real(wk(1),dp))); allocate(work(lwork)); t=a
        call zgees('V','N',select_none,n,t,n,sdim,w,q,n,work,lwork,rwork,bwork,info)
        if(info/=0) error stop "sqrtm: Schur decomposition failed"
        allocate(r(n,n)); r=cmplx(0.0_dp,0.0_dp,dp)
        do i=1,n; r(i,i)=sqrt(t(i,i)); end do
        do p=1,n-1
            do i=1,n-p
                j=i+p; sm=cmplx(0.0_dp,0.0_dp,dp)
                do k=i+1,j-1; sm=sm+r(i,k)*r(k,j); end do
                denom=r(i,i)+r(j,j)
                if(abs(denom)<=100.0_dp*epsilon(1.0_dp)) error stop "sqrtm: singular Sylvester denominator"
                r(i,j)=(t(i,j)-sm)/denom
            end do
        end do
        x=matmul(q,matmul(r,transpose(conjg(q))))
    end function sqrtm_complex

    function logm_real(a) result(l)
        real(dp), intent(in) :: a(:,:)
        complex(dp), allocatable :: l(:,:),ac(:,:)
        ac=cmplx(a,0.0_dp,dp); l=logm_complex(ac)
    end function logm_real

    function logm_complex(a) result(l)
        complex(dp), intent(in) :: a(:,:)
        complex(dp), allocatable :: l(:,:),x(:,:),id(:,:),z(:,:),z2(:,:),term(:,:),den(:,:),rhs(:,:)
        real(dp) :: nz,nt
        integer :: n,s,k,info
        n=size(a,1); if(size(a,2)/=n) error stop "logm: matrix must be square"
        id=eye_complex(n); x=a; s=0
        do
            den=x+id; rhs=x-id
            z=transpose(solve_complex(transpose(den),transpose(rhs),info))
            if(info/=0) error stop "logm: A + I became singular"
            nz=norm1_complex(z)
            if(nz<0.45_dp .or. s>=30) exit
            x=sqrtm_complex(x); s=s+1
        end do
        if(nz>=1.0_dp) error stop "logm: inverse scaling failed to reach convergence region"
        z2=matmul(z,z); term=z; l=term
        do k=1,800
            term=matmul(term,z2); l=l+term/real(2*k+1,dp)
            nt=norm1_complex(term)/real(2*k+1,dp)
            if(nt<epsilon(1.0_dp)*max(1.0_dp,norm1_complex(l))) exit
        end do
        l=(2.0_dp**(s+1))*l
    end function logm_complex

    logical function select_none(z)
        complex(dp), intent(in) :: z
        select_none=.false.
        if(real(z,dp)>huge(1.0_dp)) select_none=.true.
    end function select_none
end module expm_logsqrt
