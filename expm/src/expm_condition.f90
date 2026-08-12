! Exponential condition number computations translated from expm.
! GPL-3.0-or-later; see LICENSE and LICENSES.md.
module expm_condition
    use expm_kinds, only : dp
    use expm_linalg, only : norm1_real, normf_real, spectral_norm_real
    use expm_frechet, only : expm_frechet_sps
    implicit none
    private
    public :: expm_cond_exact, expm_cond_1_est, expm_cond_f_est
contains
    subroutine expm_cond_exact(a,cond1,condf,expa)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: cond1,condf
        real(dp), allocatable, intent(out), optional :: expa(:,:)
        real(dp), allocatable :: k(:,:),e(:,:),x(:,:),l(:,:)
        real(dp) :: nk1,nkf
        integer :: n,i,j,col
        n=size(a,1); if(size(a,2)/=n .or. n<2) error stop "expm_cond_exact: square n>=2 required"
        allocate(k(n*n,n*n),e(n,n)); k=0.0_dp; e=0.0_dp; col=0
        do j=1,n
            do i=1,n
                e=0.0_dp; e(i,j)=1.0_dp; call expm_frechet_sps(a,e,x,l); col=col+1
                k(:,col)=reshape(l,[n*n])
            end do
        end do
        nkf=spectral_norm_real(k); nk1=norm1_real(k)
        condf=nkf*normf_real(a)/normf_real(x)
        cond1=nk1*norm1_real(a)/(norm1_real(x)*real(n,dp))
        if(present(expa)) expa=x
    end subroutine expm_cond_exact

    subroutine expm_cond_1_est(a,cond1,expa)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: cond1
        real(dp), allocatable, intent(out), optional :: expa(:,:)
        real(dp), allocatable :: e(:,:),v(:,:),z(:,:),x(:,:),ea(:,:),lt(:,:)
        real(dp) :: g,g2
        integer :: n,k,j,lidx,ii,jj
        n=size(a,1); if(size(a,2)/=n .or. n<2) error stop "expm_cond_1_est: square n>=2 required"
        allocate(e(n,n),z(n,n)); e=1.0_dp/real(n*n,dp)
        call expm_frechet_sps(a,e,ea,v); g=sum(abs(v)); z=sign(1.0_dp,v)
        call expm_frechet_sps(transpose(a),z,x,lt); x=lt
        k=2
        do
            j=maxloc(abs(reshape(x,[n*n])),dim=1)
            ii=mod(j-1,n)+1; jj=(j-1)/n+1; e=0.0_dp; e(ii,jj)=1.0_dp
            call expm_frechet_sps(a,e,lt,v); g=sum(abs(v))
            if(all(((v>=0.0_dp) .and. (z>=0.0_dp)) .or. ((v<0.0_dp) .and. (z<0.0_dp))) .or. &
               all(((v>=0.0_dp) .and. (z<0.0_dp)) .or. ((v<0.0_dp) .and. (z>=0.0_dp)))) exit
            z=sign(1.0_dp,v); call expm_frechet_sps(transpose(a),z,lt,x); x=x
            k=k+1; if(k>5) exit
        end do
        do lidx=1,n*n
            ii=mod(lidx-1,n)+1; jj=(lidx-1)/n+1
            x(ii,jj)=(-1.0_dp)**(lidx+1)*(1.0_dp+real(lidx-1,dp)/real(n*n-1,dp))
        end do
        call expm_frechet_sps(a,x,lt,v); g2=2.0_dp*sum(abs(v))/(3.0_dp*real(n*n,dp)); g=max(g,g2)
        cond1=g*norm1_real(a)/(norm1_real(ea)*real(n,dp)); if(present(expa)) expa=ea
    end subroutine expm_cond_1_est

    subroutine expm_cond_f_est(a,condf,iterations,expa,abstol,reltol,maxiter)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: condf
        integer, intent(out), optional :: iterations
        real(dp), allocatable, intent(out), optional :: expa(:,:)
        real(dp), intent(in), optional :: abstol,reltol
        integer, intent(in), optional :: maxiter
        real(dp), allocatable :: z1(:,:),w1(:,:),z2(:,:),w2(:,:),ea(:,:),tmp(:,:)
        real(dp) :: g1,g2,at,rt,dg
        integer :: n,it,mi,i,j
        n=size(a,1); at=0.1_dp; rt=1.0e-6_dp; mi=100
        if(present(abstol)) at=abstol; if(present(reltol)) rt=reltol; if(present(maxiter)) mi=maxiter
        allocate(z1(n,n)); do j=1,n; do i=1,n; z1(i,j)=sin(real(i+17*j,dp)); end do; end do
        call expm_frechet_sps(a,z1,ea,w1); call expm_frechet_sps(transpose(a),w1,tmp,z1)
        g2=normf_real(z1)/normf_real(w1); it=0
        do
            g1=g2; call expm_frechet_sps(a,z1,tmp,w2); call expm_frechet_sps(transpose(a),w2,tmp,z2)
            g2=normf_real(z2)/normf_real(w2); z1=z2; dg=abs(g1-g2); it=it+1
            if(it>mi .or. (dg<at .and. dg<rt*g2)) exit
        end do
        condf=g2*normf_real(a)/normf_real(ea); if(present(iterations)) iterations=it; if(present(expa)) expa=ea
    end subroutine expm_cond_f_est
end module expm_condition
