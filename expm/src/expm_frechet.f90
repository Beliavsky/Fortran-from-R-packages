! Frechet derivative algorithms translated from expm's Higham 2008 implementation.
! GPL-3.0-or-later; see LICENSE and LICENSES.md.
module expm_frechet
    use expm_kinds, only : dp
    use expm_linalg, only : eye_real, norm1_real, solve_real
    implicit none
    private
    public :: expm_frechet_sps, expm_frechet_block
contains
    subroutine expm_frechet_sps(a,e,x,l)
        real(dp), intent(in) :: a(:,:),e(:,:)
        real(dp), allocatable, intent(out) :: x(:,:),l(:,:)
        real(dp), allocatable :: id(:,:),p(:,:),u(:,:),v(:,:),a2(:,:),m2(:,:),m(:,:),lu(:,:),lv(:,:)
        real(dp), allocatable :: b(:,:),d(:,:),b2(:,:),b4(:,:),b6(:,:),m4(:,:),m6(:,:)
        real(dp), allocatable :: w1(:,:),w2(:,:),z1(:,:),z2(:,:),w(:,:),lw1(:,:),lw2(:,:),lz1(:,:),lz2(:,:),lw(:,:)
        real(dp), parameter :: c13(14)=[ &
            64764752532480000.0_dp,32382376266240000.0_dp,7771770303897600.0_dp, &
            1187353796428800.0_dp,129060195264000.0_dp,10559470521600.0_dp, &
            670442572800.0_dp,33522128640.0_dp,1323241920.0_dp,40840800.0_dp, &
            960960.0_dp,16380.0_dp,182.0_dp,1.0_dp]
        real(dp) :: na,cc(10)
        integer :: n,ll,k,oc,s,t
        n=size(a,1)
        if(size(a,2)/=n .or. any(shape(e)/=shape(a))) error stop "expm_frechet: incompatible matrices"
        if(n==1) then
            allocate(x(1,1),l(1,1)); x(1,1)=exp(a(1,1)); l(1,1)=e(1,1)*x(1,1); return
        end if
        id=eye_real(n); na=norm1_real(a)
        if(na<=1.78_dp) then
            if(na<=0.0108_dp) then
                ll=1; cc=0.0_dp; cc(1:4)=[120.0_dp,60.0_dp,12.0_dp,1.0_dp]
            else if(na<=0.2_dp) then
                ll=2; cc=0.0_dp; cc(1:6)=[30240.0_dp,15120.0_dp,3360.0_dp,420.0_dp,30.0_dp,1.0_dp]
            else if(na<=0.783_dp) then
                ll=3; cc=0.0_dp; cc(1:8)=[17297280.0_dp,8648640.0_dp,1995840.0_dp,277200.0_dp, &
                    25200.0_dp,1512.0_dp,56.0_dp,1.0_dp]
            else
                ll=4; cc=[17643225600.0_dp,8821612800.0_dp,2075673600.0_dp,302702400.0_dp, &
                    30270240.0_dp,2162160.0_dp,110880.0_dp,3960.0_dp,90.0_dp,1.0_dp]
            end if
            p=id; u=cc(2)*id; v=cc(1)*id
            a2=matmul(a,a); m2=matmul(a,e)+matmul(e,a); m=m2
            lu=cc(4)*m; lv=cc(3)*m; oc=2
            do k=1,ll-1
                p=matmul(p,a2); u=u+cc(oc+2)*p; v=v+cc(oc+1)*p
                m=matmul(a2,m)+matmul(m2,p)
                lu=lu+cc(oc+4)*m; lv=lv+cc(oc+3)*m; oc=oc+2
            end do
            p=matmul(p,a2); u=u+cc(oc+2)*p
            lu=matmul(a,lu)+matmul(e,u); u=matmul(a,u); v=v+cc(oc+1)*p
            x=solve_real(v-u,v+u)
            l=solve_real(v-u,lu+lv+matmul(lu-lv,x))
        else
            s=max(0,ceiling(log(na/4.74_dp)/log(2.0_dp)))
            b=a/(2.0_dp**s); d=e/(2.0_dp**s)
            b2=matmul(b,b); b4=matmul(b2,b2); b6=matmul(b2,b4)
            w1=c13(14)*b6+c13(12)*b4+c13(10)*b2
            w2=c13(8)*b6+c13(6)*b4+c13(4)*b2+c13(2)*id
            z1=c13(13)*b6+c13(11)*b4+c13(9)*b2
            z2=c13(7)*b6+c13(5)*b4+c13(3)*b2+c13(1)*id
            w=matmul(b6,w1)+w2; u=matmul(b,w); v=matmul(b6,z1)+z2
            m2=matmul(b,d)+matmul(d,b)
            m4=matmul(b2,m2)+matmul(m2,b2)
            m6=matmul(b4,m2)+matmul(m4,b2)
            lw1=c13(14)*m6+c13(12)*m4+c13(10)*m2
            lw2=c13(8)*m6+c13(6)*m4+c13(4)*m2
            lz1=c13(13)*m6+c13(11)*m4+c13(9)*m2
            lz2=c13(7)*m6+c13(5)*m4+c13(3)*m2
            lw=matmul(b6,lw1)+matmul(m6,w1)+lw2
            lu=matmul(b,lw)+matmul(d,w)
            lv=matmul(b6,lz1)+matmul(m6,z1)+lz2
            x=solve_real(v-u,v+u)
            l=solve_real(v-u,lu+lv+matmul(lu-lv,x))
            do t=1,s
                l=matmul(l,x)+matmul(x,l)
                x=matmul(x,x)
            end do
        end if
    end subroutine expm_frechet_sps

    subroutine expm_frechet_block(a,e,x,l)
        use expm_matrix_functions, only : expm
        real(dp), intent(in) :: a(:,:),e(:,:)
        real(dp), allocatable, intent(out) :: x(:,:),l(:,:)
        real(dp), allocatable :: b(:,:),fb(:,:)
        integer :: n
        n=size(a,1); if(size(a,2)/=n .or. any(shape(e)/=shape(a))) error stop "expm_frechet_block: incompatible matrices"
        allocate(b(2*n,2*n)); b=0.0_dp
        b(1:n,1:n)=a; b(1:n,n+1:2*n)=e; b(n+1:2*n,n+1:2*n)=a
        fb=expm(b,balancing=.false.)
        x=fb(1:n,1:n); l=fb(1:n,n+1:2*n)
    end subroutine expm_frechet_block
end module expm_frechet
