! Matrix exponential and integer powers translated from expm.
! GPL-3.0-or-later; see LICENSE and LICENSES.md.
module expm_matrix_functions
    use expm_kinds, only : dp
    use expm_linalg, only : eye_real, eye_complex, norm1_real, norm1_complex, &
        norminf_real, solve_real, solve_complex, balance_real, balance_complex, &
        balance_real_result, balance_complex_result, reverse_balance_real, &
        reverse_balance_complex
    implicit none
    private

    public :: expm, expm_higham08, expm_pade, expm_taylor, expm_almohy09, expm_rbs
    public :: expm_ward77, matrix_power

    interface expm
        module procedure expm_higham08_real
        module procedure expm_higham08_complex
    end interface expm
    interface expm_higham08
        module procedure expm_higham08_real
        module procedure expm_higham08_complex
    end interface expm_higham08
    interface expm_pade
        module procedure expm_pade_real
        module procedure expm_pade_complex
    end interface expm_pade
    interface expm_taylor
        module procedure expm_taylor_real
        module procedure expm_taylor_complex
    end interface expm_taylor
    interface matrix_power
        module procedure matrix_power_real
        module procedure matrix_power_complex
    end interface matrix_power
contains
    function expm_higham08_real(a,balancing) result(x)
        real(dp), intent(in) :: a(:,:)
        logical, intent(in), optional :: balancing
        real(dp), allocatable :: x(:,:),b(:,:)
        type(balance_real_result) :: bp,bs
        logical :: dobal
        integer :: n
        n=size(a,1); if(size(a,2)/=n) error stop "expm: matrix must be square"
        allocate(b(n,n))
        if(n==0) then; allocate(x(0,0)); return; end if
        if(n==1) then; allocate(x(1,1)); x(1,1)=exp(a(1,1)); return; end if
        dobal=.true.; if(present(balancing)) dobal=balancing
        if(dobal) then
            bp=balance_real(a,'P'); if(bp%info/=0) error stop "expm: dgebal(P) failed"
            bs=balance_real(bp%z,'S'); if(bs%info/=0) error stop "expm: dgebal(S) failed"
            b=bs%z
        else
            b=a
        end if
        x=pade_higham_real(b)
        if(dobal) call reverse_balance_real(x,bp,bs)
    end function expm_higham08_real

    function expm_higham08_complex(a,balancing) result(x)
        complex(dp), intent(in) :: a(:,:)
        logical, intent(in), optional :: balancing
        complex(dp), allocatable :: x(:,:),b(:,:)
        type(balance_complex_result) :: bp,bs
        logical :: dobal
        integer :: n
        n=size(a,1); if(size(a,2)/=n) error stop "expm: matrix must be square"
        allocate(b(n,n))
        if(n==0) then; allocate(x(0,0)); return; end if
        if(n==1) then; allocate(x(1,1)); x(1,1)=exp(a(1,1)); return; end if
        dobal=.true.; if(present(balancing)) dobal=balancing
        if(dobal) then
            bp=balance_complex(a,'P'); if(bp%info/=0) error stop "expm: zgebal(P) failed"
            bs=balance_complex(bp%z,'S'); if(bs%info/=0) error stop "expm: zgebal(S) failed"
            b=bs%z
        else
            b=a
        end if
        x=pade_higham_complex(b)
        if(dobal) call reverse_balance_complex(x,bp,bs)
    end function expm_higham08_complex

    function pade_higham_real(a) result(x)
        real(dp), intent(in) :: a(:,:)
        real(dp), allocatable :: x(:,:),aa(:,:),a2(:,:),a4(:,:),a6(:,:),u(:,:),v(:,:),p(:,:),id(:,:)
        real(dp) :: na
        real(dp), parameter :: c13(14)=[ &
            64764752532480000.0_dp,32382376266240000.0_dp,7771770303897600.0_dp, &
            1187353796428800.0_dp,129060195264000.0_dp,10559470521600.0_dp, &
            670442572800.0_dp,33522128640.0_dp,1323241920.0_dp,40840800.0_dp, &
            960960.0_dp,16380.0_dp,182.0_dp,1.0_dp]
        real(dp), parameter :: theta13=5.4_dp
        real(dp) :: cc(10)
        integer :: n,l,k,s
        n=size(a,1); id=eye_real(n); na=norm1_real(a); aa=a
        if(na<=2.1_dp) then
            if(na<=0.015_dp) then
                l=1; cc=0.0_dp; cc(1:4)=[120.0_dp,60.0_dp,12.0_dp,1.0_dp]
            else if(na<=0.25_dp) then
                l=2; cc=0.0_dp; cc(1:6)=[30240.0_dp,15120.0_dp,3360.0_dp,420.0_dp,30.0_dp,1.0_dp]
            else if(na<=0.95_dp) then
                l=3; cc=0.0_dp; cc(1:8)=[17297280.0_dp,8648640.0_dp,1995840.0_dp,277200.0_dp, &
                    25200.0_dp,1512.0_dp,56.0_dp,1.0_dp]
            else
                l=4; cc=[17643225600.0_dp,8821612800.0_dp,2075673600.0_dp,302702400.0_dp, &
                    30270240.0_dp,2162160.0_dp,110880.0_dp,3960.0_dp,90.0_dp,1.0_dp]
            end if
            a2=matmul(aa,aa); p=id; u=cc(2)*id; v=cc(1)*id
            do k=1,l
                p=matmul(p,a2)
                u=u+cc(2*k+2)*p
                v=v+cc(2*k+1)*p
            end do
            u=matmul(aa,u)
            x=solve_real(v-u,v+u)
            return
        end if
        s=max(0,ceiling(log(na/theta13)/log(2.0_dp)))
        if(s>0) aa=aa/(2.0_dp**s)
        a2=matmul(aa,aa); a4=matmul(a2,a2); a6=matmul(a2,a4)
        u=matmul(aa,matmul(a6,c13(14)*a6+c13(12)*a4+c13(10)*a2) + &
            c13(8)*a6+c13(6)*a4+c13(4)*a2+c13(2)*id)
        v=matmul(a6,c13(13)*a6+c13(11)*a4+c13(9)*a2) + &
            c13(7)*a6+c13(5)*a4+c13(3)*a2+c13(1)*id
        x=solve_real(v-u,v+u)
        do k=1,s; x=matmul(x,x); end do
    end function pade_higham_real

    function pade_higham_complex(a) result(x)
        complex(dp), intent(in) :: a(:,:)
        complex(dp), allocatable :: x(:,:),aa(:,:),a2(:,:),a4(:,:),a6(:,:),u(:,:),v(:,:),p(:,:),id(:,:)
        real(dp) :: na
        real(dp), parameter :: c13(14)=[ &
            64764752532480000.0_dp,32382376266240000.0_dp,7771770303897600.0_dp, &
            1187353796428800.0_dp,129060195264000.0_dp,10559470521600.0_dp, &
            670442572800.0_dp,33522128640.0_dp,1323241920.0_dp,40840800.0_dp, &
            960960.0_dp,16380.0_dp,182.0_dp,1.0_dp]
        real(dp), parameter :: theta13=5.4_dp
        real(dp) :: cc(10)
        integer :: n,l,k,s
        n=size(a,1); id=eye_complex(n); na=norm1_complex(a); aa=a
        if(na<=2.1_dp) then
            if(na<=0.015_dp) then
                l=1; cc=0.0_dp; cc(1:4)=[120.0_dp,60.0_dp,12.0_dp,1.0_dp]
            else if(na<=0.25_dp) then
                l=2; cc=0.0_dp; cc(1:6)=[30240.0_dp,15120.0_dp,3360.0_dp,420.0_dp,30.0_dp,1.0_dp]
            else if(na<=0.95_dp) then
                l=3; cc=0.0_dp; cc(1:8)=[17297280.0_dp,8648640.0_dp,1995840.0_dp,277200.0_dp, &
                    25200.0_dp,1512.0_dp,56.0_dp,1.0_dp]
            else
                l=4; cc=[17643225600.0_dp,8821612800.0_dp,2075673600.0_dp,302702400.0_dp, &
                    30270240.0_dp,2162160.0_dp,110880.0_dp,3960.0_dp,90.0_dp,1.0_dp]
            end if
            a2=matmul(aa,aa); p=id; u=cc(2)*id; v=cc(1)*id
            do k=1,l
                p=matmul(p,a2); u=u+cc(2*k+2)*p; v=v+cc(2*k+1)*p
            end do
            u=matmul(aa,u); x=solve_complex(v-u,v+u); return
        end if
        s=max(0,ceiling(log(na/theta13)/log(2.0_dp)))
        if(s>0) aa=aa/(2.0_dp**s)
        a2=matmul(aa,aa); a4=matmul(a2,a2); a6=matmul(a2,a4)
        u=matmul(aa,matmul(a6,c13(14)*a6+c13(12)*a4+c13(10)*a2) + &
            c13(8)*a6+c13(6)*a4+c13(4)*a2+c13(2)*id)
        v=matmul(a6,c13(13)*a6+c13(11)*a4+c13(9)*a2) + &
            c13(7)*a6+c13(5)*a4+c13(3)*a2+c13(1)*id
        x=solve_complex(v-u,v+u)
        do k=1,s; x=matmul(x,x); end do
    end function pade_higham_complex

    function expm_pade_real(a,order) result(x)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in), optional :: order
        real(dp), allocatable :: x(:,:),as(:,:),num(:,:),den(:,:),pow(:,:),id(:,:)
        real(dp) :: c,rows
        integer :: n,m,k,s
        n=size(a,1); if(size(a,2)/=n) error stop "expm_pade: matrix must be square"
        m=8; if(present(order)) m=order; if(m<2) error stop "expm_pade: order must be >=2"
        rows=0.0_dp; do k=1,n; rows=max(rows,sum(abs(a(k,:)))); end do
        if(rows>0.0_dp) then; s=max(0,ceiling(log(rows)/log(2.0_dp))+1); else; s=0; end if
        as=a/(2.0_dp**s); id=eye_real(n); c=0.5_dp; num=id+c*as; den=id-c*as; pow=as
        do k=2,m
            c=c*real(m-k+1,dp)/(real(k,dp)*real(m-k+1+m,dp))
            pow=matmul(as,pow); num=num+c*pow
            if(mod(k,2)==0) then; den=den+c*pow; else; den=den-c*pow; end if
        end do
        x=solve_real(den,num); do k=1,s; x=matmul(x,x); end do
    end function expm_pade_real

    function expm_pade_complex(a,order) result(x)
        complex(dp), intent(in) :: a(:,:)
        integer, intent(in), optional :: order
        complex(dp), allocatable :: x(:,:),as(:,:),num(:,:),den(:,:),pow(:,:),id(:,:)
        real(dp) :: c,rows
        integer :: n,m,k,s
        n=size(a,1); if(size(a,2)/=n) error stop "expm_pade: matrix must be square"
        m=8; if(present(order)) m=order; if(m<2) error stop "expm_pade: order must be >=2"
        rows=0.0_dp; do k=1,n; rows=max(rows,sum(abs(a(k,:)))); end do
        if(rows>0.0_dp) then; s=max(0,ceiling(log(rows)/log(2.0_dp))+1); else; s=0; end if
        as=a/(2.0_dp**s); id=eye_complex(n); c=0.5_dp; num=id+c*as; den=id-c*as; pow=as
        do k=2,m
            c=c*real(m-k+1,dp)/(real(k,dp)*real(m-k+1+m,dp))
            pow=matmul(as,pow); num=num+c*pow
            if(mod(k,2)==0) then; den=den+c*pow; else; den=den-c*pow; end if
        end do
        x=solve_complex(den,num); do k=1,s; x=matmul(x,x); end do
    end function expm_pade_complex

    function expm_taylor_real(a,order) result(x)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in), optional :: order
        real(dp), allocatable :: x(:,:),as(:,:),term(:,:)
        real(dp) :: na
        integer :: n,m,k,s
        n=size(a,1); m=18; if(present(order)) m=order; if(m<1) error stop "expm_taylor: order must be >=1"
        na=sqrt(norm1_real(a)*norminf_real(a));
        if(na>0.0_dp) then; s=max(0,int(log(na)/log(2.0_dp))+4); else; s=0; end if
        as=a/(2.0_dp**s); x=eye_real(n); term=x
        do k=1,m; term=matmul(term,as)/real(k,dp); x=x+term; end do
        do k=1,s; x=matmul(x,x); end do
    end function expm_taylor_real

    function expm_taylor_complex(a,order) result(x)
        complex(dp), intent(in) :: a(:,:)
        integer, intent(in), optional :: order
        complex(dp), allocatable :: x(:,:),as(:,:),term(:,:)
        real(dp) :: na,ni
        integer :: n,m,k,s,i
        n=size(a,1); m=18; if(present(order)) m=order; if(m<1) error stop "expm_taylor: order must be >=1"
        ni=0.0_dp; do i=1,n; ni=max(ni,sum(abs(a(i,:)))); end do
        na=sqrt(norm1_complex(a)*ni)
        if(na>0.0_dp) then; s=max(0,int(log(na)/log(2.0_dp))+4); else; s=0; end if
        as=a/(2.0_dp**s); x=eye_complex(n); term=x
        do k=1,m; term=matmul(term,as)/real(k,dp); x=x+term; end do
        do k=1,s; x=matmul(x,x); end do
    end function expm_taylor_complex

    function expm_almohy09(a,p) result(x)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in), optional :: p
        real(dp), allocatable :: x(:,:),as(:,:)
        integer :: degree,sfactor,k,s
        real(dp) :: na
        degree=6; if(present(p)) degree=p
        if(degree<1 .or. degree>13) error stop "expm_almohy09: p must be 1..13"
        na=norm1_real(a); sfactor=1
        if(na>5.4_dp) then
            s=ceiling(log(na/5.4_dp)/log(2.0_dp)); sfactor=2**s
        end if
        as=a/real(sfactor,dp); x=pade_fixed_real(as,degree)
        if(sfactor>1) then
            s=nint(log(real(sfactor,dp))/log(2.0_dp))
            do k=1,s; x=matmul(x,x); end do
        end if
    end function expm_almohy09

    function pade_fixed_real(a,p) result(x)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in) :: p
        real(dp), parameter :: coef(0:13)=[1.0_dp,0.5_dp,0.12_dp, &
            1.8333333333333333e-2_dp,1.9927536231884058e-3_dp,1.6304347826086957e-4_dp, &
            1.0351966873706004e-5_dp,5.1759834368530021e-7_dp,2.0431513566525008e-8_dp, &
            6.3060227057175951e-10_dp,1.4837700484041400e-11_dp,2.5291534915979660e-13_dp, &
            2.8101705462199622e-15_dp,1.5440497506703089e-17_dp]
        real(dp), allocatable :: x(:,:),num(:,:),den(:,:),pow(:,:),id(:,:)
        integer :: k,n
        n=size(a,1); id=eye_real(n); num=id; den=id; pow=id
        do k=1,p
            pow=matmul(pow,a); num=num+coef(k)*pow
            if(mod(k,2)==0) then; den=den+coef(k)*pow; else; den=den-coef(k)*pow; end if
        end do
        x=solve_real(den,num)
    end function pade_fixed_real

    function expm_ward77(a,order) result(x)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in), optional :: order
        real(dp), allocatable :: x(:,:),b(:,:)
        real(dp) :: shift
        integer :: n,m,i
        type(balance_real_result) :: bp,bs
        n=size(a,1); m=8; if(present(order)) m=order
        shift=0.0_dp; do i=1,n; shift=shift+a(i,i); end do; shift=shift/real(n,dp)
        b=a; do i=1,n; b(i,i)=b(i,i)-shift; end do
        bp=balance_real(b,'P'); bs=balance_real(bp%z,'S')
        x=expm_pade_real(bs%z,m); call reverse_balance_real(x,bp,bs); x=exp(shift)*x
    end function expm_ward77


    function expm_rbs(a,degree,t) result(x)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in), optional :: degree
        real(dp), intent(in), optional :: t
        real(dp), allocatable :: x(:,:),as(:,:),num(:,:),den(:,:),pow(:,:),id(:,:),coef(:)
        real(dp) :: tt,hnorm,scale
        integer :: n,m,ns,k
        n=size(a,1); if(size(a,2)/=n) error stop "expm_rbs: matrix must be square"
        m=6; if(present(degree)) m=degree; if(m<1) error stop "expm_rbs: degree must be >=1"
        tt=1.0_dp; if(present(t)) tt=t
        hnorm=abs(tt)*norm1_real(a)
        if(hnorm<=tiny(1.0_dp)) then; x=eye_real(n); return; end if
        ns=max(0,int(log(hnorm)/log(2.0_dp))+2); scale=tt/(2.0_dp**ns); as=scale*a
        allocate(coef(0:m)); coef(0)=1.0_dp
        do k=1,m
            coef(k)=coef(k-1)*real(m+1-k,dp)/real(k*(2*m+1-k),dp)
        end do
        id=eye_real(n); num=id; den=id; pow=id
        do k=1,m
            pow=matmul(pow,as); num=num+coef(k)*pow
            if(mod(k,2)==0) then; den=den+coef(k)*pow; else; den=den-coef(k)*pow; end if
        end do
        x=solve_real(den,num); do k=1,ns; x=matmul(x,x); end do
    end function expm_rbs

    function matrix_power_real(a,k) result(x)
        real(dp), intent(in) :: a(:,:)
        integer, intent(in) :: k
        real(dp), allocatable :: x(:,:),b(:,:)
        integer :: e,n
        if(k<0) error stop "matrix_power: exponent must be nonnegative"
        n=size(a,1); if(size(a,2)/=n) error stop "matrix_power: matrix must be square"
        x=eye_real(n); if(k==0) return; b=a; e=k
        do while(e>0)
            if(iand(e,1)==1) x=matmul(x,b)
            e=ishft(e,-1); if(e>0) b=matmul(b,b)
        end do
    end function matrix_power_real

    function matrix_power_complex(a,k) result(x)
        complex(dp), intent(in) :: a(:,:)
        integer, intent(in) :: k
        complex(dp), allocatable :: x(:,:),b(:,:)
        integer :: e,n
        if(k<0) error stop "matrix_power: exponent must be nonnegative"
        n=size(a,1); if(size(a,2)/=n) error stop "matrix_power: matrix must be square"
        x=eye_complex(n); if(k==0) return; b=a; e=k
        do while(e>0)
            if(iand(e,1)==1) x=matmul(x,b)
            e=ishft(e,-1); if(e>0) b=matmul(b,b)
        end do
    end function matrix_power_complex
end module expm_matrix_functions
