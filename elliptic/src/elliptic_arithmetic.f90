! SPDX-License-Identifier: GPL-2.0-only
module elliptic_arithmetic
    use elliptic_kinds, only : dp, i8
    implicit none
    private
    public :: primes_upto, factorize, divisor_sigma, totient, mobius, liouville
    public :: congruence_solutions
contains
    function primes_upto(n) result(p)
        integer, intent(in) :: n
        integer, allocatable :: p(:)
        logical, allocatable :: sieve(:)
        integer :: i, j, k
        if (n < 2) then
            allocate(p(0)); return
        end if
        allocate(sieve(n)); sieve=.true.; sieve(1)=.false.
        do i=2,int(sqrt(real(n,dp)))
            if(sieve(i)) then
                do j=i*i,n,i
                    sieve(j)=.false.
                end do
            end if
        end do
        k=count(sieve); allocate(p(k)); k=0
        do i=2,n
            if(sieve(i)) then; k=k+1; p(k)=i; end if
        end do
    end function primes_upto

    function factorize(n) result(fac)
        integer(i8), intent(in) :: n
        integer(i8), allocatable :: fac(:)
        integer(i8), allocatable :: tmp(:)
        integer(i8) :: x, d
        integer :: k
        if(n<2_i8) then; allocate(fac(0)); return; end if
        allocate(tmp(64)); k=0; x=n; d=2
        do while(d*d<=x)
            do while(mod(x,d)==0_i8)
                k=k+1; tmp(k)=d; x=x/d
            end do
            if(d==2) then; d=3; else; d=d+2; end if
        end do
        if(x>1) then; k=k+1; tmp(k)=x; end if
        allocate(fac(k)); if(k>0) fac=tmp(:k)
    end function factorize

    function divisor_sigma(n,kpower) result(s)
        integer(i8), intent(in) :: n
        integer, intent(in), optional :: kpower
        integer(i8) :: s
        integer :: kp, i, j, cnt
        integer(i8), allocatable :: f(:)
        integer(i8) :: p, term, block
        kp=1; if(present(kpower)) kp=kpower
        if(n==1_i8) then; s=1_i8; return; end if
        f=factorize(n); s=1_i8; i=1
        do while(i<=size(f))
            p=f(i); cnt=1; j=i+1
            do while(j<=size(f))
                if(f(j)/=p) exit
                cnt=cnt+1; j=j+1
            end do
            if(kp==0) then
                block=int(cnt+1,i8)
            else
                block=1_i8; term=1_i8
                do j=1,cnt
                    term=term*p**kp; block=block+term
                end do
            end if
            s=s*block; i=i+cnt
        end do
    end function divisor_sigma

    function totient(n) result(phi)
        integer(i8),intent(in)::n
        integer(i8)::phi,p
        integer(i8),allocatable::f(:)
        integer::i
        if(n<=1_i8) then; phi=max(n,0_i8); return; end if
        f=factorize(n); phi=n; i=1
        do while(i<=size(f))
            p=f(i); phi=phi/p*(p-1)
            do while(i<=size(f))
                if(f(i)/=p) exit
                i=i+1
            end do
        end do
    end function totient

    function mobius(n) result(mu)
        integer(i8),intent(in)::n
        integer::mu,i
        integer(i8),allocatable::f(:)
        if(n==1_i8) then;mu=1;return;end if
        f=factorize(n)
        do i=2,size(f)
            if(f(i)==f(i-1)) then;mu=0;return;end if
        end do
        if(mod(size(f),2)==0) then;mu=1;else;mu=-1;end if
    end function mobius

    function liouville(n) result(lambda_l)
        integer(i8),intent(in)::n
        integer::lambda_l
        integer(i8),allocatable::f(:)
        if(n==1_i8) then;lambda_l=1;return;end if
        f=factorize(n)
        if(mod(size(f),2)==0) then;lambda_l=1;else;lambda_l=-1;end if
    end function liouville

    subroutine congruence_solutions(a,l,out,nrow)
        integer,intent(in)::a(2)
        integer,intent(in),optional::l
        integer,allocatable,intent(out)::out(:,:)
        integer,intent(out)::nrow
        integer::ell,m,n,i,q1,q2
        ell=1;if(present(l))ell=l;m=a(1);n=a(2)
        allocate(out(3,2)); out=0; nrow=1; out(1,:)=a
        if(m==0.and.n==1) then;nrow=0;deallocate(out);allocate(out(0,2));return;end if
        if(m==1.and.n==0) then;nrow=1;out(1,:)=[huge(1),1];return;end if
        if(m==1) then;nrow=3;out(2,:)=[1,n+ell];out(3,:)=[0,ell];return;end if
        if(n==1) then;nrow=3;out(2,:)=[m-ell,1];out(3,:)=[ell,0];return;end if
        q1=huge(1); q2=huge(1)
        if(m>0) then
            do i=1,m;if(mod(ell+i*n,m)==0) then;q1=i;exit;end if;end do
        end if
        if(n>0) then
            do i=1,n;if(mod(-ell+i*m,n)==0) then;q2=i;exit;end if;end do
        end if
        nrow=2; out(2,:)=[q1,q2]
    end subroutine congruence_solutions
end module elliptic_arithmetic
