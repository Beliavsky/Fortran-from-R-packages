! Sidje/EXPOKIT exp(t*A)*v action algorithm translated from expm::expAtv.
! GPL-3.0-or-later; see LICENSE and LICENSES.md.
module expm_action
    use expm_kinds, only : dp
    use expm_linalg, only : norminf_real
    use expm_matrix_functions, only : expm
    implicit none
    private
    type, public :: exp_action_result
        real(dp), allocatable :: value(:)
        real(dp) :: error = 0.0_dp
        integer :: nstep = 0
        integer :: nreject = 0
        integer :: info = 0
    end type exp_action_result
    public :: exp_at_v
contains
    function exp_at_v(a,v,t,tol,btol,m_max,mxrej,rescale_below) result(r)
        real(dp), intent(in) :: a(:,:),v(:)
        real(dp), intent(in), optional :: t,tol,btol,rescale_below
        integer, intent(in), optional :: m_max,mxrej
        type(exp_action_result) :: r
        real(dp), allocatable :: aa(:,:),vv(:,:),h(:,:),w(:),p(:),av(:),f(:,:)
        real(dp) :: tt,tl,bl,resc,gamma,delta,na,rndoff,t1,sgn,t_now,s_error,xm,beta,fact,t_new
        real(dp) :: t_step,s,avnorm,phi1,phi2,err_loc
        integer :: n,m,mm,mr,k1,mb,j,i,irej,mx
        n=size(a,1)
        if(size(a,2)/=n .or. size(v)/=n .or. n<1) error stop "exp_at_v: incompatible arguments"
        allocate(r%value(n)); r%value=0.0_dp
        tt=1.0_dp; tl=1.0e-7_dp; bl=1.0e-7_dp; mm=30; mr=10; resc=1.0e-6_dp
        if(present(t)) tt=t; if(present(tol)) tl=tol; if(present(btol)) bl=btol
        if(present(m_max)) mm=m_max; if(present(mxrej)) mr=mxrej; if(present(rescale_below)) resc=rescale_below
        if(n==1) then
            r%value(1)=exp(a(1,1)*tt)*v(1); return
        end if
        m=min(n,mm); if(m<2) error stop "exp_at_v: m_max must be >=2"
        gamma=0.9_dp; delta=1.2_dp; aa=a; na=norminf_real(aa)
        if(na>0.0_dp .and. na<resc) then; aa=aa/na; tt=tt*na; na=1.0_dp; end if
        rndoff=na*epsilon(1.0_dp); t1=abs(tt); sgn=sign(1.0_dp,tt); t_now=0.0_dp; s_error=0.0_dp
        k1=2; mb=m; xm=1.0_dp/real(m,dp); beta=sqrt(sum(v*v)); allocate(w(n)); w=v
        if(beta<=tiny(1.0_dp)) then; r%value=w; return; end if
        fact=((real(m+1,dp)/exp(1.0_dp))**(m+1))*sqrt(2.0_dp*acos(-1.0_dp)*real(m+1,dp))
        t_new=round_step((fact*tl/(4.0_dp*beta*na))**xm/na)
        allocate(vv(n,m+1),h(m+2,m+2),p(n),av(n)); vv=0.0_dp; h=0.0_dp
        do while(t_now<t1)
            r%nstep=r%nstep+1; t_step=min(t1-t_now,t_new); vv(:,1)=w/beta; k1=2; mb=m
            do j=1,m
                p=matmul(aa,vv(:,j))
                do i=1,j
                    h(i,j)=dot_product(vv(:,i),p); p=p-h(i,j)*vv(:,i)
                end do
                s=sqrt(dot_product(p,p))
                if(s<bl) then; k1=0; mb=j; t_step=t1-t_now; exit; end if
                h(j+1,j)=s; vv(:,j+1)=p/s
            end do
            if(k1/=0) then; h(m+2,m+1)=1.0_dp; av=matmul(aa,vv(:,m+1)); avnorm=sqrt(dot_product(av,av)); end if
            irej=0
            do while(irej<=mr)
                mx=mb+k1; f=expm(sgn*t_step*h(1:mx,1:mx),balancing=.false.)
                if(k1==0) then
                    err_loc=bl; exit
                end if
                phi1=abs(beta*f(m+1,1)); phi2=abs(beta*f(m+2,1)*avnorm)
                if(phi1>10.0_dp*phi2) then; err_loc=phi2; xm=1.0_dp/real(m,dp)
                else if(phi1>phi2) then; err_loc=(phi1*phi2)/(phi1-phi2); xm=1.0_dp/real(m,dp)
                else; err_loc=phi1; xm=1.0_dp/real(m-1,dp); end if
                if(err_loc<=delta*t_step*tl) exit
                if(irej==mr) then; r%info=1; r%value=w; return; end if
                t_step=round_step(gamma*t_step*(t_step*tl/err_loc)**xm); irej=irej+1
            end do
            r%nreject=r%nreject+irej; mx=mb+max(0,k1-1)
            w=matmul(vv(:,1:mx),beta*f(1:mx,1)); beta=sqrt(dot_product(w,w)); t_now=t_now+t_step
            if(err_loc>0.0_dp) t_new=round_step(gamma*t_step*(t_step*tl/err_loc)**xm)
            err_loc=max(err_loc,rndoff); s_error=s_error+err_loc
        end do
        r%value=w; r%error=s_error
    contains
        real(dp) function round_step(x)
            real(dp), intent(in) :: x
            real(dp) :: sc
            if(x<=0.0_dp) then; round_step=0.0_dp; return; end if
            sc=10.0_dp**(floor(log10(x))-1.0_dp); round_step=ceiling(x/sc)*sc
        end function round_step
    end function exp_at_v
end module expm_action
