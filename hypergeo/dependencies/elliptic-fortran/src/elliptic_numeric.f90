! SPDX-License-Identifier: GPL-2.0-only
module elliptic_numeric
    use elliptic_kinds, only : dp, pi, ci, finite_complex
    implicit none
    private
    public :: complex_function, integrate_real_path, integrate_segments, residue
    public :: newton_raphson_complex

    abstract interface
        function complex_function(z) result(fz)
            import dp
            complex(dp), intent(in) :: z
            complex(dp) :: fz
        end function complex_function
    end interface
contains
    function integrate_real_path(f, a, b, tol, maxdepth) result(val)
        procedure(complex_function) :: f
        real(dp),intent(in)::a,b
        real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxdepth
        complex(dp)::val,fa,fb,fm,s
        real(dp)::t
        integer::md
        t=1.0e-10_dp;if(present(tol))t=tol
        md=20;if(present(maxdepth))md=maxdepth
        fa=f(cmplx(a,0.0_dp,dp));fb=f(cmplx(b,0.0_dp,dp));fm=f(cmplx(0.5_dp*(a+b),0.0_dp,dp))
        s=(b-a)*(fa+4.0_dp*fm+fb)/6.0_dp
        val=adapt(f,a,b,fa,fm,fb,s,t,md)
    end function integrate_real_path

    recursive function adapt(f,a,b,fa,fm,fb,s,tol,depth) result(v)
        procedure(complex_function) :: f
        real(dp),intent(in)::a,b,tol
        complex(dp),intent(in)::fa,fm,fb,s
        integer,intent(in)::depth
        complex(dp)::v,fl,fr,sl,sr,s2
        real(dp)::m,l,r
        m=0.5_dp*(a+b);l=0.5_dp*(a+m);r=0.5_dp*(m+b)
        fl=f(cmplx(l,0.0_dp,dp));fr=f(cmplx(r,0.0_dp,dp))
        sl=(m-a)*(fa+4.0_dp*fl+fm)/6.0_dp
        sr=(b-m)*(fm+4.0_dp*fr+fb)/6.0_dp
        s2=sl+sr
        if(depth<=0 .or. abs(s2-s)<=15.0_dp*tol) then
            v=s2+(s2-s)/15.0_dp
        else
            v=adapt(f,a,m,fa,fl,fm,sl,0.5_dp*tol,depth-1)+ &
              adapt(f,m,b,fm,fr,fb,sr,0.5_dp*tol,depth-1)
        end if
    end function adapt

    function integrate_segments(f, points, close_path, nsub) result(val)
        procedure(complex_function) :: f
        complex(dp),intent(in)::points(:)
        logical,intent(in),optional::close_path
        integer,intent(in),optional::nsub
        complex(dp)::val,z0,z1,z,dz,segsum
        logical::cl
        integer::i,j,nseg,ns
        real(dp)::t,h,w
        cl=.true.;if(present(close_path))cl=close_path
        ns=200;if(present(nsub))ns=max(2,nsub);if(mod(ns,2)==1)ns=ns+1
        val=(0.0_dp,0.0_dp);nseg=size(points)-1
        if(cl)nseg=size(points)
        do i=1,nseg
            z0=points(i);if(i<size(points))then;z1=points(i+1);else;z1=points(1);end if
            dz=z1-z0
            h=1.0_dp/real(ns,dp)
            segsum=(0.0_dp,0.0_dp)
            do j=0,ns
                t=real(j,dp)*h
                z=z0+dz*t
                if(j==0.or.j==ns)then
                    w=1.0_dp
                else if(mod(j,2)==0)then
                    w=2.0_dp
                else
                    w=4.0_dp
                end if
                segsum=segsum+w*f(z)*dz
            end do
            val=val+segsum*h/3.0_dp
        end do
    end function integrate_segments

    function residue(f,z0,r,center,nsub) result(res)
        procedure(complex_function)::f
        complex(dp),intent(in)::z0
        real(dp),intent(in)::r
        complex(dp),intent(in),optional::center
        integer,intent(in),optional::nsub
        complex(dp)::res,o,z,dz,integrand
        integer::j,ns
        real(dp)::t,h,w
        o=z0;if(present(center))o=center
        ns=1000;if(present(nsub))ns=max(20,nsub);if(mod(ns,2)==1)ns=ns+1
        h=1.0_dp/real(ns,dp);res=(0.0_dp,0.0_dp)
        do j=0,ns
            t=real(j,dp)*h
            z=o+r*exp(2.0_dp*pi*ci*t)
            dz=r*2.0_dp*pi*ci*exp(2.0_dp*pi*ci*t)
            integrand=f(z)/(z-z0)*dz
            if(j==0.or.j==ns)then;w=1.0_dp;else if(mod(j,2)==0)then;w=2.0_dp;else;w=4.0_dp;end if
            res=res+w*integrand
        end do
        res=res*h/3.0_dp/(2.0_dp*pi*ci)
    end function residue

    function newton_raphson_complex(initial,f,fdash,maxiter,tol,converged,iterations) result(root)
        complex(dp),intent(in)::initial
        procedure(complex_function)::f,fdash
        integer,intent(in),optional::maxiter
        real(dp),intent(in),optional::tol
        logical,intent(out),optional::converged
        integer,intent(out),optional::iterations
        complex(dp)::root,old,newv,fv
        integer::i,nmax
        real(dp)::t
        nmax=100;if(present(maxiter))nmax=maxiter
        t=sqrt(epsilon(1.0_dp));if(present(tol))t=tol
        old=initial
        do i=1,nmax
            newv=old-f(old)/fdash(old);fv=f(newv)
            if(.not.finite_complex(fv))exit
            if(abs(newv-old)<=t*max(1.0_dp,abs(newv)).or.abs(fv)<=t)then
                root=newv;if(present(converged))converged=.true.;if(present(iterations))iterations=i;return
            end if
            old=newv
        end do
        root=old;if(present(converged))converged=.false.;if(present(iterations))iterations=i-1
    end function newton_raphson_complex
end module elliptic_numeric
