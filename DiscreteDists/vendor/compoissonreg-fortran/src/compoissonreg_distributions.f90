module compoissonreg_distributions
   use compoissonreg_kinds, only : dp
   use compoissonreg_types, only : cmp_control_t
   use compoissonreg_numerics, only : logadd
   use compoissonreg_normalizer, only : cmp_logz_hybrid, cmp_logz_trunc, cmp_y_trunc
   implicit none
   private
   public :: dcmp, pcmp, qcmp, rcmp, rcmp_vec, ecmp, vcmp, ncmp, tcmp
   public :: dzicmp, pzicmp, qzicmp, rzicmp, rzicmp_vec, ezicmp, vzicmp
   public :: dzip, pzip, qzip, rzip, ezip, vzip
   public :: loglik_cmp, loglik_zicmp, cmp_stats

contains

   subroutine cmp_stats(lambda,nu,control,logz,mean,var,elogfact)
      real(dp),intent(in)::lambda,nu
      type(cmp_control_t),intent(in),optional::control
      real(dp),intent(out)::logz,mean,var,elogfact
      type(cmp_control_t)::ctrl
      real(dp)::a,llambda,lp,prob,lf,m2
      integer::j,m
      logical::ok
      ctrl=cmp_control_t();if(present(control))ctrl=control
      if(lambda<0.0_dp.or.nu<0.0_dp)then
         logz=huge(1.0_dp);mean=huge(1.0_dp);var=huge(1.0_dp);elogfact=huge(1.0_dp);return
      end if
      if(lambda==0.0_dp)then
         logz=0.0_dp;mean=0.0_dp;var=0.0_dp;elogfact=0.0_dp;return
      end if
      if(nu==0.0_dp.and.lambda<1.0_dp)then
         logz=-log(1.0_dp-lambda);mean=lambda/(1.0_dp-lambda)
         var=lambda/(1.0_dp-lambda)**2
         ! No short closed form is used for E log(Y!); compute it by truncation.
         elogfact=0.0_dp;prob=1.0_dp-lambda;lf=0.0_dp
         do j=1,min(ctrl%ymax,100000)
            prob=prob*lambda;lf=lf+log(real(j,dp));elogfact=elogfact+prob*lf
            if(prob*max(1.0_dp,lf)<ctrl%truncate_tol)exit
         end do
         return
      end if
      if(nu>0.0_dp.and.nu/=1.0_dp.and.-log(lambda)/nu<log(ctrl%hybrid_tol))then
         llambda=log(lambda);a=exp(llambda/nu)
         logz=nu*a-(nu-1.0_dp)/(2.0_dp*nu)*llambda &
            -0.5_dp*(nu-1.0_dp)*log(2.0_dp*acos(-1.0_dp))-0.5_dp*log(nu)
         mean=a-(nu-1.0_dp)/(2.0_dp*nu)
         var=a/nu
         elogfact=-a*(1.0_dp-llambda/nu)+llambda/(2.0_dp*nu*nu) &
            +0.5_dp*log(2.0_dp*acos(-1.0_dp))+0.5_dp/nu
         return
      end if
      call cmp_logz_trunc(lambda,nu,ctrl%truncate_tol,ctrl%ymax,logz,m,ok)
      llambda=log(lambda);lp=-logz;lf=0.0_dp;mean=0.0_dp;m2=0.0_dp;elogfact=0.0_dp
      do j=1,m
         lp=lp+llambda-nu*log(real(j,dp));lf=lf+log(real(j,dp));prob=exp(lp)
         mean=mean+real(j,dp)*prob;m2=m2+real(j,dp)**2*prob;elogfact=elogfact+lf*prob
      end do
      var=max(0.0_dp,m2-mean*mean)
   end subroutine cmp_stats

   function dcmp(x,lambda,nu,log_p,control) result(out)
      integer, intent(in) :: x
      real(dp), intent(in) :: lambda,nu
      logical, intent(in), optional :: log_p
      type(cmp_control_t), intent(in), optional :: control
      real(dp) :: out, lp
      logical :: lg
      lg=.false.; if (present(log_p)) lg=log_p
      if (x < 0 .or. lambda < 0.0_dp .or. nu < 0.0_dp) then
         lp=-huge(1.0_dp)
      else if (lambda == 0.0_dp) then
         if (x == 0) then
            lp=0.0_dp
         else
            lp=-huge(1.0_dp)
         end if
      else
         lp=real(x,dp)*log(lambda)-nu*log_gamma(real(x+1,dp)) &
            - cmp_logz_hybrid(lambda,nu,control)
      end if
      if (lg) then; out=lp; else; out=exp(lp); end if
   end function dcmp

   function pcmp(x,lambda,nu,control) result(out)
      integer, intent(in) :: x
      real(dp), intent(in) :: lambda,nu
      type(cmp_control_t), intent(in), optional :: control
      real(dp) :: out, logz, lcp, lp, llambda
      type(cmp_control_t) :: ctrl
      integer :: j,m
      logical :: ok
      if (x < 0) then; out=0.0_dp; return; end if
      if (nu==0.0_dp .and. lambda<1.0_dp) then
         if (lambda==0.0_dp) then; out=1.0_dp; else; out=1.0_dp-lambda**real(x+1,dp); end if
         return
      end if
      ctrl=cmp_control_t(); if (present(control)) ctrl=control
      call cmp_logz_trunc(lambda,nu,ctrl%truncate_tol,ctrl%ymax,logz,m,ok)
      lcp=-logz
      if (lambda > 0.0_dp) then; llambda=log(lambda); else; llambda=-huge(1.0_dp); end if
      lp=-logz
      do j=1,min(x,m)
         if (lambda == 0.0_dp) exit
         lp=lp+llambda-nu*log(real(j,dp))
         lcp=logadd(lcp,lp)
      end do
      out=min(1.0_dp,exp(lcp))
      if (x >= m .and. ok) out=1.0_dp
   end function pcmp

   function qcmp(p,lambda,nu,log_p,control) result(q)
      real(dp), intent(in) :: p,lambda,nu
      logical, intent(in), optional :: log_p
      type(cmp_control_t), intent(in), optional :: control
      integer :: q
      real(dp) :: target, logz, lcp, lp, llambda
      type(cmp_control_t) :: ctrl
      integer :: j,m
      logical :: ok,lg
      lg=.false.; if (present(log_p)) lg=log_p
      if (lg) then; target=p; else
         if (p <= 0.0_dp) then; q=0; return; end if
         if (p >= 1.0_dp) then; q=huge(1); return; end if
         target=log(p)
      end if
      if (nu==0.0_dp .and. lambda<1.0_dp) then
         if (lambda==0.0_dp) then; q=0; return; end if
         q=max(0,ceiling(log(1.0_dp-exp(target))/log(lambda)-1.0_dp)); return
      end if
      ctrl=cmp_control_t(); if (present(control)) ctrl=control
      call cmp_logz_trunc(lambda,nu,ctrl%truncate_tol,ctrl%ymax,logz,m,ok)
      lcp=-logz
      if (target <= lcp) then; q=0; return; end if
      llambda=log(lambda);lp=-logz
      do j=1,m
         lp=lp+llambda-nu*log(real(j,dp))
         lcp=logadd(lcp,lp)
         if (target <= lcp) then; q=j; return; end if
      end do
      q=m
   end function qcmp

   subroutine rcmp(n,lambda,nu,out,control)
      integer, intent(in) :: n
      real(dp), intent(in) :: lambda,nu
      integer, intent(out) :: out(n)
      type(cmp_control_t), intent(in), optional :: control
      real(dp) :: u
      integer :: i
      do i=1,n
         call random_number(u)
         out(i)=qcmp(max(u,tiny(1.0_dp)),lambda,nu,control=control)
      end do
   end subroutine rcmp

   subroutine rcmp_vec(lambda,nu,out,control)
      real(dp), intent(in) :: lambda(:),nu(:)
      integer, intent(out) :: out(:)
      type(cmp_control_t), intent(in), optional :: control
      real(dp) :: u
      integer :: i
      if (size(lambda)/=size(nu) .or. size(out)/=size(lambda)) error stop 'rcmp_vec: size mismatch'
      do i=1,size(out)
         call random_number(u)
         out(i)=qcmp(max(u,tiny(1.0_dp)),lambda(i),nu(i),control=control)
      end do
   end subroutine rcmp_vec

   function ecmp(lambda,nu,control) result(mu)
      real(dp), intent(in) :: lambda,nu
      type(cmp_control_t), intent(in), optional :: control
      real(dp) :: mu, logz, lp, lsum, llambda
      type(cmp_control_t) :: ctrl
      integer :: j,m
      logical :: ok
      ctrl=cmp_control_t(); if (present(control)) ctrl=control
      if (lambda==0.0_dp) then; mu=0.0_dp; return; end if
      if (nu==0.0_dp .and. lambda<1.0_dp) then; mu=lambda/(1.0_dp-lambda); return; end if
      if (nu==1.0_dp) then; mu=lambda; return; end if
      if (nu>0.0_dp .and. -log(lambda)/nu < log(ctrl%hybrid_tol)) then
         mu=exp(log(lambda)/nu)-(nu-1.0_dp)/(2.0_dp*nu)
         return
      end if
      call cmp_logz_trunc(lambda,nu,ctrl%truncate_tol,ctrl%ymax,logz,m,ok)
      lsum=-huge(1.0_dp);llambda=log(lambda);lp=-logz
      do j=1,m
         lp=lp+llambda-nu*log(real(j,dp))
         lsum=logadd(lsum,log(real(j,dp))+lp)
      end do
      mu=exp(lsum)
   end function ecmp

   function vcmp(lambda,nu,control) result(vr)
      real(dp), intent(in) :: lambda,nu
      type(cmp_control_t), intent(in), optional :: control
      real(dp) :: vr,mu,m2,logz,lp,lsum,llambda
      type(cmp_control_t) :: ctrl
      integer :: j,m
      logical :: ok
      ctrl=cmp_control_t(); if (present(control)) ctrl=control
      if (lambda==0.0_dp) then; vr=0.0_dp; return; end if
      if (nu==0.0_dp .and. lambda<1.0_dp) then; vr=lambda/(1.0_dp-lambda)**2; return; end if
      if (nu==1.0_dp) then; vr=lambda; return; end if
      if (nu>0.0_dp .and. -log(lambda)/nu < log(ctrl%hybrid_tol)) then
         vr=exp(log(lambda)/nu)/nu
         return
      end if
      mu=ecmp(lambda,nu,ctrl)
      call cmp_logz_trunc(lambda,nu,ctrl%truncate_tol,ctrl%ymax,logz,m,ok)
      lsum=-huge(1.0_dp);llambda=log(lambda);lp=-logz
      do j=1,m
         lp=lp+llambda-nu*log(real(j,dp))
         lsum=logadd(lsum,2.0_dp*log(real(j,dp))+lp)
      end do
      m2=exp(lsum); vr=max(0.0_dp,m2-mu*mu)
   end function vcmp

   function ncmp(lambda,nu,log_p,control) result(z)
      real(dp), intent(in) :: lambda,nu
      logical, intent(in), optional :: log_p
      type(cmp_control_t), intent(in), optional :: control
      real(dp) :: z,lz
      logical :: lg
      lg=.false.; if (present(log_p)) lg=log_p
      lz=cmp_logz_hybrid(lambda,nu,control)
      if (lg) then; z=lz; else; z=exp(lz); end if
   end function ncmp

   function tcmp(lambda,nu,control) result(t)
      real(dp), intent(in) :: lambda,nu
      type(cmp_control_t), intent(in), optional :: control
      integer :: t
      t=cmp_y_trunc(lambda,nu,control)
   end function tcmp

   function dzicmp(x,lambda,nu,p,log_p,control) result(out)
      integer, intent(in) :: x
      real(dp), intent(in) :: lambda,nu,p
      logical, intent(in), optional :: log_p
      type(cmp_control_t), intent(in), optional :: control
      real(dp) :: out,lp,la,lb
      logical :: lg
      lg=.false.; if (present(log_p)) lg=log_p
      if (x<0 .or. p<0.0_dp .or. p>1.0_dp) then
         lp=-huge(1.0_dp)
      else if (p==1.0_dp) then
         if (x==0) then; lp=0.0_dp; else; lp=-huge(1.0_dp); end if
      else
         la=log(1.0_dp-p)+dcmp(x,lambda,nu,.true.,control)
         if (x==0 .and. p>0.0_dp) then; lb=log(p); else; lb=-huge(1.0_dp); end if
         lp=logadd(la,lb)
      end if
      if (lg) then; out=lp; else; out=exp(lp); end if
   end function dzicmp

   function pzicmp(x,lambda,nu,p,control) result(out)
      integer, intent(in) :: x
      real(dp), intent(in) :: lambda,nu,p
      type(cmp_control_t), intent(in), optional :: control
      real(dp) :: out
      if (x<0) then; out=0.0_dp; else; out=p+(1.0_dp-p)*pcmp(x,lambda,nu,control); end if
   end function pzicmp

   function qzicmp(prob,lambda,nu,p,log_p,control) result(q)
      real(dp), intent(in) :: prob,lambda,nu,p
      logical, intent(in), optional :: log_p
      type(cmp_control_t), intent(in), optional :: control
      integer :: q
      real(dp) :: pp,adj
      logical :: lg
      lg=.false.; if (present(log_p)) lg=log_p
      if (lg) then; pp=exp(prob); else; pp=prob; end if
      if (pp <= p+(1.0_dp-p)*dcmp(0,lambda,nu,control=control)) then
         q=0
      else
         adj=(pp-p)/(1.0_dp-p)
         q=qcmp(adj,lambda,nu,control=control)
      end if
   end function qzicmp

   subroutine rzicmp(n,lambda,nu,p,out,control)
      integer, intent(in) :: n
      real(dp), intent(in) :: lambda,nu,p
      integer, intent(out) :: out(n)
      type(cmp_control_t), intent(in), optional :: control
      real(dp) :: u
      integer :: i
      do i=1,n
         call random_number(u)
         if (u < p) then; out(i)=0; else
            call random_number(u); out(i)=qcmp(max(u,tiny(1.0_dp)),lambda,nu,control=control)
         end if
      end do
   end subroutine rzicmp

   subroutine rzicmp_vec(lambda,nu,p,out,control)
      real(dp), intent(in) :: lambda(:),nu(:),p(:)
      integer, intent(out) :: out(:)
      type(cmp_control_t), intent(in), optional :: control
      real(dp) :: u
      integer :: i
      if (size(out)/=size(lambda) .or. size(nu)/=size(out) .or. size(p)/=size(out)) &
         error stop 'rzicmp_vec: size mismatch'
      do i=1,size(out)
         call random_number(u)
         if (u<p(i)) then; out(i)=0; else
            call random_number(u)
            out(i)=qcmp(max(u,tiny(1.0_dp)),lambda(i),nu(i),control=control)
         end if
      end do
   end subroutine rzicmp_vec

   real(dp) function ezicmp(lambda,nu,p,control)
      real(dp), intent(in) :: lambda,nu,p
      type(cmp_control_t), intent(in), optional :: control
      ezicmp=(1.0_dp-p)*ecmp(lambda,nu,control)
   end function ezicmp

   real(dp) function vzicmp(lambda,nu,p,control)
      real(dp), intent(in) :: lambda,nu,p
      type(cmp_control_t), intent(in), optional :: control
      real(dp) :: ee,vv
      ee=ecmp(lambda,nu,control); vv=vcmp(lambda,nu,control)
      vzicmp=(1.0_dp-p)*(vv+p*ee*ee)
   end function vzicmp

   pure real(dp) function poisson_logpmf(x,lambda)
      integer, intent(in) :: x
      real(dp), intent(in) :: lambda
      if (x<0 .or. lambda<0.0_dp) then; poisson_logpmf=-huge(1.0_dp)
      else if (lambda==0.0_dp) then
         if (x==0) then; poisson_logpmf=0.0_dp; else; poisson_logpmf=-huge(1.0_dp); end if
      else
         poisson_logpmf=-lambda+real(x,dp)*log(lambda)-log_gamma(real(x+1,dp))
      end if
   end function poisson_logpmf

   function dzip(x,lambda,p,log_p) result(out)
      integer,intent(in)::x; real(dp),intent(in)::lambda,p
      logical,intent(in),optional::log_p
      real(dp)::out,lp,la,lb; logical::lg
      lg=.false.; if(present(log_p)) lg=log_p
      if(p==1.0_dp)then
         if(x==0)then;lp=0.0_dp;else;lp=-huge(1.0_dp);end if
      else
         la=log(1.0_dp-p)+poisson_logpmf(x,lambda)
         if(x==0.and.p>0.0_dp)then;lb=log(p);else;lb=-huge(1.0_dp);end if
         lp=logadd(la,lb)
      end if
      if(lg)then;out=lp;else;out=exp(lp);end if
   end function dzip

   function pzip(x,lambda,p) result(out)
      integer,intent(in)::x; real(dp),intent(in)::lambda,p; real(dp)::out,term,s
      integer::j
      if(x<0)then;out=0.0_dp;return;end if
      if(lambda==0.0_dp)then;out=1.0_dp;return;end if
      term=exp(-lambda);s=term
      do j=1,x;term=term*lambda/real(j,dp);s=s+term;end do
      out=min(1.0_dp,p+(1.0_dp-p)*s)
   end function pzip

   function qzip(prob,lambda,p,log_p) result(q)
      real(dp),intent(in)::prob,lambda,p;logical,intent(in),optional::log_p
      integer::q;real(dp)::pp,adj,term,s;logical::lg
      lg=.false.;if(present(log_p))lg=log_p
      if(lg)then;pp=exp(prob);else;pp=prob;end if
      if(pp<=p+(1.0_dp-p)*exp(-lambda))then;q=0;return;end if
      adj=(pp-p)/(1.0_dp-p);term=exp(-lambda);s=term;q=0
      do while(s<adj.and.q<10000000)
         q=q+1;term=term*lambda/real(q,dp);s=s+term
      end do
   end function qzip

   subroutine rzip(n,lambda,p,out)
      integer,intent(in)::n;real(dp),intent(in)::lambda,p;integer,intent(out)::out(n)
      real(dp)::u;integer::i
      do i=1,n;call random_number(u);out(i)=qzip(u,lambda,p);end do
   end subroutine rzip

   pure real(dp) function ezip(lambda,p)
      real(dp),intent(in)::lambda,p;ezip=(1.0_dp-p)*lambda
   end function ezip

   pure real(dp) function vzip(lambda,p)
      real(dp),intent(in)::lambda,p;vzip=(1.0_dp-p)*(lambda+p*lambda*lambda)
   end function vzip

   function loglik_cmp(y,lambda,nu,control) result(ll)
      integer,intent(in)::y(:);real(dp),intent(in)::lambda(:),nu(:)
      type(cmp_control_t),intent(in),optional::control
      real(dp)::ll;integer::i
      if(size(y)/=size(lambda).or.size(y)/=size(nu)) error stop 'loglik_cmp: size mismatch'
      ll=0.0_dp;do i=1,size(y);ll=ll+dcmp(y(i),lambda(i),nu(i),.true.,control);end do
   end function loglik_cmp

   function loglik_zicmp(y,lambda,nu,p,control) result(ll)
      integer,intent(in)::y(:);real(dp),intent(in)::lambda(:),nu(:),p(:)
      type(cmp_control_t),intent(in),optional::control
      real(dp)::ll;integer::i
      if(size(y)/=size(lambda).or.size(y)/=size(nu).or.size(y)/=size(p)) &
         error stop 'loglik_zicmp: size mismatch'
      ll=0.0_dp;do i=1,size(y);ll=ll+dzicmp(y(i),lambda(i),nu(i),p(i),.true.,control);end do
   end function loglik_zicmp

end module compoissonreg_distributions
