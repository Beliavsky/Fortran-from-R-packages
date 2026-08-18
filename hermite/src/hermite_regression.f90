module hermite_regression
   use hermite_kinds, only : dp, i64
   use hermite_numerics, only : invert_matrix, chi_square1_upper
   use hermite_distribution, only : hermite_logpmf_exact
   implicit none
   private

   integer, parameter, public :: HERMITE_LINK_LOG=1
   integer, parameter, public :: HERMITE_LINK_IDENTITY=2

   type, public :: hermite_glm_result
      real(dp), allocatable :: beta(:)
      real(dp) :: dispersion = 1.0_dp
      integer :: order = 1
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: lr_stat = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: aic = huge(1.0_dp)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: hessian(:,:)
      real(dp), allocatable :: vcov(:,:)
      integer :: status = 0
      integer :: iterations = 0
   end type hermite_glm_result

   public :: hermite_prob_mu_d, hermite_glm_loglik
   public :: fit_glm_hermite

contains

   pure elemental real(dp) function logistic(z) result(p)
      real(dp), intent(in) :: z
      if (z >= 0.0_dp) then
         p=1.0_dp/(1.0_dp+exp(-z))
      else
         p=exp(z)/(1.0_dp+exp(z))
      end if
   end function logistic

   real(dp) function hermite_prob_mu_d(y,mu,d,m,log_p) result(p)
      integer(i64), intent(in) :: y
      real(dp), intent(in) :: mu,d
      integer, intent(in) :: m
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: a,b,lpr

      lp=.false.
      if (present(log_p)) lp=log_p
      if (m == 1) then
         if (mu <= 0.0_dp .or. y < 0_i64) then
            lpr=-huge(1.0_dp)
         else if (mu <= tiny(1.0_dp)) then
            if (y == 0_i64) then
               lpr=0.0_dp
            else
               lpr=-huge(1.0_dp)
            end if
         else
            lpr=-mu+real(y,dp)*log(mu)-log_gamma(real(y+1_i64,dp))
         end if
      else
         if (mu <= 0.0_dp .or. d < 1.0_dp .or. d > real(m,dp)) then
            lpr=-huge(1.0_dp)
         else
            b=mu*(d-1.0_dp)/real(m*(m-1),dp)
            a=mu-real(m,dp)*b
            lpr=hermite_logpmf_exact(y,a,b,m)
         end if
      end if
      if (lp) then
         p=lpr
      else if (lpr <= -0.5_dp*huge(1.0_dp)) then
         p=0.0_dp
      else
         p=exp(lpr)
      end if
   end function hermite_prob_mu_d

   real(dp) function hermite_glm_loglik(y,x,beta,d,m,link) result(ll)
      integer(i64), intent(in) :: y(:)
      real(dp), intent(in) :: x(:,:),beta(:),d
      integer, intent(in) :: m,link
      real(dp) :: eta,mu,lp
      integer :: i

      if (size(x,1) /= size(y) .or. size(x,2) /= size(beta)) then
         ll=-huge(1.0_dp)
         return
      end if
      ll=0.0_dp
      do i=1,size(y)
         eta=dot_product(x(i,:),beta)
         if (link==HERMITE_LINK_LOG) then
            if (eta > log(huge(1.0_dp))-2.0_dp) then
               ll=-huge(1.0_dp)
         return
            end if
            mu=exp(eta)
         else if (link==HERMITE_LINK_IDENTITY) then
            mu=eta
         else
            ll=-huge(1.0_dp)
         return
         end if
         if (mu <= 0.0_dp) then
            ll=-huge(1.0_dp)
         return
         end if
         lp=hermite_prob_mu_d(y(i),mu,d,m,.true.)
         if (lp <= -0.5_dp*huge(1.0_dp)) then
            ll=-huge(1.0_dp)
         return
         end if
         ll=ll+lp
      end do
   end function hermite_glm_loglik

   subroutine initial_beta(y,x,link,beta)
      integer(i64), intent(in) :: y(:)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: link
      real(dp), intent(out) :: beta(size(x,2))
      real(dp) :: mean_y

      beta=0.0_dp
      mean_y=max(sum(real(y,dp))/real(size(y),dp),1.0e-3_dp)
      if (all(abs(x(:,1)-1.0_dp) <= 100.0_dp*epsilon(1.0_dp))) then
         if (link==HERMITE_LINK_LOG) then
            beta(1)=log(mean_y)
         else
            beta(1)=mean_y
         end if
      else if (link==HERMITE_LINK_IDENTITY) then
         beta=0.1_dp
      end if
   end subroutine initial_beta

   real(dp) function nll_trans(y,x,par,m,link) result(v)
      integer(i64), intent(in) :: y(:)
      real(dp), intent(in) :: x(:,:),par(:)
      integer, intent(in) :: m,link
      integer :: p
      real(dp) :: d,ll
      p=size(x,2)
      if (m==1) then
         d=1.0_dp
      else
         d=1.0_dp+(real(m,dp)-1.0_dp)*logistic(par(p+1))
      end if
      ll=hermite_glm_loglik(y,x,par(1:p),d,m,link)
      if (ll <= -0.5_dp*huge(1.0_dp)) then
         v=huge(1.0_dp)/100.0_dp
      else
         v=-ll
      end if
   end function nll_trans

   real(dp) function nll_natural(y,x,par,m,link) result(v)
      integer(i64), intent(in) :: y(:)
      real(dp), intent(in) :: x(:,:),par(:)
      integer, intent(in) :: m,link
      integer :: p
      real(dp) :: d,ll
      p=size(x,2)
      if (m==1) then
         d=1.0_dp
      else
         d=par(p+1)
         if (d < 1.0_dp .or. d > real(m,dp)) then
            v=huge(1.0_dp)/100.0_dp
            return
         end if
      end if
      ll=hermite_glm_loglik(y,x,par(1:p),d,m,link)
      if (ll <= -0.5_dp*huge(1.0_dp)) then
         v=huge(1.0_dp)/100.0_dp
      else
         v=-ll
      end if
   end function nll_natural

   subroutine order_simplex(f,idx)
      real(dp), intent(in) :: f(:)
      integer, intent(out) :: idx(size(f))
      integer :: i,j,t
      idx=[(i,i=1,size(f))]
      do i=1,size(f)-1
         do j=i+1,size(f)
            if (f(idx(j)) < f(idx(i))) then
               t=idx(i)
               idx(i)=idx(j)
               idx(j)=t
            end if
         end do
      end do
   end subroutine order_simplex

   subroutine optimize_order(y,x,m,link,start,opt,fbest,iterations,status)
      integer(i64), intent(in) :: y(:)
      real(dp), intent(in) :: x(:,:),start(:)
      integer, intent(in) :: m,link
      real(dp), allocatable, intent(out) :: opt(:)
      real(dp), intent(out) :: fbest
      integer, intent(out) :: iterations,status
      real(dp), allocatable :: simp(:,:),f(:),cent(:),xr(:),xe(:),xc(:)
      integer, allocatable :: idx(:)
      real(dp) :: fr,fe,fc,spread,step
      integer :: n,j,best,second,worst,imax

      n=size(start)
      imax=max(1800,400*n)
      step=0.10_dp
      allocate(simp(n,n+1),f(n+1),cent(n),xr(n),xe(n),xc(n),idx(n+1))
      simp(:,1)=start
      do j=2,n+1
         simp(:,j)=start
         simp(j-1,j)=simp(j-1,j)+step*max(1.0_dp,abs(start(j-1)))
      end do
      do j=1,n+1
         f(j)=nll_trans(y,x,simp(:,j),m,link)
      end do
      status=1
      do iterations=1,imax
         call order_simplex(f,idx)
         best=idx(1)
         second=idx(n)
         worst=idx(n+1)
         cent=0.0_dp
         do j=1,n
            cent=cent+simp(:,idx(j))
         end do
         cent=cent/real(n,dp)
         xr=cent+(cent-simp(:,worst))
         fr=nll_trans(y,x,xr,m,link)
         if (fr < f(best)) then
            xe=cent+2.0_dp*(xr-cent)
            fe=nll_trans(y,x,xe,m,link)
            if (fe < fr) then
               simp(:,worst)=xe
               f(worst)=fe
            else
               simp(:,worst)=xr
               f(worst)=fr
            end if
         else if (fr < f(second)) then
            simp(:,worst)=xr
            f(worst)=fr
         else
            if (fr < f(worst)) then
               xc=cent+0.5_dp*(xr-cent)
            else
               xc=cent+0.5_dp*(simp(:,worst)-cent)
            end if
            fc=nll_trans(y,x,xc,m,link)
            if (fc < min(fr,f(worst))) then
               simp(:,worst)=xc
               f(worst)=fc
            else
               do j=2,n+1
                  simp(:,idx(j))=simp(:,best)+0.5_dp*(simp(:,idx(j))-simp(:,best))
                  f(idx(j))=nll_trans(y,x,simp(:,idx(j)),m,link)
               end do
            end if
         end if
         spread=0.0_dp
         do j=2,n+1
            spread=max(spread,maxval(abs(simp(:,idx(j))-simp(:,best))))
         end do
         if (spread <= 2.0e-8_dp*(1.0_dp+maxval(abs(simp(:,best))))) then
            status=0
            exit
         end if
      end do
      call order_simplex(f,idx)
      allocate(opt(n))
      opt=simp(:,idx(1))
      fbest=f(idx(1))
   end subroutine optimize_order

   subroutine natural_hessian(y,x,par,m,link,hess,status)
      integer(i64), intent(in) :: y(:)
      real(dp), intent(in) :: x(:,:),par(:)
      integer, intent(in) :: m,link
      real(dp), allocatable, intent(out) :: hess(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: xp(:),xm(:),xpp(:),xpm(:),xmp(:),xmm(:)
      real(dp) :: f0,hi,hj
      integer :: n,i,j

      n=size(par)
      allocate(hess(n,n),xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n))
      f0=nll_natural(y,x,par,m,link)
      if (f0 >= 0.1_dp*huge(1.0_dp)) then
         hess=huge(1.0_dp)
         status=1
         return
      end if
      do i=1,n
         hi=3.0e-4_dp*max(1.0_dp,abs(par(i)))
         if (m>1 .and. i==n) then
            hi=min(hi,0.2_dp*max(1.0e-6_dp,par(i)-1.0_dp), &
                        0.2_dp*max(1.0e-6_dp,real(m,dp)-par(i)))
         end if
         xp=par
         xm=par
         xp(i)=xp(i)+hi
         xm(i)=xm(i)-hi
         hess(i,i)=(nll_natural(y,x,xp,m,link)-2.0_dp*f0+ &
                    nll_natural(y,x,xm,m,link))/(hi*hi)
         do j=i+1,n
            hj=3.0e-4_dp*max(1.0_dp,abs(par(j)))
            if (m>1 .and. j==n) then
               hj=min(hj,0.2_dp*max(1.0e-6_dp,par(j)-1.0_dp), &
                           0.2_dp*max(1.0e-6_dp,real(m,dp)-par(j)))
            end if
            xpp=par
            xpm=par
            xmp=par
            xmm=par
            xpp(i)=xpp(i)+hi
            xpp(j)=xpp(j)+hj
            xpm(i)=xpm(i)+hi
            xpm(j)=xpm(j)-hj
            xmp(i)=xmp(i)-hi
            xmp(j)=xmp(j)+hj
            xmm(i)=xmm(i)-hi
            xmm(j)=xmm(j)-hj
            hess(i,j)=(nll_natural(y,x,xpp,m,link)-nll_natural(y,x,xpm,m,link)- &
                       nll_natural(y,x,xmp,m,link)+nll_natural(y,x,xmm,m,link))/ &
                       (4.0_dp*hi*hj)
            hess(j,i)=hess(i,j)
         end do
      end do
      status=0
   end subroutine natural_hessian

   function fit_one_order(y,x,m,link,start_beta,start_d) result(res)
      integer(i64), intent(in) :: y(:)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: m,link
      real(dp), intent(in), optional :: start_beta(:),start_d
      type(hermite_glm_result) :: res
      real(dp), allocatable :: start(:),opt(:),beta0(:),nat(:),invh(:,:)
      real(dp) :: fbest,d,delta
      integer :: p,npar,hstat,istat,i

      p=size(x,2)
      npar=p+merge(0,1,m==1)
      allocate(start(npar),beta0(p))
      call initial_beta(y,x,link,beta0)
      if (present(start_beta)) then
         if (size(start_beta)==p) beta0=start_beta
      end if
      start(1:p)=beta0
      if (m>1) then
         d=1.1_dp
         if (present(start_d)) d=start_d
         d=min(real(m,dp)-1.0e-6_dp,max(1.0_dp+1.0e-6_dp,d))
         delta=log((d-1.0_dp)/(real(m,dp)-d))
         start(p+1)=delta
      end if

      call optimize_order(y,x,m,link,start,opt,fbest,res%iterations,res%status)
      allocate(res%beta(p))
      res%beta=opt(1:p)
      res%order=m
      if (m==1) then
         res%dispersion=1.0_dp
      else
         res%dispersion=1.0_dp+(real(m,dp)-1.0_dp)*logistic(opt(p+1))
      end if
      res%loglik=-fbest
      allocate(res%fitted(size(y)))
      do i=1,size(y)
         if (link==HERMITE_LINK_LOG) then
            res%fitted(i)=exp(dot_product(x(i,:),res%beta))
         else
            res%fitted(i)=dot_product(x(i,:),res%beta)
         end if
      end do

      allocate(nat(npar))
      nat(1:p)=res%beta
      if (m>1) nat(p+1)=res%dispersion
      call natural_hessian(y,x,nat,m,link,res%hessian,hstat)
      if (hstat==0) then
         call invert_matrix(res%hessian,invh,istat)
         if (istat==0) res%vcov=invh
      end if
   end function fit_one_order

   function fit_glm_hermite(y,x,link,m,start_beta,start_d,max_order) result(best)
      integer(i64), intent(in) :: y(:)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in), optional :: link,m,max_order
      real(dp), intent(in), optional :: start_beta(:),start_d
      type(hermite_glm_result) :: best
      type(hermite_glm_result) :: cur,pois
      integer :: lnk,mm,maxm,j,kpars
      real(dp) :: w

      lnk=HERMITE_LINK_LOG
      if (present(link)) lnk=link
      if (size(x,1) /= size(y) .or. size(y)==0 .or. any(y<0_i64)) then
         best%status=10
         return
      end if

      pois=fit_one_order(y,x,1,lnk,start_beta,1.0_dp)
      best=pois
      if (present(m)) then
         mm=m
         if (mm < 1) then
            best%status=11
            return
         end if
         if (mm>1) best=fit_one_order(y,x,mm,lnk,start_beta,start_d)
      else
         maxm=min(10,max(2,int(maxval(y))))
         if (present(max_order)) maxm=min(maxm,max_order)
         do j=2,maxm
            cur=fit_one_order(y,x,j,lnk,start_beta,start_d)
            if (cur%loglik > best%loglik .and. cur%status <= 1) best=cur
         end do
      end if

      if (best%order==1) then
         best%lr_stat=0.0_dp
         best%p_value=1.0_dp
         kpars=size(best%beta)
      else
         w=max(0.0_dp,2.0_dp*(best%loglik-pois%loglik))
         best%lr_stat=w
         best%p_value=0.5_dp*chi_square1_upper(w)
         kpars=size(best%beta)+1
      end if
      best%aic=2.0_dp*real(kpars,dp)-2.0_dp*best%loglik
   end function fit_glm_hermite

end module hermite_regression
