module discrete_laplace
   use, intrinsic :: iso_fortran_env, only : real64, int64
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
   implicit none
   private

   integer, parameter, public :: dp = real64

   type, public :: edlaplace_result
      real(dp) :: e1 = 0.0_dp
      real(dp) :: e1a = 0.0_dp
      real(dp) :: v = 0.0_dp
   end type

   type, public :: edlaplace2_result
      real(dp) :: e1 = 0.0_dp
      real(dp) :: e2 = 0.0_dp
   end type

   type, public :: estdlaplace_result
      real(dp) :: hatp = 0.0_dp
      real(dp) :: hatq = 0.0_dp
      real(dp) :: hatsigma(2,2) = 0.0_dp
      integer :: status = 0
   end type

   public :: ddlaplace, pdlaplace, qdlaplace, rdlaplace
   public :: ddlaplace2, palaplace2, pdlaplace2, qdlaplace2, rdlaplace2
   public :: edlaplace, edlaplace2, ifi, ifi2, iofi2
   public :: estdlaplace, estdlaplace2, dlaplacelike2, loss
   public :: set_rng_seed

contains

   pure elemental real(dp) function log1p_local(x)
      real(dp), intent(in) :: x
      real(dp) :: term, sumv
      integer :: k
      if (abs(x) > 1.0e-4_dp) then
         log1p_local=log(1.0_dp+x)
      else
         term=x
         sumv=term
         do k=2,12
            term=term*(-x)*real(k-1,dp)/real(k,dp)
            sumv=sumv+term
         end do
         log1p_local=sumv
      end if
   end function

   pure logical function valid_pq(p,q,allow_zero)
      real(dp), intent(in) :: p,q
      logical, intent(in), optional :: allow_zero
      logical :: az
      az=.true.
      if (present(allow_zero)) az=allow_zero
      if (az) then
         valid_pq = p >= 0.0_dp .and. p < 1.0_dp .and. &
                    q >= 0.0_dp .and. q < 1.0_dp
      else
         valid_pq = p > 0.0_dp .and. p < 1.0_dp .and. &
                    q > 0.0_dp .and. q < 1.0_dp
      end if
   end function

   pure elemental real(dp) function ddlaplace(x,p,q)
      integer, intent(in) :: x
      real(dp), intent(in) :: p,q
      real(dp) :: c
      if (.not. valid_pq(p,q)) then
         ddlaplace=ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      c=(1.0_dp-p)*(1.0_dp-q)/(1.0_dp-p*q)
      if (x >= 0) then
         ddlaplace=c*p**x
      else
         ddlaplace=c*q**abs(x)
      end if
   end function

   pure elemental real(dp) function pdlaplace(x,p,q)
      real(dp), intent(in) :: x,p,q
      integer :: k
      if (.not. valid_pq(p,q)) then
         pdlaplace=ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      k=floor(x)
      if (x >= 0.0_dp) then
         pdlaplace=1.0_dp-(1.0_dp-q)*p**(k+1)/(1.0_dp-p*q)
      else
         pdlaplace=(1.0_dp-p)*q**(-k)/(1.0_dp-p*q)
      end if
      pdlaplace=min(1.0_dp,max(0.0_dp,pdlaplace))
   end function

   integer(int64) function qdlaplace(prob,p,q)
      real(dp), intent(in) :: prob,p,q
      integer(int64) :: k
      real(dp) :: pk
      if (.not. valid_pq(p,q) .or. prob < 0.0_dp .or. prob > 1.0_dp) then
         qdlaplace=-huge(0_int64)
         return
      end if
      if (prob <= 0.0_dp) then
         qdlaplace=-huge(0_int64)
         return
      end if
      if (prob >= 1.0_dp) then
         qdlaplace=huge(0_int64)
         return
      end if
      k=0_int64
      pk=pdlaplace(0.0_dp,p,q)
      if (prob >= pk) then
         do while (prob >= pk)
            k=k+1_int64
            pk=pdlaplace(real(k,dp),p,q)
            if (k > 100000000_int64) exit
         end do
      else
         do while (prob < pk)
            k=k-1_int64
            pk=pdlaplace(real(k,dp),p,q)
            if (k < -100000000_int64) exit
         end do
         k=k+1_int64
      end if
      qdlaplace=k
   end function

   subroutine rdlaplace(x,p,q)
      integer(int64), intent(out) :: x(:)
      real(dp), intent(in) :: p,q
      real(dp) :: u
      integer :: i
      do i=1,size(x)
         call random_number(u)
         x(i)=qdlaplace(u,p,q)
      end do
   end subroutine

   pure elemental real(dp) function palaplace2(x,p,q)
      real(dp), intent(in) :: x,p,q
      if (.not. valid_pq(p,q,.false.)) then
         palaplace2=ieee_value(0.0_dp,ieee_quiet_nan)
      else if (x > 0.0_dp) then
         palaplace2=1.0_dp-log(q)/log(p*q)*p**x
      else if (x < 0.0_dp) then
         palaplace2=log(p)/log(p*q)*q**(-x)
      else
         palaplace2=log(p)/log(p*q)
      end if
   end function

   pure elemental real(dp) function ddlaplace2(x,p,q)
      integer, intent(in) :: x
      real(dp), intent(in) :: p,q
      ddlaplace2=palaplace2(real(x+1,dp),p,q)-palaplace2(real(x,dp),p,q)
   end function

   pure elemental real(dp) function pdlaplace2(x,p,q)
      real(dp), intent(in) :: x,p,q
      pdlaplace2=palaplace2(real(floor(x)+1,dp),p,q)
   end function

   integer(int64) function qdlaplace2(prob,p,q)
      real(dp), intent(in) :: prob,p,q
      real(dp) :: h,v
      if (.not. valid_pq(p,q,.false.) .or. prob < 0.0_dp .or. prob > 1.0_dp) then
         qdlaplace2=-huge(0_int64)
         return
      end if
      if (prob <= 0.0_dp) then
         qdlaplace2=-huge(0_int64)
         return
      end if
      if (prob >= 1.0_dp) then
         qdlaplace2=huge(0_int64)
         return
      end if
      h=1.0_dp-log(q)/log(p*q)*p
      if (prob >= h) then
         v=(log1p_local(-prob)+log(log(p*q)/log(q)))/log(p)
         qdlaplace2=int(ceiling(v),int64)-1_int64
      else
         v=(log(log(p)/log(p*q))-log(prob))/log(q)
         qdlaplace2=int(ceiling(v),int64)-1_int64
      end if
   end function

   subroutine rdlaplace2(x,p,q)
      integer(int64), intent(out) :: x(:)
      real(dp), intent(in) :: p,q
      real(dp) :: u
      integer :: i
      do i=1,size(x)
         call random_number(u)
         x(i)=qdlaplace2(u,p,q)
      end do
   end subroutine

   pure function edlaplace(p,q) result(r)
      real(dp), intent(in) :: p,q
      type(edlaplace_result) :: r
      r%e1=1.0_dp/(1.0_dp-p)-1.0_dp/(1.0_dp-q)
      r%e1a=(q*(1.0_dp-p)**2+p*(1.0_dp-q)**2) / &
            ((1.0_dp-q*p)*(1.0_dp-q)*(1.0_dp-p))
      r%v=((q*(1.0_dp-p)**3*(1.0_dp+q)+p*(1.0_dp-q)**3*(1.0_dp+p)) / &
           (1.0_dp-p*q)-(p-q)**2) / ((1.0_dp-p)**2*(1.0_dp-q)**2)
   end function

   pure function edlaplace2(p,q) result(r)
      real(dp), intent(in) :: p,q
      type(edlaplace2_result) :: r
      r%e1=log(q)/log(p*q)*p/(1.0_dp-p) - &
           log(p)/log(p*q)/(1.0_dp-q)
      r%e2=log(q)/log(p*q)*p*(1.0_dp+p)/(1.0_dp-p)**2 + &
           log(p)/log(p*q)*(1.0_dp+q)/(1.0_dp-q)**2
   end function

   pure function ifi(p,q) result(m)
      real(dp), intent(in) :: p,q
      real(dp) :: m(2,2),c
      c=p*q*(1.0_dp-p)*(1.0_dp-q)/(1.0_dp+p*q)
      m=c
      m(1,1)=m(1,1)*(1.0_dp-p)*(1.0_dp-p*q*q)/q/(1.0_dp-q)**2
      m(2,2)=m(2,2)*(1.0_dp-q)*(1.0_dp-q*p*p)/p/(1.0_dp-p)**2
   end function

   pure subroutine invert2(a,ainv,status)
      real(dp), intent(in) :: a(2,2)
      real(dp), intent(out) :: ainv(2,2)
      integer, intent(out) :: status
      real(dp) :: det
      det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      if (abs(det) <= tiny(1.0_dp)*max(1.0_dp,maxval(abs(a)))) then
         ainv=ieee_value(0.0_dp,ieee_quiet_nan)
         status=1
         return
      end if
      ainv(1,1)= a(2,2)/det
      ainv(1,2)=-a(1,2)/det
      ainv(2,1)=-a(2,1)/det
      ainv(2,2)= a(1,1)/det
      status=0
   end subroutine

   function ifi2(p,q,status) result(inv)
      real(dp), intent(in) :: p,q
      integer, intent(out), optional :: status
      real(dp) :: inv(2,2),a(2,2),edp2,edpq,edq2
      integer :: st
      edp2=log(q)/log(p)*(-(1.0_dp-p)**2-p*log(p)*log(p*q)) / &
           ((log(p*q))**2*p**2*(1.0_dp-p)**2)
      edpq=1.0_dp/(p*q*(log(p*q))**2)
      edq2=log(p)/log(q)*(-(1.0_dp-q)**2-q*log(q)*log(p*q)) / &
           ((log(p*q))**2*q**2*(1.0_dp-q)**2)
      a=reshape([-edp2,-edpq,-edpq,-edq2],[2,2])
      call invert2(a,inv,st)
      if (present(status)) status=st
   end function

   pure real(dp) function dlaplacelike2(par,x)
      real(dp), intent(in) :: par(2)
      integer, intent(in) :: x(:)
      integer :: i
      real(dp) :: d
      if (.not. valid_pq(par(1),par(2),.false.)) then
         dlaplacelike2=huge(1.0_dp)/100.0_dp
         return
      end if
      dlaplacelike2=0.0_dp
      do i=1,size(x)
         d=ddlaplace2(x(i),par(1),par(2))
         if (d <= 0.0_dp .or. ieee_is_nan(d)) then
            dlaplacelike2=huge(1.0_dp)/100.0_dp
            return
         end if
         dlaplacelike2=dlaplacelike2-log(d)
      end do
   end function

   pure real(dp) function loss(par,x)
      real(dp), intent(in) :: par(2)
      integer, intent(in) :: x(:)
      type(edlaplace2_result) :: m
      real(dp) :: mean1,mean2
      if (.not. valid_pq(par(1),par(2),.false.)) then
         loss=huge(1.0_dp)/100.0_dp
         return
      end if
      mean1=sum(real(x,dp))/real(size(x),dp)
      mean2=sum(real(x,dp)**2)/real(size(x),dp)
      m=edlaplace2(par(1),par(2))
      loss=(mean1-m%e1)**2+(mean2-m%e2)**2
   end function

   function estdlaplace(x) result(r)
      integer, intent(in) :: x(:)
      type(estdlaplace_result) :: r
      real(dp) :: xplus,xminus,mx,p,q,c
      if (size(x) == 0) then
         r%status=1
         r%hatp=ieee_value(0.0_dp,ieee_quiet_nan)
         r%hatq=r%hatp
         r%hatsigma=r%hatp
         return
      end if
      xplus=sum(real(x,dp),mask=x>0)/real(size(x),dp)
      xminus=sum(-real(x,dp),mask=x<0)/real(size(x),dp)
      mx=sum(real(x,dp))/real(size(x),dp)
      if (mx >= 0.0_dp) then
         q=2.0_dp*xminus*(1.0_dp+mx) / &
           (1.0_dp+2.0_dp*xminus*mx+sqrt(1.0_dp+4.0_dp*xminus*xplus))
         p=(q+mx*(1.0_dp-q))/(1.0_dp+mx*(1.0_dp-q))
      else
         p=2.0_dp*xplus*(1.0_dp-mx) / &
           (1.0_dp-2.0_dp*xplus*mx+sqrt(1.0_dp+4.0_dp*xminus*xplus))
         q=(p-mx*(1.0_dp-p))/(1.0_dp-mx*(1.0_dp-p))
      end if
      r%hatp=p; r%hatq=q
      if (p <= 0.0_dp .or. q <= 0.0_dp .or. p >= 1.0_dp .or. q >= 1.0_dp) then
         r%hatsigma=ieee_value(0.0_dp,ieee_quiet_nan)
         r%status=2
         return
      end if
      c=p*q*(1.0_dp-p)*(1.0_dp-q)/(1.0_dp+p*q)
      r%hatsigma=c
      r%hatsigma(1,1)=r%hatsigma(1,1)*(1.0_dp-p)*(1.0_dp-p*q*q) / &
                      q/(1.0_dp-q)**2
      r%hatsigma(2,2)=r%hatsigma(2,2)*(1.0_dp-q)*(1.0_dp-q*p*p) / &
                      p/(1.0_dp-p)**2
   end function

   pure real(dp) function logistic(z)
      real(dp), intent(in) :: z
      if (z >= 0.0_dp) then
         logistic=1.0_dp/(1.0_dp+exp(-z))
      else
         logistic=exp(z)/(1.0_dp+exp(z))
      end if
   end function

   pure real(dp) function logit(p)
      real(dp), intent(in) :: p
      logit=log(p)-log1p_local(-p)
   end function

   pure real(dp) function fit_objective(z,x,which,err)
      real(dp), intent(in) :: z(2),err
      integer, intent(in) :: x(:),which
      real(dp) :: par(2),lo,hi
      if (which == 1) then
         lo=err; hi=1.0_dp-err
         par=lo+(hi-lo)*[logistic(z(1)),logistic(z(2))]
         fit_objective=loss(par,x)
      else
         par=[logistic(z(1)),logistic(z(2))]
         fit_objective=dlaplacelike2(par,x)
      end if
   end function

   subroutine sort3(f,idx)
      real(dp), intent(in) :: f(3)
      integer, intent(out) :: idx(3)
      integer :: i,j,t
      idx=[1,2,3]
      do i=1,2
         do j=i+1,3
            if (f(idx(j)) < f(idx(i))) then
               t=idx(i); idx(i)=idx(j); idx(j)=t
            end if
         end do
      end do
   end subroutine

   subroutine nelder_mead2(x,z0,which,err,zbest,status)
      integer, intent(in) :: x(:),which
      real(dp), intent(in) :: z0(2),err
      real(dp), intent(out) :: zbest(2)
      integer, intent(out) :: status
      real(dp) :: z(2,3),f(3),cent(2),zr(2),ze(2),zc(2),fr,fe,fc
      integer :: i,it,idx(3),b,m,w
      z(:,1)=z0
      z(:,2)=z0+[0.2_dp,0.0_dp]
      z(:,3)=z0+[0.0_dp,0.2_dp]
      do i=1,3
         f(i)=fit_objective(z(:,i),x,which,err)
      end do
      status=1
      do it=1,2000
         call sort3(f,idx); b=idx(1); m=idx(2); w=idx(3)
         if (max(maxval(abs(z(:,b)-z(:,m))),maxval(abs(z(:,b)-z(:,w)))) < 1.0e-10_dp) then
            status=0; exit
         end if
         cent=0.5_dp*(z(:,b)+z(:,m))
         zr=cent+(cent-z(:,w)); fr=fit_objective(zr,x,which,err)
         if (fr < f(b)) then
            ze=cent+2.0_dp*(zr-cent); fe=fit_objective(ze,x,which,err)
            if (fe < fr) then
               z(:,w)=ze; f(w)=fe
            else
               z(:,w)=zr; f(w)=fr
            end if
         else if (fr < f(m)) then
            z(:,w)=zr; f(w)=fr
         else
            if (fr < f(w)) then
               zc=cent+0.5_dp*(zr-cent)
            else
               zc=cent+0.5_dp*(z(:,w)-cent)
            end if
            fc=fit_objective(zc,x,which,err)
            if (fc < min(fr,f(w))) then
               z(:,w)=zc; f(w)=fc
            else
               do i=1,3
                  if (i /= b) then
                     z(:,i)=z(:,b)+0.5_dp*(z(:,i)-z(:,b))
                     f(i)=fit_objective(z(:,i),x,which,err)
                  end if
               end do
            end if
         end if
      end do
      call sort3(f,idx)
      zbest=z(:,idx(1))
   end subroutine

   function estdlaplace2(x,method,err,parml,status) result(par)
      integer, intent(in) :: x(:)
      character(len=*), intent(in), optional :: method
      real(dp), intent(in), optional :: err,parml(2)
      integer, intent(out), optional :: status
      real(dp) :: par(2),er,start(2),zbest(2),r0,rplus,mp,mn,nanv
      character(len=8) :: meth
      integer :: st
      meth='M'; if (present(method)) meth=adjustl(method)
      er=0.001_dp; if (present(err)) er=err
      nanv=ieee_value(0.0_dp,ieee_quiet_nan)
      st=0
      select case(trim(meth))
      case('M','m')
         start=[exp(-1.0_dp),exp(-1.0_dp)]
         start=(start-er)/(1.0_dp-2.0_dp*er)
         start=[logit(start(1)),logit(start(2))]
         call nelder_mead2(x,start,1,er,zbest,st)
         par=er+(1.0_dp-2.0_dp*er)*[logistic(zbest(1)),logistic(zbest(2))]
      case('ML','ml','Ml','mL')
         if (count(x>=0) == size(x) .or. count(x<0) == size(x)) then
            par=nanv; st=2
         else
            if (present(parml)) then
               start=[logit(parml(1)),logit(parml(2))]
            else
               start=[logit(exp(-1.0_dp)),logit(exp(-1.0_dp))]
            end if
            call nelder_mead2(x,start,2,er,zbest,st)
            par=[logistic(zbest(1)),logistic(zbest(2))]
         end if
      case('P','p')
         r0=real(count(x==0),dp)/real(size(x),dp)
         rplus=real(count(x>=0),dp)/real(size(x),dp)
         if (count(x>=0) == count(x==0) .or. count(x==0) == 0 .or. count(x>=0) == size(x)) then
            par=nanv; st=3
         else
            par(1)=1.0_dp-r0/rplus
            par(2)=par(1)**(rplus/(1.0_dp-rplus))
         end if
      case('MM','mm','Mm','mM')
         if (count(x>=0) > 0) then
            mp=sum(real(x,dp),mask=x>=0)/real(count(x>=0),dp)
            par(1)=mp/(1.0_dp+mp)
         else
            par(1)=nanv; st=4
         end if
         if (count(x<0) > 0) then
            mn=sum(real(x,dp),mask=x<0)/real(count(x<0),dp)
            par(2)=1.0_dp-1.0_dp/abs(mn)
         else
            par(2)=nanv; st=max(st,4)
         end if
      case default
         par=nanv; st=5
      end select
      if (present(status)) status=st
   end function

   function iofi2(x,status) result(inv)
      integer, intent(in) :: x(:)
      integer, intent(out), optional :: status
      real(dp) :: inv(2,2),a(2,2),par(2),p,q,dp2,dpq,dq2
      integer :: st
      par=estdlaplace2(x,'ML',status=st)
      if (st /= 0) then
         inv=ieee_value(0.0_dp,ieee_quiet_nan)
         if (present(status)) status=st
         return
      end if
      p=par(1); q=par(2)
      dp2=real(count(x<0),dp)*(-1.0_dp-log(p))/(p*log(p))**2 + &
          real(size(x),dp)*(log(p*q)+1.0_dp)/(p*log(p*q))**2 - &
          sum(real(x,dp),mask=x>=0)/p**2 - real(count(x>=0),dp)/(1.0_dp-p)**2
      dpq=real(size(x),dp)/(p*q*(log(p*q))**2)
      dq2=real(count(x>=0),dp)*(-1.0_dp-log(q))/(q*log(q))**2 + &
          real(size(x),dp)*(log(p*q)+1.0_dp)/(q*log(p*q))**2 + &
          sum(real(x,dp),mask=x<0)/q**2 + &
          real(count(x<0),dp)*(1.0_dp-2.0_dp*q)/(q*(1.0_dp-q))**2
      a=reshape([-dp2,-dpq,-dpq,-dq2],[2,2])
      call invert2(a,inv,st)
      if (present(status)) status=st
   end function

   subroutine set_rng_seed(seed)
      integer, intent(in) :: seed
      integer :: n,i
      integer, allocatable :: s(:)
      call random_seed(size=n)
      allocate(s(n))
      do i=1,n
         s(i)=mod(seed+104729*i,2147483646)+1
      end do
      call random_seed(put=s)
   end subroutine

end module discrete_laplace
