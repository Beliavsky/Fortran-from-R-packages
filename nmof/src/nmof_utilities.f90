! SPDX-License-Identifier: GPL-3.0-only
module nmof_utilities
   use nmof_kinds, only: dp, pi
   use nmof_math, only: quantile_type7
   use nmof_linalg, only: eigen_symmetric, matrix_rank_subset
   use nmof_types, only: quadrature_rule, drawdown_summary, qtable_result, pbo_result, &
                         nmof_ok, nmof_invalid_input, nmof_numerical_failure
   implicit none
   private
   public :: moving_average, partial_moment, conditional_moment
   public :: drawdown_series, drawdown_info, diversification_ratio
   public :: ns_curve, nss_curve, ns_factors, nss_factors
   public :: change_interval, xw_gauss, bracketing
   public :: repair_matrix, column_subset, qtable_statistics, probability_backtest_overfitting
   public :: marginal_risk_contributions, marginal_risk_contributions_fd, test_ackley, test_griewank, test_rastrigin
   public :: test_rosenbrock, test_schwefel, test_trefethen, test_eggholder
   public :: metric_function, scalar_function_context

   abstract interface
      function metric_function(x, context) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         class(*), intent(in), optional :: context
         real(dp) :: value
      end function metric_function
      function scalar_function_context(x, context) result(value)
         import dp
         real(dp), intent(in) :: x
         class(*), intent(in), optional :: context
         real(dp) :: value
      end function scalar_function_context
   end interface
contains
   function moving_average(y, order, pad, use_pad, status) result(ma)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: order
      real(dp), intent(in), optional :: pad
      logical, intent(in), optional :: use_pad
      integer, intent(out), optional :: status
      real(dp) :: ma(size(y))
      real(dp), allocatable :: cs(:)
      integer :: n,i
      logical :: dopad
      n=size(y); ma=0.0_dp; dopad=.false.; if(present(use_pad)) dopad=use_pad
      if(order<1 .or. order>max(1,n)) then
         if(present(status)) status=nmof_invalid_input
         return
      end if
      allocate(cs(n)); cs(1)=y(1)
      do i=2,n; cs(i)=cs(i-1)+y(i); end do
      do i=1,n
         if(i<order) then
            ma(i)=cs(i)/real(order,dp)
            if(dopad.and.present(pad)) ma(i)=pad
         else if(i==order) then
            ma(i)=cs(i)/real(order,dp)
         else
            ma(i)=(cs(i)-cs(i-order))/real(order,dp)
         end if
      end do
      if(present(status)) status=nmof_ok
   end function moving_average

   pure real(dp) function partial_moment(x, xp, threshold, lower, normalise) result(ans)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: xp, threshold
      logical, intent(in), optional :: lower, normalise
      real(dp) :: p,t,z
      logical :: lo,norm
      integer :: i
      p=2.0_dp; if(present(xp)) p=xp
      t=0.0_dp; if(present(threshold)) t=threshold
      lo=.true.; if(present(lower)) lo=lower
      norm=.false.; if(present(normalise)) norm=normalise
      ans=0.0_dp
      do i=1,size(x)
         if(lo) then; z=max(t-x(i),0.0_dp); else; z=max(x(i)-t,0.0_dp); end if
         ans=ans+z**p
      end do
      if(size(x)>0) ans=ans/real(size(x),dp)
      if(norm.and.ans>=0.0_dp) ans=ans**(1.0_dp/p)
   end function partial_moment

   pure real(dp) function conditional_moment(x, xp, threshold, lower, normalise) result(ans)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: xp, threshold
      logical, intent(in), optional :: lower, normalise
      real(dp) :: p,t,z
      logical :: lo,norm
      integer :: i,nsel
      p=2.0_dp; if(present(xp)) p=xp
      t=0.0_dp; if(present(threshold)) t=threshold
      lo=.true.; if(present(lower)) lo=lower
      norm=.false.; if(present(normalise)) norm=normalise
      ans=0.0_dp; nsel=0
      do i=1,size(x)
         if(lo) then; z=max(t-x(i),0.0_dp); else; z=max(x(i)-t,0.0_dp); end if
         if(z>0.0_dp) then; ans=ans+z**p; nsel=nsel+1; end if
      end do
      if(nsel>0) ans=ans/real(nsel,dp)
      if(norm.and.ans>=0.0_dp) ans=ans**(1.0_dp/p)
   end function conditional_moment

   function drawdown_series(v, relative) result(d)
      real(dp), intent(in) :: v(:)
      logical, intent(in), optional :: relative
      real(dp) :: d(size(v)), high
      logical :: rel
      integer :: i
      rel=.true.; if(present(relative)) rel=relative
      if(size(v)==0) return
      high=v(1)
      do i=1,size(v)
         high=max(high,v(i)); d(i)=high-v(i)
         if(rel) then
            if(abs(high)>tiny(1.0_dp)) d(i)=d(i)/high
         end if
      end do
   end function drawdown_series

   function drawdown_info(v, relative) result(info)
      real(dp), intent(in) :: v(:)
      logical, intent(in), optional :: relative
      type(drawdown_summary) :: info
      real(dp), allocatable :: d(:)
      integer :: trough
      if(size(v)==0) return
      allocate(d(size(v))); d=drawdown_series(v,relative); trough=maxloc(d,dim=1)
      info%maximum=d(trough); info%low=v(trough); info%low_position=trough
      info%high_position=maxloc(v(1:trough),dim=1); info%high=v(info%high_position)
   end function drawdown_info

   pure real(dp) function diversification_ratio(w, covariance) result(ratio)
      real(dp), intent(in) :: w(:), covariance(:, :)
      real(dp) :: denom
      integer :: i
      ratio=0.0_dp
      denom=sqrt(max(0.0_dp,dot_product(w,matmul(covariance,w))))
      if(denom>0.0_dp) then
         do i=1,size(w); ratio=ratio+w(i)*sqrt(max(0.0_dp,covariance(i,i))); end do
         ratio=ratio/denom
      end if
   end function diversification_ratio

   function ns_curve(param, tm, status) result(y)
      real(dp), intent(in) :: param(:), tm(:)
      integer, intent(out), optional :: status
      real(dp) :: y(size(tm)), aux(size(tm)), loading(size(tm))
      if(size(param)/=4 .or. any(tm<=0.0_dp) .or. abs(param(4))<=tiny(1.0_dp)) then
         y=0.0_dp; if(present(status)) status=nmof_invalid_input; return
      end if
      aux=tm/param(4); loading=(1.0_dp-exp(-aux))/aux
      y=param(1)+param(2)*loading+param(3)*(loading-exp(-aux))
      if(present(status)) status=nmof_ok
   end function ns_curve

   function nss_curve(param, tm, status) result(y)
      real(dp), intent(in) :: param(:), tm(:)
      integer, intent(out), optional :: status
      real(dp) :: y(size(tm)),g1(size(tm)),g2(size(tm)),a1(size(tm)),a2(size(tm))
      if(size(param)/=6.or.any(tm<=0.0_dp).or.abs(param(5))<=tiny(1.0_dp).or.abs(param(6))<=tiny(1.0_dp)) then
         y=0.0_dp; if(present(status)) status=nmof_invalid_input; return
      end if
      g1=tm/param(5); g2=tm/param(6); a1=1.0_dp-exp(-g1); a2=1.0_dp-exp(-g2)
      y=param(1)+param(2)*(a1/g1)+param(3)*(a1/g1+a1-1.0_dp)+param(4)*(a2/g2+a2-1.0_dp)
      if(present(status)) status=nmof_ok
   end function nss_curve

   function ns_factors(lambda, tm) result(x)
      real(dp),intent(in)::lambda,tm(:)
      real(dp)::x(size(tm),3),a(size(tm))
      a=tm/lambda; x(:,1)=1.0_dp; x(:,2)=(1.0_dp-exp(-a))/a; x(:,3)=x(:,2)-exp(-a)
   end function ns_factors

   function nss_factors(lambda1,lambda2,tm) result(x)
      real(dp),intent(in)::lambda1,lambda2,tm(:)
      real(dp)::x(size(tm),4),a1(size(tm)),a2(size(tm))
      a1=tm/lambda1; a2=tm/lambda2; x(:,1)=1.0_dp
      x(:,2)=(1.0_dp-exp(-a1))/a1; x(:,3)=x(:,2)-exp(-a1); x(:,4)=(1.0_dp-exp(-a2))/a2-exp(-a2)
   end function nss_factors

   subroutine change_interval(nodes,weights,oldmin,oldmax,newmin,newmax)
      real(dp),intent(inout)::nodes(:),weights(:)
      real(dp),intent(in)::oldmin,oldmax,newmin,newmax
      real(dp)::nr,orr
      nr=newmax-newmin; orr=oldmax-oldmin
      nodes=(nr*nodes+newmin*oldmax-newmax*oldmin)/orr; weights=weights*nr/orr
   end subroutine change_interval

   function xw_gauss(n, method) result(rule)
      integer,intent(in)::n
      character(len=*),intent(in),optional::method
      type(quadrature_rule)::rule
      real(dp),allocatable::a(:,:),values(:),vectors(:,:)
      real(dp)::scale
      integer::i,info
      character(len=16)::m
      if(n<1) then; rule%status=nmof_invalid_input; return; end if
      m='legendre'; if(present(method)) m=lowercase(method)
      allocate(a(n,n),values(n),vectors(n,n)); a=0.0_dp
      select case(trim(m))
      case('legendre')
         do i=1,n-1; a(i,i+1)=1.0_dp/sqrt(4.0_dp-real(i,dp)**(-2)); a(i+1,i)=a(i,i+1); end do
         scale=2.0_dp
      case('laguerre')
         do i=1,n; a(i,i)=2.0_dp*real(i,dp)-1.0_dp; end do
         do i=1,n-1; a(i,i+1)=real(i,dp); a(i+1,i)=a(i,i+1); end do
         scale=1.0_dp
      case('hermite')
         do i=1,n-1; a(i,i+1)=sqrt(real(i,dp)/2.0_dp); a(i+1,i)=a(i,i+1); end do
         scale=sqrt(pi)
      case default
         rule%status=nmof_invalid_input; return
      end select
      call eigen_symmetric(a,values,vectors,info)
      if(info/=0) then; rule%status=nmof_numerical_failure; return; end if
      allocate(rule%nodes(n),rule%weights(n)); rule%nodes=values; rule%weights=scale*vectors(1,:)**2
      rule%status=nmof_ok
   end function xw_gauss

   function bracketing(fun,lower,upper,n,context) result(intervals)
      procedure(scalar_function_context)::fun
      real(dp),intent(in)::lower,upper
      integer,intent(in),optional::n
      class(*),intent(in),optional::context
      real(dp),allocatable::intervals(:,:)
      real(dp),allocatable::xs(:),fx(:)
      integer::nn,i,k,countx
      nn=20; if(present(n)) nn=n
      if(nn<2.or.lower>=upper) then; allocate(intervals(0,2)); return; end if
      allocate(xs(nn),fx(nn)); do i=1,nn; xs(i)=lower+(upper-lower)*real(i-1,dp)/real(nn-1,dp)
         if(present(context)) then; fx(i)=fun(xs(i),context); else; fx(i)=fun(xs(i)); end if
      end do
      countx=0; do i=1,nn-1; if(opposite_sign(fx(i),fx(i+1))) countx=countx+1; end do
      allocate(intervals(countx,2)); k=0
      do i=1,nn-1
         if(opposite_sign(fx(i),fx(i+1))) then; k=k+1; intervals(k,:)=[xs(i),xs(i+1)]; end if
      end do
   end function bracketing

   function repair_matrix(c,eps,status) result(bb)
      real(dp),intent(in)::c(:,:)
      real(dp),intent(in),optional::eps
      integer,intent(out),optional::status
      real(dp)::bb(size(c,1),size(c,2)),floor_value
      real(dp),allocatable::values(:),vectors(:,:),scaled(:)
      integer::info,i,j,n
      n=size(c,1); floor_value=0.0_dp; if(present(eps)) floor_value=eps
      if(size(c,2)/=n .or. n<1 .or. floor_value<0.0_dp) then
         bb=0.0_dp; if(present(status)) status=nmof_invalid_input; return
      end if
      allocate(values(n),vectors(n,n),scaled(n))
      call eigen_symmetric(0.5_dp*(c+transpose(c)),values,vectors,info)
      if(info/=0) then
         bb=0.0_dp; if(present(status)) status=nmof_numerical_failure; return
      end if
      values=max(values,floor_value)
      bb=matmul(vectors,matmul(diagonal_matrix(values),transpose(vectors)))
      do i=1,n
         scaled(i)=1.0_dp/sqrt(max(bb(i,i),tiny(1.0_dp)))
      end do
      do j=1,n
         do i=1,n
            bb(i,j)=bb(i,j)*scaled(i)*scaled(j)
         end do
      end do
      bb=0.5_dp*(bb+transpose(bb))
      do i=1,n
         bb(i,i)=1.0_dp
      end do
      if(present(status)) status=nmof_ok
   contains
      pure function diagonal_matrix(vals) result(d)
         real(dp),intent(in)::vals(:)
         real(dp)::d(size(vals),size(vals))
         integer::k
         d=0.0_dp
         do k=1,size(vals)
            d(k,k)=vals(k)
         end do
      end function diagonal_matrix
   end function repair_matrix

   subroutine column_subset(x,columns,multiplier,rank,status,tol)
      real(dp),intent(in)::x(:,:)
      integer,allocatable,intent(out)::columns(:)
      real(dp),allocatable,intent(out)::multiplier(:,:)
      integer,intent(out)::rank,status
      real(dp),intent(in),optional::tol
      call matrix_rank_subset(x,columns,multiplier,rank,status,tol)
   end subroutine column_subset

   function qtable_statistics(x) result(tab)
      real(dp),intent(in)::x(:,:)
      type(qtable_result)::tab
      real(dp)::q1,q2,q3,iqr
      integer::j
      allocate(tab%whiskers(5,size(x,2)),tab%median(size(x,2)),tab%minimum(size(x,2)),tab%maximum(size(x,2)))
      do j=1,size(x,2)
         q1=quantile_type7(x(:,j),0.25_dp); q2=quantile_type7(x(:,j),0.5_dp); q3=quantile_type7(x(:,j),0.75_dp)
         iqr=abs(q3-q1); tab%minimum(j)=minval(x(:,j)); tab%maximum(j)=maxval(x(:,j)); tab%median(j)=q2
         tab%whiskers(:,j)=[max(q1-1.5_dp*iqr,tab%minimum(j)),q1,q2,q3,min(q3+1.5_dp*iqr,tab%maximum(j))]
      end do
      tab%status=nmof_ok
   end function qtable_statistics

   function probability_backtest_overfitting(m,s,metric,threshold,context) result(ans)
      real(dp),intent(in)::m(:,:)
      integer,intent(in),optional::s
      procedure(metric_function),optional::metric
      real(dp),intent(in),optional::threshold
      class(*),intent(in),optional::context
      type(pbo_result)::ans
      integer::ss,t,n,half,ncomb,j,k,step,best,nless,nequal
      integer,allocatable::comb(:),starts(:),ends(:),idx_in(:),idx_out(:)
      real(dp),allocatable::fi(:),fo(:)
      real(dp)::th,omega
      ss=12; if(present(s)) ss=nint(real(s,dp)); if(mod(ss,2)/=0) ss=ss+1
      th=0.0_dp; if(present(threshold)) th=threshold
      t=size(m,1); n=size(m,2); half=ss/2
      if(ss<2.or.t<ss.or.n<1) then; ans%status=nmof_invalid_input; return; end if
      ncomb=binomial_count(ss,half); allocate(ans%lambda(ncomb),ans%in_sample(ncomb),ans%out_of_sample(ncomb))
      allocate(comb(half),starts(ss),ends(ss),fi(n),fo(n)); comb=[(j,j=1,half)]
      step=max(1,nint(real(t,dp)/real(ss,dp)))
      do j=1,ss; starts(j)=1+(j-1)*step; end do
      do j=1,ss-1; ends(j)=starts(j+1)-1; end do; ends(ss)=t
      do k=1,ncomb
         call build_indices(comb,.true.,starts,ends,idx_in); call build_indices(comb,.false.,starts,ends,idx_out)
         do j=1,n
            fi(j)=apply_metric(m(idx_in,j)); fo(j)=apply_metric(m(idx_out,j))
         end do
         best=maxloc(fi,dim=1)
         nless=count(fo<fo(best)-sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(fo(best))))
         nequal=count(abs(fo-fo(best))<=sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(fo(best))))
         omega=(real(nless,dp)+0.5_dp*real(nequal+1,dp))/real(n+1,dp)
         ans%lambda(k)=log(omega/(1.0_dp-omega)); ans%in_sample(k)=fi(best); ans%out_of_sample(k)=fo(best)
         call next_combination(comb,ss)
      end do
      ans%pbo=real(count(ans%lambda<=th),dp)/real(ncomb,dp); ans%status=nmof_ok
   contains
      function apply_metric(x) result(v)
         real(dp),intent(in)::x(:); real(dp)::v
         if(present(metric)) then
            if(present(context)) then; v=metric(x,context); else; v=metric(x); end if
         else
            v=sum(x)/real(size(x),dp)
         end if
      end function apply_metric
      subroutine build_indices(c,use_selected,st,en,idx)
         integer,intent(in)::c(:),st(:),en(:); logical,intent(in)::use_selected
         integer,allocatable,intent(out)::idx(:)
         logical::sel(size(st)); integer::ii,jj,p,total
         sel=.false.; sel(c)=.true.; if(.not.use_selected) sel=.not.sel
         total=0; do ii=1,size(st); if(sel(ii)) total=total+en(ii)-st(ii)+1; end do
         allocate(idx(total)); p=0
         do ii=1,size(st); if(sel(ii)) then; do jj=st(ii),en(ii); p=p+1; idx(p)=jj; end do; end if; end do
      end subroutine build_indices
      subroutine next_combination(c,nmax)
         integer,intent(inout)::c(:); integer,intent(in)::nmax; integer::ii,jj
         ii=size(c)
         do while(ii>=1)
            if(c(ii)/=nmax-size(c)+ii) exit
            ii=ii-1
         end do
         if(ii>=1) then
            c(ii)=c(ii)+1
            do jj=ii+1,size(c)
               c(jj)=c(jj-1)+1
            end do
         end if
      end subroutine next_combination
      pure integer function binomial_count(nn,kk) result(v)
         integer,intent(in)::nn,kk; integer::ii,k2
         k2=min(kk,nn-kk); v=1
         do ii=1,k2; v=v*(nn-k2+ii)/ii; end do
      end function binomial_count
   end function probability_backtest_overfitting

   function marginal_risk_contributions(w,covariance,scale) result(ans)
      real(dp),intent(in)::w(:),covariance(:,:)
      logical,intent(in),optional::scale
      real(dp)::ans(size(w)),vol
      logical::sc
      sc=.true.; if(present(scale)) sc=scale
      ans=w*matmul(covariance,w)
      if(sc) then; vol=sqrt(max(0.0_dp,dot_product(w,matmul(covariance,w)))); if(vol>0) ans=ans/vol; end if
   end function marginal_risk_contributions

   function marginal_risk_contributions_fd(w,returns,risk,portfolio_returns,h,scale,context) result(ans)
      real(dp),intent(in)::w(:),returns(:,:)
      procedure(metric_function)::risk
      real(dp),intent(in),optional::portfolio_returns(:),h
      logical,intent(in),optional::scale
      class(*),intent(in),optional::context
      real(dp)::ans(size(w)),step,base,shifted
      real(dp),allocatable::rw(:)
      logical::sc
      integer::i
      step=1.0e-8_dp; if(present(h)) step=h
      sc=.true.; if(present(scale)) sc=scale
      if(size(returns,2)/=size(w).or.step<=0.0_dp) then
         ans=0.0_dp; return
      end if
      allocate(rw(size(returns,1)))
      if(present(portfolio_returns)) then
         if(size(portfolio_returns)/=size(rw)) then; ans=0.0_dp; return; end if
         rw=portfolio_returns
      else
         rw=matmul(returns,w)
      end if
      if(present(context)) then; base=risk(rw,context); else; base=risk(rw); end if
      do i=1,size(w)
         if(present(context)) then
            shifted=risk(rw+step*returns(:,i),context)
         else
            shifted=risk(rw+step*returns(:,i))
         end if
         ans(i)=w(i)*(shifted-base)/step
      end do
      if(sc.and.abs(base)>tiny(1.0_dp)) ans=ans/base
   end function marginal_risk_contributions_fd

   pure real(dp) function test_ackley(x) result(f)
      real(dp),intent(in)::x(:); integer::n
      n=size(x); f=20.0_dp+exp(1.0_dp)-20.0_dp*exp(-0.2_dp*sqrt(dot_product(x,x)/real(n,dp)))- &
         exp(sum(cos(2.0_dp*pi*x))/real(n,dp))
   end function test_ackley
   pure real(dp) function test_griewank(x) result(f)
      real(dp),intent(in)::x(:)
      real(dp)::prod
      integer::i
      prod=1.0_dp
      do i=1,size(x)
         prod=prod*cos(x(i)/sqrt(real(i,dp)))
      end do
      f=dot_product(x,x)/4000.0_dp-prod+1.0_dp
   end function test_griewank
   pure real(dp) function test_rastrigin(x) result(f)
      real(dp),intent(in)::x(:); f=10.0_dp*size(x)+sum(x*x-10.0_dp*cos(2.0_dp*pi*x))
   end function test_rastrigin
   pure real(dp) function test_rosenbrock(x) result(f)
      real(dp),intent(in)::x(:); integer::i; f=0.0_dp
      do i=1,size(x)-1; f=f+100.0_dp*(x(i+1)-x(i)**2)**2+(1.0_dp-x(i))**2; end do
   end function test_rosenbrock
   pure real(dp) function test_schwefel(x) result(f)
      real(dp),intent(in)::x(:); f=sum(-x*sin(sqrt(abs(x))))
   end function test_schwefel
   pure real(dp) function test_trefethen(x) result(f)
      real(dp),intent(in)::x(:); real(dp)::a,b
      a=x(1); b=x(2); f=exp(sin(50*a))+sin(60*exp(b))+sin(70*sin(a))+sin(sin(80*b))-sin(10*(a+b))+(a*a+b*b)/4
   end function test_trefethen
   pure real(dp) function test_eggholder(x) result(f)
      real(dp),intent(in)::x(:); real(dp)::a,b
      a=x(1); b=x(2); f=-(b+47)*sin(sqrt(abs(b+a/2+47)))-a*sin(sqrt(abs(a-(b+47))))
   end function test_eggholder

   pure logical function opposite_sign(a,b) result(opposite)
      real(dp),intent(in)::a,b
      opposite=(a<0.0_dp .and. b>=0.0_dp) .or. (a>=0.0_dp .and. b<0.0_dp)
   end function opposite_sign

   pure function lowercase(s) result(t)
      character(len=*),intent(in)::s; character(len=len(s))::t; integer::i,c
      do i=1,len(s); c=iachar(s(i:i)); if(c>=65.and.c<=90) then; t(i:i)=achar(c+32); else; t(i:i)=s(i:i); end if; end do
   end function lowercase
end module nmof_utilities
