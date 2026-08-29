module forecast_nnetar
   use forecast_kinds, only : dp
   use forecast_types, only : forecast_result
   use forecast_stats, only : mean_value, variance_value
   use nnet, only : nnet_model_t, nnet_fit, nnet_predict
   implicit none
   private
   public :: nnetar_model, nnetar_fit, nnetar_forecast
   type :: nnetar_model
      integer :: p=1,Pseason=0,m=1,size_hidden=1,repeats=1
      integer,allocatable :: lags(:)
      real(dp) :: center=0.0_dp,scale=1.0_dp,sigma2=0.0_dp
      real(dp),allocatable :: x(:),fitted(:),residuals(:)
      type(nnet_model_t),allocatable :: networks(:)
   end type
contains
   function nnetar_fit(y,m,p,Pseason,size_hidden,repeats,decay,maxit) result(model)
      real(dp),intent(in)::y(:)
      integer,intent(in),optional::m,p,Pseason,size_hidden,repeats,maxit
      real(dp),intent(in),optional::decay
      type(nnetar_model)::model
      integer::mm,pp,ps,nh,nr,mi,maxlag,nlag,i,j,k
      real(dp)::dec
      real(dp),allocatable::z(:),xx(:,:),yy(:,:),pred(:,:),row(:,:)
      integer,allocatable::lags(:),tmp(:)
      mm=1
      if(present(m))mm=max(1,m)
      pp=max(1,min(5,size(y)/10))
      if(present(p))pp=p
      ps=merge(1,0,mm>1)
      if(present(Pseason))ps=Pseason
      if(ps>0)then
         allocate(tmp(pp+ps))
         tmp(1:pp)=[(i,i=1,pp)]
         tmp(pp+1:)=[(i*mm,i=1,ps)]
         call unique_sorted(tmp,lags)
      else
      allocate(lags(pp))
      lags=[(i,i=1,pp)]
      end if
      maxlag=maxval(lags)
      nlag=size(lags)
      nh=nint(real(nlag+1,dp)/2.0_dp)
      if(present(size_hidden))nh=size_hidden
      nr=10
      if(present(repeats))nr=max(1,repeats)
      mi=200
      if(present(maxit))mi=maxit
      dec=0.0_dp
      if(present(decay))dec=decay
      model%center=mean_value(y)
      model%scale=sqrt(max(variance_value(y,.true.),1.0e-12_dp))
      z=(y-model%center)/model%scale
      allocate(xx(size(y)-maxlag,nlag),yy(size(y)-maxlag,1))
      do i=1,size(xx,1)
      k=maxlag+i
      do j=1,nlag
      xx(i,j)=z(k-lags(j))
      end do
      yy(i,1)=z(k)
      end do
      allocate(model%networks(nr))
      do i=1,nr
      call nnet_fit(model%networks(i),xx,yy,nh,linout=.true.,decay=[dec],maxit=mi)
      end do
      model%p=pp
      model%Pseason=ps
      model%m=mm
      model%size_hidden=nh
      model%repeats=nr
      model%lags=lags
      model%x=y
      allocate(model%fitted(size(y)),model%residuals(size(y)))
      model%fitted=y
      model%residuals=0.0_dp
      allocate(row(1,nlag))
      do i=maxlag+1,size(y)
         do j=1,nlag
         row(1,j)=z(i-lags(j))
         end do
         model%fitted(i)=0.0_dp
         do k=1,nr
         pred=nnet_predict(model%networks(k),row)
         model%fitted(i)=model%fitted(i)+pred(1,1)/real(nr,dp)
         end do
         model%fitted(i)=model%center+model%scale*model%fitted(i)
         model%residuals(i)=y(i)-model%fitted(i)
      end do
      if(size(y)>maxlag)model%sigma2=sum(model%residuals(maxlag+1:)**2)/real(size(y)-maxlag,dp)
   contains
      subroutine unique_sorted(a,b)
         integer,intent(in)::a(:)
         integer,allocatable,intent(out)::b(:)
         integer,allocatable::c(:)
         integer::ii,jj,n
         c=a
         do ii=2,size(c)
         jj=ii
         do while(jj>1.and.c(jj)<c(jj-1))
         n=c(jj)
         c(jj)=c(jj-1)
         c(jj-1)=n
         jj=jj-1
         end do
         end do
         n=1
         do ii=2,size(c)
         if(c(ii)/=c(n))then
         n=n+1
         c(n)=c(ii)
         end if
         end do
         allocate(b(n))
         b=c(1:n)
      end subroutine
   end function
   function nnetar_forecast(model,h) result(fc)
      type(nnetar_model),intent(in)::model
      integer,intent(in)::h
      type(forecast_result)::fc
      real(dp),allocatable::hist(:),row(:,:),pred(:,:)
      integer::i,j,k,n0,maxlag
      n0=size(model%x)
      maxlag=maxval(model%lags)
      allocate(hist(n0+h))
      hist(1:n0)=(model%x-model%center)/model%scale
      allocate(row(1,size(model%lags)),fc%mean(h))
      do i=1,h
         do j=1,size(model%lags)
         row(1,j)=hist(n0+i-model%lags(j))
         end do
         hist(n0+i)=0.0_dp
         do k=1,model%repeats
         pred=nnet_predict(model%networks(k),row)
         hist(n0+i)=hist(n0+i)+pred(1,1)/real(model%repeats,dp)
         end do
         fc%mean(i)=model%center+model%scale*hist(n0+i)
      end do
   end function
end module forecast_nnetar
