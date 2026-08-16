module rfast_naive_bayes
   use rfast_special, only : dp
   use rfast_mle, only : mle_result, gamma_mle
   implicit none
   private

   type, public :: gaussian_nb_model
      real(dp), allocatable :: mean(:,:), variance(:,:), prior(:)
   end type gaussian_nb_model

   type, public :: poisson_nb_model
      real(dp), allocatable :: rate(:,:), prior(:)
   end type poisson_nb_model

   type, public :: geometric_nb_model
      real(dp), allocatable :: prob(:,:), prior(:)
      integer :: support_type = 1
   end type geometric_nb_model

   type, public :: multinomial_nb_model
      real(dp), allocatable :: prob(:,:), prior(:)
   end type multinomial_nb_model

   type, public :: gamma_nb_model
      real(dp), allocatable :: shape(:,:), rate(:,:), prior(:)
   end type gamma_nb_model

   public :: fit_gaussian_nb, predict_gaussian_nb
   public :: fit_poisson_nb, predict_poisson_nb
   public :: fit_geometric_nb, predict_geometric_nb
   public :: fit_multinomial_nb, predict_multinomial_nb
   public :: fit_gamma_nb, predict_gamma_nb

contains

   pure integer function class_count(y) result(k)
      integer, intent(in) :: y(:)
      if (size(y) == 0) then
         k = 0
      else
         k = maxval(y)
      end if
   end function class_count

   subroutine class_priors(y,k,prior,count)
      integer, intent(in) :: y(:), k
      real(dp), intent(out) :: prior(k)
      integer, intent(out) :: count(k)
      integer :: i
      count = 0
      do i=1,size(y)
         if (y(i)>=1 .and. y(i)<=k) count(y(i))=count(y(i))+1
      end do
      prior = real(count,dp)/real(max(1,size(y)),dp)
   end subroutine class_priors

   function fit_gaussian_nb(x,y) result(model)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: y(:)
      type(gaussian_nb_model) :: model
      integer :: k,p,c,j,i,nc
      integer, allocatable :: count(:)
      real(dp) :: s
      k=class_count(y);p=size(x,2)
      allocate(model%mean(k,p),model%variance(k,p),model%prior(k),count(k))
      call class_priors(y,k,model%prior,count)
      model%mean=0.0_dp;model%variance=0.0_dp
      do c=1,k
         nc=count(c)
         if(nc<=0) cycle
         do j=1,p
            s=0.0_dp
            do i=1,size(x,1)
               if(y(i)==c) s=s+x(i,j)
            end do
            model%mean(c,j)=s/real(nc,dp)
            s=0.0_dp
            do i=1,size(x,1)
               if(y(i)==c) s=s+(x(i,j)-model%mean(c,j))**2
            end do
            if(nc>1) then
               model%variance(c,j)=max(s/real(nc-1,dp),tiny(1.0_dp))
            else
               model%variance(c,j)=tiny(1.0_dp)
            end if
         end do
      end do
   end function fit_gaussian_nb

   function predict_gaussian_nb(model,x) result(cls)
      type(gaussian_nb_model), intent(in) :: model
      real(dp), intent(in) :: x(:,:)
      integer :: cls(size(x,1))
      integer :: i,c,k,best
      real(dp) :: score,bscore
      k=size(model%prior)
      do i=1,size(x,1)
         best=1;bscore=-huge(1.0_dp)
         do c=1,k
            if(model%prior(c)<=0.0_dp) cycle
            score=2.0_dp*log(model%prior(c))-sum(log(model%variance(c,:))) &
                  -sum((x(i,:)-model%mean(c,:))**2/model%variance(c,:))
            if(score>bscore)then;bscore=score;best=c;end if
         end do
         cls(i)=best
      end do
   end function predict_gaussian_nb

   function fit_poisson_nb(x,y) result(model)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: y(:)
      type(poisson_nb_model) :: model
      integer :: k,p,c,j,i
      integer, allocatable :: count(:)
      k=class_count(y);p=size(x,2)
      allocate(model%rate(k,p),model%prior(k),count(k))
      call class_priors(y,k,model%prior,count);model%rate=0.0_dp
      do c=1,k
         if(count(c)<=0) cycle
         do j=1,p
            do i=1,size(x,1)
               if(y(i)==c) model%rate(c,j)=model%rate(c,j)+x(i,j)
            end do
            model%rate(c,j)=max(model%rate(c,j)/real(count(c),dp),tiny(1.0_dp))
         end do
      end do
   end function fit_poisson_nb

   function predict_poisson_nb(model,x) result(cls)
      type(poisson_nb_model), intent(in) :: model
      real(dp), intent(in) :: x(:,:)
      integer :: cls(size(x,1))
      integer :: i,c,best
      real(dp) :: score,bscore
      do i=1,size(x,1)
         best=1;bscore=-huge(1.0_dp)
         do c=1,size(model%prior)
            if(model%prior(c)<=0.0_dp) cycle
            score=log(model%prior(c))+sum(x(i,:)*log(model%rate(c,:))-model%rate(c,:))
            if(score>bscore)then;bscore=score;best=c;end if
         end do
         cls(i)=best
      end do
   end function predict_poisson_nb

   function fit_geometric_nb(x,y,support_type) result(model)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: y(:)
      integer, intent(in), optional :: support_type
      type(geometric_nb_model) :: model
      integer :: k,p,c,j,i,t
      integer, allocatable :: count(:)
      real(dp) :: sx
      t=1;if(present(support_type))t=support_type
      k=class_count(y);p=size(x,2)
      allocate(model%prob(k,p),model%prior(k),count(k));model%support_type=t
      call class_priors(y,k,model%prior,count);model%prob=0.5_dp
      do c=1,k
         if(count(c)<=0) cycle
         do j=1,p
            sx=0.0_dp
            do i=1,size(x,1)
               if(y(i)==c) sx=sx+x(i,j)
            end do
            if(t==1)then
               model%prob(c,j)=1.0_dp/(1.0_dp+sx/real(count(c),dp))
            else
               model%prob(c,j)=real(count(c),dp)/max(sx,real(count(c),dp)+tiny(1.0_dp))
            end if
            model%prob(c,j)=min(1.0_dp-tiny(1.0_dp),max(tiny(1.0_dp),model%prob(c,j)))
         end do
      end do
   end function fit_geometric_nb

   function predict_geometric_nb(model,x) result(cls)
      type(geometric_nb_model), intent(in) :: model
      real(dp), intent(in) :: x(:,:)
      integer :: cls(size(x,1))
      integer :: i,c,best
      real(dp) :: score,bscore
      do i=1,size(x,1)
         best=1;bscore=-huge(1.0_dp)
         do c=1,size(model%prior)
            if(model%prior(c)<=0.0_dp) cycle
            score=log(model%prior(c))+sum(log(model%prob(c,:))+x(i,:)*log(1.0_dp-model%prob(c,:)))
            if(score>bscore)then;bscore=score;best=c;end if
         end do
         cls(i)=best
      end do
   end function predict_geometric_nb

   function fit_multinomial_nb(x,y) result(model)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: y(:)
      type(multinomial_nb_model) :: model
      integer :: k,p,c,i
      integer, allocatable :: count(:)
      real(dp) :: sx
      k=class_count(y);p=size(x,2)
      allocate(model%prob(k,p),model%prior(k),count(k));model%prob=0.0_dp
      call class_priors(y,k,model%prior,count)
      do c=1,k
         do i=1,size(x,1)
            if(y(i)/=c) cycle
            sx=sum(x(i,:))
            if(sx>0.0_dp) model%prob(c,:)=model%prob(c,:)+x(i,:)/sx
         end do
         if(count(c)>0) model%prob(c,:)=model%prob(c,:)/real(count(c),dp)
         model%prob(c,:)=max(model%prob(c,:),tiny(1.0_dp))
         model%prob(c,:)=model%prob(c,:)/sum(model%prob(c,:))
      end do
   end function fit_multinomial_nb

   function predict_multinomial_nb(model,x) result(cls)
      type(multinomial_nb_model), intent(in) :: model
      real(dp), intent(in) :: x(:,:)
      integer :: cls(size(x,1))
      integer :: i,c,best
      real(dp) :: score,bscore
      do i=1,size(x,1)
         best=1;bscore=-huge(1.0_dp)
         do c=1,size(model%prior)
            if(model%prior(c)<=0.0_dp) cycle
            score=log(model%prior(c))+sum(x(i,:)*log(model%prob(c,:)))
            if(score>bscore)then;bscore=score;best=c;end if
         end do
         cls(i)=best
      end do
   end function predict_multinomial_nb

   function fit_gamma_nb(x,y) result(model)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: y(:)
      type(gamma_nb_model) :: model
      type(mle_result) :: fit
      integer :: k,p,c,j,nc,ii,i
      integer, allocatable :: count(:)
      real(dp), allocatable :: z(:)
      k=class_count(y);p=size(x,2)
      allocate(model%shape(k,p),model%rate(k,p),model%prior(k),count(k))
      call class_priors(y,k,model%prior,count)
      model%shape=1.0_dp;model%rate=1.0_dp
      do c=1,k
         nc=count(c);if(nc<=0)cycle
         allocate(z(nc))
         do j=1,p
            ii=0
            do i=1,size(x,1)
               if(y(i)==c)then;ii=ii+1;z(ii)=x(i,j);end if
            end do
            fit=gamma_mle(z)
            if(allocated(fit%param).and.size(fit%param)>=2)then
               model%shape(c,j)=fit%param(1);model%rate(c,j)=fit%param(2)
            end if
         end do
         deallocate(z)
      end do
   end function fit_gamma_nb

   function predict_gamma_nb(model,x) result(cls)
      type(gamma_nb_model), intent(in) :: model
      real(dp), intent(in) :: x(:,:)
      integer :: cls(size(x,1))
      integer :: i,c,j,best
      real(dp) :: score,bscore,a,b
      do i=1,size(x,1)
         best=1;bscore=-huge(1.0_dp)
         do c=1,size(model%prior)
            if(model%prior(c)<=0.0_dp.or.any(x(i,:)<=0.0_dp)) cycle
            score=log(model%prior(c))
            do j=1,size(x,2)
               a=model%shape(c,j);b=model%rate(c,j)
               score=score+a*log(b)-log_gamma(a)+(a-1.0_dp)*log(x(i,j))-b*x(i,j)
            end do
            if(score>bscore)then;bscore=score;best=c;end if
         end do
         cls(i)=best
      end do
   end function predict_gamma_nb

end module rfast_naive_bayes
