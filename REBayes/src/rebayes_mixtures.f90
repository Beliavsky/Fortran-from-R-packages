module rebayes_mixtures
   use rebayes_kinds, only : dp
   use rebayes_math
   use rebayes_kw, only : kw_result, kw_control, kw_fit
   implicit none
   private
   public :: mixture_fit, bivariate_mixture_fit
   public :: glmix, bmix, b2mix, bpmix, pmix, gvmix, gammamix, weibullmix
   public :: gompertzmix, tlmix, tncpmix, umix, hlmix, npmix, cosslett
   public :: glvmix

   type :: mixture_fit
      real(dp), allocatable :: grid(:), mass(:), g(:), posterior_mean(:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: kkt_gap = huge(1.0_dp)
      integer :: status = 1
      integer :: iterations = 0
   end type mixture_fit

   type :: bivariate_mixture_fit
      real(dp), allocatable :: u(:), v(:), mass(:), g(:), post_u(:), post_v(:), post_product(:)
      real(dp), allocatable :: a(:,:)
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: status = 1
   end type bivariate_mixture_fit
contains
   subroutine make_weights(n, win, w)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: win(:)
      real(dp), allocatable, intent(out) :: w(:)
      allocate(w(n))
      if (present(win)) then
         if (size(win)/=n) error stop "observation weight size mismatch"
         w=max(win,0.0_dp); call normalize_prob(w)
      else
         w=1.0_dp/real(n,dp)
      end if
   end subroutine make_weights

   subroutine set_result(grid,a,w,norig,kw,fit)
      real(dp), intent(in) :: grid(:), a(:,:), w(:)
      integer, intent(in) :: norig
      type(kw_result), intent(in) :: kw
      type(mixture_fit), intent(out) :: fit
      allocate(fit%grid(size(grid)),fit%mass(size(grid)),fit%g(size(a,1)),fit%posterior_mean(size(a,1)))
      fit%grid=grid; fit%mass=kw%f; fit%g=kw%g
      fit%posterior_mean=matmul(a,kw%f*grid)/max(kw%g,tiny(1.0_dp))
      fit%loglik=real(norig,dp)*sum(w*safe_log(kw%g))
      fit%status=kw%status; fit%iterations=kw%iterations; fit%kkt_gap=kw%kkt_gap
   end subroutine set_result

   subroutine glmix(x,v,sigma,fit,weights,control)
      real(dp), intent(in) :: x(:),v(:),sigma(:)
      type(mixture_fit), intent(out) :: fit
      real(dp), intent(in), optional :: weights(:)
      type(kw_control), intent(in), optional :: control
      real(dp), allocatable :: a(:,:),w(:),d(:)
      type(kw_result) :: kw
      integer :: i,j,n,m
      n=size(x);m=size(v)
      if(size(sigma)/=1 .and. size(sigma)/=n) error stop "glmix: sigma size"
      allocate(a(n,m),d(m)); d=1.0_dp; call make_weights(n,weights,w)
      do j=1,m; do i=1,n
         a(i,j)=normal_pdf(x(i),v(j),sigma(merge(1,i,size(sigma)==1)))
      end do; end do
      call kw_fit(a,d,w,kw,control); call set_result(v,a,w,n,kw,fit)
   end subroutine glmix

   subroutine bmix(x,k,v,fit,weights,control)
      integer, intent(in) :: x(:),k(:)
      real(dp), intent(in) :: v(:)
      type(mixture_fit), intent(out) :: fit
      real(dp), intent(in), optional :: weights(:)
      type(kw_control), intent(in), optional :: control
      real(dp), allocatable :: a(:,:),w(:),d(:)
      type(kw_result) :: kw
      integer :: i,j,n,m
      n=size(x);m=size(v); if(size(k)/=n) error stop "bmix: k size"
      allocate(a(n,m),d(m)); d=1.0_dp; call make_weights(n,weights,w)
      do j=1,m; do i=1,n; a(i,j)=binomial_pmf(x(i),k(i),v(j)); end do; end do
      call kw_fit(a,d,w,kw,control); call set_result(v,a,w,n,kw,fit)
   end subroutine bmix

   subroutine b2mix(x1,x2,k1,k2,u,v,fit,weights,control)
      integer,intent(in)::x1(:),x2(:),k1(:),k2(:)
      real(dp),intent(in)::u(:),v(:)
      type(bivariate_mixture_fit),intent(out)::fit
      real(dp),intent(in),optional::weights(:)
      type(kw_control),intent(in),optional::control
      real(dp),allocatable::a(:,:),w(:),d(:),guv(:),ug(:),vg(:)
      type(kw_result)::kw
      integer::i,j,l,c,n,mu,mv
      n=size(x1);mu=size(u);mv=size(v)
      if(any([size(x2),size(k1),size(k2)]/=n)) error stop "b2mix: size"
      allocate(a(n,mu*mv),d(mu*mv),guv(mu*mv),ug(mu*mv),vg(mu*mv));d=1.0_dp
      call make_weights(n,weights,w); c=0
      do l=1,mv; do j=1,mu; c=c+1; ug(c)=u(j);vg(c)=v(l)
         do i=1,n
            a(i,c)=binomial_pmf(x1(i),k1(i),u(j))*binomial_pmf(x2(i),k2(i),v(l))
         end do
      end do; end do
      call kw_fit(a,d,w,kw,control)
      allocate(fit%u(mu),fit%v(mv),fit%mass(mu*mv),fit%g(n),fit%post_u(n),fit%post_v(n),fit%post_product(n),fit%a(n,mu*mv))
      fit%u=u;fit%v=v;fit%mass=kw%f;fit%g=kw%g;fit%a=a
      fit%post_u=matmul(a,kw%f*ug)/max(kw%g,tiny(1.0_dp))
      fit%post_v=matmul(a,kw%f*vg)/max(kw%g,tiny(1.0_dp))
      fit%post_product=matmul(a,kw%f*ug*vg)/max(kw%g,tiny(1.0_dp))
      fit%loglik=real(n,dp)*sum(w*safe_log(kw%g));fit%status=kw%status
   end subroutine b2mix

   subroutine bpmix(x,trials,v,u,fit,weights,control)
      integer,intent(in)::x(:),trials(:)
      real(dp),intent(in)::v(:),u(:)
      type(bivariate_mixture_fit),intent(out)::fit
      real(dp),intent(in),optional::weights(:)
      type(kw_control),intent(in),optional::control
      real(dp),allocatable::a(:,:),w(:),d(:),vg(:),ug(:)
      type(kw_result)::kw
      integer::i,j,l,c,n
      n=size(x); if(size(trials)/=n) error stop "bpmix: trials"
      allocate(a(n,size(v)*size(u)),d(size(v)*size(u)),vg(size(v)*size(u)),ug(size(v)*size(u)));d=1
      call make_weights(n,weights,w);c=0
      do l=1,size(u);do j=1,size(v);c=c+1;vg(c)=v(j);ug(c)=u(l)
         do i=1,n
            a(i,c)=binomial_pmf(x(i),trials(i),v(j))*poisson_pmf(trials(i),u(l))
         end do
      end do;end do
      call kw_fit(a,d,w,kw,control)
      allocate(fit%u(size(v)),fit%v(size(u)),fit%mass(size(d)))
      allocate(fit%g(n),fit%post_u(n),fit%post_v(n),fit%post_product(n),fit%a(n,size(d)))
      fit%u=v;fit%v=u;fit%mass=kw%f;fit%g=kw%g;fit%a=a
      fit%post_u=matmul(a,kw%f*vg)/max(kw%g,tiny(1.0_dp));fit%post_v=matmul(a,kw%f*ug)/max(kw%g,tiny(1.0_dp))
      fit%post_product=matmul(a,kw%f*vg*ug)/max(kw%g,tiny(1.0_dp))
      fit%loglik=real(n,dp)*sum(w*safe_log(kw%g));fit%status=kw%status
   end subroutine bpmix

   subroutine pmix(x,v,fit,exposure,weights,support,control)
      integer,intent(in)::x(:);real(dp),intent(in)::v(:)
      type(mixture_fit),intent(out)::fit
      real(dp),intent(in),optional::exposure(:),weights(:)
      integer,intent(in),optional::support(2)
      type(kw_control),intent(in),optional::control
      real(dp),allocatable::a(:,:),w(:),d(:),expo(:)
      type(kw_result)::kw
      integer::i,j,n,m
      n=size(x);m=size(v);allocate(a(n,m),d(m),expo(n));d=1
      if(present(exposure)) then
         if(size(exposure)==1) then;expo=exposure(1)
         else if(size(exposure)==n) then;expo=exposure
         else;error stop "pmix: exposure";end if
      else;expo=1;end if
      call make_weights(n,weights,w)
      do j=1,m;do i=1,n
         a(i,j)=poisson_pmf(x(i),v(j)*expo(i))
         if(present(support)) a(i,j)=a(i,j)/max(poisson_cdf(support(2),v(j)*expo(i))- &
            poisson_cdf(support(1),v(j)*expo(i)),tiny(1.0_dp))
      end do;end do
      call kw_fit(a,d,w,kw,control);call set_result(v,a,w,n,kw,fit)
   end subroutine pmix

   subroutine gvmix(x,mobs,v,fit,weights,control)
      real(dp),intent(in)::x(:),mobs(:),v(:);type(mixture_fit),intent(out)::fit
      real(dp),intent(in),optional::weights(:);type(kw_control),intent(in),optional::control
      real(dp),allocatable::a(:,:),w(:),d(:);type(kw_result)::kw
      integer::i,j,n
      n=size(x);if(size(mobs)/=n)error stop "gvmix: m"
      allocate(a(n,size(v)),d(size(v)));d=1;call make_weights(n,weights,w)
      do j=1,size(v);do i=1,n
         a(i,j)=gamma_pdf(x(i),0.5_dp*(mobs(i)-1.0_dp),v(j)/(0.5_dp*(mobs(i)-1.0_dp)))
      end do;end do
      call kw_fit(a,d,w,kw,control);call set_result(v,a,w,n,kw,fit)
   end subroutine gvmix

   subroutine gammamix(x,v,shape,fit,weights,control)
      real(dp),intent(in)::x(:),v(:),shape;type(mixture_fit),intent(out)::fit
      real(dp),intent(in),optional::weights(:);type(kw_control),intent(in),optional::control
      real(dp),allocatable::a(:,:),w(:),d(:);type(kw_result)::kw;integer::i,j,n
      n=size(x);allocate(a(n,size(v)),d(size(v)));d=1;call make_weights(n,weights,w)
      do j=1,size(v);do i=1,n;a(i,j)=gamma_pdf(x(i),shape,v(j));end do;end do
      call kw_fit(a,d,w,kw,control);call set_result(v,a,w,n,kw,fit)
   end subroutine gammamix

   subroutine weibullmix(x,v,alpha,lambda,fit,event,weights,control)
      real(dp),intent(in)::x(:),v(:),alpha,lambda(:);type(mixture_fit),intent(out)::fit
      integer,intent(in),optional::event(:);real(dp),intent(in),optional::weights(:)
      type(kw_control),intent(in),optional::control
      real(dp),allocatable::a(:,:),w(:),d(:);type(kw_result)::kw;real(dp)::scale
      integer::i,j,n,e
      n=size(x);if(size(lambda)/=1.and.size(lambda)/=n)error stop "weibullmix: lambda"
      if(present(event))then;if(size(event)/=n)error stop "weibullmix: event";end if
      allocate(a(n,size(v)),d(size(v)));d=1;call make_weights(n,weights,w)
      do j=1,size(v);do i=1,n
         scale=(1.0_dp/lambda(merge(1,i,size(lambda)==1)))*exp(-v(j)/alpha)
         e=1;if(present(event))e=event(i)
         if(e/=0)then;a(i,j)=weibull_pdf(x(i),alpha,scale)
         else;a(i,j)=1.0_dp-weibull_cdf(x(i),alpha,scale);end if
      end do;end do
      call kw_fit(a,d,w,kw,control);call set_result(v,a,w,n,kw,fit)
   end subroutine weibullmix

   subroutine gompertzmix(x,v,alpha,theta,fit,weights,control)
      real(dp),intent(in)::x(:),v(:),alpha,theta;type(mixture_fit),intent(out)::fit
      real(dp),intent(in),optional::weights(:);type(kw_control),intent(in),optional::control
      real(dp),allocatable::a(:,:),w(:),d(:);type(kw_result)::kw;integer::i,j,n
      n=size(x);allocate(a(n,size(v)),d(size(v)));d=1;call make_weights(n,weights,w)
      do j=1,size(v);do i=1,n;a(i,j)=gompertz_pdf(x(i),alpha,theta*exp(v(j)));end do;end do
      call kw_fit(a,d,w,kw,control);call set_result(v,a,w,n,kw,fit)
   end subroutine gompertzmix

   subroutine tlmix(x,v,df,fit,weights,control)
      real(dp),intent(in)::x(:),v(:),df;type(mixture_fit),intent(out)::fit
      real(dp),intent(in),optional::weights(:);type(kw_control),intent(in),optional::control
      real(dp),allocatable::a(:,:),w(:),d(:);type(kw_result)::kw;integer::i,j,n
      n=size(x);allocate(a(n,size(v)),d(size(v)));d=1;call make_weights(n,weights,w)
      do j=1,size(v);do i=1,n;a(i,j)=student_pdf(x(i)-v(j),df);end do;end do
      call kw_fit(a,d,w,kw,control);call set_result(v,a,w,n,kw,fit)
   end subroutine tlmix

   subroutine tncpmix(x,v,df,fit,weights,control)
      real(dp),intent(in)::x(:),v(:),df(:);type(mixture_fit),intent(out)::fit
      real(dp),intent(in),optional::weights(:);type(kw_control),intent(in),optional::control
      real(dp),allocatable::a(:,:),w(:),d(:);type(kw_result)::kw;integer::i,j,n
      n=size(x);if(size(df)/=1.and.size(df)/=n)error stop "tncpmix: df"
      allocate(a(n,size(v)),d(size(v)));d=1;call make_weights(n,weights,w)
      do j=1,size(v);do i=1,n
         a(i,j)=noncentral_t_pdf(x(i),df(merge(1,i,size(df)==1)),v(j))
      end do;end do
      call kw_fit(a,d,w,kw,control);call set_result(v,a,w,n,kw,fit)
   end subroutine tncpmix

   subroutine umix(x,fit,control)
      real(dp),intent(in)::x(:);type(mixture_fit),intent(out)::fit;type(kw_control),intent(in),optional::control
      real(dp),allocatable::v(:),a(:,:),w(:),d(:);type(kw_result)::kw;integer::i,j,n
      n=size(x);allocate(v(n));v=x;call sort_real(v);allocate(a(n,n),w(n),d(n));w=1.0_dp/real(n,dp);d=1
      do j=1,n;do i=1,n;a(i,j)=merge(1.0_dp/v(j),0.0_dp,x(i)<=v(j).and.v(j)>0);end do;end do
      call kw_fit(a,d,w,kw,control);call set_result(v,a,w,n,kw,fit)
   end subroutine umix

   subroutine hlmix(x,v,sigma,k,fit,weights,control)
      real(dp),intent(in)::x(:),v(:),sigma,k;type(mixture_fit),intent(out)::fit
      real(dp),intent(in),optional::weights(:);type(kw_control),intent(in),optional::control
      real(dp),allocatable::a(:,:),w(:),d(:);type(kw_result)::kw;real(dp)::heps;integer::i,j,n
      n=size(x);heps=huber_eps(k);allocate(a(n,size(v)),d(size(v)));d=1;call make_weights(n,weights,w)
      do j=1,size(v);do i=1,n;a(i,j)=huber_pdf(x(i)-v(j),sigma,k,heps);end do;end do
      call kw_fit(a,d,w,kw,control);call set_result(v,a,w,n,kw,fit)
   end subroutine hlmix

   subroutine npmix(x,trials,v,u,fit,weights,control)
      real(dp),intent(in)::x(:);integer,intent(in)::trials(:);real(dp),intent(in)::v(:),u(:)
      type(bivariate_mixture_fit),intent(out)::fit;real(dp),intent(in),optional::weights(:)
      type(kw_control),intent(in),optional::control
      real(dp),allocatable::a(:,:),w(:),d(:),vg(:),ug(:);type(kw_result)::kw;integer::i,j,l,c,n
      n=size(x);if(size(trials)/=n)error stop "npmix: trials"
      allocate(a(n,size(v)*size(u)),d(size(v)*size(u)))
      allocate(vg(size(v)*size(u)),ug(size(v)*size(u)))
      d=1;call make_weights(n,weights,w);c=0
      do l=1,size(u);do j=1,size(v);c=c+1;vg(c)=v(j);ug(c)=u(l)
         do i=1,n;a(i,c)=normal_pdf(x(i),v(j),1.0_dp/sqrt(real(trials(i),dp)))*poisson_pmf(trials(i),u(l));end do
      end do;end do
      call kw_fit(a,d,w,kw,control)
      allocate(fit%u(size(v)),fit%v(size(u)),fit%mass(size(d)))
      allocate(fit%g(n),fit%post_u(n),fit%post_v(n),fit%post_product(n),fit%a(n,size(d)))
      fit%u=v;fit%v=u;fit%mass=kw%f;fit%g=kw%g;fit%a=a;fit%post_u=matmul(a,kw%f*vg)/max(kw%g,tiny(1.0_dp))
      fit%post_v=matmul(a,kw%f*ug)/max(kw%g,tiny(1.0_dp))
      fit%post_product=matmul(a,kw%f*vg*ug)/max(kw%g,tiny(1.0_dp))
      fit%loglik=real(n,dp)*sum(w*safe_log(kw%g));fit%status=kw%status
   end subroutine npmix

   subroutine cosslett(x,y,v,fit,weights,control)
      real(dp),intent(in)::x(:),v(:);integer,intent(in)::y(:);type(mixture_fit),intent(out)::fit
      real(dp),intent(in),optional::weights(:);type(kw_control),intent(in),optional::control
      real(dp),allocatable::a(:,:),w(:),d(:);type(kw_result)::kw;integer::i,j,n
      n=size(x);if(size(y)/=n)error stop "cosslett: y";allocate(a(n,size(v)),d(size(v)));d=1;call make_weights(n,weights,w)
      do j=1,size(v);do i=1,n
         if(y(i)==1)then;a(i,j)=merge(1.0_dp,0.0_dp,x(i)>=v(j));else;a(i,j)=merge(1.0_dp,0.0_dp,x(i)<v(j));end if
      end do;end do
      call kw_fit(a,d,w,kw,control);call set_result(v,a,w,n,kw,fit)
   end subroutine cosslett

   subroutine glvmix(t,s,mobs,u,v,fit,control)
      real(dp),intent(in)::t(:),s(:),mobs(:),u(:),v(:);type(bivariate_mixture_fit),intent(out)::fit
      type(kw_control),intent(in),optional::control
      real(dp),allocatable::a(:,:),w(:),d(:),ug(:),vg(:);type(kw_result)::kw
      real(dp)::r,sd;integer::i,j,l,c,n
      n=size(t);if(size(s)/=n.or.size(mobs)/=n)error stop "glvmix: size"
      allocate(a(n,size(u)*size(v)),w(n),d(size(u)*size(v)),ug(size(u)*size(v)),vg(size(u)*size(v)));w=1.0_dp/real(n,dp);d=1;c=0
      do l=1,size(v);do j=1,size(u);c=c+1;ug(c)=u(j);vg(c)=v(l)
         do i=1,n
            r=0.5_dp*(mobs(i)-1.0_dp);sd=sqrt(v(l)/mobs(i))
            a(i,c)=gamma_pdf(s(i),r,v(l)/r)*normal_pdf(t(i),u(j),sd)
         end do
      end do;end do
      call kw_fit(a,d,w,kw,control)
      allocate(fit%u(size(u)),fit%v(size(v)),fit%mass(size(d)))
      allocate(fit%g(n),fit%post_u(n),fit%post_v(n),fit%post_product(n),fit%a(n,size(d)))
      fit%u=u;fit%v=v;fit%mass=kw%f;fit%g=kw%g;fit%a=a;fit%post_u=matmul(a,kw%f*ug)/max(kw%g,tiny(1.0_dp))
      fit%post_v=matmul(a,kw%f*vg)/max(kw%g,tiny(1.0_dp))
      fit%post_product=matmul(a,kw%f*ug*vg)/max(kw%g,tiny(1.0_dp))
      fit%loglik=real(n,dp)*sum(w*safe_log(kw%g));fit%status=kw%status
   end subroutine glvmix

   subroutine sort_real(x)
      real(dp),intent(inout)::x(:);real(dp)::key;integer::i,j
      do i=2,size(x);key=x(i);j=i-1;do while(j>=1);if(x(j)<=key)exit;x(j+1)=x(j);j=j-1;end do;x(j+1)=key;end do
   end subroutine sort_real
end module rebayes_mixtures
