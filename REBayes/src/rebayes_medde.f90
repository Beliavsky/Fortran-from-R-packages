module rebayes_medde
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use rebayes_kinds, only : dp
   implicit none
   private
   public :: medde_control, medde_result, medde_fit, medde_quantiles, medde_random
   type :: medde_control
      integer :: max_iter = 10000
      real(dp) :: tol = 1.0e-8_dp
      real(dp) :: initial_step = 1.0_dp
   end type medde_control
   type :: medde_result
      real(dp),allocatable :: x(:),y(:)
      real(dp) :: objective = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 1
   end type medde_result
contains
   subroutine medde_fit(data,grid,lambda,alpha,dorder,result,weights,mass,control)
      real(dp),intent(in)::data(:),grid(:),lambda,alpha
      integer,intent(in)::dorder
      type(medde_result),intent(out)::result
      real(dp),intent(in),optional::weights(:),mass
      type(medde_control),intent(in),optional::control
      type(medde_control)::ctl
      real(dp),allocatable::e(:),d(:),dv(:),amat(:,:),z(:),znew(:),f(:),fnew(:),grad(:),w(:)
      real(dp)::obj,objnew,step,chg,ms
      integer::p,q,it
      p=size(grid);if(p<=dorder)error stop "medde_fit: grid too small"
      ctl=medde_control();if(present(control))ctl=control
      allocate(e(p),d(p),dv(p),w(size(data)));call empirical_grid_mass(data,grid,weights,e)
      if(p>1)then
         d(1)=0.5_dp*(grid(2)-grid(1));d(p)=0.5_dp*(grid(p)-grid(p-1))
         if(p>2)d(2:p-1)=0.5_dp*(grid(3:p)-grid(1:p-2))
      else;d=1.0_dp;end if
      ms=1.0_dp;if(present(mass))ms=mass;dv=ms*d/sum(d)
      call difference_transpose(p,dorder,amat);q=size(amat,2)
      allocate(z(q),znew(q),f(p),fnew(p),grad(q));z=0.0_dp
      call medde_density(e,dv,amat,z,f);obj=medde_objective(f,dv,alpha)
      do it=1,ctl%max_iter
         call medde_gradient(f,amat,alpha,grad)
         step=ctl%initial_step/max(1.0_dp,sqrt(sum(grad*grad)))
         do
            znew=z+step*grad
            if(lambda>0.0_dp)then
               znew=max(-lambda,min(lambda,znew))
            else
               znew=max(0.0_dp,znew)
            end if
            call medde_density(e,dv,amat,znew,fnew)
            if(minval(fnew)>=-1.0e-12_dp)then
               objnew=medde_objective(max(fnew,0.0_dp),dv,alpha)
               if(objnew>=obj-1.0e-12_dp)exit
            end if
            step=0.5_dp*step
            if(step<1.0e-14_dp)exit
         end do
         if(step<1.0e-14_dp)exit
         chg=maxval(abs(znew-z));z=znew;f=max(fnew,0.0_dp);obj=objnew
         if(chg<ctl%tol*(1.0_dp+maxval(abs(z))))exit
      end do
      allocate(result%x(p),result%y(p));result%x=grid;result%y=f
      result%objective=obj;result%iterations=min(it,ctl%max_iter);result%status=merge(0,1,it<=ctl%max_iter)
   end subroutine medde_fit

   subroutine empirical_grid_mass(data,grid,weights,e)
      real(dp),intent(in)::data(:),grid(:)
      real(dp),intent(in),optional::weights(:)
      real(dp),intent(out)::e(:)
      real(dp),allocatable::w(:);real(dp)::t
      integer::i,j,n,p
      n=size(data);p=size(grid);allocate(w(n));if(present(weights))then
         if(size(weights)/=n)error stop "medde: weights";w=max(weights,0.0_dp);w=w/sum(w)
      else;w=1.0_dp/real(n,dp);end if
      e=0.0_dp
      do i=1,n
         if(data(i)<=grid(1))then;e(1)=e(1)+w(i)
         else if(data(i)>=grid(p))then;e(p)=e(p)+w(i)
         else
            j=1;do while(grid(j+1)<data(i));j=j+1;end do
            t=(data(i)-grid(j))/(grid(j+1)-grid(j));e(j)=e(j)+w(i)*(1.0_dp-t);e(j+1)=e(j+1)+w(i)*t
         end if
      end do
   end subroutine empirical_grid_mass

   subroutine difference_transpose(p,ord,a)
      integer,intent(in)::p,ord
      real(dp),allocatable,intent(out)::a(:,:)
      integer::q,j,k
      real(dp)::coef
      if(ord<1.or.ord>3)error stop "medde: dorder must be 1..3"
      q=p-ord;allocate(a(p,q));a=0.0_dp
      do j=1,q;do k=0,ord
         coef=real(binomial_int(ord,k),dp)
         if(mod(ord-k,2)==1)coef=-coef
         a(j+k,j)=coef
      end do;end do
   end subroutine difference_transpose

   integer function binomial_int(n,k) result(v)
      integer,intent(in)::n,k;integer::i
      v=1;do i=1,k;v=v*(n-k+i)/i;end do
   end function binomial_int

   subroutine medde_density(e,dv,a,z,f)
      real(dp),intent(in)::e(:),dv(:),a(:,:),z(:);real(dp),intent(out)::f(:)
      f=(e-matmul(a,z))/dv
   end subroutine medde_density

   real(dp) function medde_objective(f,dv,alpha) result(v)
      real(dp),intent(in)::f(:),dv(:),alpha
      real(dp)::beta
      integer::i
      v=0.0_dp
      if(abs(alpha-1.0_dp)<1.0e-12_dp)then
         do i=1,size(f);if(f(i)>0.0_dp)v=v-dv(i)*f(i)*log(f(i));end do
      else if(abs(alpha)<1.0e-12_dp)then
         do i=1,size(f);v=v+dv(i)*log(max(f(i),1.0e-14_dp));end do
      else
         beta=alpha/(alpha-1.0_dp)
         v=-sign(1.0_dp,beta)*sum(dv*max(f,1.0e-14_dp)**alpha)
      end if
   end function medde_objective

   subroutine medde_gradient(f,a,alpha,g)
      real(dp),intent(in)::f(:),a(:,:),alpha;real(dp),intent(out)::g(:)
      real(dp),allocatable::phi(:);real(dp)::beta
      integer::i
      allocate(phi(size(f)))
      if(abs(alpha-1.0_dp)<1.0e-12_dp)then
         phi=log(max(f,1.0e-14_dp))+1.0_dp;g=matmul(transpose(a),phi)
      else if(abs(alpha)<1.0e-12_dp)then
         phi=1.0_dp/max(f,1.0e-14_dp);g=-matmul(transpose(a),phi)
      else
         beta=alpha/(alpha-1.0_dp)
         phi=-sign(1.0_dp,beta)*alpha*max(f,1.0e-14_dp)**(alpha-1.0_dp)
         g=-matmul(transpose(a),phi)
      end if
      do i=1,size(g);if(ieee_is_nan(g(i)))g(i)=0.0_dp;end do
   end subroutine medde_gradient

   subroutine medde_quantiles(model,p,q)
      type(medde_result),intent(in)::model;real(dp),intent(in)::p(:);real(dp),intent(out)::q(:)
      real(dp),allocatable::cdf(:);real(dp)::area,t
      integer::i,j,n
      n=size(model%x);allocate(cdf(n));cdf=0.0_dp
      do i=2,n
         cdf(i)=cdf(i-1)+0.5_dp*(model%x(i)-model%x(i-1))*(model%y(i)+model%y(i-1))
      end do
      area=cdf(n);if(area>0)cdf=cdf/area
      do j=1,size(p)
         if(p(j)<=0)then;q(j)=model%x(1)
         else if(p(j)>=1)then;q(j)=model%x(n)
         else
            i=1;do while(i<n.and.cdf(i+1)<p(j));i=i+1;end do
            if(cdf(i+1)>cdf(i))then;t=(p(j)-cdf(i))/(cdf(i+1)-cdf(i));else;t=0;end if
            q(j)=model%x(i)+t*(model%x(i+1)-model%x(i))
         end if
      end do
   end subroutine medde_quantiles

   subroutine medde_random(model,z)
      type(medde_result),intent(in)::model;real(dp),intent(out)::z(:)
      real(dp),allocatable::u(:);allocate(u(size(z)));call random_number(u);call medde_quantiles(model,u,z)
   end subroutine medde_random
end module rebayes_medde
