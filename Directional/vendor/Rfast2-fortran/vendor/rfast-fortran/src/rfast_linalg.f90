module rfast_linalg
   use rfast_special, only : dp
   use rfast_arrays, only : colmeans
   implicit none
   private
   public :: crossprod, tcrossprod, covariance_matrix, correlation_matrix
   public :: cholesky_lower, solve_linear, inverse_matrix, logdet_spd, determinant_matrix
   public :: eigen_sym_jacobi, mahalanobis_sq, pooled_covariance, spd_inverse
   public :: standard_covariance, matrix_rank

contains

   pure function crossprod(x,y) result(a)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(in),optional::y(:,:)
      real(dp),allocatable::a(:,:)
      if(present(y))then
         allocate(a(size(x,2),size(y,2)));a=matmul(transpose(x),y)
      else
         allocate(a(size(x,2),size(x,2)));a=matmul(transpose(x),x)
      end if
   end function crossprod

   pure function tcrossprod(x,y) result(a)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(in),optional::y(:,:)
      real(dp),allocatable::a(:,:)
      if(present(y))then
         allocate(a(size(x,1),size(y,1)));a=matmul(x,transpose(y))
      else
         allocate(a(size(x,1),size(x,1)));a=matmul(x,transpose(x))
      end if
   end function tcrossprod

   pure function covariance_matrix(x,population) result(s)
      real(dp),intent(in)::x(:,:)
      logical,intent(in),optional::population
      real(dp)::s(size(x,2),size(x,2)),m(size(x,2)),den
      real(dp)::z(size(x,1),size(x,2))
      logical::pop
      pop=.false.;if(present(population))pop=population
      m=colmeans(x);z=x-spread(m,1,size(x,1))
      if(pop)then;den=real(size(x,1),dp);else;den=real(size(x,1)-1,dp);end if
      s=matmul(transpose(z),z)/den
   end function covariance_matrix

   pure function correlation_matrix(x) result(r)
      real(dp),intent(in)::x(:,:)
      real(dp)::r(size(x,2),size(x,2)),s(size(x,2),size(x,2)),sd(size(x,2))
      integer::i,j
      s=covariance_matrix(x);sd=sqrt(max(0.0_dp,[(s(i,i),i=1,size(x,2))]))
      do j=1,size(x,2);do i=1,size(x,2)
         if(sd(i)>0.0_dp.and.sd(j)>0.0_dp)then;r(i,j)=s(i,j)/(sd(i)*sd(j));else;r(i,j)=0.0_dp;end if
      end do;end do
      do i=1,size(x,2);r(i,i)=1.0_dp;end do
   end function correlation_matrix

   subroutine cholesky_lower(a,l,info)
      real(dp),intent(in)::a(:,:);real(dp),intent(out)::l(size(a,1),size(a,2));integer,intent(out)::info
      integer::i,j,k,n;real(dp)::s
      n=size(a,1);l=0.0_dp;info=0
      if(size(a,2)/=n)then;info=-1;return;end if
      do i=1,n
         do j=1,i
            s=a(i,j)
            do k=1,j-1;s=s-l(i,k)*l(j,k);end do
            if(i==j)then
               if(s<=0.0_dp)then;info=i;return;end if
               l(i,j)=sqrt(s)
            else
               l(i,j)=s/l(j,j)
            end if
         end do
      end do
   end subroutine cholesky_lower

   subroutine solve_linear(a,b,x,info)
      real(dp),intent(in)::a(:,:),b(:);real(dp),intent(out)::x(size(b));integer,intent(out)::info
      real(dp)::m(size(a,1),size(a,2)),rhs(size(b)),factor,tmp
      integer::n,i,j,k,p
      n=size(b);info=0
      if(size(a,1)/=n.or.size(a,2)/=n)then;info=-1;return;end if
      m=a;rhs=b
      do k=1,n-1
         p=k
         do i=k+1,n;if(abs(m(i,k))>abs(m(p,k)))p=i;end do
         if(abs(m(p,k))<=tiny(1.0_dp))then;info=k;return;end if
         if(p/=k)then
            do j=k,n;tmp=m(k,j);m(k,j)=m(p,j);m(p,j)=tmp;end do
            tmp=rhs(k);rhs(k)=rhs(p);rhs(p)=tmp
         end if
         do i=k+1,n
            factor=m(i,k)/m(k,k);m(i,k)=0.0_dp
            do j=k+1,n;m(i,j)=m(i,j)-factor*m(k,j);end do
            rhs(i)=rhs(i)-factor*rhs(k)
         end do
      end do
      if(abs(m(n,n))<=tiny(1.0_dp))then;info=n;return;end if
      x(n)=rhs(n)/m(n,n)
      do i=n-1,1,-1
         x(i)=(rhs(i)-dot_product(m(i,i+1:n),x(i+1:n)))/m(i,i)
      end do
   end subroutine solve_linear

   subroutine inverse_matrix(a,ainv,info)
      real(dp),intent(in)::a(:,:);real(dp),intent(out)::ainv(size(a,1),size(a,2));integer,intent(out)::info
      real(dp)::e(size(a,1)),x(size(a,1));integer::j,n,inf
      n=size(a,1);ainv=0.0_dp;info=0
      do j=1,n
         e=0.0_dp;e(j)=1.0_dp;call solve_linear(a,e,x,inf)
         if(inf/=0)then;info=inf;return;end if
         ainv(:,j)=x
      end do
   end subroutine inverse_matrix

   function determinant_matrix(a) result(det)
      real(dp),intent(in)::a(:,:);real(dp)::det,m(size(a,1),size(a,2)),factor,tmp
      integer::n,i,j,k,p,sgn
      n=size(a,1);if(size(a,2)/=n)then;det=0.0_dp;return;end if
      m=a;sgn=1
      do k=1,n-1
         p=k;do i=k+1,n;if(abs(m(i,k))>abs(m(p,k)))p=i;end do
         if(abs(m(p,k))<=tiny(1.0_dp))then;det=0.0_dp;return;end if
         if(p/=k)then
            do j=k,n;tmp=m(k,j);m(k,j)=m(p,j);m(p,j)=tmp;end do;sgn=-sgn
         end if
         do i=k+1,n
            factor=m(i,k)/m(k,k)
            do j=k+1,n;m(i,j)=m(i,j)-factor*m(k,j);end do
         end do
      end do
      det=real(sgn,dp);do i=1,n;det=det*m(i,i);end do
   end function determinant_matrix

   function logdet_spd(a,info) result(v)
      real(dp),intent(in)::a(:,:);integer,intent(out),optional::info;real(dp)::v,l(size(a,1),size(a,2));integer::i,inf
      call cholesky_lower(a,l,inf);if(present(info))info=inf
      if(inf/=0)then;v=-huge(1.0_dp);return;end if
      v=0.0_dp;do i=1,size(a,1);v=v+2.0_dp*log(l(i,i));end do
   end function logdet_spd

   subroutine eigen_sym_jacobi(a,values,vectors,info,tol,max_iter)
      real(dp),intent(in)::a(:,:);real(dp),intent(out)::values(size(a,1)),vectors(size(a,1),size(a,1))
      integer,intent(out)::info;real(dp),intent(in),optional::tol;integer,intent(in),optional::max_iter
      real(dp)::b(size(a,1),size(a,1)),t,c,s,tau,app,aqq,apq,bip,biq,vip,viq,thr
      integer::n,i,p,q,it,mi
      n=size(a,1);thr=1e-12_dp;if(present(tol))thr=tol;mi=100*n*n;if(present(max_iter))mi=max_iter
      b=a;vectors=0.0_dp;do i=1,n;vectors(i,i)=1.0_dp;end do
      info=1
      do it=1,mi
         apq=0.0_dp;p=1;q=min(2,n)
         do i=1,n-1
            if(maxval(abs(b(i,i+1:n)))>abs(apq))then
               q=i+maxloc(abs(b(i,i+1:n)),dim=1);p=i;apq=b(p,q)
            end if
         end do
         if(abs(apq)<=thr)then;info=0;exit;end if
         app=b(p,p);aqq=b(q,q)
         tau=(aqq-app)/(2.0_dp*apq)
         if(tau>=0.0_dp)then;t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau));else;t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau));end if
         c=1.0_dp/sqrt(1.0_dp+t*t);s=t*c
         do i=1,n
            if(i/=p.and.i/=q)then
               bip=b(i,p);biq=b(i,q);b(i,p)=c*bip-s*biq;b(p,i)=b(i,p);b(i,q)=s*bip+c*biq;b(q,i)=b(i,q)
            end if
         end do
         b(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
         b(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq;b(p,q)=0.0_dp;b(q,p)=0.0_dp
         do i=1,n
            vip=vectors(i,p);viq=vectors(i,q);vectors(i,p)=c*vip-s*viq;vectors(i,q)=s*vip+c*viq
         end do
      end do
      do i=1,n;values(i)=b(i,i);end do
      call sort_eigen_desc(values,vectors)
   end subroutine eigen_sym_jacobi

   subroutine sort_eigen_desc(values,vectors)
      real(dp),intent(inout)::values(:),vectors(:,:);integer::i,j,k;real(dp)::tv,col(size(values))
      do i=1,size(values)-1
         k=i;do j=i+1,size(values);if(values(j)>values(k))k=j;end do
         if (k /= i) then
            tv = values(i)
            values(i) = values(k)
            values(k) = tv
            col = vectors(:,i)
            vectors(:,i) = vectors(:,k)
            vectors(:,k) = col
         end if
      end do
   end subroutine sort_eigen_desc

   function mahalanobis_sq(x,center,cov) result(d)
      real(dp),intent(in)::x(:,:),center(:),cov(:,:);real(dp)::d(size(x,1)),inv(size(cov,1),size(cov,2)),z(size(center))
      integer::i,info
      call inverse_matrix(cov,inv,info)
      if(info/=0)then;d=huge(1.0_dp);return;end if
      do i=1,size(x,1);z=x(i,:)-center;d(i)=dot_product(z,matmul(inv,z));end do
   end function mahalanobis_sq

   function pooled_covariance(x,groups) result(s)
      real(dp),intent(in)::x(:,:);integer,intent(in)::groups(:)
      real(dp)::s(size(x,2),size(x,2)),m(size(x,2)),z(size(x,2));integer::g,i,k,ng,df
      integer,allocatable::ug(:)
      ! local unique group extraction preserving sorted values
      ug=unique_groups(groups);s=0.0_dp;df=0
      do k=1,size(ug)
         g=ug(k);ng=count(groups==g);if(ng<=1)cycle;m=0.0_dp
         do i=1,size(groups);if(groups(i)==g)m=m+x(i,:);end do;m=m/real(ng,dp)
         do i=1,size(groups);if(groups(i)==g)then;z=x(i,:)-m;s=s+outer_vec(z,z);end if;end do
         df=df+ng-1
      end do
      if(df>0)s=s/real(df,dp)
   end function pooled_covariance

   function unique_groups(g) result(u)
      integer,intent(in)::g(:);integer,allocatable::u(:),tmp(:);integer::i,n,j,key
      tmp=g
      do i=2,size(tmp);key=tmp(i);j=i-1;do while(j>=1);if(tmp(j)<=key)exit;tmp(j+1)=tmp(j);j=j-1;end do;tmp(j+1)=key;end do
      n=0;do i=1,size(tmp);if(i==1.or.tmp(i)/=tmp(i-1))n=n+1;end do;allocate(u(n));n=0
      do i=1,size(tmp);if(i==1.or.tmp(i)/=tmp(i-1))then;n=n+1;u(n)=tmp(i);end if;end do
   end function unique_groups

   pure function outer_vec(x,y) result(a)
      real(dp),intent(in)::x(:),y(:);real(dp)::a(size(x),size(y));integer::i
      do i=1,size(x);a(i,:)=x(i)*y;end do
   end function outer_vec

   subroutine spd_inverse(a,ainv,info)
      real(dp),intent(in)::a(:,:);real(dp),intent(out)::ainv(size(a,1),size(a,2));integer,intent(out)::info
      real(dp)::l(size(a,1),size(a,2)),y(size(a,1)),x(size(a,1)),e(size(a,1));integer::i,j,n
      n=size(a,1);call cholesky_lower(a,l,info);if(info/=0)return;ainv=0.0_dp
      do j=1,n
         e=0.0_dp;e(j)=1.0_dp;y=0.0_dp
         do i=1,n;y(i)=(e(i)-dot_product(l(i,1:i-1),y(1:i-1)))/l(i,i);end do
         x=0.0_dp
         do i=n,1,-1;x(i)=(y(i)-dot_product(l(i+1:n,i),x(i+1:n)))/l(i,i);end do
         ainv(:,j)=x
      end do
   end subroutine spd_inverse

   pure function standard_covariance(x) result(s)
      real(dp),intent(in)::x(:,:);real(dp)::s(size(x,2),size(x,2));s=covariance_matrix(x)
   end function standard_covariance

   integer function matrix_rank(a,tol) result(r)
      real(dp),intent(in)::a(:,:);real(dp),intent(in),optional::tol
      real(dp)::m(size(a,1),size(a,2)),t,fac,tmp;integer::i,j,k,p,nr,nc
      nr=size(a,1);nc=size(a,2);m=a;t=1e-10_dp;if(present(tol))t=tol;r=0;i=1
      do j=1,nc
         if(i>nr)exit;p=i
         do k=i+1,nr;if(abs(m(k,j))>abs(m(p,j)))p=k;end do
         if(abs(m(p,j))<=t)cycle
         if(p/=i)then;do k=j,nc;tmp=m(i,k);m(i,k)=m(p,k);m(p,k)=tmp;end do;end if
         do k=i+1,nr;fac=m(k,j)/m(i,j);m(k,j:nc)=m(k,j:nc)-fac*m(i,j:nc);end do
         r=r+1;i=i+1
      end do
   end function matrix_rank

end module rfast_linalg
