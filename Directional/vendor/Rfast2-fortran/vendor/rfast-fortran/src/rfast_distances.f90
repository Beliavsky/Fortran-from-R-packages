module rfast_distances
   use rfast_special, only : dp, pi, nan_r
   implicit none
   private
   integer, parameter, public :: DIST_EUCLIDEAN=1, DIST_MANHATTAN=2, DIST_MAXIMUM=3
   integer, parameter, public :: DIST_MINKOWSKI=4, DIST_CANBERRA=5, DIST_COSINE=6, DIST_HAVERSINE=7
   public :: vector_distance, dist_matrix, dista_matrix, total_dist, vecdist_upper
   public :: distance_covariance, distance_correlation, distance_variance, partial_distance_correlation
   public :: energy_distance, squareform_from_vector

contains

   pure real(dp) function vector_distance(x,y,method,p,square) result(d)
      real(dp),intent(in)::x(:),y(:)
      integer,intent(in),optional::method
      real(dp),intent(in),optional::p
      logical,intent(in),optional::square
      integer::m,i;real(dp)::pp,den,dlon,dlat,a,c;logical::sq
      m=DIST_EUCLIDEAN;if(present(method))m=method;pp=2.0_dp;if(present(p))pp=p;sq=.false.;if(present(square))sq=square
      select case(m)
      case(DIST_EUCLIDEAN)
         d=sum((x-y)**2);if(.not.sq)d=sqrt(d)
      case(DIST_MANHATTAN);d=sum(abs(x-y))
      case(DIST_MAXIMUM);d=maxval(abs(x-y))
      case(DIST_MINKOWSKI);d=sum(abs(x-y)**pp)**(1.0_dp/pp)
      case(DIST_CANBERRA)
         d=0.0_dp;do i=1,size(x);den=abs(x(i))+abs(y(i));if(den>0)d=d+abs(x(i)-y(i))/den;end do
      case(DIST_COSINE)
         den=sqrt(sum(x*x)*sum(y*y));if(den>0)then;d=1.0_dp-dot_product(x,y)/den;else;d=0.0_dp;end if
      case(DIST_HAVERSINE)
         if(size(x)<2)then;d=nan_r();return;end if
         dlat=(y(1)-x(1))*pi/180.0_dp;dlon=(y(2)-x(2))*pi/180.0_dp
         a=sin(0.5_dp*dlat)**2+cos(x(1)*pi/180.0_dp)*cos(y(1)*pi/180.0_dp)*sin(0.5_dp*dlon)**2
         c=2.0_dp*asin(min(1.0_dp,sqrt(max(0.0_dp,a))));d=6371.0088_dp*c
      case default;d=nan_r()
      end select
   end function vector_distance

   function dist_matrix(x,method,p,square) result(d)
      real(dp),intent(in)::x(:,:);integer,intent(in),optional::method;real(dp),intent(in),optional::p
      logical,intent(in),optional::square;real(dp)::d(size(x,1),size(x,1));integer::i,j,m;real(dp)::pp;logical::sq
      m=DIST_EUCLIDEAN;if(present(method))m=method;pp=2.0_dp;if(present(p))pp=p;sq=.false.;if(present(square))sq=square
      d=0.0_dp
      do j=1,size(x,1)-1;do i=j+1,size(x,1)
         d(i,j)=vector_distance(x(i,:),x(j,:),m,pp,sq);d(j,i)=d(i,j)
      end do;end do
   end function dist_matrix

   function dista_matrix(xnew,x,method,p,square) result(d)
      real(dp),intent(in)::xnew(:,:),x(:,:);integer,intent(in),optional::method;real(dp),intent(in),optional::p
      logical,intent(in),optional::square;real(dp)::d(size(xnew,1),size(x,1));integer::i,j,m;real(dp)::pp;logical::sq
      m=DIST_EUCLIDEAN;if(present(method))m=method;pp=2.0_dp;if(present(p))pp=p;sq=.false.;if(present(square))sq=square
      do j=1,size(x,1);do i=1,size(xnew,1);d(i,j)=vector_distance(xnew(i,:),x(j,:),m,pp,sq);end do;end do
   end function dista_matrix

   function total_dist(x,method,p,square) result(s)
      real(dp),intent(in)::x(:,:);integer,intent(in),optional::method;real(dp),intent(in),optional::p
      logical,intent(in),optional::square;real(dp)::s;integer::i,j,m;real(dp)::pp;logical::sq
      m=DIST_EUCLIDEAN;if(present(method))m=method;pp=2.0_dp;if(present(p))pp=p;sq=.false.;if(present(square))sq=square;s=0.0_dp
      do j=1,size(x,1)-1;do i=j+1,size(x,1);s=s+vector_distance(x(i,:),x(j,:),m,pp,sq);end do;end do
   end function total_dist

   function vecdist_upper(x,method,p,square) result(v)
      real(dp),intent(in)::x(:,:);integer,intent(in),optional::method;real(dp),intent(in),optional::p
      logical,intent(in),optional::square;real(dp),allocatable::v(:);integer::i,j,k,m,n;real(dp)::pp;logical::sq
      n=size(x,1);allocate(v(n*(n-1)/2));m=DIST_EUCLIDEAN;if(present(method))m=method
      pp=2.0_dp;if(present(p))pp=p;sq=.false.;if(present(square))sq=square;k=0
      do j=1,n-1;do i=j+1,n;k=k+1;v(k)=vector_distance(x(i,:),x(j,:),m,pp,sq);end do;end do
   end function vecdist_upper

   pure subroutine double_center(a,bc,bc_a)
      real(dp),intent(in)::a(:,:);logical,intent(in)::bc
      real(dp),intent(out)::bc_a(size(a,1),size(a,2))
      real(dp)::rm(size(a,1)),cm(size(a,2)),gm;integer::n,i,j
      n=size(a,1)
      if(.not.bc)then
         rm=sum(a,dim=2)/real(n,dp);cm=sum(a,dim=1)/real(n,dp);gm=sum(a)/real(n*n,dp)
         do j=1,n;do i=1,n;bc_a(i,j)=a(i,j)-rm(i)-cm(j)+gm;end do;end do
      else
         rm=sum(a,dim=2)/real(n-2,dp);cm=sum(a,dim=1)/real(n-2,dp);gm=sum(a)/real((n-1)*(n-2),dp)
         do j=1,n;do i=1,n
            if(i==j)then;bc_a(i,j)=0.0_dp;else;bc_a(i,j)=a(i,j)-rm(i)-cm(j)+gm;end if
         end do;end do
      end if
   end subroutine double_center

   function distance_covariance(x,y,bc) result(v)
      real(dp),intent(in)::x(:,:),y(:,:);logical,intent(in),optional::bc;real(dp)::v
      real(dp)::a(size(x,1),size(x,1)),b(size(y,1),size(y,1)),ac(size(x,1),size(x,1)),bcen(size(x,1),size(x,1))
      logical::biasc;integer::n
      biasc=.false.;if(present(bc))biasc=bc;n=size(x,1);a=dist_matrix(x);b=dist_matrix(y)
      call double_center(a,biasc,ac);call double_center(b,biasc,bcen)
      if(biasc)then;v=sum(ac*bcen)/real(n*(n-3),dp);else;v=sum(ac*bcen)/real(n*n,dp);end if
      v=sqrt(max(0.0_dp,v))
   end function distance_covariance

   function distance_variance(x,bc) result(v)
      real(dp),intent(in)::x(:,:);logical,intent(in),optional::bc;real(dp)::v
      v=distance_covariance(x,x,bc)
   end function distance_variance

   function distance_correlation(x,y,bc) result(r)
      real(dp),intent(in)::x(:,:),y(:,:);logical,intent(in),optional::bc;real(dp)::r,dc,dvx,dvy
      dc=distance_covariance(x,y,bc);dvx=distance_variance(x,bc);dvy=distance_variance(y,bc)
      if(dvx>0.and.dvy>0)then;r=dc/sqrt(dvx*dvy);else;r=0.0_dp;end if
   end function distance_correlation

   function partial_distance_correlation(x,y,z) result(r)
      real(dp),intent(in)::x(:,:),y(:,:),z(:,:);real(dp)::r,a1,a2,a3,den
      a1=distance_correlation(x,y,.true.);a2=distance_correlation(x,z,.true.);a3=distance_correlation(y,z,.true.)
      den=sqrt(max(0.0_dp,1-a2*a2))*sqrt(max(0.0_dp,1-a3*a3));if(den>0)then;r=(a1-a2*a3)/den;else;r=0;end if
   end function partial_distance_correlation

   function energy_distance(x,y) result(e)
      real(dp),intent(in)::x(:,:),y(:,:)
      real(dp)::e,mij,mii,mjj,n1,n2
      integer::i,j
      n1=real(size(x,1),dp);n2=real(size(y,1),dp)
      if(size(x,2)/=size(y,2).or.n1<=0.0_dp.or.n2<=0.0_dp)then
         e=nan_r();return
      end if
      mij=0.0_dp
      do j=1,size(y,1)
         do i=1,size(x,1)
            mij=mij+vector_distance(x(i,:),y(j,:),DIST_EUCLIDEAN)
         end do
      end do
      mii=total_dist(x,DIST_EUCLIDEAN)
      mjj=total_dist(y,DIST_EUCLIDEAN)
      ! Matches Rfast::edist: the usual two-sample energy statistic
      ! multiplied by n1*n2/(n1+n2).  total_dist stores each within-sample
      ! unordered pair once, hence the factors of two below.
      e=(2.0_dp*mij/(n1*n2)-2.0_dp*mii/(n1*n1)-2.0_dp*mjj/(n2*n2))*n1*n2/(n1+n2)
   end function energy_distance

   function squareform_from_vector(v,n) result(a)
      real(dp),intent(in)::v(:);integer,intent(in)::n;real(dp)::a(n,n);integer::i,j,k
      a=0.0_dp;k=0
      do j=1,n-1;do i=j+1,n;k=k+1;if(k<=size(v))then;a(i,j)=v(k);a(j,i)=v(k);end if;end do;end do
   end function squareform_from_vector

end module rfast_distances
