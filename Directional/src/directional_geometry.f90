module directional_geometry
   use directional_kinds, only : dp, pi
   implicit none
   private
   public :: euclid, euclid_inv, rotation_matrix, haversine_dist, normalize_rows
contains
   pure function euclid(u) result(x)
      real(dp), intent(in) :: u(:,:)
      real(dp) :: x(size(u,1),3), lat, lon
      integer :: i
      do i=1,size(u,1)
         lat=pi*u(i,1)/180.0_dp; lon=pi*u(i,2)/180.0_dp
         x(i,:)=[cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat)]
      end do
   end function euclid

   pure function euclid_inv(x) result(u)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: u(size(x,1),2)
      integer :: i
      do i=1,size(x,1)
         u(i,:) = 180.0_dp/pi*[asin(max(-1.0_dp,min(1.0_dp,x(i,3)))),atan2(x(i,2),x(i,1))]
      end do
   end function euclid_inv

   pure function rotation_matrix(a,b) result(r)
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: r(size(a),size(a)), c(size(a)), aa, theta, nrm
      integer :: i,j,p
      real(dp) :: skew(size(a),size(a))
      p=size(a); aa=max(-1.0_dp,min(1.0_dp,dot_product(a,b)))
      c=a-b*aa; nrm=sqrt(sum(c*c))
      r=0.0_dp; do i=1,p; r(i,i)=1.0_dp; end do
      if (nrm < 1.0e-14_dp) then
         if (aa > 0.0_dp) return
         r=-r; return
      end if
      c=c/nrm; skew=0.0_dp
      do i=1,p; do j=1,p; skew(i,j)=b(i)*c(j)-c(i)*b(j); end do; end do
      theta=acos(aa)
      do i=1,p; do j=1,p
         r(i,j)=r(i,j)+sin(theta)*skew(i,j)+(cos(theta)-1.0_dp)*(b(i)*b(j)+c(i)*c(j))
      end do; end do
   end function rotation_matrix

   pure function haversine_dist(x) result(d)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: d(size(x,1),size(x,1)), a
      integer :: i,j
      d=0.0_dp
      do i=1,size(x,1)-1; do j=i+1,size(x,1)
         a=sin(0.5_dp*(x(i,1)-x(j,1)))**2 + cos(x(i,1))*cos(x(j,1))*sin(0.5_dp*(x(i,2)-x(j,2)))**2
         d(i,j)=2.0_dp*asin(sqrt(max(0.0_dp,min(1.0_dp,a)))); d(j,i)=d(i,j)
      end do; end do
   end function haversine_dist

   pure function normalize_rows(x) result(y)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: y(size(x,1),size(x,2)), nrm
      integer :: i
      do i=1,size(x,1); nrm=sqrt(sum(x(i,:)**2)); if(nrm>0)then;y(i,:)=x(i,:)/nrm;else;y(i,:)=x(i,:);end if; end do
   end function normalize_rows
end module directional_geometry
