module rfast2_multivariate
   use rfast_special, only : dp
   use rfast_arrays, only : colmeans, colvars
   use rfast_linalg, only : covariance_matrix, eigen_sym_jacobi, inverse_matrix
   use rfast2_types, only : pca_result, pcr_result
   implicit none
   private

   public :: pca, pcr, depth_mahala, leverage, discriminability, item_difficulty
   public :: covariance_distance

contains

   function pca(x,center,scale,k,vectors) result(res)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: center,scale,vectors
      integer, intent(in), optional :: k
      type(pca_result) :: res
      logical :: cen,scl,want_vec
      integer :: n,p,kk,j,info
      real(dp), allocatable :: z(:,:),cov(:,:),eval(:),evec(:,:)

      n = size(x,1)
      p = size(x,2)
      if (n < 2 .or. p < 1) then
         res%status = 1
         return
      end if
      cen = .true.
      scl = .true.
      want_vec = .false.
      if (present(center)) cen = center
      if (present(scale)) scl = scale
      if (present(vectors)) want_vec = vectors
      kk = p
      if (present(k)) kk = max(1,min(k,p))
      allocate(z(n,p),res%center(p),res%scale(p),cov(p,p),eval(p),evec(p,p))
      res%center = 0.0_dp
      if (cen) res%center = colmeans(x)
      res%scale = 1.0_dp
      if (scl) then
         res%scale = sqrt(max(0.0_dp,colvars(x)))
         where (res%scale <= tiny(1.0_dp)) res%scale = 1.0_dp
      end if
      do j = 1, p
         z(:,j) = (x(:,j)-res%center(j))/res%scale(j)
      end do
      cov = matmul(transpose(z),z)/real(n-1,dp)
      call eigen_sym_jacobi(cov,eval,evec,info)
      if (info /= 0) then
         res%status = info
         return
      end if
      call sort_eigen_local(eval,evec)
      allocate(res%values(kk))
      res%values = max(0.0_dp,eval(1:kk))
      if (want_vec) then
         allocate(res%vectors(p,kk))
         res%vectors = evec(:,1:kk)
      end if
   end function pca

   function pcr(y,x,k,xnew) result(res)
      real(dp), intent(in) :: y(:),x(:,:)
      integer, intent(in) :: k(:)
      real(dp), intent(in), optional :: xnew(:,:)
      type(pcr_result) :: res
      type(pca_result) :: pc
      real(dp), allocatable :: z(:,:),yc(:),coef_pc(:),cum_beta(:,:),zn(:,:)
      real(dp) :: my,den
      integer :: n,p,kmax,j,r,kk

      n = size(x,1)
      p = size(x,2)
      if (size(y) /= n .or. size(k) == 0 .or. any(k < 1) .or. maxval(k) > p) then
         res%status = 1
         return
      end if
      kmax = maxval(k)
      pc = pca(x,center=.true.,scale=.true.,k=kmax,vectors=.true.)
      if (pc%status /= 0) then
         res%status = pc%status
         return
      end if
      allocate(z(n,kmax),yc(n),coef_pc(kmax),cum_beta(p,kmax))
      z = 0.0_dp
      do j = 1, p
         z = z + spread((x(:,j)-pc%center(j))/pc%scale(j),2,kmax) * &
              spread(pc%vectors(j,:),1,n)
      end do
      my = sum(y)/real(n,dp)
      yc = y-my
      do j = 1, kmax
         den = sum(z(:,j)*z(:,j))
         if (den > tiny(1.0_dp)) then
            coef_pc(j) = dot_product(z(:,j),yc)/den
         else
            coef_pc(j) = 0.0_dp
         end if
      end do
      cum_beta = 0.0_dp
      do j = 1, kmax
         if (j == 1) then
            cum_beta(:,j) = pc%vectors(:,j)*coef_pc(j)
         else
            cum_beta(:,j) = cum_beta(:,j-1)+pc%vectors(:,j)*coef_pc(j)
         end if
      end do
      allocate(res%beta(p,size(k)),res%proportion(kmax),res%vectors(p,kmax))
      res%proportion = cumulative(pc%values)/real(p,dp)
      res%vectors = pc%vectors
      do r = 1, size(k)
         kk = k(r)
         res%beta(:,r) = cum_beta(:,kk)/pc%scale
      end do
      if (present(xnew)) then
         if (size(xnew,2) /= p) then
            res%status = 2
            return
         end if
         allocate(res%fitted(size(xnew,1),size(k)),zn(size(xnew,1),p))
         do j = 1, p
            zn(:,j) = xnew(:,j)-pc%center(j)
         end do
         res%fitted = my+matmul(zn,res%beta)
      end if
   end function pcr

   pure function cumulative(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
      integer :: i
      if (size(x) == 0) return
      y(1) = x(1)
      do i = 2, size(x)
         y(i) = y(i-1)+x(i)
      end do
   end function cumulative

   function depth_mahala(x,data) result(depth)
      real(dp), intent(in) :: x(:,:),data(:,:)
      real(dp) :: depth(size(x,1))
      real(dp), allocatable :: mu(:),cov(:,:),inv(:,:),z(:)
      integer :: i,info

      mu = colmeans(data)
      cov = covariance_matrix(data)
      allocate(inv(size(cov,1),size(cov,2)),z(size(mu)))
      call inverse_matrix(cov,inv,info)
      if (info /= 0) then
         depth = 0.0_dp
         return
      end if
      do i = 1, size(x,1)
         z = x(i,:)-mu
         depth(i) = 1.0_dp/(1.0_dp+dot_product(z,matmul(inv,z)))
      end do
   end function depth_mahala

   function leverage(x) result(h)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: h(size(x,1))
      real(dp) :: xtx(size(x,2),size(x,2)),inv(size(x,2),size(x,2))
      integer :: i,info
      xtx = matmul(transpose(x),x)
      call inverse_matrix(xtx,inv,info)
      if (info /= 0) then
         h = huge(1.0_dp)
         return
      end if
      do i = 1, size(x,1)
         h(i) = dot_product(x(i,:),matmul(inv,x(i,:)))
      end do
   end function leverage

   function item_difficulty(x) result(v)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: v(size(x,2))
      v = colmeans(x)
   end function item_difficulty

   function discriminability(x,frac) result(v)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: frac
      real(dp) :: v(size(x,2))
      real(dp), allocatable :: score(:),work(:,:)
      integer, allocatable :: ord(:)
      real(dp) :: f
      integer :: n,m,i,j,t

      n = size(x,1)
      f = 1.0_dp/3.0_dp
      if (present(frac)) f = frac
      m = max(1,min(n/2,int(real(n,dp)*f)))
      allocate(score(n),ord(n),work(n,size(x,2)))
      score = sum(x,dim=2)/real(size(x,2),dp)
      ord = [(i,i=1,n)]
      do i = 1, n-1
         t = i
         do j = i+1, n
            if (score(ord(j)) < score(ord(t))) t = j
         end do
         if (t /= i) then
            j = ord(i)
            ord(i) = ord(t)
            ord(t) = j
         end if
      end do
      work = x(ord,:)
      v = (sum(work(n-m+1:n,:),dim=1)-sum(work(1:m,:),dim=1))/real(m,dp)
   end function discriminability

   function covariance_distance(x,y) result(d)
      real(dp), intent(in) :: x(:,:),y(:,:)
      real(dp) :: d
      real(dp), allocatable :: a(:,:),b(:,:)
      a = covariance_matrix(x)
      b = covariance_matrix(y)
      d = sqrt(sum((a-b)**2))
   end function covariance_distance

   subroutine sort_eigen_local(values,vectors)
      real(dp), intent(inout) :: values(:),vectors(:,:)
      real(dp) :: tv,col(size(values))
      integer :: i,j,k
      do i = 1, size(values)-1
         k = i
         do j = i+1, size(values)
            if (values(j) > values(k)) k = j
         end do
         if (k /= i) then
            tv = values(i)
            values(i) = values(k)
            values(k) = tv
            col = vectors(:,i)
            vectors(:,i) = vectors(:,k)
            vectors(:,k) = col
         end if
      end do
   end subroutine sort_eigen_local

end module rfast2_multivariate
