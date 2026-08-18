module misc_tools_matrix
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use misc_tools_kinds, only : dp
   implicit none
   private

   public :: insert_col, insert_row, sym_matrix
   public :: triang, vecli, vecli2m, veclipos
   public :: determinant, symmetric_eigenvalues
   public :: is_semidefinite, semidefiniteness
   public :: quasiconcavity, quasiconvexity

contains

   subroutine insert_col(m,pos,v,out)
      real(dp), intent(in) :: m(:,:),v(:)
      integer, intent(in) :: pos
      real(dp), allocatable, intent(out) :: out(:,:)
      integer :: nr,nc

      nr = size(m,1)
      nc = size(m,2)
      if (pos < 1 .or. pos > nc+1) error stop "insert_col: invalid position"
      if (size(v) /= 1 .and. size(v) /= nr) &
         error stop "insert_col: v must be scalar or have nrow(m) values"

      allocate(out(nr,nc+1))
      if (pos > 1) out(:,1:pos-1) = m(:,1:pos-1)
      if (size(v) == 1) then
         out(:,pos) = v(1)
      else
         out(:,pos) = v
      end if
      if (pos <= nc) out(:,pos+1:nc+1) = m(:,pos:nc)
   end subroutine insert_col

   subroutine insert_row(m,pos,v,out)
      real(dp), intent(in) :: m(:,:),v(:)
      integer, intent(in) :: pos
      real(dp), allocatable, intent(out) :: out(:,:)
      integer :: nr,nc

      nr = size(m,1)
      nc = size(m,2)
      if (pos < 1 .or. pos > nr+1) error stop "insert_row: invalid position"
      if (size(v) /= 1 .and. size(v) /= nc) &
         error stop "insert_row: v must be scalar or have ncol(m) values"

      allocate(out(nr+1,nc))
      if (pos > 1) out(1:pos-1,:) = m(1:pos-1,:)
      if (size(v) == 1) then
         out(pos,:) = v(1)
      else
         out(pos,:) = v
      end if
      if (pos <= nr) out(pos+1:nr+1,:) = m(pos:nr,:)
   end subroutine insert_row

   subroutine sym_matrix(data,out,nrow,byrow,upper,status)
      real(dp), intent(in) :: data(:)
      real(dp), allocatable, intent(out) :: out(:,:)
      integer, intent(in), optional :: nrow
      logical, intent(in), optional :: byrow,upper
      integer, intent(out), optional :: status
      integer :: n,need,k,i,j,idx
      logical :: br,up

      if (present(status)) status = 0
      if (size(data) == 0) then
         allocate(out(0,0))
         if (present(status)) status = 1
         return
      end if

      br = .false.
      up = .false.
      if (present(byrow)) br = byrow
      if (present(upper)) up = upper

      if (present(nrow)) then
         n = nrow
      else
         n = ceiling(-0.5_dp+sqrt(0.25_dp+2.0_dp*real(size(data),dp)) - &
                     sqrt(epsilon(1.0_dp)))
      end if
      if (n < 1) then
         allocate(out(0,0))
         if (present(status)) status = 2
         return
      end if

      need = n*(n+1)/2
      allocate(out(n,n))
      out = 0.0_dp
      k = 0

      ! R code fills upper triangle iff byrow != upper; otherwise lower.
      if (br .neqv. up) then
         ! Assignment to upper.tri in R column-major order.
         do j = 1, n
            do i = 1, j
               k = k+1
               idx = mod(k-1,size(data))+1
               if (k <= need) out(i,j) = data(idx)
            end do
         end do
         do j = 1, n
            do i = j+1, n
               out(i,j) = out(j,i)
            end do
         end do
      else
         ! Assignment to lower.tri in R column-major order.
         do j = 1, n
            do i = j, n
               k = k+1
               idx = mod(k-1,size(data))+1
               if (k <= need) out(i,j) = data(idx)
            end do
         end do
         do j = 1, n
            do i = 1, j-1
               out(i,j) = out(j,i)
            end do
         end do
      end if
   end subroutine sym_matrix

   pure integer function veclipos(i,j,n) result(pos)
      integer, intent(in) :: i,j,n
      integer :: mn,mx
      mn = min(i,j)
      mx = max(i,j)
      pos = n*(n-1)/2 - ((n-mn)*(n-mn+1))/2 + mx
   end function veclipos

   subroutine vecli(m,v)
      real(dp), intent(in) :: m(:,:)
      real(dp), allocatable, intent(out) :: v(:)
      integer :: n,i,j
      if (size(m,1) /= size(m,2)) error stop "vecli: matrix must be square"
      n = size(m,1)
      allocate(v(n*(n+1)/2))
      do i = 1, n
         do j = i, n
            v(veclipos(i,j,n)) = m(i,j)
         end do
      end do
   end subroutine vecli

   subroutine vecli2m(v,m,status)
      real(dp), intent(in) :: v(:)
      real(dp), allocatable, intent(out) :: m(:,:)
      integer, intent(out), optional :: status
      integer :: n,i,j
      real(dp) :: nr

      if (present(status)) status = 0
      nr = -0.5_dp+sqrt(0.25_dp+2.0_dp*real(size(v),dp))
      n = nint(nr)
      if (n*(n+1)/2 /= size(v)) then
         allocate(m(0,0))
         if (present(status)) status = 1
         return
      end if
      allocate(m(n,n))
      do i = 1, n
         do j = 1, n
            m(i,j) = v(veclipos(i,j,n))
         end do
      end do
   end subroutine vecli2m

   subroutine triang(v,n,m)
      real(dp), intent(in) :: v(:)
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: m(:,:)
      integer :: i,j,pos
      allocate(m(n,n))
      m = 0.0_dp
      do i = 1, n
         do j = i, n
            pos = veclipos(i,j,n)
            if (pos <= size(v)) m(i,j) = v(pos)
         end do
      end do
   end subroutine triang

   real(dp) function determinant(a,status) result(det)
      real(dp), intent(in) :: a(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: b(:,:)
      real(dp) :: pivot,tmpv,scale
      integer :: n,i,j,k,p,sign

      if (present(status)) status = 0
      if (size(a,1) /= size(a,2)) then
         det = 0.0_dp
         if (present(status)) status = 1
         return
      end if
      n = size(a,1)
      if (n == 0) then
         det = 1.0_dp
         return
      end if

      b = a
      sign = 1
      det = 1.0_dp
      scale = max(1.0_dp,maxval(abs(b)))

      do k = 1, n
         p = k
         do i = k+1, n
            if (abs(b(i,k)) > abs(b(p,k))) p = i
         end do
         pivot = b(p,k)
         if (abs(pivot) <= epsilon(1.0_dp)*scale) then
            det = 0.0_dp
            return
         end if
         if (p /= k) then
            do j = 1, n
               tmpv = b(k,j)
               b(k,j) = b(p,j)
               b(p,j) = tmpv
            end do
            sign = -sign
         end if
         pivot = b(k,k)
         det = det*pivot
         do i = k+1, n
            tmpv = b(i,k)/pivot
            if (abs(tmpv) <= tiny(1.0_dp)) cycle
            do j = k+1, n
               b(i,j) = b(i,j)-tmpv*b(k,j)
            end do
         end do
      end do
      det = real(sign,dp)*det
   end function determinant

   subroutine symmetric_eigenvalues(a,eig,status,tol,max_sweeps)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: eig(:)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: max_sweeps
      real(dp), allocatable :: b(:,:)
      real(dp) :: eps,off,app,aqq,apq,t,c,s,tau,bpk,bqk
      integer :: n,p,q,k,sweep,imax

      if (present(status)) status = 0
      if (size(a,1) /= size(a,2)) then
         allocate(eig(0))
         if (present(status)) status = 1
         return
      end if
      n = size(a,1)
      allocate(eig(n))
      if (n == 0) return

      eps = 100.0_dp*epsilon(1.0_dp)
      imax = max(30,10*n*n)
      if (present(tol)) eps = tol
      if (present(max_sweeps)) imax = max_sweeps
      b = 0.5_dp*(a+transpose(a))

      do sweep = 1, imax
         off = 0.0_dp
         do q = 2, n
            do p = 1, q-1
               off = max(off,abs(b(p,q)))
            end do
         end do
         if (off <= eps*max(1.0_dp,maxval(abs(b)))) exit

         do q = 2, n
            do p = 1, q-1
               apq = b(p,q)
               if (abs(apq) <= eps*max(1.0_dp,abs(b(p,p))+abs(b(q,q)))) cycle
               app = b(p,p)
               aqq = b(q,q)
               tau = (aqq-app)/(2.0_dp*apq)
               if (tau >= 0.0_dp) then
                  t = 1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
               else
                  t = -1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
               end if
               c = 1.0_dp/sqrt(1.0_dp+t*t)
               s = t*c

               b(p,p) = app-t*apq
               b(q,q) = aqq+t*apq
               b(p,q) = 0.0_dp
               b(q,p) = 0.0_dp
               do k = 1, n
                  if (k == p .or. k == q) cycle
                  bpk = b(p,k)
                  bqk = b(q,k)
                  b(p,k) = c*bpk-s*bqk
                  b(k,p) = b(p,k)
                  b(q,k) = s*bpk+c*bqk
                  b(k,q) = b(q,k)
               end do
            end do
         end do
      end do
      do k = 1, n
         eig(k) = b(k,k)
      end do
      if (sweep > imax .and. present(status)) status = 2
   end subroutine symmetric_eigenvalues

   logical function symmetric_enough(a,tol) result(ok)
      real(dp), intent(in) :: a(:,:),tol
      integer :: i,j,n
      if (size(a,1) /= size(a,2)) then
         ok = .false.
         return
      end if
      n = size(a,1)
      ok = .true.
      do j = 1, n
         do i = j+1, n
            if (abs(a(i,j)-a(j,i)) > tol*max(1.0_dp,abs(a(i,j)),abs(a(j,i)))) then
               ok = .false.
               return
            end if
         end do
      end do
   end function symmetric_enough

   logical function is_semidefinite(a,positive,tol,method,status) result(ok)
      real(dp), intent(in) :: a(:,:)
      logical, intent(in), optional :: positive
      real(dp), intent(in), optional :: tol
      character(len=*), intent(in), optional :: method
      integer, intent(out), optional :: status
      real(dp), allocatable :: m(:,:),eig(:),sub(:,:)
      real(dp) :: eps,det
      logical :: pos
      character(len=8) :: meth
      integer :: n,mask,k,i,j,ii,jj,istat

      if (present(status)) status = 0
      pos = .true.
      eps = 100.0_dp*epsilon(1.0_dp)
      if (present(positive)) pos = positive
      if (present(tol)) eps = tol

      if (size(a,1) /= size(a,2)) then
         ok = .false.
         if (present(status)) status = 1
         return
      end if
      if (.not. symmetric_enough(a,1000.0_dp*eps)) then
         ok = .false.
         if (present(status)) status = 2
         return
      end if
      n = size(a,1)
      m = 0.5_dp*(a+transpose(a))
      if (.not. pos) m = -m

      if (present(method)) then
         meth = adjustl(method)
      else if (n < 13) then
         meth = "det"
      else
         meth = "eigen"
      end if

      if (trim(meth) == "eigen") then
         call symmetric_eigenvalues(m,eig,istat,eps)
         if (istat /= 0) then
            ok = .false.
            if (present(status)) status = 3
            return
         end if
         ok = minval(eig) >= -eps
         return
      else if (trim(meth) /= "det") then
         ok = .false.
         if (present(status)) status = 4
         return
      end if

      if (n >= bit_size(mask)-1) then
         ok = .false.
         if (present(status)) status = 5
         return
      end if

      ok = .true.
      do mask = 1, 2**n-1
         k = popcnt(mask)
         allocate(sub(k,k))
         ii = 0
         do i = 1, n
            if (.not. btest(mask,i-1)) cycle
            ii = ii+1
            jj = 0
            do j = 1, n
               if (.not. btest(mask,j-1)) cycle
               jj = jj+1
               sub(ii,jj) = m(i,j)
            end do
         end do
         det = determinant(sub)
         deallocate(sub)
         if (det < -eps) then
            ok = .false.
            return
         end if
      end do
   end function is_semidefinite

   logical function semidefiniteness(a,positive,tol,method,status) result(ok)
      real(dp), intent(in) :: a(:,:)
      logical, intent(in), optional :: positive
      real(dp), intent(in), optional :: tol
      character(len=*), intent(in), optional :: method
      integer, intent(out), optional :: status
      ok = is_semidefinite(a,positive,tol,method,status)
   end function semidefiniteness

   logical function quasiconcavity(a,tol,status) result(ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(out), optional :: status
      real(dp) :: eps,det
      integer :: n,i

      eps = epsilon(1.0_dp)
      if (present(tol)) eps = tol
      if (present(status)) status = 0
      if (size(a,1) /= size(a,2) .or. size(a,1) < 2) then
         ok = .false.
         if (present(status)) status = 1
         return
      end if
      if (abs(a(1,1)) > eps) then
         ok = .false.
         if (present(status)) status = 2
         return
      end if

      n = size(a,1)
      ok = .true.
      do i = 2, n
         det = determinant(a(1:i,1:i))
         if (det*real((-1)**i,dp) > eps) then
            ok = .false.
            return
         end if
      end do
   end function quasiconcavity

   logical function quasiconvexity(a,tol,status) result(ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(out), optional :: status
      real(dp) :: eps,det
      integer :: n,i

      eps = epsilon(1.0_dp)
      if (present(tol)) eps = tol
      if (present(status)) status = 0
      if (size(a,1) /= size(a,2) .or. size(a,1) < 2) then
         ok = .false.
         if (present(status)) status = 1
         return
      end if
      if (abs(a(1,1)) > eps) then
         ok = .false.
         if (present(status)) status = 2
         return
      end if

      n = size(a,1)
      ok = .true.
      do i = 2, n
         det = determinant(a(1:i,1:i))
         if (det > eps) then
            ok = .false.
            return
         end if
      end do
   end function quasiconvexity

end module misc_tools_matrix
