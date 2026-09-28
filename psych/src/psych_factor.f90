module psych_factor
   use r_kinds, only : dp
   use r_linalg, only : symmetric_eigen, thin_svd
   use psych_stats, only : correlation_matrix, covariance_matrix, scale_columns, center_columns
   use psych_correlations, only : pseudoinverse, squared_multiple_correlations
   use psych_types, only : pca_result
   implicit none
   private
   public :: principal_components_from_matrix, principal_components_data, eigen_loadings
   public :: principal_axis_factor, factor_scores_regression, varimax_rotate, promax_rotate
   public :: target_rotation, kaiser_normalize, reflect_loadings, factor_to_clusters
   public :: procrustes_rotate, factor_residuals, factor_complexity

contains

   subroutine principal_components_from_matrix(r, n_factors, result)
      real(dp), intent(in) :: r(:, :) !! Correlation or covariance matrix.
      integer, intent(in) :: n_factors !! Number of leading components to return.
      type(pca_result), intent(out) :: result !! Eigenvalues and component loadings; scores are empty.
      real(dp), allocatable :: values(:), vectors(:, :)
      integer :: i, k
      call symmetric_eigen(r,values,vectors,result%info,descending=.true.)
      if (result%info /= 0) then
         allocate(result%values(0),result%loadings(0,0),result%scores(0,0))
         return
      end if
      k = min(max(1,n_factors),size(values))
      allocate(result%values(k),result%loadings(size(r,1),k),result%scores(0,0))
      result%values = values(1:k)
      do i = 1, k
         result%loadings(:,i) = vectors(:,i)*sqrt(max(0.0_dp,values(i)))
      end do
   end subroutine principal_components_from_matrix

   subroutine principal_components_data(x, n_factors, result, covariance_mode)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable data matrix.
      integer, intent(in) :: n_factors !! Number of leading components to return.
      type(pca_result), intent(out) :: result !! Eigenvalues, loadings, and standardized/centered component scores.
      logical, intent(in), optional :: covariance_mode !! Use covariance rather than correlation PCA when true.
      real(dp), allocatable :: r(:, :), transformed(:, :), means(:), sds(:), weights(:, :)
      logical :: use_covariance
      integer :: p, k, info
      p = size(x,2)
      use_covariance = .false.
      if (present(covariance_mode)) use_covariance = covariance_mode
      allocate(r(p,p),transformed(size(x,1),p),means(p),sds(p))
      if (use_covariance) then
         call covariance_matrix(x,r)
         call center_columns(x,transformed,means)
         sds = 1.0_dp
      else
         call correlation_matrix(x,r)
         call scale_columns(x,transformed,means,sds)
      end if
      call principal_components_from_matrix(r,n_factors,result)
      if (result%info /= 0) return
      k = size(result%values)
      call pseudoinverse(r,weights,info=info)
      if (info == 0) then
         weights = matmul(weights,result%loadings)
         if (allocated(result%scores)) deallocate(result%scores)
         allocate(result%scores(size(x,1),k))
         result%scores = matmul(transformed,weights)
      end if
   end subroutine principal_components_data

   pure subroutine eigen_loadings(values, vectors, loadings)
      real(dp), intent(in) :: values(:) !! Eigenvalues corresponding to columns of `vectors`.
      real(dp), intent(in) :: vectors(:, :) !! Eigenvector matrix.
      real(dp), intent(out) :: loadings(:, :) !! Principal-component loadings `vectors*sqrt(values)`.
      integer :: j
      do j = 1, min(size(values),size(vectors,2))
         loadings(:,j) = vectors(:,j)*sqrt(max(0.0_dp,values(j)))
      end do
   end subroutine eigen_loadings

   subroutine principal_axis_factor(r, n_factors, loadings, uniqueness, info, max_iter, tolerance)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix.
      integer, intent(in) :: n_factors !! Number of common factors.
      real(dp), intent(out) :: loadings(:, :) !! Principal-axis loading matrix.
      real(dp), intent(out) :: uniqueness(:) !! Uniqueness estimates `1-h2` after convergence.
      integer, intent(out) :: info !! Zero on successful eigensolutions.
      integer, intent(in), optional :: max_iter !! Maximum communality iterations; default 100.
      real(dp), intent(in), optional :: tolerance !! Maximum communality change for convergence; default `1e-7`.
      real(dp), allocatable :: work(:, :), values(:), vectors(:, :), h2(:), old_h2(:), smc(:)
      integer :: i, iter, p, k, limit, smc_info
      real(dp) :: tol
      p = size(r,1)
      k = min(max(1,n_factors),p)
      limit = 100
      if (present(max_iter)) limit = max(1,max_iter)
      tol = 1.0e-7_dp
      if (present(tolerance)) tol = max(0.0_dp,tolerance)
      allocate(work(p,p),h2(p),old_h2(p),smc(p))
      call squared_multiple_correlations(r,smc,smc_info)
      if (smc_info == 0) then
         h2 = smc
      else
         h2 = min(1.0_dp,max(0.0_dp,sum(abs(r),dim=2)-1.0_dp))
      end if
      do iter = 1, limit
         old_h2 = h2
         work = r
         do i = 1, p
            work(i,i) = h2(i)
         end do
         call symmetric_eigen(work,values,vectors,info,descending=.true.)
         if (info /= 0) return
         loadings = 0.0_dp
         do i = 1, k
            loadings(:,i) = vectors(:,i)*sqrt(max(0.0_dp,values(i)))
         end do
         h2 = min(1.0_dp,max(0.0_dp,sum(loadings*loadings,dim=2)))
         if (maxval(abs(h2-old_h2)) <= tol) exit
      end do
      uniqueness = max(0.0_dp,1.0_dp-h2)
   end subroutine principal_axis_factor

   subroutine factor_scores_regression(r, loadings, weights, info)
      real(dp), intent(in) :: r(:, :) !! Observed-variable correlation matrix.
      real(dp), intent(in) :: loadings(:, :) !! Factor structure/loading matrix.
      real(dp), intent(out) :: weights(:, :) !! Regression factor-score weights `R^+ L`.
      integer, intent(out) :: info !! Zero on successful pseudoinversion.
      real(dp), allocatable :: rinv(:, :)
      call pseudoinverse(r,rinv,info=info)
      if (info == 0) then
         weights = matmul(rinv,loadings)
      else
         weights = 0.0_dp
      end if
   end subroutine factor_scores_regression

   subroutine varimax_rotate(loadings, rotated, rotation, info, max_iter, tolerance, normalize)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix to rotate orthogonally.
      real(dp), intent(out) :: rotated(:, :) !! Varimax-rotated loading matrix.
      real(dp), intent(out) :: rotation(:, :) !! Orthogonal rotation matrix.
      integer, intent(out) :: info !! Zero when all pair rotations converge.
      integer, intent(in), optional :: max_iter !! Maximum sweeps over factor pairs; default 100.
      real(dp), intent(in), optional :: tolerance !! Rotation-angle tolerance; default `1e-8`.
      logical, intent(in), optional :: normalize !! Apply Kaiser row normalization when true; default false.
      real(dp), allocatable :: a(:, :), row_norm(:)
      real(dp) :: angle, c, s, u, v, aa, bb, cc, dd, change, tol
      integer :: i, j, row, p, k, iter, limit
      logical :: do_normalize
      p = size(loadings,1)
      k = size(loadings,2)
      limit = 100
      if (present(max_iter)) limit = max(1,max_iter)
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = max(0.0_dp,tolerance)
      do_normalize = .false.
      if (present(normalize)) do_normalize = normalize
      allocate(a(p,k),row_norm(p))
      a = loadings
      rotation = 0.0_dp
      do i = 1, k
         rotation(i,i) = 1.0_dp
      end do
      row_norm = 1.0_dp
      if (do_normalize) then
         do row = 1, p
            row_norm(row) = sqrt(sum(a(row,:)*a(row,:)))
            if (row_norm(row) > 0.0_dp) a(row,:) = a(row,:)/row_norm(row)
         end do
      end if
      info = 0
      do iter = 1, limit
         change = 0.0_dp
         do i = 1, k-1
            do j = i+1, k
               aa = 0.0_dp
               bb = 0.0_dp
               cc = 0.0_dp
               dd = 0.0_dp
               do row = 1, p
                  u = a(row,i)
                  v = a(row,j)
                  aa = aa+(u*u-v*v)
                  bb = bb+2.0_dp*u*v
                  cc = cc+(u*u-v*v)*(u*u-v*v)-4.0_dp*u*u*v*v
                  dd = dd+4.0_dp*u*v*(u*u-v*v)
               end do
               angle = 0.25_dp*atan2(dd-2.0_dp*aa*bb/real(p,dp),cc-(aa*aa-bb*bb)/real(p,dp))
               if (abs(angle) <= tol) cycle
               c = cos(angle)
               s = sin(angle)
               call rotate_columns(a,i,j,c,s)
               call rotate_columns(rotation,i,j,c,s)
               change = max(change,abs(angle))
            end do
         end do
         if (change <= tol) exit
      end do
      if (do_normalize) then
         do row = 1, p
            a(row,:) = a(row,:)*row_norm(row)
         end do
      end if
      rotated = a
   end subroutine varimax_rotate

   pure subroutine rotate_columns(a, i, j, c, s)
      real(dp), intent(inout) :: a(:, :) !! Matrix whose columns `i` and `j` are rotated in place.
      integer, intent(in) :: i !! First column index.
      integer, intent(in) :: j !! Second column index.
      real(dp), intent(in) :: c !! Cosine of the plane-rotation angle.
      real(dp), intent(in) :: s !! Sine of the plane-rotation angle.
      real(dp), allocatable :: first(:), second(:)
      allocate(first(size(a,1)),second(size(a,1)))
      first = a(:,i)
      second = a(:,j)
      a(:,i) = c*first+s*second
      a(:,j) = -s*first+c*second
   end subroutine rotate_columns

   subroutine promax_rotate(loadings, power, rotated, transformation, phi, info)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix to rotate.
      integer, intent(in), optional :: power !! Promax target power `m`; default four.
      real(dp), intent(out) :: rotated(:, :) !! Obliquely rotated loadings.
      real(dp), intent(out) :: transformation(:, :) !! Transformation from the original loading axes.
      real(dp), intent(out) :: phi(:, :) !! Factor correlation matrix.
      integer, intent(out) :: info !! Zero when pseudoinverse steps succeed.
      real(dp), allocatable :: orthogonal(:, :), rot(:, :), target(:, :), pinv(:, :), u(:, :), utu(:, :)
      real(dp), allocatable :: utu_inv(:, :), ui(:, :)
      integer :: m, local_info, j
      m = 4
      if (present(power)) m = max(1,power)
      allocate(orthogonal(size(loadings,1),size(loadings,2)),rot(size(loadings,2),size(loadings,2)))
      call varimax_rotate(loadings,orthogonal,rot,local_info)
      target = orthogonal*abs(orthogonal)**real(m-1,dp)
      call pseudoinverse(orthogonal,pinv,info=local_info)
      if (local_info /= 0) then
         info = local_info
         rotated = loadings
         transformation = 0.0_dp
         phi = 0.0_dp
         return
      end if
      u = matmul(pinv,target)
      utu = matmul(transpose(u),u)
      call pseudoinverse(utu,utu_inv,info=local_info)
      if (local_info /= 0) then
         info = local_info
         rotated = loadings
         transformation = 0.0_dp
         phi = 0.0_dp
         return
      end if
      do j = 1, size(u,2)
         u(:,j) = u(:,j)*sqrt(max(0.0_dp,utu_inv(j,j)))
      end do
      transformation = matmul(rot,u)
      rotated = matmul(loadings,transformation)
      call pseudoinverse(transformation,ui,info=info)
      if (info == 0) then
         phi = matmul(ui,transpose(ui))
         call covariance_to_correlation(phi)
      else
         phi = 0.0_dp
      end if
   end subroutine promax_rotate

   subroutine target_rotation(loadings, target, rotated, transformation, phi, info)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix to rotate.
      real(dp), intent(in) :: target(:, :) !! Target/key loading pattern with matching shape.
      real(dp), intent(out) :: rotated(:, :) !! Target-rotated loadings.
      real(dp), intent(out) :: transformation(:, :) !! Oblique transformation matrix.
      real(dp), intent(out) :: phi(:, :) !! Factor correlation matrix.
      integer, intent(out) :: info !! Zero on successful pseudoinverse steps.
      real(dp), allocatable :: pinv(:, :), u(:, :), utu(:, :), utu_inv(:, :), ui(:, :)
      integer :: j
      call pseudoinverse(loadings,pinv,info=info)
      if (info /= 0) return
      u = matmul(pinv,target)
      utu = matmul(transpose(u),u)
      call pseudoinverse(utu,utu_inv,info=info)
      if (info /= 0) return
      do j = 1, size(u,2)
         u(:,j) = u(:,j)*sqrt(max(0.0_dp,utu_inv(j,j)))
      end do
      transformation = u
      rotated = matmul(loadings,u)
      call pseudoinverse(u,ui,info=info)
      if (info == 0) then
         phi = matmul(ui,transpose(ui))
         call covariance_to_correlation(phi)
      end if
   end subroutine target_rotation

   pure subroutine covariance_to_correlation(a)
      real(dp), intent(inout) :: a(:, :) !! Symmetric covariance-like matrix standardized in place to unit diagonal.
      real(dp), allocatable :: sd(:)
      integer :: i, j, p
      p = size(a,1)
      allocate(sd(p))
      do i = 1, p
         sd(i) = sqrt(max(tiny(1.0_dp),a(i,i)))
      end do
      do j = 1, p
         do i = 1, p
            a(i,j) = a(i,j)/(sd(i)*sd(j))
         end do
      end do
   end subroutine covariance_to_correlation

   pure subroutine kaiser_normalize(loadings, normalized, row_norms)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix.
      real(dp), intent(out) :: normalized(:, :) !! Row-normalized loading matrix.
      real(dp), intent(out) :: row_norms(:) !! Communality square roots used for normalization.
      integer :: i
      do i = 1, size(loadings,1)
         row_norms(i) = sqrt(sum(loadings(i,:)*loadings(i,:)))
         if (row_norms(i) > 0.0_dp) then
            normalized(i,:) = loadings(i,:)/row_norms(i)
         else
            normalized(i,:) = loadings(i,:)
         end if
      end do
   end subroutine kaiser_normalize

   pure subroutine reflect_loadings(loadings, flipped, reflected)
      real(dp), intent(in) :: loadings(:, :) !! Loading matrix.
      logical, intent(in) :: flipped(:) !! True for factor columns whose signs should be reversed.
      real(dp), intent(out) :: reflected(:, :) !! Loading matrix after requested sign reflections.
      integer :: j
      reflected = loadings
      do j = 1, min(size(loadings,2),size(flipped))
         if (flipped(j)) reflected(:,j) = -reflected(:,j)
      end do
   end subroutine reflect_loadings

   pure subroutine factor_to_clusters(loadings, clusters)
      real(dp), intent(in) :: loadings(:, :) !! Item-by-factor loading matrix.
      integer, intent(out) :: clusters(:) !! Index of the largest absolute loading for each item.
      integer :: i
      do i = 1, min(size(loadings,1),size(clusters))
         clusters(i) = maxloc(abs(loadings(i,:)),dim=1)
      end do
   end subroutine factor_to_clusters

   subroutine procrustes_rotate(x, target, rotated, rotation, info)
      real(dp), intent(in) :: x(:, :) !! Loading/configuration matrix to align.
      real(dp), intent(in) :: target(:, :) !! Target configuration with matching row dimension.
      real(dp), intent(out) :: rotated(:, :) !! Orthogonally aligned version of `x`.
      real(dp), intent(out) :: rotation(:, :) !! Orthogonal Procrustes rotation matrix.
      integer, intent(out) :: info !! SVD status.
      real(dp), allocatable :: u(:, :), singular(:), vt(:, :), cross(:, :)
      cross = matmul(transpose(x),target)
      call thin_svd(cross,u,singular,vt,info)
      if (info /= 0) then
         rotation = 0.0_dp
         rotated = x
         return
      end if
      rotation = matmul(u,vt)
      rotated = matmul(x,rotation)
   end subroutine procrustes_rotate

   pure subroutine factor_residuals(r, loadings, residual)
      real(dp), intent(in) :: r(:, :) !! Observed correlation/covariance matrix.
      real(dp), intent(in) :: loadings(:, :) !! Factor loading matrix.
      real(dp), intent(out) :: residual(:, :) !! Off-diagonal residual matrix after reproducing `L L'`.
      integer :: i
      residual = r-matmul(loadings,transpose(loadings))
      do i = 1, min(size(residual,1),size(residual,2))
         residual(i,i) = 0.0_dp
      end do
   end subroutine factor_residuals

   pure subroutine factor_complexity(loadings, complexity)
      real(dp), intent(in) :: loadings(:, :) !! Item-by-factor loading matrix.
      real(dp), intent(out) :: complexity(:) !! Hofmann-style item complexity `(sum l^2)^2/sum l^4`.
      real(dp) :: s2, s4
      integer :: i
      do i = 1, size(loadings,1)
         s2 = sum(loadings(i,:)*loadings(i,:))
         s4 = sum(loadings(i,:)**4)
         if (s4 > 0.0_dp) then
            complexity(i) = s2*s2/s4
         else
            complexity(i) = 0.0_dp
         end if
      end do
   end subroutine factor_complexity

end module psych_factor
