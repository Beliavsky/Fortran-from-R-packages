! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on MarkowitzR, copyright 2014-2020 Steven E. Pav.
module markowitzr_portfolio
   use markowitzr_kinds, only: dp
   use markowitzr_types, only: theta_result, markowitz_result
   use markowitzr_linalg, only: invert_matrix, identity_matrix, symmetric_ivech
   use markowitzr_linalg, only: symmetric_vech, kronecker_product
   use markowitzr_linalg, only: duplication_matrix, lower_vector_indices
   use markowitzr_linalg, only: vech_index, matrix_rank, symmetrize_matrix
   use markowitzr_moments, only: theta_vcov, covariance_empirical
   use markowitzr_moments, only: moment_covariance_callback
   implicit none
   private

   integer, parameter, public :: weights_upstream = 0
   integer, parameter, public :: weights_all_columns = 1

   public :: mp_vcov, markowitz_weights

contains

   function mp_vcov(x, feat, fit_intercept, weights, jmat, gmat, &
                    covariance_method, hac_lags, weight_mode, &
                    covariance_callback) result(out)
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in), optional :: feat(:, :), weights(:)
      logical, intent(in), optional :: fit_intercept
      real(dp), intent(in), optional :: jmat(:, :), gmat(:, :)
      integer, intent(in), optional :: covariance_method, hac_lags, weight_mode
      procedure(moment_covariance_callback), optional :: covariance_callback
      type(markowitz_result) :: out
      type(theta_result) :: asymv
      real(dp), allocatable :: xy(:, :), theta(:, :), pj(:, :), pg(:, :)
      real(dp), allocatable :: ptheta(:, :), preh(:, :), d(:, :)
      real(dp), allocatable :: selected(:, :), h(:, :), mj(:, :), mg(:, :)
      real(dp), allocatable :: stacked(:, :)
      integer, allocatable :: lower_indices(:)
      logical :: use_intercept
      integer :: n, p, f, ff, pp, q, method, lags, mode, status
      integer :: i, a, b, k, rank_j, rank_stacked

      use_intercept = .true.
      if (present(fit_intercept)) use_intercept = fit_intercept
      method = covariance_empirical
      if (present(covariance_method)) method = covariance_method
      lags = -1
      if (present(hac_lags)) lags = hac_lags
      mode = weights_upstream
      if (present(weight_mode)) mode = weight_mode

      n = size(x,1)
      p = size(x,2)
      f = 0
      if (present(feat)) f = size(feat,2)
      ff = f+merge(1,0,use_intercept)
      out%p = p
      out%ff = ff

      if (n < 2 .or. p < 1) then
         out%status = 1
         out%message = 'x must contain at least two rows and one column'
         return
      end if
      if (present(feat)) then
         if (size(feat,1) /= n) then
            out%status = 2
            out%message = 'feat and x must have the same number of rows'
            return
         end if
      end if
      if (ff < 1) then
         out%status = 3
         out%message = 'an intercept or at least one feature is required'
         return
      end if
      if (present(weights)) then
         if (size(weights) /= n) then
            out%status = 4
            out%message = 'weights must have one value per observation'
            return
         end if
      end if
      if (mode /= weights_upstream .and. mode /= weights_all_columns) then
         out%status = 5
         out%message = 'unknown weight mode'
         return
      end if

      allocate(xy(n,f+p))
      if (present(feat)) xy(:,1:f) = feat
      xy(:,f+1:f+p) = x
      if (present(weights)) then
         select case (mode)
         case (weights_upstream)
            if (f > 0) xy(:,1:f) = xy(:,1:f)*spread(weights,2,f)
         case (weights_all_columns)
            xy = xy*spread(weights,2,f+p)
         end select
      end if

      if (present(covariance_callback)) then
         asymv = theta_vcov(xy,use_intercept,method,lags,covariance_callback)
      else
         asymv = theta_vcov(xy,use_intercept,method,lags)
      end if
      out%n = asymv%n
      if (asymv%status /= 0) then
         out%status = asymv%status
         out%message = asymv%message
         return
      end if

      pp = asymv%pp
      q = pp*(pp+1)/2
      theta = symmetric_ivech(asymv%mu,status)
      if (status /= 0) then
         out%status = 10
         out%message = 'invalid packed second moment'
         return
      end if

      if (present(jmat)) then
         if (size(jmat,2) /= p .or. size(jmat,1) < 1) then
            out%status = 11
            out%message = 'jmat must have p columns and at least one row'
            return
         end if
         mj = twidlize(jmat,ff,p)
         call projected_precision(theta,mj,pj,status)
      else
         allocate(pj(pp,pp))
         call invert_matrix(theta,pj,status)
      end if
      if (status /= 0) then
         out%status = 12
         out%message = 'constrained second moment is singular'
         return
      end if

      allocate(pg(pp,pp))
      pg = 0.0_dp
      if (present(gmat)) then
         if (size(gmat,2) /= p .or. size(gmat,1) < 1) then
            out%status = 13
            out%message = 'gmat must have p columns and at least one row'
            return
         end if
         if (present(jmat)) then
            rank_j = matrix_rank(jmat)
            allocate(stacked(size(jmat,1)+size(gmat,1),p))
            stacked(1:size(jmat,1),:) = jmat
            stacked(size(jmat,1)+1:,:) = gmat
            rank_stacked = matrix_rank(stacked)
            if (rank_stacked > rank_j) then
               out%status = 14
               out%message = 'rows of gmat must lie in the row space of jmat'
               return
            end if
         end if
         mg = twidlize(gmat,ff,p)
         call projected_precision(theta,mg,pg,status)
         if (status /= 0) then
            out%status = 15
            out%message = 'hedging second moment is singular'
            return
         end if
      end if

      allocate(ptheta(pp,pp))
      ptheta = pj-pg
      allocate(out%mu(q),out%covariance(q,q))
      out%mu = symmetric_vech(ptheta)

      preh = kronecker_product(pj,pj)-kronecker_product(pg,pg)
      d = duplication_matrix(pp)
      lower_indices = lower_vector_indices(pp)
      allocate(selected(q,pp*pp),h(q,q))
      do i = 1, q
         selected(i,:) = preh(lower_indices(i),:)
      end do
      h = -matmul(selected,d)
      out%covariance = matmul(h,matmul(asymv%covariance,transpose(h)))
      call symmetrize_matrix(out%covariance)

      allocate(out%w(p,ff),out%w_indices(p*ff),out%w_covariance(p*ff,p*ff))
      out%w = -ptheta(ff+1:ff+p,1:ff)
      k = 0
      do b = 1, ff
         do a = 1, p
            k = k+1
            out%w_indices(k) = vech_index(pp,ff+a,b)
         end do
      end do
      do b = 1, p*ff
         do a = 1, p*ff
            out%w_covariance(a,b) = &
               out%covariance(out%w_indices(a),out%w_indices(b))
         end do
      end do
      call symmetrize_matrix(out%w_covariance)
   end function mp_vcov

   function markowitz_weights(mean_returns, covariance, status) result(weights)
      real(dp), intent(in) :: mean_returns(:), covariance(:, :)
      integer, intent(out), optional :: status
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: inverse_covariance(:, :)
      integer :: local_status, p

      p = size(mean_returns)
      allocate(weights(p),inverse_covariance(p,p))
      if (any(shape(covariance) /= [p,p])) then
         weights = 0.0_dp
         if (present(status)) status = 1
         return
      end if
      call invert_matrix(covariance,inverse_covariance,local_status)
      if (local_status == 0) then
         weights = matmul(inverse_covariance,mean_returns)
      else
         weights = 0.0_dp
      end if
      if (present(status)) status = local_status
   end function markowitz_weights

   function twidlize(m, ff, p) result(twid)
      real(dp), intent(in) :: m(:, :)
      integer, intent(in) :: ff, p
      real(dp) :: twid(ff+size(m,1),ff+p)

      twid = 0.0_dp
      twid(1:ff,1:ff) = identity_matrix(ff)
      twid(ff+1:,ff+1:) = m
   end function twidlize

   subroutine projected_precision(theta, m, projected, status)
      real(dp), intent(in) :: theta(:, :), m(:, :)
      real(dp), allocatable, intent(out) :: projected(:, :)
      integer, intent(out) :: status
      real(dp), allocatable :: middle(:, :), inverse_middle(:, :)
      integer :: rows, pp

      rows = size(m,1)
      pp = size(theta,1)
      allocate(middle(rows,rows),inverse_middle(rows,rows),projected(pp,pp))
      middle = matmul(m,matmul(theta,transpose(m)))
      call invert_matrix(middle,inverse_middle,status)
      if (status /= 0) then
         projected = 0.0_dp
         return
      end if
      projected = matmul(transpose(m),matmul(inverse_middle,m))
      call symmetrize_matrix(projected)
   end subroutine projected_precision

end module markowitzr_portfolio
