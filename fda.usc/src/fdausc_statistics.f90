module fdausc_statistics
    use r_kinds, only : dp
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    public :: adot_matrix_vector, pcvm_statistic_value, rp_projection_statistics
    public :: fdr_reject, fdr_adjusted_pvalue
    public :: pred_mse, pred_rmse, pred_mae, prediction_measures
    public :: confusion_accuracy, confusion_kappa, confusion_iou, confusion_iou_by_class
    public :: confusion_recall, confusion_precision, confusion_f1, class_weights
    public :: functional_mean, functional_variance, center_curves, difference_derivative

contains

    pure subroutine adot_matrix_vector(n, inprod, adot_vec)
        integer, intent(in) :: n !! Number of residual observations defining the symmetric Adot matrix.
        real(dp), intent(in) :: inprod(:) !! Packed lower-triangle inner products in column-within-row indexing i(i-1)/2+j.
        real(dp), intent(out) :: adot_vec(:) !! Packed Adot values: common diagonal first, then the strict lower triangle.
        real(dp) :: num
        real(dp) :: den
        real(dp) :: quo
        real(dp) :: aux1
        real(dp) :: aux2
        real(dp) :: sumr
        integer :: i
        integer :: j
        integer :: r
        integer :: ij
        integer :: ir
        integer :: ii
        integer :: jr
        integer :: jj
        integer :: rj
        integer :: rr
        integer :: auxi
        integer :: auxj
        integer :: auxr

        adot_vec = 0.0_dp
        if (size(adot_vec) < 1 .or. n < 1) return
        adot_vec(1) = pi*real(n + 1, dp)
        do i = 2, n
            do j = 1, i - 1
                sumr = 0.0_dp
                do r = 1, n
                    if (i == r .or. j == r) then
                        sumr = sumr + pi
                    else
                        auxi = i*(i - 1)/2
                        auxj = j*(j - 1)/2
                        auxr = r*(r - 1)/2
                        ij = auxi + j
                        ii = auxi + i
                        jj = auxj + j
                        rr = auxr + r
                        if (i > r) then
                            ir = auxi + r
                        else
                            ir = auxr + i
                        end if
                        if (r > j) then
                            rj = auxr + j
                        else
                            rj = auxj + r
                        end if
                        jr = rj
                        num = inprod(ij) - inprod(ir) - inprod(rj) + inprod(rr)
                        aux1 = sqrt(max(0.0_dp, inprod(ii) - 2.0_dp*inprod(ir) + inprod(rr)))
                        aux2 = sqrt(max(0.0_dp, inprod(jj) - 2.0_dp*inprod(jr) + inprod(rr)))
                        den = aux1*aux2
                        if (den <= tiny(1.0_dp)) then
                            quo = 0.0_dp
                        else
                            quo = num/den
                        end if
                        quo = max(-1.0_dp, min(1.0_dp, quo))
                        sumr = sumr + abs(pi - acos(quo))
                    end if
                end do
                adot_vec(1 + (i - 1)*(i - 2)/2 + j) = sumr
            end do
        end do
    end subroutine adot_matrix_vector

    pure real(dp) function pcvm_statistic_value(adot_vec, residuals) result(statistic)
        real(dp), intent(in) :: adot_vec(:) !! Packed Adot representation from adot_matrix_vector.
        real(dp), intent(in) :: residuals(:) !! Regression residual vector entering the projected Cramer-von Mises quadratic form.
        real(dp) :: sums
        integer :: i
        integer :: j
        integer :: n

        n = size(residuals)
        sums = 0.0_dp
        do i = 2, n
            do j = 1, i - 1
                sums = sums + residuals(i)*adot_vec(1 + (i - 1)*(i - 2)/2 + j)*residuals(j)
            end do
        end do
        statistic = adot_vec(1)*dot_product(residuals, residuals) + 2.0_dp*sums
    end function pcvm_statistic_value

    pure subroutine rp_projection_statistics(proj_order, residuals, statistics)
        integer, intent(in) :: proj_order(:, :) !! One-based observation indices in projected order, one projection per column.
        real(dp), intent(in) :: residuals(:, :) !! Projection-specific residuals with shape [n_projection,n_observation].
        real(dp), intent(out) :: statistics(:, :) !! Per-projection [Cramer-von Mises, Kolmogorov-Smirnov] statistics.
        real(dp), allocatable :: cumulative(:)
        integer :: i
        integer :: j
        integer :: idx
        integer :: n

        n = size(proj_order, 1)
        allocate(cumulative(n))
        do i = 1, size(proj_order, 2)
            cumulative = 0.0_dp
            do j = 1, n
                idx = proj_order(j, i)
                if (j == 1) then
                    cumulative(j) = residuals(i, idx)
                else
                    cumulative(j) = cumulative(j - 1) + residuals(i, idx)
                end if
            end do
            statistics(i, 1) = sum(cumulative*cumulative)/real(n*n, dp)
            statistics(i, 2) = maxval(abs(cumulative))/sqrt(real(n, dp))
        end do
    end subroutine rp_projection_statistics

    pure logical function fdr_reject(pvalues, alpha, dep) result(reject)
        real(dp), intent(in) :: pvalues(:) !! Raw p-values to test by the Benjamini-Hochberg or dependency-adjusted rule.
        real(dp), intent(in), optional :: alpha !! Upstream confidence level; FDR level is 1-alpha, default alpha is 0.95.
        integer, intent(in), optional :: dep !! Dependency flag; negative applies the harmonic Benjamini-Yekutieli correction.
        real(dp), allocatable :: p(:)
        real(dp) :: a
        real(dp) :: c
        real(dp) :: threshold
        integer :: i
        integer :: d

        a = 0.95_dp
        if (present(alpha)) a = alpha
        d = 1
        if (present(dep)) d = dep
        p = pvalues
        call sort_real(p)
        c = 1.0_dp
        if (d < 0) then
            c = 0.0_dp
            do i = 1, size(p)
                c = c + 1.0_dp/real(i, dp)
            end do
        end if
        reject = .false.
        do i = 1, size(p)
            threshold = (1.0_dp - a)*real(i, dp)/(real(size(p), dp)*c)
            if (p(i) < threshold) then
                reject = .true.
                return
            end if
        end do
    end function fdr_reject

    pure real(dp) function fdr_adjusted_pvalue(pvalues, dep) result(pvalue)
        real(dp), intent(in) :: pvalues(:) !! Raw p-values whose minimum step-up adjusted value is requested.
        integer, intent(in), optional :: dep !! Dependency flag; negative applies the harmonic Benjamini-Yekutieli correction.
        real(dp), allocatable :: p(:)
        real(dp) :: c
        integer :: i
        integer :: d

        d = 1
        if (present(dep)) d = dep
        p = pvalues
        call sort_real(p)
        c = 1.0_dp
        if (d < 0) then
            c = 0.0_dp
            do i = 1, size(p)
                c = c + 1.0_dp/real(i, dp)
            end do
        end if
        pvalue = 1.0_dp
        do i = 1, size(p)
            pvalue = min(pvalue, p(i)*real(size(p), dp)*c/real(i, dp))
        end do
        pvalue = max(0.0_dp, min(1.0_dp, pvalue))
    end function fdr_adjusted_pvalue

    pure real(dp) function pred_mse(yobs, ypred) result(value)
        real(dp), intent(in) :: yobs(:) !! Observed scalar responses.
        real(dp), intent(in) :: ypred(:) !! Predicted scalar responses paired with yobs.
        value = sum((ypred - yobs)**2)/real(size(yobs), dp)
    end function pred_mse

    pure real(dp) function pred_rmse(yobs, ypred) result(value)
        real(dp), intent(in) :: yobs(:) !! Observed scalar responses.
        real(dp), intent(in) :: ypred(:) !! Predicted scalar responses paired with yobs.
        value = sqrt(pred_mse(yobs, ypred))
    end function pred_rmse

    pure real(dp) function pred_mae(yobs, ypred) result(value)
        real(dp), intent(in) :: yobs(:) !! Observed scalar responses.
        real(dp), intent(in) :: ypred(:) !! Predicted scalar responses paired with yobs.
        value = sum(abs(ypred - yobs))/real(size(yobs), dp)
    end function pred_mae

    pure subroutine prediction_measures(yobs, ypred, measures)
        real(dp), intent(in) :: yobs(:) !! Observed scalar responses.
        real(dp), intent(in) :: ypred(:) !! Predicted scalar responses paired with yobs.
        real(dp), intent(out) :: measures(3) !! Output [RMSE, MAE, MSE] in the conventional fda.usc order.
        measures(1) = pred_rmse(yobs, ypred)
        measures(2) = pred_mae(yobs, ypred)
        measures(3) = pred_mse(yobs, ypred)
    end subroutine prediction_measures

    pure real(dp) function confusion_accuracy(tab) result(value)
        real(dp), intent(in) :: tab(:, :) !! Square confusion matrix with observed classes in rows and predicted classes in columns.
        integer :: i
        value = 0.0_dp
        do i = 1, min(size(tab, 1), size(tab, 2))
            value = value + tab(i, i)
        end do
        if (sum(tab) > 0.0_dp) value = value/sum(tab)
    end function confusion_accuracy

    pure real(dp) function confusion_kappa(tab) result(value)
        real(dp), intent(in) :: tab(:, :) !! Square confusion matrix used for unweighted Cohen kappa.
        real(dp), allocatable :: rows(:)
        real(dp), allocatable :: cols(:)
        real(dp) :: total
        real(dp) :: p0
        real(dp) :: pe
        integer :: i

        total = sum(tab)
        if (total <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        rows = sum(tab, dim=2)/total
        cols = sum(tab, dim=1)/total
        p0 = 0.0_dp
        do i = 1, min(size(tab, 1), size(tab, 2))
            p0 = p0 + tab(i, i)/total
        end do
        pe = dot_product(rows, cols)
        if (abs(1.0_dp - pe) <= tiny(1.0_dp)) then
            value = 0.0_dp
        else
            value = 1.0_dp - (1.0_dp - p0)/(1.0_dp - pe)
        end if
    end function confusion_kappa

    pure real(dp) function confusion_iou(tab) result(value)
        real(dp), intent(in) :: tab(:, :) !! Square confusion matrix for aggregate intersection-over-union.
        real(dp) :: tp_total
        real(dp) :: union_total
        real(dp) :: tp
        real(dp) :: fp
        real(dp) :: fn
        integer :: j

        tp_total = 0.0_dp
        union_total = 0.0_dp
        do j = 1, size(tab, 1)
            tp = tab(j, j)
            fp = sum(tab(j, :)) - tp
            fn = sum(tab(:, j)) - tp
            tp_total = tp_total + tp
            union_total = union_total + tp + fp + fn
        end do
        if (union_total <= 0.0_dp) then
            value = 0.0_dp
        else
            value = tp_total/union_total
        end if
    end function confusion_iou

    pure subroutine confusion_iou_by_class(tab, values)
        real(dp), intent(in) :: tab(:, :) !! Square confusion matrix for classwise intersection-over-union.
        real(dp), intent(out) :: values(:) !! Per-class intersection-over-union, size at least the number of classes.
        real(dp) :: tp
        real(dp) :: fp
        real(dp) :: fn
        real(dp) :: denom
        integer :: j

        do j = 1, size(tab, 1)
            tp = tab(j, j)
            fp = sum(tab(j, :)) - tp
            fn = sum(tab(:, j)) - tp
            denom = tp + fp + fn
            if (denom > 0.0_dp) then
                values(j) = tp/denom
            else
                values(j) = 0.0_dp
            end if
        end do
    end subroutine confusion_iou_by_class

    pure subroutine confusion_recall(tab, values)
        real(dp), intent(in) :: tab(:, :) !! Square confusion matrix with observed classes in rows.
        real(dp), intent(out) :: values(:) !! Per-class recall values.
        real(dp) :: denom
        integer :: i

        do i = 1, size(tab, 1)
            denom = sum(tab(i, :))
            if (denom > 0.0_dp) then
                values(i) = tab(i, i)/denom
            else
                values(i) = 0.0_dp
            end if
        end do
    end subroutine confusion_recall

    pure subroutine confusion_precision(tab, values)
        real(dp), intent(in) :: tab(:, :) !! Square confusion matrix with predicted classes in columns.
        real(dp), intent(out) :: values(:) !! Per-class precision values.
        real(dp) :: denom
        integer :: i

        do i = 1, size(tab, 2)
            denom = sum(tab(:, i))
            if (denom > 0.0_dp) then
                values(i) = tab(i, i)/denom
            else
                values(i) = 0.0_dp
            end if
        end do
    end subroutine confusion_precision

    pure subroutine confusion_f1(tab, values)
        real(dp), intent(in) :: tab(:, :) !! Square confusion matrix used to compute classwise F1 scores.
        real(dp), intent(out) :: values(:) !! Per-class harmonic means of precision and recall.
        real(dp), allocatable :: recall(:)
        real(dp), allocatable :: precision(:)
        integer :: i

        allocate(recall(size(tab, 1)), precision(size(tab, 1)))
        call confusion_recall(tab, recall)
        call confusion_precision(tab, precision)
        do i = 1, size(tab, 1)
            if (recall(i) + precision(i) > 0.0_dp) then
                values(i) = 2.0_dp*recall(i)*precision(i)/(recall(i) + precision(i))
            else
                values(i) = 0.0_dp
            end if
        end do
    end subroutine confusion_f1

    pure subroutine class_weights(labels, weights, inverse_frequency)
        integer, intent(in) :: labels(:) !! Positive integer class labels for each observation.
        real(dp), intent(out) :: weights(:) !! Observation weights, equal or inverse class-frequency scaled.
        logical, intent(in), optional :: inverse_frequency !! If true, use inverse class frequencies; otherwise use unit weights.
        integer, allocatable :: counts(:)
        integer :: i
        integer :: k
        integer :: nclass
        logical :: inverse

        inverse = .false.
        if (present(inverse_frequency)) inverse = inverse_frequency
        if (.not. inverse) then
            weights = 1.0_dp
            return
        end if
        nclass = maxval(labels)
        allocate(counts(nclass))
        counts = 0
        do i = 1, size(labels)
            if (labels(i) >= 1 .and. labels(i) <= nclass) counts(labels(i)) = counts(labels(i)) + 1
        end do
        k = count(counts > 0)
        do i = 1, size(labels)
            if (labels(i) >= 1 .and. labels(i) <= nclass .and. counts(labels(i)) > 0) then
                weights(i) = real(size(labels), dp)/(real(k, dp)*real(counts(labels(i)), dp))
            else
                weights(i) = 0.0_dp
            end if
        end do
    end subroutine class_weights

    pure subroutine functional_mean(curves, mean_curve)
        real(dp), intent(in) :: curves(:, :) !! Functional observations, one curve per row.
        real(dp), intent(out) :: mean_curve(:) !! Pointwise arithmetic mean across curve rows.
        mean_curve = sum(curves, dim=1)/real(size(curves, 1), dp)
    end subroutine functional_mean

    pure subroutine functional_variance(curves, variance_curve)
        real(dp), intent(in) :: curves(:, :) !! Functional observations, one curve per row.
        real(dp), intent(out) :: variance_curve(:) !! Pointwise population variance with denominator equal to number of curves.
        real(dp), allocatable :: mean_curve(:)
        integer :: i

        allocate(mean_curve(size(curves, 2)))
        call functional_mean(curves, mean_curve)
        variance_curve = 0.0_dp
        do i = 1, size(curves, 1)
            variance_curve = variance_curve + (curves(i, :) - mean_curve)**2
        end do
        variance_curve = variance_curve/real(size(curves, 1), dp)
    end subroutine functional_variance

    pure subroutine center_curves(curves, centered, mean_curve)
        real(dp), intent(in) :: curves(:, :) !! Functional observations, one curve per row.
        real(dp), intent(out) :: centered(:, :) !! Input curves after subtracting their pointwise sample mean.
        real(dp), intent(out) :: mean_curve(:) !! Pointwise arithmetic mean removed from every curve.
        integer :: i

        call functional_mean(curves, mean_curve)
        do i = 1, size(curves, 1)
            centered(i, :) = curves(i, :) - mean_curve
        end do
    end subroutine center_curves

    pure subroutine difference_derivative(t, curves, nderiv, derivative)
        real(dp), intent(in) :: t(:) !! Strictly increasing observation grid for finite differences.
        real(dp), intent(in) :: curves(:, :) !! Functional observations, curves in rows.
        integer, intent(in) :: nderiv !! Positive derivative order; repeated finite differences are used.
        real(dp), intent(out) :: derivative(:, :) !! Derivative estimates on the original grid, with endpoint replication/averaging.
        real(dp), allocatable :: work(:)
        real(dp), allocatable :: next(:)
        real(dp), allocatable :: grid(:)
        real(dp), allocatable :: next_grid(:)
        integer :: i
        integer :: d
        integer :: j
        integer :: m

        derivative = 0.0_dp
        do i = 1, size(curves, 1)
            work = curves(i, :)
            grid = t
            do d = 1, max(1, nderiv)
                m = size(work) - 1
                if (m < 1) exit
                allocate(next(m), next_grid(m))
                do j = 1, m
                    next(j) = (work(j + 1) - work(j))/(grid(j + 1) - grid(j))
                    next_grid(j) = 0.5_dp*(grid(j + 1) + grid(j))
                end do
                call move_alloc(next, work)
                call move_alloc(next_grid, grid)
            end do
            call interpolate_to_grid(grid, work, t, derivative(i, :))
        end do
    end subroutine difference_derivative

    pure subroutine interpolate_to_grid(x, y, target, values)
        real(dp), intent(in) :: x(:) !! Strictly increasing source grid.
        real(dp), intent(in) :: y(:) !! Source values on x.
        real(dp), intent(in) :: target(:) !! Target grid on which values are linearly interpolated, with endpoint replication.
        real(dp), intent(out) :: values(:) !! Interpolated values on target.
        integer :: i
        integer :: j
        real(dp) :: f

        if (size(x) == 1) then
            values = y(1)
            return
        end if
        j = 1
        do i = 1, size(target)
            if (target(i) <= x(1)) then
                values(i) = y(1)
            else if (target(i) >= x(size(x))) then
                values(i) = y(size(y))
            else
                do while (j < size(x) - 1 .and. target(i) > x(j + 1))
                    j = j + 1
                end do
                f = (target(i) - x(j))/(x(j + 1) - x(j))
                values(i) = (1.0_dp - f)*y(j) + f*y(j + 1)
            end if
        end do
    end subroutine interpolate_to_grid

    pure subroutine sort_real(x)
        real(dp), intent(inout) :: x(:) !! Values sorted in ascending order in place.
        real(dp) :: key
        integer :: i
        integer :: j

        do i = 2, size(x)
            key = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= key) exit
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = key
        end do
    end subroutine sort_real

end module fdausc_statistics
