! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm_pedigree
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use mcmcglmm_rng, only : rng_state
    use mcmcglmm_matrix, only : sample_mvn_covariance
    implicit none
    private

    public :: pedigree_relationship
    public :: pedigree_inverse
    public :: prune_pedigree_mask
    public :: breeding_values_pedigree

contains

    pure subroutine pedigree_relationship(dam, sire, relationship, inbreeding, mendelian_variance, info)
        integer, intent(in) :: dam(:) !! One-based dam row index; zero is unknown and known parents must precede offspring.
        integer, intent(in) :: sire(:) !! One-based sire row index; zero is unknown and known parents must precede offspring.
        real(dp), allocatable, intent(out) :: relationship(:, :) !! Allocated numerator relationship matrix A in pedigree row order.
        real(dp), allocatable, intent(out) :: inbreeding(:) !! Allocated inbreeding coefficients F_i = A_ii - 1.
        real(dp), allocatable, intent(out) :: mendelian_variance(:) !! Allocated diagonal D terms in the decomposition A = T D T^T.
        integer, intent(out) :: info !! Zero on success; nonzero for size mismatch or invalid forward/self parent references.
        real(dp) :: parent_sum
        integer :: i
        integer :: j
        integer :: n

        info = 0
        n = size(dam)
        if (size(sire) /= n) then
            allocate(relationship(0, 0), inbreeding(0), mendelian_variance(0))
            info = 1
            return
        end if
        do i = 1, n
            if (dam(i) < 0 .or. dam(i) >= i .or. sire(i) < 0 .or. sire(i) >= i) then
                allocate(relationship(0, 0), inbreeding(0), mendelian_variance(0))
                info = 2
                return
            end if
        end do

        allocate(relationship(n, n), inbreeding(n), mendelian_variance(n))
        relationship = 0.0_dp
        inbreeding = 0.0_dp
        mendelian_variance = 1.0_dp

        do i = 1, n
            do j = 1, i - 1
                parent_sum = 0.0_dp
                if (dam(i) > 0) parent_sum = parent_sum + relationship(dam(i), j)
                if (sire(i) > 0) parent_sum = parent_sum + relationship(sire(i), j)
                relationship(i, j) = 0.5_dp * parent_sum
                relationship(j, i) = relationship(i, j)
            end do
            if (dam(i) > 0 .and. sire(i) > 0) then
                relationship(i, i) = 1.0_dp + 0.5_dp * relationship(dam(i), sire(i))
            else
                relationship(i, i) = 1.0_dp
            end if
            inbreeding(i) = relationship(i, i) - 1.0_dp
            if (dam(i) > 0) mendelian_variance(i) = mendelian_variance(i) - 0.25_dp * (1.0_dp + inbreeding(dam(i)))
            if (sire(i) > 0) mendelian_variance(i) = mendelian_variance(i) - 0.25_dp * (1.0_dp + inbreeding(sire(i)))
        end do
    end subroutine pedigree_relationship

    pure subroutine pedigree_inverse(dam, sire, inverse_relationship, inbreeding, mendelian_variance, info)
        integer, intent(in) :: dam(:) !! One-based dam row index; zero is unknown and known parents must precede offspring.
        integer, intent(in) :: sire(:) !! One-based sire row index; zero is unknown and known parents must precede offspring.
        real(dp), allocatable, intent(out) :: inverse_relationship(:, :) !! Dense inverse numerator relationship matrix A.
        real(dp), allocatable, intent(out) :: inbreeding(:) !! Allocated Meuwissen-Luo inbreeding coefficients in pedigree order.
        real(dp), allocatable, intent(out) :: mendelian_variance(:) !! Allocated diagonal D terms used by the pedigree recursion.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid pedigree ordering or matrix inversion failure.
        real(dp), allocatable :: relationship(:, :)

        call pedigree_relationship(dam, sire, relationship, inbreeding, mendelian_variance, info)
        if (info /= 0) then
            allocate(inverse_relationship(0, 0))
            return
        end if
        call inverse_matrix(relationship, inverse_relationship, info)
    end subroutine pedigree_inverse

    pure subroutine prune_pedigree_mask(dam, sire, requested, keep_mask, info)
        integer, intent(in) :: dam(:) !! One-based dam row indices with zero for missing parents.
        integer, intent(in) :: sire(:) !! One-based sire row indices with zero for missing parents.
        integer, intent(in) :: requested(:) !! One-based pedigree rows that must be retained together with all their ancestors.
        logical, allocatable, intent(out) :: keep_mask(:) !! Mask for requested individuals and all recursively required ancestors.
        integer, intent(out) :: info !! Zero on success; nonzero for vector-size mismatch or an out-of-range requested/parent index.
        integer :: i
        integer :: n
        integer :: parent
        logical :: changed

        info = 0
        n = size(dam)
        if (size(sire) /= n) then
            allocate(keep_mask(0))
            info = 1
            return
        end if
        if (any(requested < 1) .or. any(requested > n)) then
            allocate(keep_mask(0))
            info = 2
            return
        end if
        if (any(dam < 0) .or. any(dam > n) .or. any(sire < 0) .or. any(sire > n)) then
            allocate(keep_mask(0))
            info = 3
            return
        end if

        allocate(keep_mask(n))
        keep_mask = .false.
        keep_mask(requested) = .true.
        do
            changed = .false.
            do i = 1, n
                if (.not. keep_mask(i)) cycle
                parent = dam(i)
                if (parent > 0) then
                    if (.not. keep_mask(parent)) then
                        keep_mask(parent) = .true.
                        changed = .true.
                    end if
                end if
                parent = sire(i)
                if (parent > 0) then
                    if (.not. keep_mask(parent)) then
                        keep_mask(parent) = .true.
                        changed = .true.
                    end if
                end if
            end do
            if (.not. changed) exit
        end do
    end subroutine prune_pedigree_mask

    pure subroutine breeding_values_pedigree(state, dam, sire, covariance, values, info, groups, group_means)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the Mendelian sampling draws.
        integer, intent(in) :: dam(:) !! One-based dam row indices; zero denotes an unknown dam and parents must precede offspring.
        integer, intent(in) :: sire(:) !! One-based sire indices; zero is unknown and parents must precede offspring.
        real(dp), intent(in) :: covariance(:, :) !! Positive-definite additive genetic covariance matrix among traits.
        real(dp), allocatable, intent(out) :: values(:, :) !! Allocated n-individual by n-trait breeding-value matrix.
        integer, intent(out) :: info !! Zero on success; nonzero for pedigree, covariance, group, or factorization failure.
        integer, intent(in), optional :: groups(:) !! Optional one-based genetic group per individual; default is group one.
        real(dp), intent(in), optional :: group_means(:, :) !! Optional group-by-trait founder means used for missing parents.
        real(dp), allocatable :: inbreeding(:)
        real(dp), allocatable :: mendelian_variance(:)
        real(dp), allocatable :: relationship(:, :)
        real(dp), allocatable :: innovation(:)
        real(dp), allocatable :: zero_mean(:)
        real(dp), allocatable :: scaled_covariance(:, :)
        integer :: group_index
        integer :: i
        integer :: n
        integer :: ntrait

        info = 0
        n = size(dam)
        ntrait = size(covariance, 1)
        if (size(sire) /= n .or. size(covariance, 2) /= ntrait) then
            allocate(values(0, 0))
            info = 1
            return
        end if
        if (present(groups)) then
            if (size(groups) /= n .or. any(groups < 1)) then
                allocate(values(0, 0))
                info = 2
                return
            end if
        end if
        if (present(group_means)) then
            if (size(group_means, 2) /= ntrait) then
                allocate(values(0, 0))
                info = 3
                return
            end if
            if (present(groups)) then
                if (maxval(groups) > size(group_means, 1)) then
                    allocate(values(0, 0))
                    info = 4
                    return
                end if
            end if
        end if

        call pedigree_relationship(dam, sire, relationship, inbreeding, mendelian_variance, info)
        if (info /= 0) then
            allocate(values(0, 0))
            return
        end if
        allocate(values(n, ntrait), zero_mean(ntrait), scaled_covariance(ntrait, ntrait))
        values = 0.0_dp
        zero_mean = 0.0_dp

        do i = 1, n
            scaled_covariance = mendelian_variance(i) * covariance
            call sample_mvn_covariance(state, zero_mean, scaled_covariance, innovation, info)
            if (info /= 0) return
            values(i, :) = innovation
            group_index = 1
            if (present(groups)) group_index = groups(i)
            if (dam(i) > 0) then
                values(i, :) = values(i, :) + 0.5_dp * values(dam(i), :)
            else if (present(group_means)) then
                values(i, :) = values(i, :) + 0.5_dp * group_means(group_index, :)
            end if
            if (sire(i) > 0) then
                values(i, :) = values(i, :) + 0.5_dp * values(sire(i), :)
            else if (present(group_means)) then
                values(i, :) = values(i, :) + 0.5_dp * group_means(group_index, :)
            end if
        end do
    end subroutine breeding_values_pedigree

end module mcmcglmm_pedigree
