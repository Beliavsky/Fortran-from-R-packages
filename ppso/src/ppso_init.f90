module ppso_init
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use ppso_kinds, only : dp
    use ppso_rng, only : rng_state, latin_hypercube
    implicit none
    private
    public :: validate_bounds, generate_initial_population, append_pending

contains

    subroutine validate_bounds(lower, upper)
        real(dp), intent(in) :: lower(:), upper(:)
        integer :: j

        if (size(lower) /= size(upper)) error stop "ppso: lower/upper size mismatch"
        if (size(lower) == 0) error stop "ppso: at least one parameter is required"
        do j = 1, size(lower)
            if (.not. ieee_is_finite(lower(j)) .or. .not. ieee_is_finite(upper(j))) &
                error stop "ppso: bounds must be finite"
            if (lower(j) > upper(j)) error stop "ppso: reversed parameter bounds"
        end do
    end subroutine validate_bounds

    subroutine generate_initial_population(rng, lower, upper, npart, lhs, initial, x, pending)
        type(rng_state), intent(inout) :: rng
        real(dp), intent(in) :: lower(:), upper(:)
        integer, intent(in) :: npart
        logical, intent(in) :: lhs
        real(dp), intent(in), optional :: initial(:,:)
        real(dp), allocatable, intent(out) :: x(:,:), pending(:,:)
        real(dp), allocatable :: u(:,:)
        integer :: npar, i, j, ninit, nuse, npending

        call validate_bounds(lower, upper)
        npar = size(lower)
        if (npart <= 0) error stop "ppso: number_of_particles must be positive"

        allocate(x(npar,npart), u(npar,npart))
        if (lhs) then
            call latin_hypercube(rng, u)
        else
            do i = 1, npart
                do j = 1, npar
                    u(j,i) = rng%uniform()
                end do
            end do
        end if
        do i = 1, npart
            x(:,i) = lower + (upper-lower)*u(:,i)
        end do

        ninit = 0
        if (present(initial)) then
            if (size(initial,1) /= npar) error stop "ppso: initial estimates have wrong row count"
            ninit = size(initial,2)
            do i = 1, ninit
                if (any(initial(:,i) < lower) .or. any(initial(:,i) > upper)) &
                    error stop "ppso: initial estimate is outside bounds"
            end do
        end if

        nuse = min(npart, ninit)
        if (nuse > 0) x(:,1:nuse) = initial(:,1:nuse)
        npending = max(0, ninit-nuse)
        allocate(pending(npar,npending))
        if (npending > 0) pending = initial(:,nuse+1:ninit)
    end subroutine generate_initial_population

    subroutine append_pending(pending, candidate, used)
        real(dp), allocatable, intent(inout) :: pending(:,:)
        real(dp), intent(out) :: candidate(:)
        logical, intent(out) :: used
        real(dp), allocatable :: tmp(:,:)
        integer :: npar, nleft

        used = .false.
        if (.not. allocated(pending)) return
        if (size(pending,2) == 0) return
        if (size(candidate) /= size(pending,1)) error stop "ppso: pending estimate dimension mismatch"

        candidate = pending(:,1)
        used = .true.
        npar = size(pending,1)
        nleft = size(pending,2)-1
        allocate(tmp(npar,nleft))
        if (nleft > 0) tmp = pending(:,2:)
        call move_alloc(tmp, pending)
    end subroutine append_pending

end module ppso_init
