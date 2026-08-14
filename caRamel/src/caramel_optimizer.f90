module caramel_optimizer
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
    use caramel_kinds, only: dp
    use caramel_generation, only: index_block, new_xval
    use caramel_random, only: random_permutation, random_uniform
    use caramel_utils, only: decrease_pop
    implicit none
    private

    type, public :: caramel_options
        integer :: popsize = 100
        integer :: archsize = 100
        integer :: maxrun = 1000
        integer :: repart_gene(4) = [5, 5, 5, 5]
        integer :: gpp = 0
        integer :: nout = 0
        logical :: sensitivity = .false.
        real(dp) :: sensitivity_step = 1.0e-4_dp
    end type caramel_options

    type, public :: caramel_result
        logical :: success = .false.
        character(len=:), allocatable :: message
        real(dp), allocatable :: parameters(:,:)
        real(dp), allocatable :: objectives(:,:)
        real(dp), allocatable :: derivatives(:,:,:)
        real(dp), allocatable :: save_crit(:,:)
        real(dp), allocatable :: total_pop(:,:)
        integer :: gpp = 0
        integer :: nrun = 0
        integer :: ngen = 0
    end type caramel_result

    abstract interface
        subroutine objective_function(x, values)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: values(:)
        end subroutine objective_function
    end interface

    public :: objective_function, caramel_optimize

contains

    subroutine caramel_optimize(nobj, nvar, minmax, bounds, prec, func, result, options, blocks, initial_population)
        integer, intent(in) :: nobj, nvar
        logical, intent(in) :: minmax(:)
        real(dp), intent(in) :: bounds(:,:), prec(:)
        procedure(objective_function) :: func
        type(caramel_result), intent(out) :: result
        type(caramel_options), intent(in), optional :: options
        type(index_block), intent(in), optional :: blocks(:)
        real(dp), intent(in), optional :: initial_population(:,:)
        type(caramel_options) :: opt
        real(dp), allocatable :: pop(:,:), x(:,:), projected(:,:), values(:,:), pop1(:,:), arch(:,:)
        real(dp), allocatable :: param(:,:), crit(:,:), sp(:), best(:), record(:), pert_values(:,:)
        integer, allocatable :: ind_arch(:), ind_pop(:), perm(:)
        logical, allocatable :: valid(:)
        integer :: nout, nrun, ngen, gpp, i, j, nvalid, nfront
        logical :: fireworks
        real(dp) :: dx, nanv

        call initialize_empty_result(result, nvar, nobj)
        opt = caramel_options()
        if (present(options)) opt = options

        if (nobj <= 1) then
            call fail(result, "the number of objectives must be greater than one")
            return
        end if
        if (nvar < 1) then
            call fail(result, "the number of variables must be greater than zero")
            return
        end if
        if (size(minmax) /= nobj .or. size(prec) /= nobj) then
            call fail(result, "minmax/prec has an incorrect dimension")
            return
        end if
        if (size(bounds,1) /= nvar .or. size(bounds,2) /= 2) then
            call fail(result, "bounds has an incorrect dimension")
            return
        end if
        if (any(bounds(:,2) < bounds(:,1))) then
            call fail(result, "each upper bound must be at least its lower bound")
            return
        end if
        if (any(prec <= 0.0_dp)) then
            call fail(result, "prec values must be strictly positive")
            return
        end if
        if (opt%popsize < 1 .or. opt%archsize < 1 .or. opt%maxrun < 1) then
            call fail(result, "popsize, archsize, and maxrun must be strictly positive")
            return
        end if
        if (any(opt%repart_gene <= 0)) then
            call fail(result, "all repart_gene values must be strictly positive")
            return
        end if
        if (opt%gpp < 0) then
            call fail(result, "gpp must be zero (automatic) or positive")
            return
        end if
        nout = opt%nout
        if (nout == 0) nout = nobj
        if (nout < nobj) then
            call fail(result, "nout must be at least nobj")
            return
        end if

        nanv = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(sp(nvar))
        sp = (bounds(:,2) - bounds(:,1)) / (2.0_dp * sqrt(3.0_dp))
        gpp = opt%gpp
        if (gpp == 0) then
            gpp = ceiling(real(nvar*(nobj+1)*4,dp) / real(sum(opt%repart_gene),dp))
            gpp = max(1,gpp)
        end if
        result%gpp = gpp
        nrun = 0
        ngen = 0
        if (allocated(result%save_crit)) deallocate(result%save_crit)
        allocate(result%save_crit(0,nobj+1))

        if (present(initial_population)) then
            if (size(initial_population,2) == nvar) then
                x = initial_population
                call evaluate_population(x, nout, func, values, valid)
                nrun = nrun + size(x,1)
                nvalid = count(valid)
                if (nvalid == 0) then
                    call fail(result, "initial population contains no feasible objective evaluations")
                    result%nrun = nrun
                    return
                end if
                allocate(pop(nvalid,nvar+nout))
                pop(:,1:nvar) = pack_rows(x, valid)
                pop(:,nvar+1:nvar+nout) = pack_rows(values, valid)
            else if (size(initial_population,2) >= nvar+nout) then
                pop = initial_population(:,1:nvar+nout)
            else
                call fail(result, "initial_population must contain nvar or at least nvar+nout columns")
                return
            end if
        end if

        do while (nrun < opt%maxrun)
            ngen = ngen + 1
            if (.not. allocated(pop)) then
                allocate(x(opt%popsize,nvar), projected(opt%popsize,nobj), perm(opt%popsize))
                projected = nanv
                do j = 1, nvar
                    call random_permutation(opt%popsize, perm)
                    do i = 1, opt%popsize
                        x(i,j) = bounds(j,1) + &
                            (real(perm(i)-1,dp) + random_uniform()) / real(opt%popsize,dp) * &
                            (bounds(j,2) - bounds(j,1))
                    end do
                end do
            else
                if (size(pop,1) < 4) then
                    call fail(result, "the number of feasible points is insufficient; increase the population size")
                    result%nrun = nrun
                    result%ngen = ngen - 1
                    return
                end if
                param = pop(:,1:nvar)
                crit = pop(:,nvar+1:nvar+nobj)
                fireworks = modulo(ngen,gpp) == 0
                if (present(blocks)) then
                    call new_xval(param, crit, minmax, sp, bounds, opt%repart_gene, fireworks, x, projected, blocks)
                else
                    call new_xval(param, crit, minmax, sp, bounds, opt%repart_gene, fireworks, x, projected)
                end if
                if (size(x,1) == 0) then
                    call fail(result, "no new candidate points could be generated")
                    result%nrun = nrun
                    result%ngen = ngen - 1
                    return
                end if
            end if

            call evaluate_population(x, nout, func, values, valid)
            nrun = nrun + size(x,1)
            nvalid = count(valid)
            if (nvalid == 0) then
                call fail(result, "no feasible points were returned by the objective function")
                result%nrun = nrun
                result%ngen = ngen
                return
            end if
            call append_evaluated(pop, x, values, valid, pop1)

            call decrease_pop(pop1(:,nvar+1:nvar+nobj), minmax, prec, opt%archsize, opt%popsize, &
                              ind_arch, ind_pop)
            if (allocated(arch)) deallocate(arch)
            allocate(arch(size(ind_arch),nvar+nout))
            if (size(ind_arch) > 0) arch = pop1(ind_arch,:)
            call rebuild_population(pop1, ind_arch, ind_pop, pop)

            nfront = size(arch,1)
            if (nfront == 0) then
                call fail(result, "Pareto archive became empty")
                result%nrun = nrun
                result%ngen = ngen
                return
            end if
            allocate(best(nobj), record(nobj+1))
            do j = 1, nobj
                if (minmax(j)) then
                    best(j) = maxval(arch(:,nvar+j))
                else
                    best(j) = minval(arch(:,nvar+j))
                end if
            end do
            record = [real(nrun,dp), best]
            call append_record(result%save_crit, record)
            deallocate(best, record, x, projected, values, valid, pop1, ind_arch, ind_pop)
        end do

        result%success = .true.
        result%message = "ok"
        result%nrun = nrun
        result%ngen = ngen
        result%parameters = arch(:,1:nvar)
        result%objectives = arch(:,nvar+1:nvar+nout)
        result%total_pop = pop

        if (opt%sensitivity) then
            dx = opt%sensitivity_step
            if (dx <= 0.0_dp) dx = 1.0e-4_dp
            nfront = size(arch,1)
            if (allocated(result%derivatives)) deallocate(result%derivatives)
            allocate(result%derivatives(nfront,nvar,nobj))
            result%derivatives = nanv
            allocate(x(nfront,nvar))
            x = arch(:,1:nvar)
            do j = 1, nvar
                x(:,j) = arch(:,j) + dx
                call evaluate_population(x, nout, func, pert_values, valid)
                nrun = nrun + nfront
                do i = 1, nfront
                    if (valid(i)) then
                        result%derivatives(i,j,:) = &
                            (pert_values(i,1:nobj) - arch(i,nvar+1:nvar+nobj)) / dx
                    end if
                end do
                x(:,j) = arch(:,j)
                deallocate(pert_values, valid)
            end do
            result%nrun = nrun
        else
            if (allocated(result%derivatives)) deallocate(result%derivatives)
            allocate(result%derivatives(0,0,0))
        end if
    end subroutine caramel_optimize

    subroutine evaluate_population(x, nout, func, values, valid)
        real(dp), intent(in) :: x(:,:)
        integer, intent(in) :: nout
        procedure(objective_function) :: func
        real(dp), allocatable, intent(out) :: values(:,:)
        logical, allocatable, intent(out) :: valid(:)
        integer :: i

        allocate(values(size(x,1),nout), valid(size(x,1)))
        do i = 1, size(x,1)
            call func(x(i,:), values(i,:))
            valid(i) = .not. any(ieee_is_nan(values(i,:)))
        end do
    end subroutine evaluate_population

    subroutine append_evaluated(pop, x, values, valid, combined)
        real(dp), allocatable, intent(in) :: pop(:,:)
        real(dp), intent(in) :: x(:,:), values(:,:)
        logical, intent(in) :: valid(:)
        real(dp), allocatable, intent(out) :: combined(:,:)
        real(dp), allocatable :: xv(:,:), vv(:,:)
        integer :: oldn, nvar, nout, nv

        nvar = size(x,2)
        nout = size(values,2)
        nv = count(valid)
        xv = pack_rows(x,valid)
        vv = pack_rows(values,valid)
        oldn = 0
        if (allocated(pop)) oldn = size(pop,1)
        allocate(combined(oldn+nv,nvar+nout))
        if (oldn > 0) combined(1:oldn,:) = pop
        combined(oldn+1:oldn+nv,1:nvar) = xv
        combined(oldn+1:oldn+nv,nvar+1:nvar+nout) = vv
    end subroutine append_evaluated

    function pack_rows(a, mask) result(b)
        real(dp), intent(in) :: a(:,:)
        logical, intent(in) :: mask(:)
        real(dp), allocatable :: b(:,:)
        integer :: i, j

        if (size(mask) /= size(a,1)) error stop "pack_rows: inconsistent mask"
        allocate(b(count(mask),size(a,2)))
        j = 0
        do i = 1, size(a,1)
            if (mask(i)) then
                j = j + 1
                b(j,:) = a(i,:)
            end if
        end do
    end function pack_rows

    subroutine rebuild_population(source, arch_idx, pop_idx, pop)
        real(dp), intent(in) :: source(:,:)
        integer, intent(in) :: arch_idx(:), pop_idx(:)
        real(dp), allocatable, intent(inout) :: pop(:,:)
        real(dp), allocatable :: tmp(:,:)
        integer :: na, np

        na = size(arch_idx)
        np = size(pop_idx)
        allocate(tmp(na+np,size(source,2)))
        if (na > 0) tmp(1:na,:) = source(arch_idx,:)
        if (np > 0) tmp(na+1:na+np,:) = source(pop_idx,:)
        call move_alloc(tmp,pop)
    end subroutine rebuild_population

    subroutine append_record(records, record)
        real(dp), allocatable, intent(inout) :: records(:,:)
        real(dp), intent(in) :: record(:)
        real(dp), allocatable :: tmp(:,:)
        integer :: n

        if (size(records,2) /= size(record)) error stop "append_record: inconsistent record size"
        n = size(records,1)
        allocate(tmp(n+1,size(record)))
        if (n > 0) tmp(1:n,:) = records
        tmp(n+1,:) = record
        call move_alloc(tmp,records)
    end subroutine append_record

    subroutine initialize_empty_result(result, nvar, nobj)
        type(caramel_result), intent(out) :: result
        integer, intent(in) :: nvar, nobj
        result%success = .false.
        result%message = "not run"
        allocate(result%parameters(0,max(0,nvar)))
        allocate(result%objectives(0,max(0,nobj)))
        allocate(result%derivatives(0,0,0))
        allocate(result%save_crit(0,max(0,nobj+1)))
        allocate(result%total_pop(0,max(0,nvar+nobj)))
    end subroutine initialize_empty_result

    subroutine fail(result, message)
        type(caramel_result), intent(inout) :: result
        character(len=*), intent(in) :: message
        result%success = .false.
        result%message = message
    end subroutine fail

end module caramel_optimizer
