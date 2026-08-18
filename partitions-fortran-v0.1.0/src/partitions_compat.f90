! Translation of computational code from R package partitions 1.10-9.
! Upstream authors: Robin K. S. Hankin; contributor Paul Egeler.
! Upstream DESCRIPTION declares: License: GPL.
! This Fortran translation is distributed under the same GPL terms.

module partitions_compat
    use partitions_kinds, only : i8
    use partitions_counts, only : partition_count, distinct_partition_count, restricted_partition_count, &
                                  block_partition_count
    use partitions_enumeration, only : core_parts => parts, core_diffparts => distinct_parts, &
        core_restrictedparts => restricted_parts, core_blockparts => block_parts, &
        core_firstpart => first_part, core_nextpart => next_part, core_islastpart => is_last_part, &
        core_firstdiffpart => first_distinct_part, core_nextdiffpart => next_distinct_part, &
        core_islastdiffpart => is_last_distinct_part, core_firstrestrictedpart => first_restricted_part, &
        core_nextrestrictedpart => next_restricted_part, core_islastrestrictedpart => is_last_restricted_part, &
        core_firstblockpart => first_block_part, core_nextblockpart => next_block_part, &
        core_islastblockpart => is_last_block_part, core_firstcomposition => first_composition, &
        core_nextcomposition => next_composition, core_islastcomposition => is_last_composition, &
        core_compositions => compositions, core_tobin => to_binary, core_todec => to_decimal, &
        core_comptobin => composition_to_binary, core_bintocomp => binary_to_composition
    use partitions_permutations, only : core_perms => permutations, core_plainperms => plain_permutations, &
                                      core_mset => multiset_permutations
    use partitions_sets, only : core_setparts => set_partitions, &
        core_restrictedsetparts => restricted_set_partitions, core_multiset => multiset_sequences, &
        core_multinomial => multinomial_permutations, core_allbinom => all_binomial, &
        core_genrif => generalized_riffles, core_riffle => riffles
    implicit none
    private

    public :: p, q, r, s
    public :: parts, diffparts, restrictedparts, blockparts, compositions
    public :: firstpart, nextpart, islastpart
    public :: firstdiffpart, nextdiffpart, islastdiffpart
    public :: firstrestrictedpart, nextrestrictedpart, islastrestrictedpart
    public :: firstblockpart, nextblockpart, islastblockpart
    public :: firstcomposition, nextcomposition, islastcomposition
    public :: setparts, restrictedsetparts, restrictedsetparts2
    public :: perms, plainperms, mset, multiset, multinomial, allbinom, genrif, riffle
    public :: tobin, todec, comptobin, bintocomp

contains

    function p(n) result(value)
        integer, intent(in) :: n
        integer(i8) :: value
        value = partition_count(n)
    end function p

    function q(n) result(value)
        integer, intent(in) :: n
        integer(i8) :: value
        value = distinct_partition_count(n)
    end function q

    function r(m, n, include_zero) result(value)
        integer, intent(in) :: m, n
        logical, intent(in), optional :: include_zero
        integer(i8) :: value
        if (present(include_zero)) then
            value = restricted_partition_count(m, n, include_zero)
        else
            value = restricted_partition_count(m, n)
        end if
    end function r

    function s(f, n, include_fewer) result(value)
        integer, intent(in) :: f(:), n
        logical, intent(in), optional :: include_fewer
        integer(i8) :: value
        if (present(include_fewer)) then
            value = block_partition_count(f, n, include_fewer)
        else
            value = block_partition_count(f, n)
        end if
    end function s

    function parts(n) result(out)
        integer, intent(in) :: n
        integer, allocatable :: out(:,:)
        out = core_parts(n)
    end function parts

    function diffparts(n) result(out)
        integer, intent(in) :: n
        integer, allocatable :: out(:,:)
        out = core_diffparts(n)
    end function diffparts

    function restrictedparts(n, m, include_zero, decreasing) result(out)
        integer, intent(in) :: n, m
        logical, intent(in), optional :: include_zero, decreasing
        integer, allocatable :: out(:,:)
        logical :: inc0, decr
        inc0 = .true.
        decr = .true.
        if (present(include_zero)) inc0 = include_zero
        if (present(decreasing)) decr = decreasing
        out = core_restrictedparts(n, m, inc0, decr)
    end function restrictedparts

    function blockparts(f, n, include_fewer) result(out)
        integer, intent(in) :: f(:)
        integer, intent(in), optional :: n
        logical, intent(in), optional :: include_fewer
        integer, allocatable :: out(:,:)
        logical :: fewer
        fewer = .false.
        if (present(include_fewer)) fewer = include_fewer
        if (present(n)) then
            out = core_blockparts(f, n, fewer)
        else
            out = core_blockparts(f, include_fewer=fewer)
        end if
    end function blockparts

    function compositions(n, m, include_zero) result(out)
        integer, intent(in) :: n
        integer, intent(in), optional :: m
        logical, intent(in), optional :: include_zero
        integer, allocatable :: out(:,:)
        logical :: inc0
        inc0 = .true.
        if (present(include_zero)) inc0 = include_zero
        if (present(m)) then
            out = core_compositions(n, m, inc0)
        else
            out = core_compositions(n, include_zero=inc0)
        end if
    end function compositions

    function firstpart(n) result(out)
        integer, intent(in) :: n
        integer, allocatable :: out(:)
        out = core_firstpart(n)
    end function firstpart

    function nextpart(part) result(out)
        integer, intent(in) :: part(:)
        integer, allocatable :: out(:)
        out = part
        call core_nextpart(out)
    end function nextpart

    logical function islastpart(part) result(last)
        integer, intent(in) :: part(:)
        last = core_islastpart(part)
    end function islastpart

    function firstdiffpart(n) result(out)
        integer, intent(in) :: n
        integer, allocatable :: out(:)
        out = core_firstdiffpart(n)
    end function firstdiffpart

    function nextdiffpart(part) result(out)
        integer, intent(in) :: part(:)
        integer, allocatable :: out(:)
        out = part
        call core_nextdiffpart(out)
    end function nextdiffpart

    logical function islastdiffpart(part) result(last)
        integer, intent(in) :: part(:)
        last = core_islastdiffpart(part)
    end function islastdiffpart

    function firstrestrictedpart(n, m, include_zero) result(out)
        integer, intent(in) :: n, m
        logical, intent(in), optional :: include_zero
        integer, allocatable :: out(:)
        logical :: inc0
        inc0 = .true.
        if (present(include_zero)) inc0 = include_zero
        out = core_firstrestrictedpart(n, m, inc0)
    end function firstrestrictedpart

    function nextrestrictedpart(part) result(out)
        integer, intent(in) :: part(:)
        integer, allocatable :: out(:)
        out = part
        call core_nextrestrictedpart(out)
    end function nextrestrictedpart

    logical function islastrestrictedpart(part) result(last)
        integer, intent(in) :: part(:)
        last = core_islastrestrictedpart(part)
    end function islastrestrictedpart

    function firstblockpart(f, n, include_fewer) result(out)
        integer, intent(in) :: f(:)
        integer, intent(in), optional :: n
        logical, intent(in), optional :: include_fewer
        integer, allocatable :: out(:)
        logical :: fewer
        fewer = .false.
        if (present(include_fewer)) fewer = include_fewer
        if (present(n)) then
            out = core_firstblockpart(f, n, fewer)
        else
            out = core_firstblockpart(f, include_fewer=fewer)
        end if
    end function firstblockpart

    function nextblockpart(part, f, n, include_fewer) result(out)
        integer, intent(in) :: part(:), f(:)
        integer, intent(in), optional :: n
        logical, intent(in), optional :: include_fewer
        integer, allocatable :: out(:)
        logical :: fewer
        fewer = .false.
        if (present(include_fewer)) fewer = include_fewer
        out = part
        if (present(n)) then
            call core_nextblockpart(out, f, n, fewer)
        else
            call core_nextblockpart(out, f, include_fewer=fewer)
        end if
    end function nextblockpart

    logical function islastblockpart(part, f, n, include_fewer) result(last)
        integer, intent(in) :: part(:), f(:)
        integer, intent(in), optional :: n
        logical, intent(in), optional :: include_fewer
        logical :: fewer
        fewer = .false.
        if (present(include_fewer)) fewer = include_fewer
        if (present(n)) then
            last = core_islastblockpart(part, f, n, fewer)
        else
            last = core_islastblockpart(part, f, include_fewer=fewer)
        end if
    end function islastblockpart

    function firstcomposition(n, m, include_zero) result(out)
        integer, intent(in) :: n
        integer, intent(in), optional :: m
        logical, intent(in), optional :: include_zero
        integer, allocatable :: out(:)
        logical :: inc0
        inc0 = .true.
        if (present(include_zero)) inc0 = include_zero
        if (present(m)) then
            out = core_firstcomposition(n, m, inc0)
        else
            out = core_firstcomposition(n, include_zero=inc0)
        end if
    end function firstcomposition

    function nextcomposition(comp, restricted, include_zero) result(out)
        integer, intent(in) :: comp(:)
        logical, intent(in) :: restricted
        logical, intent(in), optional :: include_zero
        integer, allocatable :: out(:)
        logical :: inc0
        inc0 = .true.
        if (present(include_zero)) inc0 = include_zero
        out = core_nextcomposition(comp, restricted, inc0)
    end function nextcomposition

    logical function islastcomposition(comp, restricted, include_zero) result(last)
        integer, intent(in) :: comp(:)
        logical, intent(in) :: restricted
        logical, intent(in), optional :: include_zero
        logical :: inc0
        inc0 = .true.
        if (present(include_zero)) inc0 = include_zero
        last = core_islastcomposition(comp, restricted, inc0)
    end function islastcomposition

    function setparts(block_sizes) result(out)
        integer, intent(in) :: block_sizes(:)
        integer, allocatable :: out(:,:)
        out = core_setparts(block_sizes)
    end function setparts

    function restrictedsetparts(block_sizes) result(out)
        integer, intent(in) :: block_sizes(:)
        integer, allocatable :: out(:,:)
        out = core_restrictedsetparts(block_sizes)
    end function restrictedsetparts

    function restrictedsetparts2(block_sizes) result(out)
        integer, intent(in) :: block_sizes(:)
        integer, allocatable :: out(:,:)
        out = core_restrictedsetparts(block_sizes)
    end function restrictedsetparts2

    function perms(n) result(out)
        integer, intent(in) :: n
        integer, allocatable :: out(:,:)
        out = core_perms(n)
    end function perms

    function plainperms(n) result(out)
        integer, intent(in) :: n
        integer, allocatable :: out(:,:)
        out = core_plainperms(n)
    end function plainperms

    function mset(v) result(out)
        integer, intent(in) :: v(:)
        integer, allocatable :: out(:,:)
        out = core_mset(v)
    end function mset

    function multiset(v, n) result(out)
        integer, intent(in) :: v(:)
        integer, intent(in), optional :: n
        integer, allocatable :: out(:,:)
        if (present(n)) then
            out = core_multiset(v, n)
        else
            out = core_multiset(v)
        end if
    end function multiset

    function multinomial(v) result(out)
        integer, intent(in) :: v(:)
        integer, allocatable :: out(:,:)
        out = core_multinomial(v)
    end function multinomial

    function allbinom(n, k) result(out)
        integer, intent(in) :: n, k
        integer, allocatable :: out(:,:)
        out = core_allbinom(n, k)
    end function allbinom

    function genrif(v) result(out)
        integer, intent(in) :: v(:)
        integer, allocatable :: out(:,:)
        out = core_genrif(v)
    end function genrif

    function riffle(p, q) result(out)
        integer, intent(in) :: p
        integer, intent(in), optional :: q
        integer, allocatable :: out(:,:)
        if (present(q)) then
            out = core_riffle(p, q)
        else
            out = core_riffle(p)
        end if
    end function riffle

    function tobin(n, len) result(out)
        integer, intent(in) :: n, len
        integer, allocatable :: out(:)
        out = core_tobin(int(n, i8), len)
    end function tobin

    function todec(bin) result(value)
        integer, intent(in) :: bin(:)
        integer(i8) :: value
        value = core_todec(bin)
    end function todec

    function comptobin(comp) result(out)
        integer, intent(in) :: comp(:)
        integer, allocatable :: out(:)
        out = core_comptobin(comp)
    end function comptobin

    function bintocomp(bin) result(out)
        integer, intent(in) :: bin(:)
        integer, allocatable :: out(:)
        out = core_bintocomp(bin)
    end function bintocomp

end module partitions_compat
