! Translation of computational code from R package partitions 1.10-9.
! Upstream authors: Robin K. S. Hankin; contributor Paul Egeler.
! Upstream DESCRIPTION declares: License: GPL.
! This Fortran translation is distributed under the same GPL terms.

program demo
    use partitions
    implicit none

    integer, allocatable :: a(:,:), sp(:,:)
    integer :: j

    print '(a,i0)', 'P(100) = ', p(100)
    print '(a,i0)', 'Q(100) = ', q(100)
    print '(a,i0)', 'R(5,12) = ', r(5,12)

    a = parts(5)
    print '(a)', 'Partitions of 5 (one partition per column):'
    do j = 1, size(a,2)
        write(*,'(*(i0,1x))') a(:,j)
    end do

    sp = setparts([2,1,1])
    print '(a)', 'Set partitions with block sizes 2,1,1 (membership labels by column):'
    do j = 1, size(sp,2)
        write(*,'(*(i0,1x))') sp(:,j)
    end do
end program demo
