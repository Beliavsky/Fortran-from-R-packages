program test_names
    use goftest, only : recognise_cdf
    implicit none
    character(len=:), allocatable :: s

    s = recognise_cdf('punif')
    if (s /= 'uniform distribution') error stop 1
    s = recognise_cdf('AD')
    if (index(s, 'Anderson-Darling') == 0) error stop 1
    s = recognise_cdf('not_a_cdf')
    if (len(s) /= 0) error stop 1

    print '(a)', 'test_names: PASS'
end program test_names
