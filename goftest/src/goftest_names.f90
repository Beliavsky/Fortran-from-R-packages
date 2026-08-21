! SPDX-License-Identifier: GPL-2.0-or-later
module goftest_names
    implicit none
    private
    public :: recognise_cdf

contains

    function recognise_cdf(name) result(description)
        character(len=*), intent(in) :: name
        character(len=:), allocatable :: description
        character(len=:), allocatable :: root

        root = trim(adjustl(name))
        if (len(root) > 1) then
            if (root(1:1) == 'p') root = root(2:)
        end if
        select case (root)
        case ('beta')
            description = 'beta distribution'
        case ('binom')
            description = 'binomial distribution'
        case ('birthday')
            description = 'birthday coincidence distribution'
        case ('cauchy')
            description = 'Cauchy distribution'
        case ('chisq')
            description = 'chi-squared distribution'
        case ('exp')
            description = 'exponential distribution'
        case ('f')
            description = 'F distribution'
        case ('gamma')
            description = 'Gamma distribution'
        case ('geom')
            description = 'geometric distribution'
        case ('hyper')
            description = 'hypergeometric distribution'
        case ('lnorm')
            description = 'log-normal distribution'
        case ('logis')
            description = 'logistic distribution'
        case ('nbinom')
            description = 'negative binomial distribution'
        case ('norm')
            description = 'Normal distribution'
        case ('pois')
            description = 'Poisson distribution'
        case ('t')
            description = "Student's t distribution"
        case ('tukey')
            description = 'Tukey (Studentized range) distribution'
        case ('unif')
            description = 'uniform distribution'
        case ('weibull')
            description = 'Weibull distribution'
        case ('AD')
            description = 'null distribution of Anderson-Darling Test Statistic'
        case ('CvM')
            description = 'null distribution of Cramer-von Mises Test Statistic'
        case ('wilcox')
            description = 'null distribution of Wilcoxon Rank Sum Test Statistic'
        case default
            description = ''
        end select
    end function recognise_cdf

end module goftest_names
