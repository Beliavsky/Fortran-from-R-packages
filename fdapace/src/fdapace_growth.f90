module fdapace_growth
    use fdapace_kinds, only : dp
    use fdapace_math, only : fmm_spline_eval
    implicit none
    private

    public :: make_bw_to_zscore_02y
    public :: make_hc_to_zscore_02y
    public :: make_ln_to_zscore_02y

contains

    function make_hc_to_zscore_02y(sex, age, hc) result(z)
        character(len=*), intent(in) :: sex !! Child sex code, exactly M or F ignoring surrounding blanks.
        real(dp), intent(in) :: age(:) !! Ages in months on the documented closed interval [0,24].
        real(dp), intent(in) :: hc(:) !! Head-circumference measurements in centimeters corresponding to age.
        real(dp), allocatable :: z(:)
        real(dp), parameter :: mu_f(25) = [ &
            33.8787_dp, 36.5463_dp, 38.2521_dp, 39.5328_dp, 40.5817_dp, 41.4590_dp, 42.1995_dp, &
            42.8290_dp, 43.3671_dp, 43.8300_dp, 44.2319_dp, 44.5844_dp, 44.8965_dp, 45.1752_dp, &
            45.4265_dp, 45.6551_dp, 45.8650_dp, 46.0598_dp, 46.2424_dp, 46.4152_dp, 46.5801_dp, &
            46.7384_dp, 46.8913_dp, 47.0391_dp, 47.1822_dp ]
        real(dp), parameter :: sd_f(25) = [ &
            1.18440_dp, 1.17314_dp, 1.21183_dp, 1.24133_dp, 1.26574_dp, 1.28606_dp, 1.30270_dp, &
            1.31699_dp, 1.32833_dp, 1.33813_dp, 1.34642_dp, 1.35314_dp, 1.35902_dp, 1.36384_dp, &
            1.36825_dp, 1.37239_dp, 1.37549_dp, 1.37857_dp, 1.38126_dp, 1.38410_dp, 1.38669_dp, &
            1.38907_dp, 1.39126_dp, 1.39330_dp, 1.39518_dp ]
        real(dp), parameter :: mu_m(25) = [ &
            34.4618_dp, 37.2759_dp, 39.1285_dp, 40.5135_dp, 41.6317_dp, 42.5576_dp, 43.3306_dp, &
            43.9803_dp, 44.5300_dp, 44.9998_dp, 45.4051_dp, 45.7573_dp, 46.0661_dp, 46.3395_dp, &
            46.5844_dp, 46.8060_dp, 47.0088_dp, 47.1962_dp, 47.3711_dp, 47.5357_dp, 47.6919_dp, &
            47.8408_dp, 47.9833_dp, 48.1201_dp, 48.2515_dp ]
        real(dp), parameter :: sd_m(25) = [ &
            1.27026_dp, 1.16785_dp, 1.17268_dp, 1.18218_dp, 1.19400_dp, 1.20736_dp, 1.22062_dp, &
            1.23321_dp, 1.24506_dp, 1.25639_dp, 1.26680_dp, 1.27617_dp, 1.28478_dp, 1.29241_dp, &
            1.30017_dp, 1.30682_dp, 1.31390_dp, 1.32008_dp, 1.32639_dp, 1.33243_dp, 1.33823_dp, &
            1.34433_dp, 1.34977_dp, 1.35554_dp, 1.36667_dp ]
        real(dp) :: time(25)
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: sd(:)
        integer :: i
        integer :: info

        call validate_growth_inputs(age, hc)
        do i = 1, 25
            time(i) = real(i - 1, dp)
        end do
        allocate(mu(size(age)), sd(size(age)), z(size(age)))
        select case (trim(sex))
        case ("F")
            call fmm_spline_eval(time, mu_f, age, mu, info)
            if (info /= 0) error stop "make_hc_to_zscore_02y: spline interpolation failed"
            call fmm_spline_eval(time, sd_f, age, sd, info)
        case ("M")
            call fmm_spline_eval(time, mu_m, age, mu, info)
            if (info /= 0) error stop "make_hc_to_zscore_02y: spline interpolation failed"
            call fmm_spline_eval(time, sd_m, age, sd, info)
        case default
            error stop "make_hc_to_zscore_02y: sex must be M or F"
        end select
        if (info /= 0) error stop "make_hc_to_zscore_02y: spline interpolation failed"
        z = (hc - mu) / sd
    end function make_hc_to_zscore_02y

    function make_ln_to_zscore_02y(sex, age, ln) result(z)
        character(len=*), intent(in) :: sex !! Child sex code, exactly M or F ignoring surrounding blanks.
        real(dp), intent(in) :: age(:) !! Ages in months on the documented closed interval [0,24].
        real(dp), intent(in) :: ln(:) !! Body-length measurements in centimeters corresponding to age.
        real(dp), allocatable :: z(:)
        real(dp), parameter :: mu_f(25) = [ &
            49.1477_dp, 53.6872_dp, 57.0673_dp, 59.8029_dp, 62.0899_dp, 64.0301_dp, 65.7311_dp, &
            67.2873_dp, 68.7498_dp, 70.1435_dp, 71.4818_dp, 72.7710_dp, 74.0150_dp, 75.2176_dp, &
            76.3817_dp, 77.5099_dp, 78.6055_dp, 79.6710_dp, 80.7079_dp, 81.7182_dp, 82.7036_dp, &
            83.6654_dp, 84.6040_dp, 85.5202_dp, 86.4153_dp ]
        real(dp), parameter :: sd_f(25) = [ &
            1.8627_dp, 1.9542_dp, 2.0362_dp, 2.1051_dp, 2.1645_dp, 2.2174_dp, 2.2664_dp, &
            2.3154_dp, 2.3650_dp, 2.4157_dp, 2.4676_dp, 2.5208_dp, 2.5750_dp, 2.6296_dp, &
            2.6841_dp, 2.7392_dp, 2.7944_dp, 2.8490_dp, 2.9039_dp, 2.9582_dp, 3.0129_dp, &
            3.0672_dp, 3.1202_dp, 3.1737_dp, 3.2267_dp ]
        real(dp), parameter :: mu_m(25) = [ &
            49.8842_dp, 54.7244_dp, 58.4249_dp, 61.4292_dp, 63.8860_dp, 65.9026_dp, 67.6236_dp, &
            69.1645_dp, 70.5994_dp, 71.9687_dp, 73.2810_dp, 74.5388_dp, 75.7488_dp, 76.9186_dp, &
            78.0497_dp, 79.1458_dp, 80.2110_dp, 81.2487_dp, 82.2587_dp, 83.2418_dp, 84.1996_dp, &
            85.1348_dp, 6.0477_dp, 86.9410_dp, 87.8161_dp ]
        real(dp), parameter :: sd_m(25) = [ &
            1.8931_dp, 1.9465_dp, 2.0005_dp, 2.0444_dp, 2.0808_dp, 2.1115_dp, 2.1403_dp, &
            2.1711_dp, 2.2055_dp, 2.2433_dp, 2.2849_dp, 2.3293_dp, 2.3762_dp, 2.4260_dp, &
            2.4773_dp, 2.5303_dp, 2.5844_dp, 2.6406_dp, 2.6973_dp, 2.7553_dp, 2.8140_dp, &
            2.8742_dp, 2.9342_dp, 2.9951_dp, 3.0551_dp ]
        real(dp) :: time(25)
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: sd(:)
        integer :: i
        integer :: info

        call validate_growth_inputs(age, ln)
        do i = 1, 25
            time(i) = real(i - 1, dp)
        end do
        allocate(mu(size(age)), sd(size(age)), z(size(age)))
        select case (trim(sex))
        case ("F")
            call fmm_spline_eval(time, mu_f, age, mu, info)
            if (info /= 0) error stop "make_ln_to_zscore_02y: spline interpolation failed"
            call fmm_spline_eval(time, sd_f, age, sd, info)
        case ("M")
            call fmm_spline_eval(time, mu_m, age, mu, info)
            if (info /= 0) error stop "make_ln_to_zscore_02y: spline interpolation failed"
            call fmm_spline_eval(time, sd_m, age, sd, info)
        case default
            error stop "make_ln_to_zscore_02y: sex must be M or F"
        end select
        if (info /= 0) error stop "make_ln_to_zscore_02y: spline interpolation failed"
        z = (ln - mu) / sd
    end function make_ln_to_zscore_02y

    function make_bw_to_zscore_02y(sex, age, bw) result(z)
        character(len=*), intent(in) :: sex !! Child sex code, exactly M or F ignoring surrounding blanks.
        real(dp), intent(in) :: age(:) !! Ages in months on the documented closed interval [0,24].
        real(dp), intent(in) :: bw(:) !! Body-weight measurements corresponding one-to-one with age.
        real(dp), allocatable :: z(:)
        real(dp), parameter :: s_f(25) = [ &
            0.14171_dp, 0.13724_dp, 0.13000_dp, 0.12619_dp, 0.12402_dp, 0.12274_dp, 0.12204_dp, &
            0.12178_dp, 0.12181_dp, 0.12199_dp, 0.12223_dp, 0.12247_dp, 0.12268_dp, 0.12283_dp, &
            0.12294_dp, 0.12299_dp, 0.12303_dp, 0.12306_dp, 0.12309_dp, 0.12315_dp, 0.12323_dp, &
            0.12335_dp, 0.12350_dp, 0.12369_dp, 0.12390_dp ]
        real(dp), parameter :: m_f(25) = [ &
            3.2322_dp, 4.1873_dp, 5.1282_dp, 5.8458_dp, 6.4237_dp, 6.8985_dp, 7.2970_dp, &
            7.6422_dp, 7.9487_dp, 8.2254_dp, 8.4800_dp, 8.7192_dp, 8.9481_dp, 9.1699_dp, &
            9.3870_dp, 9.6008_dp, 9.8124_dp, 10.0226_dp, 10.2315_dp, 10.4393_dp, 10.6464_dp, &
            10.8534_dp, 11.0608_dp, 11.2688_dp, 11.4775_dp ]
        real(dp), parameter :: l_f(25) = [ &
            0.3809_dp, 0.1714_dp, 0.0962_dp, 0.0402_dp, -0.0050_dp, -0.0430_dp, -0.0756_dp, &
            -0.1039_dp, -0.1288_dp, -0.1507_dp, -0.1700_dp, -0.1872_dp, -0.2024_dp, -0.2158_dp, &
            -0.2278_dp, -0.2384_dp, -0.2478_dp, -0.2562_dp, -0.2637_dp, -0.2703_dp, -0.2762_dp, &
            -0.2815_dp, -0.2862_dp, -0.2903_dp, -0.2941_dp ]
        real(dp), parameter :: s_m(25) = [ &
            0.14602_dp, 0.13395_dp, 0.12385_dp, 0.11727_dp, 0.11316_dp, 0.11080_dp, 0.10958_dp, &
            0.10902_dp, 0.10882_dp, 0.10881_dp, 0.10891_dp, 0.10906_dp, 0.10925_dp, 0.10949_dp, &
            0.10976_dp, 0.11007_dp, 0.11041_dp, 0.11079_dp, 0.11119_dp, 0.11164_dp, 0.11211_dp, &
            0.11261_dp, 0.11314_dp, 0.11369_dp, 0.11426_dp ]
        real(dp), parameter :: m_m(25) = [ &
            3.3464_dp, 4.4709_dp, 5.5675_dp, 6.3762_dp, 7.0023_dp, 7.5105_dp, 7.9340_dp, &
            8.2970_dp, 8.6151_dp, 8.9014_dp, 9.1649_dp, 9.4122_dp, 9.6479_dp, 9.8749_dp, &
            10.0953_dp, 10.3108_dp, 10.5228_dp, 10.7319_dp, 10.9385_dp, 11.1430_dp, 11.3462_dp, &
            11.5486_dp, 11.7504_dp, 11.9514_dp, 12.1515_dp ]
        real(dp), parameter :: l_m(25) = [ &
            0.3487_dp, 0.2297_dp, 0.1970_dp, 0.1738_dp, 0.1553_dp, 0.1395_dp, 0.1257_dp, &
            0.1134_dp, 0.1021_dp, 0.0917_dp, 0.0820_dp, 0.0730_dp, 0.0644_dp, 0.0563_dp, &
            0.0487_dp, 0.0413_dp, 0.0343_dp, 0.0275_dp, 0.0211_dp, 0.0148_dp, 0.0087_dp, &
            0.0029_dp, -0.0028_dp, -0.0083_dp, -0.0137_dp ]
        real(dp) :: time(25)
        real(dp), allocatable :: lval(:)
        real(dp), allocatable :: mval(:)
        real(dp), allocatable :: sval(:)
        integer :: i
        integer :: info

        call validate_growth_inputs(age, bw)
        do i = 1, 25
            time(i) = real(i - 1, dp)
        end do
        allocate(lval(size(age)), mval(size(age)), sval(size(age)), z(size(age)))
        select case (trim(sex))
        case ("F")
            call fmm_spline_eval(time, l_f, age, lval, info)
            if (info /= 0) error stop "make_bw_to_zscore_02y: spline interpolation failed"
            call fmm_spline_eval(time, m_f, age, mval, info)
            if (info /= 0) error stop "make_bw_to_zscore_02y: spline interpolation failed"
            call fmm_spline_eval(time, s_f, age, sval, info)
        case ("M")
            call fmm_spline_eval(time, l_m, age, lval, info)
            if (info /= 0) error stop "make_bw_to_zscore_02y: spline interpolation failed"
            call fmm_spline_eval(time, m_m, age, mval, info)
            if (info /= 0) error stop "make_bw_to_zscore_02y: spline interpolation failed"
            call fmm_spline_eval(time, s_m, age, sval, info)
        case default
            error stop "make_bw_to_zscore_02y: sex must be M or F"
        end select
        if (info /= 0) error stop "make_bw_to_zscore_02y: spline interpolation failed"
        do i = 1, size(age)
            z(i) = ((bw(i) / mval(i))**lval(i) - 1.0_dp) / (sval(i) * lval(i))
        end do
    end function make_bw_to_zscore_02y

    pure subroutine validate_growth_inputs(age, measurement)
        real(dp), intent(in) :: age(:) !! Ages to validate for equal length and the documented [0,24] month support.
        real(dp), intent(in) :: measurement(:) !! Measurements that must correspond one-to-one with age.

        if (size(age) /= size(measurement)) error stop "growth z-score: age and measurement lengths differ"
        if (size(age) > 0) then
            if (minval(age) < 0.0_dp .or. maxval(age) > 24.0_dp) error stop "growth z-score: age is outside [0,24] months"
        end if
    end subroutine validate_growth_inputs

end module fdapace_growth
