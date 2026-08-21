! SPDX-License-Identifier: GPL-2.0-or-later
module desolve_rk
  use desolve_kinds, only : dp
  use desolve_types, only : ode_rhs, ode_result
  implicit none
  private

  type, public :: rk_method
    character(len=16) :: id = ''
    integer :: stages = 0
    integer :: qerr = 1
    logical :: adaptive = .false.
    logical :: implicit = .false.
    real(dp), allocatable :: a(:,:), b_low(:), b_high(:), c(:)
  end type rk_method

  public :: rk_method_by_name, rk_integrate, euler, rk4
contains

  function rk_method_by_name(name) result(m)
    character(len=*), intent(in) :: name
    type(rk_method) :: m
    character(len=:), allocatable :: key
    key = lower(trim(name))
    if (key == 'ode23') key = 'rk23bs'
    if (key == 'ode45') key = 'rk45dp7'
    select case(key)
    case ('euler')
      call allocate_method(m, 'euler', 1, .false., .false., 1)
      m%b_low = [1.0_dp]
      m%b_high = m%b_low
      m%c = [0.0_dp]
    case ('rk2')
      call allocate_method(m, 'rk2', 2, .false., .false., 1)
      m%A(2,1) = 1.0_dp
      m%b_low = [0.5_dp, 0.5_dp]
      m%b_high = m%b_low
      m%c = [0.0_dp, 1.0_dp]
    case ('rk4')
      call allocate_method(m, 'rk4', 4, .false., .false., 4)
      m%A(2,1) = 0.5_dp
      m%A(3,2) = 0.5_dp
      m%A(4,3) = 1.0_dp
      m%b_low = [ &
          0.16666666666666666_dp, 0.33333333333333331_dp, &
          0.33333333333333331_dp, 0.16666666666666666_dp ]
      m%b_high = m%b_low
      m%c = [ &
          0.0_dp, 0.5_dp, &
          0.5_dp, 1.0_dp ]
    case ('rk23')
      call allocate_method(m, 'rk23', 3, .true., .false., 2)
      call set_rowmajor(m%A, 3, 3, [ &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.5_dp, &
          0.0_dp, 0.0_dp, &
          -1.0_dp, 2.0_dp, &
          0.0_dp ])
      m%b_low = [ &
          0.0_dp, 1.0_dp, &
          0.0_dp ]
      m%b_high = [ &
          0.16666666666666666_dp, 0.66666666666666663_dp, &
          0.16666666666666666_dp ]
      m%c = [ &
          0.0_dp, 0.5_dp, &
          2.0_dp ]
    case ('rk23bs')
      call allocate_method(m, 'rk23bs', 4, .true., .false., 2)
      call set_rowmajor(m%A, 4, 4, [ &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.5_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.75_dp, &
          0.0_dp, 0.0_dp, &
          0.22222222222222221_dp, 0.33333333333333331_dp, &
          0.44444444444444442_dp, 0.0_dp ])
      m%b_low = [ &
          0.29166666666666669_dp, 0.25_dp, &
          0.33333333333333331_dp, 0.125_dp ]
      m%b_high = [ &
          0.22222222222222221_dp, 0.33333333333333331_dp, &
          0.44444444444444442_dp, 0.0_dp ]
      m%c = [ &
          0.0_dp, 0.5_dp, &
          0.75_dp, 1.0_dp ]
    case ('rk34f')
      call allocate_method(m, 'rk34f', 5, .true., .false., 3)
      call set_rowmajor(m%A, 5, 4, [ &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.2857142857142857_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.085555555555555551_dp, 0.38111111111111112_dp, &
          0.0_dp, 0.0_dp, &
          0.55747922437673125_dp, -1.4064550225980463_dp, &
          1.7700284298002624_dp, 0.0_dp, &
          0.16122448979591836_dp, 0.0_dp, &
          0.59983452840595697_dp, 0.23894098179812465_dp ])
      m%b_low = [ &
          0.16122448979591836_dp, 0.0_dp, &
          0.59983452840595697_dp, 0.23894098179812465_dp, &
          0.0_dp ]
      m%b_high = [ &
          0.15578231292517006_dp, 0.0_dp, &
          0.62051847766133483_dp, 0.16814365385793958_dp, &
          0.055555555555555552_dp ]
      m%c = [ &
          0.0_dp, 0.2857142857142857_dp, &
          0.46666666666666667_dp, 0.92105263157894735_dp, &
          1.0_dp ]
    case ('rk45f')
      call allocate_method(m, 'rk45f', 6, .true., .false., 4)
      call set_rowmajor(m%A, 6, 5, [ &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.25_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.09375_dp, 0.28125_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.87938097405553028_dp, &
          -3.2771961766044608_dp, 3.3208921256258535_dp, &
          0.0_dp, 0.0_dp, &
          2.0324074074074074_dp, -8.0_dp, &
          7.1734892787524362_dp, -0.20589668615984405_dp, &
          0.0_dp, -0.29629629629629628_dp, &
          2.0_dp, -1.3816764132553607_dp, &
          0.45297270955165692_dp, -0.27500000000000002_dp ])
      m%b_low = [ &
          0.11574074074074074_dp, 0.0_dp, &
          0.54892787524366471_dp, 0.53533138401559455_dp, &
          -0.20000000000000001_dp, 0.0_dp ]
      m%b_high = [ &
          0.11851851851851852_dp, 0.0_dp, &
          0.51898635477582844_dp, 0.50613149034201665_dp, &
          -0.17999999999999999_dp, 0.036363636363636362_dp ]
      m%c = [ &
          0.0_dp, 0.25_dp, &
          0.375_dp, 0.92307692307692313_dp, &
          1.0_dp, 0.5_dp ]
    case ('rk45ck')
      call allocate_method(m, 'rk45ck', 6, .true., .false., 4)
      call set_rowmajor(m%A, 6, 5, [ &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.20000000000000001_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.074999999999999997_dp, 0.22500000000000001_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.29999999999999999_dp, &
          -0.90000000000000002_dp, 1.2_dp, &
          0.0_dp, 0.0_dp, &
          -0.20370370370370369_dp, 2.5_dp, &
          -2.5925925925925926_dp, 1.2962962962962963_dp, &
          0.0_dp, 0.029495804398148147_dp, &
          0.341796875_dp, 0.041594328703703706_dp, &
          0.40034541377314814_dp, 0.061767578125_dp ])
      m%b_low = [ &
          0.10217737268518519_dp, 0.0_dp, &
          0.38390790343915343_dp, 0.24459273726851852_dp, &
          0.019321986607142856_dp, 0.25_dp ]
      m%b_high = [ &
          0.097883597883597878_dp, 0.0_dp, &
          0.40257648953301128_dp, 0.21043771043771045_dp, &
          0.0_dp, 0.28910220214568039_dp ]
      m%c = [ &
          0.0_dp, 0.20000000000000001_dp, &
          0.29999999999999999_dp, 0.59999999999999998_dp, &
          1.0_dp, 0.875_dp ]
    case ('rk45e')
      call allocate_method(m, 'rk45e', 6, .true., .false., 4)
      call set_rowmajor(m%A, 6, 5, [ &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.5_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.25_dp, 0.25_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          -1.0_dp, 2.0_dp, &
          0.0_dp, 0.0_dp, &
          0.25925925925925924_dp, 0.37037037037037035_dp, &
          0.0_dp, 0.037037037037037035_dp, &
          0.0_dp, 0.0448_dp, &
          -0.20000000000000001_dp, 0.87360000000000004_dp, &
          0.086400000000000005_dp, -0.6048_dp ])
      m%b_low = [ &
          0.16666666666666666_dp, 0.0_dp, &
          0.66666666666666663_dp, 0.16666666666666666_dp, &
          0.0_dp, 0.0_dp ]
      m%b_high = [ &
          0.041666666666666664_dp, 0.0_dp, &
          0.0_dp, 0.10416666666666667_dp, &
          0.48214285714285715_dp, 0.37202380952380953_dp ]
      m%c = [ &
          0.0_dp, 0.5_dp, &
          0.5_dp, 1.0_dp, &
          0.66666666666666663_dp, 0.20000000000000001_dp ]
    case ('rk45dp6')
      call allocate_method(m, 'rk45dp6', 6, .true., .false., 4)
      call set_rowmajor(m%A, 6, 5, [ &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.20000000000000001_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.074999999999999997_dp, 0.22500000000000001_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.29999999999999999_dp, &
          -0.90000000000000002_dp, 1.2_dp, &
          0.0_dp, 0.0_dp, &
          0.31001371742112482_dp, -0.92592592592592593_dp, &
          1.2071330589849107_dp, 0.075445816186556922_dp, &
          0.0_dp, -0.67037037037037039_dp, &
          2.5_dp, -0.89562289562289565_dp, &
          -3.3703703703703702_dp, 3.4363636363636365_dp ])
      m%b_low = [ &
          0.057407407407407407_dp, 0.0_dp, &
          0.63973063973063971_dp, -1.3425925925925926_dp, &
          1.5954545454545455_dp, 0.050000000000000003_dp ]
      m%b_high = [ &
          0.087962962962962965_dp, 0.0_dp, &
          0.48100048100048098_dp, -0.57870370370370372_dp, &
          0.92045454545454541_dp, 0.089285714285714288_dp ]
      m%c = [ &
          0.0_dp, 0.20000000000000001_dp, &
          0.29999999999999999_dp, 0.59999999999999998_dp, &
          0.66666666666666663_dp, 1.0_dp ]
    case ('rk45dp7')
      call allocate_method(m, 'rk45dp7', 7, .true., .false., 4)
      call set_rowmajor(m%A, 7, 6, [ &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.20000000000000001_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.074999999999999997_dp, 0.22500000000000001_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.97777777777777775_dp, -3.7333333333333334_dp, &
          3.5555555555555554_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          2.9525986892242035_dp, -11.595793324188385_dp, &
          9.8228928516994358_dp, -0.29080932784636487_dp, &
          0.0_dp, 0.0_dp, &
          2.8462752525252526_dp, -10.757575757575758_dp, &
          8.9064227177434727_dp, 0.27840909090909088_dp, &
          -0.2735313036020583_dp, 0.0_dp, &
          0.091145833333333329_dp, 0.0_dp, &
          0.44923629829290207_dp, 0.65104166666666663_dp, &
          -0.322376179245283_dp, 0.13095238095238096_dp ])
      m%b_low = [ &
          0.089913194444444441_dp, 0.0_dp, &
          0.45348906858340821_dp, 0.61406249999999996_dp, &
          -0.27151238207547168_dp, 0.089047619047619042_dp, &
          0.025000000000000001_dp ]
      m%b_high = [ &
          0.091145833333333329_dp, 0.0_dp, &
          0.44923629829290207_dp, 0.65104166666666663_dp, &
          -0.322376179245283_dp, 0.13095238095238096_dp, &
          0.0_dp ]
      m%c = [ &
          0.0_dp, 0.20000000000000001_dp, &
          0.29999999999999999_dp, 0.80000000000000004_dp, &
          0.88888888888888884_dp, 1.0_dp, &
          1.0_dp ]
    case ('rk78dp')
      call allocate_method(m, 'rk78dp', 13, .true., .false., 7)
      call set_rowmajor(m%A, 13, 12, [ &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.055555555555555552_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.020833333333333332_dp, 0.0625_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.03125_dp, 0.0_dp, &
          0.09375_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.3125_dp, 0.0_dp, &
          -1.171875_dp, 1.171875_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.037499999999999999_dp, 0.0_dp, &
          0.0_dp, 0.1875_dp, &
          0.14999999999999999_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.047910137111111112_dp, 0.0_dp, &
          0.0_dp, 0.11224871277777777_dp, &
          -0.025505673777777779_dp, 0.012846823888888888_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.016917989787292281_dp, 0.0_dp, &
          0.0_dp, 0.3878482784860432_dp, &
          0.035977369851500331_dp, 0.19697021421566607_dp, &
          -0.17271385234050185_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.069095753359192297_dp, 0.0_dp, &
          0.0_dp, -0.63424797672885413_dp, &
          -0.16119757522460407_dp, 0.13865030945882525_dp, &
          0.94092861403575623_dp, 0.21163632648194397_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.18355699683904539_dp, 0.0_dp, &
          0.0_dp, -2.4687680843155926_dp, &
          -0.29128688781630047_dp, -0.026473020233117376_dp, &
          2.8478387641928005_dp, 0.28138733146984979_dp, &
          0.12374489986331466_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          -1.2154248173958881_dp, 0.0_dp, &
          0.0_dp, 16.672608665945774_dp, &
          0.91574182841681795_dp, -6.0566058043574706_dp, &
          -16.00357359415618_dp, 14.849303086297663_dp, &
          -13.371575735289849_dp, 5.134182648179638_dp, &
          0.0_dp, 0.0_dp, &
          0.25886091643826425_dp, 0.0_dp, &
          0.0_dp, -4.7744857854892047_dp, &
          -0.43509301377703252_dp, -3.0494833320722416_dp, &
          5.5779200399360995_dp, 6.1558315898610401_dp, &
          -5.0621045867369387_dp, 2.193926173180679_dp, &
          0.13462799865933495_dp, 0.0_dp, &
          0.82242759962650747_dp, 0.0_dp, &
          0.0_dp, -11.658673257277664_dp, &
          -0.75762211669093615_dp, 0.71397358815958156_dp, &
          12.075774986890057_dp, -2.1276591139204029_dp, &
          1.9901662070489554_dp, -0.23428647154404028_dp, &
          0.17589857770794226_dp, 0.0_dp ])
      m%b_low = [ &
          0.029553213676353499_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, -0.82860627648779706_dp, &
          0.31124090005111832_dp, 2.4673451905998869_dp, &
          -2.5469416518419088_dp, 1.4435485836767752_dp, &
          0.079415595881127288_dp, 0.044444444444444446_dp, &
          0.0_dp ]
      m%b_high = [ &
          0.041747491141530244_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, -0.055452328611239311_dp, &
          0.23931280720118009_dp, 0.70351066940344298_dp, &
          -0.75975961381446089_dp, 0.6605630309222863_dp, &
          0.15818748251012332_dp, -0.23810953875286281_dp, &
          0.25_dp ]
      m%c = [ &
          0.0_dp, 0.055555555555555552_dp, &
          0.083333333333333329_dp, 0.125_dp, &
          0.3125_dp, 0.375_dp, &
          0.14749999999999999_dp, 0.46500000000000002_dp, &
          0.56486545138225952_dp, 0.65000000000000002_dp, &
          0.9246562776405044_dp, 1.0_dp, &
          1.0_dp ]
    case ('rk78f')
      call allocate_method(m, 'rk78f', 13, .true., .false., 7)
      call set_rowmajor(m%A, 13, 12, [ &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.07407407407407407_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.027777777777777776_dp, 0.083333333333333329_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.041666666666666664_dp, 0.0_dp, &
          0.125_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.41666666666666669_dp, 0.0_dp, &
          -1.5625_dp, 1.5625_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.050000000000000003_dp, 0.0_dp, &
          0.0_dp, 0.25_dp, &
          0.20000000000000001_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          -0.23148148148148148_dp, 0.0_dp, &
          0.0_dp, 1.1574074074074074_dp, &
          -2.4074074074074074_dp, 2.3148148148148149_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.10333333333333333_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.27111111111111114_dp, -0.22222222222222221_dp, &
          0.014444444444444444_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          2.0_dp, 0.0_dp, &
          0.0_dp, -8.8333333333333339_dp, &
          15.644444444444444_dp, -11.888888888888889_dp, &
          0.74444444444444446_dp, 3.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          -0.84259259259259256_dp, 0.0_dp, &
          0.0_dp, 0.21296296296296297_dp, &
          -7.2296296296296294_dp, 5.7592592592592595_dp, &
          -0.31666666666666665_dp, 2.8333333333333335_dp, &
          -0.083333333333333329_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.58121951219512191_dp, 0.0_dp, &
          0.0_dp, -2.0792682926829267_dp, &
          4.3863414634146345_dp, -3.6707317073170733_dp, &
          0.52024390243902441_dp, 0.54878048780487809_dp, &
          0.27439024390243905_dp, 0.43902439024390244_dp, &
          0.0_dp, 0.0_dp, &
          0.014634146341463415_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, -0.14634146341463414_dp, &
          -0.014634146341463415_dp, -0.073170731707317069_dp, &
          0.073170731707317069_dp, 0.14634146341463414_dp, &
          0.0_dp, 0.0_dp, &
          -0.43341463414634146_dp, 0.0_dp, &
          0.0_dp, -2.0792682926829267_dp, &
          4.3863414634146345_dp, -3.524390243902439_dp, &
          0.53487804878048784_dp, 0.62195121951219512_dp, &
          0.20121951219512196_dp, 0.29268292682926828_dp, &
          0.0_dp, 1.0_dp ])
      m%b_low = [ &
          0.04880952380952381_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.32380952380952382_dp, &
          0.25714285714285712_dp, 0.25714285714285712_dp, &
          0.03214285714285714_dp, 0.03214285714285714_dp, &
          0.04880952380952381_dp, 0.0_dp, &
          0.0_dp ]
      m%b_high = [ &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.32380952380952382_dp, &
          0.25714285714285712_dp, 0.25714285714285712_dp, &
          0.03214285714285714_dp, 0.03214285714285714_dp, &
          0.0_dp, 0.04880952380952381_dp, &
          0.04880952380952381_dp ]
      m%c = [ &
          0.0_dp, 0.07407407407407407_dp, &
          0.1111111111111111_dp, 0.16666666666666666_dp, &
          0.41666666666666669_dp, 0.5_dp, &
          0.83333333333333337_dp, 0.16666666666666666_dp, &
          0.66666666666666663_dp, 0.33333333333333331_dp, &
          1.0_dp, 0.0_dp, &
          1.0_dp ]
    case ('irk3r')
      call allocate_method(m, 'irk3r', 2, .false., .true., 3)
      call set_rowmajor(m%A, 2, 2, [ &
          0.41666666666666669_dp, -0.083333333333333329_dp, &
          0.75_dp, 0.25_dp ])
      m%b_low = [0.75_dp, 0.25_dp]
      m%b_high = m%b_low
      m%c = [0.33333333333333331_dp, 0.25_dp]
    case ('irk5r')
      call allocate_method(m, 'irk5r', 3, .false., .true., 5)
      call set_rowmajor(m%A, 3, 3, [ &
          0.19681547722366044_dp, -0.065535425850198378_dp, &
          0.023770974348220151_dp, 0.39442431473908729_dp, &
          0.29207341166522843_dp, -0.041548752125997922_dp, &
          0.37640306270046725_dp, 0.51248582618842164_dp, &
          0.1111111111111111_dp ])
      m%b_low = [ &
          0.37640306270046725_dp, 0.51248582618842164_dp, &
          0.1111111111111111_dp ]
      m%b_high = m%b_low
      m%c = [ &
          0.15505102572168222_dp, 0.64494897427831788_dp, &
          1.0_dp ]
    case ('irk4hh')
      call allocate_method(m, 'irk4hh', 2, .false., .true., 4)
      call set_rowmajor(m%A, 2, 2, [ &
          0.25_dp, -0.038675134594812866_dp, &
          0.53867513459481287_dp, 0.25_dp ])
      m%b_low = [0.5_dp, 0.5_dp]
      m%b_high = m%b_low
      m%c = [0.21132486540518713_dp, 0.78867513459481287_dp]
    case ('irk6kb')
      call allocate_method(m, 'irk6kb', 3, .false., .true., 6)
      call set_rowmajor(m%A, 3, 3, [ &
          0.1388888888888889_dp, -0.035976667524938943_dp, &
          0.0097894440153083184_dp, 0.30026319498086462_dp, &
          0.22222222222222221_dp, -0.022485417203086805_dp, &
          0.26798833376246944_dp, 0.48042111196938336_dp, &
          0.1388888888888889_dp ])
      m%b_low = [ &
          0.27777777777777779_dp, 0.44444444444444442_dp, &
          0.27777777777777779_dp ]
      m%b_high = m%b_low
      m%c = [ &
          0.1127016653792583_dp, 0.5_dp, &
          0.8872983346207417_dp ]
    case ('irk4l')
      call allocate_method(m, 'irk4l', 3, .false., .true., 4)
      call set_rowmajor(m%A, 3, 3, [ &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.25_dp, &
          0.25_dp, 0.0_dp, &
          0.0_dp, 1.0_dp, &
          0.0_dp ])
      m%b_low = [ &
          0.16666666666666666_dp, 0.66666666666666663_dp, &
          0.16666666666666666_dp ]
      m%b_high = m%b_low
      m%c = [ &
          0.0_dp, 0.5_dp, &
          1.0_dp ]
    case ('irk6l')
      call allocate_method(m, 'irk6l', 4, .false., .true., 6)
      call set_rowmajor(m%A, 4, 4, [ &
          0.0_dp, 0.0_dp, &
          0.0_dp, 0.0_dp, &
          0.12060113295832983_dp, 0.16666666666666666_dp, &
          -0.010874597374975477_dp, 0.0_dp, &
          0.046065533708336839_dp, 0.51087459737497543_dp, &
          0.16666666666666666_dp, 0.0_dp, &
          0.16666666666666666_dp, 0.23032766854168418_dp, &
          0.60300566479164919_dp, 0.0_dp ])
      m%b_low = [ &
          0.083333333333333329_dp, 0.41666666666666669_dp, &
          0.41666666666666669_dp, 0.083333333333333329_dp ]
      m%b_high = m%b_low
      m%c = [ &
          0.0_dp, 0.27639320225002101_dp, &
          0.72360679774997894_dp, 1.0_dp ]
    case default
      error stop 'deSolve rk_method_by_name: unknown Runge-Kutta method'
    end select
  end function rk_method_by_name

  function rk_integrate(rhs, y0, times, method, h, rtol, atol, max_steps) result(sol)
    procedure(ode_rhs) :: rhs
    real(dp), intent(in) :: y0(:), times(:)
    type(rk_method), intent(in) :: method
    real(dp), intent(in), optional :: h, rtol, atol
    integer, intent(in), optional :: max_steps
    type(ode_result) :: sol
    real(dp), allocatable :: y(:), ynext(:)
    real(dp) :: t, target, hs, hbase, rt, at, direction
    integer :: i, nstep, maxst, n
    logical :: ok

    n=size(y0)
    allocate(sol%t(size(times)),sol%y(n,size(times)),y(n),ynext(n))
    sol%t=times; sol%y=0.0_dp; sol%y(:,1)=y0; y=y0
    sol%status=2; sol%message='success'; sol%stats=solver_stats_zero()
    if (.not.valid_times(times)) then
      sol%status=-100;sol%message='times must be strictly monotone';return
    end if
    if (method%stages < 1) then
      sol%status=-101;sol%message='invalid Runge-Kutta method';return
    end if
    rt=1e-6_dp;if(present(rtol))rt=rtol
    at=1e-8_dp;if(present(atol))at=atol
    maxst=100000;if(present(max_steps))maxst=max_steps
    direction=merge(1.0_dp,-1.0_dp,times(size(times))>times(1))
    if(present(h))then;hbase=direction*abs(h)
    else;hbase=(times(size(times))-times(1))/max(100.0_dp,10.0_dp*real(size(times)-1,dp));end if
    if(abs(hbase)<=tiny(1.0_dp))hbase=direction*1e-6_dp
    t=times(1);hs=hbase;nstep=0
    do i=2,size(times)
      target=times(i)
      do while(direction*(target-t)>100.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(target)))
        if(nstep>=maxst)then;sol%status=-1;sol%message='maximum Runge-Kutta steps exceeded';return;end if
        if(direction*(t+hs-target)>0.0_dp)hs=target-t
        if(method%implicit)then
          call implicit_step(rhs,t,y,hs,method,ynext,ok)
          if (ok) y=ynext
          if(.not.ok)then;sol%status=-5;sol%message='implicit Runge-Kutta Newton iteration failed';return;end if
          t=t+hs;nstep=nstep+1;sol%stats%n_rhs=sol%stats%n_rhs+method%stages
          hs=sign(min(abs(hbase),abs(target-t)+abs(hbase)),direction)
        else if(method%adaptive)then
          call adaptive_step(rhs,t,y,hs,method,rt,at,target,direction,ok,sol%stats%n_rhs)
          nstep=nstep+1
          if(.not.ok)cycle
        else
          call explicit_step(rhs,t,y,hs,method,ynext,sol%stats%n_rhs)
          y=ynext
          t=t+hs;nstep=nstep+1;hs=hbase
        end if
      end do
      t=target;sol%y(:,i)=y
    end do
    sol%stats%n_steps=nstep
    sol%stats%step_last=hs
  end function rk_integrate

  function euler(rhs,y0,times,h) result(sol)
    procedure(ode_rhs)::rhs
    real(dp),intent(in)::y0(:),times(:)
    real(dp),intent(in),optional::h
    type(ode_result)::sol
    type(rk_method)::m
    m=rk_method_by_name('euler')
    if(present(h))then;sol=rk_integrate(rhs,y0,times,m,h=h)
    else;sol=rk_integrate(rhs,y0,times,m)
    end if
  end function euler

  function rk4(rhs,y0,times,h) result(sol)
    procedure(ode_rhs)::rhs
    real(dp),intent(in)::y0(:),times(:)
    real(dp),intent(in),optional::h
    type(ode_result)::sol
    type(rk_method)::m
    m=rk_method_by_name('rk4')
    if(present(h))then;sol=rk_integrate(rhs,y0,times,m,h=h)
    else;sol=rk_integrate(rhs,y0,times,m)
    end if
  end function rk4

  subroutine explicit_step(rhs,t,y,h,m,yout,nrhs)
    procedure(ode_rhs)::rhs
    real(dp),intent(in)::t,y(:),h
    type(rk_method),intent(in)::m
    real(dp),intent(out)::yout(:)
    integer,intent(inout)::nrhs
    real(dp),allocatable::k(:,:),yt(:)
    integer::i,j
    allocate(k(size(y),m%stages),yt(size(y)));k=0.0_dp
    do i=1,m%stages
      yt=y
      do j=1,i-1;yt=yt+h*m%a(i,j)*k(:,j);end do
      call rhs(t+h*m%c(i),yt,k(:,i));nrhs=nrhs+1
    end do
    yout=y
    do i=1,m%stages;yout=yout+h*m%b_high(i)*k(:,i);end do
  end subroutine explicit_step

  subroutine embedded_pair(rhs,t,y,h,m,ylow,yhigh,nrhs)
    procedure(ode_rhs)::rhs
    real(dp),intent(in)::t,y(:),h
    type(rk_method),intent(in)::m
    real(dp),intent(out)::ylow(:),yhigh(:)
    integer,intent(inout)::nrhs
    real(dp),allocatable::k(:,:),yt(:)
    integer::i,j
    allocate(k(size(y),m%stages),yt(size(y)));k=0.0_dp
    do i=1,m%stages
      yt=y
      do j=1,i-1;yt=yt+h*m%a(i,j)*k(:,j);end do
      call rhs(t+h*m%c(i),yt,k(:,i));nrhs=nrhs+1
    end do
    ylow=y;yhigh=y
    do i=1,m%stages
      ylow=ylow+h*m%b_low(i)*k(:,i);yhigh=yhigh+h*m%b_high(i)*k(:,i)
    end do
  end subroutine embedded_pair

  subroutine adaptive_step(rhs,t,y,h,m,rtol,atol,target,direction,accepted,nrhs)
    procedure(ode_rhs)::rhs
    real(dp),intent(inout)::t,h
    real(dp),intent(inout)::y(:)
    type(rk_method),intent(in)::m
    real(dp),intent(in)::rtol,atol,target,direction
    logical,intent(out)::accepted
    integer,intent(inout)::nrhs
    real(dp),allocatable::yl(:),yh(:),scale(:)
    real(dp)::err,fac,hnew
    allocate(yl(size(y)),yh(size(y)),scale(size(y)))
    call embedded_pair(rhs,t,y,h,m,yl,yh,nrhs)
    scale=atol+rtol*max(abs(y),abs(yh))
    err=maxval(abs(yh-yl)/scale)
    if(err<=1.0_dp)then
      t=t+h;y=yh;accepted=.true.
      if(err<=tiny(1.0_dp))then;fac=10.0_dp
      else;fac=min(10.0_dp,max(0.2_dp,0.9_dp*err**(-1.0_dp/real(m%qerr,dp))))
      end if
      hnew=h*fac
    else
      accepted=.false.;fac=max(0.2_dp,0.9_dp*err**(-1.0_dp/real(m%qerr,dp)));hnew=h*fac
    end if
    if(direction*(t+hnew-target)>0.0_dp)hnew=target-t
    if(abs(hnew)<=100.0_dp*tiny(1.0_dp)*max(1.0_dp,abs(t)))then
      hnew=direction*100.0_dp*tiny(1.0_dp)*max(1.0_dp,abs(t))
    end if
    h=hnew
  end subroutine adaptive_step

  subroutine implicit_step(rhs,t,y,h,m,yout,ok)
    procedure(ode_rhs)::rhs
    real(dp),intent(in)::t,y(:),h
    type(rk_method),intent(in)::m
    real(dp),intent(out)::yout(:)
    logical,intent(out)::ok
    integer::n,s,ns,iter,j,info
    real(dp),allocatable::k(:),f(:),jac(:,:),delta(:),kp(:)
    real(dp)::epsj
    n=size(y);s=m%stages;ns=n*s
    allocate(k(ns),f(ns),jac(ns,ns),delta(ns),kp(ns));k=0.0_dp
    call stage_initial_guess(rhs,t,y,h,m,k)
    ok=.false.
    do iter=1,30
      call implicit_residual(rhs,t,y,h,m,k,f)
      if(maxval(abs(f))<1e-11_dp)then;ok=.true.;exit;end if
      do j=1,ns
        kp=k;epsj=sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(k(j)));kp(j)=kp(j)+epsj
        call implicit_residual(rhs,t,y,h,m,kp,delta)
        jac(:,j)=(delta-f)/epsj
      end do
      delta=-f;call dense_solve(jac,delta,info)
      if(info/=0)return
      k=k+delta
      if(maxval(abs(delta))<1e-11_dp*(1.0_dp+maxval(abs(k))))then;ok=.true.;exit;end if
    end do
    yout=y
    do j=1,s;yout=yout+h*m%b_high(j)*k((j-1)*n+1:j*n);end do
  end subroutine implicit_step

  subroutine stage_initial_guess(rhs,t,y,h,m,k)
    procedure(ode_rhs)::rhs
    real(dp),intent(in)::t,y(:),h
    type(rk_method),intent(in)::m
    real(dp),intent(out)::k(:)
    real(dp)::dy(size(y));integer::j,n
    n=size(y);call rhs(t,y,dy)
    do j=1,m%stages;k((j-1)*n+1:j*n)=dy;end do
    if(h < -huge(1.0_dp))stop
  end subroutine stage_initial_guess

  subroutine implicit_residual(rhs,t,y,h,m,k,res)
    procedure(ode_rhs)::rhs
    real(dp),intent(in)::t,y(:),h,k(:)
    type(rk_method),intent(in)::m
    real(dp),intent(out)::res(:)
    integer::i,j,n;real(dp)::yt(size(y)),dy(size(y))
    n=size(y)
    do i=1,m%stages
      yt=y
      do j=1,m%stages;yt=yt+h*m%a(i,j)*k((j-1)*n+1:j*n);end do
      call rhs(t+h*m%c(i),yt,dy)
      res((i-1)*n+1:i*n)=k((i-1)*n+1:i*n)-dy
    end do
  end subroutine implicit_residual

  subroutine dense_solve(a,b,info)
    real(dp),intent(inout)::a(:,:)
    real(dp),intent(inout)::b(:)
    integer,intent(out)::info
    integer::n,i,j,k,p;real(dp)::mx,tmp,f
    n=size(b);info=0
    do k=1,n-1
      p=k;mx=abs(a(k,k))
      do i=k+1,n;if(abs(a(i,k))>mx)then;mx=abs(a(i,k));p=i;end if;end do
      if(mx<=tiny(1.0_dp))then;info=k;return;end if
      if(p/=k)then
        do j=k,n;tmp=a(k,j);a(k,j)=a(p,j);a(p,j)=tmp;end do
        tmp=b(k);b(k)=b(p);b(p)=tmp
      end if
      do i=k+1,n
        f=a(i,k)/a(k,k);a(i,k)=0.0_dp;a(i,k+1:n)=a(i,k+1:n)-f*a(k,k+1:n);b(i)=b(i)-f*b(k)
      end do
    end do
    if(abs(a(n,n))<=tiny(1.0_dp))then;info=n;return;end if
    do i=n,1,-1;b(i)=(b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i);end do
  end subroutine dense_solve

  subroutine allocate_method(m,id,stages,adaptive,implicit,qerr)
    type(rk_method),intent(out)::m
    character(len=*),intent(in)::id
    integer,intent(in)::stages,qerr
    logical,intent(in)::adaptive,implicit
    m%id=id;m%stages=stages;m%adaptive=adaptive;m%implicit=implicit;m%qerr=qerr
    allocate(m%a(stages,stages),m%b_low(stages),m%b_high(stages),m%c(stages))
    m%a=0.0_dp;m%b_low=0.0_dp;m%b_high=0.0_dp;m%c=0.0_dp
  end subroutine allocate_method

  subroutine set_rowmajor(a,nrow,ncol,vals)
    real(dp),intent(inout)::a(:,:)
    integer,intent(in)::nrow,ncol
    real(dp),intent(in)::vals(:)
    integer::i,j,k
    k=0
    do i=1,nrow
      do j=1,ncol;k=k+1;a(i,j)=vals(k);end do
    end do
  end subroutine set_rowmajor

  pure logical function valid_times(t) result(ok)
    real(dp),intent(in)::t(:);integer::i
    ok=size(t)>=1;if(size(t)<2)return
    if(t(2)>t(1))then
      do i=2,size(t);if(t(i)<=t(i-1))then;ok=.false.;return;end if;end do
    else if(t(2)<t(1))then
      do i=2,size(t);if(t(i)>=t(i-1))then;ok=.false.;return;end if;end do
    else;ok=.false.;end if
  end function valid_times

  pure function lower(s) result(out)
    character(len=*),intent(in)::s;character(len=len(s))::out;integer::i,k
    out=s
    do i=1,len(s);k=iachar(out(i:i));if(k>=iachar('A').and.k<=iachar('Z'))out(i:i)=achar(k+32);end do
  end function lower

  pure function solver_stats_zero() result(x)
    use desolve_types, only : solver_stats
    type(solver_stats)::x
    x=solver_stats()
  end function solver_stats_zero
end module desolve_rk
