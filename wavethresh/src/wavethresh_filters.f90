! SPDX-License-Identifier: GPL-2.0-or-later
! Filter tables translated from wavethresh 4.7.3 filter.select().
module wavethresh_filters
   use wavethresh_types, only : dp, wt_filter_t, multiple_filter_t
   implicit none
   private
   public :: filter_select, compare_filters, qmf, mfilter_select
contains

   pure function mfilter_select(filter_type) result(filter)
      character(len=*), intent(in), optional :: filter_type !! Multiple-wavelet family, Geronimo or Donovan3; default is Geronimo.
      type(multiple_filter_t) :: filter
      character(len=24) :: kind
      real(dp) :: root5

      kind = "Geronimo"
      if (present(filter_type)) kind = filter_type
      filter%filter_type = trim(kind)
      select case (trim(kind))
      case ("Geronimo")
         filter%name = "Geronimo Multiwavelets"
         filter%nphi = 2
         filter%npsi = 2
         filter%nh = 4
         filter%ndecim = 2
         allocate(filter%h(16), source=0.0_dp)
         allocate(filter%g(16), source=0.0_dp)
         filter%h(1) = 0.42426406871193001_dp
         filter%h(2) = 0.80000000000000004_dp
         filter%h(3) = -0.050000000000000003_dp
         filter%h(4) = -0.21213203435596001_dp
         filter%h(5) = 0.42426406871193001_dp
         filter%h(7) = 0.45000000000000001_dp
         filter%h(8) = 0.70710678118655002_dp
         filter%h(11) = 0.45000000000000001_dp
         filter%h(12) = -0.21213203435596001_dp
         filter%h(15) = -0.050000000000000003_dp
         filter%g(1) = -0.050000000000000003_dp
         filter%g(2) = -0.21213203435596401_dp
         filter%g(3) = 0.070710678118654793_dp
         filter%g(4) = 0.29999999999999999_dp
         filter%g(5) = 0.45000000000000001_dp
         filter%g(6) = -0.70710678118654802_dp
         filter%g(7) = -0.63639610306789296_dp
         filter%g(9) = 0.45000000000000001_dp
         filter%g(10) = -0.21213203435596401_dp
         filter%g(11) = 0.63639610306789296_dp
         filter%g(12) = -0.29999999999999999_dp
         filter%g(13) = -0.050000000000000003_dp
         filter%g(15) = -0.070710678118654793_dp
      case ("Donovan3")
         filter%name = "Donovan Multiwavelets, 3 functions"
         filter%nphi = 3
         filter%npsi = 3
         filter%nh = 4
         filter%ndecim = 2
         allocate(filter%h(36), source=0.0_dp)
         allocate(filter%g(36), source=0.0_dp)
         root5 = sqrt(5.0_dp)
         filter%h(2) = -sqrt(154.0_dp) * (3.0_dp + 2.0_dp * root5) / 3696.0_dp
         filter%h(3) = sqrt(14.0_dp) * (2.0_dp + 5.0_dp * root5) / 1232.0_dp
         filter%h(10) = -sqrt(2.0_dp) * (3.0_dp + 2.0_dp * root5) / 44.0_dp
         filter%h(11) = sqrt(154.0_dp) * (67.0_dp + 30.0_dp * root5) / 3696.0_dp
         filter%h(12) = sqrt(14.0_dp) * (-10.0_dp + root5) / 112.0_dp
         filter%h(19) = 1.0_dp / sqrt(2.0_dp)
         filter%h(20) = sqrt(154.0_dp) * (67.0_dp - 30.0_dp * root5) / 3696.0_dp
         filter%h(21) = sqrt(14.0_dp) * (10.0_dp + root5) / 112.0_dp
         filter%h(23) = 3.0_dp * sqrt(2.0_dp) / 8.0_dp
         filter%h(24) = sqrt(22.0_dp) * (-4.0_dp + root5) / 88.0_dp
         filter%h(26) = sqrt(22.0_dp) * (32.0_dp + 7.0_dp * root5) / 264.0_dp
         filter%h(27) = sqrt(2.0_dp) * (-5.0_dp + 4.0_dp * root5) / 88.0_dp
         filter%h(28) = sqrt(2.0_dp) * (-3.0_dp + 2.0_dp * root5) / 44.0_dp
         filter%h(29) = sqrt(154.0_dp) * (-3.0_dp + 2.0_dp * root5) / 3696.0_dp
         filter%h(30) = sqrt(14.0_dp) * (-2.0_dp + 5.0_dp * root5) / 1232.0_dp
         filter%h(31) = sqrt(154.0_dp) / 22.0_dp
         filter%h(32) = 3.0_dp * sqrt(2.0_dp) / 8.0_dp
         filter%h(33) = sqrt(22.0_dp) * (4.0_dp + root5) / 88.0_dp
         filter%h(34) = -sqrt(70.0_dp) / 22.0_dp
         filter%h(35) = sqrt(22.0_dp) * (-32.0_dp + 7.0_dp * root5) / 264.0_dp
         filter%h(36) = -sqrt(2.0_dp) * (5.0_dp + 4.0_dp * root5) / 88.0_dp
         filter%g(5) = sqrt(154.0_dp) * (3.0_dp + 2.0_dp * root5) / 3696.0_dp
         filter%g(6) = -sqrt(14.0_dp) * (2.0_dp + 5.0_dp * root5) / 1232.0_dp
         filter%g(8) = -sqrt(7.0_dp) * (1.0_dp + root5) / 336.0_dp
         filter%g(9) = sqrt(77.0_dp) * (-1.0_dp + 3.0_dp * root5) / 1232.0_dp
         filter%g(13) = sqrt(2.0_dp) * (3.0_dp + 2.0_dp * root5) / 44.0_dp
         filter%g(14) = -sqrt(154.0_dp) * (67.0_dp + 30.0_dp * root5) / 3696.0_dp
         filter%g(15) = sqrt(14.0_dp) * (10.0_dp - root5) / 112.0_dp
         filter%g(16) = -sqrt(11.0_dp) * (1.0_dp + root5) / 44.0_dp
         filter%g(17) = sqrt(7.0_dp) * (29.0_dp + 13.0_dp * root5) / 336.0_dp
         filter%g(18) = sqrt(77.0_dp) * (-75.0_dp + 17.0_dp * root5) / 1232.0_dp
         filter%g(20) = sqrt(77.0_dp) * (-2.0_dp + root5) / 264.0_dp
         filter%g(21) = sqrt(7.0_dp) * (13.0_dp - 6.0_dp * root5) / 88.0_dp
         filter%g(22) = 1.0_dp / sqrt(2.0_dp)
         filter%g(23) = sqrt(154.0_dp) * (-67.0_dp + 30.0_dp * root5) / 3696.0_dp
         filter%g(24) = -sqrt(14.0_dp) * (10.0_dp + root5) / 112.0_dp
         filter%g(26) = sqrt(7.0_dp) * (-29.0_dp + 13.0_dp * root5) / 336.0_dp
         filter%g(27) = -sqrt(77.0_dp) * (75.0_dp + 17.0_dp * root5) / 1232.0_dp
         filter%g(28) = 13.0_dp / 22.0_dp
         filter%g(29) = -sqrt(77.0_dp) * (2.0_dp + root5) / 264.0_dp
         filter%g(30) = -sqrt(7.0_dp) * (13.0_dp + 6.0_dp * root5) / 88.0_dp
         filter%g(31) = sqrt(2.0_dp) * (3.0_dp - 2.0_dp * root5) / 44.0_dp
         filter%g(32) = sqrt(154.0_dp) * (3.0_dp - 2.0_dp * root5) / 3696.0_dp
         filter%g(33) = sqrt(14.0_dp) * (2.0_dp - 5.0_dp * root5) / 1232.0_dp
         filter%g(34) = sqrt(11.0_dp) * (1.0_dp - root5) / 44.0_dp
         filter%g(35) = sqrt(7.0_dp) * (1.0_dp - root5) / 336.0_dp
         filter%g(36) = -sqrt(77.0_dp) * (1.0_dp + 3.0_dp * root5) / 1232.0_dp
      case default
         filter%message = "bad multiple-wavelet filter specified"
         return
      end select
      filter%ok = .true.
      filter%message = "ok"
   end function mfilter_select

   pure function qmf(low) result(high)
      real(dp), intent(in) :: low(:) !! Normalized low-pass quadrature-mirror filter coefficients.
      real(dp), allocatable :: high(:)
      integer :: i
      integer :: n
      n = size(low)
      allocate(high(n))
      do i = 1, n
         high(i) = (-1.0_dp)**real(i - 1, dp) * low(n - i + 1)
      end do
   end function qmf

   pure function compare_filters(first, second) result(equal)
      type(wt_filter_t), intent(in) :: first !! First wavelet filter descriptor to compare.
      type(wt_filter_t), intent(in) :: second !! Second wavelet filter descriptor to compare.
      logical :: equal
      equal = trim(first%family) == trim(second%family) .and. &
         abs(first%filter_number - second%filter_number) <= 16.0_dp * epsilon(1.0_dp)
   end function compare_filters

   pure function filter_select(filter_number, family, constant) result(filter)
      real(dp), intent(in) :: filter_number !! Wavethresh filter number; integer for real families.
      character(len=*), intent(in) :: family !! Wavethresh family name, such as DaubLeAsymm or Coiflets.
      real(dp), intent(in), optional :: constant !! Optional divisor applied to all selected coefficients.
      type(wt_filter_t) :: filter
      real(dp) :: scale
      integer :: n
      integer :: i
      scale = 1.0_dp
      if (present(constant)) scale = constant
      if (abs(scale) <= tiny(1.0_dp)) then
         filter%message = "constant must be nonzero"
         return
      end if
      filter%family = family
      filter%filter_number = filter_number
      select case (trim(family))
      case ("DaubExPhase")
         n = nint(filter_number)
         select case (n)
         case (1)
            filter%low = [ &

               0.707106781186548_dp, 0.707106781186548_dp ]
            filter%name = "Daubechies extremal phase 1"
         case (2)
            filter%low = [ &

               0.482962913145_dp, 0.836516303738_dp, 0.224143868042_dp, -0.129409522551_dp ]
            filter%name = "Daubechies extremal phase 2"
         case (3)
            filter%low = [ &

               0.33267055295_dp, 0.806891509311_dp, 0.459877502118_dp, -0.13501102001_dp, -0.085441273882_dp, &

               0.035226291882_dp ]
            filter%name = "Daubechies extremal phase 3"
         case (4)
            filter%low = [ &

               0.230377813309_dp, 0.714846570553_dp, 0.63088076793_dp, -0.027983769417_dp, -0.187034811719_dp, &

               0.030841381836_dp, 0.032883011667_dp, -0.010597401785_dp ]
            filter%name = "Daubechies extremal phase 4"
         case (5)
            filter%low = [ &

               0.160102397974_dp, 0.603829269797_dp, 0.724308528438_dp, 0.138428145901_dp, -0.242294887066_dp, &

               -0.032244869585_dp, 0.07757149384_dp, -0.006241490213_dp, -0.012580752_dp, 0.003335725285_dp ]
            filter%name = "Daubechies extremal phase 5"
         case (6)
            filter%low = [ &

               0.11154074335_dp, 0.494623890398_dp, 0.751133908021_dp, 0.315250351709_dp, -0.226264693965_dp, &

               -0.129766867567_dp, 0.097501605587_dp, 0.02752286553_dp, -0.031582039318_dp, 0.000553842201_dp, &

               0.004777257511_dp, -0.001077301085_dp ]
            filter%name = "Daubechies extremal phase 6"
         case (7)
            filter%low = [ &

               0.077852054085_dp, 0.396539319482_dp, 0.729132090846_dp, 0.469782287405_dp, -0.143906003929_dp, &

               -0.224036184994_dp, 0.071309219267_dp, 0.080612609151_dp, -0.038029936935_dp, &

               -0.016574541631_dp, 0.012550998556_dp, 0.000429577973_dp, -0.001801640704_dp, 0.0003537138_dp ]
            filter%name = "Daubechies extremal phase 7"
         case (8)
            filter%low = [ &

               0.054415842243_dp, 0.312871590914_dp, 0.675630736297_dp, 0.585354683654_dp, -0.015829105256_dp, &

               -0.284015542962_dp, 0.000472484574_dp, 0.12874742662_dp, -0.017369301002_dp, &

               -0.044088253931_dp, 0.013981027917_dp, 0.008746094047_dp, -0.004870352993_dp, &

               -0.000391740373_dp, 0.000675449406_dp, -0.000117476784_dp ]
            filter%name = "Daubechies extremal phase 8"
         case (9)
            filter%low = [ &

               0.038077947364_dp, 0.243834674613_dp, 0.60482312369_dp, 0.657288078051_dp, 0.133197385825_dp, &

               -0.293273783279_dp, -0.096840783223_dp, 0.148540749338_dp, 0.030725681479_dp, &

               -0.067632829061_dp, 0.000250947115_dp, 0.022361662124_dp, -0.004723204758_dp, &

               -0.004281503682_dp, 0.001847646883_dp, 0.000230385764_dp, -0.000251963189_dp, 3.934732e-05_dp ]
            filter%name = "Daubechies extremal phase 9"
         case (10)
            filter%low = [ &

               0.026670057901_dp, 0.188176800078_dp, 0.527201188932_dp, 0.688459039454_dp, 0.281172343661_dp, &

               -0.249846424327_dp, -0.195946274377_dp, 0.127369340336_dp, 0.093057364604_dp, &

               -0.071394147166_dp, -0.029457536822_dp, 0.033212674059_dp, 0.003606553567_dp, &

               -0.010733175483_dp, 0.001395351747_dp, 0.001992405295_dp, -0.000685856695_dp, &

               -0.000116466855_dp, 9.358867e-05_dp, -1.3264203e-05_dp ]
            filter%name = "Daubechies extremal phase 10"
         case default
            filter%message = "DaubExPhase filter number must be 1 through 10"
            return
         end select
      case ("DaubLeAsymm")
         n = nint(filter_number)
         select case (n)
         case (4)
            filter%low = [ &

               -0.0757657147893567_dp, -0.0296355276459604_dp, 0.497618667632563_dp, 0.803738751805386_dp, &

               0.297857795605605_dp, -0.0992195435769564_dp, -0.0126039672622638_dp, 0.0322231006040782_dp ]
            filter%name = "Daubechies least asymmetric 4"
         case (5)
            filter%low = [ &

               0.0273330683451628_dp, 0.0295194909260716_dp, -0.0391342493025811_dp, 0.199397533976983_dp, &

               0.723407690403764_dp, 0.633978963456911_dp, 0.0166021057644233_dp, -0.175328089908097_dp, &

               -0.0211018340249286_dp, 0.0195388827353857_dp ]
            filter%name = "Daubechies least asymmetric 5"
         case (6)
            filter%low = [ &

               0.0154041093273385_dp, 0.00349071208433061_dp, -0.117990111148417_dp, -0.0483117425860007_dp, &

               0.491055941927666_dp, 0.787641141028836_dp, 0.337929421728258_dp, -0.0726375227866039_dp, &

               -0.0210602925126965_dp, 0.0447249017707506_dp, 0.00176771186439837_dp, -0.00780070832476545_dp ]
            filter%name = "Daubechies least asymmetric 6"
         case (7)
            filter%low = [ &

               0.00268181456811643_dp, -0.00104738488896575_dp, -0.0126363034031528_dp, 0.0305155131659067_dp, &

               0.0678926935015983_dp, -0.0495528349370419_dp, 0.0174412550871102_dp, 0.536101917090782_dp, &

               0.767764317004585_dp, 0.288629631750988_dp, -0.140047240442706_dp, -0.107808237703619_dp, &

               0.00401024487170333_dp, 0.010268176708497_dp ]
            filter%name = "Daubechies least asymmetric 7"
         case (8)
            filter%low = [ &

               0.0018899503329009_dp, -0.000302920514551664_dp, -0.0149522583367938_dp, &

               0.00380875201406036_dp, 0.0491371796734807_dp, -0.0272190299168159_dp, -0.0519458381078792_dp, &

               0.364441894835986_dp, 0.77718575169981_dp, 0.48135965125924_dp, -0.0612733590679137_dp, &

               -0.143294238351066_dp, 0.00760748732528538_dp, 0.0316950878103478_dp, -0.00054213233163559_dp, &

               -0.00338241595135971_dp ]
            filter%name = "Daubechies least asymmetric 8"
         case (9)
            filter%low = [ &

               0.00106949003265249_dp, -0.000473154498587299_dp, -0.0102640640276723_dp, &

               0.00885926749350085_dp, 0.0620777893026874_dp, -0.0182337707798032_dp, -0.191550831296252_dp, &

               0.0352724880358911_dp, 0.617338449140593_dp, 0.717897082763343_dp, 0.238760914607125_dp, &

               -0.0545689584305094_dp, 0.000583462746330467_dp, 0.0302248788579523_dp, -0.0115282102079706_dp, &

               -0.0132719677815169_dp, 0.000619780889054126_dp, 0.00140091552556991_dp ]
            filter%name = "Daubechies least asymmetric 9"
         case (10)
            filter%low = [ &

               0.000770159808941683_dp, 9.56326707637102e-05_dp, -0.0086412992741304_dp, &

               -0.00146538258303966_dp, 0.0459272392141469_dp, 0.0116098939105411_dp, -0.15949427882413_dp, &

               -0.0708805357960178_dp, 0.471690666743878_dp, 0.769510036853189_dp, 0.383826761145002_dp, &

               -0.0355367402980268_dp, -0.0319900568214638_dp, 0.0499949720686861_dp, 0.0057649120443445_dp, &

               -0.0203549397996833_dp, -0.000804358934368569_dp, 0.00459317358270837_dp, &

               5.70360843270715e-05_dp, -0.000459329420451864_dp ]
            filter%name = "Daubechies least asymmetric 10"
         case default
            filter%message = "DaubLeAsymm filter number must be 4 through 10"
            return
         end select
      case ("Coiflets")
         n = nint(filter_number)
         select case (n)
         case (1)
            filter%low = [ &

               -0.0727326195128539_dp, 0.337897662457809_dp, 0.852572020212255_dp, 0.384864846864203_dp, &

               -0.0727329651127074_dp, -0.0156557281354645_dp ]
            filter%name = "Coiflet 1"
         case (2)
            filter%low = [ &

               0.0163873410753545_dp, -0.0414649396386779_dp, -0.0673725542838937_dp, 0.386110001012665_dp, &

               0.81272364413712_dp, 0.41700519333898_dp, -0.0764886031912219_dp, -0.0594344179948016_dp, &

               0.0236801717159357_dp, 0.00561143536672321_dp, -0.00182320836725208_dp, &

               -0.000720549446782329_dp ]
            filter%name = "Coiflet 2"
         case (3)
            filter%low = [ &

               -0.00379351332976728_dp, 0.00778259683886156_dp, 0.023452695464428_dp, -0.0657719049475929_dp, &

               -0.0611233849680726_dp, 0.405176852524648_dp, 0.793777283620651_dp, 0.428483516296625_dp, &

               -0.0717998205515808_dp, -0.0823019260292552_dp, 0.0345550214622448_dp, 0.0158805435031425_dp, &

               -0.00900797612110523_dp, -0.00257451780754416_dp, 0.00111751876947639_dp, &

               0.000466216996882439_dp, -7.09832960670734e-05_dp, -3.45997671793583e-05_dp ]
            filter%name = "Coiflet 3"
         case (4)
            filter%low = [ &

               0.00089231360352849_dp, -0.00162949222190534_dp, -0.00734616629087254_dp, &

               0.0160689450339218_dp, 0.0266823066925482_dp, -0.0812666934173269_dp, -0.0560773154140335_dp, &

               0.415308419906845_dp, 0.782238998007962_dp, 0.434386019441902_dp, -0.0666274758685574_dp, &

               -0.0962204462045038_dp, 0.0393344269069628_dp, 0.0250822654680792_dp, -0.015211733625225_dp, &

               -0.00565828684783107_dp, 0.00375143623364599_dp, 0.00126656188828675_dp, &

               -0.000589020797256532_dp, -0.000259974596328334_dp, 6.23390288041529e-05_dp, &

               3.12298801079863e-05_dp, -3.25968023688337e-06_dp, -1.78498500308826e-06_dp ]
            filter%name = "Coiflet 4"
         case (5)
            filter%low = [ &

               -0.000212080839825006_dp, 0.000358589687931597_dp, 0.00217823678259091_dp, &

               -0.0041593590646447_dp, -0.0101311176622843_dp, 0.0234081618793488_dp, 0.0281680228926348_dp, &

               -0.0919200066090903_dp, -0.0520431580902793_dp, 0.421566173498989_dp, 0.774289562152506_dp, &

               0.437991556919172_dp, -0.0620359628377363_dp, -0.105574210269531_dp, 0.0412892094614174_dp, &

               0.0326835785130209_dp, -0.0197617819143973_dp, -0.00916423116339827_dp, 0.00676418488304408_dp, &

               0.00243337363716508_dp, -0.00166286341933666_dp, -0.000638131343108925_dp, &

               0.000302259581843289_dp, 0.000140541149716088_dp, -4.13404322766466e-05_dp, &

               -2.13150268120873e-05_dp, 3.73465517551487e-06_dp, 2.06376185157106e-06_dp, &

               -1.67442885784974e-07_dp, -9.51765727477093e-08_dp ]
            filter%name = "Coiflet 5"
         case default
            filter%message = "Coiflets filter number must be 1 through 5"
            return
         end select
      case ("Yates")
         if (nint(filter_number) /= 1) then
            filter%message = "Yates only provides filter number 1"
            return
         end if
         filter%low = [ -1.0_dp / sqrt(2.0_dp), 1.0_dp / sqrt(2.0_dp) ]
         filter%name = "Yates"
      case ("LittlewoodPaley")
         n = nint(filter_number)
         if (n < 0) then
            filter%message = "LittlewoodPaley filter number must be nonnegative"
            return
         end if
         allocate(filter%low(2 * n + 1))
         filter%low(n + 1) = 1.0_dp / sqrt(2.0_dp)
         do i = 1, n
            filter%low(n + 1 - i) = sin(0.5_dp * acos(-1.0_dp) * real(i, dp)) / &
               (0.5_dp * acos(-1.0_dp) * real(i, dp) * sqrt(2.0_dp))
            filter%low(n + 1 + i) = filter%low(n + 1 - i)
         end do
         filter%name = "Littlewood-Paley"
      case ("MagKing")
         if (nint(filter_number) /= 4) then
            filter%message = "MagKing only provides filter number 4"
            return
         end if
         filter%is_complex = .true.
         filter%low_complex = [cmplx(0.1_dp, -0.1_dp, dp), cmplx(0.4_dp, -0.1_dp, dp), &
            cmplx(0.4_dp, 0.1_dp, dp), cmplx(0.1_dp, 0.1_dp, dp)]
         filter%high_complex = [cmplx(-1.0_dp, -2.0_dp, dp) / cmplx(14.0_dp, 0.0_dp, dp), &
            cmplx(5.0_dp, 2.0_dp, dp) / cmplx(14.0_dp, 0.0_dp, dp), cmplx(-5.0_dp, 2.0_dp, dp) / cmplx(14.0_dp, 0.0_dp, dp), &
            cmplx(1.0_dp, -2.0_dp, dp) / cmplx(14.0_dp, 0.0_dp, dp)]
         filter%name = "Magarey-Kingsbury 4-tap"
      case ("Nason")
         if (nint(filter_number) /= 3) then
            filter%message = "Nason only provides filter number 3"
            return
         end if
         filter%is_complex = .true.
         filter%low_complex = [cmplx(-0.066291_dp, 0.085581_dp, dp), &
            cmplx(0.110485_dp, 0.085558_dp, dp), cmplx(0.662912_dp, -0.171163_dp, dp), &
            cmplx(0.662912_dp, -0.171163_dp, dp), cmplx(0.110485_dp, 0.085558_dp, dp), &
            cmplx(-0.066291_dp, 0.085581_dp, dp)]
         filter%high_complex = [cmplx(-0.066291_dp, 0.085581_dp, dp), &
            cmplx(-0.110485_dp, -0.085558_dp, dp), cmplx(0.662912_dp, -0.171163_dp, dp), &
            cmplx(-0.662912_dp, 0.171163_dp, dp), cmplx(0.110485_dp, 0.085558_dp, dp), &
            cmplx(0.066291_dp, -0.085581_dp, dp)]
         filter%name = "Nason complex 6-tap"
      case ("Lawton")
         if (nint(filter_number) /= 3) then
            filter%message = "Lawton only provides filter number 3"
            return
         end if
         filter%is_complex = .true.
         filter%low_complex = [cmplx(-0.066291_dp, 0.085581_dp, dp), &
            cmplx(0.110485_dp, 0.085558_dp, dp), cmplx(0.662912_dp, -0.171163_dp, dp), &
            cmplx(0.662912_dp, -0.171163_dp, dp), cmplx(0.110485_dp, 0.0_dp, dp), &
            cmplx(-0.066291_dp, 0.085581_dp, dp)]
         filter%high_complex = [cmplx(-0.066291_dp, -0.085581_dp, dp), &
            cmplx(-0.110485_dp, 0.085558_dp, dp), cmplx(0.662912_dp, 0.171163_dp, dp), &
            cmplx(-0.662912_dp, -0.171163_dp, dp), cmplx(0.110485_dp, -0.085558_dp, dp), &
            cmplx(0.066291_dp, 0.085581_dp, dp)]
         filter%name = "Lawton complex 6-tap"
      case ("LinaMayrand")
         filter%message = "LinaMayrand matrix filters are not translated in this release"
         return
      case default
         filter%message = "unknown wavethresh filter family"
         return
      end select
      if (.not. filter%is_complex) then
         filter%low = filter%low / scale
         filter%high = qmf(filter%low)
      else
         filter%low_complex = filter%low_complex / cmplx(scale, 0.0_dp, dp)
         filter%high_complex = filter%high_complex / cmplx(scale, 0.0_dp, dp)
      end if
      filter%ok = .true.
      filter%message = "ok"
   end function filter_select

end module wavethresh_filters
