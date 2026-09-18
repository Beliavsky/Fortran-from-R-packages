! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the R package wavelets 0.3-0.2 by Eric Aldrich.
! Modern Fortran translation and modifications: 2026-09-10.
module wavelets_filters
   use wavelets_kinds, only : dp
   use wavelets_types, only : wt_filter_type
   implicit none
   private

   public :: wt_filter
   public :: wt_filter_named
   public :: wt_filter_coefficients
   public :: wt_filter_qmf
   public :: wt_filter_equivalent
   public :: wt_filter_shift
   public :: waveletshift_dwt
   public :: scalingshift_dwt

   interface wt_filter
      module procedure wt_filter_named
      module procedure wt_filter_coefficients
   end interface wt_filter

contains

   pure subroutine wt_filter_named(name, filter, modwt, level, ierr)
      character(len=*), intent(in) :: name !! Named upstream filter, such as "la8", "d4", or "haar".
      type(wt_filter_type), intent(out) :: filter !! Constructed wavelet/scaling filter and metadata.
      logical, intent(in), optional :: modwt !! If true, normalize coefficients for the MODWT; default is false.
      integer, intent(in), optional :: level !! Equivalent-filter level to construct; must be positive and defaults to 1.
      integer, intent(out), optional :: ierr !! Status: 0 success, 1 unknown filter name, 2 invalid level.

      type(wt_filter_type) :: base_filter
      real(dp), allocatable :: g(:)
      character(len=:), allocatable :: filter_class
      character(len=:), allocatable :: canonical_name
      logical :: use_modwt
      integer :: lev
      integer :: status

      use_modwt = .false.
      if (present(modwt)) use_modwt = modwt
      lev = 1
      if (present(level)) lev = level
      status = 0
      if (lev < 1) then
         status = 2
         call clear_filter(filter)
         if (present(ierr)) ierr = status
         return
      end if

      call named_scaling_coefficients(trim(name), g, filter_class, canonical_name, status)
      if (status /= 0) then
         call clear_filter(filter)
         if (present(ierr)) ierr = status
         return
      end if

      if (use_modwt) g = g / sqrt(2.0_dp)
      base_filter%l = size(g)
      base_filter%level = 1
      base_filter%g = g
      call wt_filter_qmf(g, base_filter%h, inverse=.true.)
      base_filter%wt_class = filter_class
      base_filter%wt_name = canonical_name
      if (use_modwt) then
         base_filter%transform = "modwt"
      else
         base_filter%transform = "dwt"
      end if

      if (lev > 1) then
         call wt_filter_equivalent(base_filter, lev, filter, status)
      else
         filter = base_filter
      end if
      if (present(ierr)) ierr = status
   end subroutine wt_filter_named

   pure subroutine wt_filter_coefficients(coefficients, filter, modwt, level, ierr)
      real(dp), intent(in) :: coefficients(:) !! User-supplied wavelet filter coefficients; length must be positive and even.
      type(wt_filter_type), intent(out) :: filter !! Constructed custom filter with quadrature-mirror scaling coefficients.
      logical, intent(in), optional :: modwt !! If true, mark the custom filter for MODWT use; default is false.
      integer, intent(in), optional :: level !! Equivalent-filter level to construct; must be positive and defaults to 1.
      integer, intent(out), optional :: ierr !! Status: 0 success, 1 invalid coefficient length, 2 invalid level.

      type(wt_filter_type) :: base_filter
      logical :: use_modwt
      integer :: lev
      integer :: status

      use_modwt = .false.
      if (present(modwt)) use_modwt = modwt
      lev = 1
      if (present(level)) lev = level
      status = 0

      if (size(coefficients) < 2 .or. mod(size(coefficients), 2) /= 0) then
         status = 1
         call clear_filter(filter)
         if (present(ierr)) ierr = status
         return
      end if
      if (lev < 1) then
         status = 2
         call clear_filter(filter)
         if (present(ierr)) ierr = status
         return
      end if

      base_filter%l = size(coefficients)
      base_filter%level = 1
      base_filter%h = coefficients
      call wt_filter_qmf(coefficients, base_filter%g)
      base_filter%wt_class = "none"
      base_filter%wt_name = "none"
      if (use_modwt) then
         base_filter%transform = "modwt"
      else
         base_filter%transform = "dwt"
      end if

      if (lev > 1) then
         call wt_filter_equivalent(base_filter, lev, filter, status)
      else
         filter = base_filter
      end if
      if (present(ierr)) ierr = status
   end subroutine wt_filter_coefficients

   pure subroutine wt_filter_qmf(x, y, inverse)
      real(dp), intent(in) :: x(:) !! Input filter coefficients in upstream wavelets ordering.
      real(dp), allocatable, intent(out) :: y(:) !! Quadrature-mirror coefficients with reversed order and alternating signs.
      logical, intent(in), optional :: inverse !! If true, use the inverse sign convention used to form wavelet coefficients.

      logical :: inv
      integer :: i
      integer :: l

      inv = .false.
      if (present(inverse)) inv = inverse
      l = size(x)
      allocate(y(l))
      do i = 1, l
         if (inv) then
            y(i) = x(l - i + 1) * alternating_sign(i - 1)
         else
            y(i) = x(l - i + 1) * alternating_sign(i)
         end if
      end do
   end subroutine wt_filter_qmf

   pure subroutine wt_filter_equivalent(filter, j, equivalent, ierr)
      type(wt_filter_type), intent(in) :: filter !! Base level-one wavelet filter.
      integer, intent(in) :: j !! Positive equivalent-filter level.
      type(wt_filter_type), intent(out) :: equivalent !! Equivalent level-j filter using the upstream cascade construction.
      integer, intent(out), optional :: ierr !! Status: 0 success, 1 invalid level or empty filter.

      real(dp), allocatable :: h_last(:)
      real(dp), allocatable :: g_last(:)
      real(dp), allocatable :: hj(:)
      real(dp), allocatable :: gj(:)
      real(dp) :: g_mult
      integer :: base_l
      integer :: l_last
      integer :: l_new
      integer :: lev
      integer :: ell
      integer :: k
      integer :: u
      integer :: status

      status = 0
      if (j < 1 .or. filter%l < 1 .or. .not. allocated(filter%h) .or. .not. allocated(filter%g)) then
         status = 1
         call clear_filter(equivalent)
         if (present(ierr)) ierr = status
         return
      end if

      equivalent = filter
      equivalent%level = j
      if (j == 1) then
         if (present(ierr)) ierr = status
         return
      end if

      base_l = filter%l
      l_last = base_l
      h_last = filter%h
      g_last = filter%g

      do lev = 2, j
         l_new = (2**lev - 1) * (base_l - 1) + 1
         allocate(hj(l_new), gj(l_new))
         do ell = 0, l_new - 1
            u = ell
            if (u >= base_l) then
               g_mult = 0.0_dp
            else
               g_mult = filter%g(u + 1)
            end if
            hj(ell + 1) = g_mult * h_last(1)
            gj(ell + 1) = g_mult * g_last(1)
            do k = 1, l_last - 1
               u = u - 2
               if (u < 0 .or. u >= base_l) then
                  g_mult = 0.0_dp
               else
                  g_mult = filter%g(u + 1)
               end if
               hj(ell + 1) = hj(ell + 1) + g_mult * h_last(k + 1)
               gj(ell + 1) = gj(ell + 1) + g_mult * g_last(k + 1)
            end do
         end do
         call move_alloc(hj, h_last)
         call move_alloc(gj, g_last)
         l_last = l_new
      end do

      equivalent%l = l_last
      equivalent%h = h_last
      equivalent%g = g_last
      if (present(ierr)) ierr = status
   end subroutine wt_filter_equivalent

   pure subroutine wt_filter_shift(filter, levels, shifts, wavelet, coe, modwt, ierr)
      type(wt_filter_type), intent(in) :: filter !! Filter whose phase shift is requested.
      integer, intent(in) :: levels(:) !! Positive decomposition levels at which to compute shifts.
      real(dp), allocatable, intent(out) :: shifts(:) !! Shift for each requested level, in coefficient samples.
      logical, intent(in), optional :: wavelet !! If true compute wavelet-filter shifts; false selects scaling shifts.
      logical, intent(in), optional :: coe !! If true use center-of-energy shifts instead of named-filter phase rules.
      logical, intent(in), optional :: modwt !! Haar special-case selector matching the R argument; default is false.
      integer, intent(out), optional :: ierr !! Status: 0 success, 1 invalid level, 2 invalid filter coefficients.

      logical :: use_wavelet
      logical :: use_coe
      logical :: use_modwt
      real(dp) :: nu
      real(dp) :: delta
      real(dp) :: raw_shift
      real(dp) :: denom
      integer :: i
      integer :: j
      integer :: status

      status = 0
      allocate(shifts(size(levels)))
      shifts = 0.0_dp
      if (any(levels <= 0)) then
         status = 1
         if (present(ierr)) ierr = status
         return
      end if
      if (filter%l < 1 .or. .not. allocated(filter%g) .or. .not. allocated(filter%wt_class) .or. &
          .not. allocated(filter%wt_name) .or. .not. allocated(filter%transform)) then
         status = 2
         if (present(ierr)) ierr = status
         return
      end if

      use_wavelet = .true.
      if (present(wavelet)) use_wavelet = wavelet
      use_coe = .false.
      if (present(coe)) use_coe = coe
      use_modwt = .false.
      if (present(modwt)) use_modwt = modwt

      denom = sum(filter%g**2)
      if (denom <= 0.0_dp) then
         status = 2
         if (present(ierr)) ierr = status
         return
      end if

      if (use_coe) then
         nu = center_of_energy(filter%g)
      else
         select case (filter%wt_class)
         case ("Daubechies")
            if (filter%l == 2) then
               nu = 0.0_dp
            else if (filter%l == 4) then
               nu = 1.0_dp
            else
               nu = center_of_energy(filter%g)
            end if
         case ("Least Asymmetric")
            delta = 0.0_dp
            if (any(filter%l == [8, 12, 16, 20])) delta = 1.0_dp
            if (filter%l == 14) delta = 2.0_dp
            nu = abs(-real(filter%l, dp) / 2.0_dp + delta)
         case ("Coiflet")
            nu = abs(-2.0_dp * real(filter%l, dp) / 3.0_dp + 1.0_dp)
         case ("Best Localized")
            select case (filter%l)
            case (14)
               nu = 5.0_dp
            case (18)
               nu = 11.0_dp
            case (20)
               nu = 9.0_dp
            case default
               nu = center_of_energy(filter%g)
            end select
         case default
            nu = center_of_energy(filter%g)
         end select
      end if

      do i = 1, size(levels)
         j = levels(i)
         if (filter%wt_name == "haar") then
            if (use_modwt) then
               raw_shift = real(2**(j - 1), dp)
            else
               raw_shift = 0.0_dp
            end if
         else if (j > 1) then
            if (use_wavelet) then
               raw_shift = real(2**(j - 1) * (filter%l - 1), dp) - nu
            else
               raw_shift = real(2**j - 1, dp) * nu
            end if
         else
            if (use_wavelet) then
               raw_shift = real(filter%l - 1, dp) - nu
            else
               raw_shift = nu
            end if
         end if

         if (filter%transform == "dwt") then
            shifts(i) = real(ceiling((raw_shift + 1.0_dp) / real(2**j, dp) - 1.0_dp), dp)
         else
            shifts(i) = raw_shift
         end if
      end do
      if (present(ierr)) ierr = status
   end subroutine wt_filter_shift

   pure subroutine waveletshift_dwt(l, j, shift, n, ierr)
      integer, intent(in) :: l !! Even filter length.
      integer, intent(in) :: j !! Positive decomposition level.
      integer, intent(out) :: shift !! Nonnegative wavelet shift, optionally reduced modulo n.
      integer, intent(in), optional :: n !! Optional positive series length used for modulo reduction.
      integer, intent(out), optional :: ierr !! Status: 0 success, 1 invalid filter length or level, 2 invalid n.

      integer :: lj
      integer :: vjh
      integer :: status

      status = 0
      shift = 0
      if (l <= 0 .or. mod(l, 2) /= 0 .or. j <= 0) then
         status = 1
      else
         lj = (2**j - 1) * (l - 1) + 1
         if (l == 10 .or. l == 18) then
            vjh = -lj / 2 + 1
         else if (l == 14) then
            vjh = -lj / 2 - 1
         else
            vjh = -lj / 2
         end if
         shift = abs(vjh)
         if (present(n)) then
            if (n <= 0) then
               status = 2
               shift = 0
            else
               shift = mod(shift, n)
            end if
         end if
      end if
      if (present(ierr)) ierr = status
   end subroutine waveletshift_dwt

   pure subroutine scalingshift_dwt(l, j, shift, n, ierr)
      integer, intent(in) :: l !! Even filter length.
      integer, intent(in) :: j !! Positive decomposition level.
      integer, intent(out) :: shift !! Nonnegative scaling shift, optionally reduced modulo n.
      integer, intent(in), optional :: n !! Optional positive series length used for modulo reduction.
      integer, intent(out), optional :: ierr !! Status: 0 success, 1 invalid filter length or level, 2 invalid n.

      integer :: lj
      integer :: vjg
      integer :: status

      status = 0
      shift = 0
      if (l <= 0 .or. mod(l, 2) /= 0 .or. j <= 0) then
         status = 1
      else
         lj = (2**j - 1) * (l - 1) + 1
         if (l == 10 .or. l == 18) then
            vjg = -((lj - 1) * l) / (2 * (l - 1))
         else if (l == 14) then
            vjg = -((lj - 1) * (l - 4)) / (2 * (l - 1))
         else
            vjg = -((lj - 1) * (l - 2)) / (2 * (l - 1))
         end if
         shift = abs(vjg)
         if (present(n)) then
            if (n <= 0) then
               status = 2
               shift = 0
            else
               shift = mod(shift, n)
            end if
         end if
      end if
      if (present(ierr)) ierr = status
   end subroutine scalingshift_dwt

   pure elemental real(dp) function alternating_sign(power) result(value)
      integer, intent(in) :: power !! Integer exponent in (-1)^power.

      if (mod(power, 2) == 0) then
         value = 1.0_dp
      else
         value = -1.0_dp
      end if
   end function alternating_sign

   pure real(dp) function center_of_energy(g) result(nu)
      real(dp), intent(in) :: g(:) !! Scaling-filter coefficients used to compute the energy centroid.

      integer :: i
      real(dp) :: denom

      denom = sum(g**2)
      nu = 0.0_dp
      do i = 2, size(g)
         nu = nu + real(i - 1, dp) * g(i)**2
      end do
      nu = nu / denom
   end function center_of_energy

   pure subroutine clear_filter(filter)
      type(wt_filter_type), intent(out) :: filter !! Filter object reset to an empty, valid state after an input error.

      filter%l = 0
      filter%level = 1
      filter%wt_class = ""
      filter%wt_name = ""
      filter%transform = ""
   end subroutine clear_filter

   pure subroutine named_scaling_coefficients(name, g, filter_class, canonical_name, ierr)
      character(len=*), intent(in) :: name !! Exact upstream filter name or the d2 alias for Haar.
      real(dp), allocatable, intent(out) :: g(:) !! Level-one scaling-filter coefficients in upstream ordering.
      character(len=:), allocatable, intent(out) :: filter_class !! Upstream filter family name.
      character(len=:), allocatable, intent(out) :: canonical_name !! Canonical upstream filter name.
      integer, intent(out) :: ierr !! Status: 0 success, 1 unknown filter name.

      ierr = 0
      select case (name)
      case ("haar", "d2")
         allocate(g(2))
         g = [ &
            0.7071067811865475_dp, &
            0.7071067811865475_dp ]
         filter_class = "Daubechies"
         canonical_name = "haar"
      case ("d4")
         allocate(g(4))
         g = [ &
            0.4829629131445341_dp, &
            0.8365163037378077_dp, &
            0.2241438680420134_dp, &
            -0.1294095225512603_dp ]
         filter_class = "Daubechies"
         canonical_name = "d4"
      case ("d6")
         allocate(g(6))
         g = [ &
            0.3326705529500827_dp, &
            0.8068915093110928_dp, &
            0.4598775021184915_dp, &
            -0.1350110200102546_dp, &
            -0.0854412738820267_dp, &
            0.0352262918857096_dp ]
         filter_class = "Daubechies"
         canonical_name = "d6"
      case ("d8")
         allocate(g(8))
         g = [ &
            0.2303778133074431_dp, &
            0.7148465705484058_dp, &
            0.6308807679358788_dp, &
            -0.0279837694166834_dp, &
            -0.1870348117179132_dp, &
            0.0308413818353661_dp, &
            0.0328830116666778_dp, &
            -0.0105974017850021_dp ]
         filter_class = "Daubechies"
         canonical_name = "d8"
      case ("d10")
         allocate(g(10))
         g = [ &
            0.1601023979741930_dp, &
            0.6038292697971898_dp, &
            0.7243085284377729_dp, &
            0.1384281459013204_dp, &
            -0.2422948870663824_dp, &
            -0.0322448695846381_dp, &
            0.0775714938400459_dp, &
            -0.0062414902127983_dp, &
            -0.0125807519990820_dp, &
            0.0033357252854738_dp ]
         filter_class = "Daubechies"
         canonical_name = "d10"
      case ("d12")
         allocate(g(12))
         g = [ &
            0.1115407433501094_dp, &
            0.4946238903984530_dp, &
            0.7511339080210954_dp, &
            0.3152503517091980_dp, &
            -0.2262646939654399_dp, &
            -0.1297668675672624_dp, &
            0.0975016055873224_dp, &
            0.0275228655303053_dp, &
            -0.0315820393174862_dp, &
            0.0005538422011614_dp, &
            0.0047772575109455_dp, &
            -0.0010773010853085_dp ]
         filter_class = "Daubechies"
         canonical_name = "d12"
      case ("d14")
         allocate(g(14))
         g = [ &
            0.0778520540850081_dp, &
            0.3965393194819136_dp, &
            0.7291320908462368_dp, &
            0.4697822874052154_dp, &
            -0.1439060039285293_dp, &
            -0.2240361849938538_dp, &
            0.0713092192668312_dp, &
            0.0806126091510820_dp, &
            -0.0380299369350125_dp, &
            -0.0165745416306664_dp, &
            0.0125509985560993_dp, &
            0.0004295779729214_dp, &
            -0.0018016407040474_dp, &
            0.0003537137999745_dp ]
         filter_class = "Daubechies"
         canonical_name = "d14"
      case ("d16")
         allocate(g(16))
         g = [ &
            0.0544158422431049_dp, &
            0.3128715909143031_dp, &
            0.6756307362972904_dp, &
            0.5853546836541907_dp, &
            -0.0158291052563816_dp, &
            -0.2840155429615702_dp, &
            0.0004724845739124_dp, &
            0.1287474266204837_dp, &
            -0.0173693010018083_dp, &
            -0.0440882539307952_dp, &
            0.0139810279173995_dp, &
            0.0087460940474061_dp, &
            -0.0048703529934518_dp, &
            -0.0003917403733770_dp, &
            0.0006754494064506_dp, &
            -0.0001174767841248_dp ]
         filter_class = "Daubechies"
         canonical_name = "d16"
      case ("d18")
         allocate(g(18))
         g = [ &
            0.0380779473638791_dp, &
            0.2438346746125939_dp, &
            0.6048231236901156_dp, &
            0.6572880780512955_dp, &
            0.1331973858249927_dp, &
            -0.2932737832791761_dp, &
            -0.0968407832229524_dp, &
            0.1485407493381306_dp, &
            0.0307256814793395_dp, &
            -0.0676328290613302_dp, &
            0.0002509471148340_dp, &
            0.0223616621236805_dp, &
            -0.0047232047577520_dp, &
            -0.0042815036824636_dp, &
            0.0018476468830564_dp, &
            0.0002303857635232_dp, &
            -0.0002519631889427_dp, &
            0.0000393473203163_dp ]
         filter_class = "Daubechies"
         canonical_name = "d18"
      case ("d20")
         allocate(g(20))
         g = [ &
            0.0266700579005546_dp, &
            0.1881768000776863_dp, &
            0.5272011889317202_dp, &
            0.6884590394536250_dp, &
            0.2811723436606485_dp, &
            -0.2498464243272283_dp, &
            -0.1959462743773399_dp, &
            0.1273693403357890_dp, &
            0.0930573646035802_dp, &
            -0.0713941471663697_dp, &
            -0.0294575368218480_dp, &
            0.0332126740593703_dp, &
            0.0036065535669880_dp, &
            -0.0107331754833036_dp, &
            0.0013953517470692_dp, &
            0.0019924052951930_dp, &
            -0.0006858566949566_dp, &
            -0.0001164668551285_dp, &
            0.0000935886703202_dp, &
            -0.0000132642028945_dp ]
         filter_class = "Daubechies"
         canonical_name = "d20"
      case ("la8")
         allocate(g(8))
         g = [ &
            -0.0757657147893407_dp, &
            -0.0296355276459541_dp, &
            0.4976186676324578_dp, &
            0.8037387518052163_dp, &
            0.2978577956055422_dp, &
            -0.0992195435769354_dp, &
            -0.0126039672622612_dp, &
            0.0322231006040713_dp ]
         filter_class = "Least Asymmetric"
         canonical_name = "la8"
      case ("la10")
         allocate(g(10))
         g = [ &
            0.0195388827353869_dp, &
            -0.0211018340249298_dp, &
            -0.1753280899081075_dp, &
            0.0166021057644243_dp, &
            0.6339789634569490_dp, &
            0.7234076904038076_dp, &
            0.1993975339769955_dp, &
            -0.0391342493025834_dp, &
            0.0295194909260734_dp, &
            0.0273330683451645_dp ]
         filter_class = "Least Asymmetric"
         canonical_name = "la10"
      case ("la12")
         allocate(g(12))
         g = [ &
            0.0154041093273377_dp, &
            0.0034907120843304_dp, &
            -0.1179901111484105_dp, &
            -0.0483117425859981_dp, &
            0.4910559419276396_dp, &
            0.7876411410287941_dp, &
            0.3379294217282401_dp, &
            -0.0726375227866000_dp, &
            -0.0210602925126954_dp, &
            0.0447249017707482_dp, &
            0.0017677118643983_dp, &
            -0.0078007083247650_dp ]
         filter_class = "Least Asymmetric"
         canonical_name = "la12"
      case ("la14")
         allocate(g(14))
         g = [ &
            0.0102681767084968_dp, &
            0.0040102448717033_dp, &
            -0.1078082377036168_dp, &
            -0.1400472404427030_dp, &
            0.2886296317509833_dp, &
            0.7677643170045710_dp, &
            0.5361019170907720_dp, &
            0.0174412550871099_dp, &
            -0.0495528349370410_dp, &
            0.0678926935015971_dp, &
            0.0305155131659062_dp, &
            -0.0126363034031526_dp, &
            -0.0010473848889657_dp, &
            0.0026818145681164_dp ]
         filter_class = "Least Asymmetric"
         canonical_name = "la14"
      case ("la16")
         allocate(g(16))
         g = [ &
            -0.0033824159513594_dp, &
            -0.0005421323316355_dp, &
            0.0316950878103452_dp, &
            0.0076074873252848_dp, &
            -0.1432942383510542_dp, &
            -0.0612733590679088_dp, &
            0.4813596512592012_dp, &
            0.7771857516997478_dp, &
            0.3644418948359564_dp, &
            -0.0519458381078751_dp, &
            -0.0272190299168137_dp, &
            0.0491371796734768_dp, &
            0.0038087520140601_dp, &
            -0.0149522583367926_dp, &
            -0.0003029205145516_dp, &
            0.0018899503329007_dp ]
         filter_class = "Least Asymmetric"
         canonical_name = "la16"
      case ("la18")
         allocate(g(18))
         g = [ &
            0.0010694900326538_dp, &
            -0.0004731544985879_dp, &
            -0.0102640640276849_dp, &
            0.0088592674935117_dp, &
            0.0620777893027638_dp, &
            -0.0182337707798257_dp, &
            -0.1915508312964873_dp, &
            0.0352724880359345_dp, &
            0.6173384491413523_dp, &
            0.7178970827642257_dp, &
            0.2387609146074182_dp, &
            -0.0545689584305765_dp, &
            0.0005834627463312_dp, &
            0.0302248788579895_dp, &
            -0.0115282102079848_dp, &
            -0.0132719677815332_dp, &
            0.0006197808890549_dp, &
            0.0014009155255716_dp ]
         filter_class = "Least Asymmetric"
         canonical_name = "la18"
      case ("la20")
         allocate(g(20))
         g = [ &
            0.0007701598091030_dp, &
            0.0000956326707837_dp, &
            -0.0086412992759401_dp, &
            -0.0014653825833465_dp, &
            0.0459272392237649_dp, &
            0.0116098939129724_dp, &
            -0.1594942788575307_dp, &
            -0.0708805358108615_dp, &
            0.4716906668426588_dp, &
            0.7695100370143388_dp, &
            0.3838267612253823_dp, &
            -0.0355367403054689_dp, &
            -0.0319900568281631_dp, &
            0.0499949720791560_dp, &
            0.0057649120455518_dp, &
            -0.0203549398039460_dp, &
            -0.0008043589345370_dp, &
            0.0045931735836703_dp, &
            0.0000570360843390_dp, &
            -0.0004593294205481_dp ]
         filter_class = "Least Asymmetric"
         canonical_name = "la20"
      case ("bl14")
         allocate(g(14))
         g = [ &
            0.0120154192834842_dp, &
            0.0172133762994439_dp, &
            -0.0649080035533744_dp, &
            -0.0641312898189170_dp, &
            0.3602184608985549_dp, &
            0.7819215932965554_dp, &
            0.4836109156937821_dp, &
            -0.0568044768822707_dp, &
            -0.1010109208664125_dp, &
            0.0447423494687405_dp, &
            0.0204642075778225_dp, &
            -0.0181266051311065_dp, &
            -0.0032832978473081_dp, &
            0.0022918339541009_dp ]
         filter_class = "Best Localized"
         canonical_name = "bl14"
      case ("bl18")
         allocate(g(18))
         g = [ &
            0.0002594576266544_dp, &
            -0.0006273974067728_dp, &
            -0.0019161070047557_dp, &
            0.0059845525181721_dp, &
            0.0040676562965785_dp, &
            -0.0295361433733604_dp, &
            -0.0002189514157348_dp, &
            0.0856124017265279_dp, &
            -0.0211480310688774_dp, &
            -0.1432929759396520_dp, &
            0.2337782900224977_dp, &
            0.7374707619933686_dp, &
            0.5926551374433956_dp, &
            0.0805670008868546_dp, &
            -0.1143343069619310_dp, &
            -0.0348460237698368_dp, &
            0.0139636362487191_dp, &
            0.0057746045512475_dp ]
         filter_class = "Best Localized"
         canonical_name = "bl18"
      case ("bl20")
         allocate(g(20))
         g = [ &
            0.0008625782242896_dp, &
            0.0007154205305517_dp, &
            -0.0070567640909701_dp, &
            0.0005956827305406_dp, &
            0.0496861265075979_dp, &
            0.0262403647054251_dp, &
            -0.1215521061578162_dp, &
            -0.0150192395413644_dp, &
            0.5137098728334054_dp, &
            0.7669548365010849_dp, &
            0.3402160135110789_dp, &
            -0.0878787107378667_dp, &
            -0.0670899071680668_dp, &
            0.0338423550064691_dp, &
            -0.0008687519578684_dp, &
            -0.0230054612862905_dp, &
            -0.0011404297773324_dp, &
            0.0050716491945793_dp, &
            0.0003401492622332_dp, &
            -0.0004101159165852_dp ]
         filter_class = "Best Localized"
         canonical_name = "bl20"
      case ("c6")
         allocate(g(6))
         g = [ &
            -0.0156557285289848_dp, &
            -0.0727326213410511_dp, &
            0.3848648565381134_dp, &
            0.8525720416423900_dp, &
            0.3378976709511590_dp, &
            -0.0727322757411889_dp ]
         filter_class = "Coiflet"
         canonical_name = "c6"
      case ("c12")
         allocate(g(12))
         g = [ &
            -0.0007205494453679_dp, &
            -0.0018232088707116_dp, &
            0.0056114348194211_dp, &
            0.0236801719464464_dp, &
            -0.0594344186467388_dp, &
            -0.0764885990786692_dp, &
            0.4170051844236707_dp, &
            0.8127236354493977_dp, &
            0.3861100668229939_dp, &
            -0.0673725547222826_dp, &
            -0.0414649367819558_dp, &
            0.0163873364635998_dp ]
         filter_class = "Coiflet"
         canonical_name = "c12"
      case ("c18")
         allocate(g(18))
         g = [ &
            -0.0000345997728362_dp, &
            -0.0000709833031381_dp, &
            0.0004662169601129_dp, &
            0.0011175187708906_dp, &
            -0.0025745176887502_dp, &
            -0.0090079761366615_dp, &
            0.0158805448636158_dp, &
            0.0345550275730615_dp, &
            -0.0823019271068856_dp, &
            -0.0717998216193117_dp, &
            0.4284834763776168_dp, &
            0.7937772226256169_dp, &
            0.4051769024096150_dp, &
            -0.0611233900026726_dp, &
            -0.0657719112818552_dp, &
            0.0234526961418362_dp, &
            0.0077825964273254_dp, &
            -0.0037935128644910_dp ]
         filter_class = "Coiflet"
         canonical_name = "c18"
      case ("c24")
         allocate(g(24))
         g = [ &
            -0.0000017849850031_dp, &
            -0.0000032596802369_dp, &
            0.0000312298758654_dp, &
            0.0000623390344610_dp, &
            -0.0002599745524878_dp, &
            -0.0005890207562444_dp, &
            0.0012665619292991_dp, &
            0.0037514361572790_dp, &
            -0.0056582866866115_dp, &
            -0.0152117315279485_dp, &
            0.0250822618448678_dp, &
            0.0393344271233433_dp, &
            -0.0962204420340021_dp, &
            -0.0666274742634348_dp, &
            0.4343860564915321_dp, &
            0.7822389309206135_dp, &
            0.4153084070304910_dp, &
            -0.0560773133167630_dp, &
            -0.0812666996808907_dp, &
            0.0266823001560570_dp, &
            0.0160689439647787_dp, &
            -0.0073461663276432_dp, &
            -0.0016294920126020_dp, &
            0.0008923136685824_dp ]
         filter_class = "Coiflet"
         canonical_name = "c24"
      case ("c30")
         allocate(g(30))
         g = [ &
            -0.0000000951765727_dp, &
            -0.0000001674428858_dp, &
            0.0000020637618516_dp, &
            0.0000037346551755_dp, &
            -0.0000213150268122_dp, &
            -0.0000413404322768_dp, &
            0.0001405411497166_dp, &
            0.0003022595818445_dp, &
            -0.0006381313431115_dp, &
            -0.0016628637021860_dp, &
            0.0024333732129107_dp, &
            0.0067641854487565_dp, &
            -0.0091642311634348_dp, &
            -0.0197617789446276_dp, &
            0.0326835742705106_dp, &
            0.0412892087544753_dp, &
            -0.1055742087143175_dp, &
            -0.0620359639693546_dp, &
            0.4379916262173834_dp, &
            0.7742896037334738_dp, &
            0.4215662067346898_dp, &
            -0.0520431631816557_dp, &
            -0.0919200105692549_dp, &
            0.0281680289738655_dp, &
            0.0234081567882734_dp, &
            -0.0101311175209033_dp, &
            -0.0041593587818186_dp, &
            0.0021782363583355_dp, &
            0.0003585896879330_dp, &
            -0.0002120808398259_dp ]
         filter_class = "Coiflet"
         canonical_name = "c30"
      case default
         allocate(g(0))
         filter_class = ""
         canonical_name = ""
         ierr = 1
      end select
   end subroutine named_scaling_coefficients

end module wavelets_filters
