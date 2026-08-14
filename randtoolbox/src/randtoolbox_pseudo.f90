module randtoolbox_pseudo
   use, intrinsic :: iso_fortran_env, only : int64, real64
   use randtoolbox_base, only : system_seed
   use randtoolbox_congruential, only : congru_rng
   use randtoolbox_sfmt, only : sfmt_rng, sfmt_random
   use randtoolbox_mt19937, only : mt19937_rng
   use randtoolbox_knuth, only : knuth_rng
   use rngwell, only : well_rng, well_from_options
   implicit none
   private
   integer(int64), save :: package_seed=1_int64
   logical, save :: seed_set=.false.
   public :: set_seed, congru_rand, sfmt_generate, well_generate, knuth_taocp, mt19937_generate
contains
   subroutine set_seed(seed)
      integer(int64), intent(in) :: seed
      package_seed=seed; seed_set=.true.
   end subroutine set_seed

   integer(int64) function choose_seed(seed) result(s)
      integer(int64), intent(in), optional :: seed
      if(present(seed)) then
         s=seed
      else if(seed_set) then
         s=package_seed; seed_set=.false.
      else
         s=system_seed()
      end if
   end function choose_seed

   function congru_rand(n,dim,modulus,multiplier,increment,seed) result(x)
      integer, intent(in) :: n,dim
      integer(int64), intent(in), optional :: modulus,multiplier,increment,seed
      real(real64), allocatable :: x(:,:)
      type(congru_rng) :: g
      integer(int64) :: m,a,c,s
      if(n<0 .or. dim<1) error stop 'randtoolbox: invalid congruRand dimensions'
      m=2147483647_int64; a=16807_int64; c=0_int64
      if(present(modulus))m=modulus; if(present(multiplier))a=multiplier; if(present(increment))c=increment
      s=choose_seed(seed); call g%init(s,m,a,c); allocate(x(n,dim)); call g%fill_matrix(x)
   end function congru_rand

   function sfmt_generate(n,dim,mexp,seed,use_parameter_sets) result(x)
      integer, intent(in) :: n,dim
      integer, intent(in), optional :: mexp
      integer(int64), intent(in), optional :: seed
      logical, intent(in), optional :: use_parameter_sets
      real(real64), allocatable :: x(:,:)
      integer :: me; integer(int64)::s; logical::ups
      me=19937; if(present(mexp))me=mexp; s=choose_seed(seed); ups=.true.; if(present(use_parameter_sets))ups=use_parameter_sets
      x=sfmt_random(n,dim,me,s,ups)
   end function sfmt_generate

   function well_generate(n,dim,order,temper,version,seed) result(x)
      integer, intent(in) :: n,dim
      integer, intent(in), optional :: order
      logical, intent(in), optional :: temper
      character(len=*), intent(in), optional :: version
      integer(int64), intent(in), optional :: seed
      real(real64), allocatable :: x(:,:)
      type(well_rng)::g
      integer::ord; logical::t; character(len=1)::v; integer(int64)::s
      if(n<0 .or. dim<1) error stop 'randtoolbox: invalid WELL dimensions'
      ord=512; if(present(order))ord=order; t=.false.; if(present(temper))t=temper; v='a'; if(present(version))v=version(1:1)
      s=choose_seed(seed); g=well_from_options(ord,v,t,s); allocate(x(n,dim)); call g%fill_matrix(x)
   end function well_generate

   function knuth_taocp(n,dim,seed) result(x)
      integer,intent(in)::n,dim
      integer(int64),intent(in),optional::seed
      real(real64),allocatable::x(:,:)
      type(knuth_rng)::g; integer(int64)::s
      if(n<0 .or. dim<1) error stop 'randtoolbox: invalid Knuth dimensions'
      s=choose_seed(seed); call g%seed(s); allocate(x(n,dim)); call g%fill_matrix(x)
   end function knuth_taocp

   function mt19937_generate(n,dim,seed,resolution,array_seed) result(x)
      integer,intent(in)::n,dim
      integer(int64),intent(in),optional::seed,array_seed(:)
      integer,intent(in),optional::resolution
      real(real64),allocatable::x(:,:)
      type(mt19937_rng)::g; integer(int64)::s; integer::r
      if(n<0 .or. dim<1) error stop 'randtoolbox: invalid MT19937 dimensions'
      r=32; if(present(resolution))r=resolution
      if(present(array_seed)) then; call g%init_array(array_seed); else; s=choose_seed(seed); call g%init(s); end if
      allocate(x(n,dim)); call g%fill_matrix(x,r)
   end function mt19937_generate
end module randtoolbox_pseudo
