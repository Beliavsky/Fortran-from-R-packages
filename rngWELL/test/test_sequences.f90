program test_sequences
   use, intrinsic :: iso_fortran_env, only : int32, int64, real64
   use rngwell, only : well_rng, well_from_options, well_state_size
   implicit none
   type(well_rng) :: a, b
   integer(int32), allocatable :: st(:)
   integer(int32) :: xa, xb
   real(real64) :: m(4,3), v(12)
   integer :: i

   call check('512a', [ &
      int(z'3e099644',int32), int(z'1ed900fc',int32), &
      int(z'95171007',int32), int(z'0941445e',int32), &
      int(z'c9b16d1f',int32), int(z'e05ebe34',int32), &
      int(z'13d573b7',int32), int(z'2667c974',int32)])
   call check('521a', [ &
      int(z'ee9a28c9',int32), int(z'471ed5fe',int32), &
      int(z'caee6e48',int32), int(z'40316a4d',int32), &
      int(z'a5ebc9c8',int32), int(z'cefab14b',int32), &
      int(z'5462f5d7',int32), int(z'39e6851d',int32)])
   call check('521b', [ &
      int(z'5928b748',int32), int(z'b7d92ee2',int32), &
      int(z'686e761d',int32), int(z'e508dca2',int32), &
      int(z'c086d916',int32), int(z'1ecd1bf8',int32), &
      int(z'abad3370',int32), int(z'bf9c43a2',int32)])
   call check('607a', [ &
      int(z'cc193cdf',int32), int(z'fa48863a',int32), &
      int(z'85497bfe',int32), int(z'9f60db8e',int32), &
      int(z'c2073d64',int32), int(z'5b565d92',int32), &
      int(z'517b1155',int32), int(z'bb65129d',int32)])
   call check('607b', [ &
      int(z'24c7560e',int32), int(z'8bade9d5',int32), &
      int(z'412f5616',int32), int(z'99f3b20a',int32), &
      int(z'2e175b96',int32), int(z'112fb2ea',int32), &
      int(z'2f689890',int32), int(z'e1200c9a',int32)])
   call check('800a', [ &
      int(z'c953f217',int32), int(z'da4956e6',int32), &
      int(z'18c78543',int32), int(z'4789f006',int32), &
      int(z'f57fa4ce',int32), int(z'e08f2233',int32), &
      int(z'728a2d0e',int32), int(z'21200b25',int32)])
   call check('800b', [ &
      int(z'b493e85a',int32), int(z'd5f430e8',int32), &
      int(z'e43a4128',int32), int(z'34edcc45',int32), &
      int(z'9c852846',int32), int(z'4ec349e0',int32), &
      int(z'c7c7e35c',int32), int(z'f7a5623f',int32)])
   call check('1024a', [ &
      int(z'51443936',int32), int(z'b13fde12',int32), &
      int(z'3f52c5f6',int32), int(z'fe7c6aa9',int32), &
      int(z'3ac732e3',int32), int(z'3bf073a8',int32), &
      int(z'b9690b1a',int32), int(z'b9c33f1d',int32)])
   call check('1024b', [ &
      int(z'6a0fb81c',int32), int(z'4e2ba1f3',int32), &
      int(z'4bdb2bb2',int32), int(z'900df389',int32), &
      int(z'f22408fb',int32), int(z'ba764d47',int32), &
      int(z'2e9e1122',int32), int(z'2af417fc',int32)])
   call check('19937a', [ &
      int(z'3010eeb0',int32), int(z'55e4fcf6',int32), &
      int(z'e76d64a0',int32), int(z'ab403f42',int32), &
      int(z'2d0679b9',int32), int(z'5c7bba29',int32), &
      int(z'37b4da9f',int32), int(z'70653dc7',int32)])
   call check('19937b', [ &
      int(z'0c657908',int32), int(z'3582a1e2',int32), &
      int(z'd45ce4d3',int32), int(z'8ee8a4a8',int32), &
      int(z'899f037c',int32), int(z'a8d43bfb',int32), &
      int(z'8dd570e4',int32), int(z'4da3aa9b',int32)])
   call check('19937c', [ &
      int(z'2b76feb0',int32), int(z'a608eff6',int32), &
      int(z'd94f74a0',int32), int(z'104e3e42',int32), &
      int(z'bfaeedb9',int32), int(z'eb332e29',int32), &
      int(z'fd5e5d9f',int32), int(z'cb69bec7',int32)])
   call check('21701a', [ &
      int(z'567a235e',int32), int(z'64a22fb3',int32), &
      int(z'9e4a7ee9',int32), int(z'38406cec',int32), &
      int(z'd5d48cb8',int32), int(z'8f9b79ae',int32), &
      int(z'992b4be3',int32), int(z'e89fd572',int32)])
   call check('23209a', [ &
      int(z'22101679',int32), int(z'1506e42f',int32), &
      int(z'060eb2be',int32), int(z'60372ac9',int32), &
      int(z'a78fe909',int32), int(z'59f2a57f',int32), &
      int(z'1d12d358',int32), int(z'91f2cc9c',int32)])
   call check('23209b', [ &
      int(z'15624142',int32), int(z'6f6306ed',int32), &
      int(z'401f7fae',int32), int(z'5c25fb89',int32), &
      int(z'53b1128d',int32), int(z'43030b9e',int32), &
      int(z'fc857645',int32), int(z'7f0b6484',int32)])
   call check('44497a', [ &
      int(z'88c1af99',int32), int(z'e997ae96',int32), &
      int(z'0afe95b7',int32), int(z'e7ff851d',int32), &
      int(z'bfc8a9b3',int32), int(z'fac6897c',int32), &
      int(z'38cb8aba',int32), int(z'99cbe765',int32)])
   call check('44497b', [ &
      int(z'd8142b99',int32), int(z'3843ae96',int32), &
      int(z'5ba705b7',int32), int(z'b43f011d',int32), &
      int(z'678d39b3',int32), int(z'b3929d7c',int32), &
      int(z'731f9eba',int32), int(z'620a7765',int32)])


   ! Canonical state export/import must continue the exact stream.
   call a%init('19937c', seed=7_int64)
   do i=1,1000
      xa=a%next_uint32()
   end do
   allocate(st(well_state_size('19937c')))
   call a%get_state(st)
   call b%init('19937c', state=st)
   do i=1,200
      xa=a%next_uint32()
      xb=b%next_uint32()
      if (xa /= xb) error stop 'state round-trip failed'
   end do
   deallocate(st)

   ! R-compatible option mapping for the special tempered variants.
   a=well_from_options(19937,'a',.true.,123_int64)
   if (trim(a%get_variant()) /= '19937c') error stop '19937 temper mapping failed'
   a=well_from_options(44497,'a',.true.,123_int64)
   if (trim(a%get_variant()) /= '44497b') error stop '44497 temper mapping failed'
   a=well_from_options(800,'b',.true.,123_int64)
   if (trim(a%get_variant()) /= '800b') error stop '800 mapping failed'

   ! Matrix fill follows Fortran/R column-major order.
   call a%init('512a', seed=19_int64)
   call a%fill_matrix(m)
   call b%init('512a', seed=19_int64)
   call b%fill(v)
   if (maxval(abs(reshape(m,[12])-v)) > epsilon(1.0_real64)) error stop 'matrix fill order failed'

   print '(a)', 'rngWELL sequence tests: PASS'

contains

   subroutine check(name, expected)
      character(len=*), intent(in) :: name
      integer(int32), intent(in) :: expected(:)
      type(well_rng) :: r
      integer(int32) :: got
      integer :: j
      call r%init(name, seed=123456789_int64)
      do j=1,size(expected)
         got=r%next_uint32()
         if (got /= expected(j)) then
            write(*,'(a,1x,a,1x,i0,2(1x,z8.8))') 'sequence mismatch',trim(name),j,got,expected(j)
            error stop 'sequence mismatch'
         end if
      end do
   end subroutine check

end program test_sequences
