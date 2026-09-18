module randompack_engines
   use iso_fortran_env, only : int64, int32
   use randompack_uint, only : uadd, usub, umul, umul_wide, urotl, add128, mul128_mod
   use randompack_uint, only : u64_from_u32, lo32, hi32, unsigned_lt
   implicit none
   private

   integer, parameter, public :: eng_x256ppsimd = 1
   integer, parameter, public :: eng_x256sssimd = 2
   integer, parameter, public :: eng_sfc64simd = 3
   integer, parameter, public :: eng_x256pp = 4
   integer, parameter, public :: eng_x256ss = 5
   integer, parameter, public :: eng_x128p = 6
   integer, parameter, public :: eng_xoro = 7
   integer, parameter, public :: eng_pcg64 = 8
   integer, parameter, public :: eng_sfc64 = 9
   integer, parameter, public :: eng_squares = 10
   integer, parameter, public :: eng_philox = 11
   integer, parameter, public :: eng_cwg128 = 12
   integer, parameter, public :: eng_ranluxpp = 13
   integer, parameter, public :: eng_chacha20 = 14

   integer(int64), parameter :: mask32 = int(z'00000000FFFFFFFF', int64)
   integer(int64), parameter :: mask16 = int(z'000000000000FFFF', int64)

   type, public :: engine_state_type
      integer :: engine = eng_x256ppsimd
      integer(int64) :: state(9) = 0_int64
      integer(int64) :: simd(4, 8) = 0_int64
      integer :: simd_lane = 1
      integer(int64) :: chacha_cache(8) = 0_int64
      integer :: chacha_pos = 9
   end type engine_state_type

   public :: engine_from_name, engine_name, engine_description, engine_state_words
   public :: seed_engine, randomize_engine, next_u64, jump_engine, advance_pcg
   public :: set_engine_state, set_pcg_inc, set_cwg_weyl, set_sfc_abc
   public :: set_chacha_nonce, set_philox_key, set_squares_key

contains
   pure integer function engine_from_name(name) result(engine)
      character(len=*), intent(in) :: name !! Engine identifier; matching is case-insensitive and treats '-' as '_'.
      character(len=:), allocatable :: key
      integer :: i, c
      key = adjustl(trim(name))
      do i = 1, len(key)
         c = iachar(key(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) key(i:i) = achar(c + 32)
         if (key(i:i) == '-') key(i:i) = '_'
      end do
      select case (key)
      case ('', 'x256++simd')
         engine = eng_x256ppsimd
      case ('x256**simd')
         engine = eng_x256sssimd
      case ('sfc64simd')
         engine = eng_sfc64simd
      case ('x256++')
         engine = eng_x256pp
      case ('x256**')
         engine = eng_x256ss
      case ('x128+')
         engine = eng_x128p
      case ('xoro++')
         engine = eng_xoro
      case ('pcg64')
         engine = eng_pcg64
      case ('sfc64')
         engine = eng_sfc64
      case ('squares')
         engine = eng_squares
      case ('philox')
         engine = eng_philox
      case ('cwg128')
         engine = eng_cwg128
      case ('ranlux++')
         engine = eng_ranluxpp
      case ('chacha20')
         engine = eng_chacha20
      case default
         engine = 0
      end select
   end function engine_from_name

   pure function engine_name(engine) result(name)
      integer, intent(in) :: engine !! Internal engine identifier in 1..14.
      character(len=12) :: name
      select case (engine)
      case (eng_x256ppsimd); name = 'x256++simd'
      case (eng_x256sssimd); name = 'x256**simd'
      case (eng_sfc64simd); name = 'sfc64simd'
      case (eng_x256pp); name = 'x256++'
      case (eng_x256ss); name = 'x256**'
      case (eng_x128p); name = 'x128+'
      case (eng_xoro); name = 'xoro++'
      case (eng_pcg64); name = 'pcg64'
      case (eng_sfc64); name = 'sfc64'
      case (eng_squares); name = 'squares'
      case (eng_philox); name = 'philox'
      case (eng_cwg128); name = 'cwg128'
      case (eng_ranluxpp); name = 'ranlux++'
      case (eng_chacha20); name = 'chacha20'
      case default; name = 'invalid'
      end select
   end function engine_name

   pure function engine_description(engine) result(description)
      integer, intent(in) :: engine !! Internal engine identifier in 1..14.
      character(len=72) :: description
      select case (engine)
      case (eng_x256ppsimd); description = 'xoshiro256++, portable scalar emulation of 8 streams'
      case (eng_x256sssimd); description = 'xoshiro256**, portable scalar emulation of 8 streams'
      case (eng_sfc64simd); description = 'sfc64, portable scalar emulation of 8 streams'
      case (eng_x256pp); description = 'xoshiro256++, Vigna and Blackman, 2019'
      case (eng_x256ss); description = 'xoshiro256**, Vigna and Blackman, 2019'
      case (eng_x128p); description = 'xorshift128+, Vigna, 2014'
      case (eng_xoro); description = 'xoroshiro128++, Vigna and Blackman, 2016'
      case (eng_pcg64); description = 'PCG64-DXSM, O''Neill, 2014'
      case (eng_sfc64); description = 'sfc64, Chris Doty-Humphrey, 2013'
      case (eng_squares); description = 'squares64, Widynski, 2021'
      case (eng_philox); description = 'Philox-4x64-10, Salmon and Moraes, 2011'
      case (eng_cwg128); description = 'cwg128, Dziala, 2022'
      case (eng_ranluxpp); description = 'ranlux++, Sibidanov, 2017'
      case (eng_chacha20); description = 'ChaCha20, Bernstein, 2008'
      case default; description = 'invalid engine'
      end select
   end function engine_description

   pure integer function engine_state_words(engine) result(n)
      integer, intent(in) :: engine !! Internal engine identifier in 1..14.
      select case (engine)
      case (eng_x256ppsimd, eng_x256sssimd, eng_sfc64simd, eng_x256pp, eng_x256ss, eng_pcg64, eng_sfc64)
         n = 4
      case (eng_x128p, eng_xoro, eng_squares)
         n = 2
      case (eng_philox, eng_chacha20)
         n = 6
      case (eng_cwg128)
         n = 8
      case (eng_ranluxpp)
         n = 9
      case default
         n = 0
      end select
   end function engine_state_words

   pure elemental function u32_mul(a, b) result(c)
      integer(int64), intent(in) :: a !! First unsigned 32-bit factor in an int64 container.
      integer(int64), intent(in) :: b !! Second unsigned 32-bit factor in an int64 container.
      integer(int64) :: c
      c = lo32(umul(a, b))
   end function u32_mul

   pure elemental function u32_add(a, b) result(c)
      integer(int64), intent(in) :: a !! First unsigned 32-bit addend in an int64 container.
      integer(int64), intent(in) :: b !! Second unsigned 32-bit addend in an int64 container.
      integer(int64) :: c
      c = lo32(a + b)
   end function u32_add

   pure elemental function u32_sub(a, b) result(c)
      integer(int64), intent(in) :: a !! Unsigned 32-bit minuend in an int64 container.
      integer(int64), intent(in) :: b !! Unsigned 32-bit subtrahend in an int64 container.
      integer(int64) :: c
      c = lo32(a - b)
   end function u32_sub

   pure subroutine seed_hash(value, hash_const, hashed)
      integer(int64), intent(in) :: value !! Unsigned 32-bit entropy word.
      integer(int64), intent(inout) :: hash_const !! Running unsigned 32-bit hash multiplier.
      integer(int64), intent(out) :: hashed !! Mixed unsigned 32-bit output word.
      integer(int64), parameter :: mult_a = int(z'931E8875', int64)
      integer(int64) :: v
      v = ieor(lo32(value), lo32(hash_const))
      hash_const = u32_mul(hash_const, mult_a)
      v = u32_mul(v, hash_const)
      hashed = ieor(v, shiftr(v, 16))
      hashed = lo32(hashed)
   end subroutine seed_hash

   pure function seed_mix(x, y) result(r)
      integer(int64), intent(in) :: x !! First unsigned 32-bit seed-mixing word.
      integer(int64), intent(in) :: y !! Second unsigned 32-bit seed-mixing word.
      integer(int64) :: r
      integer(int64), parameter :: ml = int(z'CA01F9DD', int64)
      integer(int64), parameter :: mr = int(z'4973F715', int64)
      r = u32_sub(u32_mul(ml, x), u32_mul(mr, y))
      r = lo32(ieor(r, shiftr(r, 16)))
   end function seed_mix

   pure subroutine seed_generate(seed32, spawn_key, out)
      integer(int64), intent(in) :: seed32 !! Deterministic 32-bit seed bit pattern.
      integer(int64), intent(in) :: spawn_key(:) !! Optional spawn-key words; only low 32 bits are used.
      integer(int64), intent(out) :: out(:) !! Generated unsigned 32-bit words, one per element.
      integer(int64), parameter :: init_a = int(z'43B0D7E5', int64)
      integer(int64), parameter :: init_b = int(z'8B51F9DD', int64)
      integer(int64), parameter :: mult_b = int(z'58F38DED', int64)
      integer(int64), allocatable :: seq(:)
      integer(int64) :: mixer(4), hash_const, h, v
      integer :: current, di, i, si, src
      allocate(seq(size(spawn_key) + 4))
      seq(1) = 314159265_int64
      seq(2) = lo32(seed32)
      seq(3) = int(size(spawn_key) + 1, int64)
      seq(4) = 271828182_int64
      if (size(spawn_key) > 0) seq(5:) = lo32(spawn_key)
      hash_const = init_a
      current = 1
      do i = 1, 4
         v = 0_int64
         if (current <= size(seq)) then
            v = seq(current)
            current = current + 1
         end if
         call seed_hash(v, hash_const, mixer(i))
      end do
      do si = 1, 4
         v = mixer(si)
         do di = 1, 4
            if (di == si) cycle
            call seed_hash(v, hash_const, h)
            mixer(di) = seed_mix(mixer(di), h)
         end do
      end do
      do while (current <= size(seq))
         v = seq(current)
         current = current + 1
         do di = 1, 4
            call seed_hash(v, hash_const, h)
            mixer(di) = seed_mix(mixer(di), h)
         end do
      end do
      hash_const = init_b
      src = 1
      do i = 1, size(out)
         v = ieor(mixer(src), hash_const)
         hash_const = u32_mul(hash_const, mult_b)
         v = u32_mul(v, hash_const)
         out(i) = lo32(ieor(v, shiftr(v, 16)))
         src = src + 1
         if (src == 5) src = 1
      end do
   end subroutine seed_generate

   pure subroutine seed_engine(state, seed, spawn_key, info)
      type(engine_state_type), intent(inout) :: state !! Engine state to deterministically initialize.
      integer, intent(in) :: seed !! Signed Fortran seed whose low 32 bits match the upstream R integer seed.
      integer(int64), intent(in), optional :: spawn_key(:) !! Optional unsigned 32-bit spawn-key words.
      integer, intent(out) :: info !! Zero on success; nonzero for an invalid engine.
      integer(int64), allocatable :: words32(:), key(:)
      integer :: i, n
      n = engine_state_words(state%engine)
      if (n <= 0) then
         info = 1
         return
      end if
      if (present(spawn_key)) then
         allocate(key(size(spawn_key)))
         key = spawn_key
      else
         allocate(key(0))
      end if
      allocate(words32(2 * n))
      call seed_generate(iand(int(seed, int64), mask32), key, words32)
      state%state = 0_int64
      do i = 1, n
         state%state(i) = u64_from_u32(words32(2 * i - 1), words32(2 * i))
      end do
      state%simd_lane = 1
      state%chacha_pos = 9
      if (state%engine == eng_ranluxpp) state%chacha_pos = 10
      if (state%engine == eng_x256ppsimd .or. state%engine == eng_x256sssimd) then
         call setup_x256_streams(state)
      else if (state%engine == eng_sfc64simd) then
         call setup_sfc_streams(state)
      end if
      info = 0
   end subroutine seed_engine

   subroutine randomize_engine(state)
      type(engine_state_type), intent(inout) :: state !! Engine state to initialize from portable clock/date entropy.
      integer :: values(8), count, rate, i, n, info
      integer(int64) :: x, z, key(4)
      call date_and_time(values=values)
      call system_clock(count, rate)
      x = int(count, int64)
      do i = 1, 8
         x = ieor(x, shiftl(int(values(i), int64), modulo(7 * i, 48)))
      end do
      do i = 1, 4
         x = uadd(x, int(z'9E3779B97F4A7C15', int64))
         z = x
         z = umul(ieor(z, shiftr(z, 30)), int(z'BF58476D1CE4E5B9', int64))
         z = umul(ieor(z, shiftr(z, 27)), int(z'94D049BB133111EB', int64))
         key(i) = ieor(z, shiftr(z, 31))
      end do
      n = engine_state_words(state%engine)
      call seed_engine(state, int(lo32(key(1)), kind=kind(1)), key(2:min(4, n + 1)), info)
   end subroutine randomize_engine

   pure subroutine setup_x256_streams(st)
      type(engine_state_type), intent(inout) :: st !! SIMD-compatible xoshiro state whose stream 0 is in `state(1:4)`.
      integer :: lane
      st%simd(:, 1) = st%state(1:4)
      do lane = 2, 8
         st%simd(:, lane) = st%simd(:, lane - 1)
         call xoshiro_jump_state(st%simd(:, lane), 253)
      end do
      st%simd_lane = 1
   end subroutine setup_x256_streams

   pure subroutine setup_sfc_streams(st)
      type(engine_state_type), intent(inout) :: st !! SIMD-compatible SFC state whose base words are in `state(1:4)`.
      integer :: lane
      integer(int64), parameter :: delta = int(z'2000000000000000', int64)
      do lane = 1, 8
         st%simd(1:3, lane) = st%state(1:3)
         st%simd(4, lane) = uadd(st%state(4), umul(int(lane - 1, int64), delta))
      end do
      st%simd_lane = 1
   end subroutine setup_sfc_streams

   pure subroutine xoshiro_step(s)
      integer(int64), intent(inout) :: s(4) !! Four-word xoshiro256 state advanced by one transition.
      integer(int64) :: t
      t = shiftl(s(2), 17)
      s(3) = ieor(s(3), s(1))
      s(4) = ieor(s(4), s(2))
      s(2) = ieor(s(2), s(3))
      s(1) = ieor(s(1), s(4))
      s(3) = ieor(s(3), t)
      s(4) = urotl(s(4), 45)
   end subroutine xoshiro_step

   pure subroutine xoshiro_next(s, starstar, value)
      integer(int64), intent(inout) :: s(4) !! Four-word xoshiro256 state advanced in place.
      logical, intent(in) :: starstar !! Use xoshiro256** output when true, xoshiro256++ otherwise.
      integer(int64), intent(out) :: value !! Next 64-bit output word.
      if (starstar) then
         value = umul(urotl(umul(s(2), 5_int64), 7), 9_int64)
      else
         value = uadd(urotl(uadd(s(1), s(4)), 23), s(1))
      end if
      call xoshiro_step(s)
   end subroutine xoshiro_next

   pure subroutine x128_next(s, value)
      integer(int64), intent(inout) :: s(2) !! Two-word xorshift128+ state advanced in place.
      integer(int64), intent(out) :: value !! Next 64-bit output word.
      integer(int64) :: t, old0
      old0 = s(1)
      t = ieor(s(2), shiftl(s(2), 23))
      value = uadd(old0, t)
      s(2) = old0
      s(1) = ieor(ieor(t, old0), ieor(shiftr(t, 18), shiftr(old0, 5)))
   end subroutine x128_next

   pure subroutine xoro_next(s, value)
      integer(int64), intent(inout) :: s(2) !! Two-word xoroshiro128++ state advanced in place.
      integer(int64), intent(out) :: value !! Next 64-bit output word.
      integer(int64) :: s0, s1
      s0 = s(1)
      s1 = s(2)
      value = uadd(urotl(uadd(s0, s1), 17), s0)
      s1 = ieor(s1, s0)
      s(1) = ieor(ieor(urotl(s0, 49), s1), shiftl(s1, 21))
      s(2) = urotl(s1, 28)
   end subroutine xoro_next

   pure subroutine pcg_next(s, value)
      integer(int64), intent(inout) :: s(4) !! PCG64 state `(state_lo,state_hi,inc_lo,inc_hi)`.
      integer(int64), intent(out) :: value !! Next PCG64-DXSM output word.
      integer(int64), parameter :: multiplier = int(z'DA942042E4DD58B5', int64)
      integer(int64) :: state_lo, state_hi, prod_lo, prod_hi, next_lo, next_hi, hi, lo, hh
      state_lo = s(1)
      state_hi = s(2)
      call umul_wide(state_lo, multiplier, prod_hi, prod_lo)
      prod_hi = uadd(prod_hi, umul(state_hi, multiplier))
      call add128(prod_lo, prod_hi, s(3), s(4), next_lo, next_hi)
      hi = ieor(state_hi, shiftr(state_hi, 32))
      hi = umul(hi, multiplier)
      hi = ieor(hi, shiftr(hi, 48))
      lo = ior(state_lo, 1_int64)
      hh = umul(hi, lo)
      value = hh
      s(1) = next_lo
      s(2) = next_hi
   end subroutine pcg_next

   pure subroutine sfc_next(s, value)
      integer(int64), intent(inout) :: s(4) !! SFC64 state `(a,b,c,counter)` advanced in place.
      integer(int64), intent(out) :: value !! Next SFC64 output word.
      integer(int64) :: a, b, c, counter
      a = s(1)
      b = s(2)
      c = s(3)
      counter = s(4)
      value = uadd(uadd(a, b), counter)
      s(4) = uadd(counter, 1_int64)
      s(1) = ieor(b, shiftr(b, 11))
      s(2) = uadd(c, shiftl(c, 3))
      s(3) = uadd(urotl(c, 24), value)
   end subroutine sfc_next

   pure function squares_word(counter, key) result(value)
      integer(int64), intent(in) :: counter !! Squares64 counter word.
      integer(int64), intent(in) :: key !! Squares64 key word.
      integer(int64) :: value
      integer(int64) :: x, y, z, t
      y = umul(counter, key)
      x = y
      z = uadd(y, key)
      x = uadd(umul(x, x), y)
      x = urotl(x, 32)
      x = uadd(umul(x, x), z)
      x = urotl(x, 32)
      x = uadd(umul(x, x), y)
      x = urotl(x, 32)
      x = uadd(umul(x, x), z)
      t = x
      x = urotl(x, 32)
      value = ieor(t, shiftr(uadd(umul(x, x), y), 32))
   end function squares_word

   pure subroutine squares_next(s, value)
      integer(int64), intent(inout) :: s(2) !! Squares64 `(counter,key)` state advanced in place.
      integer(int64), intent(out) :: value !! Next Squares64 output word.
      value = squares_word(s(1), s(2))
      s(1) = uadd(s(1), 1_int64)
   end subroutine squares_next

   pure subroutine philox_round(counter, key, out)
      integer(int64), intent(in) :: counter(4) !! Four-word Philox counter for one round.
      integer(int64), intent(in) :: key(2) !! Two-word Philox key for one round.
      integer(int64), intent(out) :: out(4) !! Four-word output of one Philox round.
      integer(int64) :: hi0, hi1, lo0, lo1
      call umul_wide(int(z'D2E7470EE14C6C93', int64), counter(1), hi0, lo0)
      call umul_wide(int(z'CA5A826395121157', int64), counter(3), hi1, lo1)
      out(1) = ieor(ieor(hi1, counter(2)), key(1))
      out(2) = lo1
      out(3) = ieor(ieor(hi0, counter(4)), key(2))
      out(4) = lo0
   end subroutine philox_round

   pure subroutine philox_block(counter, key0, block)
      integer(int64), intent(in) :: counter(4) !! Four-word Philox block counter.
      integer(int64), intent(in) :: key0(2) !! Initial two-word Philox key.
      integer(int64), intent(out) :: block(4) !! Four-word Philox-4x64-10 block.
      integer(int64) :: ctr(4), key(2), tmp(4)
      integer :: round
      ctr = counter
      key = key0
      do round = 1, 10
         call philox_round(ctr, key, tmp)
         ctr = tmp
         if (round < 10) then
            key(1) = uadd(key(1), int(z'9E3779B97F4A7C15', int64))
            key(2) = uadd(key(2), int(z'BB67AE8584CAA73B', int64))
         end if
      end do
      block = ctr
   end subroutine philox_block

   pure subroutine philox_next(st, value)
      type(engine_state_type), intent(inout) :: st !! Philox state and four-word output cache position.
      integer(int64), intent(out) :: value !! Next Philox output word.
      integer :: pos
      ! state(7) is a local output position 1..4 encoded as a small integer; zero means refill.
      pos = int(iand(st%state(7), 7_int64))
      if (pos < 1 .or. pos > 4) then
         call philox_block(st%state(1:4), st%state(5:6), st%chacha_cache(1:4))
         st%state(1) = uadd(st%state(1), 1_int64)
         pos = 1
      end if
      value = st%chacha_cache(pos)
      pos = pos + 1
      if (pos > 4) pos = 0
      st%state(7) = int(pos, int64)
   end subroutine philox_next

   pure subroutine cwg_next_pair(s, out1, out2)
      integer(int64), intent(inout) :: s(8) !! CWG128 state as four low/high 128-bit pairs.
      integer(int64), intent(out) :: out1 !! Low 64 bits of the next CWG128 output pair.
      integer(int64), intent(out) :: out2 !! High 64 bits of the next CWG128 output pair.
      integer(int64) :: c0l, c0h, c1l, c1h, c2l, c2h, c3l, c3h
      integer(int64) :: t1l, t1h, t2l, t2h, zl, zh, cy
      c0l = s(1)
      c0h = s(2)
      c1l = s(3)
      c1h = s(4)
      c2l = s(5)
      c2h = s(6)
      c3l = s(7)
      c3h = s(8)
      call add128(c2l, c2h, c1l, c1h, t1l, t1h)
      c2l = t1l
      c2h = t1h
      t1l = ior(shiftr(c1l, 1), shiftl(c1h, 63))
      t1h = shiftr(c1h, 1)
      call mul128_mod(t1l, t1h, ior(c2l, 1_int64), c2h, t2l, t2h)
      call add128(c3l, c3h, c0l, c0h, zl, zh)
      c3l = zl
      c3h = zh
      c1l = ieor(t2l, c3l)
      c1h = ieor(t2h, c3h)
      cy = shiftr(c2h, 32)
      zl = ieor(c1l, cy)
      zh = c1h
      out1 = zl
      out2 = zh
      s = [c0l, c0h, c1l, c1h, c2l, c2h, c3l, c3h]
   end subroutine cwg_next_pair

   pure subroutine cwg_next(st, value)
      type(engine_state_type), intent(inout) :: st !! CWG128 state with a one-word pair cache.
      integer(int64), intent(out) :: value !! Next CWG128 64-bit word.
      integer(int64) :: a, b
      if (st%chacha_pos == 1) then
         value = st%chacha_cache(1)
         st%chacha_pos = 9
      else
         call cwg_next_pair(st%state(1:8), a, b)
         value = a
         st%chacha_cache(1) = b
         st%chacha_pos = 1
      end if
   end subroutine cwg_next

   pure subroutine addc64(a, b, cin, result, carry)
      integer(int64), intent(in) :: a !! First unsigned 64-bit addend.
      integer(int64), intent(in) :: b !! Second unsigned 64-bit addend.
      integer(int64), intent(in) :: cin !! Carry-in, constrained to zero or one.
      integer(int64), intent(out) :: result !! Unsigned sum modulo 2^64.
      integer(int64), intent(out) :: carry !! Carry-out, zero or one.
      integer(int64) :: t
      logical :: c0, c1
      t = uadd(a, b)
      c0 = unsigned_lt(t, a)
      result = uadd(t, cin)
      c1 = (cin /= 0_int64 .and. result == 0_int64)
      carry = merge(1_int64, 0_int64, c0 .or. c1)
   end subroutine addc64

   pure subroutine subb64(a, b, bin, result, borrow)
      integer(int64), intent(in) :: a !! Unsigned 64-bit minuend.
      integer(int64), intent(in) :: b !! Unsigned 64-bit subtrahend.
      integer(int64), intent(in) :: bin !! Borrow-in, constrained to zero or one.
      integer(int64), intent(out) :: result !! Difference modulo 2^64.
      integer(int64), intent(out) :: borrow !! Borrow-out, zero or one.
      integer(int64) :: t
      logical :: b0, b1
      b0 = unsigned_lt(a, b)
      t = usub(a, b)
      b1 = (bin /= 0_int64 .and. t == 0_int64)
      result = usub(t, bin)
      borrow = merge(1_int64, 0_int64, b0 .or. b1)
   end subroutine subb64

   pure subroutine mul9x9(z, x, y)
      integer(int64), intent(out) :: z(18) !! Eighteen-word unsigned product workspace.
      integer(int64), intent(in) :: x(9) !! First nine-word little-endian factor.
      integer(int64), intent(in) :: y(9) !! Second nine-word little-endian factor.
      integer(int64) :: lo, hi, c, cy, tmp, cy2
      integer :: i, j
      z = 0_int64
      do i = 1, 9
         c = 0_int64
         cy = 0_int64
         do j = 1, 9
            call umul_wide(x(i), y(j), hi, lo)
            call addc64(lo, c, 0_int64, tmp, cy2)
            lo = tmp
            hi = uadd(hi, cy2)
            c = hi
            call addc64(z(i + j - 1), lo, cy, tmp, cy2)
            z(i + j - 1) = tmp
            cy = cy2
         end do
         z(i + 9) = uadd(c, cy)
      end do
   end subroutine mul9x9

   pure subroutine sub_bw_n(x, y, n, borrow_out)
      integer(int64), intent(inout) :: x(9) !! Nine-word accumulator modified by subtraction.
      integer(int64), intent(in) :: y(9) !! Nine-word subtrahend; only first `n` words are used directly.
      integer, intent(in) :: n !! Number of low words from `y` participating in the subtraction.
      integer, intent(out) :: borrow_out !! Final unsigned borrow bit, returned as zero or one.
      integer(int64) :: borrow, next_borrow, tmp
      integer :: i
      borrow = 0_int64
      do i = 1, n
         call subb64(x(i), y(i), borrow, tmp, next_borrow)
         x(i) = tmp
         borrow = next_borrow
      end do
      if (borrow /= 0_int64) then
         do i = n + 1, 9
            call subb64(x(i), 0_int64, borrow, tmp, next_borrow)
            x(i) = tmp
            borrow = next_borrow
            if (borrow == 0_int64) exit
         end do
      end if
      borrow_out = int(borrow)
   end subroutine sub_bw_n

   pure subroutine add_cy_n(x, y, n, carry_out)
      integer(int64), intent(inout) :: x(9) !! Nine-word accumulator modified by addition.
      integer(int64), intent(in) :: y(9) !! Nine-word addend; only first `n` words are used directly.
      integer, intent(in) :: n !! Number of low words from `y` participating in the addition.
      integer, intent(out) :: carry_out !! Final unsigned carry bit, returned as zero or one.
      integer(int64) :: carry, next_carry, tmp
      integer :: i
      carry = 0_int64
      do i = 1, n
         call addc64(x(i), y(i), carry, tmp, next_carry)
         x(i) = tmp
         carry = next_carry
      end do
      if (carry /= 0_int64) then
         do i = n + 1, 9
            call addc64(x(i), 0_int64, carry, tmp, next_carry)
            x(i) = tmp
            carry = next_carry
            if (carry == 0_int64) exit
         end do
      end if
      carry_out = int(carry)
   end subroutine add_cy_n

   pure subroutine add1_words(x, first, value)
      integer(int64), intent(inout) :: x(9) !! Nine-word accumulator receiving a scalar addition.
      integer, intent(in) :: first !! One-based first word at which to add `value`.
      integer(int64), intent(in) :: value !! Unsigned scalar addend.
      integer(int64) :: carry, next_carry, tmp
      integer :: i
      call addc64(x(first), value, 0_int64, tmp, carry)
      x(first) = tmp
      do i = first + 1, 9
         if (carry == 0_int64) exit
         call addc64(x(i), 0_int64, carry, tmp, next_carry)
         x(i) = tmp
         carry = next_carry
      end do
   end subroutine add1_words

   pure subroutine sub1_words(x, first, value)
      integer(int64), intent(inout) :: x(9) !! Nine-word accumulator receiving a scalar subtraction.
      integer, intent(in) :: first !! One-based first word at which to subtract `value`.
      integer(int64), intent(in) :: value !! Unsigned scalar subtrahend.
      integer(int64) :: borrow, next_borrow, tmp
      integer :: i
      call subb64(x(first), value, 0_int64, tmp, borrow)
      x(first) = tmp
      do i = first + 1, 9
         if (borrow == 0_int64) exit
         call subb64(x(i), 0_int64, borrow, tmp, next_borrow)
         x(i) = tmp
         borrow = next_borrow
      end do
   end subroutine sub1_words

   pure subroutine ranlux_mod(x, z)
      integer(int64), intent(out) :: x(9) !! Reduced nine-word ranlux++ state.
      integer(int64), intent(in) :: z(18) !! Eighteen-word product to reduce modulo 2^576-2^240+1.
      integer(int64) :: t1(9), t2(9), t3(9), c0, c3
      integer :: c
      t1 = z(10:18)
      t2 = 0_int64
      t2(1) = ior(shiftr(t1(6), 16), shiftl(t1(7), 48))
      t2(2) = ior(shiftr(t1(7), 16), shiftl(t1(8), 48))
      t2(3) = ior(shiftr(t1(8), 16), shiftl(t1(9), 48))
      t2(4) = shiftr(t1(9), 16)
      t3 = 0_int64
      t3(1:6) = t1(1:6)
      t3(6) = iand(t3(6), mask16)
      x = z(1:9)
      c = 0
      block
         integer :: bit_result
         call sub_bw_n(x, t1, 9, bit_result)
         c = c - bit_result
         call sub_bw_n(x, t2, 5, bit_result)
         c = c - bit_result
         call add_cy_n(t3, t2, 5, bit_result)
         c = c + bit_result
      end block
      ! Shift t3 left by 240 bits = 3 words + 48 bits.
      t3(9) = ior(shiftl(t3(6), 48), shiftr(t3(5), 16))
      t3(8) = ior(shiftl(t3(5), 48), shiftr(t3(4), 16))
      t3(7) = ior(shiftl(t3(4), 48), shiftr(t3(3), 16))
      t3(6) = ior(shiftl(t3(3), 48), shiftr(t3(2), 16))
      t3(5) = ior(shiftl(t3(2), 48), shiftr(t3(1), 16))
      t3(4) = shiftl(t3(1), 48)
      t3(1:3) = 0_int64
      block
         integer :: bit_result
         call add_cy_n(x, t3, 9, bit_result)
         c = c + bit_result
      end block
      c0 = int(abs(c), int64)
      c3 = shiftl(c0, 48)
      if (c > 0) then
         call sub1_words(x, 1, c0)
         call add1_words(x, 4, c3)
      else if (c < 0) then
         call add1_words(x, 1, c0)
         call sub1_words(x, 4, c3)
      end if
   end subroutine ranlux_mod

   pure subroutine ranlux_step(s)
      integer(int64), intent(inout) :: s(9) !! Nine-word ranlux++ state advanced by one modular multiplication.
      integer(int64), parameter :: a(9) = [ &
         int(z'ED7FAA90747AAAD9', int64), int(z'4CEC2C78AF55C101', int64), &
         int(z'E64DCB31C48228EC', int64), int(z'6D8A15A13BEE7CB0', int64), &
         int(z'20B2CA60CB78C509', int64), int(z'256C3D3C662EA36C', int64), &
         int(z'FF74E54107684ED2', int64), int(z'492EDFCC0CC8E753', int64), &
         int(z'B48C187CF5B22097', int64)]
      integer(int64) :: z(18), out(9)
      call mul9x9(z, s, a)
      call ranlux_mod(out, z)
      s = out
   end subroutine ranlux_step

   pure subroutine ranlux_next(st, value)
      type(engine_state_type), intent(inout) :: st !! ranlux++ state with a nine-word block output position.
      integer(int64), intent(out) :: value !! Next ranlux++ 64-bit word.
      integer :: pos
      pos = st%chacha_pos
      if (pos < 1 .or. pos > 9) then
         call ranlux_step(st%state(1:9))
         st%chacha_cache(1:8) = st%state(1:8)
         pos = 1
      end if
      if (pos <= 8) then
         value = st%chacha_cache(pos)
      else
         value = st%state(9)
      end if
      pos = pos + 1
      if (pos > 9) pos = 10
      st%chacha_pos = pos
   end subroutine ranlux_next

   pure elemental function rotl32(x, k) result(y)
      integer(int64), intent(in) :: x !! Unsigned 32-bit word stored in an int64.
      integer, intent(in) :: k !! Left rotation distance in bits, in 0..31.
      integer(int64) :: y
      integer :: r
      r = modulo(k, 32)
      if (r == 0) then
         y = lo32(x)
      else
         y = lo32(ior(shiftl(lo32(x), r), shiftr(lo32(x), 32 - r)))
      end if
   end function rotl32

   pure subroutine chacha_qr(a, b, c, d)
      integer(int64), intent(inout) :: a !! First unsigned 32-bit ChaCha quarter-round word.
      integer(int64), intent(inout) :: b !! Second unsigned 32-bit ChaCha quarter-round word.
      integer(int64), intent(inout) :: c !! Third unsigned 32-bit ChaCha quarter-round word.
      integer(int64), intent(inout) :: d !! Fourth unsigned 32-bit ChaCha quarter-round word.
      a = u32_add(a, b)
      d = rotl32(ieor(d, a), 16)
      c = u32_add(c, d)
      b = rotl32(ieor(b, c), 12)
      a = u32_add(a, b)
      d = rotl32(ieor(d, a), 8)
      c = u32_add(c, d)
      b = rotl32(ieor(b, c), 7)
   end subroutine chacha_qr

   pure subroutine get_state_u32(s, w)
      integer(int64), intent(in) :: s(6) !! Six packed 64-bit ChaCha state words.
      integer(int64), intent(out) :: w(12) !! Twelve unpacked unsigned 32-bit state words.
      integer :: i
      do i = 1, 6
         w(2 * i - 1) = lo32(s(i))
         w(2 * i) = hi32(s(i))
      end do
   end subroutine get_state_u32

   pure subroutine set_state_u32(s, w)
      integer(int64), intent(inout) :: s(6) !! Six packed 64-bit ChaCha state words to replace.
      integer(int64), intent(in) :: w(12) !! Twelve unsigned 32-bit state words to pack.
      integer :: i
      do i = 1, 6
         s(i) = u64_from_u32(w(2 * i - 1), w(2 * i))
      end do
   end subroutine set_state_u32

   pure subroutine chacha_block(st)
      type(engine_state_type), intent(inout) :: st !! ChaCha20 engine state receiving an eight-word output block.
      integer(int64), parameter :: constants(4) = [int(z'61707865', int64), int(z'3320646E', int64), &
         int(z'79622D32', int64), int(z'6B206574', int64)]
      integer(int64) :: sw(12), base(16), x(16), out32(16)
      integer :: i, round
      call get_state_u32(st%state(1:6), sw)
      base(1:4) = constants
      base(5:12) = sw(1:8)
      base(13) = sw(12)
      base(14:16) = sw(9:11)
      x = base
      do round = 1, 10
         call chacha_qr(x(1), x(5), x(9), x(13))
         call chacha_qr(x(2), x(6), x(10), x(14))
         call chacha_qr(x(3), x(7), x(11), x(15))
         call chacha_qr(x(4), x(8), x(12), x(16))
         call chacha_qr(x(1), x(6), x(11), x(16))
         call chacha_qr(x(2), x(7), x(12), x(13))
         call chacha_qr(x(3), x(8), x(9), x(14))
         call chacha_qr(x(4), x(5), x(10), x(15))
      end do
      do i = 1, 16
         out32(i) = u32_add(x(i), base(i))
      end do
      do i = 1, 8
         st%chacha_cache(i) = u64_from_u32(out32(2 * i - 1), out32(2 * i))
      end do
      sw(12) = u32_add(sw(12), 1_int64)
      call set_state_u32(st%state(1:6), sw)
      st%chacha_pos = 1
   end subroutine chacha_block

   pure subroutine chacha_next(st, value)
      type(engine_state_type), intent(inout) :: st !! ChaCha20 engine state and block cache position.
      integer(int64), intent(out) :: value !! Next ChaCha20 output word in little-endian byte order.
      if (st%chacha_pos < 1 .or. st%chacha_pos > 8) call chacha_block(st)
      value = st%chacha_cache(st%chacha_pos)
      st%chacha_pos = st%chacha_pos + 1
   end subroutine chacha_next

   pure subroutine next_u64(st, value)
      type(engine_state_type), intent(inout) :: st !! RNG engine state advanced by one 64-bit output.
      integer(int64), intent(out) :: value !! Next raw 64-bit engine word.
      integer :: lane
      select case (st%engine)
      case (eng_x256ppsimd, eng_x256sssimd)
         lane = st%simd_lane
         call xoshiro_next(st%simd(:, lane), st%engine == eng_x256sssimd, value)
         st%simd_lane = lane + 1
         if (st%simd_lane > 8) st%simd_lane = 1
      case (eng_sfc64simd)
         lane = st%simd_lane
         call sfc_next(st%simd(:, lane), value)
         st%simd_lane = lane + 1
         if (st%simd_lane > 8) st%simd_lane = 1
      case (eng_x256pp)
         call xoshiro_next(st%state(1:4), .false., value)
      case (eng_x256ss)
         call xoshiro_next(st%state(1:4), .true., value)
      case (eng_x128p)
         call x128_next(st%state(1:2), value)
      case (eng_xoro)
         call xoro_next(st%state(1:2), value)
      case (eng_pcg64)
         call pcg_next(st%state(1:4), value)
      case (eng_sfc64)
         call sfc_next(st%state(1:4), value)
      case (eng_squares)
         call squares_next(st%state(1:2), value)
      case (eng_philox)
         call philox_next(st, value)
      case (eng_cwg128)
         call cwg_next(st, value)
      case (eng_ranluxpp)
         call ranlux_next(st, value)
      case (eng_chacha20)
         call chacha_next(st, value)
      case default
         value = 0_int64
      end select
   end subroutine next_u64

   pure subroutine xoshiro_jump_state(s, p)
      integer(int64), intent(inout) :: s(4) !! xoshiro256 state advanced by a supported power-of-two jump.
      integer, intent(in) :: p !! Jump exponent: 32, 64, 96, 128, 192, or internal 253.
      integer(int64) :: jump(4), acc(4)
      integer :: i, b
      select case (p)
      case (32)
         jump = [int(z'58120D583C112F69', int64), int(z'7D8D0632BD08E6AC', int64), &
            int(z'214FAFC0FBDBC208', int64), int(z'0E055D3520FDB9D7', int64)]
      case (64)
         jump = [int(z'B13C16E8096F0754', int64), int(z'B60D6C5B8C78F106', int64), &
            int(z'34FAFF184785C20A', int64), int(z'12E4A2FBFC19BFF9', int64)]
      case (96)
         jump = [int(z'148C356C3114B7A9', int64), int(z'CDB45D7DEF42C317', int64), &
            int(z'B27C05962EA56A13', int64), int(z'31EEBB6C82A9615F', int64)]
      case (128)
         jump = [int(z'180EC6D33CFD0ABA', int64), int(z'D5A61266F0C9392C', int64), &
            int(z'A9582618E03FC9AA', int64), int(z'39ABDC4529B1661C', int64)]
      case (192)
         jump = [int(z'76E15D3EFEFDCBBF', int64), int(z'C5004E441C522FB3', int64), &
            int(z'77710069854EE241', int64), int(z'39109BB02ACBE635', int64)]
      case (253)
         jump = [int(z'DFCA68648B28C5AF', int64), int(z'B56437FB2B753802', int64), &
            int(z'EBB82AACDF6CA80D', int64), int(z'A170E108788DB093', int64)]
      case default
         return
      end select
      acc = 0_int64
      do i = 1, 4
         do b = 0, 63
            if (btest(jump(i), b)) acc = ieor(acc, s)
            call xoshiro_step(s)
         end do
      end do
      s = acc
   end subroutine xoshiro_jump_state

   pure subroutine x128_jump_state(s, p)
      integer(int64), intent(inout) :: s(2) !! xorshift128+ state advanced by a supported jump.
      integer, intent(in) :: p !! Jump exponent: 32, 64, or 96.
      integer(int64) :: jump(2), acc(2), dummy
      integer :: i, b
      select case (p)
      case (32); jump = [int(z'55FCB2D81DA9D3CA', int64), int(z'43B6C9723B2AA348', int64)]
      case (64); jump = [int(z'8A5CD789635D2DFF', int64), int(z'121FD2155C472F96', int64)]
      case (96); jump = [int(z'EA61C9F1F13962AE', int64), int(z'A1FE50EF79CFAFB2', int64)]
      case default; return
      end select
      acc = 0_int64
      do i = 1, 2
         do b = 0, 63
            if (btest(jump(i), b)) acc = ieor(acc, s)
            call x128_next(s, dummy)
         end do
      end do
      s = acc
   end subroutine x128_jump_state

   pure subroutine xoro_jump_state(s, p)
      integer(int64), intent(inout) :: s(2) !! xoroshiro128++ state advanced by a supported jump.
      integer, intent(in) :: p !! Jump exponent: 32, 64, or 96.
      integer(int64) :: jump(2), acc(2), dummy
      integer :: i, b
      select case (p)
      case (32); jump = [int(z'FCCEEC21D5C306D9', int64), int(z'2E1BCF52F1051044', int64)]
      case (64); jump = [int(z'2BD7A6A6E99C2DDC', int64), int(z'0992CCAF6A6FCA05', int64)]
      case (96); jump = [int(z'360FD5F2CF8D5D99', int64), int(z'9C6E6877736C46E3', int64)]
      case default; return
      end select
      acc = 0_int64
      do i = 1, 2
         do b = 0, 63
            if (btest(jump(i), b)) acc = ieor(acc, s)
            call xoro_next(s, dummy)
         end do
      end do
      s = acc
   end subroutine xoro_jump_state

   pure subroutine ranlux_jump_state(s, p)
      integer(int64), intent(inout) :: s(9) !! ranlux++ state advanced by a supported precomputed jump.
      integer, intent(in) :: p !! Jump exponent: 32, 64, 96, 128, or 192.
      integer(int64) :: jump(9), z(18), out(9)
      select case (p)
      case (32)
         jump = [int(z'225D447B6009CA81',int64),int(z'203038DD5C7C7AF4',int64),int(z'24431478EE6B2A34',int64), &
            int(z'1C5BA10CF5DE7FF4',int64),int(z'3DA7537AAB4D5C21',int64),int(z'B9BF984F07B9D6EA',int64), &
            int(z'58E772114B65B524',int64),int(z'0F535213BAABDAFA',int64),int(z'6F0CF98612E25C09',int64)]
      case (64)
         jump = [int(z'B51C1DFE8C39D3CE',int64),int(z'C5ADF3554A5FFF69',int64),int(z'0166E1F247D5ECCA',int64), &
            int(z'1E24936423C97857',int64),int(z'9FECD738B2DD0192',int64),int(z'10D8C6D28F09F108',int64), &
            int(z'02DE9777A3C7EEC7',int64),int(z'5D30DD7D398EB561',int64),int(z'A9EA143361405140',int64)]
      case (96)
         jump = [int(z'9F1C67142C84C502',int64),int(z'024D94E3C4B490E8',int64),int(z'E9D460859F0659B6',int64), &
            int(z'D697D9321E8373B1',int64),int(z'1164275F61142884',int64),int(z'D644D1BD1837C737',int64), &
            int(z'AD4191BCF0926C6B',int64),int(z'2624A1B9EF2C42C0',int64),int(z'F671BBCEE85222AB',int64)]
      case (128)
         jump = [int(z'E5E5397B5C20CB13',int64),int(z'0D2DC1DFB1338F79',int64),int(z'9C80C974601AD658',int64), &
            int(z'90FDE743F57480FA',int64),int(z'144FEA51E6AFE769',int64),int(z'C6B1C11018312892',int64), &
            int(z'C9E664CCABD5195C',int64),int(z'7640D9F7D294B533',int64),int(z'48DECE251DE941E4',int64)]
      case (192)
         jump = [int(z'CDD87E6F8FF89DCE',int64),int(z'CF580FE4528CECF6',int64),int(z'64EBB69A5327B45F',int64), &
            int(z'F254C1E441D1568A',int64),int(z'4E448EFC9F9F97B1',int64),int(z'1443A786DAA27CFB',int64), &
            int(z'A4F6B958BE0E7E54',int64),int(z'EBC258B2F86E436A',int64),int(z'E8331FCF219A53F9',int64)]
      case default
         return
      end select
      call mul9x9(z, jump, s)
      call ranlux_mod(out, z)
      s = out
   end subroutine ranlux_jump_state

   pure subroutine pcg_jump_state(s, p)
      integer(int64), intent(inout) :: s(4) !! PCG64 state advanced by exactly 2^p LCG steps.
      integer, intent(in) :: p !! Jump exponent in 0..127.
      integer(int64), parameter :: multiplier = int(z'DA942042E4DD58B5', int64)
      integer(int64) :: mlo, mhi, ilo, ihi, mp1lo, mp1hi, lo, hi
      integer :: i
      mlo = multiplier
      mhi = 0_int64
      ilo = s(3)
      ihi = s(4)
      do i = 1, p
         call add128(mlo, mhi, 1_int64, 0_int64, mp1lo, mp1hi)
         call mul128_mod(mp1lo, mp1hi, ilo, ihi, lo, hi)
         ilo = lo
         ihi = hi
         call mul128_mod(mlo, mhi, mlo, mhi, lo, hi)
         mlo = lo
         mhi = hi
      end do
      call mul128_mod(mlo, mhi, s(1), s(2), lo, hi)
      call add128(lo, hi, ilo, ihi, s(1), s(2))
   end subroutine pcg_jump_state

   pure subroutine advance_pcg(st, delta, info)
      type(engine_state_type), intent(inout) :: st !! RNG state; must select the PCG64 engine.
      integer(int64), intent(in) :: delta(2) !! Unsigned 128-bit advance as low/high 64-bit words.
      integer, intent(out) :: info !! Zero on success; nonzero if the selected engine is not PCG64.
      integer(int64), parameter :: multiplier = int(z'DA942042E4DD58B5', int64)
      integer(int64) :: mlo, mhi, ilo, ihi, dlo, dhi, malo, mahi, ialo, iahi
      integer(int64) :: mp1lo, mp1hi, lo, hi, state_lo, state_hi
      if (st%engine /= eng_pcg64) then
         info = 1
         return
      end if
      state_lo = st%state(1)
      state_hi = st%state(2)
      mlo = multiplier
      mhi = 0_int64
      ilo = st%state(3)
      ihi = st%state(4)
      dlo = delta(1)
      dhi = delta(2)
      malo = 1_int64
      mahi = 0_int64
      ialo = 0_int64
      iahi = 0_int64
      do while (dlo /= 0_int64 .or. dhi /= 0_int64)
         if (btest(dlo, 0)) then
            call mul128_mod(malo, mahi, mlo, mhi, lo, hi)
            malo = lo
            mahi = hi
            call mul128_mod(ialo, iahi, mlo, mhi, lo, hi)
            call add128(lo, hi, ilo, ihi, ialo, iahi)
         end if
         call add128(mlo, mhi, 1_int64, 0_int64, mp1lo, mp1hi)
         call mul128_mod(mp1lo, mp1hi, ilo, ihi, lo, hi)
         ilo = lo
         ihi = hi
         call mul128_mod(mlo, mhi, mlo, mhi, lo, hi)
         mlo = lo
         mhi = hi
         dlo = ior(shiftr(dlo, 1), shiftl(dhi, 63))
         dhi = shiftr(dhi, 1)
      end do
      call mul128_mod(malo, mahi, state_lo, state_hi, lo, hi)
      call add128(lo, hi, ialo, iahi, st%state(1), st%state(2))
      st%chacha_pos = 9
      info = 0
   end subroutine advance_pcg

   pure subroutine jump_engine(st, p, info)
      type(engine_state_type), intent(inout) :: st !! RNG state to advance by a supported power-of-two jump.
      integer, intent(in) :: p !! Jump exponent requested by the upstream-compatible API.
      integer, intent(out) :: info !! Zero on success; nonzero for unsupported engine/exponent combinations.
      integer :: lane
      info = 0
      select case (st%engine)
      case (eng_x256pp, eng_x256ss)
         if (.not. any(p == [32,64,96,128,192])) then
         info = 2
         return
         end if
         call xoshiro_jump_state(st%state(1:4), p)
      case (eng_x256ppsimd, eng_x256sssimd)
         if (.not. any(p == [32,64,96,128,192])) then
         info = 2
         return
         end if
         do lane = 1, 8
            call xoshiro_jump_state(st%simd(:, lane), p)
         end do
      case (eng_x128p)
         if (.not. any(p == [32,64,96])) then
         info = 2
         return
         end if
         call x128_jump_state(st%state(1:2), p)
      case (eng_xoro)
         if (.not. any(p == [32,64,96])) then
         info = 2
         return
         end if
         call xoro_jump_state(st%state(1:2), p)
      case (eng_pcg64)
         if (p < 0 .or. p > 127) then
         info = 2
         return
         end if
         call pcg_jump_state(st%state(1:4), p)
      case (eng_ranluxpp)
         if (.not. any(p == [32,64,96,128,192])) then
         info = 2
         return
         end if
         call ranlux_jump_state(st%state(1:9), p)
      case default
         info = 1
         return
      end select
      st%chacha_pos = 9
      if (st%engine == eng_ranluxpp) st%chacha_pos = 10
   end subroutine jump_engine

   pure subroutine set_engine_state(st, state_words, info)
      type(engine_state_type), intent(inout) :: st !! RNG state whose primary engine words are to be replaced.
      integer(int64), intent(in) :: state_words(:) !! Packed unsigned 64-bit state words; exact required length depends on engine.
      integer, intent(out) :: info !! Zero on success; nonzero if the state-word count is invalid.
      integer :: n
      n = engine_state_words(st%engine)
      if (size(state_words) /= n) then
         info = 1
         return
      end if
      st%state = 0_int64
      st%state(1:n) = state_words
      st%chacha_pos = 9
      st%simd_lane = 1
      if (st%engine == eng_x256ppsimd .or. st%engine == eng_x256sssimd) call setup_x256_streams(st)
      if (st%engine == eng_sfc64simd) call setup_sfc_streams(st)
      info = 0
   end subroutine set_engine_state

   pure subroutine set_pcg_inc(st, inc, info)
      type(engine_state_type), intent(inout) :: st !! RNG state; must select PCG64.
      integer(int64), intent(in) :: inc(2) !! Packed low/high 64-bit PCG increment words.
      integer, intent(out) :: info !! Zero on success; nonzero for a non-PCG engine.
      if (st%engine /= eng_pcg64) then
      info = 1
      return
      end if
      st%state(3:4) = inc
      info = 0
   end subroutine set_pcg_inc

   pure subroutine set_cwg_weyl(st, weyl, info)
      type(engine_state_type), intent(inout) :: st !! RNG state; must select CWG128.
      integer(int64), intent(in) :: weyl(2) !! Packed low/high 64-bit Weyl increment words.
      integer, intent(out) :: info !! Zero on success; nonzero for a non-CWG engine.
      integer(int64) :: a, b
      integer :: i
      if (st%engine /= eng_cwg128) then
      info = 1
      return
      end if
      st%state(1:2) = weyl
      st%state(1) = ior(st%state(1), 1_int64)
      do i = 1, 96
         call cwg_next_pair(st%state(1:8), a, b)
      end do
      st%chacha_pos = 9
      info = 0
   end subroutine set_cwg_weyl

   pure subroutine set_sfc_abc(st, abc, info)
      type(engine_state_type), intent(inout) :: st !! RNG state; must select scalar SFC64.
      integer(int64), intent(in) :: abc(3) !! Replacement SFC64 a/b/c words; counter is retained.
      integer, intent(out) :: info !! Zero on success; nonzero for a non-SFC64 engine.
      integer(int64) :: dummy
      integer :: i
      if (st%engine /= eng_sfc64) then
      info = 1
      return
      end if
      st%state(1:3) = abc
      do i = 1, 18
         call sfc_next(st%state(1:4), dummy)
      end do
      info = 0
   end subroutine set_sfc_abc

   pure subroutine set_chacha_nonce(st, nonce, info)
      type(engine_state_type), intent(inout) :: st !! RNG state; must select ChaCha20.
      integer(int64), intent(in) :: nonce(3) !! Three unsigned 32-bit nonce words in int64 containers.
      integer, intent(out) :: info !! Zero on success; nonzero for a non-ChaCha engine.
      integer(int64) :: sw(12)
      if (st%engine /= eng_chacha20) then
      info = 1
      return
      end if
      call get_state_u32(st%state(1:6), sw)
      sw(9:11) = lo32(nonce)
      call set_state_u32(st%state(1:6), sw)
      st%chacha_pos = 9
      info = 0
   end subroutine set_chacha_nonce

   pure subroutine set_philox_key(st, key, info)
      type(engine_state_type), intent(inout) :: st !! RNG state; must select Philox.
      integer(int64), intent(in) :: key(2) !! Two packed 64-bit Philox key words.
      integer, intent(out) :: info !! Zero on success; nonzero for a non-Philox engine.
      if (st%engine /= eng_philox) then
      info = 1
      return
      end if
      st%state(5:6) = key
      st%state(7) = 0_int64
      info = 0
   end subroutine set_philox_key

   pure subroutine set_squares_key(st, key, info)
      type(engine_state_type), intent(inout) :: st !! RNG state; must select Squares64.
      integer(int64), intent(in) :: key !! Packed 64-bit Squares64 key.
      integer, intent(out) :: info !! Zero on success; nonzero for a non-Squares engine.
      if (st%engine /= eng_squares) then
      info = 1
      return
      end if
      st%state(2) = key
      info = 0
   end subroutine set_squares_key
end module randompack_engines
