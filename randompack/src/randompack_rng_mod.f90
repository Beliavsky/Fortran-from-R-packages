module randompack_rng_mod
   use iso_fortran_env, only : int8, int64
   use randompack_kinds, only : dp
   use randompack_uint, only : lo32, hi32, umul_wide
   use randompack_ziggurat, only : ki_double, wi_double, fi_double, ke_double, we_double, fe_double
   use randompack_ziggurat, only : dnorm, dexp, dzig52, ziggurat_nor_r, ziggurat_nor_inv_r, ziggurat_exp_r
   use randompack_openlibm, only : openlibm_exp, openlibm_log, openlibm_log1p
   use randompack_engines, only : engine_state_type, engine_from_name, engine_name, engine_description
   use randompack_engines, only : seed_engine, randomize_engine, next_u64, jump_engine, advance_pcg
   use randompack_engines, only : set_engine_state, set_pcg_inc, set_cwg_weyl, set_sfc_abc
   use randompack_engines, only : set_chacha_nonce, set_philox_key, set_squares_key
   use ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use r_linalg, only : cholesky_factor
   implicit none
   private
   real(dp), parameter :: pi_dp = 3.141592653589793238462643383279502884197_dp
   integer(int64), parameter :: mask32 = int(z'00000000FFFFFFFF', int64)

   type, public :: randompack_snapshot
      type(engine_state_type) :: core
      logical :: bitexact = .false.
      logical :: full_mantissa = .false.
      logical :: have_spare = .false.
      real(dp) :: spare_normal = 0.0_dp
      integer(int64) :: draw_word = 0_int64
      integer :: draw_byte = 8
   end type randompack_snapshot

   type, public :: randompack_rng_type
      private
      type(engine_state_type) :: core
      logical :: bitexact = .false.
      logical :: full_mantissa = .false.
      logical :: have_spare = .false.
      real(dp) :: spare_normal = 0.0_dp
      integer(int64) :: draw_word = 0_int64
      integer :: draw_byte = 8
   contains
      procedure :: seed => rng_seed
      procedure :: randomize => rng_randomize
      procedure :: jump => rng_jump
      procedure :: advance => rng_advance
      procedure :: set_state => rng_set_state
      procedure :: pcg64_set_inc => rng_pcg64_set_inc
      procedure :: cwg128_set_weyl => rng_cwg128_set_weyl
      procedure :: sfc64_set_abc => rng_sfc64_set_abc
      procedure :: chacha_set_nonce => rng_chacha_set_nonce
      procedure :: philox_set_key => rng_philox_set_key
      procedure :: squares_set_key => rng_squares_set_key
      procedure :: duplicate => rng_duplicate
      procedure :: serialize => rng_serialize
      procedure :: deserialize => rng_deserialize
      procedure :: engine_name => rng_engine_name
      procedure :: unif => rng_unif
      procedure :: normal => rng_normal
      procedure :: skew_normal => rng_skew_normal
      procedure :: lognormal => rng_lognormal
      procedure :: gumbel => rng_gumbel
      procedure :: pareto => rng_pareto
      procedure :: exp => rng_exp
      procedure :: gamma => rng_gamma
      procedure :: chi2 => rng_chi2
      procedure :: beta => rng_beta
      procedure :: t => rng_t
      procedure :: f => rng_f
      procedure :: weibull => rng_weibull
      procedure :: mvn => rng_mvn
      procedure :: int => rng_int
      procedure :: perm => rng_perm
      procedure :: sample => rng_sample
      procedure :: raw => rng_raw
   end type randompack_rng_type

   public :: randompack_rng, randompack_engines

contains
   function randompack_rng(engine, seed, bitexact, full_mantissa) result(rng)
      character(len=*), intent(in), optional :: engine !! Engine identifier; defaults to `x256++simd`.
      integer, intent(in), optional :: seed !! Optional deterministic signed 32-bit-style seed.
      logical, intent(in), optional :: bitexact !! Compatibility flag; exact support is method- and path-dependent.
      logical, intent(in), optional :: full_mantissa !! Use 53 random mantissa bits for uniform doubles when true.
      type(randompack_rng_type) :: rng
      character(len=12) :: selected
      integer :: info
      selected = 'x256++simd'
      if (present(engine)) selected = trim(engine)
      rng%core%engine = engine_from_name(selected)
      if (rng%core%engine == 0) then
         rng%core%engine = engine_from_name('x256++simd')
      end if
      if (present(bitexact)) rng%bitexact = bitexact
      if (present(full_mantissa)) rng%full_mantissa = full_mantissa
      if (present(seed)) then
         call seed_engine(rng%core, seed, info=info)
      else
         call randomize_engine(rng%core)
      end if
      rng%draw_byte = 8
      rng%have_spare = .false.
   end function randompack_rng

   subroutine randompack_engines(names, descriptions)
      character(len=12), intent(out) :: names(14) !! Engine identifiers in the same table order as upstream `randompack_engines()`.
      character(len=72), intent(out) :: descriptions(14) !! Human-readable descriptions corresponding to `names`.
      integer :: i
      do i = 1, 14
         names(i) = engine_name(i)
         descriptions(i) = engine_description(i)
      end do
   end subroutine randompack_engines

   pure subroutine set_info(info, value)
      integer, intent(out), optional :: info !! Optional status code receiving zero on success or a method-specific nonzero code.
      integer, intent(in) :: value !! Status code to store when `info` is present.
      if (present(info)) info = value
   end subroutine set_info

   pure subroutine reset_draw_cache(self)
      class(randompack_rng_type), intent(inout) :: self !! RNG object whose buffered byte and normal caches are invalidated.
      self%draw_byte = 8
      self%have_spare = .false.
   end subroutine reset_draw_cache

   subroutine rng_seed(self, seed, spawn_key, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object to deterministically reinitialize.
      integer, intent(in) :: seed !! Signed seed; low 32 bits reproduce the upstream integer seed convention.
      integer(int64), intent(in), optional :: spawn_key(:) !! Optional spawn-key words; only low 32 bits are used.
      integer, intent(out), optional :: info !! Zero on success; nonzero if engine initialization fails.
      integer :: ierr
      if (present(spawn_key)) then
         call seed_engine(self%core, seed, spawn_key, ierr)
      else
         call seed_engine(self%core, seed, info=ierr)
      end if
      call reset_draw_cache(self)
      call set_info(info, ierr)
   end subroutine rng_seed

   subroutine rng_randomize(self, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object to initialize from portable clock/date entropy.
      integer, intent(out), optional :: info !! Zero on success.
      call randomize_engine(self%core)
      call reset_draw_cache(self)
      call set_info(info, 0)
   end subroutine rng_randomize

   pure subroutine rng_jump(self, p, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object to jump ahead by a supported power of two.
      integer, intent(in) :: p !! Jump exponent accepted by the selected engine.
      integer, intent(out), optional :: info !! Zero on success; nonzero for unsupported engine/exponent combinations.
      integer :: ierr
      call jump_engine(self%core, p, ierr)
      if (ierr == 0) call reset_draw_cache(self)
      call set_info(info, ierr)
   end subroutine rng_jump

   pure subroutine rng_advance(self, delta, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object; selected engine must be PCG64.
      integer(int64), intent(in) :: delta(:) !! One or two packed 64-bit words forming an unsigned 128-bit advance.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid size or a non-PCG engine.
      integer(int64) :: d(2)
      integer :: ierr
      if (size(delta) < 1 .or. size(delta) > 2) then
         call set_info(info, 2)
         return
      end if
      d = 0_int64
      d(1:size(delta)) = delta
      call advance_pcg(self%core, d, ierr)
      if (ierr == 0) call reset_draw_cache(self)
      call set_info(info, ierr)
   end subroutine rng_advance

   pure subroutine rng_set_state(self, state, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object whose packed primary engine state is replaced.
      integer(int64), intent(in) :: state(:) !! Packed 64-bit state words; exact length depends on engine.
      integer, intent(out), optional :: info !! Zero on success; nonzero if the state length is invalid.
      integer :: ierr
      call set_engine_state(self%core, state, ierr)
      if (ierr == 0) call reset_draw_cache(self)
      call set_info(info, ierr)
   end subroutine rng_set_state

   pure subroutine rng_pcg64_set_inc(self, inc, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object; selected engine must be PCG64.
      integer(int64), intent(in) :: inc(:) !! One or two packed 64-bit words forming the PCG increment.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid size or engine.
      integer(int64) :: packed(2)
      integer :: ierr
      if (size(inc) < 1 .or. size(inc) > 2) then
      call set_info(info, 2)
      return
      end if
      packed = 0_int64
      packed(1:size(inc)) = inc
      call set_pcg_inc(self%core, packed, ierr)
      if (ierr == 0) call reset_draw_cache(self)
      call set_info(info, ierr)
   end subroutine rng_pcg64_set_inc

   pure subroutine rng_cwg128_set_weyl(self, weyl, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object; selected engine must be CWG128.
      integer(int64), intent(in) :: weyl(:) !! One or two packed 64-bit words forming the Weyl increment.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid size or engine.
      integer(int64) :: packed(2)
      integer :: ierr
      if (size(weyl) < 1 .or. size(weyl) > 2) then
      call set_info(info, 2)
      return
      end if
      packed = 0_int64
      packed(1:size(weyl)) = weyl
      call set_cwg_weyl(self%core, packed, ierr)
      if (ierr == 0) call reset_draw_cache(self)
      call set_info(info, ierr)
   end subroutine rng_cwg128_set_weyl

   pure subroutine rng_sfc64_set_abc(self, abc, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object; selected engine must be scalar SFC64.
      integer(int64), intent(in) :: abc(:) !! One to three packed 64-bit words replacing SFC64 a/b/c state.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid size or engine.
      integer(int64) :: packed(3)
      integer :: ierr
      if (size(abc) < 1 .or. size(abc) > 3) then
      call set_info(info, 2)
      return
      end if
      packed = 0_int64
      packed(1:size(abc)) = abc
      call set_sfc_abc(self%core, packed, ierr)
      if (ierr == 0) call reset_draw_cache(self)
      call set_info(info, ierr)
   end subroutine rng_sfc64_set_abc

   pure subroutine rng_chacha_set_nonce(self, nonce, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object; selected engine must be ChaCha20.
      integer(int64), intent(in) :: nonce(:) !! One to three unsigned 32-bit nonce words held in int64 containers.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid size or engine.
      integer(int64) :: packed(3)
      integer :: ierr
      if (size(nonce) < 1 .or. size(nonce) > 3) then
      call set_info(info, 2)
      return
      end if
      packed = 0_int64
      packed(1:size(nonce)) = nonce
      call set_chacha_nonce(self%core, packed, ierr)
      if (ierr == 0) call reset_draw_cache(self)
      call set_info(info, ierr)
   end subroutine rng_chacha_set_nonce

   pure subroutine rng_philox_set_key(self, key, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object; selected engine must be Philox.
      integer(int64), intent(in) :: key(:) !! One or two packed 64-bit Philox key words.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid size or engine.
      integer(int64) :: packed(2)
      integer :: ierr
      if (size(key) < 1 .or. size(key) > 2) then
      call set_info(info, 2)
      return
      end if
      packed = 0_int64
      packed(1:size(key)) = key
      call set_philox_key(self%core, packed, ierr)
      if (ierr == 0) call reset_draw_cache(self)
      call set_info(info, ierr)
   end subroutine rng_philox_set_key

   pure subroutine rng_squares_set_key(self, key, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object; selected engine must be Squares64.
      integer(int64), intent(in) :: key !! Packed 64-bit Squares64 key.
      integer, intent(out), optional :: info !! Zero on success; nonzero for a non-Squares engine.
      integer :: ierr
      call set_squares_key(self%core, key, ierr)
      if (ierr == 0) call reset_draw_cache(self)
      call set_info(info, ierr)
   end subroutine rng_squares_set_key

   pure function rng_duplicate(self) result(copy)
      class(randompack_rng_type), intent(in) :: self !! RNG object to copy including all cached stream state.
      type(randompack_rng_type) :: copy
      copy%core = self%core
      copy%bitexact = self%bitexact
      copy%full_mantissa = self%full_mantissa
      copy%have_spare = self%have_spare
      copy%spare_normal = self%spare_normal
      copy%draw_word = self%draw_word
      copy%draw_byte = self%draw_byte
   end function rng_duplicate

   pure function rng_serialize(self) result(snapshot)
      class(randompack_rng_type), intent(in) :: self !! RNG object whose complete Fortran state is captured.
      type(randompack_snapshot) :: snapshot
      snapshot%core = self%core
      snapshot%bitexact = self%bitexact
      snapshot%full_mantissa = self%full_mantissa
      snapshot%have_spare = self%have_spare
      snapshot%spare_normal = self%spare_normal
      snapshot%draw_word = self%draw_word
      snapshot%draw_byte = self%draw_byte
   end function rng_serialize

   pure subroutine rng_deserialize(self, snapshot)
      class(randompack_rng_type), intent(inout) :: self !! RNG object receiving a previously captured Fortran snapshot.
      type(randompack_snapshot), intent(in) :: snapshot !! Snapshot produced by `serialize` in this Fortran translation.
      self%core = snapshot%core
      self%bitexact = snapshot%bitexact
      self%full_mantissa = snapshot%full_mantissa
      self%have_spare = snapshot%have_spare
      self%spare_normal = snapshot%spare_normal
      self%draw_word = snapshot%draw_word
      self%draw_byte = snapshot%draw_byte
   end subroutine rng_deserialize

   pure function rng_engine_name(self) result(name)
      class(randompack_rng_type), intent(in) :: self !! RNG object whose selected engine identifier is requested.
      character(len=12) :: name
      name = engine_name(self%core%engine)
   end function rng_engine_name

   pure subroutine align_draw(self, bytes)
      class(randompack_rng_type), intent(inout) :: self !! RNG object whose subword buffer is aligned forward.
      integer, intent(in) :: bytes !! Requested alignment in bytes: 2, 4, or 8.
      integer :: pos
      if (self%draw_byte >= 8) return
      pos = ((self%draw_byte + bytes - 1) / bytes) * bytes
      if (pos >= 8) then
         self%draw_byte = 8
      else
         self%draw_byte = pos
      end if
   end subroutine align_draw

   pure subroutine ensure_draw_word(self)
      class(randompack_rng_type), intent(inout) :: self !! RNG object whose byte buffer is filled if exhausted.
      if (self%draw_byte >= 8) then
         call next_u64(self%core, self%draw_word)
         self%draw_byte = 0
      end if
   end subroutine ensure_draw_word

   pure subroutine draw_u64(self, word)
      class(randompack_rng_type), intent(inout) :: self !! RNG object advanced to the next 64-bit-aligned output word.
      integer(int64), intent(out) :: word !! Next raw 64-bit word after alignment.
      call align_draw(self, 8)
      call next_u64(self%core, word)
      self%draw_byte = 8
   end subroutine draw_u64

   pure subroutine draw_u16(self, word)
      class(randompack_rng_type), intent(inout) :: self !! RNG object advanced by one aligned 16-bit subword.
      integer(int64), intent(out) :: word !! Unsigned 16-bit result held in an int64 container.
      call align_draw(self, 2)
      call ensure_draw_word(self)
      word = iand(shiftr(self%draw_word, 8 * self%draw_byte), int(z'FFFF', int64))
      self%draw_byte = self%draw_byte + 2
   end subroutine draw_u16

   pure subroutine draw_u32(self, word)
      class(randompack_rng_type), intent(inout) :: self !! RNG object advanced by one aligned 32-bit subword.
      integer(int64), intent(out) :: word !! Unsigned 32-bit result held in an int64 container.
      call align_draw(self, 4)
      call ensure_draw_word(self)
      if (self%draw_byte == 0) then
         word = lo32(self%draw_word)
      else
         word = hi32(self%draw_word)
      end if
      self%draw_byte = self%draw_byte + 4
   end subroutine draw_u32

   pure subroutine draw_byte(self, byte)
      class(randompack_rng_type), intent(inout) :: self !! RNG object advanced by one raw byte.
      integer(int64), intent(out) :: byte !! Unsigned byte in the range 0..255.
      call ensure_draw_word(self)
      byte = iand(shiftr(self%draw_word, 8 * self%draw_byte), 255_int64)
      self%draw_byte = self%draw_byte + 1
   end subroutine draw_byte

   pure subroutine uniform_scalar(self, value)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying one raw 64-bit word.
      real(dp), intent(out) :: value !! Uniform variate on [0,1).
      integer(int64) :: word, k
      call draw_u64(self, word)
      if (self%full_mantissa) then
         k = shiftr(word, 11)
         value = real(k, dp) * 2.0_dp**(-53)
      else
         k = shiftr(word, 12)
         value = real(k, dp) * 2.0_dp**(-52)
      end if
   end subroutine uniform_scalar

   pure elemental function u64_to_unit(word) result(value)
      integer(int64), intent(in) :: word !! Unsigned 64-bit word whose upper 53 bits define a uniform variate.
      real(dp) :: value
      value = real(shiftr(word, 11), dp) * 2.0_dp**(-53)
   end function u64_to_unit

   pure subroutine norm_decode(word, idx, ki, rabs, sign, value)
      integer(int64), intent(in) :: word !! Raw Ziggurat proposal word.
      integer, intent(out) :: idx !! Zero-based Ziggurat strip index in 0..255.
      integer(int64), intent(out) :: ki !! Fast-acceptance threshold for the selected strip.
      integer(int64), intent(out) :: rabs !! Unsigned 52-bit proposal magnitude.
      logical, intent(out) :: sign !! True when the proposal is negative.
      real(dp), intent(out) :: value !! Decoded normal proposal.
      idx = int(iand(word, 255_int64))
      ki = ki_double(idx + 1)
      rabs = iand(shiftr(word, 9), int(z'000FFFFFFFFFFFFF', int64))
      sign = btest(word, 8)
      value = real(rabs, dp) * wi_double(idx + 1)
      if (sign) value = -value
   end subroutine norm_decode

   pure subroutine norm_tail_draw(self, sign, y0, value)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying tail-rejection uniforms.
      logical, intent(in) :: sign !! Sign of the original tail proposal.
      real(dp), intent(in) :: y0 !! Initial independent uniform variate on [0,1).
      real(dp), intent(out) :: value !! Accepted normal-tail variate.
      real(dp) :: x, xx, y, yy
      y = y0
      do
         call uniform_scalar(self, x)
         xx = -ziggurat_nor_inv_r * openlibm_log1p(-x)
         yy = -openlibm_log1p(-y)
         if (yy + yy > xx * xx) then
            value = ziggurat_nor_r + xx
            if (sign) value = -value
            return
         end if
         call uniform_scalar(self, y)
      end do
   end subroutine norm_tail_draw

   pure subroutine norm_from_word(self, word, value)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying rare-path Ziggurat words.
      integer(int64), intent(in) :: word !! Initial raw proposal word.
      real(dp), intent(out) :: value !! Accepted standard-normal variate.
      integer(int64), parameter :: two52 = 4503599627370496_int64
      integer(int64) :: current, high, ki, low, lval, rabs, rval, yword
      integer :: idx
      logical :: sign
      real(dp) :: u, x
      current = word
      do
         call norm_decode(current, idx, ki, rabs, sign, x)
         if (rabs < ki) then
            value = x
            return
         end if
         call draw_u64(self, yword)
         if (idx == 0) then
            u = u64_to_unit(yword)
            call norm_tail_draw(self, sign, u, value)
            return
         end if
         lval = two52 - ki
         rval = two52 - rabs
         call umul_wide(yword, lval, high, low)
         if (idx > 52) then
            if (high <= rval) then
               if (high + dnorm(idx + 1) < rval) then
                  value = x
                  return
               end if
            else
               call draw_u64(self, current)
               cycle
            end if
         else if (idx < 52) then
            if (high < rval) then
               value = x
               return
            end if
            if (high > rval + dnorm(idx + 1)) then
               call draw_u64(self, current)
               cycle
            end if
         else
            if (high > rval + dzig52) then
               call draw_u64(self, current)
               cycle
            end if
            if (high + dnorm(idx + 1) < rval) then
               value = x
               return
            end if
         end if
         u = u64_to_unit(yword)
         if (fi_double(idx + 1) + u * (fi_double(idx) - fi_double(idx + 1)) < openlibm_exp(-0.5_dp * x * x)) then
            value = x
            return
         end if
         call draw_u64(self, current)
      end do
   end subroutine norm_from_word

   pure subroutine normal_scalar(self, value)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying one Ziggurat normal proposal stream.
      real(dp), intent(out) :: value !! Standard-normal variate.
      integer(int64) :: word
      call draw_u64(self, word)
      call norm_from_word(self, word, value)
   end subroutine normal_scalar

   pure subroutine exp_tail_draw(self, value)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying the exponential-tail uniform.
      real(dp), intent(out) :: value !! Accepted unit-exponential tail variate.
      real(dp) :: u
      call uniform_scalar(self, u)
      value = ziggurat_exp_r - openlibm_log1p(-u)
   end subroutine exp_tail_draw

   pure subroutine exp_from_word(self, word, value)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying rare-path Ziggurat words.
      integer(int64), intent(in) :: word !! Initial raw exponential proposal word.
      real(dp), intent(out) :: value !! Accepted unit-exponential variate.
      integer(int64), parameter :: two53 = 9007199254740992_int64
      integer(int64) :: current, high, ke, low, lval, rval, rmag, yword
      integer :: idx
      real(dp) :: u, x
      current = word
      do
         idx = int(iand(current, 255_int64))
         rmag = iand(shiftr(current, 8), int(z'001FFFFFFFFFFFFF', int64))
         ke = ke_double(idx + 1)
         x = real(rmag, dp) * we_double(idx + 1)
         if (rmag < ke) then
            value = x
            return
         end if
         if (idx == 0) then
            call exp_tail_draw(self, value)
            return
         end if
         call draw_u64(self, yword)
         lval = two53 - ke
         rval = two53 - rmag
         call umul_wide(yword, lval, high, low)
         if (high <= rval) then
            if (high + dexp(idx + 1) < rval) then
               value = x
               return
            end if
            u = u64_to_unit(yword)
            if (fe_double(idx + 1) + u * (fe_double(idx) - fe_double(idx + 1)) < openlibm_exp(-x)) then
               value = x
               return
            end if
         end if
         call draw_u64(self, current)
      end do
   end subroutine exp_from_word

   pure subroutine u01_53_scalar(self, value)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying one raw word for a 53-bit [0,1) draw.
      real(dp), intent(out) :: value !! Uniform variate formed from the upper 53 bits, independent of `full_mantissa`.
      integer(int64) :: word
      call draw_u64(self, word)
      value = u64_to_unit(word)
   end subroutine u01_53_scalar

   pure recursive subroutine gamma_scalar(self, shape, scale, value, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying gamma acceptance-rejection draws.
      real(dp), intent(in) :: shape !! Gamma shape parameter; must be strictly positive.
      real(dp), intent(in) :: scale !! Gamma scale parameter; must be strictly positive.
      real(dp), intent(out) :: value !! Generated gamma variate on success.
      integer, intent(out) :: info !! Zero on success; one for invalid shape/scale.
      real(dp) :: c, d, u, v, x
      if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         value = 0.0_dp
         info = 1
         return
      end if
      if (shape < 1.0_dp) then
         call gamma_scalar(self, shape + 1.0_dp, scale, value, info)
         if (info /= 0) return
         call uniform_scalar(self, u)
         if (self%bitexact) then
            value = value * openlibm_exp(openlibm_log(u) / shape)
         else
            value = value * exp(log(u) / shape)
         end if
         return
      end if
      d = shape - 1.0_dp / 3.0_dp
      c = 1.0_dp / sqrt(9.0_dp * d)
      do
         call normal_scalar(self, x)
         v = 1.0_dp + c * x
         if (v <= 0.0_dp) cycle
         v = v * v * v
         call u01_53_scalar(self, u)
         if (u < 1.0_dp - 0.0331_dp * x**4) exit
         if (openlibm_log(u) < 0.5_dp * x * x + d * (1.0_dp - v + openlibm_log(v))) exit
      end do
      value = scale * d * v
      info = 0
   end subroutine gamma_scalar

   pure subroutine gamma_shape1_vector(self, x, shape, scale)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying the upstream vector gamma proposal stream.
      real(dp), intent(out) :: x(:) !! Output gamma variates for a shape parameter at least one.
      real(dp), intent(in) :: shape !! Gamma shape parameter; must be at least one.
      real(dp), intent(in) :: scale !! Positive gamma scale parameter.
      real(dp) :: c, d, dscale, logu, logv3, u, v, v3, z, z2, z4
      integer :: i, ierr
      d = shape - 1.0_dp / 3.0_dp
      dscale = d * scale
      c = 1.0_dp / sqrt(9.0_dp * d)
      call rng_normal(self, x, info=ierr)
      do i = 1, size(x)
         do
            z = x(i)
            v = 1.0_dp + c * z
            if (v > 0.0_dp) then
               z2 = z * z
               z4 = z2 * z2
               v3 = v * v * v
               call u01_53_scalar(self, u)
               if (u < 1.0_dp - 0.0331_dp * z4) exit
               logu = openlibm_log(u)
               logv3 = openlibm_log(v3)
               if (logu < 0.5_dp * z2 + d * (1.0_dp - v3 + logv3)) exit
            end if
            call normal_scalar(self, x(i))
         end do
         x(i) = dscale * v3
      end do
   end subroutine gamma_shape1_vector

   pure subroutine rng_unif(self, x, a, b, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying uniform variates.
      real(dp), intent(out) :: x(:) !! Output variates, one per element, on [a,b).
      real(dp), intent(in), optional :: a !! Lower endpoint; defaults to zero and must be smaller than `b`.
      real(dp), intent(in), optional :: b !! Upper endpoint; defaults to one and must exceed `a`.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid endpoints.
      real(dp) :: lower, upper, u
      integer :: i
      lower = 0.0_dp
      upper = 1.0_dp
      if (present(a)) lower = a
      if (present(b)) upper = b
      if (.not. (lower < upper)) then
      x = 0.0_dp
      call set_info(info, 1)
      return
      end if
      do i = 1, size(x)
         call uniform_scalar(self, u)
         x(i) = lower + (upper - lower) * u
      end do
      call set_info(info, 0)
   end subroutine rng_unif

   pure subroutine rng_normal(self, x, mu, sigma, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying upstream-order Ziggurat normal variates.
      real(dp), intent(out) :: x(:) !! Output normal variates.
      real(dp), intent(in), optional :: mu !! Mean; defaults to zero.
      real(dp), intent(in), optional :: sigma !! Standard deviation; defaults to one and must be positive.
      integer, intent(out), optional :: info !! Zero on success; nonzero for nonpositive sigma.
      integer(int64), allocatable :: raw(:)
      real(dp) :: center, sd, z
      integer :: i
      center = 0.0_dp
      sd = 1.0_dp
      if (present(mu)) center = mu
      if (present(sigma)) sd = sigma
      if (sd <= 0.0_dp) then
         x = 0.0_dp
         call set_info(info, 1)
         return
      end if
      allocate(raw(size(x)))
      do i = 1, size(x)
         call draw_u64(self, raw(i))
      end do
      do i = size(x), 1, -1
         call norm_from_word(self, raw(i), z)
         x(i) = center + sd * z
      end do
      call set_info(info, 0)
   end subroutine rng_normal

   pure subroutine rng_skew_normal(self, x, mu, sigma, alpha, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying skew-normal component variates.
      real(dp), intent(out) :: x(:) !! Output skew-normal variates.
      real(dp), intent(in), optional :: mu !! Location parameter; defaults to zero.
      real(dp), intent(in), optional :: sigma !! Positive scale parameter; defaults to one.
      real(dp), intent(in) :: alpha !! Skew-normal shape parameter.
      integer, intent(out), optional :: info !! Zero on success; nonzero for nonpositive scale.
      real(dp) :: center, sd, denom, delta, scale0, u(128), v
      integer :: i, j, n, ierr
      center = 0.0_dp
      sd = 1.0_dp
      if (present(mu)) center = mu
      if (present(sigma)) sd = sigma
      if (sd <= 0.0_dp) then
         x = 0.0_dp
         call set_info(info, 1)
         return
      end if
      denom = sqrt(1.0_dp + alpha * alpha)
      delta = alpha / denom
      scale0 = 1.0_dp / denom
      call rng_normal(self, x, info=ierr)
      if (ierr /= 0) then
         call set_info(info, ierr)
         return
      end if
      i = 1
      do while (i <= size(x))
         n = min(128, size(x) - i + 1)
         call rng_normal(self, u(1:n), info=ierr)
         if (ierr /= 0) then
            call set_info(info, ierr)
            return
         end if
         do j = 1, n
            v = delta * abs(u(j)) + scale0 * x(i + j - 1)
            x(i + j - 1) = center + sd * v
         end do
         i = i + n
      end do
      call set_info(info, 0)
   end subroutine rng_skew_normal

   pure subroutine rng_lognormal(self, x, mu, sigma, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying lognormal variates.
      real(dp), intent(out) :: x(:) !! Output lognormal variates.
      real(dp), intent(in), optional :: mu !! Mean of the underlying normal; defaults to zero.
      real(dp), intent(in), optional :: sigma !! Positive standard deviation of the underlying normal; defaults to one.
      integer, intent(out), optional :: info !! Zero on success; nonzero for nonpositive sigma.
      real(dp) :: center, sd
      integer :: i, ierr
      center = 0.0_dp
      sd = 1.0_dp
      if (present(mu)) center = mu
      if (present(sigma)) sd = sigma
      if (sd <= 0.0_dp) then
         x = 0.0_dp
         call set_info(info, 1)
         return
      end if
      call rng_normal(self, x, info=ierr)
      if (ierr /= 0) then
         call set_info(info, ierr)
         return
      end if
      if (center /= 0.0_dp .or. sd /= 1.0_dp) then
         do i = 1, size(x)
            x(i) = center + sd * x(i)
         end do
      end if
      if (self%bitexact) then
         do i = 1, size(x)
            x(i) = openlibm_exp(x(i))
         end do
      else
         do i = 1, size(x)
            x(i) = exp(x(i))
         end do
      end if
      call set_info(info, 0)
   end subroutine rng_lognormal

   pure subroutine rng_gumbel(self, x, mu, beta, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying inverse-transform Gumbel variates.
      real(dp), intent(out) :: x(:) !! Output Gumbel variates.
      real(dp), intent(in), optional :: mu !! Location parameter; defaults to zero.
      real(dp), intent(in), optional :: beta !! Positive scale parameter; defaults to one.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid scale.
      real(dp) :: center, scale
      integer :: i
      center = 0.0_dp
      scale = 1.0_dp
      if (present(mu)) center = mu
      if (present(beta)) scale = beta
      if (scale <= 0.0_dp) then
         x = 0.0_dp
         call set_info(info, 1)
         return
      end if
      do i = 1, size(x)
         call uniform_scalar(self, x(i))
      end do
      if (self%bitexact) then
         do i = 1, size(x)
            x(i) = openlibm_log(x(i))
            x(i) = -x(i)
            x(i) = openlibm_log(x(i))
            x(i) = center + (-scale) * x(i)
         end do
      else
         do i = 1, size(x)
            x(i) = log(x(i))
            x(i) = -x(i)
            x(i) = log(x(i))
            x(i) = center + (-scale) * x(i)
         end do
      end if
      call set_info(info, 0)
   end subroutine rng_gumbel

   pure subroutine rng_pareto(self, x, xm, alpha, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying Pareto variates.
      real(dp), intent(out) :: x(:) !! Output Pareto variates with lower endpoint `xm`.
      real(dp), intent(in) :: xm !! Positive Pareto minimum value.
      real(dp), intent(in) :: alpha !! Positive Pareto shape parameter.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid parameters.
      integer :: i, ierr
      if (xm <= 0.0_dp .or. alpha <= 0.0_dp) then
         x = 0.0_dp
         call set_info(info, 1)
         return
      end if
      call rng_exp(self, x, info=ierr)
      if (ierr /= 0) then
         call set_info(info, ierr)
         return
      end if
      if (self%bitexact) then
         do i = 1, size(x)
            x(i) = x(i) / alpha
            x(i) = openlibm_exp(x(i))
            x(i) = xm * x(i)
         end do
      else
         do i = 1, size(x)
            x(i) = x(i) / alpha
            x(i) = exp(x(i))
            x(i) = xm * x(i)
         end do
      end if
      call set_info(info, 0)
   end subroutine rng_pareto

   pure subroutine rng_exp(self, x, scale, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying upstream-order Ziggurat exponential variates.
      real(dp), intent(out) :: x(:) !! Output exponential variates.
      real(dp), intent(in), optional :: scale !! Positive exponential scale; defaults to one.
      integer, intent(out), optional :: info !! Zero on success; nonzero for nonpositive scale.
      integer(int64), allocatable :: raw(:)
      real(dp) :: s, value
      integer :: i
      s = 1.0_dp
      if (present(scale)) s = scale
      if (s <= 0.0_dp) then
         x = 0.0_dp
         call set_info(info, 1)
         return
      end if
      allocate(raw(size(x)))
      do i = 1, size(x)
         call draw_u64(self, raw(i))
      end do
      do i = size(x), 1, -1
         call exp_from_word(self, raw(i), value)
         x(i) = s * value
      end do
      call set_info(info, 0)
   end subroutine rng_exp

   pure subroutine rng_gamma(self, x, shape, scale, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying upstream-order gamma variates.
      real(dp), intent(out) :: x(:) !! Output gamma variates.
      real(dp), intent(in) :: shape !! Positive gamma shape parameter.
      real(dp), intent(in), optional :: scale !! Positive gamma scale parameter; defaults to one.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid parameters.
      real(dp), allocatable :: boost(:)
      real(dp) :: inv_shape, s
      integer :: i
      s = 1.0_dp
      if (present(scale)) s = scale
      if (shape <= 0.0_dp .or. s <= 0.0_dp) then
         x = 0.0_dp
         call set_info(info, 1)
         return
      end if
      if (shape >= 1.0_dp) then
         call gamma_shape1_vector(self, x, shape, s)
      else
         call gamma_shape1_vector(self, x, shape + 1.0_dp, s)
         allocate(boost(size(x)))
         do i = 1, size(boost)
            call uniform_scalar(self, boost(i))
         end do
         inv_shape = 1.0_dp / shape
         if (self%bitexact) then
            do i = 1, size(boost)
               boost(i) = openlibm_log(boost(i))
               boost(i) = boost(i) * inv_shape
               boost(i) = openlibm_exp(boost(i))
            end do
         else
            do i = 1, size(boost)
               boost(i) = log(boost(i))
               boost(i) = boost(i) * inv_shape
               boost(i) = exp(boost(i))
            end do
         end if
         x = x * boost
      end if
      call set_info(info, 0)
   end subroutine rng_gamma

   pure subroutine rng_chi2(self, x, nu, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying chi-square variates.
      real(dp), intent(out) :: x(:) !! Output chi-square variates.
      real(dp), intent(in) :: nu !! Positive degrees of freedom.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid degrees of freedom.
      integer :: ierr
      if (nu <= 0.0_dp) then
      x = 0.0_dp
      call set_info(info, 1)
      return
      end if
      call rng_gamma(self, x, 0.5_dp * nu, 2.0_dp, ierr)
      call set_info(info, ierr)
   end subroutine rng_chi2

   pure subroutine rng_beta(self, x, a, b, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying beta variates.
      real(dp), intent(out) :: x(:) !! Output beta variates constrained to the open interval (0,1).
      real(dp), intent(in) :: a !! Positive first beta shape parameter.
      real(dp), intent(in) :: b !! Positive second beta shape parameter.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid shapes.
      real(dp) :: gb(128), ga, gbj, hi, y
      integer :: i, j, n, ierr
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         x = 0.0_dp
         call set_info(info, 1)
         return
      end if
      call rng_gamma(self, x, a, 1.0_dp, ierr)
      if (ierr /= 0) then
         call set_info(info, ierr)
         return
      end if
      hi = 1.0_dp - epsilon(1.0_dp)
      i = 1
      do while (i <= size(x))
         n = min(128, size(x) - i + 1)
         call rng_gamma(self, gb(1:n), b, 1.0_dp, ierr)
         if (ierr /= 0) then
            call set_info(info, ierr)
            return
         end if
         do j = 1, n
            ga = x(i + j - 1)
            gbj = gb(j)
            y = ga / (ga + gbj)
            if (y < tiny(1.0_dp)) y = tiny(1.0_dp)
            if (y > hi) y = hi
            x(i + j - 1) = y
         end do
         i = i + n
      end do
      call set_info(info, 0)
   end subroutine rng_beta

   pure subroutine rng_t(self, x, nu, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying Student-t variates.
      real(dp), intent(out) :: x(:) !! Output Student-t variates.
      real(dp), intent(in) :: nu !! Positive degrees of freedom.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid degrees of freedom.
      real(dp) :: u(128)
      integer :: i, j, n, ierr
      if (nu <= 0.0_dp) then
         x = 0.0_dp
         call set_info(info, 1)
         return
      end if
      call rng_normal(self, x, info=ierr)
      if (ierr /= 0) then
         call set_info(info, ierr)
         return
      end if
      i = 1
      do while (i <= size(x))
         n = min(128, size(x) - i + 1)
         call rng_gamma(self, u(1:n), 0.5_dp * nu, 2.0_dp, ierr)
         if (ierr /= 0) then
            call set_info(info, ierr)
            return
         end if
         do j = 1, n
            x(i + j - 1) = x(i + j - 1) / sqrt(u(j) / nu)
         end do
         i = i + n
      end do
      call set_info(info, 0)
   end subroutine rng_t

   pure subroutine rng_f(self, x, nu1, nu2, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying F-distribution variates.
      real(dp), intent(out) :: x(:) !! Output F variates.
      real(dp), intent(in) :: nu1 !! Positive numerator degrees of freedom.
      real(dp), intent(in) :: nu2 !! Positive denominator degrees of freedom.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid degrees of freedom.
      real(dp) :: x2(128)
      integer :: i, j, n, ierr
      if (nu1 <= 0.0_dp .or. nu2 <= 0.0_dp) then
         x = 0.0_dp
         call set_info(info, 1)
         return
      end if
      call rng_gamma(self, x, 0.5_dp * nu1, 1.0_dp, ierr)
      if (ierr /= 0) then
         call set_info(info, ierr)
         return
      end if
      i = 1
      do while (i <= size(x))
         n = min(128, size(x) - i + 1)
         call rng_gamma(self, x2(1:n), 0.5_dp * nu2, 1.0_dp, ierr)
         if (ierr /= 0) then
            call set_info(info, ierr)
            return
         end if
         do j = 1, n
            x(i + j - 1) = (x(i + j - 1) * nu2) / (x2(j) * nu1)
         end do
         i = i + n
      end do
      call set_info(info, 0)
   end subroutine rng_f

   pure subroutine rng_weibull(self, x, shape, scale, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying Weibull variates.
      real(dp), intent(out) :: x(:) !! Output Weibull variates.
      real(dp), intent(in) :: shape !! Positive Weibull shape parameter.
      real(dp), intent(in), optional :: scale !! Positive Weibull scale; defaults to one.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid parameters.
      real(dp) :: s, inv_shape
      integer :: i, ierr
      s = 1.0_dp
      if (present(scale)) s = scale
      if (shape <= 0.0_dp .or. s <= 0.0_dp) then
         x = 0.0_dp
         call set_info(info, 1)
         return
      end if
      call rng_exp(self, x, info=ierr)
      if (ierr /= 0) then
         call set_info(info, ierr)
         return
      end if
      if (shape == 1.0_dp) then
         if (s /= 1.0_dp) then
            do i = 1, size(x)
               x(i) = s * x(i)
            end do
         end if
         call set_info(info, 0)
         return
      end if
      if (shape == 2.0_dp) then
         do i = 1, size(x)
            x(i) = sqrt(x(i))
         end do
         if (s /= 1.0_dp) then
            do i = 1, size(x)
               x(i) = s * x(i)
            end do
         end if
         call set_info(info, 0)
         return
      end if
      inv_shape = 1.0_dp / shape
      if (self%bitexact) then
         do i = 1, size(x)
            x(i) = openlibm_log(x(i))
            x(i) = inv_shape * x(i)
            x(i) = openlibm_exp(x(i))
            if (s /= 1.0_dp) x(i) = s * x(i)
         end do
      else
         do i = 1, size(x)
            x(i) = log(x(i))
            x(i) = inv_shape * x(i)
            x(i) = exp(x(i))
            if (s /= 1.0_dp) x(i) = s * x(i)
         end do
      end if
      call set_info(info, 0)
   end subroutine rng_weibull

   pure subroutine pivoted_cholesky_lower(a, factor, rank, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric covariance matrix to factor using complete diagonal pivoting.
      real(dp), allocatable, intent(out) :: factor(:, :) !! Unpivoted lower factor with zero trailing columns after rank.
      integer, intent(out) :: rank !! Numerical rank determined using the upstream fixed tolerance `1.0e-14`.
      integer, intent(out) :: info !! Zero for full rank; one when the positive-semidefinite factor stops early.
      real(dp), allocatable :: work_a(:, :), work(:), residual(:)
      real(dp) :: ajj, dstop, s, temp
      integer, allocatable :: piv(:)
      integer :: i, j, k, n, pvt, itemp

      n = size(a, 1)
      allocate(factor(n, n))
      factor = 0.0_dp
      rank = 0
      info = 1
      if (n == 0 .or. size(a, 2) /= n) return

      allocate(work_a(n, n), work(n), residual(n), piv(n))
      work_a = 0.0_dp
      do j = 1, n
         do i = j, n
            work_a(i, j) = a(i, j)
         end do
      end do
      do i = 1, n
         piv(i) = i
      end do

      pvt = 1
      ajj = work_a(1, 1)
      do i = 2, n
         if (work_a(i, i) > ajj) then
            pvt = i
            ajj = work_a(i, i)
         end if
      end do
      if (ajj <= 0.0_dp .or. ieee_is_nan(ajj)) return

      dstop = 1.0e-14_dp
      work = 0.0_dp
      info = 0
      rank = n
      do j = 1, n
         do i = j, n
            if (j > 1) work(i) = work(i) + work_a(i, j - 1) * work_a(i, j - 1)
            residual(i) = work_a(i, i) - work(i)
         end do
         if (j > 1) then
            pvt = j
            ajj = residual(j)
            do i = j + 1, n
               if (residual(i) > ajj) then
                  pvt = i
                  ajj = residual(i)
               end if
            end do
            if (ajj <= dstop .or. ieee_is_nan(ajj)) then
               work_a(j, j) = ajj
               rank = j - 1
               info = 1
               exit
            end if
         end if

         if (j /= pvt) then
            work_a(pvt, pvt) = work_a(j, j)
            do k = 1, j - 1
               temp = work_a(j, k)
               work_a(j, k) = work_a(pvt, k)
               work_a(pvt, k) = temp
            end do
            do i = pvt + 1, n
               temp = work_a(i, j)
               work_a(i, j) = work_a(i, pvt)
               work_a(i, pvt) = temp
            end do
            do i = j + 1, pvt - 1
               temp = work_a(i, j)
               work_a(i, j) = work_a(pvt, i)
               work_a(pvt, i) = temp
            end do
            temp = work(j)
            work(j) = work(pvt)
            work(pvt) = temp
            itemp = piv(pvt)
            piv(pvt) = piv(j)
            piv(j) = itemp
         end if

         ajj = sqrt(ajj)
         work_a(j, j) = ajj
         if (j < n) then
            do i = j + 1, n
               s = work_a(i, j)
               do k = 1, j - 1
                  s = s - work_a(i, k) * work_a(j, k)
               end do
               work_a(i, j) = s / ajj
            end do
         end if
      end do

      factor = 0.0_dp
      do j = 1, n
         do k = 1, rank
            factor(piv(j), k) = work_a(j, k)
         end do
      end do
   end subroutine pivoted_cholesky_lower

   subroutine rng_mvn(self, x, sigma, mu, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying independent standard-normal components.
      real(dp), intent(out) :: x(:, :) !! Output matrix shaped `(n,d)`, matching the R method's row-per-draw layout.
      real(dp), intent(in) :: sigma(:, :) !! Finite symmetric positive-semidefinite covariance matrix shaped `(d,d)`.
      real(dp), intent(in), optional :: mu(:) !! Optional finite mean vector of length `d`; defaults to zero.
      integer, intent(out), optional :: info !! Zero on success; nonzero for shape, finiteness, symmetry, or factorization failure.
      real(dp), allocatable :: factor(:, :), center(:), zflat(:)
      real(dp) :: s
      integer :: d, i, ierr, j, k, n, rank

      n = size(x, 1)
      d = size(sigma, 1)
      if (d <= 0 .or. size(sigma, 2) /= d .or. size(x, 2) /= d) then
         x = 0.0_dp
         call set_info(info, 1)
         return
      end if
      if (.not. all(ieee_is_finite(sigma))) then
         x = 0.0_dp
         call set_info(info, 2)
         return
      end if
      if (.not. all(sigma == transpose(sigma))) then
         x = 0.0_dp
         call set_info(info, 3)
         return
      end if

      allocate(center(d))
      center = 0.0_dp
      if (present(mu)) then
         if (size(mu) /= d) then
            x = 0.0_dp
            call set_info(info, 4)
            return
         end if
         if (.not. all(ieee_is_finite(mu))) then
            x = 0.0_dp
            call set_info(info, 5)
            return
         end if
         center = mu
      end if

      call cholesky_factor(sigma, factor, ierr)
      if (ierr == 0) then
         do j = 1, d
            call rng_normal(self, x(:, j), info=ierr)
            if (ierr /= 0) then
               x = 0.0_dp
               call set_info(info, 6)
               return
            end if
         end do
         do i = 1, n
            do j = d, 1, -1
               s = factor(j, j) * x(i, j)
               do k = 1, j - 1
                  s = s + factor(j, k) * x(i, k)
               end do
               x(i, j) = s
            end do
         end do
      else
         call pivoted_cholesky_lower(sigma, factor, rank, ierr)
         if (rank > 0) then
            allocate(zflat(n * rank))
            call rng_normal(self, zflat, info=ierr)
            if (ierr /= 0) then
               x = 0.0_dp
               call set_info(info, 7)
               return
            end if
            do j = 1, d
               do i = 1, n
                  s = 0.0_dp
                  do k = 1, rank
                     s = s + factor(j, k) * zflat(i + (k - 1) * n)
                  end do
                  x(i, j) = s
               end do
            end do
         else
            x = 0.0_dp
         end if
      end if

      do i = 1, n
         do j = 1, d
            x(i, j) = x(i, j) + center(j)
         end do
      end do
      call set_info(info, 0)
   end subroutine rng_mvn

   pure subroutine bounded_index(self, bound, value, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying upstream-compatible bounded integers.
      integer, intent(in) :: bound !! Number of equally likely outcomes; must be positive and fit default integer range.
      integer, intent(out) :: value !! Uniform integer in `0..bound-1` on success.
      integer, intent(out) :: info !! Zero on success; one for nonpositive bound.
      integer(int64) :: r, prod, low, threshold
      if (bound <= 0) then
         value = 0
         info = 1
         return
      end if
      if (bound <= 1000) then
         threshold = modulo(65536_int64, int(bound, int64))
         do
            call draw_u16(self, r)
            prod = r * int(bound, int64)
            low = iand(prod, int(z'FFFF', int64))
            if (low >= threshold) exit
         end do
         value = int(shiftr(prod, 16))
      else
         threshold = modulo(4294967296_int64, int(bound, int64))
         do
            call draw_u32(self, r)
            prod = r * int(bound, int64)
            low = iand(prod, mask32)
            if (low >= threshold) exit
         end do
         value = int(shiftr(prod, 32))
      end if
      info = 0
   end subroutine bounded_index

   pure subroutine rng_int(self, x, min_value, max_value, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying discrete-uniform integer draws.
      integer, intent(out) :: x(:) !! Output integers on the inclusive interval `[min_value,max_value]`.
      integer, intent(in) :: min_value !! Inclusive lower integer endpoint.
      integer, intent(in) :: max_value !! Inclusive upper integer endpoint.
      integer, intent(out), optional :: info !! Zero on success; nonzero for an invalid/too-wide interval.
      integer(int64) :: width64
      integer :: i, ierr, offset
      width64 = int(max_value, int64) - int(min_value, int64) + 1_int64
      if (width64 <= 0_int64 .or. width64 > int(huge(1), int64)) then
         x = 0
         call set_info(info, 1)
         return
      end if
      do i = 1, size(x)
         call bounded_index(self, int(width64), offset, ierr)
         x(i) = min_value + offset
      end do
      call set_info(info, 0)
   end subroutine rng_int

   pure subroutine rng_perm(self, p, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying incremental Fisher-Yates choices.
      integer, intent(out) :: p(:) !! Random permutation of `1:size(p)`.
      integer, intent(out), optional :: info !! Zero on success.
      integer :: i, ierr, j
      do i = 1, size(p)
         call bounded_index(self, i, j, ierr)
         j = j + 1
         if (j /= i) p(i) = p(j)
         p(j) = i
      end do
      call set_info(info, 0)
   end subroutine rng_perm

   pure subroutine rng_sample(self, n, p, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying sample-without-replacement choices.
      integer, intent(in) :: n !! Population size; population labels are `1..n`.
      integer, intent(out) :: p(:) !! Sample of `size(p)` distinct population labels.
      integer, intent(out), optional :: info !! Zero on success; nonzero when sample size is outside `0..n`.
      integer :: i, j, k, pick, ierr
      logical :: seen
      k = size(p)
      if (n < 0 .or. k > n) then
      p = 0
      call set_info(info, 1)
      return
      end if
      if (k < n / 2) then
         do i = 1, k
            j = n - k + i
            call bounded_index(self, j, pick, ierr)
            pick = pick + 1
            seen = .false.
            if (i > 1) seen = any(p(1:i - 1) == pick)
            if (seen) then
               p(i) = j
            else
               p(i) = pick
            end if
         end do
      else
         if (k > 0) p(1) = 1
         do i = 1, k
            call bounded_index(self, i, pick, ierr)
            pick = pick + 1
            if (pick /= i) p(i) = p(pick)
            p(pick) = i
         end do
         do i = k + 1, n
            call bounded_index(self, i, pick, ierr)
            pick = pick + 1
            if (pick <= k) p(pick) = i
         end do
      end if
      call set_info(info, 0)
   end subroutine rng_sample

   pure subroutine rng_raw(self, bytes, info)
      class(randompack_rng_type), intent(inout) :: self !! RNG object supplying raw little-endian bytes from engine output words.
      integer(int8), intent(out) :: bytes(:) !! Output raw bytes; bit patterns match unsigned values 0..255.
      integer, intent(out), optional :: info !! Zero on success.
      integer(int64) :: b
      integer :: i
      do i = 1, size(bytes)
         call draw_byte(self, b)
         bytes(i) = int(iand(b, 127_int64), int8)
         if (b >= 128_int64) bytes(i) = int(b - 256_int64, int8)
      end do
      call set_info(info, 0)
   end subroutine rng_raw
end module randompack_rng_mod
