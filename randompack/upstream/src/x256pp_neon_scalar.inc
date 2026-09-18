// -*- C -*-

#if HAVE_VEC
#define FAST_STEP_VEC(s0,s1,s2,s3,outv) do { \
  VEC_T r_; \
  VEC_T t_; \
  r_ = VEC_ADD(VEC_ROTL23(VEC_ADD((s0),(s3))),(s0)); \
  t_ = VEC_SHL17((s1)); \
  (s2) = VEC_XOR((s2),(s0)); \
  (s3) = VEC_XOR((s3),(s1)); \
  (s1) = VEC_XOR((s1),(s2)); \
  (s0) = VEC_XOR((s0),(s3)); \
  (s2) = VEC_XOR((s2),t_); \
  (s3) = VEC_ROTL45((s3)); \
  (outv) = r_; \
} while (0)
#endif

static inline ALIGN64 void fill_x256ppsimd(uint64_t *buf, size_t len,
  randompack_state *state) {
  uint64_t *out = buf;
  xo256 *st = &state->xo;
#if HAVE_NEON
  VEC_T a0 = VEC_LOAD(&st->s0[0]);
  VEC_T a1 = VEC_LOAD(&st->s1[0]);
  VEC_T a2 = VEC_LOAD(&st->s2[0]);
  VEC_T a3 = VEC_LOAD(&st->s3[0]);
  VEC_T b0 = VEC_LOAD(&st->s0[2]);
  VEC_T b1 = VEC_LOAD(&st->s1[2]);
  VEC_T b2 = VEC_LOAD(&st->s2[2]);
  VEC_T b3 = VEC_LOAD(&st->s3[2]);
  VEC_T c0 = VEC_LOAD(&st->s0[4]);
  VEC_T c1 = VEC_LOAD(&st->s1[4]);
  VEC_T c2 = VEC_LOAD(&st->s2[4]);
  VEC_T c3 = VEC_LOAD(&st->s3[4]);
  VEC_T d0 = VEC_LOAD(&st->s0[6]);
  VEC_T d1 = VEC_LOAD(&st->s1[6]);
  VEC_T d2 = VEC_LOAD(&st->s2[6]);
  VEC_T d3 = VEC_LOAD(&st->s3[6]);
  for (size_t i = 0; i < len; i += 8) {
    VEC_T ra, rb, rc, rd;
    FAST_STEP_VEC(a0, a1, a2, a3, ra);
    FAST_STEP_VEC(b0, b1, b2, b3, rb);
    FAST_STEP_VEC(c0, c1, c2, c3, rc);
    FAST_STEP_VEC(d0, d1, d2, d3, rd);
    VEC_STORE(out + i, ra);
    VEC_STORE(out + i + 2, rb);
    VEC_STORE(out + i + 4, rc);
    VEC_STORE(out + i + 6, rd);
  }
  VEC_STORE(&st->s0[0], a0);
  VEC_STORE(&st->s1[0], a1);
  VEC_STORE(&st->s2[0], a2);
  VEC_STORE(&st->s3[0], a3);
  VEC_STORE(&st->s0[2], b0);
  VEC_STORE(&st->s1[2], b1);
  VEC_STORE(&st->s2[2], b2);
  VEC_STORE(&st->s3[2], b3);
  VEC_STORE(&st->s0[4], c0);
  VEC_STORE(&st->s1[4], c1);
  VEC_STORE(&st->s2[4], c2);
  VEC_STORE(&st->s3[4], c3);
  VEC_STORE(&st->s0[6], d0);
  VEC_STORE(&st->s1[6], d1);
  VEC_STORE(&st->s2[6], d2);
  VEC_STORE(&st->s3[6], d3);
#else
  uint64_t a0 = st->s0[0], a1 = st->s1[0], a2 = st->s2[0], a3 = st->s3[0];
  uint64_t b0 = st->s0[1], b1 = st->s1[1], b2 = st->s2[1], b3 = st->s3[1];
  uint64_t c0 = st->s0[2], c1 = st->s1[2], c2 = st->s2[2], c3 = st->s3[2];
  uint64_t d0 = st->s0[3], d1 = st->s1[3], d2 = st->s2[3], d3 = st->s3[3];
  uint64_t e0 = st->s0[4], e1 = st->s1[4], e2 = st->s2[4], e3 = st->s3[4];
  uint64_t f0 = st->s0[5], f1 = st->s1[5], f2 = st->s2[5], f3 = st->s3[5];
  uint64_t g0 = st->s0[6], g1 = st->s1[6], g2 = st->s2[6], g3 = st->s3[6];
  uint64_t h0 = st->s0[7], h1 = st->s1[7], h2 = st->s2[7], h3 = st->s3[7];
  for (size_t i = 0; i < len; i += 8) {
    FAST_STEP_SCALAR(a0,a1,a2,a3,out[i]);
    FAST_STEP_SCALAR(b0,b1,b2,b3,out[i+1]);
    FAST_STEP_SCALAR(c0,c1,c2,c3,out[i+2]);
    FAST_STEP_SCALAR(d0,d1,d2,d3,out[i+3]);
    FAST_STEP_SCALAR(e0,e1,e2,e3,out[i+4]);
    FAST_STEP_SCALAR(f0,f1,f2,f3,out[i+5]);
    FAST_STEP_SCALAR(g0,g1,g2,g3,out[i+6]);
    FAST_STEP_SCALAR(h0,h1,h2,h3,out[i+7]);
  }
  st->s0[0]=a0; st->s1[0]=a1; st->s2[0]=a2; st->s3[0]=a3;
  st->s0[1]=b0; st->s1[1]=b1; st->s2[1]=b2; st->s3[1]=b3;
  st->s0[2]=c0; st->s1[2]=c1; st->s2[2]=c2; st->s3[2]=c3;
  st->s0[3]=d0; st->s1[3]=d1; st->s2[3]=d2; st->s3[3]=d3;
  st->s0[4]=e0; st->s1[4]=e1; st->s2[4]=e2; st->s3[4]=e3;
  st->s0[5]=f0; st->s1[5]=f1; st->s2[5]=f2; st->s3[5]=f3;
  st->s0[6]=g0; st->s1[6]=g1; st->s2[6]=g2; st->s3[6]=g3;
  st->s0[7]=h0; st->s1[7]=h1; st->s2[7]=h2; st->s3[7]=h3;
#endif
}
