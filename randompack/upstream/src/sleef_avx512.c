// -*- C -*-
// AVX-512 SLEEF-based in-place exp/log helpers.

#if !defined(__x86_64__) && !defined(_M_X64)
// No AVX-512 on this architecture.
#else

#include <immintrin.h>
#include <math.h>
#include <stddef.h>

#define STATINLINE static inline __attribute__((always_inline))
#define CONST __attribute__((const))

typedef __m512i vmask;
typedef __mmask8 vopmask;
typedef __m512d vdouble;
typedef __m256i vint;
typedef __mmask16 vopmaskf;
typedef __m512 vfloat;
typedef __m512i vint2;
typedef struct {
  vfloat x, y;
} vfloat2;

STATINLINE vdouble vcast_vd_d(double d) { return _mm512_set1_pd(d); }
STATINLINE vmask vreinterpret_vm_vd(vdouble vd) { return _mm512_castpd_si512(vd); }
STATINLINE vdouble vreinterpret_vd_vm(vmask vm) { return _mm512_castsi512_pd(vm); }
STATINLINE vint vcast_vi_i(int i) { return _mm256_set1_epi32(i); }
STATINLINE vmask vcastu_vm_vi(vint vi) { return _mm512_slli_epi64(_mm512_cvtepi32_epi64(vi), 32); }
STATINLINE vint vcastu_vi_vm(vmask vm) { return _mm512_cvtepi64_epi32(_mm512_srli_epi64(vm, 32)); }
STATINLINE vdouble vcast_vd_vi(vint vi) { return _mm512_cvtepi32_pd(vi); }
STATINLINE vint vrint_vi_vd(vdouble vd) {
  return _mm512_cvt_roundpd_epi32(vd, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
}
STATINLINE vdouble vrint_vd_vd(vdouble vd) {
  return _mm512_roundscale_pd(vd, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
}
STATINLINE vint vadd_vi_vi_vi(vint x, vint y) { return _mm256_add_epi32(x, y); }
STATINLINE vint vsub_vi_vi_vi(vint x, vint y) { return _mm256_sub_epi32(x, y); }
STATINLINE vint vneg_vi_vi(vint e) { return vsub_vi_vi_vi(vcast_vi_i(0), e); }
STATINLINE vint vand_vi_vi_vi(vint x, vint y) { return _mm256_and_si256(x, y); }
STATINLINE vint vsll_vi_vi_i(vint x, int c) { return _mm256_slli_epi32(x, c); }
STATINLINE vint vsra_vi_vi_i(vint x, int c) { return _mm256_srai_epi32(x, c); }
STATINLINE vint vsrl_vi_vi_i(vint x, int c) { return _mm256_srli_epi32(x, c); }
STATINLINE vmask vadd64_vm_vm_vm(vmask x, vmask y) { return _mm512_add_epi64(x, y); }
STATINLINE vdouble vadd_vd_vd_vd(vdouble x, vdouble y) { return _mm512_add_pd(x, y); }
STATINLINE vdouble vsub_vd_vd_vd(vdouble x, vdouble y) { return _mm512_sub_pd(x, y); }
STATINLINE vdouble vmul_vd_vd_vd(vdouble x, vdouble y) { return _mm512_mul_pd(x, y); }
STATINLINE vdouble vdiv_vd_vd_vd(vdouble x, vdouble y) { return _mm512_div_pd(x, y); }
STATINLINE vdouble vrec_vd_vd(vdouble x) { return _mm512_div_pd(vcast_vd_d(1), x); }
STATINLINE vdouble vmla_vd_vd_vd_vd(vdouble x, vdouble y, vdouble z) { return _mm512_fmadd_pd(x, y, z); }
STATINLINE vdouble vfma_vd_vd_vd_vd(vdouble x, vdouble y, vdouble z) { return _mm512_fmadd_pd(x, y, z); }
STATINLINE vdouble vfmapn_vd_vd_vd_vd(vdouble x, vdouble y, vdouble z) { return _mm512_fmsub_pd(x, y, z); }
STATINLINE vdouble vfmanp_vd_vd_vd_vd(vdouble x, vdouble y, vdouble z) { return _mm512_fnmadd_pd(x, y, z); }
STATINLINE vopmask vlt_vo_vd_vd(vdouble x, vdouble y) { return _mm512_cmp_pd_mask(x, y, _CMP_LT_OQ); }
STATINLINE vopmask vgt_vo_vd_vd(vdouble x, vdouble y) { return _mm512_cmp_pd_mask(x, y, _CMP_GT_OQ); }
STATINLINE vopmask veq_vo_vd_vd(vdouble x, vdouble y) { return _mm512_cmp_pd_mask(x, y, _CMP_EQ_OQ); }
STATINLINE vopmask visnan_vo_vd(vdouble d) { return _mm512_cmp_pd_mask(d, d, _CMP_NEQ_UQ); }
STATINLINE vopmask vispinf_vo_vd(vdouble d) { return _mm512_cmp_pd_mask(d, vcast_vd_d(__builtin_inf()), _CMP_EQ_OQ); }
STATINLINE vopmask vor_vo_vo_vo(vopmask x, vopmask y) { return x | y; }
STATINLINE vdouble vsel_vd_vo_vd_vd(vopmask o, vdouble x, vdouble y) { return _mm512_mask_blend_pd(o, y, x); }
STATINLINE vint vsel_vi_vo_vi_vi(vopmask m, vint x, vint y) {
  return _mm512_castsi512_si256(_mm512_mask_blend_epi32(m, _mm512_castsi256_si512(y), _mm512_castsi256_si512(x)));
}
STATINLINE vmask vandnot_vm_vo64_vm(vopmask o, vmask m) { return _mm512_mask_blend_epi64(o, m, _mm512_setzero_si512()); }

STATINLINE CONST vdouble vpow2i_vd_vi(vint q) {
  q = vadd_vi_vi_vi(vcast_vi_i(0x3ff), q);
  return vreinterpret_vd_vm(vcastu_vm_vi(vsll_vi_vi_i(q, 20)));
}

STATINLINE CONST vdouble vldexp2_vd_vd_vi(vdouble d, vint e) {
  vint e2 = vsra_vi_vi_i(e, 1);
  return vmul_vd_vd_vd(vmul_vd_vd_vd(d, vpow2i_vd_vi(e2)), vpow2i_vd_vi(vsub_vi_vi_vi(e, e2)));
}

STATINLINE CONST vint vilogb2k_vi_vd(vdouble d) {
  vint q = vcastu_vi_vm(vreinterpret_vm_vd(d));
  q = vsrl_vi_vi_i(q, 20);
  q = vand_vi_vi_vi(q, vcast_vi_i(0x7ff));
  q = vsub_vi_vi_vi(q, vcast_vi_i(0x3ff));
  return q;
}

STATINLINE CONST vdouble vldexp3_vd_vd_vi(vdouble d, vint q) {
  return vreinterpret_vd_vm(vadd64_vm_vm_vm(vreinterpret_vm_vd(d), vcastu_vm_vi(vsll_vi_vi_i(q, 20))));
}

CONST vdouble Sleef_expd8_u10avx512f(vdouble d) {
  vdouble u = vrint_vd_vd(vmul_vd_vd_vd(d, vcast_vd_d(1.442695040888963407359924681001892137426645954152985934135449406931)));
  vint q = vrint_vi_vd(u);
  vdouble s = vmla_vd_vd_vd_vd(u, vcast_vd_d(-.69314718055966295651160180568695068359375), d);
  s = vmla_vd_vd_vd_vd(u, vcast_vd_d(-.28235290563031577122588448175013436025525412068e-12), s);
  vdouble s2 = vmul_vd_vd_vd(s, s);
  vdouble s4 = vmul_vd_vd_vd(s2, s2);
  vdouble s8 = vmul_vd_vd_vd(s4, s4);
  u = vmla_vd_vd_vd_vd(s8, vmla_vd_vd_vd_vd(s, vcast_vd_d(+0.2081276378237164457e-8), vcast_vd_d(+0.2511210703042288022e-7)),
      vmla_vd_vd_vd_vd(s4, vmla_vd_vd_vd_vd(s2, vmla_vd_vd_vd_vd(s, vcast_vd_d(+0.2755762628169491192e-6), vcast_vd_d(+0.2755723402025388239e-5)),
      vmla_vd_vd_vd_vd(s, vcast_vd_d(+0.2480158687479686264e-4), vcast_vd_d(+0.1984126989855865850e-3))),
      vmla_vd_vd_vd_vd(s2, vmla_vd_vd_vd_vd(s, vcast_vd_d(+0.1388888888914497797e-2), vcast_vd_d(+0.8333333333314938210e-2)),
      vmla_vd_vd_vd_vd(s, vcast_vd_d(+0.4166666666666602598e-1), vcast_vd_d(+0.1666666666666669072)))));
  u = vfma_vd_vd_vd_vd(u, s, vcast_vd_d(+0.5));
  u = vfma_vd_vd_vd_vd(u, s, vcast_vd_d(+1));
  u = vfma_vd_vd_vd_vd(u, s, vcast_vd_d(+1));
  u = vldexp2_vd_vd_vi(u, q);
  vopmask o = vgt_vo_vd_vd(d, vcast_vd_d(0x1.62e42fefa39efp+9));
  u = vsel_vd_vo_vd_vd(o, vcast_vd_d(__builtin_inf()), u);
  u = vreinterpret_vd_vm(vandnot_vm_vo64_vm(vlt_vo_vd_vd(d, vcast_vd_d(-1000)), vreinterpret_vm_vd(u)));
  return u;
}

CONST vdouble Sleef_logd8_u35avx512f(vdouble d) {
  vopmask o = vlt_vo_vd_vd(d, vcast_vd_d(0x1p-1022));
  d = vsel_vd_vo_vd_vd(o, vmul_vd_vd_vd(d, vcast_vd_d(4294967296.0 * 4294967296.0)), d);
  vint e = vilogb2k_vi_vd(vmul_vd_vd_vd(d, vcast_vd_d(1.0/0.75)));
  vdouble m = vldexp3_vd_vd_vi(d, vneg_vi_vi(e));
  e = vsel_vi_vo_vi_vi(o, vsub_vi_vi_vi(e, vcast_vi_i(64)), e);
  vdouble x = vdiv_vd_vd_vd(vsub_vd_vd_vd(m, vcast_vd_d(1)), vadd_vd_vd_vd(vcast_vd_d(1), m));
  vdouble x2 = vmul_vd_vd_vd(x, x);
  vdouble x4 = vmul_vd_vd_vd(x2, x2);
  vdouble x8 = vmul_vd_vd_vd(x4, x4);
  vdouble x3 = vmul_vd_vd_vd(x, x2);
  vdouble t = vmla_vd_vd_vd_vd(x8, vmla_vd_vd_vd_vd(x4, vcast_vd_d(0.153487338491425068243146),
      vmla_vd_vd_vd_vd(x2, vcast_vd_d(0.152519917006351951593857), vcast_vd_d(0.181863266251982985677316))),
      vmla_vd_vd_vd_vd(x4, vmla_vd_vd_vd_vd(x2, vcast_vd_d(0.222221366518767365905163), vcast_vd_d(0.285714294746548025383248)),
      vmla_vd_vd_vd_vd(x2, vcast_vd_d(0.399999999950799600689777), vcast_vd_d(0.6666666666667778740063))));
  x = vmla_vd_vd_vd_vd(x, vcast_vd_d(2), vmul_vd_vd_vd(vcast_vd_d(0.693147180559945286226764), vcast_vd_vi(e)));
  x = vmla_vd_vd_vd_vd(x3, t, x);
  x = vsel_vd_vo_vd_vd(vispinf_vo_vd(d), vcast_vd_d(__builtin_inf()), x);
  x = vsel_vd_vo_vd_vd(vor_vo_vo_vo(vlt_vo_vd_vd(d, vcast_vd_d(0)), visnan_vo_vd(d)), vcast_vd_d(__builtin_nan("")), x);
  x = vsel_vd_vo_vd_vd(veq_vo_vd_vd(d, vcast_vd_d(0)), vcast_vd_d(-__builtin_inf()), x);
  return x;
}

void sleef_exp_inplace_avx512(double *x, size_t len) {
  size_t i = 0;
  for (; i + 8 <= len; i += 8) {
    __m512d d = _mm512_loadu_pd(x + i);
    d = Sleef_expd8_u10avx512f(d);
    _mm512_storeu_pd(x + i, d);
  }
  for (; i < len; i++) x[i] = exp(x[i]);
}

void sleef_log_inplace_avx512(double *x, size_t len) {
  size_t i = 0;
  for (; i + 8 <= len; i += 8) {
    __m512d d = _mm512_loadu_pd(x + i);
    d = Sleef_logd8_u35avx512f(d);
    _mm512_storeu_pd(x + i, d);
  }
  for (; i < len; i++) x[i] = log(x[i]);
}

STATINLINE vmask vreinterpret_vm_vf_sp(vfloat vf) { return _mm512_castps_si512(vf); }
STATINLINE vfloat vreinterpret_vf_vm_sp(vmask vm) { return _mm512_castsi512_ps(vm); }
STATINLINE vmask vandnot_vm_vo32_vm_sp(vopmaskf x, vmask y) {
  return _mm512_mask_blend_epi32(x, y, _mm512_setzero_si512());
}
STATINLINE vint2 vcast_vi2_vm_sp(vmask vm) { return vm; }
STATINLINE vmask vcast_vm_vi2_sp(vint2 vi) { return vi; }
STATINLINE vint2 vrint_vi2_vf_sp(vfloat vf) {
  return _mm512_cvt_roundps_epi32(vf, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
}
STATINLINE vfloat vcast_vf_vi2_sp(vint2 vi) { return _mm512_cvtepi32_ps(vcast_vm_vi2_sp(vi)); }
STATINLINE vfloat vcast_vf_f_sp(float f) { return _mm512_set1_ps(f); }
STATINLINE vint2 vcast_vi2_i_sp(int i) { return _mm512_set1_epi32(i); }
STATINLINE vfloat vreinterpret_vf_vi2_sp(vint2 vi) { return vreinterpret_vf_vm_sp(vcast_vm_vi2_sp(vi)); }
STATINLINE vfloat vadd_vf_vf_vf_sp(vfloat x, vfloat y) { return _mm512_add_ps(x, y); }
STATINLINE vfloat vmul_vf_vf_vf_sp(vfloat x, vfloat y) { return _mm512_mul_ps(x, y); }
STATINLINE vfloat vmla_vf_vf_vf_vf_sp(vfloat x, vfloat y, vfloat z) {
  return _mm512_fmadd_ps(x, y, z);
}
STATINLINE vopmaskf vlt_vo_vf_vf_sp(vfloat x, vfloat y) {
  return _mm512_cmp_ps_mask(x, y, _CMP_LT_OQ);
}
STATINLINE vint2 vadd_vi2_vi2_vi2_sp(vint2 x, vint2 y) { return _mm512_add_epi32(x, y); }
STATINLINE vint2 vsub_vi2_vi2_vi2_sp(vint2 x, vint2 y) { return _mm512_sub_epi32(x, y); }
STATINLINE vint2 vsll_vi2_vi2_i_sp(vint2 x, int c) { return _mm512_slli_epi32(x, c); }
STATINLINE vint2 vsra_vi2_vi2_i_sp(vint2 x, int c) { return _mm512_srai_epi32(x, c); }
STATINLINE vfloat vsel_vf_vo_vf_vf_sp(vopmaskf o, vfloat x, vfloat y) {
  return _mm512_mask_blend_ps(o, y, x);
}
STATINLINE CONST vfloat vpow2i_vf_vi2_sp(vint2 q) {
  return vreinterpret_vf_vi2_sp(vsll_vi2_vi2_i_sp(
    vadd_vi2_vi2_vi2_sp(q, vcast_vi2_i_sp(0x7f)), 23));
}
STATINLINE CONST vfloat vldexp2_vf_vf_vi2_sp(vfloat d, vint2 e) {
  return vmul_vf_vf_vf_sp(vmul_vf_vf_vf_sp(d, vpow2i_vf_vi2_sp(vsra_vi2_vi2_i_sp(e, 1))),
    vpow2i_vf_vi2_sp(vsub_vi2_vi2_vi2_sp(e, vsra_vi2_vi2_i_sp(e, 1))));
}

CONST vfloat Sleef_expf16_u10avx512f(vfloat d) {
  vint2 q = vrint_vi2_vf_sp(vmul_vf_vf_vf_sp(d,
    vcast_vf_f_sp(1.442695040888963407359924681001892137426645954152985934135449406931f)));
  vfloat s = vmla_vf_vf_vf_vf_sp(vcast_vf_vi2_sp(q), vcast_vf_f_sp(-0.693145751953125f), d);
  s = vmla_vf_vf_vf_vf_sp(vcast_vf_vi2_sp(q), vcast_vf_f_sp(-1.428606765330187045e-06f), s);
  vfloat u = vcast_vf_f_sp(0.000198527617612853646278381f);
  u = vmla_vf_vf_vf_vf_sp(u, s, vcast_vf_f_sp(0.00139304355252534151077271f));
  u = vmla_vf_vf_vf_vf_sp(u, s, vcast_vf_f_sp(0.00833336077630519866943359f));
  u = vmla_vf_vf_vf_vf_sp(u, s, vcast_vf_f_sp(0.0416664853692054748535156f));
  u = vmla_vf_vf_vf_vf_sp(u, s, vcast_vf_f_sp(0.166666671633720397949219f));
  u = vmla_vf_vf_vf_vf_sp(u, s, vcast_vf_f_sp(0.5f));
  u = vadd_vf_vf_vf_sp(vcast_vf_f_sp(1), vmla_vf_vf_vf_vf_sp(vmul_vf_vf_vf_sp(s, s), u, s));
  u = vldexp2_vf_vf_vi2_sp(u, q);
  u = vreinterpret_vf_vm_sp(vandnot_vm_vo32_vm_sp(vlt_vo_vf_vf_sp(d, vcast_vf_f_sp(-104)),
    vreinterpret_vm_vf_sp(u)));
  u = vsel_vf_vo_vf_vf_sp(vlt_vo_vf_vf_sp(vcast_vf_f_sp(100), d),
    vcast_vf_f_sp(__builtin_inff()), u);
  return u;
}

STATINLINE vopmaskf vor_vo_vo_vo_sp(vopmaskf x, vopmaskf y) { return x | y; }
STATINLINE vint2 vreinterpret_vi2_vf_sp(vfloat vf) { return vcast_vi2_vm_sp(vreinterpret_vm_vf_sp(vf)); }
STATINLINE vfloat vsub_vf_vf_vf_sp(vfloat x, vfloat y) { return _mm512_sub_ps(x, y); }
STATINLINE vfloat vdiv_vf_vf_vf_sp(vfloat x, vfloat y) { return _mm512_div_ps(x, y); }
STATINLINE vfloat vrec_vf_vf_sp(vfloat x) { return vdiv_vf_vf_vf_sp(vcast_vf_f_sp(1), x); }
STATINLINE vfloat vfma_vf_vf_vf_vf_sp(vfloat x, vfloat y, vfloat z) {
  return _mm512_fmadd_ps(x, y, z);
}
STATINLINE vfloat vfmapn_vf_vf_vf_vf_sp(vfloat x, vfloat y, vfloat z) {
  return _mm512_fmsub_ps(x, y, z);
}
STATINLINE vfloat vfmanp_vf_vf_vf_vf_sp(vfloat x, vfloat y, vfloat z) {
  return _mm512_fnmadd_ps(x, y, z);
}
STATINLINE vopmaskf veq_vo_vf_vf_sp(vfloat x, vfloat y) {
  return _mm512_cmp_ps_mask(x, y, _CMP_EQ_OQ);
}
STATINLINE vopmaskf vneq_vo_vf_vf_sp(vfloat x, vfloat y) {
  return _mm512_cmp_ps_mask(x, y, _CMP_NEQ_UQ);
}
STATINLINE vint2 vneg_vi2_vi2_sp(vint2 e) { return vsub_vi2_vi2_vi2_sp(vcast_vi2_i_sp(0), e); }
STATINLINE vint2 vand_vi2_vi2_vi2_sp(vint2 x, vint2 y) { return _mm512_and_si512(x, y); }
STATINLINE vint2 vsrl_vi2_vi2_i_sp(vint2 x, int c) { return _mm512_srli_epi32(x, c); }
STATINLINE vint2 vsel_vi2_vo_vi2_vi2_sp(vopmaskf m, vint2 x, vint2 y) {
  return _mm512_mask_blend_epi32(m, y, x);
}
STATINLINE vopmaskf vispinf_vo_vf_sp(vfloat d) {
  return veq_vo_vf_vf_sp(d, vcast_vf_f_sp(__builtin_inff()));
}
STATINLINE vopmaskf visnan_vo_vf_sp(vfloat d) { return vneq_vo_vf_vf_sp(d, d); }
STATINLINE CONST vfloat vf2getx_vf_vf2_sp(vfloat2 v) { return v.x; }
STATINLINE CONST vfloat vf2gety_vf_vf2_sp(vfloat2 v) { return v.y; }
STATINLINE CONST vfloat2 vf2setxy_vf2_vf_vf_sp(vfloat x, vfloat y) {
  vfloat2 v;
  v.x = x;
  v.y = y;
  return v;
}
STATINLINE CONST vfloat2 vcast_vf2_f_f_sp(float h, float l) {
  return vf2setxy_vf2_vf_vf_sp(vcast_vf_f_sp(h), vcast_vf_f_sp(l));
}
STATINLINE CONST vfloat vadd_vf_3vf_sp(vfloat v0, vfloat v1, vfloat v2) {
  return vadd_vf_vf_vf_sp(vadd_vf_vf_vf_sp(v0, v1), v2);
}
STATINLINE CONST vfloat vadd_vf_4vf_sp(vfloat v0, vfloat v1, vfloat v2, vfloat v3) {
  return vadd_vf_3vf_sp(vadd_vf_vf_vf_sp(v0, v1), v2, v3);
}
STATINLINE CONST vfloat2 dfscale_vf2_vf2_vf_sp(vfloat2 d, vfloat s) {
  return vf2setxy_vf2_vf_vf_sp(vmul_vf_vf_vf_sp(vf2getx_vf_vf2_sp(d), s),
    vmul_vf_vf_vf_sp(vf2gety_vf_vf2_sp(d), s));
}
STATINLINE CONST vfloat2 dfadd2_vf2_vf_vf_sp(vfloat x, vfloat y) {
  vfloat s = vadd_vf_vf_vf_sp(x, y);
  vfloat v = vsub_vf_vf_vf_sp(s, x);
  return vf2setxy_vf2_vf_vf_sp(s, vadd_vf_vf_vf_sp(vsub_vf_vf_vf_sp(x,
    vsub_vf_vf_vf_sp(s, v)), vsub_vf_vf_vf_sp(y, v)));
}
STATINLINE CONST vfloat2 dfadd_vf2_vf2_vf_sp(vfloat2 x, vfloat y) {
  vfloat s = vadd_vf_vf_vf_sp(vf2getx_vf_vf2_sp(x), y);
  return vf2setxy_vf2_vf_vf_sp(s, vadd_vf_3vf_sp(vsub_vf_vf_vf_sp(vf2getx_vf_vf2_sp(x), s),
    y, vf2gety_vf_vf2_sp(x)));
}
STATINLINE CONST vfloat2 dfadd_vf2_vf2_vf2_sp(vfloat2 x, vfloat2 y) {
  vfloat s = vadd_vf_vf_vf_sp(vf2getx_vf_vf2_sp(x), vf2getx_vf_vf2_sp(y));
  return vf2setxy_vf2_vf_vf_sp(s, vadd_vf_4vf_sp(vsub_vf_vf_vf_sp(vf2getx_vf_vf2_sp(x), s),
    vf2getx_vf_vf2_sp(y), vf2gety_vf_vf2_sp(x), vf2gety_vf_vf2_sp(y)));
}
STATINLINE CONST vfloat2 dfdiv_vf2_vf2_vf2_sp(vfloat2 n, vfloat2 d) {
  vfloat t = vrec_vf_vf_sp(vf2getx_vf_vf2_sp(d));
  vfloat s = vmul_vf_vf_vf_sp(vf2getx_vf_vf2_sp(n), t);
  vfloat u = vfmapn_vf_vf_vf_vf_sp(t, vf2getx_vf_vf2_sp(n), s);
  vfloat v = vfmanp_vf_vf_vf_vf_sp(vf2gety_vf_vf2_sp(d), t,
    vfmanp_vf_vf_vf_vf_sp(vf2getx_vf_vf2_sp(d), t, vcast_vf_f_sp(1)));
  return vf2setxy_vf2_vf_vf_sp(s, vfma_vf_vf_vf_vf_sp(s, v,
    vfma_vf_vf_vf_vf_sp(vf2gety_vf_vf2_sp(n), t, u)));
}
STATINLINE CONST vfloat2 dfmul_vf2_vf2_vf_sp(vfloat2 x, vfloat y) {
  vfloat s = vmul_vf_vf_vf_sp(vf2getx_vf_vf2_sp(x), y);
  return vf2setxy_vf2_vf_vf_sp(s, vfma_vf_vf_vf_vf_sp(vf2gety_vf_vf2_sp(x), y,
    vfmapn_vf_vf_vf_vf_sp(vf2getx_vf_vf2_sp(x), y, s)));
}
STATINLINE CONST vint2 vilogb2k_vi2_vf_sp(vfloat d) {
  vint2 q = vreinterpret_vi2_vf_sp(d);
  q = vsrl_vi2_vi2_i_sp(q, 23);
  q = vand_vi2_vi2_vi2_sp(q, vcast_vi2_i_sp(0xff));
  q = vsub_vi2_vi2_vi2_sp(q, vcast_vi2_i_sp(0x7f));
  return q;
}
STATINLINE CONST vfloat vldexp3_vf_vf_vi2_sp(vfloat d, vint2 q) {
  return vreinterpret_vf_vi2_sp(vadd_vi2_vi2_vi2_sp(vreinterpret_vi2_vf_sp(d),
    vsll_vi2_vi2_i_sp(q, 23)));
}

CONST vfloat Sleef_logf16_u10avx512f(vfloat d) {
  vopmaskf o = vlt_vo_vf_vf_sp(d, vcast_vf_f_sp(0x1p-126));
  d = vsel_vf_vo_vf_vf_sp(o, vmul_vf_vf_vf_sp(d,
    vcast_vf_f_sp(4294967296.0f * 4294967296.0f)), d);
  vint2 e = vilogb2k_vi2_vf_sp(vmul_vf_vf_vf_sp(d, vcast_vf_f_sp(1.0f/0.75f)));
  vfloat m = vldexp3_vf_vf_vi2_sp(d, vneg_vi2_vi2_sp(e));
  e = vsel_vi2_vo_vi2_vi2_sp(o, vsub_vi2_vi2_vi2_sp(e, vcast_vi2_i_sp(64)), e);
  vfloat2 s = dfmul_vf2_vf2_vf_sp(vcast_vf2_f_f_sp(0.69314718246459960938f,
    -1.904654323148236017e-09f), vcast_vf_vi2_sp(e));
  vfloat2 x = dfdiv_vf2_vf2_vf2_sp(dfadd2_vf2_vf_vf_sp(vcast_vf_f_sp(-1), m),
    dfadd2_vf2_vf_vf_sp(vcast_vf_f_sp(1), m));
  vfloat x2 = vmul_vf_vf_vf_sp(vf2getx_vf_vf2_sp(x), vf2getx_vf_vf2_sp(x));
  vfloat t = vcast_vf_f_sp(+0.3027294874e+0f);
  t = vmla_vf_vf_vf_vf_sp(t, x2, vcast_vf_f_sp(+0.3996108174e+0f));
  t = vmla_vf_vf_vf_vf_sp(t, x2, vcast_vf_f_sp(+0.6666694880e+0f));
  s = dfadd_vf2_vf2_vf2_sp(s, dfscale_vf2_vf2_vf_sp(x, vcast_vf_f_sp(2)));
  s = dfadd_vf2_vf2_vf_sp(s, vmul_vf_vf_vf_sp(vmul_vf_vf_vf_sp(x2,
    vf2getx_vf_vf2_sp(x)), t));
  vfloat r = vadd_vf_vf_vf_sp(vf2getx_vf_vf2_sp(s), vf2gety_vf_vf2_sp(s));
  r = vsel_vf_vo_vf_vf_sp(vispinf_vo_vf_sp(d), vcast_vf_f_sp(__builtin_inff()), r);
  r = vsel_vf_vo_vf_vf_sp(vor_vo_vo_vo_sp(vlt_vo_vf_vf_sp(d, vcast_vf_f_sp(0)),
    visnan_vo_vf_sp(d)), vcast_vf_f_sp(__builtin_nanf("")), r);
  r = vsel_vf_vo_vf_vf_sp(veq_vo_vf_vf_sp(d, vcast_vf_f_sp(0)),
    vcast_vf_f_sp(-__builtin_inff()), r);
  return r;
}

void sleef_expf_inplace_avx512(float *x, size_t len) {
  size_t i = 0;
  for (; i + 16 <= len; i += 16) {
    __m512 d = _mm512_loadu_ps(x + i);
    d = Sleef_expf16_u10avx512f(d);
    _mm512_storeu_ps(x + i, d);
  }
  for (; i < len; i++) x[i] = expf(x[i]);
}

void sleef_logf_inplace_avx512(float *x, size_t len) {
  size_t i = 0;
  for (; i + 16 <= len; i += 16) {
    __m512 d = _mm512_loadu_ps(x + i);
    d = Sleef_logf16_u10avx512f(d);
    _mm512_storeu_ps(x + i, d);
  }
  for (; i < len; i++) x[i] = logf(x[i]);
}

#endif
