// -*- C -*-
// Combined SLEEF functions (AVX2).

#if !defined(__x86_64__) && !defined(_M_X64)
// No AVX2 on this architecture.
#else

#include <immintrin.h>
#include <math.h>
#include <stddef.h>
#include "randompack_config.h"

#define STATINLINE ALWAYS_INLINE // MinGW crashes on the out-of-line AVX2 vector helper form.
#define STATINLINECONST ALWAYS_INLINE CONST_ATTR
#define SLEEF_CONST CONST_ATTR
#define SLEEF_INF INF_VALUE
#define SLEEF_INFF INFF_VALUE
#define SLEEF_NAN NAN_VALUE
#define SLEEF_NANF NANF_VALUE

typedef __m256i vmask;
typedef __m256i vopmask;
typedef __m256d vdouble;
typedef __m128i vint;
typedef __m256 vfloat;
typedef __m256i vint2;
typedef struct {
  vdouble x, y;
} vdouble2;
typedef struct {
  vfloat x, y;
} vfloat2;

STATINLINE vdouble vcast_vd_d(double d) { return _mm256_set1_pd(d); }
STATINLINE vmask vreinterpret_vm_vd(vdouble vd) { return _mm256_castpd_si256(vd); }
STATINLINE vdouble vreinterpret_vd_vm(vmask vm) { return _mm256_castsi256_pd(vm); }
STATINLINE vmask vandnot_vm_vo64_vm(vopmask x, vmask y) { return
    vreinterpret_vm_vd(_mm256_andnot_pd(vreinterpret_vd_vm(x), vreinterpret_vd_vm(y))); }
STATINLINE vint vrint_vi_vd(vdouble vd) { return _mm256_cvtpd_epi32(vd); }
STATINLINE vdouble vrint_vd_vd(vdouble vd) { return _mm256_round_pd(vd, 0x00 |0x08); }
STATINLINE vint vcast_vi_i(int i) { return _mm_set1_epi32(i); }
STATINLINE vmask vcastu_vm_vi(vint vi) {
  return _mm256_slli_epi64(_mm256_cvtepi32_epi64(vi), 32);
}
STATINLINE vdouble vmul_vd_vd_vd(vdouble x, vdouble y) { return _mm256_mul_pd(x, y); }
STATINLINE vdouble vmla_vd_vd_vd_vd(vdouble x, vdouble y, vdouble z) { return _mm256_fmadd_pd(x, y, z); }
STATINLINE vdouble vfma_vd_vd_vd_vd(vdouble x, vdouble y, vdouble z) { return _mm256_fmadd_pd(x, y, z); }
STATINLINE vopmask vlt_vo_vd_vd(vdouble x, vdouble y) { return vreinterpret_vm_vd(_mm256_cmp_pd(x, y, 0x11)); }
STATINLINE vopmask vgt_vo_vd_vd(vdouble x, vdouble y) { return vreinterpret_vm_vd(_mm256_cmp_pd(x, y, 0x1e)); }
STATINLINE vint vadd_vi_vi_vi(vint x, vint y) { return _mm_add_epi32(x, y); }
STATINLINE vint vsub_vi_vi_vi(vint x, vint y) { return _mm_sub_epi32(x, y); }
STATINLINE vint vsll_vi_vi_i(vint x, int c) { return _mm_slli_epi32(x, c); }
STATINLINE vint vsra_vi_vi_i(vint x, int c) { return _mm_srai_epi32(x, c); }
STATINLINE vdouble vsel_vd_vo_vd_vd(vopmask o, vdouble x, vdouble y) { return _mm256_blendv_pd(y, x, _mm256_castsi256_pd(o)); }

STATINLINE SLEEF_CONST vdouble vpow2i_vd_vi(vint q) {
  q = vadd_vi_vi_vi(vcast_vi_i(0x3ff), q);
  vmask r = vcastu_vm_vi(vsll_vi_vi_i(q, 20));
  return vreinterpret_vd_vm(r);
}
STATINLINE SLEEF_CONST vdouble vldexp2_vd_vd_vi(vdouble d, vint e) {
  return vmul_vd_vd_vd(vmul_vd_vd_vd(d, vpow2i_vd_vi(vsra_vi_vi_i(e, 1))), vpow2i_vd_vi(vsub_vi_vi_vi(e, vsra_vi_vi_i(e, 1))));
}

STATINLINECONST vdouble Sleef_expd4_u10avx2(vdouble d) {
  vdouble u = vrint_vd_vd(vmul_vd_vd_vd(d, vcast_vd_d(1.442695040888963407359924681001892137426645954152985934135449406931))), s;
  vint q = vrint_vi_vd(u);
  s = vmla_vd_vd_vd_vd(u, vcast_vd_d(-.69314718055966295651160180568695068359375), d);
  s = vmla_vd_vd_vd_vd(u, vcast_vd_d(-.28235290563031577122588448175013436025525412068e-12), s);
  vdouble s2 = vmul_vd_vd_vd(s, s), s4 = vmul_vd_vd_vd(s2, s2), s8 = vmul_vd_vd_vd(s4, s4);
    u = vmla_vd_vd_vd_vd((s8), (vmla_vd_vd_vd_vd((s),
    (vcast_vd_d(+0.2081276378237164457e-8)), (vcast_vd_d(+0.2511210703042288022e-7)))),
    (vmla_vd_vd_vd_vd((s4), (vmla_vd_vd_vd_vd((s2), (vmla_vd_vd_vd_vd((s),
    (vcast_vd_d(+0.2755762628169491192e-6)), (vcast_vd_d(+0.2755723402025388239e-5)))),
    (vmla_vd_vd_vd_vd((s), (vcast_vd_d(+0.2480158687479686264e-4)),
    (vcast_vd_d(+0.1984126989855865850e-3)))))), (vmla_vd_vd_vd_vd((s2),
    (vmla_vd_vd_vd_vd((s), (vcast_vd_d(+0.1388888888914497797e-2)),
    (vcast_vd_d(+0.8333333333314938210e-2)))), (vmla_vd_vd_vd_vd((s),
    (vcast_vd_d(+0.4166666666666602598e-1)), (vcast_vd_d(+0.1666666666666669072e+0)))))))));
  u = vfma_vd_vd_vd_vd(u, s, vcast_vd_d(+0.5000000000000000000e+0));
  u = vfma_vd_vd_vd_vd(u, s, vcast_vd_d(+0.1000000000000000000e+1));
  u = vfma_vd_vd_vd_vd(u, s, vcast_vd_d(+0.1000000000000000000e+1));
  u = vldexp2_vd_vd_vi(u, q);
  vopmask o = vgt_vo_vd_vd(d, vcast_vd_d(0x1.62e42fefa39efp+9));
  u = vsel_vd_vo_vd_vd(o, vcast_vd_d(SLEEF_INF), u);
  u = vreinterpret_vd_vm(vandnot_vm_vo64_vm(vlt_vo_vd_vd(d, vcast_vd_d(-1000)), vreinterpret_vm_vd(u)));
  return u;
}

STATINLINE vmask vreinterpret_vm_vd_sp(vdouble vd) { return _mm256_castpd_si256(vd); }
STATINLINE vdouble vreinterpret_vd_vm_sp(vmask vm) { return _mm256_castsi256_pd(vm); }
STATINLINE vmask vandnot_vm_vo32_vm_sp(vopmask x, vmask y) { return
    vreinterpret_vm_vd_sp(_mm256_andnot_pd(vreinterpret_vd_vm_sp(x),
    vreinterpret_vd_vm_sp(y))); }
STATINLINE vint2 vcast_vi2_vm_sp(vmask vm) { return vm; }
STATINLINE vmask vcast_vm_vi2_sp(vint2 vi) { return vi; }
STATINLINE vint2 vrint_vi2_vf_sp(vfloat vf) { return vcast_vi2_vm_sp(_mm256_cvtps_epi32(vf)); }
STATINLINE vfloat vcast_vf_vi2_sp(vint2 vi) { return _mm256_cvtepi32_ps(vcast_vm_vi2_sp(vi)); }
STATINLINE vfloat vcast_vf_f_sp(float f) { return _mm256_set1_ps(f); }
STATINLINE vint2 vcast_vi2_i_sp(int i) { return _mm256_set1_epi32(i); }
STATINLINE vmask vreinterpret_vm_vf_sp(vfloat vf) { return _mm256_castps_si256(vf); }
STATINLINE vfloat vreinterpret_vf_vm_sp(vmask vm) { return _mm256_castsi256_ps(vm); }
STATINLINE vfloat vreinterpret_vf_vi2_sp(vint2 vi) { return vreinterpret_vf_vm_sp(vcast_vm_vi2_sp(vi)); }
STATINLINE vfloat vadd_vf_vf_vf_sp(vfloat x, vfloat y) { return _mm256_add_ps(x, y); }
STATINLINE vfloat vmul_vf_vf_vf_sp(vfloat x, vfloat y) { return _mm256_mul_ps(x, y); }
STATINLINE vfloat vmla_vf_vf_vf_vf_sp(vfloat x, vfloat y, vfloat z) { return _mm256_fmadd_ps(x, y, z); }
STATINLINE vopmask vlt_vo_vf_vf_sp(vfloat x, vfloat y) { return
    vreinterpret_vm_vf_sp(_mm256_cmp_ps(x, y, 0x11 )); }
STATINLINE vint2 vadd_vi2_vi2_vi2_sp(vint2 x, vint2 y) { return _mm256_add_epi32(x, y); }
STATINLINE vint2 vsub_vi2_vi2_vi2_sp(vint2 x, vint2 y) { return _mm256_sub_epi32(x, y); }
STATINLINE vint2 vsll_vi2_vi2_i_sp(vint2 x, int c) { return _mm256_slli_epi32(x, c); }
STATINLINE vint2 vsra_vi2_vi2_i_sp(vint2 x, int c) { return _mm256_srai_epi32(x, c); }
STATINLINE vfloat vsel_vf_vo_vf_vf_sp(vopmask o, vfloat x, vfloat y) { return _mm256_blendv_ps(y, x, _mm256_castsi256_ps(o)); }
STATINLINE SLEEF_CONST vfloat vpow2i_vf_vi2_sp(vint2 q) {
  return vreinterpret_vf_vi2_sp(vsll_vi2_vi2_i_sp(vadd_vi2_vi2_vi2_sp(q, vcast_vi2_i_sp(0x7f)), 23));
}
STATINLINE SLEEF_CONST vfloat vldexp2_vf_vf_vi2_sp(vfloat d, vint2 e) {
  return vmul_vf_vf_vf_sp(vmul_vf_vf_vf_sp(d, vpow2i_vf_vi2_sp(vsra_vi2_vi2_i_sp(e, 1))),
  vpow2i_vf_vi2_sp(vsub_vi2_vi2_vi2_sp(e, vsra_vi2_vi2_i_sp(e, 1))));
}

STATINLINECONST vfloat Sleef_expf8_u10avx2(vfloat d) {
  vint2 q = vrint_vi2_vf_sp(vmul_vf_vf_vf_sp(d, vcast_vf_f_sp(1.442695040888963407359924681001892137426645954152985934135449406931f)));
  vfloat s, u;
  s = vmla_vf_vf_vf_vf_sp(vcast_vf_vi2_sp(q), vcast_vf_f_sp(-0.693145751953125f), d);
  s = vmla_vf_vf_vf_vf_sp(vcast_vf_vi2_sp(q), vcast_vf_f_sp(-1.428606765330187045e-06f), s);
  u = vcast_vf_f_sp(0.000198527617612853646278381);
  u = vmla_vf_vf_vf_vf_sp(u, s, vcast_vf_f_sp(0.00139304355252534151077271));
  u = vmla_vf_vf_vf_vf_sp(u, s, vcast_vf_f_sp(0.00833336077630519866943359));
  u = vmla_vf_vf_vf_vf_sp(u, s, vcast_vf_f_sp(0.0416664853692054748535156));
  u = vmla_vf_vf_vf_vf_sp(u, s, vcast_vf_f_sp(0.166666671633720397949219));
  u = vmla_vf_vf_vf_vf_sp(u, s, vcast_vf_f_sp(0.5));
  u = vadd_vf_vf_vf_sp(vcast_vf_f_sp(1.0f), vmla_vf_vf_vf_vf_sp(vmul_vf_vf_vf_sp(s, s), u, s));
  u = vldexp2_vf_vf_vi2_sp(u, q);
  u = vreinterpret_vf_vm_sp(vandnot_vm_vo32_vm_sp(vlt_vo_vf_vf_sp(d, vcast_vf_f_sp(-104)), vreinterpret_vm_vf_sp(u)));
  u = vsel_vf_vo_vf_vf_sp(vlt_vo_vf_vf_sp(vcast_vf_f_sp(100), d), vcast_vf_f_sp(SLEEF_INFF), u);
  return u;
}

STATINLINE vopmask vor_vo_vo_vo(vopmask x, vopmask y) { return vreinterpret_vm_vd(_mm256_or_pd(vreinterpret_vd_vm(x), vreinterpret_vd_vm(y))); }

STATINLINE vopmask vcast_vo32_vo64(vopmask o) {
  return _mm256_permutevar8x32_epi32(o, _mm256_set_epi32(0, 0, 0, 0, 6, 4, 2, 0));
}

STATINLINE vdouble vcast_vd_vi(vint vi) { return _mm256_cvtepi32_pd(vi); }

STATINLINE vint vcastu_vi_vm(vmask vi) {
  return _mm_or_si128(_mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(_mm256_castsi256_si128(vi)), _mm_set1_ps(0), 0x0d)),
          _mm_castps_si128(_mm_shuffle_ps(_mm_set1_ps(0), _mm_castsi128_ps(_mm256_extractf128_si256(vi, 1)), 0xd0)));
}

STATINLINE vmask vadd64_vm_vm_vm(vmask x, vmask y) { return _mm256_add_epi64(x, y); }
STATINLINE vdouble vadd_vd_vd_vd(vdouble x, vdouble y) { return _mm256_add_pd(x, y); }
STATINLINE vdouble vsub_vd_vd_vd(vdouble x, vdouble y) { return _mm256_sub_pd(x, y); }
STATINLINE vdouble vrec_vd_vd(vdouble x) { return _mm256_div_pd(_mm256_set1_pd(1), x); }
STATINLINE vdouble vfmapn_vd_vd_vd_vd(vdouble x, vdouble y, vdouble z) { return _mm256_fmsub_pd(x, y, z); }
STATINLINE vdouble vfmanp_vd_vd_vd_vd(vdouble x, vdouble y, vdouble z) { return _mm256_fnmadd_pd(x, y, z); }
STATINLINE vopmask veq_vo_vd_vd(vdouble x, vdouble y) { return vreinterpret_vm_vd(_mm256_cmp_pd(x, y, 0x00)); }
STATINLINE vint vneg_vi_vi(vint e) { return vsub_vi_vi_vi(vcast_vi_i(0), e); }
STATINLINE vint vand_vi_vi_vi(vint x, vint y) { return _mm_and_si128(x, y); }
STATINLINE vint vsrl_vi_vi_i(vint x, int c) { return _mm_srli_epi32(x, c); }
STATINLINE vint vsel_vi_vo_vi_vi(vopmask m, vint x, vint y) { return _mm_blendv_epi8(y, x, _mm256_castsi256_si128(m)); }
STATINLINE vopmask vispinf_vo_vd(vdouble d) {
  return vreinterpret_vm_vd(_mm256_cmp_pd(d, _mm256_set1_pd(SLEEF_INF), 0x00));
}
STATINLINE vopmask visnan_vo_vd(vdouble d) {
  return vreinterpret_vm_vd(_mm256_cmp_pd(d, d, 0x04));
}
STATINLINE SLEEF_CONST vdouble vd2getx_vd_vd2(vdouble2 v) { return v.x; }
STATINLINE SLEEF_CONST vdouble vd2gety_vd_vd2(vdouble2 v) { return v.y; }
STATINLINE SLEEF_CONST vdouble2 vd2setxy_vd2_vd_vd(vdouble x, vdouble y) { vdouble2 v; v.x = x; v.y = y; return v; }
STATINLINE SLEEF_CONST vdouble2 vcast_vd2_d_d(double h, double l) {
  return vd2setxy_vd2_vd_vd(vcast_vd_d(h), vcast_vd_d(l));
}
STATINLINE SLEEF_CONST vdouble vadd_vd_3vd(vdouble v0, vdouble v1, vdouble v2) {
  return vadd_vd_vd_vd(vadd_vd_vd_vd(v0, v1), v2);
}

STATINLINE SLEEF_CONST vdouble vadd_vd_4vd(vdouble v0, vdouble v1, vdouble v2, vdouble v3) {
  return vadd_vd_3vd(vadd_vd_vd_vd(v0, v1), v2, v3);
}

STATINLINE SLEEF_CONST vdouble2 ddscale_vd2_vd2_vd(vdouble2 d, vdouble s) {
  return vd2setxy_vd2_vd_vd(vmul_vd_vd_vd(vd2getx_vd_vd2(d), s), vmul_vd_vd_vd(vd2gety_vd_vd2(d), s));
}

STATINLINE SLEEF_CONST vdouble2 ddadd2_vd2_vd_vd(vdouble x, vdouble y) {
  vdouble s = vadd_vd_vd_vd(x, y);
  vdouble v = vsub_vd_vd_vd(s, x);
  return vd2setxy_vd2_vd_vd(s, vadd_vd_vd_vd(vsub_vd_vd_vd(x, vsub_vd_vd_vd(s, v)), vsub_vd_vd_vd(y, v)));
}

STATINLINE SLEEF_CONST vdouble2 ddadd_vd2_vd2_vd(vdouble2 x, vdouble y) {
  vdouble s = vadd_vd_vd_vd(vd2getx_vd_vd2(x), y);
  return vd2setxy_vd2_vd_vd(s, vadd_vd_3vd(vsub_vd_vd_vd(vd2getx_vd_vd2(x), s), y, vd2gety_vd_vd2(x)));
}

STATINLINE SLEEF_CONST vdouble2 ddadd_vd2_vd2_vd2(vdouble2 x, vdouble2 y) {
  vdouble s = vadd_vd_vd_vd(vd2getx_vd_vd2(x), vd2getx_vd_vd2(y));
  return vd2setxy_vd2_vd_vd(s, vadd_vd_4vd(vsub_vd_vd_vd(vd2getx_vd_vd2(x), s), vd2getx_vd_vd2(y), vd2gety_vd_vd2(x), vd2gety_vd_vd2(y)));
}

STATINLINE SLEEF_CONST vdouble2 dddiv_vd2_vd2_vd2(vdouble2 n, vdouble2 d) {
  vdouble t = vrec_vd_vd(vd2getx_vd_vd2(d));
  vdouble s = vmul_vd_vd_vd(vd2getx_vd_vd2(n), t);
  vdouble u = vfmapn_vd_vd_vd_vd(t, vd2getx_vd_vd2(n), s);
  vdouble v = vfmanp_vd_vd_vd_vd(vd2gety_vd_vd2(d), t, vfmanp_vd_vd_vd_vd(vd2getx_vd_vd2(d), t, vcast_vd_d(1)));
  return vd2setxy_vd2_vd_vd(s, vfma_vd_vd_vd_vd(s, v, vfma_vd_vd_vd_vd(vd2gety_vd_vd2(n), t, u)));
}

STATINLINE SLEEF_CONST vdouble2 ddmul_vd2_vd2_vd(vdouble2 x, vdouble y) {
  vdouble s = vmul_vd_vd_vd(vd2getx_vd_vd2(x), y);
  return vd2setxy_vd2_vd_vd(s, vfma_vd_vd_vd_vd(vd2gety_vd_vd2(x), y, vfmapn_vd_vd_vd_vd(vd2getx_vd_vd2(x), y, s)));
}

STATINLINE SLEEF_CONST vint vilogb2k_vi_vd(vdouble d) {
  vint q = vcastu_vi_vm(vreinterpret_vm_vd(d));
  q = vsrl_vi_vi_i(q, 20);
  q = vand_vi_vi_vi(q, vcast_vi_i(0x7ff));
  q = vsub_vi_vi_vi(q, vcast_vi_i(0x3ff));
  return q;
}

STATINLINE SLEEF_CONST vdouble vldexp3_vd_vd_vi(vdouble d, vint q) {
  return vreinterpret_vd_vm(vadd64_vm_vm_vm(vreinterpret_vm_vd(d),
    vcastu_vm_vi(vsll_vi_vi_i(q, 20))));
}

STATINLINE vdouble vdiv_vd_vd_vd(vdouble x, vdouble y) { return _mm256_div_pd(x, y); }

STATINLINECONST vdouble Sleef_logd4_u35avx2(vdouble d) {
  vdouble x, x2;
  vdouble t, m;
  vopmask o = vlt_vo_vd_vd(d, vcast_vd_d(0x1p-1022));
  d = vsel_vd_vo_vd_vd(o, vmul_vd_vd_vd(d, vcast_vd_d(4294967296.0 * 4294967296.0)), d);
  vint e = vilogb2k_vi_vd(vmul_vd_vd_vd(d, vcast_vd_d(1.0/0.75)));
  m = vldexp3_vd_vd_vi(d, vneg_vi_vi(e));
  e = vsel_vi_vo_vi_vi(vcast_vo32_vo64(o), vsub_vi_vi_vi(e, vcast_vi_i(64)), e);
  x = vdiv_vd_vd_vd(vsub_vd_vd_vd(m, vcast_vd_d(1)), vadd_vd_vd_vd(vcast_vd_d(1), m));
  x2 = vmul_vd_vd_vd(x, x);
  vdouble x4 = vmul_vd_vd_vd(x2, x2), x8 = vmul_vd_vd_vd(x4, x4), x3 = vmul_vd_vd_vd(x, x2);
  t = vmla_vd_vd_vd_vd((x8), (vmla_vd_vd_vd_vd((x4),
    (vcast_vd_d(0.153487338491425068243146)), (vmla_vd_vd_vd_vd((x2),
    (vcast_vd_d(0.152519917006351951593857)), (vcast_vd_d(0.181863266251982985677316)))))),
    (vmla_vd_vd_vd_vd((x4), (vmla_vd_vd_vd_vd((x2),
    (vcast_vd_d(0.222221366518767365905163)), (vcast_vd_d(0.285714294746548025383248)))),
    (vmla_vd_vd_vd_vd((x2), (vcast_vd_d(0.399999999950799600689777)),
    (vcast_vd_d(0.6666666666667778740063)))))));
  x = vmla_vd_vd_vd_vd(x, vcast_vd_d(2), vmul_vd_vd_vd(vcast_vd_d(0.693147180559945286226764), vcast_vd_vi(e)));
  x = vmla_vd_vd_vd_vd(x3, t, x);
  x = vsel_vd_vo_vd_vd(vispinf_vo_vd(d), vcast_vd_d(SLEEF_INF), x);
  x = vsel_vd_vo_vd_vd(vor_vo_vo_vo(vlt_vo_vd_vd(d, vcast_vd_d(0)), visnan_vo_vd(d)), vcast_vd_d(SLEEF_NAN), x);
  x = vsel_vd_vo_vd_vd(veq_vo_vd_vd(d, vcast_vd_d(0)), vcast_vd_d(-SLEEF_INF), x);
  return x;
}

STATINLINE vopmask vor_vo_vo_vo_sp(vopmask x, vopmask y) { return
    vreinterpret_vm_vd_sp(_mm256_or_pd(vreinterpret_vd_vm_sp(x), vreinterpret_vd_vm_sp(y))); }

STATINLINE vint2 vreinterpret_vi2_vf_sp(vfloat vf) { return vcast_vi2_vm_sp(vreinterpret_vm_vf_sp(vf)); }
STATINLINE vfloat vsub_vf_vf_vf_sp(vfloat x, vfloat y) { return _mm256_sub_ps(x, y); }
STATINLINE vfloat vdiv_vf_vf_vf_sp(vfloat x, vfloat y) { return _mm256_div_ps(x, y); }
STATINLINE vfloat vrec_vf_vf_sp(vfloat x) { return vdiv_vf_vf_vf_sp(vcast_vf_f_sp(1.0f), x); }
STATINLINE vfloat vfma_vf_vf_vf_vf_sp(vfloat x, vfloat y, vfloat z) { return _mm256_fmadd_ps(x, y, z); }
STATINLINE vfloat vfmapn_vf_vf_vf_vf_sp(vfloat x, vfloat y, vfloat z) { return _mm256_fmsub_ps(x, y, z); }
STATINLINE vfloat vfmanp_vf_vf_vf_vf_sp(vfloat x, vfloat y, vfloat z) { return _mm256_fnmadd_ps(x, y, z); }
STATINLINE vopmask veq_vo_vf_vf_sp(vfloat x, vfloat y) { return
    vreinterpret_vm_vf_sp(_mm256_cmp_ps(x, y, 0x00 )); }
STATINLINE vopmask vneq_vo_vf_vf_sp(vfloat x, vfloat y) { return
    vreinterpret_vm_vf_sp(_mm256_cmp_ps(x, y, 0x04 )); }
STATINLINE vint2 vneg_vi2_vi2_sp(vint2 e) { return vsub_vi2_vi2_vi2_sp(vcast_vi2_i_sp(0), e); }
STATINLINE vint2 vand_vi2_vi2_vi2_sp(vint2 x, vint2 y) { return _mm256_and_si256(x, y); }
STATINLINE vint2 vsrl_vi2_vi2_i_sp(vint2 x, int c) { return _mm256_srli_epi32(x, c); }
STATINLINE vint2 vsel_vi2_vo_vi2_vi2_sp(vopmask m, vint2 x, vint2 y) {
  return _mm256_blendv_epi8(y, x, m);
}
STATINLINE vopmask vispinf_vo_vf_sp(vfloat d) { return veq_vo_vf_vf_sp(d, vcast_vf_f_sp(SLEEF_INFF)); }
STATINLINE vopmask visnan_vo_vf_sp(vfloat d) { return vneq_vo_vf_vf_sp(d, d); }
STATINLINE SLEEF_CONST vfloat vf2getx_vf_vf2_sp(vfloat2 v) { return v.x; }
STATINLINE SLEEF_CONST vfloat vf2gety_vf_vf2_sp(vfloat2 v) { return v.y; }
STATINLINE SLEEF_CONST vfloat2 vf2setxy_vf2_vf_vf_sp(vfloat x, vfloat y) { vfloat2 v; v.x = x; v.y = y; return v; }
STATINLINE SLEEF_CONST vfloat2 vcast_vf2_f_f_sp(float h, float l) {
  return vf2setxy_vf2_vf_vf_sp(vcast_vf_f_sp(h), vcast_vf_f_sp(l));
}
STATINLINE SLEEF_CONST vfloat vadd_vf_3vf_sp(vfloat v0, vfloat v1, vfloat v2) {
  return vadd_vf_vf_vf_sp(vadd_vf_vf_vf_sp(v0, v1), v2);
}
STATINLINE SLEEF_CONST vfloat vadd_vf_4vf_sp(vfloat v0, vfloat v1, vfloat v2, vfloat v3) {
  return vadd_vf_3vf_sp(vadd_vf_vf_vf_sp(v0, v1), v2, v3);
}
STATINLINE SLEEF_CONST vfloat2 dfscale_vf2_vf2_vf_sp(vfloat2 d, vfloat s) {
  return vf2setxy_vf2_vf_vf_sp(vmul_vf_vf_vf_sp(vf2getx_vf_vf2_sp(d), s), vmul_vf_vf_vf_sp(vf2gety_vf_vf2_sp(d), s));
}
STATINLINE SLEEF_CONST vfloat2 dfadd2_vf2_vf_vf_sp(vfloat x, vfloat y) {
  vfloat s = vadd_vf_vf_vf_sp(x, y);
  vfloat v = vsub_vf_vf_vf_sp(s, x);
  return vf2setxy_vf2_vf_vf_sp(s, vadd_vf_vf_vf_sp(vsub_vf_vf_vf_sp(x, vsub_vf_vf_vf_sp(s, v)), vsub_vf_vf_vf_sp(y, v)));
}
STATINLINE SLEEF_CONST vfloat2 dfadd_vf2_vf2_vf_sp(vfloat2 x, vfloat y) {
  vfloat s = vadd_vf_vf_vf_sp(vf2getx_vf_vf2_sp(x), y);
  return vf2setxy_vf2_vf_vf_sp(s, vadd_vf_3vf_sp(vsub_vf_vf_vf_sp(vf2getx_vf_vf2_sp(x), s), y, vf2gety_vf_vf2_sp(x)));
}
STATINLINE SLEEF_CONST vfloat2 dfadd_vf2_vf2_vf2_sp(vfloat2 x, vfloat2 y) {
  vfloat s = vadd_vf_vf_vf_sp(vf2getx_vf_vf2_sp(x), vf2getx_vf_vf2_sp(y));
  return vf2setxy_vf2_vf_vf_sp(s, vadd_vf_4vf_sp(vsub_vf_vf_vf_sp(vf2getx_vf_vf2_sp(x), s), vf2getx_vf_vf2_sp(y), vf2gety_vf_vf2_sp(x), vf2gety_vf_vf2_sp(y)));
}

STATINLINE SLEEF_CONST vfloat2 dfdiv_vf2_vf2_vf2_sp(vfloat2 n, vfloat2 d) {
  vfloat t = vrec_vf_vf_sp(vf2getx_vf_vf2_sp(d));
  vfloat s = vmul_vf_vf_vf_sp(vf2getx_vf_vf2_sp(n), t);
  vfloat u = vfmapn_vf_vf_vf_vf_sp(t, vf2getx_vf_vf2_sp(n), s);
  vfloat v = vfmanp_vf_vf_vf_vf_sp(vf2gety_vf_vf2_sp(d), t, vfmanp_vf_vf_vf_vf_sp(vf2getx_vf_vf2_sp(d), t, vcast_vf_f_sp(1)));
  return vf2setxy_vf2_vf_vf_sp(s, vfma_vf_vf_vf_vf_sp(s, v, vfma_vf_vf_vf_vf_sp(vf2gety_vf_vf2_sp(n), t, u)));
}

STATINLINE SLEEF_CONST vfloat2 dfmul_vf2_vf2_vf_sp(vfloat2 x, vfloat y) {
  vfloat s = vmul_vf_vf_vf_sp(vf2getx_vf_vf2_sp(x), y);
  return vf2setxy_vf2_vf_vf_sp(s, vfma_vf_vf_vf_vf_sp(vf2gety_vf_vf2_sp(x), y, vfmapn_vf_vf_vf_vf_sp(vf2getx_vf_vf2_sp(x), y, s)));
}

STATINLINE SLEEF_CONST vint2 vilogb2k_vi2_vf_sp(vfloat d) {
  vint2 q = vreinterpret_vi2_vf_sp(d);
  q = vsrl_vi2_vi2_i_sp(q, 23);
  q = vand_vi2_vi2_vi2_sp(q, vcast_vi2_i_sp(0xff));
  q = vsub_vi2_vi2_vi2_sp(q, vcast_vi2_i_sp(0x7f));
  return q;
}

STATINLINE SLEEF_CONST vfloat vldexp3_vf_vf_vi2_sp(vfloat d, vint2 q) {
  return vreinterpret_vf_vi2_sp(vadd_vi2_vi2_vi2_sp(vreinterpret_vi2_vf_sp(d), vsll_vi2_vi2_i_sp(q, 23)));
}

STATINLINECONST vfloat Sleef_logf8_u10avx2(vfloat d) {
  vfloat2 x;
  vfloat t, m, x2;
  vopmask o = vlt_vo_vf_vf_sp(d, vcast_vf_f_sp(0x1p-126));
  d = vsel_vf_vo_vf_vf_sp(o, vmul_vf_vf_vf_sp(d, vcast_vf_f_sp(4294967296.0f * 4294967296.0f)), d);
  vint2 e = vilogb2k_vi2_vf_sp(vmul_vf_vf_vf_sp(d, vcast_vf_f_sp(1.0f/0.75f)));
  m = vldexp3_vf_vf_vi2_sp(d, vneg_vi2_vi2_sp(e));
  e = vsel_vi2_vo_vi2_vi2_sp(o, vsub_vi2_vi2_vi2_sp(e, vcast_vi2_i_sp(64)), e);
  vfloat2 s = dfmul_vf2_vf2_vf_sp(vcast_vf2_f_f_sp(0.69314718246459960938f, -1.904654323148236017e-09f), vcast_vf_vi2_sp(e));
  x = dfdiv_vf2_vf2_vf2_sp(dfadd2_vf2_vf_vf_sp(vcast_vf_f_sp(-1), m), dfadd2_vf2_vf_vf_sp(vcast_vf_f_sp(1), m));
  x2 = vmul_vf_vf_vf_sp(vf2getx_vf_vf2_sp(x), vf2getx_vf_vf2_sp(x));
  t = vcast_vf_f_sp(+0.3027294874e+0f);
  t = vmla_vf_vf_vf_vf_sp(t, x2, vcast_vf_f_sp(+0.3996108174e+0f));
  t = vmla_vf_vf_vf_vf_sp(t, x2, vcast_vf_f_sp(+0.6666694880e+0f));
  s = dfadd_vf2_vf2_vf2_sp(s, dfscale_vf2_vf2_vf_sp(x, vcast_vf_f_sp(2)));
  s = dfadd_vf2_vf2_vf_sp(s, vmul_vf_vf_vf_sp(vmul_vf_vf_vf_sp(x2, vf2getx_vf_vf2_sp(x)), t));
  vfloat r = vadd_vf_vf_vf_sp(vf2getx_vf_vf2_sp(s), vf2gety_vf_vf2_sp(s));
  r = vsel_vf_vo_vf_vf_sp(vispinf_vo_vf_sp(d), vcast_vf_f_sp(SLEEF_INFF), r);
  r = vsel_vf_vo_vf_vf_sp(vor_vo_vo_vo_sp(vlt_vo_vf_vf_sp(d, vcast_vf_f_sp(0)), visnan_vo_vf_sp(d)), vcast_vf_f_sp(SLEEF_NANF), r);
  r = vsel_vf_vo_vf_vf_sp(veq_vo_vf_vf_sp(d, vcast_vf_f_sp(0)), vcast_vf_f_sp(-SLEEF_INFF), r);
  return r;
}

HIDDEN void sleef_expd4_u10avx2_array(double *x, size_t len) {
  size_t i = 0;
  for (; i + 4 <= len; i += 4) {
    __m256d d = _mm256_loadu_pd(x + i);
    d = Sleef_expd4_u10avx2(d);
    _mm256_storeu_pd(x + i, d);
  }
  for (; i < len; i++) x[i] = exp(x[i]);
}

HIDDEN void sleef_logd4_u35avx2_array(double *x, size_t len) {
  size_t i = 0;
  for (; i + 4 <= len; i += 4) {
    __m256d d = _mm256_loadu_pd(x + i);
    d = Sleef_logd4_u35avx2(d);
    _mm256_storeu_pd(x + i, d);
  }
  for (; i < len; i++) x[i] = log(x[i]);
}

HIDDEN void sleef_expf8_u10avx2_array(float *x, size_t len) {
  size_t i = 0;
  for (; i + 8 <= len; i += 8) {
    __m256 d = _mm256_loadu_ps(x + i);
    d = Sleef_expf8_u10avx2(d);
    _mm256_storeu_ps(x + i, d);
  }
  for (; i < len; i++) x[i] = expf(x[i]);
}

HIDDEN void sleef_logf8_u10avx2_array(float *x, size_t len) {
  size_t i = 0;
  for (; i + 8 <= len; i += 8) {
    __m256 d = _mm256_loadu_ps(x + i);
    d = Sleef_logf8_u10avx2(d);
    _mm256_storeu_ps(x + i, d);
  }
  for (; i < len; i++) x[i] = logf(x[i]);
}
#endif
