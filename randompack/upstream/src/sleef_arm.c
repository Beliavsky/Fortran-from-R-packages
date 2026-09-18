// -*- C -*-
// Combined SLEEF functions (AArch64 AdvSIMD).

#if !defined(__aarch64__) && !defined(_M_ARM64)
// No AdvSIMD on this architecture.
#else

#include <arm_neon.h>
#include <stdint.h>

#if defined(__NVCOMPILER)
  #define INLINE inline
#else
  #define INLINE inline __attribute__((always_inline))
#endif
#define CONST __attribute__((const))

typedef uint32x4_t vmask;
typedef uint32x4_t vopmask;
typedef float32x4_t vfloat;
typedef int32x4_t vint2;
typedef float64x2_t vdouble;
typedef int32x2_t vint;

static INLINE vfloat vcast_vf_f(float f) { return vdupq_n_f32(f); }
static INLINE vdouble vcast_vd_d(double f) { return vdupq_n_f64(f); }
static INLINE vmask vreinterpret_vm_vf(vfloat vf) {
  return vreinterpretq_u32_f32(vf);
}
static INLINE vfloat vreinterpret_vf_vm(vmask vm) {
  return vreinterpretq_f32_u32(vm);
}
static INLINE vmask vreinterpret_vm_vd(vdouble vd) {
  return vreinterpretq_u32_f64(vd);
}
static INLINE vdouble vreinterpret_vd_vm(vmask vm) {
  return vreinterpretq_f64_u32(vm);
}
static INLINE vint2 vreinterpret_vi2_vf(vfloat vf) {
  return vreinterpretq_s32_f32(vf);
}
static INLINE vfloat vreinterpret_vf_vi2(vint2 vi) {
  return vreinterpretq_f32_s32(vi);
}
static INLINE vfloat vadd_vf_vf_vf(vfloat x, vfloat y) { return vaddq_f32(x, y); }
static INLINE vfloat vsub_vf_vf_vf(vfloat x, vfloat y) { return vsubq_f32(x, y); }
static INLINE vfloat vmul_vf_vf_vf(vfloat x, vfloat y) { return vmulq_f32(x, y); }
static INLINE vfloat vdiv_vf_vf_vf(vfloat n, vfloat d) { return vdivq_f32(n, d); }
static INLINE vfloat vrec_vf_vf(vfloat d) {
  return vdiv_vf_vf_vf(vcast_vf_f(1), d);
}
static INLINE vfloat vmla_vf_vf_vf_vf(vfloat x, vfloat y, vfloat z) {
  return vfmaq_f32(z, x, y);
}
static INLINE vfloat vfma_vf_vf_vf_vf(vfloat x, vfloat y, vfloat z) {
  return vfmaq_f32(z, x, y);
}
static INLINE vfloat vfmanp_vf_vf_vf_vf(vfloat x, vfloat y, vfloat z) {
  return vfmsq_f32(z, x, y);
}
static INLINE vfloat vfmapn_vf_vf_vf_vf(vfloat x, vfloat y, vfloat z) {
  return vfmaq_f32(vnegq_f32(z), x, y);
}
static INLINE vdouble vadd_vd_vd_vd(vdouble x, vdouble y) { return vaddq_f64(x, y); }
static INLINE vdouble vsub_vd_vd_vd(vdouble x, vdouble y) { return vsubq_f64(x, y); }
static INLINE vdouble vmul_vd_vd_vd(vdouble x, vdouble y) { return vmulq_f64(x, y); }
static INLINE vdouble vdiv_vd_vd_vd(vdouble n, vdouble d) { return vdivq_f64(n, d); }
static INLINE vdouble vmla_vd_vd_vd_vd(vdouble x, vdouble y, vdouble z) {
  return vfmaq_f64(z, x, y);
}
static INLINE vdouble vfma_vd_vd_vd_vd(vdouble x, vdouble y, vdouble z) {
  return vfmaq_f64(z, x, y);
}
static INLINE vopmask vlt_vo_vf_vf(vfloat x, vfloat y) { return vcltq_f32(x, y); }
static INLINE vopmask vlt_vo_vd_vd(vdouble x, vdouble y) {
  return vreinterpretq_u32_u64(vcltq_f64(x, y));
}
static INLINE vopmask vgt_vo_vd_vd(vdouble x, vdouble y) {
  return vreinterpretq_u32_u64(vcgtq_f64(x, y));
}
static INLINE vopmask veq_vo_vf_vf(vfloat x, vfloat y) { return vceqq_f32(x, y); }
static INLINE vopmask veq_vo_vd_vd(vdouble x, vdouble y) {
  return vreinterpretq_u32_u64(vceqq_f64(x, y));
}
static INLINE vfloat vsel_vf_vo_vf_vf(vopmask mask, vfloat x, vfloat y) {
  return vbslq_f32(mask, x, y);
}
static INLINE vdouble vsel_vd_vo_vd_vd(vopmask mask, vdouble x, vdouble y) {
  return vbslq_f64(vreinterpretq_u64_u32(mask), x, y);
}
static INLINE vint2 vcast_vi2_i(int i) { return vdupq_n_s32(i); }
static INLINE vint vcast_vi_i(int i) { return vdup_n_s32(i); }
static INLINE vfloat vcast_vf_vi2(vint2 vi) { return vcvtq_f32_s32(vi); }
static INLINE vdouble vcast_vd_vi(vint vi) { return vcvtq_f64_s64(vmovl_s32(vi)); }
static INLINE vint2 vrint_vi2_vf(vfloat d) {
  return vcvtq_s32_f32(vrndnq_f32(d));
}
static INLINE vdouble vrint_vd_vd(vdouble d) { return vrndnq_f64(d); }
static INLINE vint vrint_vi_vd(vdouble d) {
  return vqmovn_s64(vcvtq_s64_f64(vrndnq_f64(d)));
}
static INLINE vint2 vadd_vi2_vi2_vi2(vint2 x, vint2 y) { return vaddq_s32(x, y); }
static INLINE vint2 vsub_vi2_vi2_vi2(vint2 x, vint2 y) { return vsubq_s32(x, y); }
static INLINE vint2 vneg_vi2_vi2(vint2 e) { return vnegq_s32(e); }
static INLINE vint2 vand_vi2_vi2_vi2(vint2 x, vint2 y) { return vandq_s32(x, y); }
static INLINE vint2 vsel_vi2_vo_vi2_vi2(vopmask m, vint2 x, vint2 y) {
  return vbslq_s32(m, x, y);
}
static INLINE vint vadd_vi_vi_vi(vint x, vint y) { return vadd_s32(x, y); }
static INLINE vint vsub_vi_vi_vi(vint x, vint y) { return vsub_s32(x, y); }
static INLINE vint vneg_vi_vi(vint e) { return vneg_s32(e); }
static INLINE vint vand_vi_vi_vi(vint x, vint y) { return vand_s32(x, y); }
static INLINE vint vsel_vi_vo_vi_vi(vopmask m, vint x, vint y) {
  return vbsl_s32(vget_low_u32(m), x, y);
}
static INLINE vmask vcastu_vm_vi(vint vi) {
  return vrev64q_u32(vreinterpretq_u32_u64(vmovl_u32(vreinterpret_u32_s32(vi))));
}
static INLINE vint vcastu_vi_vm(vmask vm) {
  return vreinterpret_s32_u32(vmovn_u64(vreinterpretq_u64_u32(vrev64q_u32(vm))));
}
static INLINE vopmask vcast_vo32_vo64(vopmask m) { return vuzpq_u32(m, m).val[0]; }
static INLINE vopmask vor_vo_vo_vo(vopmask x, vopmask y) { return vorrq_u32(x, y); }
static INLINE vmask vandnot_vm_vo32_vm(vopmask x, vmask y) { return vbicq_u32(y, x); }
static INLINE vmask vandnot_vm_vo64_vm(vopmask x, vmask y) { return vbicq_u32(y, x); }
static INLINE vmask vadd64_vm_vm_vm(vmask x, vmask y) {
  return vreinterpretq_u32_s64(vaddq_s64(vreinterpretq_s64_u32(x),
    vreinterpretq_s64_u32(y)));
}
static INLINE vopmask visnan_vo_vd(vdouble d) {
  return vmvnq_u32(vreinterpretq_u32_u64(vceqq_f64(d, d)));
}
static INLINE vopmask vispinf_vo_vd(vdouble d) {
  return vreinterpretq_u32_u64(vceqq_f64(d, vdupq_n_f64(__builtin_inf())));
}
static INLINE vopmask visnan_vo_vf(vfloat d) { return vmvnq_u32(vceqq_f32(d, d)); }
static INLINE vopmask vispinf_vo_vf(vfloat d) {
  return vceqq_f32(d, vdupq_n_f32(__builtin_inff()));
}

typedef struct {
  vfloat x, y;
} vfloat2;

static INLINE CONST vfloat vf2getx_vf_vf2(vfloat2 v) { return v.x; }
static INLINE CONST vfloat vf2gety_vf_vf2(vfloat2 v) { return v.y; }
static INLINE CONST vfloat2 vf2setxy_vf2_vf_vf(vfloat x, vfloat y) {
  vfloat2 v;
  v.x = x;
  v.y = y;
  return v;
}
static INLINE CONST vfloat2 vcast_vf2_f_f(float h, float l) {
  return vf2setxy_vf2_vf_vf(vcast_vf_f(h), vcast_vf_f(l));
}
static INLINE CONST vfloat2 dfscale_vf2_vf2_vf(vfloat2 d, vfloat s) {
  return vf2setxy_vf2_vf_vf(vmul_vf_vf_vf(vf2getx_vf_vf2(d), s),
    vmul_vf_vf_vf(vf2gety_vf_vf2(d), s));
}
static INLINE CONST vfloat2 dfadd2_vf2_vf_vf(vfloat x, vfloat y) {
  vfloat s = vadd_vf_vf_vf(x, y);
  vfloat v = vsub_vf_vf_vf(s, x);
  return vf2setxy_vf2_vf_vf(s, vadd_vf_vf_vf(vsub_vf_vf_vf(x,
    vsub_vf_vf_vf(s, v)), vsub_vf_vf_vf(y, v)));
}
static INLINE CONST vfloat2 dfadd_vf2_vf2_vf(vfloat2 x, vfloat y) {
  vfloat s = vadd_vf_vf_vf(vf2getx_vf_vf2(x), y);
  return vf2setxy_vf2_vf_vf(s, vadd_vf_vf_vf(vadd_vf_vf_vf(
    vsub_vf_vf_vf(vf2getx_vf_vf2(x), s), y), vf2gety_vf_vf2(x)));
}
static INLINE CONST vfloat2 dfadd_vf2_vf2_vf2(vfloat2 x, vfloat2 y) {
  vfloat s = vadd_vf_vf_vf(vf2getx_vf_vf2(x), vf2getx_vf_vf2(y));
  return vf2setxy_vf2_vf_vf(s, vadd_vf_vf_vf(vadd_vf_vf_vf(vadd_vf_vf_vf(
    vsub_vf_vf_vf(vf2getx_vf_vf2(x), s), vf2getx_vf_vf2(y)),
    vf2gety_vf_vf2(x)), vf2gety_vf_vf2(y)));
}
static INLINE CONST vfloat2 dfdiv_vf2_vf2_vf2(vfloat2 n, vfloat2 d) {
  vfloat t = vrec_vf_vf(vf2getx_vf_vf2(d));
  vfloat s = vmul_vf_vf_vf(vf2getx_vf_vf2(n), t);
  vfloat u = vfmapn_vf_vf_vf_vf(t, vf2getx_vf_vf2(n), s);
  vfloat v = vfmanp_vf_vf_vf_vf(vf2gety_vf_vf2(d), t,
    vfmanp_vf_vf_vf_vf(vf2getx_vf_vf2(d), t, vcast_vf_f(1)));
  return vf2setxy_vf2_vf_vf(s, vfma_vf_vf_vf_vf(s, v,
    vfma_vf_vf_vf_vf(vf2gety_vf_vf2(n), t, u)));
}
static INLINE CONST vfloat2 dfmul_vf2_vf2_vf(vfloat2 x, vfloat y) {
  vfloat s = vmul_vf_vf_vf(vf2getx_vf_vf2(x), y);
  return vf2setxy_vf2_vf_vf(s, vfma_vf_vf_vf_vf(vf2gety_vf_vf2(x), y,
    vfmapn_vf_vf_vf_vf(vf2getx_vf_vf2(x), y, s)));
}
static INLINE CONST vint2 vilogb2k_vi2_vf(vfloat d) {
  vint2 q = vreinterpret_vi2_vf(d);
  q = vreinterpretq_s32_u32(vshrq_n_u32(vreinterpretq_u32_s32(q), 23));
  q = vand_vi2_vi2_vi2(q, vcast_vi2_i(0xff));
  q = vsub_vi2_vi2_vi2(q, vcast_vi2_i(0x7f));
  return q;
}
static INLINE CONST vfloat vpow2i_vf_vi2(vint2 q) {
  return vreinterpret_vf_vi2(vshlq_n_s32(vadd_vi2_vi2_vi2(q,
    vcast_vi2_i(0x7f)), 23));
}
static INLINE CONST vfloat vldexp2_vf_vf_vi2(vfloat d, vint2 e) {
  vint2 e2 = vshrq_n_s32(e, 1);
  return vmul_vf_vf_vf(vmul_vf_vf_vf(d, vpow2i_vf_vi2(e2)),
    vpow2i_vf_vi2(vsub_vi2_vi2_vi2(e, e2)));
}
static INLINE CONST vfloat vldexp3_vf_vf_vi2(vfloat d, vint2 q) {
  return vreinterpret_vf_vi2(vadd_vi2_vi2_vi2(vreinterpret_vi2_vf(d),
    vshlq_n_s32(q, 23)));
}

static INLINE CONST vint vilogb2k_vi_vd(vdouble d) {
  vint q = vcastu_vi_vm(vreinterpret_vm_vd(d));
  q = vreinterpret_s32_u32(vshr_n_u32(vreinterpret_u32_s32(q), 20));
  q = vand_vi_vi_vi(q, vcast_vi_i(0x7ff));
  q = vsub_vi_vi_vi(q, vcast_vi_i(0x3ff));
  return q;
}
static INLINE CONST vdouble vpow2i_vd_vi(vint q) {
  q = vadd_vi_vi_vi(vcast_vi_i(0x3ff), q);
  return vreinterpret_vd_vm(vcastu_vm_vi(vshl_n_s32(q, 20)));
}
static INLINE CONST vdouble vldexp2_vd_vd_vi(vdouble d, vint e) {
  vint e2 = vshr_n_s32(e, 1);
  return vmul_vd_vd_vd(vmul_vd_vd_vd(d, vpow2i_vd_vi(e2)),
    vpow2i_vd_vi(vsub_vi_vi_vi(e, e2)));
}
static INLINE CONST vdouble vldexp3_vd_vd_vi(vdouble d, vint q) {
  return vreinterpret_vd_vm(vadd64_vm_vm_vm(vreinterpret_vm_vd(d),
    vcastu_vm_vi(vshl_n_s32(q, 20))));
}

CONST vdouble Sleef_logd2_u35advsimd(vdouble d) {
  vdouble x, x2;
  vdouble t, m;
  vopmask o = vlt_vo_vd_vd(d, vcast_vd_d(0x1p-1022));
  d = vsel_vd_vo_vd_vd(o, vmul_vd_vd_vd(d,
    vcast_vd_d((double)((1LL) << 32) * (double)((1LL) << 32))), d);
  vint e = vilogb2k_vi_vd(vmul_vd_vd_vd(d, vcast_vd_d(1.0/0.75)));
  m = vldexp3_vd_vd_vi(d, vneg_vi_vi(e));
  e = vsel_vi_vo_vi_vi(vcast_vo32_vo64(o),
    vsub_vi_vi_vi(e, vcast_vi_i(64)), e);
  x = vdiv_vd_vd_vd(vsub_vd_vd_vd(m, vcast_vd_d(1)),
    vadd_vd_vd_vd(vcast_vd_d(1), m));
  x2 = vmul_vd_vd_vd(x, x);
  vdouble x4 = vmul_vd_vd_vd(x2, x2);
  vdouble x8 = vmul_vd_vd_vd(x4, x4);
  vdouble x3 = vmul_vd_vd_vd(x, x2);
  t = vmla_vd_vd_vd_vd((x8), (vmla_vd_vd_vd_vd((x4),
    (vcast_vd_d(0.153487338491425068243146)), (vmla_vd_vd_vd_vd((x2),
    (vcast_vd_d(0.152519917006351951593857)),
    (vcast_vd_d(0.181863266251982985677316)))))), (vmla_vd_vd_vd_vd((x4),
    (vmla_vd_vd_vd_vd((x2), (vcast_vd_d(0.222221366518767365905163)),
    (vcast_vd_d(0.285714294746548025383248)))), (vmla_vd_vd_vd_vd((x2),
    (vcast_vd_d(0.399999999950799600689777)),
    (vcast_vd_d(0.6666666666667778740063)))))));
  x = vmla_vd_vd_vd_vd(x, vcast_vd_d(2), vmul_vd_vd_vd(
    vcast_vd_d(0.693147180559945286226764), vcast_vd_vi(e)));
  x = vmla_vd_vd_vd_vd(x3, t, x);
  x = vsel_vd_vo_vd_vd(vispinf_vo_vd(d), vcast_vd_d(__builtin_inf()), x);
  x = vsel_vd_vo_vd_vd(vor_vo_vo_vo(vlt_vo_vd_vd(d, vcast_vd_d(0)),
    visnan_vo_vd(d)), vcast_vd_d(__builtin_nan("")), x);
  x = vsel_vd_vo_vd_vd(veq_vo_vd_vd(d, vcast_vd_d(0)),
    vcast_vd_d(-__builtin_inf()), x);
  return x;
}

CONST vdouble Sleef_expd2_u10advsimd(vdouble d) {
  vdouble u = vrint_vd_vd(vmul_vd_vd_vd(d,
    vcast_vd_d(1.442695040888963407359924681001892137426645954152985934135449406931)));
  vdouble s;
  vint q = vrint_vi_vd(u);
  s = vmla_vd_vd_vd_vd(u, vcast_vd_d(-.69314718055966295651160180568695068359375), d);
  s = vmla_vd_vd_vd_vd(u,
    vcast_vd_d(-.28235290563031577122588448175013436025525412068e-12), s);
  vdouble s2 = vmul_vd_vd_vd(s, s);
  vdouble s4 = vmul_vd_vd_vd(s2, s2);
  vdouble s8 = vmul_vd_vd_vd(s4, s4);
  u = vmla_vd_vd_vd_vd((s8), (vmla_vd_vd_vd_vd((s),
    (vcast_vd_d(+0.2081276378237164457e-8)),
    (vcast_vd_d(+0.2511210703042288022e-7)))), (vmla_vd_vd_vd_vd((s4),
    (vmla_vd_vd_vd_vd((s2), (vmla_vd_vd_vd_vd((s),
    (vcast_vd_d(+0.2755762628169491192e-6)),
    (vcast_vd_d(+0.2755723402025388239e-5)))), (vmla_vd_vd_vd_vd((s),
    (vcast_vd_d(+0.2480158687479686264e-4)),
    (vcast_vd_d(+0.1984126989855865850e-3)))))), (vmla_vd_vd_vd_vd((s2),
    (vmla_vd_vd_vd_vd((s), (vcast_vd_d(+0.1388888888914497797e-2)),
    (vcast_vd_d(+0.8333333333314938210e-2)))), (vmla_vd_vd_vd_vd((s),
    (vcast_vd_d(+0.4166666666666602598e-1)),
    (vcast_vd_d(+0.1666666666666669072e+0)))))))));
  u = vfma_vd_vd_vd_vd(u, s, vcast_vd_d(+0.5000000000000000000e+0));
  u = vfma_vd_vd_vd_vd(u, s, vcast_vd_d(+0.1000000000000000000e+1));
  u = vfma_vd_vd_vd_vd(u, s, vcast_vd_d(+0.1000000000000000000e+1));
  u = vldexp2_vd_vd_vi(u, q);
  vopmask o = vgt_vo_vd_vd(d, vcast_vd_d(0x1.62e42fefa39efp+9));
  u = vsel_vd_vo_vd_vd(o, vcast_vd_d(__builtin_inf()), u);
  u = vreinterpret_vd_vm(vandnot_vm_vo64_vm(vlt_vo_vd_vd(d,
    vcast_vd_d(-1000)), vreinterpret_vm_vd(u)));
  return u;
}

CONST vfloat Sleef_expf4_u10advsimd(vfloat d) {
  vint2 q = vrint_vi2_vf(vmul_vf_vf_vf(d,
    vcast_vf_f(1.442695040888963407359924681001892137426645954152985934135449406931f)));
  vfloat s, u;
  s = vmla_vf_vf_vf_vf(vcast_vf_vi2(q), vcast_vf_f(-0.693145751953125f), d);
  s = vmla_vf_vf_vf_vf(vcast_vf_vi2(q), vcast_vf_f(-1.428606765330187045e-06f), s);
  u = vcast_vf_f(0.000198527617612853646278381);
  u = vmla_vf_vf_vf_vf(u, s, vcast_vf_f(0.00139304355252534151077271));
  u = vmla_vf_vf_vf_vf(u, s, vcast_vf_f(0.00833336077630519866943359));
  u = vmla_vf_vf_vf_vf(u, s, vcast_vf_f(0.0416664853692054748535156));
  u = vmla_vf_vf_vf_vf(u, s, vcast_vf_f(0.166666671633720397949219));
  u = vmla_vf_vf_vf_vf(u, s, vcast_vf_f(0.5));
  u = vadd_vf_vf_vf(vcast_vf_f(1.0f), vmla_vf_vf_vf_vf(vmul_vf_vf_vf(s, s), u, s));
  u = vldexp2_vf_vf_vi2(u, q);
  u = vreinterpret_vf_vm(vandnot_vm_vo32_vm(vlt_vo_vf_vf(d, vcast_vf_f(-104)),
    vreinterpret_vm_vf(u)));
  u = vsel_vf_vo_vf_vf(vlt_vo_vf_vf(vcast_vf_f(100), d),
    vcast_vf_f(__builtin_inff()), u);
  return u;
}

CONST vfloat Sleef_logf4_u10advsimd(vfloat d) {
  vfloat2 x;
  vfloat t, m, x2;
  vopmask o = vlt_vo_vf_vf(d, vcast_vf_f(0x1p-126));
  d = vsel_vf_vo_vf_vf(o, vmul_vf_vf_vf(d,
    vcast_vf_f((float)((1LL) << 32) * (float)((1LL) << 32))), d);
  vint2 e = vilogb2k_vi2_vf(vmul_vf_vf_vf(d, vcast_vf_f(1.0f/0.75f)));
  m = vldexp3_vf_vf_vi2(d, vneg_vi2_vi2(e));
  e = vsel_vi2_vo_vi2_vi2(o, vsub_vi2_vi2_vi2(e, vcast_vi2_i(64)), e);
  vfloat2 s = dfmul_vf2_vf2_vf(vcast_vf2_f_f(0.69314718246459960938f,
    -1.904654323148236017e-09f), vcast_vf_vi2(e));
  x = dfdiv_vf2_vf2_vf2(dfadd2_vf2_vf_vf(vcast_vf_f(-1), m),
    dfadd2_vf2_vf_vf(vcast_vf_f(1), m));
  x2 = vmul_vf_vf_vf(vf2getx_vf_vf2(x), vf2getx_vf_vf2(x));
  t = vcast_vf_f(+0.3027294874e+0f);
  t = vmla_vf_vf_vf_vf(t, x2, vcast_vf_f(+0.3996108174e+0f));
  t = vmla_vf_vf_vf_vf(t, x2, vcast_vf_f(+0.6666694880e+0f));
  s = dfadd_vf2_vf2_vf2(s, dfscale_vf2_vf2_vf(x, vcast_vf_f(2)));
  s = dfadd_vf2_vf2_vf(s, vmul_vf_vf_vf(vmul_vf_vf_vf(x2, vf2getx_vf_vf2(x)), t));
  vfloat r = vadd_vf_vf_vf(vf2getx_vf_vf2(s), vf2gety_vf_vf2(s));
  r = vsel_vf_vo_vf_vf(vispinf_vo_vf(d), vcast_vf_f(__builtin_inff()), r);
  r = vsel_vf_vo_vf_vf(vor_vo_vo_vo(vlt_vo_vf_vf(d, vcast_vf_f(0)),
    visnan_vo_vf(d)), vcast_vf_f(__builtin_nanf("")), r);
  r = vsel_vf_vo_vf_vf(veq_vo_vf_vf(d, vcast_vf_f(0)),
    vcast_vf_f(-__builtin_inff()), r);
  return r;
}

#endif
