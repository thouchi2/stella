#ifndef STELLA_GEN_H
#define STELLA_GEN_H

/**********************************************************
DO NOT EDIT!
I AM AUTO GENERATED FROM utils/electric/gen/gen_metrics.py.
I WAS ADDED TO VERSION CONTROL AT LUKE'S REQUEST.
***********************************************************/

#include <float.h>
#define FDIV(num, den) (fabs(den) < DBL_MIN ? 0.0 : num / den)
#define INTERIOR_METRICS_3D {\
enum {c=0,xh=1,yh=2,zh=3};\
PetscScalar ***vec__x_r;\
double x_r[4];\
vec__x_r = jac_v[met->t3map[0][0]];\
PetscScalar ***vec__x_s;\
double x_s[4];\
vec__x_s = jac_v[met->t3map[0][1]];\
PetscScalar ***vec__x_t;\
double x_t[4];\
vec__x_t = jac_v[met->t3map[0][2]];\
PetscScalar ***vec__y_r;\
double y_r[4];\
vec__y_r = jac_v[met->t3map[1][0]];\
PetscScalar ***vec__y_s;\
double y_s[4];\
vec__y_s = jac_v[met->t3map[1][1]];\
PetscScalar ***vec__y_t;\
double y_t[4];\
vec__y_t = jac_v[met->t3map[1][2]];\
PetscScalar ***vec__z_r;\
double z_r[4];\
vec__z_r = jac_v[met->t3map[2][0]];\
PetscScalar ***vec__z_s;\
double z_s[4];\
vec__z_s = jac_v[met->t3map[2][1]];\
PetscScalar ***vec__z_t;\
double z_t[4];\
vec__z_t = jac_v[met->t3map[2][2]];\
double g_11;\
double g11;\
double g_12;\
double g12;\
double g_13;\
double g13;\
double g_22;\
double g22;\
double g_23;\
double g23;\
double g_33;\
double g33;\
double g;\
double pt0[3], pt1[3], pt2[3];\
\
for (k = kbeg; k < kend; k++)\
{\
    for (j = jbeg; j < jend; j++)\
    {\
        for (i = ibeg; i < iend; i++)\
        {\
            x_r[c] = si[-1] * coors[k][j][i-1].x + si[1] * coors[k][j][i+1].x;\
            y_r[c] = si[-1] * coors[k][j][i-1].y + si[1] * coors[k][j][i+1].y;\
            z_r[c] = si[-1] * coors[k][j][i-1].z + si[1] * coors[k][j][i+1].z;\
            x_s[c] = si[-1] * coors[k][j-1][i].x + si[1] * coors[k][j+1][i].x;\
            y_s[c] = si[-1] * coors[k][j-1][i].y + si[1] * coors[k][j+1][i].y;\
            z_s[c] = si[-1] * coors[k][j-1][i].z + si[1] * coors[k][j+1][i].z;\
            x_t[c] = si[-1] * coors[k-1][j][i].x + si[1] * coors[k+1][j][i].x;\
            y_t[c] = si[-1] * coors[k-1][j][i].y + si[1] * coors[k+1][j][i].y;\
            z_t[c] = si[-1] * coors[k-1][j][i].z + si[1] * coors[k+1][j][i].z;\
            g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
            g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
            g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
            g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
            g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
            g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            jac[k][j][i] = sqrt(g);\
            vec__x_r[k][j][i] = x_r[c];\
            vec__y_r[k][j][i] = y_r[c];\
            vec__z_r[k][j][i] = z_r[c];\
            vec__x_s[k][j][i] = x_s[c];\
            vec__y_s[k][j][i] = y_s[c];\
            vec__z_s[k][j][i] = z_s[c];\
            vec__x_t[k][j][i] = x_t[c];\
            vec__y_t[k][j][i] = y_t[c];\
            vec__z_t[k][j][i] = z_t[c];\
            g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
            g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
            g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
            coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
            coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
            coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
            coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
            coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
            coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
            coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
            coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
            coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
            coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
            coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
            coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
\
        }\
\
    }\
\
}\
\
\
for (k = kbeg; k < kend; k++)\
{\
    for (j = jbeg; j < jend; j++)\
    {\
        for (i = ibeg; i < xs + xm; i++)\
        {\
            x_r[xh] = coors[k][j][i].x - coors[k][j][i-1].x;\
            y_r[xh] = coors[k][j][i].y - coors[k][j][i-1].y;\
            z_r[xh] = coors[k][j][i].z - coors[k][j][i-1].z;\
\
            pt0[0] = .5 * coors[k][j-1][i-1].x + .5 * coors[k][j-1][i].x;\
            pt1[0] = .5 * coors[k][j+1][i-1].x + .5 * coors[k][j+1][i].x;\
            pt0[1] = .5 * coors[k][j-1][i-1].y + .5 * coors[k][j-1][i].y;\
            pt1[1] = .5 * coors[k][j+1][i-1].y + .5 * coors[k][j+1][i].y;\
            pt0[2] = .5 * coors[k][j-1][i-1].z + .5 * coors[k][j-1][i].z;\
            pt1[2] = .5 * coors[k][j+1][i-1].z + .5 * coors[k][j+1][i].z;\
            x_s[xh] = (pt1[0] - pt0[0]) * .5;\
            y_s[xh] = (pt1[1] - pt0[1]) * .5;\
            z_s[xh] = (pt1[2] - pt0[2]) * .5;\
\
            pt0[0] = .5 * coors[k-1][j][i-1].x + .5 * coors[k-1][j][i].x;\
            pt1[0] = .5 * coors[k+1][j][i-1].x + .5 * coors[k+1][j][i].x;\
            pt0[1] = .5 * coors[k-1][j][i-1].y + .5 * coors[k-1][j][i].y;\
            pt1[1] = .5 * coors[k+1][j][i-1].y + .5 * coors[k+1][j][i].y;\
            pt0[2] = .5 * coors[k-1][j][i-1].z + .5 * coors[k-1][j][i].z;\
            pt1[2] = .5 * coors[k+1][j][i-1].z + .5 * coors[k+1][j][i].z;\
            x_t[xh] = (pt1[0] - pt0[0]) * .5;\
            y_t[xh] = (pt1[1] - pt0[1]) * .5;\
            z_t[xh] = (pt1[2] - pt0[2]) * .5;\
\
            g_11 = x_r[xh]*x_r[xh] + y_r[xh]*y_r[xh] + z_r[xh]*z_r[xh];\
            g_12 = x_r[xh]*x_s[xh] + y_r[xh]*y_s[xh] + z_r[xh]*z_s[xh];\
            g_13 = x_r[xh]*x_t[xh] + y_r[xh]*y_t[xh] + z_r[xh]*z_t[xh];\
            g_22 = x_s[xh]*x_s[xh] + y_s[xh]*y_s[xh] + z_s[xh]*z_s[xh];\
            g_23 = x_s[xh]*x_t[xh] + y_s[xh]*y_t[xh] + z_s[xh]*z_t[xh];\
            g_33 = x_t[xh]*x_t[xh] + y_t[xh]*y_t[xh] + z_t[xh]*z_t[xh];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            g11 = FDIV(1.0,g) * (g_22*g_33 - g_23*g_23);\
            coef[MET_W][k][j][i] = sqrt(g) * g11;\
        }\
\
    }\
\
}\
\
\
for (k = kbeg; k < kend; k++)\
{\
    for (j = jbeg; j < ys + ym; j++)\
    {\
        for (i = ibeg; i < iend; i++)\
        {\
            pt0[0] = .5 * coors[k][j-1][i-1].x + .5 * coors[k][j][i-1].x;\
            pt1[0] = .5 * coors[k][j-1][i+1].x + .5 * coors[k][j][i+1].x;\
            pt0[1] = .5 * coors[k][j-1][i-1].y + .5 * coors[k][j][i-1].y;\
            pt1[1] = .5 * coors[k][j-1][i+1].y + .5 * coors[k][j][i+1].y;\
            pt0[2] = .5 * coors[k][j-1][i-1].z + .5 * coors[k][j][i-1].z;\
            pt1[2] = .5 * coors[k][j-1][i+1].z + .5 * coors[k][j][i+1].z;\
            x_r[yh] = (pt1[0] - pt0[0]) * .5;\
            y_r[yh] = (pt1[1] - pt0[1]) * .5;\
            z_r[yh] = (pt1[2] - pt0[2]) * .5;\
\
            x_s[yh] = coors[k][j][i].x - coors[k][j-1][i].x;\
            y_s[yh] = coors[k][j][i].y - coors[k][j-1][i].y;\
            z_s[yh] = coors[k][j][i].z - coors[k][j-1][i].z;\
\
            pt0[0] = .5 * coors[k-1][j-1][i].x + .5 * coors[k-1][j][i].x;\
            pt1[0] = .5 * coors[k+1][j-1][i].x + .5 * coors[k+1][j][i].x;\
            pt0[1] = .5 * coors[k-1][j-1][i].y + .5 * coors[k-1][j][i].y;\
            pt1[1] = .5 * coors[k+1][j-1][i].y + .5 * coors[k+1][j][i].y;\
            pt0[2] = .5 * coors[k-1][j-1][i].z + .5 * coors[k-1][j][i].z;\
            pt1[2] = .5 * coors[k+1][j-1][i].z + .5 * coors[k+1][j][i].z;\
            x_t[yh] = (pt1[0] - pt0[0]) * .5;\
            y_t[yh] = (pt1[1] - pt0[1]) * .5;\
            z_t[yh] = (pt1[2] - pt0[2]) * .5;\
\
            g_11 = x_r[yh]*x_r[yh] + y_r[yh]*y_r[yh] + z_r[yh]*z_r[yh];\
            g_12 = x_r[yh]*x_s[yh] + y_r[yh]*y_s[yh] + z_r[yh]*z_s[yh];\
            g_13 = x_r[yh]*x_t[yh] + y_r[yh]*y_t[yh] + z_r[yh]*z_t[yh];\
            g_22 = x_s[yh]*x_s[yh] + y_s[yh]*y_s[yh] + z_s[yh]*z_s[yh];\
            g_23 = x_s[yh]*x_t[yh] + y_s[yh]*y_t[yh] + z_s[yh]*z_t[yh];\
            g_33 = x_t[yh]*x_t[yh] + y_t[yh]*y_t[yh] + z_t[yh]*z_t[yh];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            g22 = FDIV(1.0,g) * (g_11*g_33 - g_13*g_13);\
            coef[MET_S][k][j][i] = sqrt(g) * g22;\
        }\
\
    }\
\
}\
\
\
for (k = kbeg; k < zs + zm; k++)\
{\
    for (j = jbeg; j < jend; j++)\
    {\
        for (i = ibeg; i < iend; i++)\
        {\
            pt0[0] = .5 * coors[k-1][j][i-1].x + .5 * coors[k][j][i-1].x;\
            pt1[0] = .5 * coors[k-1][j][i+1].x + .5 * coors[k][j][i+1].x;\
            pt0[1] = .5 * coors[k-1][j][i-1].y + .5 * coors[k][j][i-1].y;\
            pt1[1] = .5 * coors[k-1][j][i+1].y + .5 * coors[k][j][i+1].y;\
            pt0[2] = .5 * coors[k-1][j][i-1].z + .5 * coors[k][j][i-1].z;\
            pt1[2] = .5 * coors[k-1][j][i+1].z + .5 * coors[k][j][i+1].z;\
            x_r[zh] = (pt1[0] - pt0[0]) * .5;\
            y_r[zh] = (pt1[1] - pt0[1]) * .5;\
            z_r[zh] = (pt1[2] - pt0[2]) * .5;\
\
            pt0[0] = .5 * coors[k-1][j-1][i].x + .5 * coors[k][j-1][i].x;\
            pt1[0] = .5 * coors[k-1][j+1][i].x + .5 * coors[k][j+1][i].x;\
            pt0[1] = .5 * coors[k-1][j-1][i].y + .5 * coors[k][j-1][i].y;\
            pt1[1] = .5 * coors[k-1][j+1][i].y + .5 * coors[k][j+1][i].y;\
            pt0[2] = .5 * coors[k-1][j-1][i].z + .5 * coors[k][j-1][i].z;\
            pt1[2] = .5 * coors[k-1][j+1][i].z + .5 * coors[k][j+1][i].z;\
            x_s[zh] = (pt1[0] - pt0[0]) * .5;\
            y_s[zh] = (pt1[1] - pt0[1]) * .5;\
            z_s[zh] = (pt1[2] - pt0[2]) * .5;\
\
            x_t[zh] = coors[k][j][i].x - coors[k-1][j][i].x;\
            y_t[zh] = coors[k][j][i].y - coors[k-1][j][i].y;\
            z_t[zh] = coors[k][j][i].z - coors[k-1][j][i].z;\
\
            g_11 = x_r[zh]*x_r[zh] + y_r[zh]*y_r[zh] + z_r[zh]*z_r[zh];\
            g_12 = x_r[zh]*x_s[zh] + y_r[zh]*y_s[zh] + z_r[zh]*z_s[zh];\
            g_13 = x_r[zh]*x_t[zh] + y_r[zh]*y_t[zh] + z_r[zh]*z_t[zh];\
            g_22 = x_s[zh]*x_s[zh] + y_s[zh]*y_s[zh] + z_s[zh]*z_s[zh];\
            g_23 = x_s[zh]*x_t[zh] + y_s[zh]*y_t[zh] + z_s[zh]*z_t[zh];\
            g_33 = x_t[zh]*x_t[zh] + y_t[zh]*y_t[zh] + z_t[zh]*z_t[zh];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            g33 = FDIV(1.0,g) * (g_11*g_22 - g_12*g_12);\
            coef[MET_B][k][j][i] = sqrt(g) * g33;\
        }\
\
    }\
\
}\
\
\
}\


#define BOUNDARY_METRICS_3D {\
enum {c=0,xh=1,yh=2,zh=3};\
PetscScalar ***vec__x_r;\
double x_r[4];\
vec__x_r = jac_v[met->t3map[0][0]];\
PetscScalar ***vec__x_s;\
double x_s[4];\
vec__x_s = jac_v[met->t3map[0][1]];\
PetscScalar ***vec__x_t;\
double x_t[4];\
vec__x_t = jac_v[met->t3map[0][2]];\
PetscScalar ***vec__y_r;\
double y_r[4];\
vec__y_r = jac_v[met->t3map[1][0]];\
PetscScalar ***vec__y_s;\
double y_s[4];\
vec__y_s = jac_v[met->t3map[1][1]];\
PetscScalar ***vec__y_t;\
double y_t[4];\
vec__y_t = jac_v[met->t3map[1][2]];\
PetscScalar ***vec__z_r;\
double z_r[4];\
vec__z_r = jac_v[met->t3map[2][0]];\
PetscScalar ***vec__z_s;\
double z_s[4];\
vec__z_s = jac_v[met->t3map[2][1]];\
PetscScalar ***vec__z_t;\
double z_t[4];\
vec__z_t = jac_v[met->t3map[2][2]];\
double g_11;\
double g11;\
double g_12;\
double g12;\
double g_13;\
double g13;\
double g_22;\
double g22;\
double g_23;\
double g23;\
double g_33;\
double g33;\
double g;\
double pt0[3], pt1[3], pt2[3];\
ibeg = xs; iend = xs + xm; jbeg = ys; jend = ys + ym; kbeg = zs; kend = zs + zm;\
if (xs == 0 && !per_x) ibeg++;\
if (xs + xm == ngx && !per_x) iend--;\
if (ys == 0 && !per_y) jbeg++;\
if (ys + ym == ngy && !per_y) jend--;	\
if (zs == 0 && !per_z) kbeg++;\
if ((zs + zm == ngz && !per_z) && !per_z) kend--;\
if (xs + xm == ngx && !per_x) {\
    i = ngx - 1;\
    for (k = kbeg; k < kend; k++)\
    {\
        for (j = jbeg; j < ys + ym; j++)\
        {\
            pt0[0] = .5 * coors[k][j-1][i].x + .5 * coors[k][j][i].x;\
            pt1[0] = .5 * coors[k][j-1][i-1].x + .5 * coors[k][j][i-1].x;\
            pt2[0] = .5 * coors[k][j-1][i-2].x + .5 * coors[k][j][i-2].x;\
            pt0[1] = .5 * coors[k][j-1][i].y + .5 * coors[k][j][i].y;\
            pt1[1] = .5 * coors[k][j-1][i-1].y + .5 * coors[k][j][i-1].y;\
            pt2[1] = .5 * coors[k][j-1][i-2].y + .5 * coors[k][j][i-2].y;\
            pt0[2] = .5 * coors[k][j-1][i].z + .5 * coors[k][j][i].z;\
            pt1[2] = .5 * coors[k][j-1][i-1].z + .5 * coors[k][j][i-1].z;\
            pt2[2] = .5 * coors[k][j-1][i-2].z + .5 * coors[k][j][i-2].z;\
            x_r[yh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
            y_r[yh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
            z_r[yh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
            x_s[yh] = coors[k][j][i].x - coors[k][j-1][i].x;\
            y_s[yh] = coors[k][j][i].y - coors[k][j-1][i].y;\
            z_s[yh] = coors[k][j][i].z - coors[k][j-1][i].z;\
            pt0[0] = .5 * coors[k-1][j-1][i].x + .5 * coors[k-1][j][i].x;\
            pt1[0] = .5 * coors[k+1][j-1][i].x + .5 * coors[k+1][j][i].x;\
            pt0[1] = .5 * coors[k-1][j-1][i].y + .5 * coors[k-1][j][i].y;\
            pt1[1] = .5 * coors[k+1][j-1][i].y + .5 * coors[k+1][j][i].y;\
            pt0[2] = .5 * coors[k-1][j-1][i].z + .5 * coors[k-1][j][i].z;\
            pt1[2] = .5 * coors[k+1][j-1][i].z + .5 * coors[k+1][j][i].z;\
            x_t[yh] = (pt1[0] - pt0[0]) * .5;\
            y_t[yh] = (pt1[1] - pt0[1]) * .5;\
            z_t[yh] = (pt1[2] - pt0[2]) * .5;\
            g_11 = x_r[yh]*x_r[yh] + y_r[yh]*y_r[yh] + z_r[yh]*z_r[yh];\
            g_12 = x_r[yh]*x_s[yh] + y_r[yh]*y_s[yh] + z_r[yh]*z_s[yh];\
            g_13 = x_r[yh]*x_t[yh] + y_r[yh]*y_t[yh] + z_r[yh]*z_t[yh];\
            g_22 = x_s[yh]*x_s[yh] + y_s[yh]*y_s[yh] + z_s[yh]*z_s[yh];\
            g_23 = x_s[yh]*x_t[yh] + y_s[yh]*y_t[yh] + z_s[yh]*z_t[yh];\
            g_33 = x_t[yh]*x_t[yh] + y_t[yh]*y_t[yh] + z_t[yh]*z_t[yh];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            g22 = FDIV(1.0,g) * (g_11*g_33 - g_13*g_13);\
            coef[MET_S][k][j][i] = sqrt(g) * g22;\
\
        }\
\
    }\
    for (k = kbeg; k < zs + zm; k++)\
    {\
        for (j = jbeg; j < jend; j++)\
        {\
            pt0[0] = .5 * coors[k-1][j][i].x + .5 * coors[k][j][i].x;\
            pt1[0] = .5 * coors[k-1][j][i-1].x + .5 * coors[k][j][i-1].x;\
            pt2[0] = .5 * coors[k-1][j][i-2].x + .5 * coors[k][j][i-2].x;\
            pt0[1] = .5 * coors[k-1][j][i].y + .5 * coors[k][j][i].y;\
            pt1[1] = .5 * coors[k-1][j][i-1].y + .5 * coors[k][j][i-1].y;\
            pt2[1] = .5 * coors[k-1][j][i-2].y + .5 * coors[k][j][i-2].y;\
            pt0[2] = .5 * coors[k-1][j][i].z + .5 * coors[k][j][i].z;\
            pt1[2] = .5 * coors[k-1][j][i-1].z + .5 * coors[k][j][i-1].z;\
            pt2[2] = .5 * coors[k-1][j][i-2].z + .5 * coors[k][j][i-2].z;\
            x_r[zh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
            y_r[zh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
            z_r[zh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
            pt0[0] = .5 * coors[k-1][j-1][i].x + .5 * coors[k][j-1][i].x;\
            pt1[0] = .5 * coors[k-1][j+1][i].x + .5 * coors[k][j+1][i].x;\
            pt0[1] = .5 * coors[k-1][j-1][i].y + .5 * coors[k][j-1][i].y;\
            pt1[1] = .5 * coors[k-1][j+1][i].y + .5 * coors[k][j+1][i].y;\
            pt0[2] = .5 * coors[k-1][j-1][i].z + .5 * coors[k][j-1][i].z;\
            pt1[2] = .5 * coors[k-1][j+1][i].z + .5 * coors[k][j+1][i].z;\
            x_s[zh] = (pt1[0] - pt0[0]) * .5;\
            y_s[zh] = (pt1[1] - pt0[1]) * .5;\
            z_s[zh] = (pt1[2] - pt0[2]) * .5;\
            x_t[zh] = coors[k][j][i].x - coors[k-1][j][i].x;\
            y_t[zh] = coors[k][j][i].y - coors[k-1][j][i].y;\
            z_t[zh] = coors[k][j][i].z - coors[k-1][j][i].z;\
            g_11 = x_r[zh]*x_r[zh] + y_r[zh]*y_r[zh] + z_r[zh]*z_r[zh];\
            g_12 = x_r[zh]*x_s[zh] + y_r[zh]*y_s[zh] + z_r[zh]*z_s[zh];\
            g_13 = x_r[zh]*x_t[zh] + y_r[zh]*y_t[zh] + z_r[zh]*z_t[zh];\
            g_22 = x_s[zh]*x_s[zh] + y_s[zh]*y_s[zh] + z_s[zh]*z_s[zh];\
            g_23 = x_s[zh]*x_t[zh] + y_s[zh]*y_t[zh] + z_s[zh]*z_t[zh];\
            g_33 = x_t[zh]*x_t[zh] + y_t[zh]*y_t[zh] + z_t[zh]*z_t[zh];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            g33 = FDIV(1.0,g) * (g_11*g_22 - g_12*g_12);\
            coef[MET_B][k][j][i] = sqrt(g) * g33;\
\
        }\
\
    }\
    for (k = kbeg; k < kend; k++)\
    {\
        for (j = jbeg; j < jend; j++)\
        {\
            x_r[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j][i-1].x + sp[-2] * coors[k][j][i-2].x;\
            y_r[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j][i-1].y + sp[-2] * coors[k][j][i-2].y;\
            z_r[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j][i-1].z + sp[-2] * coors[k][j][i-2].z;\
            x_s[c] = si[-1] * coors[k][j-1][i].x + si[1] * coors[k][j+1][i].x;\
            y_s[c] = si[-1] * coors[k][j-1][i].y + si[1] * coors[k][j+1][i].y;\
            z_s[c] = si[-1] * coors[k][j-1][i].z + si[1] * coors[k][j+1][i].z;\
            x_t[c] = si[-1] * coors[k-1][j][i].x + si[1] * coors[k+1][j][i].x;\
            y_t[c] = si[-1] * coors[k-1][j][i].y + si[1] * coors[k+1][j][i].y;\
            z_t[c] = si[-1] * coors[k-1][j][i].z + si[1] * coors[k+1][j][i].z;\
            g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
            g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
            g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
            g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
            g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
            g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            jac[k][j][i] = sqrt(g);\
            vec__x_r[k][j][i] = x_r[c];\
            vec__y_r[k][j][i] = y_r[c];\
            vec__z_r[k][j][i] = z_r[c];\
            vec__x_s[k][j][i] = x_s[c];\
            vec__y_s[k][j][i] = y_s[c];\
            vec__z_s[k][j][i] = z_s[c];\
            vec__x_t[k][j][i] = x_t[c];\
            vec__y_t[k][j][i] = y_t[c];\
            vec__z_t[k][j][i] = z_t[c];\
            g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
            g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
            g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
            coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
            coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
            coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
            coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
            coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
            coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
            coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
            coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
            coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
            coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
\
        }\
\
    }\
\
}\
\
\
if (xs == 0 && !per_x) {\
    i = 0;\
    for (k = kbeg; k < kend; k++)\
    {\
        for (j = jbeg; j < ys + ym; j++)\
        {\
            pt0[0] = .5 * coors[k][j-1][i].x + .5 * coors[k][j][i].x;\
            pt1[0] = .5 * coors[k][j-1][i+1].x + .5 * coors[k][j][i+1].x;\
            pt2[0] = .5 * coors[k][j-1][i+2].x + .5 * coors[k][j][i+2].x;\
            pt0[1] = .5 * coors[k][j-1][i].y + .5 * coors[k][j][i].y;\
            pt1[1] = .5 * coors[k][j-1][i+1].y + .5 * coors[k][j][i+1].y;\
            pt2[1] = .5 * coors[k][j-1][i+2].y + .5 * coors[k][j][i+2].y;\
            pt0[2] = .5 * coors[k][j-1][i].z + .5 * coors[k][j][i].z;\
            pt1[2] = .5 * coors[k][j-1][i+1].z + .5 * coors[k][j][i+1].z;\
            pt2[2] = .5 * coors[k][j-1][i+2].z + .5 * coors[k][j][i+2].z;\
            x_r[yh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
            y_r[yh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
            z_r[yh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
            x_s[yh] = coors[k][j][i].x - coors[k][j-1][i].x;\
            y_s[yh] = coors[k][j][i].y - coors[k][j-1][i].y;\
            z_s[yh] = coors[k][j][i].z - coors[k][j-1][i].z;\
            pt0[0] = .5 * coors[k-1][j-1][i].x + .5 * coors[k-1][j][i].x;\
            pt1[0] = .5 * coors[k+1][j-1][i].x + .5 * coors[k+1][j][i].x;\
            pt0[1] = .5 * coors[k-1][j-1][i].y + .5 * coors[k-1][j][i].y;\
            pt1[1] = .5 * coors[k+1][j-1][i].y + .5 * coors[k+1][j][i].y;\
            pt0[2] = .5 * coors[k-1][j-1][i].z + .5 * coors[k-1][j][i].z;\
            pt1[2] = .5 * coors[k+1][j-1][i].z + .5 * coors[k+1][j][i].z;\
            x_t[yh] = (pt1[0] - pt0[0]) * .5;\
            y_t[yh] = (pt1[1] - pt0[1]) * .5;\
            z_t[yh] = (pt1[2] - pt0[2]) * .5;\
            g_11 = x_r[yh]*x_r[yh] + y_r[yh]*y_r[yh] + z_r[yh]*z_r[yh];\
            g_12 = x_r[yh]*x_s[yh] + y_r[yh]*y_s[yh] + z_r[yh]*z_s[yh];\
            g_13 = x_r[yh]*x_t[yh] + y_r[yh]*y_t[yh] + z_r[yh]*z_t[yh];\
            g_22 = x_s[yh]*x_s[yh] + y_s[yh]*y_s[yh] + z_s[yh]*z_s[yh];\
            g_23 = x_s[yh]*x_t[yh] + y_s[yh]*y_t[yh] + z_s[yh]*z_t[yh];\
            g_33 = x_t[yh]*x_t[yh] + y_t[yh]*y_t[yh] + z_t[yh]*z_t[yh];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            g22 = FDIV(1.0,g) * (g_11*g_33 - g_13*g_13);\
            coef[MET_S][k][j][i] = sqrt(g) * g22;\
\
        }\
\
    }\
    for (k = kbeg; k < zs + zm; k++)\
    {\
        for (j = jbeg; j < jend; j++)\
        {\
            pt0[0] = .5 * coors[k-1][j][i].x + .5 * coors[k][j][i].x;\
            pt1[0] = .5 * coors[k-1][j][i+1].x + .5 * coors[k][j][i+1].x;\
            pt2[0] = .5 * coors[k-1][j][i+2].x + .5 * coors[k][j][i+2].x;\
            pt0[1] = .5 * coors[k-1][j][i].y + .5 * coors[k][j][i].y;\
            pt1[1] = .5 * coors[k-1][j][i+1].y + .5 * coors[k][j][i+1].y;\
            pt2[1] = .5 * coors[k-1][j][i+2].y + .5 * coors[k][j][i+2].y;\
            pt0[2] = .5 * coors[k-1][j][i].z + .5 * coors[k][j][i].z;\
            pt1[2] = .5 * coors[k-1][j][i+1].z + .5 * coors[k][j][i+1].z;\
            pt2[2] = .5 * coors[k-1][j][i+2].z + .5 * coors[k][j][i+2].z;\
            x_r[zh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
            y_r[zh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
            z_r[zh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
            pt0[0] = .5 * coors[k-1][j-1][i].x + .5 * coors[k][j-1][i].x;\
            pt1[0] = .5 * coors[k-1][j+1][i].x + .5 * coors[k][j+1][i].x;\
            pt0[1] = .5 * coors[k-1][j-1][i].y + .5 * coors[k][j-1][i].y;\
            pt1[1] = .5 * coors[k-1][j+1][i].y + .5 * coors[k][j+1][i].y;\
            pt0[2] = .5 * coors[k-1][j-1][i].z + .5 * coors[k][j-1][i].z;\
            pt1[2] = .5 * coors[k-1][j+1][i].z + .5 * coors[k][j+1][i].z;\
            x_s[zh] = (pt1[0] - pt0[0]) * .5;\
            y_s[zh] = (pt1[1] - pt0[1]) * .5;\
            z_s[zh] = (pt1[2] - pt0[2]) * .5;\
            x_t[zh] = coors[k][j][i].x - coors[k-1][j][i].x;\
            y_t[zh] = coors[k][j][i].y - coors[k-1][j][i].y;\
            z_t[zh] = coors[k][j][i].z - coors[k-1][j][i].z;\
            g_11 = x_r[zh]*x_r[zh] + y_r[zh]*y_r[zh] + z_r[zh]*z_r[zh];\
            g_12 = x_r[zh]*x_s[zh] + y_r[zh]*y_s[zh] + z_r[zh]*z_s[zh];\
            g_13 = x_r[zh]*x_t[zh] + y_r[zh]*y_t[zh] + z_r[zh]*z_t[zh];\
            g_22 = x_s[zh]*x_s[zh] + y_s[zh]*y_s[zh] + z_s[zh]*z_s[zh];\
            g_23 = x_s[zh]*x_t[zh] + y_s[zh]*y_t[zh] + z_s[zh]*z_t[zh];\
            g_33 = x_t[zh]*x_t[zh] + y_t[zh]*y_t[zh] + z_t[zh]*z_t[zh];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            g33 = FDIV(1.0,g) * (g_11*g_22 - g_12*g_12);\
            coef[MET_B][k][j][i] = sqrt(g) * g33;\
\
        }\
\
    }\
    for (k = kbeg; k < kend; k++)\
    {\
        for (j = jbeg; j < jend; j++)\
        {\
            x_r[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j][i+1].x + sm[2] * coors[k][j][i+2].x;\
            y_r[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j][i+1].y + sm[2] * coors[k][j][i+2].y;\
            z_r[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j][i+1].z + sm[2] * coors[k][j][i+2].z;\
            x_s[c] = si[-1] * coors[k][j-1][i].x + si[1] * coors[k][j+1][i].x;\
            y_s[c] = si[-1] * coors[k][j-1][i].y + si[1] * coors[k][j+1][i].y;\
            z_s[c] = si[-1] * coors[k][j-1][i].z + si[1] * coors[k][j+1][i].z;\
            x_t[c] = si[-1] * coors[k-1][j][i].x + si[1] * coors[k+1][j][i].x;\
            y_t[c] = si[-1] * coors[k-1][j][i].y + si[1] * coors[k+1][j][i].y;\
            z_t[c] = si[-1] * coors[k-1][j][i].z + si[1] * coors[k+1][j][i].z;\
            g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
            g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
            g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
            g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
            g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
            g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            jac[k][j][i] = sqrt(g);\
            vec__x_r[k][j][i] = x_r[c];\
            vec__y_r[k][j][i] = y_r[c];\
            vec__z_r[k][j][i] = z_r[c];\
            vec__x_s[k][j][i] = x_s[c];\
            vec__y_s[k][j][i] = y_s[c];\
            vec__z_s[k][j][i] = z_s[c];\
            vec__x_t[k][j][i] = x_t[c];\
            vec__y_t[k][j][i] = y_t[c];\
            vec__z_t[k][j][i] = z_t[c];\
            g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
            g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
            g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
            coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
            coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
            coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
            coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
            coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
            coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
            coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
            coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
            coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
            coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
\
        }\
\
    }\
\
}\
\
\
if (ys + ym == ngy && !per_y) {\
    j = ngy - 1;\
    for (k = kbeg; k < kend; k++)\
    {\
        for (i = ibeg; i < xs + xm; i++)\
        {\
            x_r[xh] = coors[k][j][i].x - coors[k][j][i-1].x;\
            y_r[xh] = coors[k][j][i].y - coors[k][j][i-1].y;\
            z_r[xh] = coors[k][j][i].z - coors[k][j][i-1].z;\
            pt0[0] = .5 * coors[k][j][i-1].x + .5 * coors[k][j][i].x;\
            pt1[0] = .5 * coors[k][j-1][i-1].x + .5 * coors[k][j-1][i].x;\
            pt2[0] = .5 * coors[k][j-2][i-1].x + .5 * coors[k][j-2][i].x;\
            pt0[1] = .5 * coors[k][j][i-1].y + .5 * coors[k][j][i].y;\
            pt1[1] = .5 * coors[k][j-1][i-1].y + .5 * coors[k][j-1][i].y;\
            pt2[1] = .5 * coors[k][j-2][i-1].y + .5 * coors[k][j-2][i].y;\
            pt0[2] = .5 * coors[k][j][i-1].z + .5 * coors[k][j][i].z;\
            pt1[2] = .5 * coors[k][j-1][i-1].z + .5 * coors[k][j-1][i].z;\
            pt2[2] = .5 * coors[k][j-2][i-1].z + .5 * coors[k][j-2][i].z;\
            x_s[xh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
            y_s[xh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
            z_s[xh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
            pt0[0] = .5 * coors[k-1][j][i-1].x + .5 * coors[k-1][j][i].x;\
            pt1[0] = .5 * coors[k+1][j][i-1].x + .5 * coors[k+1][j][i].x;\
            pt0[1] = .5 * coors[k-1][j][i-1].y + .5 * coors[k-1][j][i].y;\
            pt1[1] = .5 * coors[k+1][j][i-1].y + .5 * coors[k+1][j][i].y;\
            pt0[2] = .5 * coors[k-1][j][i-1].z + .5 * coors[k-1][j][i].z;\
            pt1[2] = .5 * coors[k+1][j][i-1].z + .5 * coors[k+1][j][i].z;\
            x_t[xh] = (pt1[0] - pt0[0]) * .5;\
            y_t[xh] = (pt1[1] - pt0[1]) * .5;\
            z_t[xh] = (pt1[2] - pt0[2]) * .5;\
            g_11 = x_r[xh]*x_r[xh] + y_r[xh]*y_r[xh] + z_r[xh]*z_r[xh];\
            g_12 = x_r[xh]*x_s[xh] + y_r[xh]*y_s[xh] + z_r[xh]*z_s[xh];\
            g_13 = x_r[xh]*x_t[xh] + y_r[xh]*y_t[xh] + z_r[xh]*z_t[xh];\
            g_22 = x_s[xh]*x_s[xh] + y_s[xh]*y_s[xh] + z_s[xh]*z_s[xh];\
            g_23 = x_s[xh]*x_t[xh] + y_s[xh]*y_t[xh] + z_s[xh]*z_t[xh];\
            g_33 = x_t[xh]*x_t[xh] + y_t[xh]*y_t[xh] + z_t[xh]*z_t[xh];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            g11 = FDIV(1.0,g) * (g_22*g_33 - g_23*g_23);\
            coef[MET_W][k][j][i] = sqrt(g) * g11;\
\
        }\
\
    }\
    for (k = kbeg; k < zs + zm; k++)\
    {\
        for (i = ibeg; i < iend; i++)\
        {\
            pt0[0] = .5 * coors[k-1][j][i-1].x + .5 * coors[k][j][i-1].x;\
            pt1[0] = .5 * coors[k-1][j][i+1].x + .5 * coors[k][j][i+1].x;\
            pt0[1] = .5 * coors[k-1][j][i-1].y + .5 * coors[k][j][i-1].y;\
            pt1[1] = .5 * coors[k-1][j][i+1].y + .5 * coors[k][j][i+1].y;\
            pt0[2] = .5 * coors[k-1][j][i-1].z + .5 * coors[k][j][i-1].z;\
            pt1[2] = .5 * coors[k-1][j][i+1].z + .5 * coors[k][j][i+1].z;\
            x_r[zh] = (pt1[0] - pt0[0]) * .5;\
            y_r[zh] = (pt1[1] - pt0[1]) * .5;\
            z_r[zh] = (pt1[2] - pt0[2]) * .5;\
            pt0[0] = .5 * coors[k-1][j][i].x + .5 * coors[k][j][i].x;\
            pt1[0] = .5 * coors[k-1][j-1][i].x + .5 * coors[k][j-1][i].x;\
            pt2[0] = .5 * coors[k-1][j-2][i].x + .5 * coors[k][j-2][i].x;\
            pt0[1] = .5 * coors[k-1][j][i].y + .5 * coors[k][j][i].y;\
            pt1[1] = .5 * coors[k-1][j-1][i].y + .5 * coors[k][j-1][i].y;\
            pt2[1] = .5 * coors[k-1][j-2][i].y + .5 * coors[k][j-2][i].y;\
            pt0[2] = .5 * coors[k-1][j][i].z + .5 * coors[k][j][i].z;\
            pt1[2] = .5 * coors[k-1][j-1][i].z + .5 * coors[k][j-1][i].z;\
            pt2[2] = .5 * coors[k-1][j-2][i].z + .5 * coors[k][j-2][i].z;\
            x_s[zh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
            y_s[zh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
            z_s[zh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
            x_t[zh] = coors[k][j][i].x - coors[k-1][j][i].x;\
            y_t[zh] = coors[k][j][i].y - coors[k-1][j][i].y;\
            z_t[zh] = coors[k][j][i].z - coors[k-1][j][i].z;\
            g_11 = x_r[zh]*x_r[zh] + y_r[zh]*y_r[zh] + z_r[zh]*z_r[zh];\
            g_12 = x_r[zh]*x_s[zh] + y_r[zh]*y_s[zh] + z_r[zh]*z_s[zh];\
            g_13 = x_r[zh]*x_t[zh] + y_r[zh]*y_t[zh] + z_r[zh]*z_t[zh];\
            g_22 = x_s[zh]*x_s[zh] + y_s[zh]*y_s[zh] + z_s[zh]*z_s[zh];\
            g_23 = x_s[zh]*x_t[zh] + y_s[zh]*y_t[zh] + z_s[zh]*z_t[zh];\
            g_33 = x_t[zh]*x_t[zh] + y_t[zh]*y_t[zh] + z_t[zh]*z_t[zh];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            g33 = FDIV(1.0,g) * (g_11*g_22 - g_12*g_12);\
            coef[MET_B][k][j][i] = sqrt(g) * g33;\
\
        }\
\
    }\
    for (k = kbeg; k < kend; k++)\
    {\
        for (i = ibeg; i < iend; i++)\
        {\
            x_r[c] = si[-1] * coors[k][j][i-1].x + si[1] * coors[k][j][i+1].x;\
            y_r[c] = si[-1] * coors[k][j][i-1].y + si[1] * coors[k][j][i+1].y;\
            z_r[c] = si[-1] * coors[k][j][i-1].z + si[1] * coors[k][j][i+1].z;\
            x_s[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j-1][i].x + sp[-2] * coors[k][j-2][i].x;\
            y_s[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j-1][i].y + sp[-2] * coors[k][j-2][i].y;\
            z_s[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j-1][i].z + sp[-2] * coors[k][j-2][i].z;\
            x_t[c] = si[-1] * coors[k-1][j][i].x + si[1] * coors[k+1][j][i].x;\
            y_t[c] = si[-1] * coors[k-1][j][i].y + si[1] * coors[k+1][j][i].y;\
            z_t[c] = si[-1] * coors[k-1][j][i].z + si[1] * coors[k+1][j][i].z;\
            g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
            g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
            g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
            g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
            g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
            g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            jac[k][j][i] = sqrt(g);\
            vec__x_r[k][j][i] = x_r[c];\
            vec__y_r[k][j][i] = y_r[c];\
            vec__z_r[k][j][i] = z_r[c];\
            vec__x_s[k][j][i] = x_s[c];\
            vec__y_s[k][j][i] = y_s[c];\
            vec__z_s[k][j][i] = z_s[c];\
            vec__x_t[k][j][i] = x_t[c];\
            vec__y_t[k][j][i] = y_t[c];\
            vec__z_t[k][j][i] = z_t[c];\
            g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
            g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
            g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
            coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
            coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
            coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
            coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
            coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
            coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
            coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
            coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
            coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
\
        }\
\
    }\
\
}\
\
\
if (ys == 0 && !per_y) {\
    j = 0;\
    for (k = kbeg; k < kend; k++)\
    {\
        for (i = ibeg; i < xs + xm; i++)\
        {\
            x_r[xh] = coors[k][j][i].x - coors[k][j][i-1].x;\
            y_r[xh] = coors[k][j][i].y - coors[k][j][i-1].y;\
            z_r[xh] = coors[k][j][i].z - coors[k][j][i-1].z;\
            pt0[0] = .5 * coors[k][j][i-1].x + .5 * coors[k][j][i].x;\
            pt1[0] = .5 * coors[k][j+1][i-1].x + .5 * coors[k][j+1][i].x;\
            pt2[0] = .5 * coors[k][j+2][i-1].x + .5 * coors[k][j+2][i].x;\
            pt0[1] = .5 * coors[k][j][i-1].y + .5 * coors[k][j][i].y;\
            pt1[1] = .5 * coors[k][j+1][i-1].y + .5 * coors[k][j+1][i].y;\
            pt2[1] = .5 * coors[k][j+2][i-1].y + .5 * coors[k][j+2][i].y;\
            pt0[2] = .5 * coors[k][j][i-1].z + .5 * coors[k][j][i].z;\
            pt1[2] = .5 * coors[k][j+1][i-1].z + .5 * coors[k][j+1][i].z;\
            pt2[2] = .5 * coors[k][j+2][i-1].z + .5 * coors[k][j+2][i].z;\
            x_s[xh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
            y_s[xh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
            z_s[xh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
            pt0[0] = .5 * coors[k-1][j][i-1].x + .5 * coors[k-1][j][i].x;\
            pt1[0] = .5 * coors[k+1][j][i-1].x + .5 * coors[k+1][j][i].x;\
            pt0[1] = .5 * coors[k-1][j][i-1].y + .5 * coors[k-1][j][i].y;\
            pt1[1] = .5 * coors[k+1][j][i-1].y + .5 * coors[k+1][j][i].y;\
            pt0[2] = .5 * coors[k-1][j][i-1].z + .5 * coors[k-1][j][i].z;\
            pt1[2] = .5 * coors[k+1][j][i-1].z + .5 * coors[k+1][j][i].z;\
            x_t[xh] = (pt1[0] - pt0[0]) * .5;\
            y_t[xh] = (pt1[1] - pt0[1]) * .5;\
            z_t[xh] = (pt1[2] - pt0[2]) * .5;\
            g_11 = x_r[xh]*x_r[xh] + y_r[xh]*y_r[xh] + z_r[xh]*z_r[xh];\
            g_12 = x_r[xh]*x_s[xh] + y_r[xh]*y_s[xh] + z_r[xh]*z_s[xh];\
            g_13 = x_r[xh]*x_t[xh] + y_r[xh]*y_t[xh] + z_r[xh]*z_t[xh];\
            g_22 = x_s[xh]*x_s[xh] + y_s[xh]*y_s[xh] + z_s[xh]*z_s[xh];\
            g_23 = x_s[xh]*x_t[xh] + y_s[xh]*y_t[xh] + z_s[xh]*z_t[xh];\
            g_33 = x_t[xh]*x_t[xh] + y_t[xh]*y_t[xh] + z_t[xh]*z_t[xh];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            g11 = FDIV(1.0,g) * (g_22*g_33 - g_23*g_23);\
            coef[MET_W][k][j][i] = sqrt(g) * g11;\
\
        }\
\
    }\
    for (k = kbeg; k < zs + zm; k++)\
    {\
        for (i = ibeg; i < iend; i++)\
        {\
            pt0[0] = .5 * coors[k-1][j][i-1].x + .5 * coors[k][j][i-1].x;\
            pt1[0] = .5 * coors[k-1][j][i+1].x + .5 * coors[k][j][i+1].x;\
            pt0[1] = .5 * coors[k-1][j][i-1].y + .5 * coors[k][j][i-1].y;\
            pt1[1] = .5 * coors[k-1][j][i+1].y + .5 * coors[k][j][i+1].y;\
            pt0[2] = .5 * coors[k-1][j][i-1].z + .5 * coors[k][j][i-1].z;\
            pt1[2] = .5 * coors[k-1][j][i+1].z + .5 * coors[k][j][i+1].z;\
            x_r[zh] = (pt1[0] - pt0[0]) * .5;\
            y_r[zh] = (pt1[1] - pt0[1]) * .5;\
            z_r[zh] = (pt1[2] - pt0[2]) * .5;\
            pt0[0] = .5 * coors[k-1][j][i].x + .5 * coors[k][j][i].x;\
            pt1[0] = .5 * coors[k-1][j+1][i].x + .5 * coors[k][j+1][i].x;\
            pt2[0] = .5 * coors[k-1][j+2][i].x + .5 * coors[k][j+2][i].x;\
            pt0[1] = .5 * coors[k-1][j][i].y + .5 * coors[k][j][i].y;\
            pt1[1] = .5 * coors[k-1][j+1][i].y + .5 * coors[k][j+1][i].y;\
            pt2[1] = .5 * coors[k-1][j+2][i].y + .5 * coors[k][j+2][i].y;\
            pt0[2] = .5 * coors[k-1][j][i].z + .5 * coors[k][j][i].z;\
            pt1[2] = .5 * coors[k-1][j+1][i].z + .5 * coors[k][j+1][i].z;\
            pt2[2] = .5 * coors[k-1][j+2][i].z + .5 * coors[k][j+2][i].z;\
            x_s[zh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
            y_s[zh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
            z_s[zh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
            x_t[zh] = coors[k][j][i].x - coors[k-1][j][i].x;\
            y_t[zh] = coors[k][j][i].y - coors[k-1][j][i].y;\
            z_t[zh] = coors[k][j][i].z - coors[k-1][j][i].z;\
            g_11 = x_r[zh]*x_r[zh] + y_r[zh]*y_r[zh] + z_r[zh]*z_r[zh];\
            g_12 = x_r[zh]*x_s[zh] + y_r[zh]*y_s[zh] + z_r[zh]*z_s[zh];\
            g_13 = x_r[zh]*x_t[zh] + y_r[zh]*y_t[zh] + z_r[zh]*z_t[zh];\
            g_22 = x_s[zh]*x_s[zh] + y_s[zh]*y_s[zh] + z_s[zh]*z_s[zh];\
            g_23 = x_s[zh]*x_t[zh] + y_s[zh]*y_t[zh] + z_s[zh]*z_t[zh];\
            g_33 = x_t[zh]*x_t[zh] + y_t[zh]*y_t[zh] + z_t[zh]*z_t[zh];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            g33 = FDIV(1.0,g) * (g_11*g_22 - g_12*g_12);\
            coef[MET_B][k][j][i] = sqrt(g) * g33;\
\
        }\
\
    }\
    for (k = kbeg; k < kend; k++)\
    {\
        for (i = ibeg; i < iend; i++)\
        {\
            x_r[c] = si[-1] * coors[k][j][i-1].x + si[1] * coors[k][j][i+1].x;\
            y_r[c] = si[-1] * coors[k][j][i-1].y + si[1] * coors[k][j][i+1].y;\
            z_r[c] = si[-1] * coors[k][j][i-1].z + si[1] * coors[k][j][i+1].z;\
            x_s[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j+1][i].x + sm[2] * coors[k][j+2][i].x;\
            y_s[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j+1][i].y + sm[2] * coors[k][j+2][i].y;\
            z_s[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j+1][i].z + sm[2] * coors[k][j+2][i].z;\
            x_t[c] = si[-1] * coors[k-1][j][i].x + si[1] * coors[k+1][j][i].x;\
            y_t[c] = si[-1] * coors[k-1][j][i].y + si[1] * coors[k+1][j][i].y;\
            z_t[c] = si[-1] * coors[k-1][j][i].z + si[1] * coors[k+1][j][i].z;\
            g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
            g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
            g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
            g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
            g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
            g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            jac[k][j][i] = sqrt(g);\
            vec__x_r[k][j][i] = x_r[c];\
            vec__y_r[k][j][i] = y_r[c];\
            vec__z_r[k][j][i] = z_r[c];\
            vec__x_s[k][j][i] = x_s[c];\
            vec__y_s[k][j][i] = y_s[c];\
            vec__z_s[k][j][i] = z_s[c];\
            vec__x_t[k][j][i] = x_t[c];\
            vec__y_t[k][j][i] = y_t[c];\
            vec__z_t[k][j][i] = z_t[c];\
            g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
            g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
            g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
            coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
            coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
            coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
            coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
            coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
            coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
            coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
            coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
            coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
            coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
            coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
\
        }\
\
    }\
\
}\
\
\
if ((zs + zm == ngz && !per_z) && !per_z) {\
    k = ngz - 1;\
    for (j = jbeg; j < jend; j++)\
    {\
        for (i = ibeg; i < xs + xm; i++)\
        {\
            x_r[xh] = coors[k][j][i].x - coors[k][j][i-1].x;\
            y_r[xh] = coors[k][j][i].y - coors[k][j][i-1].y;\
            z_r[xh] = coors[k][j][i].z - coors[k][j][i-1].z;\
            pt0[0] = .5 * coors[k][j-1][i-1].x + .5 * coors[k][j-1][i].x;\
            pt1[0] = .5 * coors[k][j+1][i-1].x + .5 * coors[k][j+1][i].x;\
            pt0[1] = .5 * coors[k][j-1][i-1].y + .5 * coors[k][j-1][i].y;\
            pt1[1] = .5 * coors[k][j+1][i-1].y + .5 * coors[k][j+1][i].y;\
            pt0[2] = .5 * coors[k][j-1][i-1].z + .5 * coors[k][j-1][i].z;\
            pt1[2] = .5 * coors[k][j+1][i-1].z + .5 * coors[k][j+1][i].z;\
            x_s[xh] = (pt1[0] - pt0[0]) * .5;\
            y_s[xh] = (pt1[1] - pt0[1]) * .5;\
            z_s[xh] = (pt1[2] - pt0[2]) * .5;\
            pt0[0] = .5 * coors[k][j][i-1].x + .5 * coors[k][j][i].x;\
            pt1[0] = .5 * coors[k-1][j][i-1].x + .5 * coors[k-1][j][i].x;\
            pt2[0] = .5 * coors[k-2][j][i-1].x + .5 * coors[k-2][j][i].x;\
            pt0[1] = .5 * coors[k][j][i-1].y + .5 * coors[k][j][i].y;\
            pt1[1] = .5 * coors[k-1][j][i-1].y + .5 * coors[k-1][j][i].y;\
            pt2[1] = .5 * coors[k-2][j][i-1].y + .5 * coors[k-2][j][i].y;\
            pt0[2] = .5 * coors[k][j][i-1].z + .5 * coors[k][j][i].z;\
            pt1[2] = .5 * coors[k-1][j][i-1].z + .5 * coors[k-1][j][i].z;\
            pt2[2] = .5 * coors[k-2][j][i-1].z + .5 * coors[k-2][j][i].z;\
            x_t[xh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
            y_t[xh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
            z_t[xh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
            g_11 = x_r[xh]*x_r[xh] + y_r[xh]*y_r[xh] + z_r[xh]*z_r[xh];\
            g_12 = x_r[xh]*x_s[xh] + y_r[xh]*y_s[xh] + z_r[xh]*z_s[xh];\
            g_13 = x_r[xh]*x_t[xh] + y_r[xh]*y_t[xh] + z_r[xh]*z_t[xh];\
            g_22 = x_s[xh]*x_s[xh] + y_s[xh]*y_s[xh] + z_s[xh]*z_s[xh];\
            g_23 = x_s[xh]*x_t[xh] + y_s[xh]*y_t[xh] + z_s[xh]*z_t[xh];\
            g_33 = x_t[xh]*x_t[xh] + y_t[xh]*y_t[xh] + z_t[xh]*z_t[xh];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            g11 = FDIV(1.0,g) * (g_22*g_33 - g_23*g_23);\
            coef[MET_W][k][j][i] = sqrt(g) * g11;\
\
        }\
\
    }\
    for (j = jbeg; j < ys + ym; j++)\
    {\
        for (i = ibeg; i < iend; i++)\
        {\
            pt0[0] = .5 * coors[k][j-1][i-1].x + .5 * coors[k][j][i-1].x;\
            pt1[0] = .5 * coors[k][j-1][i+1].x + .5 * coors[k][j][i+1].x;\
            pt0[1] = .5 * coors[k][j-1][i-1].y + .5 * coors[k][j][i-1].y;\
            pt1[1] = .5 * coors[k][j-1][i+1].y + .5 * coors[k][j][i+1].y;\
            pt0[2] = .5 * coors[k][j-1][i-1].z + .5 * coors[k][j][i-1].z;\
            pt1[2] = .5 * coors[k][j-1][i+1].z + .5 * coors[k][j][i+1].z;\
            x_r[yh] = (pt1[0] - pt0[0]) * .5;\
            y_r[yh] = (pt1[1] - pt0[1]) * .5;\
            z_r[yh] = (pt1[2] - pt0[2]) * .5;\
            x_s[yh] = coors[k][j][i].x - coors[k][j-1][i].x;\
            y_s[yh] = coors[k][j][i].y - coors[k][j-1][i].y;\
            z_s[yh] = coors[k][j][i].z - coors[k][j-1][i].z;\
            pt0[0] = .5 * coors[k][j-1][i].x + .5 * coors[k][j][i].x;\
            pt1[0] = .5 * coors[k-1][j-1][i].x + .5 * coors[k-1][j][i].x;\
            pt2[0] = .5 * coors[k-2][j-1][i].x + .5 * coors[k-2][j][i].x;\
            pt0[1] = .5 * coors[k][j-1][i].y + .5 * coors[k][j][i].y;\
            pt1[1] = .5 * coors[k-1][j-1][i].y + .5 * coors[k-1][j][i].y;\
            pt2[1] = .5 * coors[k-2][j-1][i].y + .5 * coors[k-2][j][i].y;\
            pt0[2] = .5 * coors[k][j-1][i].z + .5 * coors[k][j][i].z;\
            pt1[2] = .5 * coors[k-1][j-1][i].z + .5 * coors[k-1][j][i].z;\
            pt2[2] = .5 * coors[k-2][j-1][i].z + .5 * coors[k-2][j][i].z;\
            x_t[yh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
            y_t[yh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
            z_t[yh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
            g_11 = x_r[yh]*x_r[yh] + y_r[yh]*y_r[yh] + z_r[yh]*z_r[yh];\
            g_12 = x_r[yh]*x_s[yh] + y_r[yh]*y_s[yh] + z_r[yh]*z_s[yh];\
            g_13 = x_r[yh]*x_t[yh] + y_r[yh]*y_t[yh] + z_r[yh]*z_t[yh];\
            g_22 = x_s[yh]*x_s[yh] + y_s[yh]*y_s[yh] + z_s[yh]*z_s[yh];\
            g_23 = x_s[yh]*x_t[yh] + y_s[yh]*y_t[yh] + z_s[yh]*z_t[yh];\
            g_33 = x_t[yh]*x_t[yh] + y_t[yh]*y_t[yh] + z_t[yh]*z_t[yh];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            g22 = FDIV(1.0,g) * (g_11*g_33 - g_13*g_13);\
            coef[MET_S][k][j][i] = sqrt(g) * g22;\
\
        }\
\
    }\
    for (j = jbeg; j < jend; j++)\
    {\
        for (i = ibeg; i < iend; i++)\
        {\
            x_r[c] = si[-1] * coors[k][j][i-1].x + si[1] * coors[k][j][i+1].x;\
            y_r[c] = si[-1] * coors[k][j][i-1].y + si[1] * coors[k][j][i+1].y;\
            z_r[c] = si[-1] * coors[k][j][i-1].z + si[1] * coors[k][j][i+1].z;\
            x_s[c] = si[-1] * coors[k][j-1][i].x + si[1] * coors[k][j+1][i].x;\
            y_s[c] = si[-1] * coors[k][j-1][i].y + si[1] * coors[k][j+1][i].y;\
            z_s[c] = si[-1] * coors[k][j-1][i].z + si[1] * coors[k][j+1][i].z;\
            x_t[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k-1][j][i].x + sp[-2] * coors[k-2][j][i].x;\
            y_t[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k-1][j][i].y + sp[-2] * coors[k-2][j][i].y;\
            z_t[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k-1][j][i].z + sp[-2] * coors[k-2][j][i].z;\
            g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
            g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
            g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
            g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
            g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
            g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            jac[k][j][i] = sqrt(g);\
            vec__x_r[k][j][i] = x_r[c];\
            vec__y_r[k][j][i] = y_r[c];\
            vec__z_r[k][j][i] = z_r[c];\
            vec__x_s[k][j][i] = x_s[c];\
            vec__y_s[k][j][i] = y_s[c];\
            vec__z_s[k][j][i] = z_s[c];\
            vec__x_t[k][j][i] = x_t[c];\
            vec__y_t[k][j][i] = y_t[c];\
            vec__z_t[k][j][i] = z_t[c];\
            g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
            g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
            g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
            coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
            coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
            coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
            coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
            coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
            coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
            coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
            coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
\
        }\
\
    }\
\
}\
\
\
if (zs == 0 && !per_z) {\
    k = 0;\
    for (j = jbeg; j < jend; j++)\
    {\
        for (i = ibeg; i < xs + xm; i++)\
        {\
            x_r[xh] = coors[k][j][i].x - coors[k][j][i-1].x;\
            y_r[xh] = coors[k][j][i].y - coors[k][j][i-1].y;\
            z_r[xh] = coors[k][j][i].z - coors[k][j][i-1].z;\
            pt0[0] = .5 * coors[k][j-1][i-1].x + .5 * coors[k][j-1][i].x;\
            pt1[0] = .5 * coors[k][j+1][i-1].x + .5 * coors[k][j+1][i].x;\
            pt0[1] = .5 * coors[k][j-1][i-1].y + .5 * coors[k][j-1][i].y;\
            pt1[1] = .5 * coors[k][j+1][i-1].y + .5 * coors[k][j+1][i].y;\
            pt0[2] = .5 * coors[k][j-1][i-1].z + .5 * coors[k][j-1][i].z;\
            pt1[2] = .5 * coors[k][j+1][i-1].z + .5 * coors[k][j+1][i].z;\
            x_s[xh] = (pt1[0] - pt0[0]) * .5;\
            y_s[xh] = (pt1[1] - pt0[1]) * .5;\
            z_s[xh] = (pt1[2] - pt0[2]) * .5;\
            pt0[0] = .5 * coors[k][j][i-1].x + .5 * coors[k][j][i].x;\
            pt1[0] = .5 * coors[k+1][j][i-1].x + .5 * coors[k+1][j][i].x;\
            pt2[0] = .5 * coors[k+2][j][i-1].x + .5 * coors[k+2][j][i].x;\
            pt0[1] = .5 * coors[k][j][i-1].y + .5 * coors[k][j][i].y;\
            pt1[1] = .5 * coors[k+1][j][i-1].y + .5 * coors[k+1][j][i].y;\
            pt2[1] = .5 * coors[k+2][j][i-1].y + .5 * coors[k+2][j][i].y;\
            pt0[2] = .5 * coors[k][j][i-1].z + .5 * coors[k][j][i].z;\
            pt1[2] = .5 * coors[k+1][j][i-1].z + .5 * coors[k+1][j][i].z;\
            pt2[2] = .5 * coors[k+2][j][i-1].z + .5 * coors[k+2][j][i].z;\
            x_t[xh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
            y_t[xh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
            z_t[xh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
            g_11 = x_r[xh]*x_r[xh] + y_r[xh]*y_r[xh] + z_r[xh]*z_r[xh];\
            g_12 = x_r[xh]*x_s[xh] + y_r[xh]*y_s[xh] + z_r[xh]*z_s[xh];\
            g_13 = x_r[xh]*x_t[xh] + y_r[xh]*y_t[xh] + z_r[xh]*z_t[xh];\
            g_22 = x_s[xh]*x_s[xh] + y_s[xh]*y_s[xh] + z_s[xh]*z_s[xh];\
            g_23 = x_s[xh]*x_t[xh] + y_s[xh]*y_t[xh] + z_s[xh]*z_t[xh];\
            g_33 = x_t[xh]*x_t[xh] + y_t[xh]*y_t[xh] + z_t[xh]*z_t[xh];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            g11 = FDIV(1.0,g) * (g_22*g_33 - g_23*g_23);\
            coef[MET_W][k][j][i] = sqrt(g) * g11;\
\
        }\
\
    }\
    for (j = jbeg; j < ys + ym; j++)\
    {\
        for (i = ibeg; i < iend; i++)\
        {\
            pt0[0] = .5 * coors[k][j-1][i-1].x + .5 * coors[k][j][i-1].x;\
            pt1[0] = .5 * coors[k][j-1][i+1].x + .5 * coors[k][j][i+1].x;\
            pt0[1] = .5 * coors[k][j-1][i-1].y + .5 * coors[k][j][i-1].y;\
            pt1[1] = .5 * coors[k][j-1][i+1].y + .5 * coors[k][j][i+1].y;\
            pt0[2] = .5 * coors[k][j-1][i-1].z + .5 * coors[k][j][i-1].z;\
            pt1[2] = .5 * coors[k][j-1][i+1].z + .5 * coors[k][j][i+1].z;\
            x_r[yh] = (pt1[0] - pt0[0]) * .5;\
            y_r[yh] = (pt1[1] - pt0[1]) * .5;\
            z_r[yh] = (pt1[2] - pt0[2]) * .5;\
            x_s[yh] = coors[k][j][i].x - coors[k][j-1][i].x;\
            y_s[yh] = coors[k][j][i].y - coors[k][j-1][i].y;\
            z_s[yh] = coors[k][j][i].z - coors[k][j-1][i].z;\
            pt0[0] = .5 * coors[k][j-1][i].x + .5 * coors[k][j][i].x;\
            pt1[0] = .5 * coors[k+1][j-1][i].x + .5 * coors[k+1][j][i].x;\
            pt2[0] = .5 * coors[k+2][j-1][i].x + .5 * coors[k+2][j][i].x;\
            pt0[1] = .5 * coors[k][j-1][i].y + .5 * coors[k][j][i].y;\
            pt1[1] = .5 * coors[k+1][j-1][i].y + .5 * coors[k+1][j][i].y;\
            pt2[1] = .5 * coors[k+2][j-1][i].y + .5 * coors[k+2][j][i].y;\
            pt0[2] = .5 * coors[k][j-1][i].z + .5 * coors[k][j][i].z;\
            pt1[2] = .5 * coors[k+1][j-1][i].z + .5 * coors[k+1][j][i].z;\
            pt2[2] = .5 * coors[k+2][j-1][i].z + .5 * coors[k+2][j][i].z;\
            x_t[yh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
            y_t[yh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
            z_t[yh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
            g_11 = x_r[yh]*x_r[yh] + y_r[yh]*y_r[yh] + z_r[yh]*z_r[yh];\
            g_12 = x_r[yh]*x_s[yh] + y_r[yh]*y_s[yh] + z_r[yh]*z_s[yh];\
            g_13 = x_r[yh]*x_t[yh] + y_r[yh]*y_t[yh] + z_r[yh]*z_t[yh];\
            g_22 = x_s[yh]*x_s[yh] + y_s[yh]*y_s[yh] + z_s[yh]*z_s[yh];\
            g_23 = x_s[yh]*x_t[yh] + y_s[yh]*y_t[yh] + z_s[yh]*z_t[yh];\
            g_33 = x_t[yh]*x_t[yh] + y_t[yh]*y_t[yh] + z_t[yh]*z_t[yh];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            g22 = FDIV(1.0,g) * (g_11*g_33 - g_13*g_13);\
            coef[MET_S][k][j][i] = sqrt(g) * g22;\
\
        }\
\
    }\
    for (j = jbeg; j < jend; j++)\
    {\
        for (i = ibeg; i < iend; i++)\
        {\
            x_r[c] = si[-1] * coors[k][j][i-1].x + si[1] * coors[k][j][i+1].x;\
            y_r[c] = si[-1] * coors[k][j][i-1].y + si[1] * coors[k][j][i+1].y;\
            z_r[c] = si[-1] * coors[k][j][i-1].z + si[1] * coors[k][j][i+1].z;\
            x_s[c] = si[-1] * coors[k][j-1][i].x + si[1] * coors[k][j+1][i].x;\
            y_s[c] = si[-1] * coors[k][j-1][i].y + si[1] * coors[k][j+1][i].y;\
            z_s[c] = si[-1] * coors[k][j-1][i].z + si[1] * coors[k][j+1][i].z;\
            x_t[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k+1][j][i].x + sm[2] * coors[k+2][j][i].x;\
            y_t[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k+1][j][i].y + sm[2] * coors[k+2][j][i].y;\
            z_t[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k+1][j][i].z + sm[2] * coors[k+2][j][i].z;\
            g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
            g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
            g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
            g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
            g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
            g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
\
            g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
            jac[k][j][i] = sqrt(g);\
            vec__x_r[k][j][i] = x_r[c];\
            vec__y_r[k][j][i] = y_r[c];\
            vec__z_r[k][j][i] = z_r[c];\
            vec__x_s[k][j][i] = x_s[c];\
            vec__y_s[k][j][i] = y_s[c];\
            vec__z_s[k][j][i] = z_s[c];\
            vec__x_t[k][j][i] = x_t[c];\
            vec__y_t[k][j][i] = y_t[c];\
            vec__z_t[k][j][i] = z_t[c];\
            g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
            g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
            g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
            coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
            coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
            coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
            coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
            coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
            coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
            coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
            coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
            coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
            coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
            coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
            coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
\
        }\
\
    }\
\
}\
\
\
if ((ys + ym == ngy && !per_y) && ((zs + zm == ngz && !per_z) && !per_z)) {	\
    j = ngy - 1;\
    k = ngz - 1;\
    for (i = ibeg; i < xs + xm; i++)\
    {\
        x_r[xh] = coors[k][j][i].x - coors[k][j][i-1].x;\
        y_r[xh] = coors[k][j][i].y - coors[k][j][i-1].y;\
        z_r[xh] = coors[k][j][i].z - coors[k][j][i-1].z;\
        pt0[0] = .5 * coors[k][j][i-1].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k][j-1][i-1].x + .5 * coors[k][j-1][i].x;\
        pt2[0] = .5 * coors[k][j-2][i-1].x + .5 * coors[k][j-2][i].x;\
        pt0[1] = .5 * coors[k][j][i-1].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k][j-1][i-1].y + .5 * coors[k][j-1][i].y;\
        pt2[1] = .5 * coors[k][j-2][i-1].y + .5 * coors[k][j-2][i].y;\
        pt0[2] = .5 * coors[k][j][i-1].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k][j-1][i-1].z + .5 * coors[k][j-1][i].z;\
        pt2[2] = .5 * coors[k][j-2][i-1].z + .5 * coors[k][j-2][i].z;\
        x_s[xh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
        y_s[xh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
        z_s[xh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
        pt0[0] = .5 * coors[k][j][i-1].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k-1][j][i-1].x + .5 * coors[k-1][j][i].x;\
        pt2[0] = .5 * coors[k-2][j][i-1].x + .5 * coors[k-2][j][i].x;\
        pt0[1] = .5 * coors[k][j][i-1].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k-1][j][i-1].y + .5 * coors[k-1][j][i].y;\
        pt2[1] = .5 * coors[k-2][j][i-1].y + .5 * coors[k-2][j][i].y;\
        pt0[2] = .5 * coors[k][j][i-1].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k-1][j][i-1].z + .5 * coors[k-1][j][i].z;\
        pt2[2] = .5 * coors[k-2][j][i-1].z + .5 * coors[k-2][j][i].z;\
        x_t[xh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
        y_t[xh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
        z_t[xh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
        g_11 = x_r[xh]*x_r[xh] + y_r[xh]*y_r[xh] + z_r[xh]*z_r[xh];\
        g_12 = x_r[xh]*x_s[xh] + y_r[xh]*y_s[xh] + z_r[xh]*z_s[xh];\
        g_13 = x_r[xh]*x_t[xh] + y_r[xh]*y_t[xh] + z_r[xh]*z_t[xh];\
        g_22 = x_s[xh]*x_s[xh] + y_s[xh]*y_s[xh] + z_s[xh]*z_s[xh];\
        g_23 = x_s[xh]*x_t[xh] + y_s[xh]*y_t[xh] + z_s[xh]*z_t[xh];\
        g_33 = x_t[xh]*x_t[xh] + y_t[xh]*y_t[xh] + z_t[xh]*z_t[xh];\
\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        g11 = FDIV(1.0,g) * (g_22*g_33 - g_23*g_23);\
        coef[MET_W][k][j][i] = sqrt(g) * g11;\
\
    }\
\
\
    for (i = ibeg; i < iend; i++)\
    {\
        x_r[c] = si[-1] * coors[k][j][i-1].x + si[1] * coors[k][j][i+1].x;\
        y_r[c] = si[-1] * coors[k][j][i-1].y + si[1] * coors[k][j][i+1].y;\
        z_r[c] = si[-1] * coors[k][j][i-1].z + si[1] * coors[k][j][i+1].z;\
        x_s[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j-1][i].x + sp[-2] * coors[k][j-2][i].x;\
        y_s[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j-1][i].y + sp[-2] * coors[k][j-2][i].y;\
        z_s[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j-1][i].z + sp[-2] * coors[k][j-2][i].z;\
        x_t[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k-1][j][i].x + sp[-2] * coors[k-2][j][i].x;\
        y_t[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k-1][j][i].y + sp[-2] * coors[k-2][j][i].y;\
        z_t[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k-1][j][i].z + sp[-2] * coors[k-2][j][i].z;\
        g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
        g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
        g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
        g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
        g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
        g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        jac[k][j][i] = sqrt(g);\
        vec__x_r[k][j][i] = x_r[c];\
        vec__y_r[k][j][i] = y_r[c];\
        vec__z_r[k][j][i] = z_r[c];\
        vec__x_s[k][j][i] = x_s[c];\
        vec__y_s[k][j][i] = y_s[c];\
        vec__z_s[k][j][i] = z_s[c];\
        vec__x_t[k][j][i] = x_t[c];\
        vec__y_t[k][j][i] = y_t[c];\
        vec__z_t[k][j][i] = z_t[c];\
        g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
        g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
        g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
        coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
        coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
        coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
        coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
        coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
\
    }\
\
}\
\
\
if ((xs + xm == ngx && !per_x) && ((zs + zm == ngz && !per_z) && !per_z)) {	\
    i = ngx - 1;\
    k = ngz - 1;\
    for (j = jbeg; j < ys + ym; j++)\
    {\
        pt0[0] = .5 * coors[k][j-1][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k][j-1][i-1].x + .5 * coors[k][j][i-1].x;\
        pt2[0] = .5 * coors[k][j-1][i-2].x + .5 * coors[k][j][i-2].x;\
        pt0[1] = .5 * coors[k][j-1][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k][j-1][i-1].y + .5 * coors[k][j][i-1].y;\
        pt2[1] = .5 * coors[k][j-1][i-2].y + .5 * coors[k][j][i-2].y;\
        pt0[2] = .5 * coors[k][j-1][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k][j-1][i-1].z + .5 * coors[k][j][i-1].z;\
        pt2[2] = .5 * coors[k][j-1][i-2].z + .5 * coors[k][j][i-2].z;\
        x_r[yh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
        y_r[yh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
        z_r[yh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
        x_s[yh] = coors[k][j][i].x - coors[k][j-1][i].x;\
        y_s[yh] = coors[k][j][i].y - coors[k][j-1][i].y;\
        z_s[yh] = coors[k][j][i].z - coors[k][j-1][i].z;\
        pt0[0] = .5 * coors[k][j-1][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k-1][j-1][i].x + .5 * coors[k-1][j][i].x;\
        pt2[0] = .5 * coors[k-2][j-1][i].x + .5 * coors[k-2][j][i].x;\
        pt0[1] = .5 * coors[k][j-1][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k-1][j-1][i].y + .5 * coors[k-1][j][i].y;\
        pt2[1] = .5 * coors[k-2][j-1][i].y + .5 * coors[k-2][j][i].y;\
        pt0[2] = .5 * coors[k][j-1][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k-1][j-1][i].z + .5 * coors[k-1][j][i].z;\
        pt2[2] = .5 * coors[k-2][j-1][i].z + .5 * coors[k-2][j][i].z;\
        x_t[yh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
        y_t[yh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
        z_t[yh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
        g_11 = x_r[yh]*x_r[yh] + y_r[yh]*y_r[yh] + z_r[yh]*z_r[yh];\
        g_12 = x_r[yh]*x_s[yh] + y_r[yh]*y_s[yh] + z_r[yh]*z_s[yh];\
        g_13 = x_r[yh]*x_t[yh] + y_r[yh]*y_t[yh] + z_r[yh]*z_t[yh];\
        g_22 = x_s[yh]*x_s[yh] + y_s[yh]*y_s[yh] + z_s[yh]*z_s[yh];\
        g_23 = x_s[yh]*x_t[yh] + y_s[yh]*y_t[yh] + z_s[yh]*z_t[yh];\
        g_33 = x_t[yh]*x_t[yh] + y_t[yh]*y_t[yh] + z_t[yh]*z_t[yh];\
\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        g22 = FDIV(1.0,g) * (g_11*g_33 - g_13*g_13);\
        coef[MET_S][k][j][i] = sqrt(g) * g22;\
\
    }\
\
\
    for (j = jbeg; j < jend; j++)\
    {\
        x_r[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j][i-1].x + sp[-2] * coors[k][j][i-2].x;\
        y_r[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j][i-1].y + sp[-2] * coors[k][j][i-2].y;\
        z_r[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j][i-1].z + sp[-2] * coors[k][j][i-2].z;\
        x_s[c] = si[-1] * coors[k][j-1][i].x + si[1] * coors[k][j+1][i].x;\
        y_s[c] = si[-1] * coors[k][j-1][i].y + si[1] * coors[k][j+1][i].y;\
        z_s[c] = si[-1] * coors[k][j-1][i].z + si[1] * coors[k][j+1][i].z;\
        x_t[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k-1][j][i].x + sp[-2] * coors[k-2][j][i].x;\
        y_t[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k-1][j][i].y + sp[-2] * coors[k-2][j][i].y;\
        z_t[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k-1][j][i].z + sp[-2] * coors[k-2][j][i].z;\
        g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
        g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
        g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
        g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
        g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
        g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        jac[k][j][i] = sqrt(g);\
        vec__x_r[k][j][i] = x_r[c];\
        vec__y_r[k][j][i] = y_r[c];\
        vec__z_r[k][j][i] = z_r[c];\
        vec__x_s[k][j][i] = x_s[c];\
        vec__y_s[k][j][i] = y_s[c];\
        vec__z_s[k][j][i] = z_s[c];\
        vec__x_t[k][j][i] = x_t[c];\
        vec__y_t[k][j][i] = y_t[c];\
        vec__z_t[k][j][i] = z_t[c];\
        g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
        g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
        g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
        coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
        coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
        coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
        coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
        coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
        coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
\
    }\
\
}\
\
\
if ((xs + xm == ngx && !per_x) && (ys + ym == ngy && !per_y)) {	\
    i = ngx - 1;\
    j = ngy - 1;\
    for (k = kbeg; k < zs + zm; k++)\
    {\
        pt0[0] = .5 * coors[k-1][j][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k-1][j][i-1].x + .5 * coors[k][j][i-1].x;\
        pt2[0] = .5 * coors[k-1][j][i-2].x + .5 * coors[k][j][i-2].x;\
        pt0[1] = .5 * coors[k-1][j][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k-1][j][i-1].y + .5 * coors[k][j][i-1].y;\
        pt2[1] = .5 * coors[k-1][j][i-2].y + .5 * coors[k][j][i-2].y;\
        pt0[2] = .5 * coors[k-1][j][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k-1][j][i-1].z + .5 * coors[k][j][i-1].z;\
        pt2[2] = .5 * coors[k-1][j][i-2].z + .5 * coors[k][j][i-2].z;\
        x_r[zh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
        y_r[zh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
        z_r[zh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
        pt0[0] = .5 * coors[k-1][j][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k-1][j-1][i].x + .5 * coors[k][j-1][i].x;\
        pt2[0] = .5 * coors[k-1][j-2][i].x + .5 * coors[k][j-2][i].x;\
        pt0[1] = .5 * coors[k-1][j][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k-1][j-1][i].y + .5 * coors[k][j-1][i].y;\
        pt2[1] = .5 * coors[k-1][j-2][i].y + .5 * coors[k][j-2][i].y;\
        pt0[2] = .5 * coors[k-1][j][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k-1][j-1][i].z + .5 * coors[k][j-1][i].z;\
        pt2[2] = .5 * coors[k-1][j-2][i].z + .5 * coors[k][j-2][i].z;\
        x_s[zh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
        y_s[zh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
        z_s[zh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
        x_t[zh] = coors[k][j][i].x - coors[k-1][j][i].x;\
        y_t[zh] = coors[k][j][i].y - coors[k-1][j][i].y;\
        z_t[zh] = coors[k][j][i].z - coors[k-1][j][i].z;\
        g_11 = x_r[zh]*x_r[zh] + y_r[zh]*y_r[zh] + z_r[zh]*z_r[zh];\
        g_12 = x_r[zh]*x_s[zh] + y_r[zh]*y_s[zh] + z_r[zh]*z_s[zh];\
        g_13 = x_r[zh]*x_t[zh] + y_r[zh]*y_t[zh] + z_r[zh]*z_t[zh];\
        g_22 = x_s[zh]*x_s[zh] + y_s[zh]*y_s[zh] + z_s[zh]*z_s[zh];\
        g_23 = x_s[zh]*x_t[zh] + y_s[zh]*y_t[zh] + z_s[zh]*z_t[zh];\
        g_33 = x_t[zh]*x_t[zh] + y_t[zh]*y_t[zh] + z_t[zh]*z_t[zh];\
\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        g33 = FDIV(1.0,g) * (g_11*g_22 - g_12*g_12);\
        coef[MET_B][k][j][i] = sqrt(g) * g33;\
\
    }\
\
\
    for (k = kbeg; k < kend; k++)\
    {\
        x_r[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j][i-1].x + sp[-2] * coors[k][j][i-2].x;\
        y_r[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j][i-1].y + sp[-2] * coors[k][j][i-2].y;\
        z_r[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j][i-1].z + sp[-2] * coors[k][j][i-2].z;\
        x_s[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j-1][i].x + sp[-2] * coors[k][j-2][i].x;\
        y_s[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j-1][i].y + sp[-2] * coors[k][j-2][i].y;\
        z_s[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j-1][i].z + sp[-2] * coors[k][j-2][i].z;\
        x_t[c] = si[-1] * coors[k-1][j][i].x + si[1] * coors[k+1][j][i].x;\
        y_t[c] = si[-1] * coors[k-1][j][i].y + si[1] * coors[k+1][j][i].y;\
        z_t[c] = si[-1] * coors[k-1][j][i].z + si[1] * coors[k+1][j][i].z;\
        g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
        g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
        g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
        g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
        g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
        g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        jac[k][j][i] = sqrt(g);\
        vec__x_r[k][j][i] = x_r[c];\
        vec__y_r[k][j][i] = y_r[c];\
        vec__z_r[k][j][i] = z_r[c];\
        vec__x_s[k][j][i] = x_s[c];\
        vec__y_s[k][j][i] = y_s[c];\
        vec__z_s[k][j][i] = z_s[c];\
        vec__x_t[k][j][i] = x_t[c];\
        vec__y_t[k][j][i] = y_t[c];\
        vec__z_t[k][j][i] = z_t[c];\
        g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
        g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
        g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
        coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
        coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
        coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
        coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
        coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
        coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
        coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
\
    }\
\
}\
\
\
if ((ys + ym == ngy && !per_y) && (zs == 0 && !per_z)) {	\
    j = ngy - 1;\
    k = 0;\
    for (i = ibeg; i < xs + xm; i++)\
    {\
        x_r[xh] = coors[k][j][i].x - coors[k][j][i-1].x;\
        y_r[xh] = coors[k][j][i].y - coors[k][j][i-1].y;\
        z_r[xh] = coors[k][j][i].z - coors[k][j][i-1].z;\
        pt0[0] = .5 * coors[k][j][i-1].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k][j-1][i-1].x + .5 * coors[k][j-1][i].x;\
        pt2[0] = .5 * coors[k][j-2][i-1].x + .5 * coors[k][j-2][i].x;\
        pt0[1] = .5 * coors[k][j][i-1].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k][j-1][i-1].y + .5 * coors[k][j-1][i].y;\
        pt2[1] = .5 * coors[k][j-2][i-1].y + .5 * coors[k][j-2][i].y;\
        pt0[2] = .5 * coors[k][j][i-1].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k][j-1][i-1].z + .5 * coors[k][j-1][i].z;\
        pt2[2] = .5 * coors[k][j-2][i-1].z + .5 * coors[k][j-2][i].z;\
        x_s[xh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
        y_s[xh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
        z_s[xh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
        pt0[0] = .5 * coors[k][j][i-1].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k+1][j][i-1].x + .5 * coors[k+1][j][i].x;\
        pt2[0] = .5 * coors[k+2][j][i-1].x + .5 * coors[k+2][j][i].x;\
        pt0[1] = .5 * coors[k][j][i-1].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k+1][j][i-1].y + .5 * coors[k+1][j][i].y;\
        pt2[1] = .5 * coors[k+2][j][i-1].y + .5 * coors[k+2][j][i].y;\
        pt0[2] = .5 * coors[k][j][i-1].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k+1][j][i-1].z + .5 * coors[k+1][j][i].z;\
        pt2[2] = .5 * coors[k+2][j][i-1].z + .5 * coors[k+2][j][i].z;\
        x_t[xh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
        y_t[xh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
        z_t[xh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
        g_11 = x_r[xh]*x_r[xh] + y_r[xh]*y_r[xh] + z_r[xh]*z_r[xh];\
        g_12 = x_r[xh]*x_s[xh] + y_r[xh]*y_s[xh] + z_r[xh]*z_s[xh];\
        g_13 = x_r[xh]*x_t[xh] + y_r[xh]*y_t[xh] + z_r[xh]*z_t[xh];\
        g_22 = x_s[xh]*x_s[xh] + y_s[xh]*y_s[xh] + z_s[xh]*z_s[xh];\
        g_23 = x_s[xh]*x_t[xh] + y_s[xh]*y_t[xh] + z_s[xh]*z_t[xh];\
        g_33 = x_t[xh]*x_t[xh] + y_t[xh]*y_t[xh] + z_t[xh]*z_t[xh];\
\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        g11 = FDIV(1.0,g) * (g_22*g_33 - g_23*g_23);\
        coef[MET_W][k][j][i] = sqrt(g) * g11;\
\
    }\
\
\
    for (i = ibeg; i < iend; i++)\
    {\
        x_r[c] = si[-1] * coors[k][j][i-1].x + si[1] * coors[k][j][i+1].x;\
        y_r[c] = si[-1] * coors[k][j][i-1].y + si[1] * coors[k][j][i+1].y;\
        z_r[c] = si[-1] * coors[k][j][i-1].z + si[1] * coors[k][j][i+1].z;\
        x_s[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j-1][i].x + sp[-2] * coors[k][j-2][i].x;\
        y_s[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j-1][i].y + sp[-2] * coors[k][j-2][i].y;\
        z_s[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j-1][i].z + sp[-2] * coors[k][j-2][i].z;\
        x_t[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k+1][j][i].x + sm[2] * coors[k+2][j][i].x;\
        y_t[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k+1][j][i].y + sm[2] * coors[k+2][j][i].y;\
        z_t[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k+1][j][i].z + sm[2] * coors[k+2][j][i].z;\
        g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
        g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
        g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
        g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
        g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
        g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        jac[k][j][i] = sqrt(g);\
        vec__x_r[k][j][i] = x_r[c];\
        vec__y_r[k][j][i] = y_r[c];\
        vec__z_r[k][j][i] = z_r[c];\
        vec__x_s[k][j][i] = x_s[c];\
        vec__y_s[k][j][i] = y_s[c];\
        vec__z_s[k][j][i] = z_s[c];\
        vec__x_t[k][j][i] = x_t[c];\
        vec__y_t[k][j][i] = y_t[c];\
        vec__z_t[k][j][i] = z_t[c];\
        g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
        g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
        g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
        coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
        coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
        coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
        coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
        coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
        coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
        coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
        coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
        coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
\
    }\
\
}\
\
\
if ((xs + xm == ngx && !per_x) && (zs == 0 && !per_z)) {	\
    i = ngx - 1;\
    k = 0;\
    for (j = jbeg; j < ys + ym; j++)\
    {\
        pt0[0] = .5 * coors[k][j-1][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k][j-1][i-1].x + .5 * coors[k][j][i-1].x;\
        pt2[0] = .5 * coors[k][j-1][i-2].x + .5 * coors[k][j][i-2].x;\
        pt0[1] = .5 * coors[k][j-1][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k][j-1][i-1].y + .5 * coors[k][j][i-1].y;\
        pt2[1] = .5 * coors[k][j-1][i-2].y + .5 * coors[k][j][i-2].y;\
        pt0[2] = .5 * coors[k][j-1][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k][j-1][i-1].z + .5 * coors[k][j][i-1].z;\
        pt2[2] = .5 * coors[k][j-1][i-2].z + .5 * coors[k][j][i-2].z;\
        x_r[yh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
        y_r[yh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
        z_r[yh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
        x_s[yh] = coors[k][j][i].x - coors[k][j-1][i].x;\
        y_s[yh] = coors[k][j][i].y - coors[k][j-1][i].y;\
        z_s[yh] = coors[k][j][i].z - coors[k][j-1][i].z;\
        pt0[0] = .5 * coors[k][j-1][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k+1][j-1][i].x + .5 * coors[k+1][j][i].x;\
        pt2[0] = .5 * coors[k+2][j-1][i].x + .5 * coors[k+2][j][i].x;\
        pt0[1] = .5 * coors[k][j-1][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k+1][j-1][i].y + .5 * coors[k+1][j][i].y;\
        pt2[1] = .5 * coors[k+2][j-1][i].y + .5 * coors[k+2][j][i].y;\
        pt0[2] = .5 * coors[k][j-1][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k+1][j-1][i].z + .5 * coors[k+1][j][i].z;\
        pt2[2] = .5 * coors[k+2][j-1][i].z + .5 * coors[k+2][j][i].z;\
        x_t[yh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
        y_t[yh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
        z_t[yh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
        g_11 = x_r[yh]*x_r[yh] + y_r[yh]*y_r[yh] + z_r[yh]*z_r[yh];\
        g_12 = x_r[yh]*x_s[yh] + y_r[yh]*y_s[yh] + z_r[yh]*z_s[yh];\
        g_13 = x_r[yh]*x_t[yh] + y_r[yh]*y_t[yh] + z_r[yh]*z_t[yh];\
        g_22 = x_s[yh]*x_s[yh] + y_s[yh]*y_s[yh] + z_s[yh]*z_s[yh];\
        g_23 = x_s[yh]*x_t[yh] + y_s[yh]*y_t[yh] + z_s[yh]*z_t[yh];\
        g_33 = x_t[yh]*x_t[yh] + y_t[yh]*y_t[yh] + z_t[yh]*z_t[yh];\
\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        g22 = FDIV(1.0,g) * (g_11*g_33 - g_13*g_13);\
        coef[MET_S][k][j][i] = sqrt(g) * g22;\
\
    }\
\
\
    for (j = jbeg; j < jend; j++)\
    {\
        x_r[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j][i-1].x + sp[-2] * coors[k][j][i-2].x;\
        y_r[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j][i-1].y + sp[-2] * coors[k][j][i-2].y;\
        z_r[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j][i-1].z + sp[-2] * coors[k][j][i-2].z;\
        x_s[c] = si[-1] * coors[k][j-1][i].x + si[1] * coors[k][j+1][i].x;\
        y_s[c] = si[-1] * coors[k][j-1][i].y + si[1] * coors[k][j+1][i].y;\
        z_s[c] = si[-1] * coors[k][j-1][i].z + si[1] * coors[k][j+1][i].z;\
        x_t[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k+1][j][i].x + sm[2] * coors[k+2][j][i].x;\
        y_t[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k+1][j][i].y + sm[2] * coors[k+2][j][i].y;\
        z_t[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k+1][j][i].z + sm[2] * coors[k+2][j][i].z;\
        g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
        g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
        g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
        g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
        g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
        g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        jac[k][j][i] = sqrt(g);\
        vec__x_r[k][j][i] = x_r[c];\
        vec__y_r[k][j][i] = y_r[c];\
        vec__z_r[k][j][i] = z_r[c];\
        vec__x_s[k][j][i] = x_s[c];\
        vec__y_s[k][j][i] = y_s[c];\
        vec__z_s[k][j][i] = z_s[c];\
        vec__x_t[k][j][i] = x_t[c];\
        vec__y_t[k][j][i] = y_t[c];\
        vec__z_t[k][j][i] = z_t[c];\
        g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
        g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
        g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
        coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
        coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
        coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
        coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
        coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
        coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
        coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
        coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
        coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
        coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
\
    }\
\
}\
\
\
if ((xs + xm == ngx && !per_x) && (ys == 0 && !per_y)) {	\
    i = ngx - 1;\
    j = 0;\
    for (k = kbeg; k < zs + zm; k++)\
    {\
        pt0[0] = .5 * coors[k-1][j][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k-1][j][i-1].x + .5 * coors[k][j][i-1].x;\
        pt2[0] = .5 * coors[k-1][j][i-2].x + .5 * coors[k][j][i-2].x;\
        pt0[1] = .5 * coors[k-1][j][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k-1][j][i-1].y + .5 * coors[k][j][i-1].y;\
        pt2[1] = .5 * coors[k-1][j][i-2].y + .5 * coors[k][j][i-2].y;\
        pt0[2] = .5 * coors[k-1][j][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k-1][j][i-1].z + .5 * coors[k][j][i-1].z;\
        pt2[2] = .5 * coors[k-1][j][i-2].z + .5 * coors[k][j][i-2].z;\
        x_r[zh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
        y_r[zh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
        z_r[zh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
        pt0[0] = .5 * coors[k-1][j][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k-1][j+1][i].x + .5 * coors[k][j+1][i].x;\
        pt2[0] = .5 * coors[k-1][j+2][i].x + .5 * coors[k][j+2][i].x;\
        pt0[1] = .5 * coors[k-1][j][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k-1][j+1][i].y + .5 * coors[k][j+1][i].y;\
        pt2[1] = .5 * coors[k-1][j+2][i].y + .5 * coors[k][j+2][i].y;\
        pt0[2] = .5 * coors[k-1][j][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k-1][j+1][i].z + .5 * coors[k][j+1][i].z;\
        pt2[2] = .5 * coors[k-1][j+2][i].z + .5 * coors[k][j+2][i].z;\
        x_s[zh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
        y_s[zh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
        z_s[zh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
        x_t[zh] = coors[k][j][i].x - coors[k-1][j][i].x;\
        y_t[zh] = coors[k][j][i].y - coors[k-1][j][i].y;\
        z_t[zh] = coors[k][j][i].z - coors[k-1][j][i].z;\
        g_11 = x_r[zh]*x_r[zh] + y_r[zh]*y_r[zh] + z_r[zh]*z_r[zh];\
        g_12 = x_r[zh]*x_s[zh] + y_r[zh]*y_s[zh] + z_r[zh]*z_s[zh];\
        g_13 = x_r[zh]*x_t[zh] + y_r[zh]*y_t[zh] + z_r[zh]*z_t[zh];\
        g_22 = x_s[zh]*x_s[zh] + y_s[zh]*y_s[zh] + z_s[zh]*z_s[zh];\
        g_23 = x_s[zh]*x_t[zh] + y_s[zh]*y_t[zh] + z_s[zh]*z_t[zh];\
        g_33 = x_t[zh]*x_t[zh] + y_t[zh]*y_t[zh] + z_t[zh]*z_t[zh];\
\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        g33 = FDIV(1.0,g) * (g_11*g_22 - g_12*g_12);\
        coef[MET_B][k][j][i] = sqrt(g) * g33;\
\
    }\
\
\
    for (k = kbeg; k < kend; k++)\
    {\
        x_r[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j][i-1].x + sp[-2] * coors[k][j][i-2].x;\
        y_r[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j][i-1].y + sp[-2] * coors[k][j][i-2].y;\
        z_r[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j][i-1].z + sp[-2] * coors[k][j][i-2].z;\
        x_s[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j+1][i].x + sm[2] * coors[k][j+2][i].x;\
        y_s[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j+1][i].y + sm[2] * coors[k][j+2][i].y;\
        z_s[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j+1][i].z + sm[2] * coors[k][j+2][i].z;\
        x_t[c] = si[-1] * coors[k-1][j][i].x + si[1] * coors[k+1][j][i].x;\
        y_t[c] = si[-1] * coors[k-1][j][i].y + si[1] * coors[k+1][j][i].y;\
        z_t[c] = si[-1] * coors[k-1][j][i].z + si[1] * coors[k+1][j][i].z;\
        g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
        g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
        g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
        g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
        g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
        g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        jac[k][j][i] = sqrt(g);\
        vec__x_r[k][j][i] = x_r[c];\
        vec__y_r[k][j][i] = y_r[c];\
        vec__z_r[k][j][i] = z_r[c];\
        vec__x_s[k][j][i] = x_s[c];\
        vec__y_s[k][j][i] = y_s[c];\
        vec__z_s[k][j][i] = z_s[c];\
        vec__x_t[k][j][i] = x_t[c];\
        vec__y_t[k][j][i] = y_t[c];\
        vec__z_t[k][j][i] = z_t[c];\
        g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
        g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
        g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
        coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
        coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
        coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
        coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
        coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
        coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
        coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
        coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
        coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
\
    }\
\
}\
\
\
if ((ys == 0 && !per_y) && ((zs + zm == ngz && !per_z) && !per_z)) {	\
    j = 0;\
    k = ngz - 1;\
    for (i = ibeg; i < xs + xm; i++)\
    {\
        x_r[xh] = coors[k][j][i].x - coors[k][j][i-1].x;\
        y_r[xh] = coors[k][j][i].y - coors[k][j][i-1].y;\
        z_r[xh] = coors[k][j][i].z - coors[k][j][i-1].z;\
        pt0[0] = .5 * coors[k][j][i-1].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k][j+1][i-1].x + .5 * coors[k][j+1][i].x;\
        pt2[0] = .5 * coors[k][j+2][i-1].x + .5 * coors[k][j+2][i].x;\
        pt0[1] = .5 * coors[k][j][i-1].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k][j+1][i-1].y + .5 * coors[k][j+1][i].y;\
        pt2[1] = .5 * coors[k][j+2][i-1].y + .5 * coors[k][j+2][i].y;\
        pt0[2] = .5 * coors[k][j][i-1].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k][j+1][i-1].z + .5 * coors[k][j+1][i].z;\
        pt2[2] = .5 * coors[k][j+2][i-1].z + .5 * coors[k][j+2][i].z;\
        x_s[xh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
        y_s[xh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
        z_s[xh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
        pt0[0] = .5 * coors[k][j][i-1].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k-1][j][i-1].x + .5 * coors[k-1][j][i].x;\
        pt2[0] = .5 * coors[k-2][j][i-1].x + .5 * coors[k-2][j][i].x;\
        pt0[1] = .5 * coors[k][j][i-1].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k-1][j][i-1].y + .5 * coors[k-1][j][i].y;\
        pt2[1] = .5 * coors[k-2][j][i-1].y + .5 * coors[k-2][j][i].y;\
        pt0[2] = .5 * coors[k][j][i-1].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k-1][j][i-1].z + .5 * coors[k-1][j][i].z;\
        pt2[2] = .5 * coors[k-2][j][i-1].z + .5 * coors[k-2][j][i].z;\
        x_t[xh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
        y_t[xh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
        z_t[xh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
        g_11 = x_r[xh]*x_r[xh] + y_r[xh]*y_r[xh] + z_r[xh]*z_r[xh];\
        g_12 = x_r[xh]*x_s[xh] + y_r[xh]*y_s[xh] + z_r[xh]*z_s[xh];\
        g_13 = x_r[xh]*x_t[xh] + y_r[xh]*y_t[xh] + z_r[xh]*z_t[xh];\
        g_22 = x_s[xh]*x_s[xh] + y_s[xh]*y_s[xh] + z_s[xh]*z_s[xh];\
        g_23 = x_s[xh]*x_t[xh] + y_s[xh]*y_t[xh] + z_s[xh]*z_t[xh];\
        g_33 = x_t[xh]*x_t[xh] + y_t[xh]*y_t[xh] + z_t[xh]*z_t[xh];\
\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        g11 = FDIV(1.0,g) * (g_22*g_33 - g_23*g_23);\
        coef[MET_W][k][j][i] = sqrt(g) * g11;\
\
    }\
\
\
    for (i = ibeg; i < iend; i++)\
    {\
        x_r[c] = si[-1] * coors[k][j][i-1].x + si[1] * coors[k][j][i+1].x;\
        y_r[c] = si[-1] * coors[k][j][i-1].y + si[1] * coors[k][j][i+1].y;\
        z_r[c] = si[-1] * coors[k][j][i-1].z + si[1] * coors[k][j][i+1].z;\
        x_s[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j+1][i].x + sm[2] * coors[k][j+2][i].x;\
        y_s[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j+1][i].y + sm[2] * coors[k][j+2][i].y;\
        z_s[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j+1][i].z + sm[2] * coors[k][j+2][i].z;\
        x_t[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k-1][j][i].x + sp[-2] * coors[k-2][j][i].x;\
        y_t[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k-1][j][i].y + sp[-2] * coors[k-2][j][i].y;\
        z_t[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k-1][j][i].z + sp[-2] * coors[k-2][j][i].z;\
        g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
        g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
        g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
        g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
        g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
        g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        jac[k][j][i] = sqrt(g);\
        vec__x_r[k][j][i] = x_r[c];\
        vec__y_r[k][j][i] = y_r[c];\
        vec__z_r[k][j][i] = z_r[c];\
        vec__x_s[k][j][i] = x_s[c];\
        vec__y_s[k][j][i] = y_s[c];\
        vec__z_s[k][j][i] = z_s[c];\
        vec__x_t[k][j][i] = x_t[c];\
        vec__y_t[k][j][i] = y_t[c];\
        vec__z_t[k][j][i] = z_t[c];\
        g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
        g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
        g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
        coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
        coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
        coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
        coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
        coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
        coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
        coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
\
    }\
\
}\
\
\
if ((xs == 0 && !per_x) && ((zs + zm == ngz && !per_z) && !per_z)) {	\
    i = 0;\
    k = ngz - 1;\
    for (j = jbeg; j < ys + ym; j++)\
    {\
        pt0[0] = .5 * coors[k][j-1][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k][j-1][i+1].x + .5 * coors[k][j][i+1].x;\
        pt2[0] = .5 * coors[k][j-1][i+2].x + .5 * coors[k][j][i+2].x;\
        pt0[1] = .5 * coors[k][j-1][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k][j-1][i+1].y + .5 * coors[k][j][i+1].y;\
        pt2[1] = .5 * coors[k][j-1][i+2].y + .5 * coors[k][j][i+2].y;\
        pt0[2] = .5 * coors[k][j-1][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k][j-1][i+1].z + .5 * coors[k][j][i+1].z;\
        pt2[2] = .5 * coors[k][j-1][i+2].z + .5 * coors[k][j][i+2].z;\
        x_r[yh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
        y_r[yh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
        z_r[yh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
        x_s[yh] = coors[k][j][i].x - coors[k][j-1][i].x;\
        y_s[yh] = coors[k][j][i].y - coors[k][j-1][i].y;\
        z_s[yh] = coors[k][j][i].z - coors[k][j-1][i].z;\
        pt0[0] = .5 * coors[k][j-1][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k-1][j-1][i].x + .5 * coors[k-1][j][i].x;\
        pt2[0] = .5 * coors[k-2][j-1][i].x + .5 * coors[k-2][j][i].x;\
        pt0[1] = .5 * coors[k][j-1][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k-1][j-1][i].y + .5 * coors[k-1][j][i].y;\
        pt2[1] = .5 * coors[k-2][j-1][i].y + .5 * coors[k-2][j][i].y;\
        pt0[2] = .5 * coors[k][j-1][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k-1][j-1][i].z + .5 * coors[k-1][j][i].z;\
        pt2[2] = .5 * coors[k-2][j-1][i].z + .5 * coors[k-2][j][i].z;\
        x_t[yh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
        y_t[yh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
        z_t[yh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
        g_11 = x_r[yh]*x_r[yh] + y_r[yh]*y_r[yh] + z_r[yh]*z_r[yh];\
        g_12 = x_r[yh]*x_s[yh] + y_r[yh]*y_s[yh] + z_r[yh]*z_s[yh];\
        g_13 = x_r[yh]*x_t[yh] + y_r[yh]*y_t[yh] + z_r[yh]*z_t[yh];\
        g_22 = x_s[yh]*x_s[yh] + y_s[yh]*y_s[yh] + z_s[yh]*z_s[yh];\
        g_23 = x_s[yh]*x_t[yh] + y_s[yh]*y_t[yh] + z_s[yh]*z_t[yh];\
        g_33 = x_t[yh]*x_t[yh] + y_t[yh]*y_t[yh] + z_t[yh]*z_t[yh];\
\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        g22 = FDIV(1.0,g) * (g_11*g_33 - g_13*g_13);\
        coef[MET_S][k][j][i] = sqrt(g) * g22;\
\
    }\
\
\
    for (j = jbeg; j < jend; j++)\
    {\
        x_r[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j][i+1].x + sm[2] * coors[k][j][i+2].x;\
        y_r[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j][i+1].y + sm[2] * coors[k][j][i+2].y;\
        z_r[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j][i+1].z + sm[2] * coors[k][j][i+2].z;\
        x_s[c] = si[-1] * coors[k][j-1][i].x + si[1] * coors[k][j+1][i].x;\
        y_s[c] = si[-1] * coors[k][j-1][i].y + si[1] * coors[k][j+1][i].y;\
        z_s[c] = si[-1] * coors[k][j-1][i].z + si[1] * coors[k][j+1][i].z;\
        x_t[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k-1][j][i].x + sp[-2] * coors[k-2][j][i].x;\
        y_t[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k-1][j][i].y + sp[-2] * coors[k-2][j][i].y;\
        z_t[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k-1][j][i].z + sp[-2] * coors[k-2][j][i].z;\
        g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
        g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
        g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
        g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
        g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
        g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        jac[k][j][i] = sqrt(g);\
        vec__x_r[k][j][i] = x_r[c];\
        vec__y_r[k][j][i] = y_r[c];\
        vec__z_r[k][j][i] = z_r[c];\
        vec__x_s[k][j][i] = x_s[c];\
        vec__y_s[k][j][i] = y_s[c];\
        vec__z_s[k][j][i] = z_s[c];\
        vec__x_t[k][j][i] = x_t[c];\
        vec__y_t[k][j][i] = y_t[c];\
        vec__z_t[k][j][i] = z_t[c];\
        g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
        g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
        g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
        coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
        coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
        coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
        coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
        coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
        coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
\
    }\
\
}\
\
\
if ((xs == 0 && !per_x) && (ys + ym == ngy && !per_y)) {	\
    i = 0;\
    j = ngy - 1;\
    for (k = kbeg; k < zs + zm; k++)\
    {\
        pt0[0] = .5 * coors[k-1][j][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k-1][j][i+1].x + .5 * coors[k][j][i+1].x;\
        pt2[0] = .5 * coors[k-1][j][i+2].x + .5 * coors[k][j][i+2].x;\
        pt0[1] = .5 * coors[k-1][j][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k-1][j][i+1].y + .5 * coors[k][j][i+1].y;\
        pt2[1] = .5 * coors[k-1][j][i+2].y + .5 * coors[k][j][i+2].y;\
        pt0[2] = .5 * coors[k-1][j][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k-1][j][i+1].z + .5 * coors[k][j][i+1].z;\
        pt2[2] = .5 * coors[k-1][j][i+2].z + .5 * coors[k][j][i+2].z;\
        x_r[zh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
        y_r[zh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
        z_r[zh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
        pt0[0] = .5 * coors[k-1][j][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k-1][j-1][i].x + .5 * coors[k][j-1][i].x;\
        pt2[0] = .5 * coors[k-1][j-2][i].x + .5 * coors[k][j-2][i].x;\
        pt0[1] = .5 * coors[k-1][j][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k-1][j-1][i].y + .5 * coors[k][j-1][i].y;\
        pt2[1] = .5 * coors[k-1][j-2][i].y + .5 * coors[k][j-2][i].y;\
        pt0[2] = .5 * coors[k-1][j][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k-1][j-1][i].z + .5 * coors[k][j-1][i].z;\
        pt2[2] = .5 * coors[k-1][j-2][i].z + .5 * coors[k][j-2][i].z;\
        x_s[zh] = sp[0] * pt0[0] + sp[-1] * pt1[0] + sp[-2] * pt2[0];\
        y_s[zh] = sp[0] * pt0[1] + sp[-1] * pt1[1] + sp[-2] * pt2[1];\
        z_s[zh] = sp[0] * pt0[2] + sp[-1] * pt1[2] + sp[-2] * pt2[2];\
        x_t[zh] = coors[k][j][i].x - coors[k-1][j][i].x;\
        y_t[zh] = coors[k][j][i].y - coors[k-1][j][i].y;\
        z_t[zh] = coors[k][j][i].z - coors[k-1][j][i].z;\
        g_11 = x_r[zh]*x_r[zh] + y_r[zh]*y_r[zh] + z_r[zh]*z_r[zh];\
        g_12 = x_r[zh]*x_s[zh] + y_r[zh]*y_s[zh] + z_r[zh]*z_s[zh];\
        g_13 = x_r[zh]*x_t[zh] + y_r[zh]*y_t[zh] + z_r[zh]*z_t[zh];\
        g_22 = x_s[zh]*x_s[zh] + y_s[zh]*y_s[zh] + z_s[zh]*z_s[zh];\
        g_23 = x_s[zh]*x_t[zh] + y_s[zh]*y_t[zh] + z_s[zh]*z_t[zh];\
        g_33 = x_t[zh]*x_t[zh] + y_t[zh]*y_t[zh] + z_t[zh]*z_t[zh];\
\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        g33 = FDIV(1.0,g) * (g_11*g_22 - g_12*g_12);\
        coef[MET_B][k][j][i] = sqrt(g) * g33;\
\
    }\
\
\
    for (k = kbeg; k < kend; k++)\
    {\
        x_r[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j][i+1].x + sm[2] * coors[k][j][i+2].x;\
        y_r[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j][i+1].y + sm[2] * coors[k][j][i+2].y;\
        z_r[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j][i+1].z + sm[2] * coors[k][j][i+2].z;\
        x_s[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j-1][i].x + sp[-2] * coors[k][j-2][i].x;\
        y_s[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j-1][i].y + sp[-2] * coors[k][j-2][i].y;\
        z_s[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j-1][i].z + sp[-2] * coors[k][j-2][i].z;\
        x_t[c] = si[-1] * coors[k-1][j][i].x + si[1] * coors[k+1][j][i].x;\
        y_t[c] = si[-1] * coors[k-1][j][i].y + si[1] * coors[k+1][j][i].y;\
        z_t[c] = si[-1] * coors[k-1][j][i].z + si[1] * coors[k+1][j][i].z;\
        g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
        g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
        g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
        g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
        g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
        g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        jac[k][j][i] = sqrt(g);\
        vec__x_r[k][j][i] = x_r[c];\
        vec__y_r[k][j][i] = y_r[c];\
        vec__z_r[k][j][i] = z_r[c];\
        vec__x_s[k][j][i] = x_s[c];\
        vec__y_s[k][j][i] = y_s[c];\
        vec__z_s[k][j][i] = z_s[c];\
        vec__x_t[k][j][i] = x_t[c];\
        vec__y_t[k][j][i] = y_t[c];\
        vec__z_t[k][j][i] = z_t[c];\
        g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
        g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
        g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
        coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
        coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
        coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
        coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
        coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
        coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
        coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
\
    }\
\
}\
\
\
if ((ys == 0 && !per_y) && (zs == 0 && !per_z)) {	\
    j = 0;\
    k = 0;\
    for (i = ibeg; i < xs + xm; i++)\
    {\
        x_r[xh] = coors[k][j][i].x - coors[k][j][i-1].x;\
        y_r[xh] = coors[k][j][i].y - coors[k][j][i-1].y;\
        z_r[xh] = coors[k][j][i].z - coors[k][j][i-1].z;\
        pt0[0] = .5 * coors[k][j][i-1].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k][j+1][i-1].x + .5 * coors[k][j+1][i].x;\
        pt2[0] = .5 * coors[k][j+2][i-1].x + .5 * coors[k][j+2][i].x;\
        pt0[1] = .5 * coors[k][j][i-1].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k][j+1][i-1].y + .5 * coors[k][j+1][i].y;\
        pt2[1] = .5 * coors[k][j+2][i-1].y + .5 * coors[k][j+2][i].y;\
        pt0[2] = .5 * coors[k][j][i-1].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k][j+1][i-1].z + .5 * coors[k][j+1][i].z;\
        pt2[2] = .5 * coors[k][j+2][i-1].z + .5 * coors[k][j+2][i].z;\
        x_s[xh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
        y_s[xh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
        z_s[xh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
        pt0[0] = .5 * coors[k][j][i-1].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k+1][j][i-1].x + .5 * coors[k+1][j][i].x;\
        pt2[0] = .5 * coors[k+2][j][i-1].x + .5 * coors[k+2][j][i].x;\
        pt0[1] = .5 * coors[k][j][i-1].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k+1][j][i-1].y + .5 * coors[k+1][j][i].y;\
        pt2[1] = .5 * coors[k+2][j][i-1].y + .5 * coors[k+2][j][i].y;\
        pt0[2] = .5 * coors[k][j][i-1].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k+1][j][i-1].z + .5 * coors[k+1][j][i].z;\
        pt2[2] = .5 * coors[k+2][j][i-1].z + .5 * coors[k+2][j][i].z;\
        x_t[xh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
        y_t[xh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
        z_t[xh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
        g_11 = x_r[xh]*x_r[xh] + y_r[xh]*y_r[xh] + z_r[xh]*z_r[xh];\
        g_12 = x_r[xh]*x_s[xh] + y_r[xh]*y_s[xh] + z_r[xh]*z_s[xh];\
        g_13 = x_r[xh]*x_t[xh] + y_r[xh]*y_t[xh] + z_r[xh]*z_t[xh];\
        g_22 = x_s[xh]*x_s[xh] + y_s[xh]*y_s[xh] + z_s[xh]*z_s[xh];\
        g_23 = x_s[xh]*x_t[xh] + y_s[xh]*y_t[xh] + z_s[xh]*z_t[xh];\
        g_33 = x_t[xh]*x_t[xh] + y_t[xh]*y_t[xh] + z_t[xh]*z_t[xh];\
\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        g11 = FDIV(1.0,g) * (g_22*g_33 - g_23*g_23);\
        coef[MET_W][k][j][i] = sqrt(g) * g11;\
\
    }\
\
\
    for (i = ibeg; i < iend; i++)\
    {\
        x_r[c] = si[-1] * coors[k][j][i-1].x + si[1] * coors[k][j][i+1].x;\
        y_r[c] = si[-1] * coors[k][j][i-1].y + si[1] * coors[k][j][i+1].y;\
        z_r[c] = si[-1] * coors[k][j][i-1].z + si[1] * coors[k][j][i+1].z;\
        x_s[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j+1][i].x + sm[2] * coors[k][j+2][i].x;\
        y_s[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j+1][i].y + sm[2] * coors[k][j+2][i].y;\
        z_s[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j+1][i].z + sm[2] * coors[k][j+2][i].z;\
        x_t[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k+1][j][i].x + sm[2] * coors[k+2][j][i].x;\
        y_t[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k+1][j][i].y + sm[2] * coors[k+2][j][i].y;\
        z_t[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k+1][j][i].z + sm[2] * coors[k+2][j][i].z;\
        g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
        g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
        g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
        g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
        g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
        g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        jac[k][j][i] = sqrt(g);\
        vec__x_r[k][j][i] = x_r[c];\
        vec__y_r[k][j][i] = y_r[c];\
        vec__z_r[k][j][i] = z_r[c];\
        vec__x_s[k][j][i] = x_s[c];\
        vec__y_s[k][j][i] = y_s[c];\
        vec__z_s[k][j][i] = z_s[c];\
        vec__x_t[k][j][i] = x_t[c];\
        vec__y_t[k][j][i] = y_t[c];\
        vec__z_t[k][j][i] = z_t[c];\
        g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
        g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
        g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
        coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
        coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
        coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
        coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
        coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
        coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
        coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
        coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
        coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
        coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
        coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
\
    }\
\
}\
\
\
if ((xs == 0 && !per_x) && (zs == 0 && !per_z)) {	\
    i = 0;\
    k = 0;\
    for (j = jbeg; j < ys + ym; j++)\
    {\
        pt0[0] = .5 * coors[k][j-1][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k][j-1][i+1].x + .5 * coors[k][j][i+1].x;\
        pt2[0] = .5 * coors[k][j-1][i+2].x + .5 * coors[k][j][i+2].x;\
        pt0[1] = .5 * coors[k][j-1][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k][j-1][i+1].y + .5 * coors[k][j][i+1].y;\
        pt2[1] = .5 * coors[k][j-1][i+2].y + .5 * coors[k][j][i+2].y;\
        pt0[2] = .5 * coors[k][j-1][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k][j-1][i+1].z + .5 * coors[k][j][i+1].z;\
        pt2[2] = .5 * coors[k][j-1][i+2].z + .5 * coors[k][j][i+2].z;\
        x_r[yh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
        y_r[yh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
        z_r[yh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
        x_s[yh] = coors[k][j][i].x - coors[k][j-1][i].x;\
        y_s[yh] = coors[k][j][i].y - coors[k][j-1][i].y;\
        z_s[yh] = coors[k][j][i].z - coors[k][j-1][i].z;\
        pt0[0] = .5 * coors[k][j-1][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k+1][j-1][i].x + .5 * coors[k+1][j][i].x;\
        pt2[0] = .5 * coors[k+2][j-1][i].x + .5 * coors[k+2][j][i].x;\
        pt0[1] = .5 * coors[k][j-1][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k+1][j-1][i].y + .5 * coors[k+1][j][i].y;\
        pt2[1] = .5 * coors[k+2][j-1][i].y + .5 * coors[k+2][j][i].y;\
        pt0[2] = .5 * coors[k][j-1][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k+1][j-1][i].z + .5 * coors[k+1][j][i].z;\
        pt2[2] = .5 * coors[k+2][j-1][i].z + .5 * coors[k+2][j][i].z;\
        x_t[yh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
        y_t[yh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
        z_t[yh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
        g_11 = x_r[yh]*x_r[yh] + y_r[yh]*y_r[yh] + z_r[yh]*z_r[yh];\
        g_12 = x_r[yh]*x_s[yh] + y_r[yh]*y_s[yh] + z_r[yh]*z_s[yh];\
        g_13 = x_r[yh]*x_t[yh] + y_r[yh]*y_t[yh] + z_r[yh]*z_t[yh];\
        g_22 = x_s[yh]*x_s[yh] + y_s[yh]*y_s[yh] + z_s[yh]*z_s[yh];\
        g_23 = x_s[yh]*x_t[yh] + y_s[yh]*y_t[yh] + z_s[yh]*z_t[yh];\
        g_33 = x_t[yh]*x_t[yh] + y_t[yh]*y_t[yh] + z_t[yh]*z_t[yh];\
\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        g22 = FDIV(1.0,g) * (g_11*g_33 - g_13*g_13);\
        coef[MET_S][k][j][i] = sqrt(g) * g22;\
\
    }\
\
\
    for (j = jbeg; j < jend; j++)\
    {\
        x_r[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j][i+1].x + sm[2] * coors[k][j][i+2].x;\
        y_r[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j][i+1].y + sm[2] * coors[k][j][i+2].y;\
        z_r[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j][i+1].z + sm[2] * coors[k][j][i+2].z;\
        x_s[c] = si[-1] * coors[k][j-1][i].x + si[1] * coors[k][j+1][i].x;\
        y_s[c] = si[-1] * coors[k][j-1][i].y + si[1] * coors[k][j+1][i].y;\
        z_s[c] = si[-1] * coors[k][j-1][i].z + si[1] * coors[k][j+1][i].z;\
        x_t[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k+1][j][i].x + sm[2] * coors[k+2][j][i].x;\
        y_t[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k+1][j][i].y + sm[2] * coors[k+2][j][i].y;\
        z_t[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k+1][j][i].z + sm[2] * coors[k+2][j][i].z;\
        g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
        g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
        g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
        g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
        g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
        g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        jac[k][j][i] = sqrt(g);\
        vec__x_r[k][j][i] = x_r[c];\
        vec__y_r[k][j][i] = y_r[c];\
        vec__z_r[k][j][i] = z_r[c];\
        vec__x_s[k][j][i] = x_s[c];\
        vec__y_s[k][j][i] = y_s[c];\
        vec__z_s[k][j][i] = z_s[c];\
        vec__x_t[k][j][i] = x_t[c];\
        vec__y_t[k][j][i] = y_t[c];\
        vec__z_t[k][j][i] = z_t[c];\
        g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
        g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
        g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
        coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
        coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
        coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
        coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
        coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
        coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
        coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
        coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
        coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
        coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
\
    }\
\
}\
\
\
if ((xs == 0 && !per_x) && (ys == 0 && !per_y)) {\
    i = 0;\
    j = 0;\
    for (k = kbeg; k < zs + zm; k++)\
    {\
        pt0[0] = .5 * coors[k-1][j][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k-1][j][i+1].x + .5 * coors[k][j][i+1].x;\
        pt2[0] = .5 * coors[k-1][j][i+2].x + .5 * coors[k][j][i+2].x;\
        pt0[1] = .5 * coors[k-1][j][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k-1][j][i+1].y + .5 * coors[k][j][i+1].y;\
        pt2[1] = .5 * coors[k-1][j][i+2].y + .5 * coors[k][j][i+2].y;\
        pt0[2] = .5 * coors[k-1][j][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k-1][j][i+1].z + .5 * coors[k][j][i+1].z;\
        pt2[2] = .5 * coors[k-1][j][i+2].z + .5 * coors[k][j][i+2].z;\
        x_r[zh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
        y_r[zh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
        z_r[zh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
        pt0[0] = .5 * coors[k-1][j][i].x + .5 * coors[k][j][i].x;\
        pt1[0] = .5 * coors[k-1][j+1][i].x + .5 * coors[k][j+1][i].x;\
        pt2[0] = .5 * coors[k-1][j+2][i].x + .5 * coors[k][j+2][i].x;\
        pt0[1] = .5 * coors[k-1][j][i].y + .5 * coors[k][j][i].y;\
        pt1[1] = .5 * coors[k-1][j+1][i].y + .5 * coors[k][j+1][i].y;\
        pt2[1] = .5 * coors[k-1][j+2][i].y + .5 * coors[k][j+2][i].y;\
        pt0[2] = .5 * coors[k-1][j][i].z + .5 * coors[k][j][i].z;\
        pt1[2] = .5 * coors[k-1][j+1][i].z + .5 * coors[k][j+1][i].z;\
        pt2[2] = .5 * coors[k-1][j+2][i].z + .5 * coors[k][j+2][i].z;\
        x_s[zh] = sm[0] * pt0[0] + sm[1] * pt1[0] + sm[2] * pt2[0];\
        y_s[zh] = sm[0] * pt0[1] + sm[1] * pt1[1] + sm[2] * pt2[1];\
        z_s[zh] = sm[0] * pt0[2] + sm[1] * pt1[2] + sm[2] * pt2[2];\
        x_t[zh] = coors[k][j][i].x - coors[k-1][j][i].x;\
        y_t[zh] = coors[k][j][i].y - coors[k-1][j][i].y;\
        z_t[zh] = coors[k][j][i].z - coors[k-1][j][i].z;\
        g_11 = x_r[zh]*x_r[zh] + y_r[zh]*y_r[zh] + z_r[zh]*z_r[zh];\
        g_12 = x_r[zh]*x_s[zh] + y_r[zh]*y_s[zh] + z_r[zh]*z_s[zh];\
        g_13 = x_r[zh]*x_t[zh] + y_r[zh]*y_t[zh] + z_r[zh]*z_t[zh];\
        g_22 = x_s[zh]*x_s[zh] + y_s[zh]*y_s[zh] + z_s[zh]*z_s[zh];\
        g_23 = x_s[zh]*x_t[zh] + y_s[zh]*y_t[zh] + z_s[zh]*z_t[zh];\
        g_33 = x_t[zh]*x_t[zh] + y_t[zh]*y_t[zh] + z_t[zh]*z_t[zh];\
\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        g33 = FDIV(1.0,g) * (g_11*g_22 - g_12*g_12);\
        coef[MET_B][k][j][i] = sqrt(g) * g33;\
\
    }\
\
\
    for (k = kbeg; k < kend; k++)\
    {\
        x_r[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j][i+1].x + sm[2] * coors[k][j][i+2].x;\
        y_r[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j][i+1].y + sm[2] * coors[k][j][i+2].y;\
        z_r[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j][i+1].z + sm[2] * coors[k][j][i+2].z;\
        x_s[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j+1][i].x + sm[2] * coors[k][j+2][i].x;\
        y_s[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j+1][i].y + sm[2] * coors[k][j+2][i].y;\
        z_s[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j+1][i].z + sm[2] * coors[k][j+2][i].z;\
        x_t[c] = si[-1] * coors[k-1][j][i].x + si[1] * coors[k+1][j][i].x;\
        y_t[c] = si[-1] * coors[k-1][j][i].y + si[1] * coors[k+1][j][i].y;\
        z_t[c] = si[-1] * coors[k-1][j][i].z + si[1] * coors[k+1][j][i].z;\
        g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
        g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
        g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
        g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
        g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
        g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
        g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
        jac[k][j][i] = sqrt(g);\
        vec__x_r[k][j][i] = x_r[c];\
        vec__y_r[k][j][i] = y_r[c];\
        vec__z_r[k][j][i] = z_r[c];\
        vec__x_s[k][j][i] = x_s[c];\
        vec__y_s[k][j][i] = y_s[c];\
        vec__z_s[k][j][i] = z_s[c];\
        vec__x_t[k][j][i] = x_t[c];\
        vec__y_t[k][j][i] = y_t[c];\
        vec__z_t[k][j][i] = z_t[c];\
        g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
        g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
        g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
        coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
        coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
        coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
        coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
        coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
        coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
        coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
        coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
        coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
\
    }\
\
}\
\
\
if ((xs == 0 && !per_x) && (ys == 0 && !per_y) && (zs == 0 && !per_z)) {\
    i = 0;\
    j = 0;\
    k = 0;\
    x_r[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j][i+1].x + sm[2] * coors[k][j][i+2].x;\
    y_r[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j][i+1].y + sm[2] * coors[k][j][i+2].y;\
    z_r[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j][i+1].z + sm[2] * coors[k][j][i+2].z;\
    x_s[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j+1][i].x + sm[2] * coors[k][j+2][i].x;\
    y_s[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j+1][i].y + sm[2] * coors[k][j+2][i].y;\
    z_s[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j+1][i].z + sm[2] * coors[k][j+2][i].z;\
    x_t[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k+1][j][i].x + sm[2] * coors[k+2][j][i].x;\
    y_t[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k+1][j][i].y + sm[2] * coors[k+2][j][i].y;\
    z_t[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k+1][j][i].z + sm[2] * coors[k+2][j][i].z;\
    g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
    g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
    g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
    g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
    g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
    g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
    g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
    jac[k][j][i] = sqrt(g);\
    vec__x_r[k][j][i] = x_r[c];\
    vec__y_r[k][j][i] = y_r[c];\
    vec__z_r[k][j][i] = z_r[c];\
    vec__x_s[k][j][i] = x_s[c];\
    vec__y_s[k][j][i] = y_s[c];\
    vec__z_s[k][j][i] = z_s[c];\
    vec__x_t[k][j][i] = x_t[c];\
    vec__y_t[k][j][i] = y_t[c];\
    vec__z_t[k][j][i] = z_t[c];\
    g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
    g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
    g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
    coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
    coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
    coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
    coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
    coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
    coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
    coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
    coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
    coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
}\
\
\
if ((xs == 0 && !per_x) && (ys == 0 && !per_y) && (zs + zm == ngz && !per_z)) {\
    i = 0;\
    j = 0;\
    k = ngz - 1;\
    x_r[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j][i+1].x + sm[2] * coors[k][j][i+2].x;\
    y_r[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j][i+1].y + sm[2] * coors[k][j][i+2].y;\
    z_r[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j][i+1].z + sm[2] * coors[k][j][i+2].z;\
    x_s[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j+1][i].x + sm[2] * coors[k][j+2][i].x;\
    y_s[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j+1][i].y + sm[2] * coors[k][j+2][i].y;\
    z_s[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j+1][i].z + sm[2] * coors[k][j+2][i].z;\
    x_t[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k-1][j][i].x + sp[-2] * coors[k-2][j][i].x;\
    y_t[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k-1][j][i].y + sp[-2] * coors[k-2][j][i].y;\
    z_t[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k-1][j][i].z + sp[-2] * coors[k-2][j][i].z;\
    g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
    g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
    g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
    g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
    g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
    g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
    g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
    jac[k][j][i] = sqrt(g);\
    vec__x_r[k][j][i] = x_r[c];\
    vec__y_r[k][j][i] = y_r[c];\
    vec__z_r[k][j][i] = z_r[c];\
    vec__x_s[k][j][i] = x_s[c];\
    vec__y_s[k][j][i] = y_s[c];\
    vec__z_s[k][j][i] = z_s[c];\
    vec__x_t[k][j][i] = x_t[c];\
    vec__y_t[k][j][i] = y_t[c];\
    vec__z_t[k][j][i] = z_t[c];\
    g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
    g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
    g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
    coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
    coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
    coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
    coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
    coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
}\
\
\
if ((xs == 0 && !per_x) && (ys + ym == ngy && !per_y) && (zs == 0 && !per_z)) {\
    i = 0;\
    j = ngy - 1;\
    k = 0;\
    x_r[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j][i+1].x + sm[2] * coors[k][j][i+2].x;\
    y_r[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j][i+1].y + sm[2] * coors[k][j][i+2].y;\
    z_r[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j][i+1].z + sm[2] * coors[k][j][i+2].z;\
    x_s[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j-1][i].x + sp[-2] * coors[k][j-2][i].x;\
    y_s[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j-1][i].y + sp[-2] * coors[k][j-2][i].y;\
    z_s[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j-1][i].z + sp[-2] * coors[k][j-2][i].z;\
    x_t[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k+1][j][i].x + sm[2] * coors[k+2][j][i].x;\
    y_t[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k+1][j][i].y + sm[2] * coors[k+2][j][i].y;\
    z_t[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k+1][j][i].z + sm[2] * coors[k+2][j][i].z;\
    g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
    g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
    g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
    g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
    g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
    g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
    g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
    jac[k][j][i] = sqrt(g);\
    vec__x_r[k][j][i] = x_r[c];\
    vec__y_r[k][j][i] = y_r[c];\
    vec__z_r[k][j][i] = z_r[c];\
    vec__x_s[k][j][i] = x_s[c];\
    vec__y_s[k][j][i] = y_s[c];\
    vec__z_s[k][j][i] = z_s[c];\
    vec__x_t[k][j][i] = x_t[c];\
    vec__y_t[k][j][i] = y_t[c];\
    vec__z_t[k][j][i] = z_t[c];\
    g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
    g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
    g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
    coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
    coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
    coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
    coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
    coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
    coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
    coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
}\
\
\
if ((xs == 0 && !per_x) && (ys + ym == ngy && !per_y) && (zs + zm == ngz && !per_z)) {\
    i = 0;\
    j = ngy - 1;\
    k = ngz - 1;\
    x_r[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j][i+1].x + sm[2] * coors[k][j][i+2].x;\
    y_r[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j][i+1].y + sm[2] * coors[k][j][i+2].y;\
    z_r[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j][i+1].z + sm[2] * coors[k][j][i+2].z;\
    x_s[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j-1][i].x + sp[-2] * coors[k][j-2][i].x;\
    y_s[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j-1][i].y + sp[-2] * coors[k][j-2][i].y;\
    z_s[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j-1][i].z + sp[-2] * coors[k][j-2][i].z;\
    x_t[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k-1][j][i].x + sp[-2] * coors[k-2][j][i].x;\
    y_t[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k-1][j][i].y + sp[-2] * coors[k-2][j][i].y;\
    z_t[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k-1][j][i].z + sp[-2] * coors[k-2][j][i].z;\
    g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
    g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
    g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
    g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
    g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
    g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
    g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
    jac[k][j][i] = sqrt(g);\
    vec__x_r[k][j][i] = x_r[c];\
    vec__y_r[k][j][i] = y_r[c];\
    vec__z_r[k][j][i] = z_r[c];\
    vec__x_s[k][j][i] = x_s[c];\
    vec__y_s[k][j][i] = y_s[c];\
    vec__z_s[k][j][i] = z_s[c];\
    vec__x_t[k][j][i] = x_t[c];\
    vec__y_t[k][j][i] = y_t[c];\
    vec__z_t[k][j][i] = z_t[c];\
    g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
    g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
    g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
    coef[MET_SW][k][j][i+1] += .25 * sqrt(g) * g12;\
    coef[MET_WB][k][j][i+1] += .25 * sqrt(g) * g13;\
    coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
}\
\
\
if ((xs + xm == ngx && !per_x) && (ys == 0 && !per_y) && (zs == 0 && !per_z)) {\
    i = ngx - 1;\
    j = 0;\
    k = 0;\
    x_r[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j][i-1].x + sp[-2] * coors[k][j][i-2].x;\
    y_r[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j][i-1].y + sp[-2] * coors[k][j][i-2].y;\
    z_r[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j][i-1].z + sp[-2] * coors[k][j][i-2].z;\
    x_s[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j+1][i].x + sm[2] * coors[k][j+2][i].x;\
    y_s[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j+1][i].y + sm[2] * coors[k][j+2][i].y;\
    z_s[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j+1][i].z + sm[2] * coors[k][j+2][i].z;\
    x_t[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k+1][j][i].x + sm[2] * coors[k+2][j][i].x;\
    y_t[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k+1][j][i].y + sm[2] * coors[k+2][j][i].y;\
    z_t[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k+1][j][i].z + sm[2] * coors[k+2][j][i].z;\
    g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
    g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
    g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
    g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
    g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
    g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
    g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
    jac[k][j][i] = sqrt(g);\
    vec__x_r[k][j][i] = x_r[c];\
    vec__y_r[k][j][i] = y_r[c];\
    vec__z_r[k][j][i] = z_r[c];\
    vec__x_s[k][j][i] = x_s[c];\
    vec__y_s[k][j][i] = y_s[c];\
    vec__z_s[k][j][i] = z_s[c];\
    vec__x_t[k][j][i] = x_t[c];\
    vec__y_t[k][j][i] = y_t[c];\
    vec__z_t[k][j][i] = z_t[c];\
    g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
    g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
    g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
    coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
    coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
    coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
    coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
    coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
    coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
    coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
    coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
    coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
}\
\
\
if ((xs + xm == ngx && !per_x) && (ys == 0 && !per_y) && (zs + zm == ngz && !per_z)) {\
    i = ngx - 1;\
    j = 0;\
    k = ngz - 1;\
    x_r[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j][i-1].x + sp[-2] * coors[k][j][i-2].x;\
    y_r[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j][i-1].y + sp[-2] * coors[k][j][i-2].y;\
    z_r[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j][i-1].z + sp[-2] * coors[k][j][i-2].z;\
    x_s[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k][j+1][i].x + sm[2] * coors[k][j+2][i].x;\
    y_s[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k][j+1][i].y + sm[2] * coors[k][j+2][i].y;\
    z_s[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k][j+1][i].z + sm[2] * coors[k][j+2][i].z;\
    x_t[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k-1][j][i].x + sp[-2] * coors[k-2][j][i].x;\
    y_t[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k-1][j][i].y + sp[-2] * coors[k-2][j][i].y;\
    z_t[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k-1][j][i].z + sp[-2] * coors[k-2][j][i].z;\
    g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
    g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
    g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
    g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
    g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
    g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
    g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
    jac[k][j][i] = sqrt(g);\
    vec__x_r[k][j][i] = x_r[c];\
    vec__y_r[k][j][i] = y_r[c];\
    vec__z_r[k][j][i] = z_r[c];\
    vec__x_s[k][j][i] = x_s[c];\
    vec__y_s[k][j][i] = y_s[c];\
    vec__z_s[k][j][i] = z_s[c];\
    vec__x_t[k][j][i] = x_t[c];\
    vec__y_t[k][j][i] = y_t[c];\
    vec__z_t[k][j][i] = z_t[c];\
    g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
    g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
    g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
    coef[MET_SW][k][j+1][i] += .25 * sqrt(g) * g12;\
    coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
    coef[MET_SE][k][j+1][i] -= .25 * sqrt(g) * g12;\
    coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
    coef[MET_SB][k][j+1][i] += .25 * sqrt(g) * g23;\
}\
\
\
if ((xs + xm == ngx && !per_x) && (ys + ym == ngy && !per_y) && (zs == 0 && !per_z)) {\
    i = ngx - 1;\
    j = ngy - 1;\
    k = 0;\
    x_r[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j][i-1].x + sp[-2] * coors[k][j][i-2].x;\
    y_r[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j][i-1].y + sp[-2] * coors[k][j][i-2].y;\
    z_r[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j][i-1].z + sp[-2] * coors[k][j][i-2].z;\
    x_s[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j-1][i].x + sp[-2] * coors[k][j-2][i].x;\
    y_s[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j-1][i].y + sp[-2] * coors[k][j-2][i].y;\
    z_s[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j-1][i].z + sp[-2] * coors[k][j-2][i].z;\
    x_t[c] = sm[0] * coors[k][j][i].x + sm[1] * coors[k+1][j][i].x + sm[2] * coors[k+2][j][i].x;\
    y_t[c] = sm[0] * coors[k][j][i].y + sm[1] * coors[k+1][j][i].y + sm[2] * coors[k+2][j][i].y;\
    z_t[c] = sm[0] * coors[k][j][i].z + sm[1] * coors[k+1][j][i].z + sm[2] * coors[k+2][j][i].z;\
    g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
    g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
    g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
    g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
    g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
    g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
    g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
    jac[k][j][i] = sqrt(g);\
    vec__x_r[k][j][i] = x_r[c];\
    vec__y_r[k][j][i] = y_r[c];\
    vec__z_r[k][j][i] = z_r[c];\
    vec__x_s[k][j][i] = x_s[c];\
    vec__y_s[k][j][i] = y_s[c];\
    vec__z_s[k][j][i] = z_s[c];\
    vec__x_t[k][j][i] = x_t[c];\
    vec__y_t[k][j][i] = y_t[c];\
    vec__z_t[k][j][i] = z_t[c];\
    g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
    g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
    g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
    coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
    coef[MET_WB][k+1][j][i] += .25 * sqrt(g) * g13;\
    coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
    coef[MET_EB][k+1][j][i] -= .25 * sqrt(g) * g13;\
    coef[MET_SB][k+1][j][i] += .25 * sqrt(g) * g23;\
    coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
    coef[MET_NB][k+1][j][i] -= .25 * sqrt(g) * g23;\
}\
\
\
if ((xs + xm == ngx && !per_x) && (ys + ym == ngy && !per_y) && (zs + zm == ngz && !per_z)) {\
    i = ngx - 1;\
    j = ngy - 1;\
    k = ngz - 1;\
    x_r[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j][i-1].x + sp[-2] * coors[k][j][i-2].x;\
    y_r[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j][i-1].y + sp[-2] * coors[k][j][i-2].y;\
    z_r[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j][i-1].z + sp[-2] * coors[k][j][i-2].z;\
    x_s[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k][j-1][i].x + sp[-2] * coors[k][j-2][i].x;\
    y_s[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k][j-1][i].y + sp[-2] * coors[k][j-2][i].y;\
    z_s[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k][j-1][i].z + sp[-2] * coors[k][j-2][i].z;\
    x_t[c] = sp[0] * coors[k][j][i].x + sp[-1] * coors[k-1][j][i].x + sp[-2] * coors[k-2][j][i].x;\
    y_t[c] = sp[0] * coors[k][j][i].y + sp[-1] * coors[k-1][j][i].y + sp[-2] * coors[k-2][j][i].y;\
    z_t[c] = sp[0] * coors[k][j][i].z + sp[-1] * coors[k-1][j][i].z + sp[-2] * coors[k-2][j][i].z;\
    g_11 = x_r[c]*x_r[c] + y_r[c]*y_r[c] + z_r[c]*z_r[c];\
    g_12 = x_r[c]*x_s[c] + y_r[c]*y_s[c] + z_r[c]*z_s[c];\
    g_13 = x_r[c]*x_t[c] + y_r[c]*y_t[c] + z_r[c]*z_t[c];\
    g_22 = x_s[c]*x_s[c] + y_s[c]*y_s[c] + z_s[c]*z_s[c];\
    g_23 = x_s[c]*x_t[c] + y_s[c]*y_t[c] + z_s[c]*z_t[c];\
    g_33 = x_t[c]*x_t[c] + y_t[c]*y_t[c] + z_t[c]*z_t[c];\
    g = (g_11*(g_22*g_33 - g_23*g_23) - g_12*(g_33*g_12 - g_23*g_13) + g_23*(g_12*g_23 - g_22*g_13));\
    jac[k][j][i] = sqrt(g);\
    vec__x_r[k][j][i] = x_r[c];\
    vec__y_r[k][j][i] = y_r[c];\
    vec__z_r[k][j][i] = z_r[c];\
    vec__x_s[k][j][i] = x_s[c];\
    vec__y_s[k][j][i] = y_s[c];\
    vec__z_s[k][j][i] = z_s[c];\
    vec__x_t[k][j][i] = x_t[c];\
    vec__y_t[k][j][i] = y_t[c];\
    vec__z_t[k][j][i] = z_t[c];\
    g12 = FDIV(1.0,g) * (g_13*g_23 - g_12*g_33);\
    g13 = FDIV(1.0,g) * (g_12*g_23 - g_13*g_22);\
    g23 = FDIV(1.0,g) * (g_13*g_12 - g_11*g_23);\
    coef[MET_SE][k][j][i-1] -= .25 * sqrt(g) * g12;\
    coef[MET_EB][k][j][i-1] -= .25 * sqrt(g) * g13;\
    coef[MET_NB][k][j-1][i] -= .25 * sqrt(g) * g23;\
}\
\
\
}\

#endif
