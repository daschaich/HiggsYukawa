// -----------------------------------------------------------------
// Defines and subroutine declarations for SO(4)
// with fermions in the DIMF-dimensional adjoint rep
// The original names now refer to objects of dimension DIMF or DIMFxDIMF
// New objects with suffix _f have dimension NCOL or NCOLxNCOL
#ifndef _SUN_H
#define _SUN_H

#include "../include/complex.h"
#include "../include/random.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
//typedef struct { float e[6]; } fso4_antisym;
typedef struct { float e[3]; } fso4_selfdual;
typedef struct { float c[4]; } fso4_vector;

//typedef struct { double e[6]; } dso4_antisym;
typedef struct { double e[3]; } dso4_selfdual;
typedef struct { double c[4]; } dso4_vector;

#if (PRECISION == 1)
#define so4_selfdual    fso4_selfdual
#define so4_vector      fso4_vector
#else
#define so4_selfdual    dso4_selfdual
#define so4_vector      dso4_vector
#endif

#define PLUS 1          // Flags for selecting M or M_adjoint
#define MINUS -1
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Subroutine definitions
// Vector operations
// In file clearvec_f.c
void clearvec_f(so4_vector_f *v);

void dumpvec(so4_vector *v);

// In file clearvec.c
void clearvec(so4_vector *v);

// In file msq_su3vec.c
Real magsq_su3vec(so4_vector *v);

// In file su3vec_copy.c
void su3vec_copy(so4_vector *a, so4_vector *b);

// In file so4_rdot.c
Real so4_rdot(so4_vector *a, so4_vector *b);

// In file so4_dot.c
complex so4_dot(so4_vector *a, so4_vector *b);

// In file addvec.c
void add_so4_vector(so4_vector *a, so4_vector *b, so4_vector *c);

// In file subvec.c
void sub_so4_vector(so4_vector *a, so4_vector *b, so4_vector *c);

// In file s_m_vec.c
void scalar_mult_so4_vector(so4_vector *src, Real scalar,
                            so4_vector *dest);

// In file s_m_sum_vec.c
void scalar_mult_sum_so4_vector(so4_vector *src1, so4_vector *src2,
                                Real scalar);

// In file s_m_a_vec.c
void scalar_mult_add_so4_vector(so4_vector *src1, so4_vector *src2,
                                Real scalar, so4_vector *dest);

// In file s_m_s_vec.c
void scalar_mult_sub_so4_vector(so4_vector *src1, so4_vector *src2,
                                Real scalar, so4_vector *dest);

// In file cs_m_vec.c
void c_scalar_mult_su3vec(so4_vector *src, complex *phase, so4_vector *dest);

// In file cs_m_a_vec.c
void c_scalar_mult_add_su3vec(so4_vector *v1, complex *phase, so4_vector *v2);

// In file cs_m_s_vec.c
void c_scalar_mult_sub_su3vec(so4_vector *v1, complex *phase, so4_vector *v2);

// In file sub4vecs.c
void sub_four_so4_vecs(so4_vector *a, so4_vector *b1, so4_vector *b2,
                       so4_vector *b3, so4_vector *b4);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Matrix operations
// In file clear_mat_f
void clear_su3mat_f(so4_matrix_f *dest);

// In file trace_so4_f.c
complex trace_so4_f(so4_matrix_f *a);

// In file realtr_f.c
Real realtrace_so4_f(so4_matrix_f *a, so4_matrix_f *b);

// In file complextr_f.c
complex complextrace_so4_f(so4_matrix_f *a, so4_matrix_f *b);

// b <-- a, in file su3mat_copy_f.c
void su3mat_copy_f(so4_matrix_f *a, so4_matrix_f *b);

// b <-- adag, in file so4_adjoint_f.c
void so4_adjoint_f(so4_matrix_f *a, so4_matrix_f *b);

// In file addmat_f.c
void add_so4_matrix_f(so4_matrix_f *a, so4_matrix_f *b, so4_matrix_f *c);

// In file submat_f.c
void sub_so4_matrix_f(so4_matrix_f *a, so4_matrix_f *b, so4_matrix_f *c);

// In file s_a_d_mat_f.c
void scalar_add_diag_so4_f(so4_matrix_f *a, Real s);

// In file s_m_mat_f.c
void scalar_mult_so4_matrix_f(so4_matrix_f *src, Real scalar,
                              so4_matrix_f *dest);

// In file s_m_a_mat_f.c
void scalar_mult_add_so4_matrix_f(so4_matrix_f *src1, so4_matrix_f *src2,
                                  Real scalar, so4_matrix_f *dest);

// In file s_m_s_mat_f.c
void scalar_mult_sub_so4_matrix_f(so4_matrix_f *src1, so4_matrix_f *src2,
                                  Real scalar, so4_matrix_f *dest);

// In file cs_a_d_mat_f.c
void c_scalar_add_diag_so4_f(so4_matrix_f *a, complex *f);

// In file cs_m_mat_f.c
void c_scalar_mult_su3mat_f(so4_matrix_f *b, complex *s, so4_matrix_f *c);

// In file cs_m_a_mat_f.c
void c_scalar_mult_add_su3mat_f(so4_matrix_f *m1, so4_matrix_f *m2,
                                complex *phase, so4_matrix_f *m3);

// In file dumpmat_f.c
void dumpmat_f(so4_matrix_f *m);

// In file m_mat_nn_f.c
void mult_so4_nn_f(so4_matrix_f *a, so4_matrix_f *b, so4_matrix_f *c);

// In file m_mat_na_f.c
void mult_so4_na_f(so4_matrix_f *a, so4_matrix_f *b, so4_matrix_f *c);

// In file m_mat_an_f.c
void mult_so4_an_f(so4_matrix_f *a, so4_matrix_f *b, so4_matrix_f *c);

// Relate fundamental matrices and anti-hermitian matrices
// In file make_ahmat.c
void make_anti_hermitian(so4_matrix_f *m, anti_hermitmat *ah);

// In file rand_ahmat.c
void random_anti_hermitian(anti_hermitmat *ah, double_prn *prn_pt);

// In file uncmp_ahmat.c
void uncompress_anti_hermitian(anti_hermitmat *ah, so4_matrix_f *m);

// In file cmp_ahmat.c
void compress_anti_hermitian(so4_matrix_f *m, anti_hermitmat *ah);
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Fermion rep matrix operations
// In file trace_su3.c
complex trace_su3(so4_matrix *a);

// In file realtr.c
Real realtrace_su3(so4_matrix *a, so4_matrix *b);

// In file complextr.c
complex complextrace_su3(so4_matrix *a, so4_matrix *b);

// In file addmat.c
void add_so4_matrix(so4_matrix *a, so4_matrix *b, so4_matrix *c);

// In file submat.c
void sub_so4_matrix(so4_matrix *a, so4_matrix *b, so4_matrix *c);

// In file s_m_mat.c
void scalar_mult_so4_matrix(so4_matrix *src, Real scalar, so4_matrix *dest);

// In file s_m_a_mat.c
void scalar_mult_add_so4_matrix(so4_matrix *src1, so4_matrix *src2,
                                Real scalar, so4_matrix *dest);

// In file s_m_s_mat.c
void scalar_mult_sub_so4_matrix(so4_matrix *src1, so4_matrix *src2,
                                Real scalar, so4_matrix *dest);

// In file cs_m_mat.c
void c_scalar_mult_su3mat(so4_matrix *src, complex *scalar,
                          so4_matrix *dest);

// In file cs_m_a_mat.c
void c_scalar_mult_add_su3mat(so4_matrix *src1, so4_matrix *src2,
                              complex *scalar, so4_matrix *dest);

// In file cs_m_s_mat.c
void c_scalar_mult_sub_su3mat(so4_matrix *src1, so4_matrix *src2,
                              complex *scalar, so4_matrix *dest);

// In file so4_adjoint.c
void so4_adjoint(so4_matrix *a, so4_matrix *b);

// In file clear_mat.c
void clear_su3mat(so4_matrix *dest);

// In file su3mat_copy.c
void su3mat_copy(so4_matrix *a, so4_matrix *b);

// In file dumpmat.c
void dumpmat(so4_matrix *m);

// In file m_mat_nn.c
void mult_so4_nn(so4_matrix *a, so4_matrix *b, so4_matrix *c);

// In file m_mat_na.c
void mult_so4_na(so4_matrix *a, so4_matrix *b, so4_matrix *c);

// In file m_mat_an.c
void mult_so4_an(so4_matrix *a, so4_matrix *b, so4_matrix *c);

// c <-- a * b, in file m_matvec.c
void mult_so4_mat_vec(so4_matrix *a, so4_vector *b, so4_vector *c);

// c <-- c + a * b, in file m_matvec_s.c
void mult_so4_mat_vec_sum(so4_matrix *a, so4_vector *b, so4_vector *c);

// c <-- c - a * b, in file m_matvec_ns.c
void mult_so4_mat_vec_nsum(so4_matrix *a, so4_vector *b, so4_vector *c);

// c <-- adag * b, in file m_amatvec.c
void mult_adj_so4_mat_vec(so4_matrix *a, so4_vector *b, so4_vector *c);

// c <-- c + adag * b, in file m_amatvec_s.c
void mult_adj_so4_mat_vec_sum(so4_matrix *a, so4_vector *b, so4_vector *c);

// c <-- c - adag * b, in file m_amatvec_ns.c
void mult_adj_so4_mat_vec_nsum(so4_matrix *a, so4_vector *b, so4_vector *c);

// c <-- b * a, in file m_vecmat.c
void mult_so4_vec_mat(so4_vector *b, so4_matrix *a, so4_vector *c);

// c <-- b * adag, in file m_vecamat.c
void mult_so4_vec_adj_mat(so4_vector *b, so4_matrix *a, so4_vector *c);

// c <-- c + b * adag, in file m_vecamat_s.c
void mult_so4_vec_adj_mat_sum(so4_vector *b, so4_matrix *a, so4_vector *c);

// In file m_amv_4dir.c
void mult_adj_so4_mat_vec_4dir(so4_matrix *a, so4_vector *b, so4_vector *c);

// In file m_amv_4vec.c
void mult_adj_so4_mat_4vec(so4_matrix *mat, so4_vector *src,
                           so4_vector *dest0, so4_vector *dest1,
                           so4_vector *dest2, so4_vector *dest3);

// In file m_mv_s_4dir.c
void mult_so4_mat_vec_sum_4dir(so4_matrix *a, so4_vector *b0,
                               so4_vector *b1, so4_vector *b2,
                               so4_vector *b3, so4_vector *c);

// In file so4_proj.c
void so4_projector(so4_vector *a, so4_vector *b, so4_matrix *c);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Miscellaneous routines
// In file gaussrand.c
Real gaussian_rand_no(double_prn *prn_pt);

#include "../include/int32type.h"
void byterevn(int32type w[], int n);
void byterevn64(int32type w[], int n);

#endif
// -----------------------------------------------------------------
