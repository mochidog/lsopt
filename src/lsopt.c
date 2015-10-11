#include "lsopt.h"

#include <stdio.h>
#include <math.h>
#include <f2c.h>
#include <clapack.h>

#define SQRT2 1.41421356237309504880168872420969807

const lsopt_lm_option_t lsopt_lm_option_default = {
    100,
    1e-6,
    1e-6,
    1.0,
    2.0,
    1e-6,
    0
};

inline void print_vector(const char *v_name, const double *vec, const int *n) {
    int i;
    printf("%s", v_name);
    for (i = 0; i < (*n); ++i) {
        printf("\t%.8f", vec[i]);
    }
    printf("\n");
}

inline void print_matrix(const double *A, const int *n, const int *m) {
    int i, j;
    printf("\n");
    for (i = 0; i < (*n); ++i) {
        printf("|\t");
        for (j = 0; j < (*m); ++j) {
            printf("%.7f\t", A[j * (*n) + i]);
        }
        printf("|\n");
    }
    printf("\n");
}

inline void max_diag_element(const double *A, const int *n, double *a_max) {
    *a_max = A[0];
    int i;
    for (i = 1; i < (*n); ++i) {
        if (*a_max < A[i * (*n) + i]) 
            *a_max = A[i * (*n) + i];
    }
}

inline void damp_by_mu(double *A, const int *n, const double *mu) {
    int i;
    for (i = 0; i < (*n); ++i) {
        A[i * (*n) + i] += (*mu);
    }
}

inline void clone_vector(const double *src, const int *n, double *des) {
    int i;
    for (i = 0; i < *n; ++i) {
        des[i] = src[i];
    }
}

inline void negate_vector(double *vec, const int *n) {
    int i;
    for (i = 0; i < *n; ++i) {
        vec[i] *= -1.0;
    }
}

inline void add_vectors(const double *v1, const double *v2, const int *n, double *v_ret) {
    int i;
    for (i = 0; i < *n; ++i) {
        v_ret[i] = v1[i] + v2[i];
    }
}

inline void subtract_vectors(const double *v1, const double *v2, const int *n, double *v_ret) {
    int i;
    for (i = 0; i < *n; ++i) {
        v_ret[i] = v1[i] - v2[i];
    }
}

inline void multiply_vectors(const double *v1, const double *v2, const int *n, double *v_ret) {
    int i;
    *v_ret = 0;
    for (i = 0; i < *n; ++i) {
        *v_ret += v1[i] * v2[i];
    }
}

inline void max_element(const double *vec, const int *n, double *e_max) {
    *e_max = vec[0];
    int i;
    for (i = 1; i < (*n); ++i) {
        if (*e_max < vec[i])
            *e_max = vec[i];
    }
}

inline void range_vector(integer * const vec, const integer *low, const integer *high) {
    int i;
    for (i = 0; i < (*high) - (*low); ++i) {
        vec[i] = (*low) + i;
    }
}

void compute_jacobian(lsopt_func_t func, const double *x, const double *delta_x,
        const int *x_dim, const int *y_dim, double *J) {
    double y[*y_dim];
    double *y_new = J;
    double x_new[*x_dim];
    double inv_delta_x = 1.0 / (*delta_x);
    clone_vector(x, x_dim, x_new);
    func(x, x_dim, y, y_dim);
    int i, j;
    for (j = 0; j < (*x_dim); j++) {
        y_new += j * (*y_dim);
        x_new[j] = x[j] + (*delta_x);
        func(x_new, x_dim, y_new, y_dim);
        x_new[j] = x[j];
        for (i = 0; i < (*y_dim); i++) {
            y_new[i] = (y_new[i] - y[i]) * inv_delta_x;
        }
    }
}

void lsopt_lm(lsopt_func_t func, const double *x_init, const int *x_dim, const int *y_dim,
        double *x_min, double *F_min, const lsopt_lm_option_t* options) {
    int k = 0, k_max, print_steps;
    double e1, e2, tau, nu, dx;
    if (!options) {
        options = &lsopt_lm_option_default;
    }
    k_max = options->lm_max_iteration;
    e1 = options->lm_epsilon_1;
    e2 = options->lm_epsilon_2;
    tau = options->lm_tau;
    nu = options->lm_nu_0;
    dx = options->lm_jacobian_delta;
    print_steps = options->lm_print_steps;
    integer n = *x_dim;
    integer m = *y_dim;
    double A[n * n];
    double J[m * n];
    double y[m];
    double g[n];
    double h[n];
    double x_new[n];
    double y_new[m];
    integer pivot[n];
    integer l = 1;
    int one = 1;
    integer n1 = n + 1;
    range_vector(pivot, &l, &n1); 
    // pivot is used in dgetrs_ lapack routine
    
    clone_vector(x_init, x_dim, x_min);
    double *x = x_min;
    // y = f(x)
    func(x, x_dim, y, y_dim);
    if (print_steps) {
        printf("step %d\n", 0);
        print_vector("x", x, x_dim);
        print_vector("y", y, y_dim);
        double F;
        multiply_vectors(y, y, y_dim, &F);
        F /= 2.0;
        print_vector("F", &F, &one);
    }
    
    compute_jacobian(func, x, &dx, x_dim, y_dim, J);
    char transA = 't';
    char transB = 'n';
    double alpha = 1.0;
    double beta = 0.0;
    // A = J' * J
    dgemm_(&transA, &transB, &n, &n, &m, &alpha, J, &n, J, &m, &beta, A, &n);
    // g = J' * y
    dgemm_(&transA, &transB, &n, &l, &m, &alpha, J, &n, y, &m, &beta, g, &n);
    
    double g_inf_norm;
    max_element(g, x_dim, &g_inf_norm);
    int found = (g_inf_norm < e1);
    double max_A_diag_element;
    max_diag_element(A, x_dim, &max_A_diag_element);
    double mu = tau * max_A_diag_element;
    
    int iterator = 0;
    while (!found && k < k_max) {
        k++;
        // A = A + mu * I
        damp_by_mu(A, x_dim, &mu);
        integer info;
        clone_vector(g, x_dim, h);
        negate_vector(h, x_dim);
        // solve linear system (A + mu * I) * h = -g for h
        transA = 'n';
        dgetrs_(&transA, &n, &l, A, &n, pivot, h, &n, &info); 
        // h holds the result
        // TODO: change to double precision symmetric positive definitive solver
        
        double x_norm = dnrm2_(&n, x, &l);
        double h_norm = dnrm2_(&n, h, &l);
        if (h_norm <= e2 * (x_norm + e2)) {
            found = 1;
        }
        else {
            // L(0) - L(h) = (mu * h' * h - h' * g) / 2
            double hg = 0.0;
            // g' * h = f'Jh
            multiply_vectors(h, g, x_dim, &hg);
            double hh;
            multiply_vectors(h, h, x_dim, &hh);
            double dL = (mu * hh - hg) / 2;
            // x_new = x + h
            add_vectors(x, h, x_dim, x_new);
            func(x_new, x_dim, y_new, y_dim);
            // dF = (||y||^2 - ||y_new||^2) / 2
            double dF = 0.0, yy = 0.0;
            multiply_vectors(y, y, y_dim, &dF);
            multiply_vectors(y_new, y_new, y_dim, &yy);
            dF -= yy;
            dF /= 2.0;
            
            double rho = dF / dL;
            if (rho > 0) {      // step acceptable
                iterator++;
                x = x_new;
                // y = f(x)
                func(x, x_dim, y, y_dim);
                if (print_steps) {
                    printf("step %d\n", iterator);
                    print_vector("x", x, x_dim);
                    print_vector("y", y, y_dim);
                    double F;
                    multiply_vectors(y, y, y_dim, &F);
                    F /= 2.0;
                    print_vector("F", &F, &one);
                }
                compute_jacobian(func, x, &dx, x_dim, y_dim, J);
                // A = J' * J
                transA = 't';
                dgemm_(&transA, &transB, &n, &n, &m, &alpha, J, &n, J, &m, &beta, A, &n);
                // g = J' * y
                dgemm_(&transA, &transB, &n, &l, &m, &alpha, J, &n, y, &m, &beta, g, &n);
                max_element(g, x_dim, &g_inf_norm);
                found = (g_inf_norm < e1);
                mu = mu * dmax(1.0 / 3.0, 1.0 - pow(2 * rho - 1, 3));
                nu = 2;
            }
            else {
                mu = mu * nu;
                nu = 2 * nu;
            }
        }
    }
    clone_vector(x, x_dim, x_min);
    multiply_vectors(y, y, y_dim, F_min);
    *F_min /= 2.0;
}

