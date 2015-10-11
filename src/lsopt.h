/*
 * Non-linear least square optimizer based on
 * K. Madsen, H.B. Nielsen, O.Tingleff, "Methods for Non-Linear Least Squares Problems", 2014
 *
 * It includes the following optimizers:
 * Gradient-Descent
 * Gauss-Newton
 * Levenberg-Marquardt
 * Secant Levenberg-Marquardt
 * Powell Dog-Leg
 * Secant Powell Dog-Leg
 *
 *
 * Abbreviations:
 * Least Sqaure OPTimization    lsopt
 * Levenberg-Marquardt          lm
 * x vector dimension           n
 * y vector dimension           m
 *
 */


typedef struct lsopt_lm_option_s {
    int     lm_max_iteration;
    double  lm_epsilon_1;
    double  lm_epsilon_2;
    double  lm_tau;
    double  lm_nu_0;
    double  lm_jacobian_delta;
    int     lm_print_steps;
} lsopt_lm_option_t;

extern const lsopt_lm_option_t lsopt_lm_option_default;

typedef void (*lsopt_func_t)(const double *x, const int *x_dim, double *y, const int *y_dim);

void lsopt_lm(lsopt_func_t func, const double *x_init, const int *x_dim, const int *y_dim,
        double *x_min, double *F_min, const lsopt_lm_option_t* options);

