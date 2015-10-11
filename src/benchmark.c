#include <stdio.h>
#include "lsopt.h"

void func(const double *x, const int *n, double *y, const int *m) {
    y[0] = x[0] - 1;
    y[1] = x[1] - 4;
}

int main() {
    printf("this is a sample\n");
    const int n = 2, m = 2;
    double x_init[n];
    x_init[0] = 5.0;
    x_init[1] = 6.0;
    double x_min[n];
    double F_min;
    lsopt_lm_option_t options = lsopt_lm_option_default;
    options.lm_print_steps = 1;
    
    lsopt_lm(func, x_init, &n, &m, x_min, &F_min, &options);
    
    printf("x_min = (%.6f, %.6f)\n", x_min[0], x_min[1]);
    printf("F_min = %.6f\n", F_min);
    
    return 0;
}

