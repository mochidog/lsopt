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
    double y_min[n];
    lsopt_lm_option_t options = lsopt_lm_option_default;
    options.lm_print_steps = 1;
    
    lsopt_lm(func, x_init, &n, &m, x_min, y_min, &options);
    
    printf("x_min = (%.6f, %.6f)\n", x_min[0], x_min[1]);
    printf("y_min = (%.6f, %.6f)\n", y_min[0], y_min[1]);
    
    return 0;
}

