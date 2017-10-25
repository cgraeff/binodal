//
//  Functions.c
//  binodal
//
//  Created by Clebson Graeff on 2017-06-20.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#include <math.h>

#include "Constants.h"
#include <stdio.h>

double F0(double mass,
          double momentum)
{
    if (mass < CONST_ZERO_MASS_TOL)
        return pow(momentum, 2.0) / 2.0;

    double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));

    double first_term = momentum * E;
    double second_term = - pow(mass, 2.0) * log(momentum + E);
    double third_term = pow(mass, 2.0) * log(mass); // This term won't ever
                                                    // make a difference in a
                                                    // definite integral

    // In the limit of mass -> 0, the logarithm may give a NaN,
    // so when that happens, just make sure the third term
    // have the correct value
    if (third_term != third_term)
        third_term = 0;

    return (1.0 / 2.0) * (first_term + second_term + third_term);
}

double F2(double mass,
          double momentum)
{
    if (mass < CONST_ZERO_MASS_TOL){
        return pow(momentum, 4.0) / 4.0;
    }

	double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));

	return (1.0 / 8.0) * (-3.0 * pow(mass, 2.0) * momentum
                          + 2.0 * pow(momentum, 3.0)) * E
           + (3.0 / 8.0) * pow(mass, 4.0) * log((momentum + E) / mass);
}

double F_E(double mass, double momentum)
{
    if (mass < CONST_ZERO_MASS_TOL)
        return pow(momentum, 4.0) / 4.0;

    double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));

    double first_term = momentum * pow(E, 3.0);
    double second_term = - 0.5 * pow(mass, 2.0) * momentum * E;
    double third_term = - 0.5 * pow(mass, 4.0) * log ((momentum + E) / mass);

    return (first_term + second_term + third_term) / 4.0;
}

