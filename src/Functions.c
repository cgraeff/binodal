
#include <math.h>

#include "Constants.h"

double F0(double mass,
          double momentum)
{
	double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));

    double first_term = momentum * E;
    double second_term = - pow(mass, 2.0) * log(momentum + E);
    double third_term = pow(mass, 2.0) * log(mass);

    // In the limit of mass -> 0, the logarithm may give a NaN,
    // so when that happens, just make sure the third term
    // have the correct value
    if (third_term != third_term)
        third_term = 0;

    return (1.0 / 2.0) * (first_term +second_term + third_term);
}

double F2(double mass,
          double momentum)
{
    if (mass == 0.0){
        return pow(momentum, 4.0) / 4.0;
    }

	double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));

	return (1.0 / 8.0) * (-3.0 * pow(mass, 2.0) * momentum
                          + 2.0 * pow(momentum, 3.0)) * E
           + (3.0 / 8.0) * pow(mass, 4.0) * log((momentum + E) / mass);
}

double F_E(double mass, double momentum)
{
    if (mass == 0)
        return pow(momentum, 4.0) / 4.0;

    double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));

    return (momentum * pow(E, 3.0)
            - 0.5 * pow(mass, 2.0) * momentum * E
            - 0.5 * pow(mass, 4.0) * log ((momentum + E) / mass))
           / 4.0;
}

