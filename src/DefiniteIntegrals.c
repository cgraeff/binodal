//
//  DefiniteIntegrals.c
//  binodal
//
//  Created by Clebson Graeff on 2017-06-20.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#include <math.h>

#include "Constants.h"
#include <stdio.h>

// F_0 calculates the integral
//
// \int \frac{k^2}{\sqrt{k^2 + m^2}} dk &=
//     \frac{1}{2}\left[k\epsilon
//                      - m^2\ln\left(\frac{\epsilon + k}{m}\right)
//               \right]
//
// in the interval [0, p]. ($\epsilon \equiv E$)
//
// Inputs: mass and momentum in MeV
// Output: the definite integral, in (MeV)^2
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

// F_2 calculates the integral
//
// \int \frac{k^4}{\sqrt{k^2 + m^2}} dk =
//  \frac{1}{4}\left[k^3\epsilon
//                   - \frac{3}{2} m^2k\epsilon
//                   + \frac{3}{2}m^4\ln\left(\frac{\epsilon + k}{m}\right)
//            \right]
//
// in the interval [0, p]. ($\epsilon \equiv E$)
// Inputs: mass and momentum in MeV
// Output: the definite integral, in (MeV)^4
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

// F_E calculates the integral
//
// \int k^2 \sqrt{k^2 + m^2} \;dk =
//     \frac{1}{4}\left[k \epsilon^3
//                      - \frac{1}{2} m^2 k \epsilon
//                      - \frac{1}{2} m^4\ln\left(\frac{k+\epsilon}{m}\right)
//               \right]
//
// in the interval [0, p]. ($\epsilon \equiv E$)
// Inputs: mass and momentum in MeV
// Output: the definite integral, in (MeV)^4
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

