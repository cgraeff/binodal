//
//  EOS.c
//  NJLv
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#include <math.h>

#include "libdatafun/libdatafun.h"

#include "Constants.h"
#include "Parameters.h"
#include "CommandlineOptions.h"
#include "FermiDiracDistributions.h"
#include "Functions.h"

typedef struct _therm_pot_free_gas_contrib_params{
    double mass;
    double temperature;
    double renormalized_chemical_potential;
} therm_pot_free_gas_contrib_params;

typedef struct _entropy_integrand_parameters{
    double mass;
    double temperature;
    double renormalized_chemical_potential;
} entropy_integrand_parameters;

double QuarkThermodynamicPotentialFreeGasContributionIntegrand(double momentum,
                                                          void * parameters);

// TODO: how is this different from the hadron case?
double QuarkBarionicDensity(double mass,
                            double renormalized_chemical_potential,
                            double temperature)
{
    double constant = NUM_Q_FLAVORS / (pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));

    if (temperature == 0){

        double fermi_momentum_3rd_power =
            pow(QPFermiMomentum(mass, renormalized_chemical_potential), 3.0);

        return constant * fermi_momentum_3rd_power / 3.0;
    }

    double integral =
        FermiDiracDistributionFromDensityIntegral(temperature,
                                                  mass,
                                                  renormalized_chemical_potential);

    return constant * integral;
}

double QuarkThermodynamicPotential(double mass,
                                   double chemical_potential,
                                   double renormalized_chemical_potential,
                                   double temperature)
{
    double first_term =
        QuarkThermodynamicPotentialFreeGasContribution(mass,
                                                  chemical_potential,
                                                  renormalized_chemical_potential,
                                                  temperature);

    double second_term = pow(mass - parameters.quark_model.bare_mass, 2.0)
                         / (4.0 * parameters.quark_model.G_S * CONST_HBAR_C);

    // If G_V == 0, we have to avoid a division by zero
    double third_term = 0.0;
    if (parameters.quark_model.G_V != 0)
        third_term = pow(chemical_potential - renormalized_chemical_potential, 2.0)
        / (4.0 * parameters.quark_model.G_V * CONST_HBAR_C);

    return first_term + second_term - third_term;
}

// TODO: Does the hadron case have a similar function?
double QuarkThermodynamicPotentialFreeGasContribution(double mass,
                                                      double chemical_potential,
                                                      double renormalized_chemical_potential,
                                                      double temperature)
{
    double constant = - NUM_Q_FLAVORS * NUM_Q_COLORS * pow(CONST_HBAR_C, -3.0)
                        / pow(M_PI, 2.0);

    if (temperature == 0){

        double fermi_momentum = QPFermiMomentum(mass, renormalized_chemical_potential);

        double F_diff = F_E(mass, parameters.quark_model.cutoff) - F_E(mass, fermi_momentum);

        return constant
               * (F_diff
                  + renormalized_chemical_potential * pow(fermi_momentum, 3.0) / 3.0);
    }

    therm_pot_free_gas_contrib_params p;
    p.mass = mass;
    p.renormalized_chemical_potential = renormalized_chemical_potential;
    p.temperature = temperature;

    gsl_function F;
    F.function = &QuarkThermodynamicPotentialFreeGasContributionIntegrand;
    F.params = &p;

    double integral = OnedimensionalIntegrator(&F,
                                               parameters.therm_pot_free_gas_integral);

    return constant * integral;
}

// TODO: Does the hadron case have a similar function?
double QuarkThermodynamicPotentialFreeGasContributionIntegrand(double momentum,
                                                               void * parameters)
{
    therm_pot_free_gas_contrib_params * p =
        (therm_pot_free_gas_contrib_params *)parameters;

    double T = p->temperature;
    double mu_r = p->renormalized_chemical_potential;

    double energy = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));

    // From docs:
    // If x is nearly zero, then the common expression log(1 + x) will
    // not be able to produce accurate results, as most (or all) of the
    // information in x will be lost by addition.  Instead, use
    // log1p(x) to perform the same computation without undue
    // loss of accuracy.
    double first_term = T * log1p(exp(-(energy - mu_r)/T));
    double second_term = T * log1p(exp(-(energy + mu_r)/T));

    return pow(momentum, 2.0) * (energy + first_term + second_term);
}

double QuarkPressure(double regularized_thermodynamic_potential, double temperature)
{
    return - regularized_thermodynamic_potential;
}

double QuarkVacuumMassDetermination()
{
    // Prepare function to be passed to the root finding algorithm.
    // No parameters are needed.
    gsl_function F;
    F.function = &QuarkVacuumMassEquation;

    double root;
    int status =
        UnidimensionalRootFinder(&F,
                                 parameters.vacuum_mass_determination,
                                 &root);

    if (status == -1){
        return 0;
    }

    return root;
}

double QuarkVacuumMassEquation(double mass, void * input)
{
    double F_diff = F0(mass, parameters.quark_model.cutoff) - F0(mass, 0.0);
    double term = 2.0 * NUM_Q_COLORS * NUM_Q_FLAVORS * pow(CONST_HBAR_C, -2.0)
                  * parameters.quark_model.G_S * mass * F_diff
                  / pow(M_PI, 2.0);

    return mass - parameters.quark_model.bare_mass - term;
}

double QuarkScalarDensity(double temperature,
                          double mass,
                          double renormalized_chemical_potential)
{
    if (mass == 0){
        return 0;
    }

    double constant = - NUM_Q_COLORS * NUM_Q_FLAVORS
                        / (pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));

    if (temperature == 0){

        double fermi_momentum = QPFermiMomentum (mass, renormalized_chemical_potential);

        return constant * mass * (F0(mass, parameters.quark_model.cutoff)
                                  - F0(mass, fermi_momentum));
    }

    double integral =
        FermiDiracDistributionIntegralFromScalarDensity(temperature,
                                                        mass,
                                                        renormalized_chemical_potential);

    return constant * mass * integral;
}

double QPFermiMomentum(double mass, double renormalized_chemical_potential)
{
    if (pow(renormalized_chemical_potential, 2.0) < pow(mass, 2.0))
        return 0.0;

    return sqrt(pow(renormalized_chemical_potential, 2.0) - pow(mass, 2.0));
}

