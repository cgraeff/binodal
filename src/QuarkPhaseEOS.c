//
//  EOS.c
//  NJLv
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#include <math.h>

#include <gsl/gsl_multimin.h>

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
    double constant = NUM_Q_COLORS / (pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));

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

double QuarkThermodynamicPotential(double up_mass,
                                   double down_mass,
                                   double up_chemical_potential,
                                   double down_chemical_potential,
                                   double up_renormalized_chemical_potential,
                                   double down_renormalized_chemical_potential,
                                   double temperature)
{
    double up_free_term =
        QuarkThermodynamicPotentialFreeGasContribution(up_mass,
                                                       up_chemical_potential,
                                                       up_renormalized_chemical_potential,
                                                       temperature);
    double down_free_term =
        QuarkThermodynamicPotentialFreeGasContribution(down_mass,
                                                       down_chemical_potential,
                                                       down_renormalized_chemical_potential,
                                                       temperature);

    double up_scalar_density =
        QuarkScalarDensity (temperature,
                            up_mass,
                            up_renormalized_chemical_potential);

    double down_scalar_density =
        QuarkScalarDensity (temperature,
                            down_mass,
                            down_renormalized_chemical_potential);

    double scalar_term = parameters.quark_model.G_S * CONST_HBAR_C
                         * pow(up_scalar_density + down_scalar_density, 2.0);

    double up_barionic_density =
        QuarkBarionicDensity (up_mass,
                              up_renormalized_chemical_potential,
                              temperature);
    double down_barionic_density =
        QuarkBarionicDensity (down_mass,
                              down_renormalized_chemical_potential,
                              temperature);

    double vector_term = - parameters.quark_model.G_V * CONST_HBAR_C
                           * pow(up_barionic_density + down_barionic_density, 2.0);

    return (up_free_term + down_free_term) + scalar_term + 0.0 * vector_term;
}

// TODO: Does the hadron case have a similar function?
double QuarkThermodynamicPotentialFreeGasContribution(double mass,
                                                      double chemical_potential,
                                                      double renormalized_chemical_potential,
                                                      double temperature)
{
    double constant = - NUM_Q_COLORS * pow(CONST_HBAR_C, -3.0)
                        / pow(M_PI, 2.0);

    if (temperature == 0){

        double fermi_momentum = QPFermiMomentum(mass, renormalized_chemical_potential);

        double upper_limit = F_E(mass, parameters.quark_model.cutoff);
        double lower_limit =  F_E(mass, fermi_momentum);
        double F_diff =  upper_limit - lower_limit;

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

typedef struct _minimizer_function_parameters{
    double up_chemical_potential;
    double down_chemical_potential;
    double up_renormalized_chemical_potential;
    double down_renormalized_chemical_potential;
    double temperature;
} minimizer_function_parameters;

double MinimizationHelperFunction(const gsl_vector * x, void * params)
{
    minimizer_function_parameters * p = (minimizer_function_parameters *)params;

    double up_mass = gsl_vector_get(x, 0);
    double down_mass = gsl_vector_get(x,1);

    return -QuarkThermodynamicPotential(up_mass,
                                        down_mass,
                                        p->up_chemical_potential,
                                        p->down_chemical_potential,
                                        p->up_renormalized_chemical_potential,
                                        p->down_renormalized_chemical_potential,
                                        p->temperature);
}

void QuarkThermodynamicPotentialMinimum(double * up_mass_at_minimum,
                                        double * down_mass_at_minimum)
{
    int max_iter = 20000;
    double tolerance = 1.0E-3;
    double up_vacuum_mass_guess = 300.0;
    double down_vacuum_mass_guess = 300.0;
    double up_mass_step = 0.5;
    double down_mass_step = 0.5;

    minimizer_function_parameters p;
    p.up_chemical_potential = 0.0;
    p.down_chemical_potential = 0.0;
    p.up_renormalized_chemical_potential = 0.0;
    p.down_renormalized_chemical_potential = 0.0;
    p.temperature = 0.0;

    const gsl_multimin_fminimizer_type * minimizer_type = gsl_multimin_fminimizer_nmsimplex2;
    const int minimizer_dimension = 2;

    gsl_multimin_fminimizer * minimizer =
        gsl_multimin_fminimizer_alloc(minimizer_type, minimizer_dimension);

    gsl_multimin_function f;
    f.f = &MinimizationHelperFunction;
    f.params = &p;
    f.n = minimizer_dimension;

    gsl_vector * guesses = gsl_vector_alloc(minimizer_dimension);
    gsl_vector * step_size = gsl_vector_alloc(minimizer_dimension);

    gsl_vector_set(guesses, 0, up_vacuum_mass_guess);
    gsl_vector_set(guesses, 1, down_vacuum_mass_guess);

    gsl_vector_set(step_size, 0, up_mass_step);
    gsl_vector_set(step_size, 1, down_mass_step);

    gsl_multimin_fminimizer_set(minimizer,
                                &f,
                                guesses,
                                step_size);

    int status;
    int iter = 0;
    do {

        if (iter > max_iter){
            printf("%s:%d: Max iterations reached.\n",
                   __FILE__,
                   __LINE__);

            abort();
        }

        status = gsl_multimin_fminimizer_iterate (minimizer);

        if (status == GSL_ENOPROG){
            printf("%s:%d: Can't improve the result beyond current one (no progress).",
                   __FILE__,
                   __LINE__);
            abort();
        }

        double size = gsl_multimin_fminimizer_size(minimizer);

        status = gsl_multimin_test_size(size, tolerance);

        iter++;

    } while (status == GSL_CONTINUE);

    gsl_vector_free(guesses);
    gsl_vector_free(step_size);

    *up_mass_at_minimum = gsl_vector_get(minimizer->x, 0);
    *down_mass_at_minimum = gsl_vector_get(minimizer->x, 1);

    gsl_multimin_fminimizer_free (minimizer);

    return;
}

void QuarkVacuumMassDetermination(double * up_vacuum_mass,
                                  double * down_vacuum_mass)
{
    // Set up parameters to be passed to helper function


    // Set dimension (number of equations|variables to solve|find)
    const int dimension = 2;

    gsl_multiroot_function f;
    f.f = &QuarkVacuumMassDeterminationEquation;
    f.n = dimension;
    f.params = NULL;

    gsl_vector * initial_guess = gsl_vector_alloc(dimension);
    gsl_vector * return_results = gsl_vector_alloc(dimension);

    gsl_vector_set(initial_guess,
                   0,
                   sqrt(parameters.vacuum_mass_determination.up_vacuum_mass_guess));
    gsl_vector_set(initial_guess,
                   1,
                   sqrt(parameters.vacuum_mass_determination.down_vacuum_mass_guess));

    int status =
        MultidimensionalRootFinder(dimension,
                                   &f,
                                   initial_guess,
                                   parameters.vacuum_mass_determination.abs_error,
                                   parameters.vacuum_mass_determination.rel_error,
                                   parameters.vacuum_mass_determination.max_iter,
                                   return_results);

    if (status != 0){
        printf("%s:%d: Something is wrong with the rootfinding.\n",
               __FILE__,
               __LINE__);
        abort();
    }

    // Save results in return variables,
    // taking care of the mappinps
    *up_vacuum_mass = pow(gsl_vector_get(return_results, 0), 2.0);
    *down_vacuum_mass = pow(gsl_vector_get(return_results, 1), 2.0);

    // Free vectors
    gsl_vector_free(initial_guess);
    gsl_vector_free(return_results);

    return;
}

int QuarkVacuumMassDeterminationEquation(const gsl_vector   *x,
                                         void *params,
                                         gsl_vector *return_values)
{

    const double up_renorm_chem_pot = 0.0;
    const double down_renorm_chem_pot = 0.0;

    const double up_vacuum_mass = pow(gsl_vector_get(x, 0), 2.0);
    const double down_vacuum_mass = pow(gsl_vector_get(x, 1), 2.0);

    double up_scalar_density =
        QuarkScalarDensity(parameters.variables.temperature,
                           up_vacuum_mass,
                           up_renorm_chem_pot);

    double down_scalar_density =
        QuarkScalarDensity(parameters.variables.temperature,
                           down_vacuum_mass,
                           down_renorm_chem_pot);

    double up_zeroed_gap_eq =
        QuarkZeroedGapEquation(up_vacuum_mass,
                               up_scalar_density,
                               down_scalar_density);

    double down_zeroed_gap_eq =
        QuarkZeroedGapEquation(down_vacuum_mass,
                               up_scalar_density,
                               down_scalar_density);

    gsl_vector_set(return_values, 0, up_zeroed_gap_eq);
    gsl_vector_set(return_values, 1, down_zeroed_gap_eq);

    return GSL_SUCCESS;
}

double QuarkVacuumMassEquation(double mass, void * input)
{
    double F_diff = F0(mass, parameters.quark_model.cutoff) - F0(mass, 0.0);

    double F_P_diff = F_P_0(mass, parameters.quark_model.cutoff) - F_P_0(mass, 0.0);

    return (2.0 * mass + parameters.quark_model.G_S) * F_diff
           + parameters.quark_model.G_S * mass * F_P_diff;
}

double QuarkScalarDensity(double temperature,
                          double mass,
                          double renormalized_chemical_potential)
{
    if (mass == 0){
        return 0;
    }

    double constant = - NUM_Q_COLORS
                        / (pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));

    if (temperature == 0){

        double fermi_momentum = QPFermiMomentum (mass, renormalized_chemical_potential);

        double upper_limit = F0(mass, parameters.quark_model.cutoff);
        double lower_limit = F0(mass, fermi_momentum);

        return constant * mass * (upper_limit - lower_limit);
    }

    double integral =
        FermiDiracDistributionIntegralFromScalarDensity(temperature,
                                                        mass,
                                                        renormalized_chemical_potential);

    return constant * mass * integral;
}

double QPFermiMomentum(double mass, double renormalized_chemical_potential)
{
    if (pow(renormalized_chemical_potential, 2.0) <= pow(mass, 2.0))
        return 0.0;

    return sqrt(pow(renormalized_chemical_potential, 2.0) - pow(mass, 2.0));
}

double QuarkZeroedGapEquation(double mass,
                              double up_scalar_density,
                              double down_scalar_density)
{

    double term = 2.0 * parameters.quark_model.G_S * CONST_HBAR_C
                  * (up_scalar_density + down_scalar_density);

    return mass - parameters.quark_model.bare_mass + term;
}

void QuarkSelfConsistentRenormalizedChemicalPotential(RenormChemPotSolutionParameters params,
                                                      double up_quark_mass,
                                                      double down_quark_mass,
                                                      double up_chemical_potential,
                                                      double down_chemical_potential,
                                                      double temperature,
                                                      double *return_up_renorm_chem_pot,
                                                      double *return_down_renorm_chem_pot)
{
    if (parameters.quark_model.G_V == 0.0){
        *return_up_renorm_chem_pot = up_chemical_potential;
        *return_down_renorm_chem_pot = down_chemical_potential;

        return;
    }

    // Set up parameters to be passed to helper function
    renorm_chem_pot_equation_input p;
    p.up_quark_mass = up_quark_mass;
    p.down_quark_mass = down_quark_mass;
    p.up_chemical_potential = up_chemical_potential;
    p.down_chemical_potential = down_chemical_potential;
    p.temperature = temperature;

    // Set dimension (number of equations|variables to solve|find)
    const int dimension = 2;

    gsl_multiroot_function f;
    f.f = &ZeroedRenormalizedQuarkChemicalPotentialEquation;
    f.n = dimension;
    f.params = (void *)&p;

    gsl_vector * initial_guess = gsl_vector_alloc(dimension);
    gsl_vector * return_results = gsl_vector_alloc(dimension);

    gsl_vector_set(initial_guess,
                   0,
                   sqrt(params.up_renorm_chem_pot_guess));
    gsl_vector_set(initial_guess,
                   1,
                   sqrt(params.down_renorm_chem_pot_guess));

    int status =
        MultidimensionalRootFinder(dimension,
                                   &f,
                                   initial_guess,
                                   params.abs_error,
                                   params.rel_error,
                                   params.max_iter,
                                   return_results);

    if (status != 0){
        printf("%s:%d: Something is wrong with the rootfinding.\n",
               __FILE__,
               __LINE__);
        abort();
    }

    // Save results in return variables,
    // taking care of the mappinps
    *return_up_renorm_chem_pot = pow(gsl_vector_get(return_results, 0), 2.0);
    *return_down_renorm_chem_pot = pow(gsl_vector_get(return_results, 1), 2.0);

    // Free vectors
    gsl_vector_free(initial_guess);
    gsl_vector_free(return_results);

    return;
}

int ZeroedRenormalizedQuarkChemicalPotentialEquation(const gsl_vector   *x,
                                                        void *params,
                                                        gsl_vector *return_values)
{
    renorm_chem_pot_equation_input * p =
        (renorm_chem_pot_equation_input *)params;

    const double up_renormalized_chemical_potential =
        pow(gsl_vector_get(x, 0), 2.0);
    const double down_renormalized_chemical_potential =
        pow(gsl_vector_get(x, 1), 2.0);

    double term = 2.0 * parameters.quark_model.G_V * CONST_HBAR_C
                  * (QuarkBarionicDensity(p->up_quark_mass,
                                          up_renormalized_chemical_potential,
                                          p->temperature)
                     + QuarkBarionicDensity(p->down_quark_mass,
                                            down_renormalized_chemical_potential,
                                            p->temperature));

    double up_mu_zeroed_eq = up_renormalized_chemical_potential
                             - p->up_chemical_potential
                             + term;

    double down_mu_zeroed_eq = down_renormalized_chemical_potential
                               - p->down_chemical_potential
                               + term;
/*
    printf(">%s:%d:\n>\tup_mu_guess: %f\n>\tdown_mu_guess: %f\n",
           __FILE__,
           __LINE__,
           up_renormalized_chemical_potential,
           down_renormalized_chemical_potential);
*/
    gsl_vector_set(return_values, 0, up_mu_zeroed_eq);
    gsl_vector_set(return_values, 1, down_mu_zeroed_eq);

    return GSL_SUCCESS;
}

typedef struct _TestMassAndRenormChemPot{
    double up_chemical_potential;
    double down_chemical_potential;

    double up_renorm_chem_pot;
    double down_renorm_chem_pot;
} TestMassAndRenormChemPot;

int TestMassAndRenormChemPotSimultaneousSolutionEquation(const gsl_vector   *x,
                                                         void *params,
                                                         gsl_vector *return_values);

int TestMassAndRenormChemPotSimultaneousSolution(double up_chemical_potential,
                                                 double down_chemical_potential,
                                                 double up_mass_guess,
                                                 double down_mass_guess,
                                                 double abs_error,
                                                 double rel_error,
                                                 int max_iter,
                                                 double * return_up_mass,
                                                 double * return_down_mass,
                                                 double * return_up_renorm_chem_pot,
                                                 double * return_down_renorm_chem_pot)
{
    // Set up parameters to be passed to helper function
    TestMassAndRenormChemPot params;
    params.up_chemical_potential = up_chemical_potential;
    params.down_chemical_potential = down_chemical_potential;

    // Set dimension (number of equations|variables to solve|find)
    const int dimension = 2;

    gsl_multiroot_function f;
    f.f = &TestMassAndRenormChemPotSimultaneousSolutionEquation;
    f.n = dimension;
    f.params = (void *)&params;

    gsl_vector * initial_guess = gsl_vector_alloc(dimension);
    gsl_vector * return_results = gsl_vector_alloc(dimension);

    gsl_vector_set(initial_guess,
                   0,
                   sqrt(up_mass_guess));
    gsl_vector_set(initial_guess,
                   1,
                   sqrt(down_mass_guess));

    int status =
        MultidimensionalRootFinder(dimension,
                                   &f,
                                   initial_guess,
                                   abs_error,
                                   rel_error,
                                   max_iter,
                                   return_results);

    if (status != 0){
        printf("%s:%d: Something is wrong with the rootfinding.\n",
               __FILE__,
               __LINE__);
        abort();
    }

    // Save results in return variables,
    // taking care of the mappinps
    *return_up_mass = pow(gsl_vector_get(return_results, 0), 2.0);
    *return_down_mass = pow(gsl_vector_get(return_results, 1), 2.0);

    *return_up_renorm_chem_pot = params.up_renorm_chem_pot;
    *return_down_renorm_chem_pot = params.down_renorm_chem_pot;

    // Free vectors
    gsl_vector_free(initial_guess);
    gsl_vector_free(return_results);

    return status;
}

int TestMassAndRenormChemPotSimultaneousSolutionEquation(const gsl_vector   *x,
                                                         void *params,
                                                         gsl_vector *return_values)
{
    const double up_mass = pow(gsl_vector_get(x, 0), 2.0);
    const double down_mass = pow(gsl_vector_get(x, 1), 2.0);

    TestMassAndRenormChemPot * p = (TestMassAndRenormChemPot *)params;

    double up_renorm_chem_pot;
    double down_renorm_chem_pot;

    QuarkSelfConsistentRenormalizedChemicalPotential(parameters.simultaneous_solution.renorm_chem_pot_solution,
                                                     up_mass,
                                                     down_mass,
                                                     p->up_chemical_potential,
                                                     p->down_chemical_potential,
                                                     parameters.variables.temperature,
                                                     &up_renorm_chem_pot,
                                                     &down_renorm_chem_pot);
    // save renormalized chemical potentials
    p->up_renorm_chem_pot = up_renorm_chem_pot;
    p->down_renorm_chem_pot = down_renorm_chem_pot;

    // Gap equations:
    double up_scalar_density =
        QuarkScalarDensity(parameters.variables.temperature,
                           up_mass,
                           up_renorm_chem_pot);

    double down_scalar_density =
        QuarkScalarDensity(parameters.variables.temperature,
                           down_mass,
                           down_renorm_chem_pot);

    double up_quark_zeroed_gap_eq =
        QuarkZeroedGapEquation(up_mass,
                               up_scalar_density,
                               down_scalar_density);

    double down_quark_zeroed_gap_eq =
        QuarkZeroedGapEquation(down_mass,
                               up_scalar_density,
                               down_scalar_density);

    gsl_vector_set(return_values, 0, up_quark_zeroed_gap_eq);
    gsl_vector_set(return_values, 1, down_quark_zeroed_gap_eq);

    return GSL_SUCCESS;
}
