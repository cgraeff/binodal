//
//  QuarkPhaseEOS.c
//  binodal
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

double QuarkFermiMomentum(double mass, double renormalized_chemical_potential)
{
    if (pow(renormalized_chemical_potential, 2.0) <= pow(mass, 2.0))
        return 0.0;

    return sqrt(pow(renormalized_chemical_potential, 2.0) - pow(mass, 2.0));
}

double QuarkDensity(double mass,
                    double renormalized_chemical_potential,
                    double temperature)
{
    double constant = NUM_Q_COLORS / (pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));

    if (temperature == 0){

        double fermi_momentum_3rd_power =
            pow(QuarkFermiMomentum(mass, renormalized_chemical_potential), 3.0);

        return constant * fermi_momentum_3rd_power / 3.0;
    }

    double integral =
    FermiDiracDistributionFromDensityIntegral(temperature,
                                              mass,
                                              renormalized_chemical_potential);

    return constant * integral;
}

double QuarkProtonFraction(double up_quark_density,
                           double down_quark_density)
{
    double proton_density = 2.0 * up_quark_density + down_quark_density;
    double neutron_density = 2.0 * down_quark_density + up_quark_density;

    return proton_density / (proton_density + neutron_density);
}

double QuarkScalarDensity(double temperature,
                          double mass,
                          double renorm_chem_pot)
{
    if (mass == 0){
        return 0;
    }

    double constant = - NUM_Q_COLORS
                        / (pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));

    if (temperature == 0){

        double fermi_momentum =
        QuarkFermiMomentum (mass, renorm_chem_pot);

        double upper_limit = F0(mass, parameters.quark.model.cutoff);
        double lower_limit = F0(mass, fermi_momentum);

        return constant * mass * (upper_limit - lower_limit);
    }

    double integral =
    FermiDiracDistributionIntegralFromScalarDensity(temperature,
                                                    mass,
                                                    renorm_chem_pot);

    return constant * mass * integral;
}

double QuarkZeroedGapEquation(double mass,
                              double up_scalar_density,
                              double down_scalar_density)
{

    double term = 2.0 * parameters.quark.model.G_S * CONST_HBAR_C
                  * (up_scalar_density + down_scalar_density);

    return mass - parameters.quark.model.bare_mass + term;
}

void QuarkVacuumMassDetermination(double * up_vacuum_mass,
                                  double * down_vacuum_mass)
{
    QuarkVacuumMassDeterminationParameters params =
    parameters.quark.vacuum_mass_determination;

    MultidimensionalRootFinderParams rootf_pars;
    rootf_pars.solver_type = gsl_multiroot_fsolver_dnewton;
    rootf_pars.abs_error = params.abs_error;
    rootf_pars.rel_error = params.rel_error;
    rootf_pars.max_iterations = params.max_iter;


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
                   sqrt(params.up_vacuum_mass_guess));
    gsl_vector_set(initial_guess,
                   1,
                   sqrt(params.down_vacuum_mass_guess));

    int status =
    MultidimensionalRootFinder(&f,
                               &rootf_pars,
                               dimension,
                               initial_guess,
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

typedef struct _therm_pot_free_gas_contrib_params{
    double mass;
    double temperature;
    double renormalized_chemical_potential;
} therm_pot_free_gas_contrib_params;

double QuarkThermodynamicPotentialFreeGasTermIntegrand(double momentum,
                                                       void * parameters);

double QuarkThermodynamicPotentialFreeGasTerm(double mass,
                                              double chemical_potential,
                                              double renorm_chem_pot,
                                              double temperature)
{
    double constant = - NUM_Q_COLORS * pow(CONST_HBAR_C, -3.0)
                        / pow(M_PI, 2.0);

    if (temperature == 0){

        double fermi_momentum = QuarkFermiMomentum(mass, renorm_chem_pot);

        double upper_limit = F_E(mass, parameters.quark.model.cutoff);
        double lower_limit =  F_E(mass, fermi_momentum);
        double F_diff =  upper_limit - lower_limit;

        return constant
               * (F_diff
                  + renorm_chem_pot * pow(fermi_momentum, 3.0) / 3.0);
    }

    therm_pot_free_gas_contrib_params p;
    p.mass = mass;
    p.renormalized_chemical_potential = renorm_chem_pot;
    p.temperature = temperature;

    gsl_function F;
    F.function = &QuarkThermodynamicPotentialFreeGasTermIntegrand;
    F.params = &p;

    double integral =
    OnedimensionalIntegrator(&F,
                             parameters.therm_pot_free_gas_integral);

    return constant * integral;
}

double QuarkThermodynamicPotentialFreeGasTermIntegrand(double momentum,
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

double QuarkThermodynamicPotential(double up_mass,
                                   double down_mass,
                                   double up_chemical_potential,
                                   double down_chemical_potential,
                                   double up_renormalized_chemical_potential,
                                   double down_renormalized_chemical_potential,
                                   double temperature)
{
    double up_free_term =
    QuarkThermodynamicPotentialFreeGasTerm(up_mass,
                                           up_chemical_potential,
                                           up_renormalized_chemical_potential,
                                           temperature);
    double down_free_term =
    QuarkThermodynamicPotentialFreeGasTerm(down_mass,
                                           down_chemical_potential,
                                           down_renormalized_chemical_potential,
                                           temperature);

    double up_scalar_density =
    QuarkScalarDensity(temperature,
                       up_mass,
                       up_renormalized_chemical_potential);

    double down_scalar_density =
    QuarkScalarDensity(temperature,
                       down_mass,
                       down_renormalized_chemical_potential);

    double scalar_term = parameters.quark.model.G_S * CONST_HBAR_C
                         * pow(up_scalar_density + down_scalar_density, 2.0);

    double up_barionic_density =
    QuarkDensity(up_mass,
                 up_renormalized_chemical_potential,
                 temperature);

    double down_barionic_density =
    QuarkDensity(down_mass,
                 down_renormalized_chemical_potential,
                 temperature);

    double vector_term =
    - parameters.quark.model.G_V * CONST_HBAR_C
    * pow(up_barionic_density + down_barionic_density, 2.0);

    return (up_free_term + down_free_term) + scalar_term + vector_term;
}

double QuarkPressure(double regularized_thermodynamic_potential,
                     double temperature)
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
    QuarkThermodynamicPotMinDetermParams params =
    parameters.quark.therm_pot_minimum;

    minimizer_function_parameters p;
    p.up_chemical_potential = 0.0;
    p.down_chemical_potential = 0.0;
    p.up_renormalized_chemical_potential = 0.0;
    p.down_renormalized_chemical_potential = 0.0;
    p.temperature = 0.0;

    const gsl_multimin_fminimizer_type * minimizer_type =
    gsl_multimin_fminimizer_nmsimplex2;
    const int minimizer_dimension = 2;

    gsl_multimin_fminimizer * minimizer =
        gsl_multimin_fminimizer_alloc(minimizer_type, minimizer_dimension);

    gsl_multimin_function f;
    f.f = &MinimizationHelperFunction;
    f.params = &p;
    f.n = minimizer_dimension;

    gsl_vector * guesses = gsl_vector_alloc(minimizer_dimension);
    gsl_vector * step_size = gsl_vector_alloc(minimizer_dimension);

    gsl_vector_set(guesses, 0, params.up_vacuum_mass_guess);
    gsl_vector_set(guesses, 1, params.down_vacuum_mass_guess);

    gsl_vector_set(step_size, 0, params.up_mass_step);
    gsl_vector_set(step_size, 1, params.down_mass_step);

    gsl_multimin_fminimizer_set(minimizer,
                                &f,
                                guesses,
                                step_size);

    int status;
    int iter = 0;
    do {

        if (iter > params.max_iter){
            printf("%s:%d: Max iterations reached.\n",
                   __FILE__,
                   __LINE__);

            abort();
        }

        status = gsl_multimin_fminimizer_iterate (minimizer);

        if (status == GSL_ENOPROG){
            printf("%s:%d: Can't improve the result beyond "
                   "current one (no progress).",
                   __FILE__,
                   __LINE__);
            abort();
        }

        double size = gsl_multimin_fminimizer_size(minimizer);

        status = gsl_multimin_test_size(size, params.tolerance);

        iter++;

    } while (status == GSL_CONTINUE);

    gsl_vector_free(guesses);
    gsl_vector_free(step_size);

    *up_mass_at_minimum = gsl_vector_get(minimizer->x, 0);
    *down_mass_at_minimum = gsl_vector_get(minimizer->x, 1);

    gsl_multimin_fminimizer_free (minimizer);

    return;
}

typedef struct _renorm_chem_pot_equation_input{
    double up_quark_mass;
    double down_quark_mass;
    double up_chemical_potential;
    double down_chemical_potential;
    double temperature;
} renorm_chem_pot_equation_input;

int QuarkSelfConsistentRenormChemPot(double up_quark_mass,
                                      double down_quark_mass,
                                      double up_chemical_potential,
                                      double down_chemical_potential,
                                      double temperature,
                                      double *return_up_renorm_chem_pot,
                                      double *return_down_renorm_chem_pot)
{
    *return_up_renorm_chem_pot = NAN;
    *return_down_renorm_chem_pot = NAN;

    if (parameters.quark.model.G_V == 0.0){
        *return_up_renorm_chem_pot = up_chemical_potential;
        *return_down_renorm_chem_pot = down_chemical_potential;

        return 0;
    }

    QuarkRenormChemPotSolutionParameters params =
    parameters.quark.renorm_chem_pot_solution;

    // Set up parameters to be passed to helper function
    renorm_chem_pot_equation_input p;
    p.up_quark_mass = up_quark_mass;
    p.down_quark_mass = down_quark_mass;
    p.up_chemical_potential = up_chemical_potential;
    p.down_chemical_potential = down_chemical_potential;
    p.temperature = temperature;

    // Set dimension (number of equations|variables to solve|find)
    const int dimension = 2;

    MultidimensionalRootFinderParams rootf_pars;
    rootf_pars.solver_type = gsl_multiroot_fsolver_dnewton;
    rootf_pars.abs_error = params.abs_error;
    rootf_pars.rel_error = params.rel_error;
    rootf_pars.max_iterations = params.max_iter;

    gsl_multiroot_function f;
    f.f = &ZeroedRenormalizedQuarkChemPotEquation;
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
    MultidimensionalRootFinder(&f,
                               &rootf_pars,
                               dimension,
                               initial_guess,
                               return_results);

    if (status != 0){
        if (options.abort_on_error){
            printf("%s:%d: Something is wrong with the rootfinding.\n",
                   __FILE__,
                   __LINE__);
            abort();
        }

        return -1;
    }

    // Save results in return variables,
    // taking care of the mappinps
    *return_up_renorm_chem_pot = pow(gsl_vector_get(return_results, 0), 2.0);
    *return_down_renorm_chem_pot = pow(gsl_vector_get(return_results, 1), 2.0);

    // Free vectors
    gsl_vector_free(initial_guess);
    gsl_vector_free(return_results);

    return 0;
}

int ZeroedRenormalizedQuarkChemPotEquation(const gsl_vector   *x,
                                           void *params,
                                           gsl_vector *return_values)
{
    renorm_chem_pot_equation_input * p =
        (renorm_chem_pot_equation_input *)params;

    const double up_renormalized_chemical_potential =
        pow(gsl_vector_get(x, 0), 2.0);
    const double down_renormalized_chemical_potential =
        pow(gsl_vector_get(x, 1), 2.0);

    double term =
    2.0 * parameters.quark.model.G_V * CONST_HBAR_C
    * (QuarkDensity(p->up_quark_mass,
                    up_renormalized_chemical_potential,
                    p->temperature)
       + QuarkDensity(p->down_quark_mass,
                      down_renormalized_chemical_potential,
                      p->temperature)
       );

    double up_mu_zeroed_eq = up_renormalized_chemical_potential
                             - p->up_chemical_potential
                             + term;

    double down_mu_zeroed_eq = down_renormalized_chemical_potential
                               - p->down_chemical_potential
                               + term;

    gsl_vector_set(return_values, 0, up_mu_zeroed_eq);
    gsl_vector_set(return_values, 1, down_mu_zeroed_eq);

    return GSL_SUCCESS;
}

int QuarkMassAndRenormChemPotSolEquation(const gsl_vector   *x,
                                         void *params,
                                         gsl_vector *return_values);

int QuarkMassAndRenormChemPotSolution(double up_chemical_potential,
                                      double down_chemical_potential,
                                      double up_mass_guess,
                                      double down_mass_guess,
                                      double * return_up_mass,
                                      double * return_down_mass,
                                      double * return_up_renorm_chem_pot,
                                      double * return_down_renorm_chem_pot)
{
    *return_up_mass = NAN;
    *return_down_mass = NAN;
    *return_up_renorm_chem_pot = NAN;
    *return_down_renorm_chem_pot = NAN;

    QuarkMassAndRenormChemPotSolParams params =
    parameters.quark.mass_and_renorm_chem_pot_solution;

    // Set up parameters to be passed to helper function
    quark_mass_and_renorm_chem_pot_input_params p;
    p.up_chemical_potential = up_chemical_potential;
    p.down_chemical_potential = down_chemical_potential;
    p.up_renorm_chem_pot = NAN;
    p.down_renorm_chem_pot = NAN;

    // Set dimension (number of equations|variables to solve|find)
    const int dimension = 2;

    MultidimensionalRootFinderParams rootf_pars;
    rootf_pars.solver_type = gsl_multiroot_fsolver_dnewton;
    rootf_pars.abs_error = params.abs_error;
    rootf_pars.rel_error = params.rel_error;
    rootf_pars.max_iterations = params.max_iter;

    gsl_multiroot_function f;
    f.f = &QuarkMassAndRenormChemPotSolEquation;
    f.n = dimension;
    f.params = (void *)&p;

    gsl_vector * initial_guess = gsl_vector_alloc(dimension);
    gsl_vector * return_results = gsl_vector_alloc(dimension);

    gsl_vector_set(initial_guess,
                   0,
                   sqrt(up_mass_guess));
    gsl_vector_set(initial_guess,
                   1,
                   sqrt(down_mass_guess));

    int status =
    MultidimensionalRootFinder(&f,
                               &rootf_pars,
                               dimension,
                               initial_guess,
                               return_results);

    // If no progress is made, we may have reached chiral restoration
    if (status == -1){
            gsl_vector_set(initial_guess,
                           0,
                           0.0);
            gsl_vector_set(initial_guess,
                           1,
                           0.0);

            status =
            MultidimensionalRootFinder(&f,
                                       &rootf_pars,
                                       dimension,
                                       initial_guess,
                                       return_results);

    }

    if (status){
        if (options.abort_on_error){
            printf("%s:%d: Something is wrong with the rootfinding (error: %d).\n",
                   __FILE__,
                   __LINE__,
                   status);
            printf("up_chem_pot: %f, down_chem_pot: %f\n",
                   up_chemical_potential,
                   down_chemical_potential);
            abort();
        }

        return -1;
    }

    // Save results in return variables,
    // taking care of the mappinps
    *return_up_mass = pow(gsl_vector_get(return_results, 0), 2.0);
    *return_down_mass = pow(gsl_vector_get(return_results, 1), 2.0);

    *return_up_renorm_chem_pot = p.up_renorm_chem_pot;
    *return_down_renorm_chem_pot = p.down_renorm_chem_pot;

    // Free vectors
    gsl_vector_free(initial_guess);
    gsl_vector_free(return_results);

    return status;
}

int QuarkMassAndRenormChemPotSolEquation(const gsl_vector   *x,
                                         void *params,
                                         gsl_vector *return_values)
{
    const double up_mass = pow(gsl_vector_get(x, 0), 2.0);
    const double down_mass = pow(gsl_vector_get(x, 1), 2.0);

    quark_mass_and_renorm_chem_pot_input_params * p =
    (quark_mass_and_renorm_chem_pot_input_params *)params;

    double up_renorm_chem_pot;
    double down_renorm_chem_pot;

    int status =
    QuarkSelfConsistentRenormChemPot(up_mass,
                                     down_mass,
                                     p->up_chemical_potential,
                                     p->down_chemical_potential,
                                     parameters.variables.temperature,
                                     &up_renorm_chem_pot,
                                     &down_renorm_chem_pot);

    if (status){
        printf("%s:%d: Problems in renormalized chemical "
               "potential determination.\n",
               __FILE__,
               __LINE__);
        abort();
    }

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

double MyAdapterFunction(double val, void *params)
{
    const int dimension = 2;

    gsl_vector * x = gsl_vector_alloc(dimension);
    gsl_vector * return_values = gsl_vector_alloc(dimension);

    const double up_mass_sq = sqrt(val);
    const double down_mass_sq = sqrt(val);

    gsl_vector_set(x, 0, up_mass_sq);
    gsl_vector_set(x, 1, down_mass_sq);

    int status =
    QuarkMassAndRenormChemPotSolEquation(x,
                                         params,
                                         return_values);

    if (status != GSL_SUCCESS){
        printf("%s:%d: Problems.\n", __FILE__, __LINE__);
        abort();
    }

    double result = gsl_vector_get(return_values, 0);

    gsl_vector_free(x);
    gsl_vector_free(return_values);

    return result;
}

int MySign(double value)
{
    int sgn;
    if (value < 0){
        sgn = -1;
    }
    else if (value > 0){
        sgn = 1;
    }
    else {
        sgn = 0;
    }

    return sgn;
}

int QuarkMassAndRenormChemPotSolutionBissection(double up_chemical_potential,
                                                double down_chemical_potential,
                                                double temperature,
                                                double * return_up_mass,
                                                double * return_down_mass,
                                                double * return_up_renorm_chem_pot,
                                                double * return_down_renorm_chem_pot)
{
    // If no solutions are found, we probably have reached
    // the chiral restoration, so we will calculate this solution
    // and set it to the return variables. If another solution
    // is found, it will be overwritten. The return status will be
    // different from zero if no other solutions are found.
    *return_up_mass = 0.0;
    *return_down_mass = 0.0;

    double up_renorm_chem_pot;
    double down_renorm_chem_pot;

    QuarkSelfConsistentRenormChemPot(*return_up_mass,
                                     *return_down_mass,
                                     up_chemical_potential,
                                     down_chemical_potential,
                                     temperature,
                                     &up_renorm_chem_pot,
                                     &down_renorm_chem_pot);

    *return_up_renorm_chem_pot = up_renorm_chem_pot;
    *return_down_renorm_chem_pot = down_renorm_chem_pot;

    //
    // Search for other solutions
    //

    QuarkMassAndRenormChemPotSolParBissec bissection_pars =
    parameters.quark.quark_mass_and_renorm_chem_pot_bissec_params;

     // Set up parameters to be passed to helper function
    quark_mass_and_renorm_chem_pot_input_params p;
    p.up_chemical_potential = up_chemical_potential;
    p.down_chemical_potential = down_chemical_potential;
    p.up_renorm_chem_pot = NAN;
    p.down_renorm_chem_pot = NAN;

    UnidimensionalRootFindingParameters params;
    params.abs_error = 1.0E-4;
    params.rel_error = 1.0E-4;
    params.max_iterations = 1000;

    gsl_function F;
    F.function = &MyAdapterFunction;
    F.params = &p;

    bool no_solutions = true;

    //
    // Loop for search of multiple solutions
    //

    double min_potential = INFINITY;

    double point = bissection_pars.min_mass;

    while (true){

        //printf("mass: %f\n", point);

        int sgn = MySign(MyAdapterFunction(point, (void *)&p));

        double lower_bound;
        int next_point_sgn;

        do {
            lower_bound = point;
            point += bissection_pars.mass_step;
            next_point_sgn = MySign(MyAdapterFunction(point, (void *)&p));

            if (point > bissection_pars.max_mass)
                return no_solutions; // Will be false if there are solutions

        }while (sgn == next_point_sgn);

        no_solutions = false;

        double result;

        if (next_point_sgn == 0){
            result = point;
        }
        else{
            double upper_bound = point;

            params.lower_bound = lower_bound;
            params.upper_bound = upper_bound;

            int status =
            UnidimensionalRootFinder(&F,
                                     params,
                                     &result);

            if (status){
                printf("%s:%d: UnidimensionalRootFinder didn't work.\n",
                       __FILE__,
                       __LINE__);
                abort();

                return -1;
            }
        }

        QuarkSelfConsistentRenormChemPot(result,
                                         result,
                                         up_chemical_potential,
                                         down_chemical_potential,
                                         temperature,
                                         &up_renorm_chem_pot,
                                         &down_renorm_chem_pot);

        double potential =
        QuarkThermodynamicPotential(result,
                                    result,
                                    up_chemical_potential,
                                    down_chemical_potential,
                                    up_renorm_chem_pot,
                                    down_renorm_chem_pot,
                                    temperature);

        if (potential < min_potential){

            min_potential = potential;

            *return_up_mass = result;
            *return_down_mass = result;
            *return_up_renorm_chem_pot = up_renorm_chem_pot;
            *return_down_renorm_chem_pot = down_renorm_chem_pot;

        }
    }

    // Never reached
    return 0;
}

double QuarkPhaseAsymmetry(double up_quark_density,
                           double down_quark_density)
{
    return 3.0 * (down_quark_density - up_quark_density)
           / (up_quark_density + down_quark_density);
}

