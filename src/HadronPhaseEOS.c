//
//  HadronPhaseEOS.c
//  binodal
//
//  Created by Clebson Graeff on 2016-02-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>

#include "libdatafun/libdatafun.h"

#include "CommandlineOptions.h"
#include "Parameters.h"
#include "Constants.h"
#include "Functions.h"

double HadronZeroedGapEquation(double mass,
                               void * params)
{
    hadron_gap_eq_input_params * p = (hadron_gap_eq_input_params *)params;

    double barionic_density = p->proton_density + p->neutron_density;
    double rho_3 = p->proton_density - p->neutron_density;

    double scalar_density =
    HadronScalarDensity(mass,
                        p->proton_fermi_momentum,
                        parameters.hadron.model.cutoff)
    + HadronScalarDensity(mass,
                          p->neutron_fermi_momentum,
                          parameters.hadron.model.cutoff);

    double gap_1st_term = parameters.hadron.model.G_S * scalar_density;
    double gap_2nd_term = - parameters.hadron.model.G_SV
                            * scalar_density
                            * pow(barionic_density, 2.0);
    double gap_3rd_term = - parameters.hadron.model.G_SRHO
                            * scalar_density
                            * pow(rho_3, 2.0);

    return mass
           + 2.0 * CONST_HBAR_C
             * (gap_1st_term + gap_2nd_term + gap_3rd_term)
           - parameters.hadron.model.bare_mass;
}

double HadronProtonFraction(double proton_barionic_density,
                            double neutron_barionic_density)
{
    return proton_barionic_density
           / (proton_barionic_density + neutron_barionic_density);
}

double HadronScalarDensity(double mass,
                           double fermi_momentum,
                           double cutoff)
{
    if (mass == 0.0){
        return 0.0;
    }

	return pow(CONST_HBAR_C, -3.0)
           * (mass / pow(M_PI, 2.0))
           * (F0(mass, fermi_momentum) - F0(mass, cutoff));
}

double HadronVacuumScalarDensity()
{
    return 2.0 * pow(CONST_HBAR_C, -3.0)
           * (parameters.hadron.model.nucleon_mass / pow(M_PI, 2.0))
           * (F0(parameters.hadron.model.nucleon_mass, 0.0)
              - F0(parameters.hadron.model.nucleon_mass,
                   parameters.hadron.model.cutoff));
}

double ProtonChemicalPotentialEquation(double proton_fermi_momentum,
                                       double scalar_density,
                                       double mass,
                                       double proton_density,
                                       double neutron_density)
{
    double E = sqrt(pow(mass, 2.0) + pow(proton_fermi_momentum, 2.0));

    double barionic_density = proton_density + neutron_density;
    double rho_3 = (proton_density - neutron_density);

    double rho_terms = (parameters.hadron.model.G_V * barionic_density
                        + parameters.hadron.model.G_SV * barionic_density
                          * pow(scalar_density, 2.0)
                        + parameters.hadron.model.G_RHO * rho_3
                        + parameters.hadron.model.G_VRHO * pow(rho_3, 2.0)
                          * barionic_density
                        + parameters.hadron.model.G_VRHO * pow(barionic_density, 2.0)
                          * rho_3
                        + parameters.hadron.model.G_SRHO * pow(scalar_density, 2.0)
                          * rho_3)
                       * 2.0 * CONST_HBAR_C;

    return E + rho_terms;
}

double NeutronChemicalPotentialEquation(double neutron_fermi_momentum,
                                        double scalar_density,
                                        double mass,
                                        double proton_density,
                                        double neutron_density)
{
    double E = sqrt(pow(mass, 2.0)
                    + pow(neutron_fermi_momentum, 2.0));

    double barionic_density = proton_density + neutron_density;
    double rho_3 = (proton_density - neutron_density);

    double rho_terms = (parameters.hadron.model.G_V * barionic_density
                        + parameters.hadron.model.G_SV * barionic_density
                          * pow(scalar_density, 2.0)
                        - parameters.hadron.model.G_RHO * rho_3
                        + parameters.hadron.model.G_VRHO * pow(rho_3, 2.0)
                          * barionic_density
                        - parameters.hadron.model.G_VRHO * pow(barionic_density, 2.0)
                          * rho_3
                        - parameters.hadron.model.G_SRHO * pow(scalar_density, 2.0)
                          * rho_3)
                       * 2.0 * CONST_HBAR_C;

    return E + rho_terms;
}


double HadronKinecticEnergyDensity(double mass,
                                   double proton_fermi_momentum,
                                   double neutron_fermi_momentum)
{
    double proton_kinectic_energy = (NUM_H_COLORS  / pow(M_PI, 2.0))
                                    * (F2(mass, proton_fermi_momentum)
                                       - F2(mass, parameters.hadron.model.cutoff));

    double neutron_kinectic_energy = (NUM_H_COLORS / pow(M_PI, 2.0))
                                     * (F2(mass, neutron_fermi_momentum)
                                        - F2(mass, parameters.hadron.model.cutoff));

    return (proton_kinectic_energy + neutron_kinectic_energy)
           / (pow(CONST_HBAR_C, 3.0));
}

double HadronVacuumKinecticEnergyDensity()
{
    return 2.0 * (NUM_H_COLORS / pow(M_PI, 2.0)) * (pow(CONST_HBAR_C, -3.0))
           * (F2(parameters.hadron.model.nucleon_mass, 0)
              - F2(parameters.hadron.model.nucleon_mass,
                   parameters.hadron.model.cutoff));
}

double HadronVacuumEnergyDensity()
{
    double scalar_density_0 = HadronVacuumScalarDensity();

    return HadronVacuumKinecticEnergyDensity()
           + parameters.hadron.model.bare_mass * scalar_density_0
           - parameters.hadron.model.G_S * pow(scalar_density_0, 2.0) * CONST_HBAR_C;
}

double HadronEnergyDensity(double pressure,
                           double proton_chemical_potential,
                           double neutron_chemical_potential,
                           double proton_density,
                           double neutron_density)
{
    return - pressure
           + proton_chemical_potential * proton_density
           + neutron_chemical_potential * neutron_density;
}

double HadronPressure(double termodynamic_potential)
{
    return - termodynamic_potential;
}

double HadronThermodynamicPotential(double total_scalar_density,
                                    double barionic_density,
                                    double proton_density,
                                    double neutron_density,
                                    double proton_chemical_potential,
                                    double neutron_chemical_potential,
                                    double kinectic_energy_density)
{
    double rho_3 = proton_density - neutron_density;

    double omega = kinectic_energy_density
                   + parameters.hadron.model.bare_mass * total_scalar_density;

    omega += - proton_chemical_potential * proton_density
             - neutron_chemical_potential * neutron_density;

    omega += (- parameters.hadron.model.G_S * pow(total_scalar_density, 2.0)
              + parameters.hadron.model.G_V * pow(barionic_density, 2.0)
              + parameters.hadron.model.G_SV
                * pow(total_scalar_density * barionic_density, 2.0)
              + parameters.hadron.model.G_RHO * pow(rho_3, 2.0)
              + parameters.hadron.model.G_VRHO
                * pow(barionic_density * rho_3, 2.0)
              + parameters.hadron.model.G_SRHO
                * pow(total_scalar_density * rho_3, 2.0))
             * CONST_HBAR_C;

    return omega;
}

double HadronFermiMomentum(double density)
{
    return CONST_HBAR_C * pow(3.0 * pow(M_PI, 2.0) * density, 1.0 / 3.0);
}

double HadronPhaseAsymmetry(double proton_density, double neutron_density)
{
    double barionic_density = proton_density + neutron_density;
    double proton_fraction = proton_density / barionic_density;

    return 1.0 - 2.0 * proton_fraction;
}


// Solve for chempot

int
HadronMassAndDensitiesSolutionEquationZeroMassCase(const gsl_vector   *x,
                                                   void *params,
                                                   gsl_vector *return_values);

int
HadronMassAndDensitiesSolutionEquationZeroProtonDensCase(const gsl_vector   *x,
                                                         void *params,
                                                         gsl_vector *return_values);

int
HadronMassAndDensitiesSolutionEquationZeroNeutronDensCase(const gsl_vector   *x,
                                                         void *params,
                                                         gsl_vector *return_values);


int HadronMassAndDensitiesSolution(double proton_chemical_potential,
                                   double neutron_chemical_potential,
                                   double hadron_mass_guess,
                                   double proton_density_guess,
                                   double neutron_density_guess,
                                   double * return_mass,
                                   double * return_proton_density,
                                   double * return_neutron_density)
{
    *return_mass = NAN;
    *return_proton_density = NAN;
    *return_neutron_density = NAN;

    double barionic_chemical_potential =
    BarionicChemicalPotential(proton_chemical_potential,
                              neutron_chemical_potential);

    double isovector_chemical_potential =
    BarionicChemicalPotential(proton_chemical_potential,
                              neutron_chemical_potential);

    double zero_dens_tolerance = 1.0E-4; // move to parameters
    bool debug = false;

    // We will try with the following order of solver, until one
    // of them works
    const gsl_multiroot_fsolver_type * solver_types[4] =
        {
            gsl_multiroot_fsolver_dnewton,
            gsl_multiroot_fsolver_hybrids,
            gsl_multiroot_fsolver_broyden,
            gsl_multiroot_fsolver_hybrid
        };

    HadronMassAndDensitiesSolutionParams params =
    parameters.hadron.mass_and_densities_solution;

    MultidimensionalRootFinderParams rootf_pars;
    rootf_pars.abs_error = params.abs_error;
    rootf_pars.rel_error = params.rel_error;
    rootf_pars.max_iterations = params.max_iter;

    hadron_mass_and_renorm_chem_pot_input_params p;
    p.proton_chemical_potential = proton_chemical_potential;
    p.neutron_chemical_potential = neutron_chemical_potential;

    /* Begin cases */

    // If all guesses are below the tolerances, there is
    // nothing to do
    if (hadron_mass_guess < params.zero_mass_tolerance
        && proton_density_guess < zero_dens_tolerance
        && neutron_density_guess < zero_dens_tolerance){

        *return_mass = 0.0;
        *return_proton_density = 0.0;
        *return_neutron_density = 0.0;

        return 0;
    }

    // Case where only mass goes to zero
    if (hadron_mass_guess < params.zero_mass_tolerance
        && proton_density_guess > zero_dens_tolerance
        && neutron_density_guess > zero_dens_tolerance){

        if (debug)
            printf("bar chem pot: %f\n"
                   "iso chem pot: %f\n"
                   "hadron mass guess: %f\n"
                   "proton density guess: %f\n"
                   "neutron density guess: %f\n",
                   barionic_chemical_potential,
                   isovector_chemical_potential,
                   hadron_mass_guess,
                   proton_density_guess,
                   neutron_density_guess);

        // Set dimension (number of equations|variables to solve|find)
        const int dimension = 2;

        rootf_pars.solver_type = solver_types[0];

        gsl_multiroot_function f;
        f.f = &HadronMassAndDensitiesSolutionEquationZeroMassCase;
        f.n = dimension;
        f.params = (void *)&p;

        gsl_vector * initial_guess = gsl_vector_alloc(dimension);
        gsl_vector * return_results = gsl_vector_alloc(dimension);

        gsl_vector_set(initial_guess, 0, sqrt(proton_density_guess));
        gsl_vector_set(initial_guess, 1, sqrt(neutron_density_guess));

        int status =
        MultidimensionalRootFinder(&f,
                                   &rootf_pars,
                                   dimension,
                                   initial_guess,
                                   return_results);

       if (status != 0){
            // Free vectors
            gsl_vector_free(initial_guess);
            gsl_vector_free(return_results);

            if (options.abort_on_error){

                printf("%s:%d:"
                       "Abort on error enabled, aborting due to error.\n",
                       __FILE__, __LINE__);

                abort();
            }

            return -1;
        }

        // Save results in return variables,
        // taking care of the mappings
        *return_mass = 0.0;
        *return_proton_density = pow(gsl_vector_get(return_results, 0), 2.0);
        *return_neutron_density = pow(gsl_vector_get(return_results, 1), 2.0);

        // Free vectors
        gsl_vector_free(initial_guess);
        gsl_vector_free(return_results);

        return 0;
    }

    // Handle case where proton density is below zero tolerance
    if (hadron_mass_guess > params.zero_mass_tolerance
        && proton_density_guess < zero_dens_tolerance
        && neutron_density_guess > zero_dens_tolerance){

        if (debug)
            printf("bar chem pot: %f\n"
                   "iso chem pot: %f\n"
                   "hadron mass guess: %f\n"
                   "proton density guess: %f\n"
                   "neutron density guess: %f\n",
                   barionic_chemical_potential,
                   isovector_chemical_potential,
                   hadron_mass_guess,
                   proton_density_guess,
                   neutron_density_guess);

        // Set dimension (number of equations|variables to solve|find)
        const int dimension = 2;

        rootf_pars.solver_type = solver_types[0];

        gsl_multiroot_function f;
        f.f = &HadronMassAndDensitiesSolutionEquationZeroProtonDensCase;
        f.n = dimension;
        f.params = (void *)&p;

        gsl_vector * initial_guess = gsl_vector_alloc(dimension);
        gsl_vector * return_results = gsl_vector_alloc(dimension);

        gsl_vector_set(initial_guess, 0, sqrt(hadron_mass_guess));
        gsl_vector_set(initial_guess, 1, sqrt(neutron_density_guess));

        int status =
        MultidimensionalRootFinder(&f,
                                   &rootf_pars,
                                   dimension,
                                   initial_guess,
                                   return_results);

       if (status != 0){
            // Free vectors
            gsl_vector_free(initial_guess);
            gsl_vector_free(return_results);

            if (options.abort_on_error){

                printf("%s:%d:"
                       "Abort on error enabled, aborting due to error.\n",
                       __FILE__, __LINE__);

                abort();
            }

            return -1;
        }

        // Save results in return variables,
        // taking care of the mappings
        *return_mass = pow(gsl_vector_get(return_results, 0), 2.0);
        *return_proton_density = 0.0;
        *return_neutron_density = pow(gsl_vector_get(return_results, 1), 2.0);

        // Free vectors
        gsl_vector_free(initial_guess);
        gsl_vector_free(return_results);

        return 0;
    }

    // Handle case where neutron density is below zero tolerance
    if (hadron_mass_guess > params.zero_mass_tolerance
        && proton_density_guess > zero_dens_tolerance
        && neutron_density_guess < zero_dens_tolerance){

        if (debug)
            printf("bar chem pot: %f\n"
                   "iso chem pot: %f\n"
                   "hadron mass guess: %f\n"
                   "proton density guess: %f\n"
                   "neutron density guess: %f\n",
                   barionic_chemical_potential,
                   isovector_chemical_potential,
                   hadron_mass_guess,
                   proton_density_guess,
                   neutron_density_guess);

        // Set dimension (number of equations|variables to solve|find)
        const int dimension = 2;

        rootf_pars.solver_type = solver_types[0];

        gsl_multiroot_function f;
        f.f = &HadronMassAndDensitiesSolutionEquationZeroNeutronDensCase;
        f.n = dimension;
        f.params = (void *)&p;

        gsl_vector * initial_guess = gsl_vector_alloc(dimension);
        gsl_vector * return_results = gsl_vector_alloc(dimension);

        gsl_vector_set(initial_guess, 0, sqrt(hadron_mass_guess));
        gsl_vector_set(initial_guess, 1, sqrt(proton_density_guess));

        int status =
        MultidimensionalRootFinder(&f,
                                   &rootf_pars,
                                   dimension,
                                   initial_guess,
                                   return_results);

       if (status != 0){
            // Free vectors
            gsl_vector_free(initial_guess);
            gsl_vector_free(return_results);

            if (options.abort_on_error){

                printf("%s:%d:"
                       "Abort on error enabled, aborting due to error.\n",
                       __FILE__, __LINE__);

                abort();
            }

            return -1;
        }

        // Save results in return variables,
        // taking care of the mappings
        *return_mass = pow(gsl_vector_get(return_results, 0), 2.0);
        *return_proton_density = pow(gsl_vector_get(return_results, 1), 2.0);
        *return_neutron_density = 0.0;

        // Free vectors
        gsl_vector_free(initial_guess);
        gsl_vector_free(return_results);

        return 0;
    }

    // Standard case where no variables are near zero
    int status = -1;
    if (hadron_mass_guess > params.zero_mass_tolerance
        && proton_density_guess > zero_dens_tolerance
        && neutron_density_guess > zero_dens_tolerance){

        // Set dimension (number of equations|variables to solve|find)
        const int dimension = 3;

        gsl_multiroot_function f;
        f.f = &HadronMassAndDensitiesSolutionEquation;
        f.n = dimension;
        f.params = (void *)&p;

        gsl_vector * initial_guess = gsl_vector_alloc(dimension);
        gsl_vector * return_results = gsl_vector_alloc(dimension);

        int selected_solver = 0;
        rootf_pars.solver_type = solver_types[selected_solver];
        do {

            gsl_vector_set(initial_guess,
                           0,
                           sqrt(hadron_mass_guess));

            gsl_vector_set(initial_guess,
                           1,
                           sqrt(proton_density_guess));
            gsl_vector_set(initial_guess,
                           2,
                           sqrt(neutron_density_guess));

            status =
            MultidimensionalRootFinder(&f,
                                       &rootf_pars,
                                       dimension,
                                       initial_guess,
                                       return_results);

            if (!status){

                // Save results in return variables,
                // taking care of the mappinps
                *return_mass = pow(gsl_vector_get(return_results, 0), 2.0);
                *return_proton_density = pow(gsl_vector_get(return_results, 1), 2.0);
                *return_neutron_density = pow(gsl_vector_get(return_results, 2), 2.0);

                if (selected_solver){
                    printf("mass: %f\n", *return_mass);
                    printf("pdens: %f\n", *return_proton_density);
                    printf("ndens: %f\n", *return_neutron_density);
                }

                // Free vectors
                gsl_vector_free(initial_guess);
                gsl_vector_free(return_results);

                return 0;

            }

            gsl_vector_set(initial_guess, 0, hadron_mass_guess);
            gsl_vector_set(initial_guess, 1, proton_density_guess);
            gsl_vector_set(initial_guess, 2, neutron_density_guess);

            selected_solver++;

        } while (status && selected_solver <= 4);

        printf("¡not good!\n");
        printf("%s:%d\n", __FILE__, __LINE__);
        printf("hadron mass guess: %f\n"
               "proton density guess: %f\n"
               "neutron density guess: %f\n",
               hadron_mass_guess,
               proton_density_guess,
               neutron_density_guess);

        // Free vectors
        gsl_vector_free(initial_guess);
        gsl_vector_free(return_results);

        return -1;
    }

    printf("%s:%d ", __FILE__, __LINE__);
    printf("Problems in management of possible calculation scenarios.\n");
    printf("bar chem pot: %f\n"
           "iso chem pot: %f\n"
           "hadron mass guess: %f\n"
           "proton density guess: %f\n"
           "neutron density guess: %f\n",
           barionic_chemical_potential,
           isovector_chemical_potential,
           hadron_mass_guess,
           proton_density_guess,
           neutron_density_guess);

    abort();

    return -1;
}

int
HadronMassAndDensitiesSolutionEquationZeroMassCase(const gsl_vector   *x,
                                                   void *params,
                                                   gsl_vector *return_values)
{
    gsl_vector * guesses = gsl_vector_alloc(x->size + 1);
    gsl_vector * results = gsl_vector_alloc(x->size + 1);

    double mass = 0.0;

    gsl_vector_set(guesses, 0, mass);
    gsl_vector_set(guesses, 1, gsl_vector_get(x, 0));
    gsl_vector_set(guesses, 2, gsl_vector_get(x, 1));

    HadronMassAndDensitiesSolutionEquation(guesses,
                                           params,
                                           results);

    gsl_vector_set(return_values, 0, gsl_vector_get(results, 1));
    gsl_vector_set(return_values, 1, gsl_vector_get(results, 2));

    gsl_vector_free(guesses);
    gsl_vector_free(results);

    return 0;
}

int
HadronMassAndDensitiesSolutionEquationZeroProtonDensCase(const gsl_vector   *x,
                                                         void *params,
                                                         gsl_vector *return_values)
{
    gsl_vector * guesses = gsl_vector_alloc(x->size + 1);
    gsl_vector * results = gsl_vector_alloc(x->size + 1);

    double proton_dens = 0.0;

    gsl_vector_set(guesses, 0, gsl_vector_get(x, 0));
    gsl_vector_set(guesses, 1, proton_dens);
    gsl_vector_set(guesses, 2, gsl_vector_get(x, 1));

    HadronMassAndDensitiesSolutionEquation(guesses,
                                           params,
                                           results);

    gsl_vector_set(return_values, 0, gsl_vector_get(results, 0));
    gsl_vector_set(return_values, 1, gsl_vector_get(results, 2));

    gsl_vector_free(guesses);
    gsl_vector_free(results);

    return 0;
}

int
HadronMassAndDensitiesSolutionEquationZeroNeutronDensCase(const gsl_vector   *x,
                                                          void *params,
                                                          gsl_vector *return_values)
{
    gsl_vector * guesses = gsl_vector_alloc(x->size + 1);
    gsl_vector * results = gsl_vector_alloc(x->size + 1);

    double neutron_dens = 0.0;

    gsl_vector_set(guesses, 0, gsl_vector_get(x, 0));
    gsl_vector_set(guesses, 1, gsl_vector_get(x, 1));
    gsl_vector_set(guesses, 2, neutron_dens);

    HadronMassAndDensitiesSolutionEquation(guesses,
                                           params,
                                           results);

    gsl_vector_set(return_values, 0, gsl_vector_get(results, 0));
    gsl_vector_set(return_values, 1, gsl_vector_get(results, 1));

    gsl_vector_free(guesses);
    gsl_vector_free(results);

    return 0;
}

int HadronMassAndDensitiesSolutionEquation(const gsl_vector   *x,
                                           void *params,
                                           gsl_vector *return_values)
{
    const double mass = pow(gsl_vector_get(x, 0), 2.0);
    const double proton_density = pow(gsl_vector_get(x, 1), 2.0);
    const double neutron_density = pow (gsl_vector_get(x, 2), 2.0);

    hadron_mass_and_renorm_chem_pot_input_params * p =
    (hadron_mass_and_renorm_chem_pot_input_params *) params;

    double proton_fermi_momentum = HadronFermiMomentum(proton_density);
    double neutron_fermi_momentum = HadronFermiMomentum(neutron_density);

    hadron_gap_eq_input_params  gap_params;
    gap_params.proton_fermi_momentum = proton_fermi_momentum;
    gap_params.neutron_fermi_momentum = neutron_fermi_momentum;
    gap_params.proton_density = proton_density;
    gap_params.neutron_density = neutron_density;

    double gap_equation = HadronZeroedGapEquation(mass,
                                                  (void *)&gap_params);

    //calculate scalar densities
    double proton_scalar_density =
    HadronScalarDensity(mass,
                        proton_fermi_momentum,
                        parameters.hadron.model.cutoff);

    double neutron_scalar_density =
    HadronScalarDensity(mass,
                        neutron_fermi_momentum,
                        parameters.hadron.model.cutoff);

    double proton_chem_pot_gap =
    ProtonChemicalPotentialEquation(proton_fermi_momentum,
                                    proton_scalar_density,
                                    mass,
                                    proton_density,
                                    neutron_density)
    - p->proton_chemical_potential;

    double neutron_chem_pot_gap =
    NeutronChemicalPotentialEquation(neutron_fermi_momentum,
                                    neutron_scalar_density,
                                    mass,
                                    proton_density,
                                    neutron_density)
    - p->neutron_chemical_potential;

    gsl_vector_set(return_values, 0, gap_equation);
    gsl_vector_set(return_values, 1, proton_chem_pot_gap);
    gsl_vector_set(return_values, 2, neutron_chem_pot_gap);

    return GSL_SUCCESS;
}

int GridRootFinder(double proton_chemical_potential,
                   double neutron_chemical_potential,
                   double * return_hadron_mass,
                   double * return_proton_density,
                   double * return_neutron_density,
                   int * return_num_solutions)
{
    GridRootFinderParameters par = parameters.hadron.grid_root_finder;

    hadron_mass_and_renorm_chem_pot_input_params p;
    p.proton_chemical_potential = proton_chemical_potential;
    p.neutron_chemical_potential = neutron_chemical_potential;

    const double hadron_mass_step = Step(par.min_mass,
                                         par.max_mass,
                                         par.num_pts_mass);
    const double density_step = Step(par.min_density,
                                     par.max_density,
                                     par.num_pts_dens);

    const int dimension = 3;

    gsl_vector * input_vec = gsl_vector_alloc(dimension);
    gsl_vector * output_vec = gsl_vector_alloc(dimension);

    const int num_pts = par.num_pts_dens * par.num_pts_dens * par.num_pts_mass;

    int solutions = 0;
    gsl_vector * sols[dimension];
    for (int index = 0; index < dimension; index++)
        sols[index] = gsl_vector_alloc(num_pts);

    double hadron_mass = par.min_mass;
    for(int i = 0; i < par.num_pts_mass; i++){

        double proton_density = par.min_density;
        for(int j = 0; j < par.num_pts_dens; j++){

            double neutron_density = par.min_density;
            for (int k = 0; k < par.num_pts_dens; k++){

                gsl_vector_set(input_vec,
                               0,
                               sqrt(hadron_mass));

                gsl_vector_set(input_vec,
                               1,
                               sqrt(proton_density));

                gsl_vector_set(input_vec,
                               2,
                               sqrt(neutron_density));

                HadronMassAndDensitiesSolutionEquation(input_vec,
                                                       (void *)&p,
                                                       output_vec);

                bool is_solution = true;
                for (int index = 0; index < dimension; index++){

                    double val = gsl_vector_get(output_vec, index);

                    if (fabs(val) > par.zero_tol){

                        is_solution = false;
                        break;
                    }
                }

                if (is_solution){
                    gsl_vector_set(sols[0], solutions, hadron_mass);
                    gsl_vector_set(sols[1], solutions, proton_density);
                    gsl_vector_set(sols[2], solutions, neutron_density);
                    solutions++;
                }

                neutron_density += density_step;
            }

            proton_density += density_step;
        }

        hadron_mass += hadron_mass_step;
    }

    if (solutions > 0){

        double sums[3] = {0,0,0};
        for (int index = 0; index < dimension; index++){

            sums[index] = 0;
            for (int n = 0; n < solutions; n++)
                sums[index] += gsl_vector_get(sols[index], n);

            gsl_vector_free(sols[index]);
        }

        double means[3];
        for (int index = 0; index < dimension; index++)
            means[index] = sums[index] / (double)solutions;

        *return_hadron_mass = means[0];
        *return_proton_density = means[1];
        *return_neutron_density = means[2];
        *return_num_solutions = solutions;

        gsl_vector_free(input_vec);
        gsl_vector_free(output_vec);

        return 0;
    }

    // no solutions found
    gsl_vector_free(input_vec);
    gsl_vector_free(output_vec);

    return -1;
}

double ProtonChemicalPotential(double barionic_chemical_potential,
                               double isovector_chemical_potential)
{
    return barionic_chemical_potential + 0.5 * isovector_chemical_potential;
}

double NeutronChemicalPotential(double barionic_chemical_potential,
                                double isovector_chemical_potential)
{
    return barionic_chemical_potential - 0.5 * isovector_chemical_potential;
}
