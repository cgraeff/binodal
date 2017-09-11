//
//  HadronPhaseEOS.c
//  binodal
//
//  Created by Clebson Graeff on 2016-02-17.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>

#include "libdatafun/libdatafun.h"

#include "Parameters.h"
#include "Constants.h"
#include "Functions.h"

typedef struct _hadron_gap_eq_input_params{
    double proton_fermi_momentum;
    double neutron_fermi_momentum;
    double proton_density;
    double neutron_density;
    double renormalized_chemical_potential; // TODO: is this necessary?
} hadron_gap_eq_input_params;

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
typedef struct _hadron_mass_and_renorm_chem_pot_input_params{

    double proton_chemical_potential;
    double neutron_chemical_potential;

} hadron_mass_and_renorm_chem_pot_input_params;

int HadronMassAndDensitiesSolutionEquation(const gsl_vector   *x,
                                           void *params,
                                           gsl_vector *return_values);
int
HadronMassAndDensitiesSolutionEquationZeroMassCase(const gsl_vector   *x,
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
    HadronMassAndDensitiesSolutionParams params =
    parameters.hadron.mass_and_densities_solution;

    hadron_mass_and_renorm_chem_pot_input_params p;
    p.proton_chemical_potential = proton_chemical_potential;
    p.neutron_chemical_potential = neutron_chemical_potential;

    // Handle mass != 0 first as it is the least common case
    // and the one we start with. If no progress is made, we
    // probably reached chiral restoration (mass tends to zero),
    // so we just try to solve with a mass = 0 zero guess. For
    // mass values below ZERO_MASS_TOL, we don't even try the
    // non zero mass solution
    if (hadron_mass_guess > params.zero_mass_tolerance){

        // Set dimension (number of equations|variables to solve|find)
        const int dimension = 3;

        gsl_multiroot_function f;
        f.f = &HadronMassAndDensitiesSolutionEquation;
        f.n = dimension;
        f.params = (void *)&p;

        gsl_vector * initial_guess = gsl_vector_alloc(dimension);
        gsl_vector * return_results = gsl_vector_alloc(dimension);

        gsl_vector_set(initial_guess,
                       0,
                       sqrt(hadron_mass_guess));

        gsl_vector_set(initial_guess,
                       1,
                       sqrt(proton_density_guess));
        gsl_vector_set(initial_guess,
                       2,
                       sqrt(neutron_density_guess));

        int status =
            MultidimensionalRootFinder(dimension,
                                       &f,
                                       initial_guess,
                                       params.abs_error,
                                       params.rel_error,
                                       params.max_iter,
                                       return_results);

        if (status == 0){

            // Save results in return variables,
            // taking care of the mappinps
            *return_mass = pow(gsl_vector_get(return_results, 0), 2.0);
            *return_proton_density = pow(gsl_vector_get(return_results, 1), 2.0);
            *return_neutron_density = pow(gsl_vector_get(return_results, 2), 2.0);

            // Free vectors
            gsl_vector_free(initial_guess);
            gsl_vector_free(return_results);

            return status;

        }
    }

    // Set dimension (number of equations|variables to solve|find)
    const int dimension = 2;

    gsl_multiroot_function f;
    f.f = &HadronMassAndDensitiesSolutionEquationZeroMassCase;
    f.n = dimension;
    f.params = (void *)&p;

    gsl_vector * initial_guess = gsl_vector_alloc(dimension);
    gsl_vector * return_results = gsl_vector_alloc(dimension);

    gsl_vector_set(initial_guess, 0, sqrt(proton_density_guess));
    gsl_vector_set(initial_guess, 1, sqrt(neutron_density_guess));

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
    *return_mass = 0.0;
    *return_proton_density = pow(gsl_vector_get(return_results, 0), 2.0);
    *return_neutron_density = pow(gsl_vector_get(return_results, 1), 2.0);

    // Free vectors
    gsl_vector_free(initial_guess);
    gsl_vector_free(return_results);

    return status;
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
