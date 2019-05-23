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
#include "DefiniteIntegrals.h"
#include "FermiDiracDistributions.h"

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

double HadronProtonFraction(double proton_barionic_density,
                            double neutron_barionic_density)
{
    return proton_barionic_density
           / (proton_barionic_density + neutron_barionic_density);
}

double HadronZeroTemperatureChemicalPotential(double mass,
                                              double fermi_momentum)
{
    return sqrt(pow(mass, 2.0) + pow(fermi_momentum, 2.0));
}

double HadronFermiMomentumFromBarionicDensity(double density)
{
    return CONST_HBAR_C * pow(3.0 * pow(M_PI, 2.0) * density, 1.0 / 3.0);
}

double HadronBarionicDensityFromFermiMomentum(double fermi_momentum)
{
    return pow(CONST_HBAR_C * fermi_momentum, 3.0) / (3.0 * pow(M_PI, 2.0));
}

double HadronPhaseAsymmetry(double proton_density, double neutron_density)
{
    double barionic_density = proton_density + neutron_density;
    double proton_fraction = proton_density / barionic_density;

    return 1.0 - 2.0 * proton_fraction;
}

// TODO: is this function used somewhere?
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

double HadronBarionicDensity(double mass,
                             double renorm_chem_pot,
                             double temperature,
                             double density,
                             double cutoff)
{
    double integral = 0.0;

    if (temperature == 0){

        double fermi_momentum = HadronFermiMomentumFromBarionicDensity(density);
        integral = pow(fermi_momentum, 3.0) / 3.0;

    }
    else{
        integral =
        FermiDiracDistributionIntegralFromBarionicDensity(temperature,
                                                          mass,
                                                          renorm_chem_pot,
                                                          cutoff);
    }

    return NUM_H_COLORS * NUM_H_FLAVORS * pow(CONST_HBAR_C, -3.0)
           * (integral / pow(M_PI, 2.0));
}

double HadronZeroTemperatureScalarDensity(double mass,
                                          double fermi_momentum,
                                          double cutoff)
{
    return pow(CONST_HBAR_C, -3.0)
           * (mass / pow(M_PI, 2.0))
           * (F0(mass, fermi_momentum) - F0(mass, cutoff));
}

double HadronScalarDensity(double mass,
                           double renorm_chem_pot,
                           double temperature,
                           double density)
{
    if (mass == 0.0){
        return 0.0;
    }

    double cutoff = parameters.hadron.model.cutoff;

    double integral = 0.0;
    if (temperature == 0.0){

        double fermi_momentum = HadronFermiMomentumFromBarionicDensity(density);

        integral =
        F0(mass, cutoff) - F0(mass, fermi_momentum);

    }
    else {

        integral =
        FermiDiracDistributionIntegralFromScalarDensity(temperature,
                                                        mass,
                                                        renorm_chem_pot,
                                                        cutoff);
    }

    return - NUM_H_COLORS * NUM_H_FLAVORS * pow(CONST_HBAR_C, -3.0)
           * mass * integral / pow(M_PI, 2.0);
}

// TODO: This should be a special case of the full HadronScalarDensity()
//       function
double HadronVacuumScalarDensity()
{
    return 2.0 * pow(CONST_HBAR_C, -3.0)
           * (parameters.hadron.model.nucleon_mass / pow(M_PI, 2.0))
           * (F0(parameters.hadron.model.nucleon_mass, 0.0)
              - F0(parameters.hadron.model.nucleon_mass,
                   parameters.hadron.model.cutoff));
}

double ProtonChemicalPotentialEquation(double proton_fermi_momentum,
                                       double proton_scalar_density,
                                       double neutron_scalar_density,
                                       double mass,
                                       double proton_density,
                                       double neutron_density)
{
    double E = sqrt(pow(mass, 2.0) + pow(proton_fermi_momentum, 2.0));

    double barionic_density = proton_density + neutron_density;
    double rho_3 = (proton_density - neutron_density);
    double scalar_density = proton_scalar_density + neutron_scalar_density;

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
                                        double proton_scalar_density,
                                        double neutron_scalar_density,
                                        double mass,
                                        double proton_density,
                                        double neutron_density)
{
    double E = sqrt(pow(mass, 2.0) + pow(neutron_fermi_momentum, 2.0));

    double barionic_density = proton_density + neutron_density;
    double rho_3 = (proton_density - neutron_density);
    double scalar_density = proton_scalar_density + neutron_scalar_density;

    double rho_terms =
    (parameters.hadron.model.G_V * barionic_density
     + parameters.hadron.model.G_SV * barionic_density
       * pow(scalar_density, 2.0)
     - parameters.hadron.model.G_RHO * rho_3
     + parameters.hadron.model.G_VRHO * pow(rho_3, 2.0)
       * barionic_density
     - parameters.hadron.model.G_VRHO * pow(barionic_density, 2.0)
       * rho_3
     - parameters.hadron.model.G_SRHO * pow(scalar_density, 2.0)
       * rho_3) * 2.0 * CONST_HBAR_C;

    return E + rho_terms;
}

double ProtonRenormalizedChemicalPotential(double proton_chemical_potential,
                                           double scalar_density,
                                           double mass,
                                           double proton_density,
                                           double neutron_density)
{
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

    return proton_chemical_potential - rho_terms;
}

double NeutronRenormalizedChemicalPotential(double neutron_chemical_potential,
                                            double scalar_density,
                                            double mass,
                                            double proton_density,
                                            double neutron_density)
{
    double barionic_density = proton_density + neutron_density;
    double rho_3 = (proton_density - neutron_density);

    double rho_terms =
    (parameters.hadron.model.G_V * barionic_density
     + parameters.hadron.model.G_SV * barionic_density
       * pow(scalar_density, 2.0)
     - parameters.hadron.model.G_RHO * rho_3
     + parameters.hadron.model.G_VRHO * pow(rho_3, 2.0)
       * barionic_density
     - parameters.hadron.model.G_VRHO * pow(barionic_density, 2.0)
       * rho_3
     - parameters.hadron.model.G_SRHO * pow(scalar_density, 2.0)
       * rho_3) * 2.0 * CONST_HBAR_C;

    return neutron_chemical_potential - rho_terms;
}

// All inputs in MeV, output in MeV/fm^3
double HadronKinecticEnergyDensity(double temperature,
                                   double mass,
                                   double proton_density,
                                   double neutron_density,
                                   double proton_renorm_chem_pot,
                                   double neutron_renorm_chem_pot)
{

    double cutoff = parameters.hadron.model.cutoff;

    double integral = 0.0;

    if (temperature == 0.0){

        double proton_fermi_momentum =
        HadronFermiMomentumFromBarionicDensity(proton_density);

        double neutron_fermi_momentum =
        HadronFermiMomentumFromBarionicDensity(neutron_density);

        double proton_kinectic_energy =
        F2(mass, cutoff) - F2(mass, proton_fermi_momentum);

        double neutron_kinectic_energy =
        F2(mass, cutoff) - F2(mass, neutron_fermi_momentum);

        integral = proton_kinectic_energy + neutron_kinectic_energy;
    }
    else{
        integral =
        FermiDiracDistributionIntegralFromHadronEnergy(temperature,
                                                       mass,
                                                       proton_renorm_chem_pot,
                                                       cutoff)
        + FermiDiracDistributionIntegralFromHadronEnergy(temperature,
                                                         mass,
                                                         neutron_renorm_chem_pot,
                                                         cutoff);
    }

    return - NUM_H_FLAVORS * NUM_H_COLORS * pow(CONST_HBAR_C, -3.0)
           * integral / pow(M_PI, 2.0);
}

// TODO: Both functions below should not exist, they should just be special
//       cases of HadronKinecticEnergyDensity()
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

/*
double HadronEntropyIntegrand(double  momentum,
                              void   *par)
{
    fermi_dirac_distrib_integrand * p = (fermi_dirac_distrib_integrand *)par;

    double E = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));

    double n_p = FermiDiracDistributionForParticles(E,
                                                    p->chemical_potential,
                                                    p->temperature);

    double n_ap = FermiDiracDistributionForAntiparticles(E,
                                                         p->chemical_potential,
                                                         p->temperature);

    return pow(momentum, 2.0) * (n_p * log(n_p) + (1.0 - n_p) * log(1.0 - n_p)
                          + n_ap * log(n_ap) + (1.0 - n_ap) * log(1.0 - n_ap));
}
 */

double HadronEntropy(double mass,
                     double renorm_chem_pot,
                     double temperature)
{
    /*
    fermi_dirac_distrib_integrand p;
    p.mass = mass;
    p.chemical_potential = renorm_chem_pot;
    p.temperature = temperature;

    gsl_function F;
    F.function = &HadronEntropyIntegrand;
    F.params = &p;

    double integral =
    OnedimensionalIntegrator(&F, parameters.fermi_dirac_integrals);
     */

    double integral =
    FermiDiracDistributionIntegralFromHadronEntropy(temperature,
                                                    mass,
                                                    renorm_chem_pot,
                                                    parameters.hadron.model.cutoff);

    return -2.0 * NUM_H_FLAVORS * NUM_H_COLORS * pow(CONST_HBAR_C, -3.0)
           * integral / pow(M_PI, 2.0);
}

double HadronThermodynamicPotential(double mass,
                                    double proton_scalar_density,
                                    double neutron_scalar_density,
                                    double proton_barionic_density,
                                    double neutron_barionic_density,
                                    double proton_chemical_potential,
                                    double neutron_chemical_potential,
                                    double proton_renorm_chem_pot,
                                    double neutron_renorm_chem_pot,
                                    double kinectic_energy_density)
{
    double barionic_density = proton_barionic_density
                              + neutron_barionic_density;
    double rho_3 = proton_barionic_density - neutron_barionic_density;

    double scalar_density = proton_scalar_density + neutron_scalar_density;

    double omega = kinectic_energy_density
                   + parameters.hadron.model.bare_mass * scalar_density;

    omega += - proton_chemical_potential * proton_barionic_density
             - neutron_chemical_potential * neutron_barionic_density;

    omega += (- parameters.hadron.model.G_S * pow(scalar_density, 2.0)
              + parameters.hadron.model.G_V * pow(barionic_density, 2.0)
              + parameters.hadron.model.G_SV
                * pow(scalar_density * barionic_density, 2.0)
              + parameters.hadron.model.G_RHO * pow(rho_3, 2.0)
              + parameters.hadron.model.G_VRHO
                * pow(barionic_density * rho_3, 2.0)
              + parameters.hadron.model.G_SRHO
                * pow(scalar_density * rho_3, 2.0))
             * CONST_HBAR_C;

    if (parameters.variables.temperature != 0){

        double proton_contrib =
        HadronEntropy(mass,
                      proton_renorm_chem_pot,
                      parameters.variables.temperature);

        double neutron_contrib =
        HadronEntropy(mass,
                      neutron_renorm_chem_pot,
                      parameters.variables.temperature);

        omega -= parameters.variables.temperature
                 * (proton_contrib + neutron_contrib);
    }

    return omega;
}

double HadronZeroedGapEquation(double mass,
                               double proton_density,
                               double neutron_density,
                               double proton_scalar_density,
                               double neutron_scalar_density)
{
    double barionic_density = proton_density + neutron_density;
    double rho_3 = proton_density - neutron_density;
    double scalar_density = proton_scalar_density + neutron_scalar_density;

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
    // TODO: TEST exp(X) instead of x^2

    const double mass = pow(gsl_vector_get(x, 0), 2.0);
    const double proton_density = pow(gsl_vector_get(x, 1), 2.0);
    const double neutron_density = pow (gsl_vector_get(x, 2), 2.0);

    // TODO: maybe there is no need to make two cases
    if (parameters.variables.temperature == 0){
        hadron_mass_and_renorm_chem_pot_input_params * p =
        (hadron_mass_and_renorm_chem_pot_input_params *) params;

        double proton_fermi_momentum =
        HadronFermiMomentumFromBarionicDensity(proton_density);

        double neutron_fermi_momentum =
        HadronFermiMomentumFromBarionicDensity(neutron_density);

        //calculate scalar densities
        double proton_scalar_density =
        HadronZeroTemperatureScalarDensity(mass,
                                           proton_fermi_momentum,
                                           parameters.hadron.model.cutoff);

        double neutron_scalar_density =
        HadronZeroTemperatureScalarDensity(mass,
                                           neutron_fermi_momentum,
                                           parameters.hadron.model.cutoff);

        double gap_equation =
        HadronZeroedGapEquation(mass,
                                proton_density,
                                neutron_density,
                                proton_scalar_density,
                                neutron_scalar_density);


        double proton_chem_pot_gap =
        ProtonChemicalPotentialEquation(proton_fermi_momentum,
                                        proton_scalar_density,
                                        neutron_scalar_density,
                                        mass,
                                        proton_density,
                                        neutron_density)
        - p->proton_chemical_potential;

        double neutron_chem_pot_gap =
        NeutronChemicalPotentialEquation(neutron_fermi_momentum,
                                         proton_scalar_density,
                                         neutron_scalar_density,
                                         mass,
                                         proton_density,
                                         neutron_density)
        - p->neutron_chemical_potential;

        gsl_vector_set(return_values, 0, gap_equation);
        gsl_vector_set(return_values, 1, proton_chem_pot_gap);
        gsl_vector_set(return_values, 2, neutron_chem_pot_gap);
    }
    else{

        hadron_mass_and_renorm_chem_pot_input_params * p =
        (hadron_mass_and_renorm_chem_pot_input_params *) params;

        double proton_scalar_density =
        HadronScalarDensity(mass,
                            p->proton_chemical_potential,
                            parameters.variables.temperature,
                            proton_density);

        double neutron_scalar_density =
        HadronScalarDensity(mass,
                            p->neutron_chemical_potential,
                            parameters.variables.temperature,
                            neutron_density);

        double proton_renorm_chem_pot =
        ProtonRenormalizedChemicalPotential(p->proton_chemical_potential,
                                            proton_scalar_density,
                                            mass,
                                            proton_density,
                                            neutron_density);

        double neutron_renorm_chem_pot =
        NeutronRenormalizedChemicalPotential(p->neutron_chemical_potential,
                                             neutron_scalar_density,
                                             mass,
                                             proton_density,
                                             neutron_density);

        double gap_equation =
        HadronZeroedGapEquation(mass,
                                proton_density,
                                neutron_density,
                                proton_scalar_density,
                                neutron_scalar_density);

        double proton_dens_eq =
        HadronBarionicDensity(mass,
                              proton_renorm_chem_pot,
                              parameters.variables.temperature,
                              proton_density,
                              parameters.hadron.model.cutoff)
        - proton_density;

        double neutron_dens_eq =
        HadronBarionicDensity(mass,
                              neutron_renorm_chem_pot,
                              parameters.variables.temperature,
                              neutron_density,
                              parameters.hadron.model.cutoff)
        - neutron_density;

        gsl_vector_set(return_values, 0, gap_equation);
        gsl_vector_set(return_values, 1, proton_dens_eq);
        gsl_vector_set(return_values, 2, neutron_dens_eq);
    }

    return GSL_SUCCESS;
}

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
    IsovectorChemicalPotential(proton_chemical_potential,
                               neutron_chemical_potential);

    const gsl_multiroot_fsolver_type * solver_type =
    gsl_multiroot_fsolver_hybrids;

    HadronMassAndDensitiesSolutionParams params =
    parameters.hadron.mass_and_densities_solution;

    MultidimensionalRootFinderParams rootf_pars;
    rootf_pars.solver_type = solver_type;
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
        && proton_density_guess < params.zero_dens_tolerance
        && neutron_density_guess < params.zero_dens_tolerance){

        *return_mass = 0.0;
        *return_proton_density = 0.0;
        *return_neutron_density = 0.0;

        return 0;
    }

    // TODO: Is there a way to shorten this function using the proper
    //       "f.f = &Hadron ..."? This is too long and look quite
    //       repetitive.

    // Case where only mass goes to zero
    if (hadron_mass_guess < params.zero_mass_tolerance
        && proton_density_guess > params.zero_dens_tolerance
        && neutron_density_guess > params.zero_dens_tolerance){

        // Set dimension (number of equations|variables to solve|find)
        const int dimension = 2;

        rootf_pars.solver_type = solver_type;

        gsl_multiroot_function f;
        f.f = &HadronMassAndDensitiesSolutionEquationZeroMassCase;
        f.n = dimension;
        f.params = (void *)&p;

        gsl_vector * initial_guess = gsl_vector_alloc(dimension);
        gsl_vector * return_results = gsl_vector_alloc(dimension);

        // TODO: TEST exp(x) instead of x^2

        gsl_vector_set(initial_guess, 0, sqrt(proton_density_guess));
        gsl_vector_set(initial_guess, 1, sqrt(neutron_density_guess));

        int status =
        MultidimensionalRootFinder(&f,
                                   &rootf_pars,
                                   dimension,
                                   initial_guess,
                                   return_results);

        if (status != 0){

            if (options.debug)
                printf("\n[%s:%d: No solution for:\n"
                           "\tbar chem pot: %f\n"
                           "\tiso chem pot: %f\n"
                           "\thadron mass guess: %f\n"
                           "\tproton density guess: %f\n"
                           "\tneutron density guess: %f]\n",
                           __FILE__,
                           __LINE__,
                           barionic_chemical_potential,
                           isovector_chemical_potential,
                           hadron_mass_guess,
                           proton_density_guess,
                           neutron_density_guess);

            // Free vectors
            gsl_vector_free(initial_guess);
            gsl_vector_free(return_results);

            if (options.abort_on_error){

                printf("%s:%d:"
                       "No solution found for hadron mass and densities.\n",
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
        && proton_density_guess < params.zero_dens_tolerance
        && neutron_density_guess > params.zero_dens_tolerance){

        // Set dimension (number of equations|variables to solve|find)
        const int dimension = 2;

        rootf_pars.solver_type = solver_type;

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

           if (options.debug)
                printf("\n[%s:%d: No solution for:\n"
                           "\tbar chem pot: %f\n"
                           "\tiso chem pot: %f\n"
                           "\thadron mass guess: %f\n"
                           "\tproton density guess: %f\n"
                           "\tneutron density guess: %f]\n",
                           __FILE__,
                           __LINE__,
                           barionic_chemical_potential,
                           isovector_chemical_potential,
                           hadron_mass_guess,
                           proton_density_guess,
                           neutron_density_guess);

            // Free vectors
            gsl_vector_free(initial_guess);
            gsl_vector_free(return_results);

            if (options.abort_on_error){

                printf("%s:%d:"
                       "No solution found for hadron mass and densities.\n",
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
        && proton_density_guess > params.zero_dens_tolerance
        && neutron_density_guess < params.zero_dens_tolerance){

        // Set dimension (number of equations|variables to solve|find)
        const int dimension = 2;

        rootf_pars.solver_type = solver_type;

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

           if (options.debug)
                printf("\n[%s:%d: No solution for:\n"
                           "\tbar chem pot: %f\n"
                           "\tiso chem pot: %f\n"
                           "\thadron mass guess: %f\n"
                           "\tproton density guess: %f\n"
                           "\tneutron density guess: %f]\n",
                           __FILE__,
                           __LINE__,
                           barionic_chemical_potential,
                           isovector_chemical_potential,
                           hadron_mass_guess,
                           proton_density_guess,
                           neutron_density_guess);

            // Free vectors
            gsl_vector_free(initial_guess);
            gsl_vector_free(return_results);

            if (options.abort_on_error){

                printf("%s:%d:"
                       "No solution found for hadron mass and densities.\n",
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

    // Standard case, where no variables are near zero
    int status = -1;
    if (hadron_mass_guess > params.zero_mass_tolerance
        && proton_density_guess > params.zero_dens_tolerance
        && neutron_density_guess > params.zero_dens_tolerance){

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

            *return_proton_density =
            pow(gsl_vector_get(return_results, 1), 2.0);

            *return_neutron_density =
            pow(gsl_vector_get(return_results, 2), 2.0);

            // Free vectors
            gsl_vector_free(initial_guess);
            gsl_vector_free(return_results);

            return 0;

        }

        if (options.debug)
            printf("\n[%s:%d: No solution for:\n"
                   "\tbar chem pot: %f\n"
                   "\tiso chem pot: %f\n"
                   "\thadron mass guess: %f\n"
                   "\tproton density guess: %f\n"
                   "\tneutron density guess: %f]\n",
                   __FILE__,
                   __LINE__,
                   barionic_chemical_potential,
                   isovector_chemical_potential,
                   hadron_mass_guess,
                   proton_density_guess,
                   neutron_density_guess);

        // Free vectors
        gsl_vector_free(initial_guess);
        gsl_vector_free(return_results);

        if (options.abort_on_error){

            printf("%s:%d:"
                   "No solution found for hadron mass and densities.\n",
                   __FILE__, __LINE__);

            abort();
        }

        return -1;
    }

    if (hadron_mass_guess > params.zero_mass_tolerance
        && proton_density_guess < params.zero_dens_tolerance
        && neutron_density_guess < params.zero_dens_tolerance){
            printf("both dens < tol\n");
            printf("%s:%d\n", __FILE__, __LINE__);
            abort();
    }

    if (hadron_mass_guess < params.zero_mass_tolerance
        && proton_density_guess < params.zero_dens_tolerance
        && neutron_density_guess > params.zero_dens_tolerance){
            printf("m, p < tol\n");
            printf("%s:%d\n", __FILE__, __LINE__);
            abort();
    }

    if (hadron_mass_guess < params.zero_mass_tolerance
        && proton_density_guess > params.zero_dens_tolerance
        && neutron_density_guess > params.zero_dens_tolerance){
            printf("m, n < tol\n");
            printf("%s:%d\n", __FILE__, __LINE__);
            abort();
    }

    // This return point should not be reached
    printf("%s:%d ", __FILE__, __LINE__);
    printf("Problems in management of possible calculation scenarios.\n");

    abort();

    return -1;
}

