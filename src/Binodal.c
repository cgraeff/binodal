//
//  Binodal.c
//  binodal
//
//  Created by Clebson Graeff on 2017-06-20.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#include <assert.h>
#include <math.h>

#include <gsl/gsl_vector.h>

#include "libdatafun/libdatafun.h"

#include "Constants.h"
#include "CommandlineOptions.h"
#include "QuarkPhaseEOS.h"
#include "HadronPhaseEOS.h"
#include "Parameters.h"

#include "Binodal.h"

typedef struct _binodal_parameters {
    double isovector_chemical_potential;
    double temperature;

    double hadron_vacuum_thermodynamic_potential;
    double quark_vacuum_thermodynamic_potential;

    double hadron_mass_guess;
    double proton_density_guess;
    double neutron_density_guess;

    double up_mass_guess;
    double down_mass_guess;

} binodal_parameters;


int BinodalPointCandidate(double barionic_chemical_potential,
                          double isovector_chemical_potential,
                          double temperature,
                          double hadron_vacuum_thermodynamic_potential,
                          double hadron_mass_guess,
                          double proton_density_guess,
                          double neutron_density_guess,
                          double quark_vacuum_thermodynamic_potential,
                          double up_mass_guess,
                          double down_mass_guess,
                          BinodalPoint *return_point){

    double proton_chemical_potential =
    ProtonChemicalPotential(barionic_chemical_potential,
                            isovector_chemical_potential);

    double neutron_chemical_potential =
    NeutronChemicalPotential(barionic_chemical_potential,
                             isovector_chemical_potential);

    double hadron_mass;
    double hadron_pressure;
    double proton_density;
    double neutron_density;

    int status_h =
    DetermineHadronPressureAndDensities(proton_chemical_potential,
                                        neutron_chemical_potential,
                                        hadron_vacuum_thermodynamic_potential,
                                        hadron_mass_guess,
                                        proton_density_guess,
                                        neutron_density_guess,
                                        &hadron_mass,
                                        &proton_density,
                                        &neutron_density,
                                        &hadron_pressure);

    if (status_h)
        return -1;

    double up_chemical_potential =
    UpChemicalPotentialFromGibbsConditions(proton_chemical_potential,
                                           neutron_chemical_potential);

    double down_chemical_potential =
    DownChemicalPotentialFromGibbsConditions(proton_chemical_potential,
                                             neutron_chemical_potential);

    double up_quark_mass;
    double down_quark_mass;
    double quark_pressure;

    int status_q =
    DetermineQuarkPressureAndMasses(up_chemical_potential,
                                    down_chemical_potential,
                                    temperature,
                                    quark_vacuum_thermodynamic_potential,
                                    up_mass_guess,
                                    down_mass_guess,
                                    &up_quark_mass,
                                    &down_quark_mass,
                                    &quark_pressure);

    if (status_q)
        return -1;

    double up_renormalized_chemical_potential;
    double down_renormalized_chemical_potential;

    status_q =
    QuarkSelfConsistentRenormChemPot(up_quark_mass,
                                     down_quark_mass,
                                     up_chemical_potential,
                                     down_chemical_potential,
                                     parameters.variables.temperature,
                                     &up_renormalized_chemical_potential,
                                     &down_renormalized_chemical_potential);

    if (status_q)
        return -1;

    double up_quark_density =
    QuarkDensity(up_quark_mass,
                 up_renormalized_chemical_potential,
                 parameters.variables.temperature);

    double down_quark_density =
    QuarkDensity(down_quark_mass,
                 down_renormalized_chemical_potential,
                 parameters.variables.temperature);

    BinodalPoint point;
    point.hadron_mass = hadron_mass;
    point.proton_density = proton_density;
    point.neutron_density = neutron_density;
    point.proton_chemical_potential = proton_chemical_potential;
    point.neutron_chemical_potential = neutron_chemical_potential;
    point.hadron_pressure = hadron_pressure;

    point.up_quark_mass = up_quark_mass;
    point.down_quark_mass = down_quark_mass;
    point.up_quark_density = up_quark_density;
    point.down_quark_density = down_quark_density;
    point.up_chemical_potential = up_chemical_potential;
    point.down_chemical_potential = down_chemical_potential;
    point.quark_pressure = quark_pressure;

    *return_point = point;

    return 0;
}

double BinodalPointEquation(double  barionic_density,
                            void   *params);

// TODO: Test, if this doesn't work, throw away
int
DetermineBinodalPointCandidateStepwise(double temperature,
                                       double isovector_chemical_potential,
                                       double hadron_vacuum_potential,
                                       double quark_vacuum_potential,
                                       double *barionic_chem_pot_lower_bound,
                                       double *barionic_chem_pot_upper_bound,
                                       BinodalPoint * return_point)
{
    double hadron_mass_guess =
    parameters.hadron.mass_and_densities_solution.initial_mass_guess;

    double proton_density_guess =
    parameters.hadron.mass_and_densities_solution.initial_proton_density_guess;

    double neutron_density_guess =
    parameters.hadron.mass_and_densities_solution.initial_neutron_density_guess;

    double up_mass_guess =
    parameters.quark.mass_and_renorm_chem_pot_solution.initial_up_mass_guess;

    double down_mass_guess =
    parameters.quark.mass_and_renorm_chem_pot_solution.initial_down_mass_guess;

    double min_barionic_chemical_potential =
    parameters.binodal_rootfinding_params.lower_bound;

    double max_barionic_chemical_potential =
    parameters.binodal_rootfinding_params.upper_bound;

    double barionic_chemical_potential_step =
    parameters.binodal_rootfinding_params.step_size;

    double barionic_chemical_potential = min_barionic_chemical_potential;
    while (barionic_chemical_potential <= max_barionic_chemical_potential){

        BinodalPoint point;

        int status =
        BinodalPointCandidate(barionic_chemical_potential,
                              isovector_chemical_potential,
                              temperature,
                              hadron_vacuum_potential,
                              hadron_mass_guess,
                              proton_density_guess,
                              neutron_density_guess,
                              quark_vacuum_potential,
                              up_mass_guess,
                              down_mass_guess,
                              &point);

        if (status){

            if (options.abort_on_error){

                printf("%s:%d:"
                       "Abort on error enabled, aborting due to error.\n",
                       __FILE__, __LINE__);

                abort();
            }
            return -1;
        }

        printf("hp: %f\nqp: %f\n", point.hadron_pressure, point.quark_pressure);
        if ((point.hadron_pressure - point.quark_pressure) > 0.0){

            *barionic_chem_pot_lower_bound = barionic_chemical_potential;

            *barionic_chem_pot_upper_bound =
            barionic_chemical_potential + barionic_chemical_potential_step;

            // For masses with values below ZERO_MASS_TOL,
            // just assume zero.
            double hadron_zero_mass_tol =
            parameters.hadron.mass_and_densities_solution.zero_mass_tolerance;

            hadron_mass_guess =
            point.hadron_mass < hadron_zero_mass_tol? 0.0 : point.hadron_mass;
            proton_density_guess = point.proton_density;
            neutron_density_guess = point.neutron_density;

            double quark_zero_mass_tol =
            parameters.quark.mass_and_renorm_chem_pot_solution.zero_mass_tolerance;

            up_mass_guess =
            point.up_quark_mass < quark_zero_mass_tol? 0.0 : point.up_quark_mass;

            down_mass_guess =
            point.down_quark_mass < quark_zero_mass_tol? 0.0 : point.down_quark_mass;

        }
        else{

            *return_point = point;

            return 0;

        }

        barionic_chemical_potential += barionic_chemical_potential_step;

    }

   if (options.abort_on_error){

        printf("%s:%d:"
               "Abort on error enabled, aborting due to error.\n",
               __FILE__, __LINE__);

        abort();
    }

    return -1;
}


// TODO: Test, if it doesn't work, throw away
int
DetermineBinodalPointByBissection(double temperature,
                                  double isovector_chemical_potential,
                                  double hadron_vacuum_thermodynamic_potential,
                                  double quark_vacuum_thermodynamic_potential,
                                  double transition_bar_chem_pot_lower_bound,
                                  double transition_bar_chem_pot_upper_bound,
                                  double hadron_mass_guess,
                                  double proton_density_guess,
                                  double neutron_density_guess,
                                  double up_mass_guess,
                                  double down_mass_guess,
                                  BinodalPoint *return_point)
{

    // Determine which value of density gives equal pressures for each phase:

    binodal_parameters params;
    params.temperature = temperature;
    params.isovector_chemical_potential = isovector_chemical_potential;
    params.hadron_vacuum_thermodynamic_potential =
    hadron_vacuum_thermodynamic_potential;
    params.quark_vacuum_thermodynamic_potential =
    quark_vacuum_thermodynamic_potential;
    params.hadron_mass_guess = hadron_mass_guess;
    params.proton_density_guess = proton_density_guess;
    params.neutron_density_guess = neutron_density_guess;
    params.up_mass_guess = up_mass_guess;
    params.down_mass_guess = down_mass_guess;

    gsl_function F;
    F.function = &BinodalPointEquation;
    F.params = &params;

    double barionic_chemical_potential;

    UnidimensionalRootFindingParameters p;
    p.lower_bound = transition_bar_chem_pot_lower_bound;
    p.upper_bound = transition_bar_chem_pot_upper_bound;
    p.max_iterations = parameters.binodal_rootfinding_params.max_iterations;
    p.abs_error = parameters.binodal_rootfinding_params.abs_error;
    p.rel_error = parameters.binodal_rootfinding_params.rel_error;

    int status = UnidimensionalRootFinder(&F,
                                          p,
                                          &barionic_chemical_potential);

    if (status){

        if (options.abort_on_error){

            printf("%s:%d:"
                   "Abort on error enabled, aborting due to error.\n",
                   __FILE__, __LINE__);

            abort();
        }
        return -1;
    }

    BinodalPoint point;

    int status_c =
    BinodalPointCandidate(barionic_chemical_potential,
                          isovector_chemical_potential,
                          temperature,
                          hadron_vacuum_thermodynamic_potential,
                          hadron_mass_guess,
                          proton_density_guess,
                          neutron_density_guess,
                          quark_vacuum_thermodynamic_potential,
                          up_mass_guess,
                          down_mass_guess,
                          &point);

    if (status_c){

        if (options.abort_on_error){

            printf("%s:%d:"
                   "Abort on error enabled, aborting due to error.\n",
                   __FILE__, __LINE__);

            abort();
        }

        return -1;
    }

    *return_point = point;

    return 0;
}

double BinodalPointEquation(double  barionic_chemical_potential,
                            void   *params)
{

    binodal_parameters * p = (binodal_parameters *)params;

    BinodalPoint point = {0};
    BinodalPointCandidate(barionic_chemical_potential,
                          p->isovector_chemical_potential,
                          p->temperature,
                          p->hadron_vacuum_thermodynamic_potential,
                          p->hadron_mass_guess,
                          p->proton_density_guess,
                          p->neutron_density_guess,
                          p->quark_vacuum_thermodynamic_potential,
                          p->up_mass_guess,
                          p->down_mass_guess,
                          &point);

    return point.hadron_pressure - point.quark_pressure;
}

int DetermineHadronPressureAndDensities(double proton_chemical_potential,
                                        double neutron_chemical_potential,
                                        double hadron_vacuum_potential,
                                        double hadron_mass_guess,
                                        double proton_density_guess,
                                        double neutron_density_guess,
                                        double *return_hadron_mass,
                                        double *return_proton_density,
                                        double *return_neutron_density,
                                        double *return_pressure)
{
    double mass = NAN;
    double proton_density = NAN;
    double neutron_density = NAN;

    int status =
    HadronMassAndDensitiesSolution(proton_chemical_potential,
                                   neutron_chemical_potential,
                                   hadron_mass_guess,
                                   proton_density_guess,
                                   neutron_density_guess,
                                   &mass,
                                   &proton_density,
                                   &neutron_density);

    if (status)
        return -1;

    *return_hadron_mass = mass;
    *return_proton_density = proton_density;
    *return_neutron_density = neutron_density;

    double proton_fermi_momentum = HadronFermiMomentum(proton_density);
    double neutron_fermi_momentum = HadronFermiMomentum(neutron_density);

    double kinectic_energy_density =
    HadronKinecticEnergyDensity(mass,
                                proton_fermi_momentum,
                                neutron_fermi_momentum);

    double scalar_density =
    HadronScalarDensity(mass,
                        proton_fermi_momentum,
                        parameters.hadron.model.cutoff)
    + HadronScalarDensity(mass,
                          neutron_fermi_momentum,
                          parameters.hadron.model.cutoff);

    double barionic_density = proton_density + neutron_density;
    double potential =
    HadronThermodynamicPotential(scalar_density,
                                 barionic_density,
                                 proton_density,
                                 neutron_density,
                                 proton_chemical_potential,
                                 neutron_chemical_potential,
                                 kinectic_energy_density)
    - hadron_vacuum_potential;

    double hadron_pressure = HadronPressure(potential);

    *return_pressure = hadron_pressure;

    return 0;
}

int DetermineQuarkPressureAndMasses(double up_chemical_potential,
                                    double down_chemical_potential,
                                    double temperature,
                                    double quark_vacuum_thermodynamic_potential,
                                    double up_mass_guess,
                                    double down_mass_guess,
                                    double *return_up_mass,
                                    double *return_down_mass,
                                    double *return_pressure)
{
    double up_mass = NAN;
    double down_mass = NAN;
    double up_renorm_chem_pot = NAN;
    double down_renorm_chem_pot = NAN;

    QuarkMassAndRenormChemPotSolutionBissection(up_chemical_potential,
                                                down_chemical_potential,
                                                temperature,
                                                &up_mass,
                                                &down_mass,
                                                &up_renorm_chem_pot,
                                                &down_renorm_chem_pot);

    double omega =
    QuarkThermodynamicPotential(up_mass,
                                down_mass,
                                up_chemical_potential,
                                down_chemical_potential,
                                up_renorm_chem_pot,
                                down_renorm_chem_pot,
                                temperature);

    double quark_pressure =
    QuarkPressure(omega - quark_vacuum_thermodynamic_potential,
                  temperature);

    *return_up_mass = up_mass;
    *return_down_mass = down_mass;
    *return_pressure = quark_pressure;

    return 0;
}

double BarionicChemicalPotential(double proton_chemical_potential,
                                 double neutron_chemical_potential)
{
    double barionic_chemical_potential =
    (proton_chemical_potential + neutron_chemical_potential) / 2.0;

    return barionic_chemical_potential;
}

double IsovectorChemicalPotential(double proton_chemical_potential,
                                  double neutron_chemical_potential)
{
    double isovector_chemical_potential =
    (proton_chemical_potential - neutron_chemical_potential);

    return isovector_chemical_potential;
}

double UpChemicalPotentialFromGibbsConditions(double proton_chemical_potential,
                                              double neutron_chemical_potential)
{
    double up_chemical_potential =
    (2.0 * proton_chemical_potential - neutron_chemical_potential) / 3.0;

    return up_chemical_potential;
}

double DownChemicalPotentialFromGibbsConditions(double proton_chemical_potential,
                                                double neutron_chemical_potential)
{
    double down_chemical_potential =
    (-proton_chemical_potential + 2.0 * neutron_chemical_potential) / 3.0;

    return down_chemical_potential;
}
