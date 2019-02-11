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
    double proton_barionic_density = NAN;
    double neutron_barionic_density = NAN;

    int status =
    HadronMassAndDensitiesSolution(proton_chemical_potential,
                                   neutron_chemical_potential,
                                   hadron_mass_guess,
                                   proton_density_guess,
                                   neutron_density_guess,
                                   &mass,
                                   &proton_barionic_density,
                                   &neutron_barionic_density);

    if (status)
        return -1;

    *return_hadron_mass = mass;
    *return_proton_density = proton_barionic_density;
    *return_neutron_density = neutron_barionic_density;

//    double proton_fermi_momentum = HadronFermiMomentum(proton_density);
//    double neutron_fermi_momentum = HadronFermiMomentum(neutron_density);

    // TODO: Maybe we will need a case for zero temperature here? or
    //       carry the results that are obtained inside the functions and not
    //       recalculate the final results here
    // TODO: I believe the chemical potentials below should be renorm
    double kinectic_energy_density =
    HadronKinecticEnergyDensity(parameters.variables.temperature,
                                mass,
                                proton_barionic_density,
                                neutron_barionic_density,
                                proton_chemical_potential,
                                neutron_chemical_potential);

    // TODO: I believe the chemical potentials below should be renorm
    double proton_scalar_density =
    HadronScalarDensity(mass,
                        proton_chemical_potential,
                        parameters.variables.temperature,
                        proton_barionic_density);
    double neutron_scalar_density =
    HadronScalarDensity(mass,
                        neutron_chemical_potential,
                        parameters.variables.temperature,
                        neutron_barionic_density);

    // TODO: Shouldn't the chemical potentials below be renorm?
    double potential =
    HadronThermodynamicPotential(mass,
                                 proton_scalar_density,
                                 neutron_scalar_density,
                                 proton_barionic_density,
                                 neutron_barionic_density,
                                 proton_chemical_potential,
                                 neutron_chemical_potential,
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
