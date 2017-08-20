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
#include "QuarkPhaseEOS.h"
#include "HadronPhaseEOS.h"
#include "Parameters.h"

#include "Binodal.h"

typedef struct _binodal_parameters {
    double isovector_chemical_potential;
    double temperature;

    double hadron_vacuum_thermodynamic_potential;
    double quark_vacuum_thermodynamic_potential;
} binodal_parameters;

double BinodalPointEquation(double  barionic_density,
                            void   *params);

BinodalPoint DetermineBinodalPoint(double temperature,
                                   double isovector_chemical_potential,
                                   double hadron_vacuum_thermodynamic_potential,
                                   double quark_vacuum_thermodynamic_potential)
{

    // Determine which value of density gives equal pressures for each phase:

    binodal_parameters params;
    params.temperature = temperature;
    params.isovector_chemical_potential = isovector_chemical_potential;
    params.hadron_vacuum_thermodynamic_potential = hadron_vacuum_thermodynamic_potential;
    params.quark_vacuum_thermodynamic_potential = quark_vacuum_thermodynamic_potential;

    gsl_function F;
    F.function = &BinodalPointEquation;
    F.params = &params;

    double barionic_chemical_potential;

    int status = UnidimensionalRootFinder(&F,
                                          parameters.binodal_rootfinding_params,
                                          &barionic_chemical_potential);

    if (status != 0){
        printf("%s:%d: Problem with rootfinding.\n", __FILE__, __LINE__);
        abort();
    }

    // Determine results to return:
    //      Note that the following calculation is the same as the one performed
    //      in BinodalPointEquation() in the last iteration of the root finding.
    //      We recalculate just to get the results that we are interested in,
    //      but in a clean way. Saving the results at BinodalPointEquation and
    //      returning would require updating a series of variables inside the
    //      parameters binodal_parameters struct.

    double proton_chemical_potential =
    ProtonChemicalPotential(barionic_chemical_potential,
                            isovector_chemical_potential);

    double neutron_chemical_potential =
    NeutronChemicalPotential(barionic_chemical_potential,
                             isovector_chemical_potential);

    double hadron_mass = NAN;
    double hadron_pressure = NAN;
    double proton_density = NAN;
    double neutron_density = NAN;

    DetermineHadronPressureAndDensities(proton_chemical_potential,
                                        neutron_chemical_potential,
                                        hadron_vacuum_thermodynamic_potential,
                                        &hadron_mass,
                                        &proton_density,
                                        &neutron_density,
                                        &hadron_pressure);

    double up_chemical_potential =
    UpChemicalPotentialFromGibbsConditions(proton_chemical_potential,
                                           neutron_chemical_potential);

    double down_chemical_potential =
    DownChemicalPotentialFromGibbsConditions(proton_chemical_potential,
                                             neutron_chemical_potential);

    double up_mass = NAN;
    double down_mass = NAN;
    double quark_pressure = NAN;
    DetermineQuarkPressure(up_chemical_potential,
                           down_chemical_potential,
                           parameters.variables.temperature,
                           quark_vacuum_thermodynamic_potential,
                           &up_mass,
                           &down_mass,
                           &quark_pressure);

    double up_renormalized_chemical_potential = NAN;
    double down_renormalized_chemical_potential = NAN;
    QuarkSelfConsistentRenormChemPot(up_mass,
                                     down_mass,
                                     up_chemical_potential,
                                     down_chemical_potential,
                                     parameters.variables.temperature,
                                     &up_renormalized_chemical_potential,
                                     &down_renormalized_chemical_potential);
    double up_quark_density =
    QuarkDensity(up_mass,
                 up_renormalized_chemical_potential,
                 parameters.variables.temperature);

    double down_quark_density =
    QuarkDensity(down_mass,
                 down_renormalized_chemical_potential,
                 parameters.variables.temperature);

    BinodalPoint point;

    point.hadron_mass = hadron_mass;
    point.proton_density = proton_density;
    point.neutron_density = neutron_density;
    point.proton_chemical_potential = proton_chemical_potential;
    point.neutron_chemical_potential = neutron_chemical_potential;

    point.up_quark_mass = up_mass;
    point.down_quark_mass = down_mass;
    point.up_quark_density = up_quark_density;
    point.down_quark_density = down_quark_density;
    point.up_chemical_potential = up_chemical_potential;
    point.down_chemical_potential = down_chemical_potential;

    point.pressure = quark_pressure;

    return point;
}

double BinodalPointEquation(double  barionic_chemical_potential,
                            void   *params)
{

    binodal_parameters * p = (binodal_parameters *)params;

    double proton_chemical_potential =
    ProtonChemicalPotential(barionic_chemical_potential,
                            p->isovector_chemical_potential);

    double neutron_chemical_potential =
    NeutronChemicalPotential(barionic_chemical_potential,
                             p->isovector_chemical_potential);

    double hadron_mass;
    double hadron_pressure;
    double proton_density;
    double neutron_density;

    DetermineHadronPressureAndDensities(proton_chemical_potential,
                                        neutron_chemical_potential,
                                        p->hadron_vacuum_thermodynamic_potential,
                                        &hadron_mass,
                                        &proton_density,
                                        &neutron_density,
                                        &hadron_pressure);

    // From Gibbs conditions
    double up_chemical_potential =
    UpChemicalPotentialFromGibbsConditions(proton_chemical_potential,
                                           neutron_chemical_potential);

    double down_chemical_potential =
    DownChemicalPotentialFromGibbsConditions(proton_chemical_potential,
                                             neutron_chemical_potential);

    double up_quark_mass;
    double down_quark_mass;
    double quark_pressure;

    DetermineQuarkPressure(up_chemical_potential,
                           down_chemical_potential,
                           p->temperature,
                           p->quark_vacuum_thermodynamic_potential,
                           &up_quark_mass,
                           &down_quark_mass,
                           &quark_pressure);

    return hadron_pressure - quark_pressure;
}

void DetermineHadronPressureAndDensities(double proton_chemical_potential,
                                         double neutron_chemical_potential,
                                         double hadron_vacuum_potential,
                                         double *return_hadron_mass,
                                         double *return_proton_density,
                                         double *return_neutron_density,
                                         double *return_pressure)
{   double mass = NAN;
    double proton_density = NAN;
    double neutron_density = NAN;

    HadronMassAndDensitiesSolution(proton_chemical_potential,
                                   neutron_chemical_potential,
                                   &mass,
                                   &proton_density,
                                   &neutron_density);

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
}

void DetermineQuarkPressure(double up_chemical_potential,
                            double down_chemical_potential,
                            double temperature,
                            double quark_vacuum_thermodynamic_potential,
                            double *return_up_mass,
                            double *return_down_mass,
                            double *return_pressure)
{
    double up_mass = NAN;
    double down_mass = NAN;
    double up_renorm_chem_pot = NAN;
    double down_renorm_chem_pot = NAN;

    QuarkMassAndRenormChemPotSolution(up_chemical_potential,
                                      down_chemical_potential,
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

    return;
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
