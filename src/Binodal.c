//
//  Binodal.c
//  binodal
//
//  Created by Clebson Graeff on 2017-06-20.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#include <gsl/gsl_vector.h>

#include "libdatafun/libdatafun.h"

#include "QuarkPhaseEOS.h"
#include "HadronPhaseEOS.h"
#include "Parameters.h"

#include "Binodal.h"

typedef struct _binodal_parameters {
    double proton_fraction;
    double temperature;

    double hadron_vacuum_thermodynamic_potential;
    double quark_vacuum_thermodynamic_potential;
} binodal_parameters;

double BinodalPointEquation(double  barionic_density,
                            void   *params);

BinodalPoint DetermineBinodalPoint(double temperature,
                                   double proton_fraction,
                                   double hadron_vacuum_thermodynamic_potential,
                                   double quark_vacuum_thermodynamic_potential)
{

    binodal_parameters params;
    params.temperature = temperature;
    params.proton_fraction = proton_fraction;
    params.hadron_vacuum_thermodynamic_potential = hadron_vacuum_thermodynamic_potential;
    params.quark_vacuum_thermodynamic_potential = quark_vacuum_thermodynamic_potential;

    gsl_function F;
    F.function = &BinodalPointEquation;
    F.params = &params;

    double barionic_density;

    int status = UnidimensionalRootFinder(&F,
                                          parameters.binodal_rootfinding_params,
                                          &barionic_density);

    if (status != 0){
        printf("%s:%d: Problem with rootfinding.\n", __FILE__, __LINE__);
        abort();
    }

    // Determine results to return
    double hadron_mass = NAN;
    double pressure = NAN;
    double proton_chemical_potential = NAN;
    double neutron_chemical_potential = NAN;

    DetermineHadronPressureAndChemPots(barionic_density,
                                       proton_fraction,
                                       hadron_vacuum_thermodynamic_potential,
                                       &hadron_mass,
                                       &proton_chemical_potential,
                                       &neutron_chemical_potential,
                                       &pressure);

    BinodalPoint point;

    point.hadron_mass = hadron_mass;
    point.barionic_density = barionic_density;
    point.proton_chemical_potential = proton_chemical_potential;

    point.neutron_chemical_potential = neutron_chemical_potential;

    point.pressure = pressure;

    return point;
}

double BinodalPointEquation(double  barionic_density,
                            void   *params)
{

    binodal_parameters * p = (binodal_parameters *)params;

    double hadron_mass;
    double hadron_pressure;
    double proton_chemical_potential;
    double neutron_chemical_potential;

    DetermineHadronPressureAndChemPots(barionic_density,
                                       p->proton_fraction,
                                       p->hadron_vacuum_thermodynamic_potential,
                                       &hadron_mass,
                                       &proton_chemical_potential,
                                       &neutron_chemical_potential,
                                       &hadron_pressure);

    // From Gibbs conditions
    double up_chemical_potential =
    (2.0 * proton_chemical_potential - neutron_chemical_potential) / 3.0;

    double down_chemical_potential =
    (-proton_chemical_potential + 2.0 * neutron_chemical_potential) / 3.0;

    double quark_pressure;

    DetermineQuarkPressure(up_chemical_potential,
                           down_chemical_potential,
                           p->temperature,
                           p->quark_vacuum_thermodynamic_potential,
                           &quark_pressure);

    return hadron_pressure - quark_pressure;
}

void DetermineHadronPressureAndChemPots(double barionic_density,
                                        double proton_fraction,
                                        double hadron_vacuum_potential,
                                        double *return_hadron_mass,
                                        double *return_proton_chem_pot,
                                        double *return_neutron_chem_pot,
                                        double *return_pressure)
{
    double proton_density = proton_fraction * barionic_density;
    double neutron_density = (1.0 - proton_fraction) * barionic_density;

    double proton_fermi_momentum = HadronFermiMomentum(proton_density);
    double neutron_fermi_momentum = HadronFermiMomentum(neutron_density);

    double mass = HadronSolveGapEquation(proton_density,
                                         neutron_density,
                                         proton_fermi_momentum,
                                         neutron_fermi_momentum);

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

    double proton_chemical_potential =
    ProtonChemicalPotential(proton_fermi_momentum,
                            scalar_density,
                            mass,
                            barionic_density,
                            proton_density,
                            neutron_density);

    double neutron_chemical_potential =
    NeutronChemicalPotential(neutron_fermi_momentum,
                             scalar_density,
                             mass,
                             barionic_density,
                             proton_density,
                             neutron_density);

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
    *return_proton_chem_pot = proton_chemical_potential;
    *return_neutron_chem_pot = neutron_chemical_potential;
}

void DetermineQuarkPressure(double up_chemical_potential,
                            double down_chemical_potential,
                            double temperature,
                            double quark_vacuum_thermodynamic_potential,
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
