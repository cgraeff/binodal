//
//  HadronPhaseEOS.h
//  binodal
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#ifndef HadronPhaseEOS_h
#define HadronPhaseEOS_h

// Hadron phase parameters
typedef struct _HadronModelParameters{
    char * parameters_set_identifier;
    char * parameters_set_origin;

    double G_S;             // scalar-isoscalar coupling (fm^2)
    double G_V;             // vector-isoscalar coupling (fm^2)
    double G_RHO;           // vector-isovector coupling (fm^2)
    double G_SV;            // scalar-isovector coupling (fm^8)
    double G_VRHO;          // (fm^8)
    double G_SRHO;          // (fm^8)

    double cutoff;          // (MeV)
    double nucleon_mass;    // (MeV)
    double bare_mass;       // (MeV)
} HadronModelParameters;

double HadronZeroedGapEquation(double mass,
                               void * params);

double HadronSolveGapEquation(double proton_density,
                              double neutron_density,
                              double proton_fermi_momentum,
                              double neutron_fermi_momentum);

double HadronScalarDensity(double mass,
                           double fermi_momentum,
                           double cutoff);

double HadronVacuumScalarDensity();

double ProtonChemicalPotential(double proton_fermi_momentum,
                               double scalar_density,
                               double mass,
                               double barionic_density,
                               double proton_density,
                               double neutron_density);

double NeutronChemicalPotential(double neutron_fermi_momentum,
                                double scalar_density,
                                double mass,
                                double barionic_density,
                                double proton_density,
                                double neutron_density);

double HadronKinecticEnergyDensity(double mass,
                                   double proton_fermi_momentum,
                                   double neutron_fermi_momentum);

double HadronVacuumKinecticEnergyDensity();

double HadronVacuumEnergyDensity();

double HadronEnergyDensity(double pressure,
                           double proton_chemical_potential,
                           double neutron_chemical_potential,
                           double proton_density,
                           double neutron_density);

double HadronPressure(double termodynamic_potential);

double HadronThermodynamicPotential(double scalar_density,
                                    double barionic_density,
                                    double proton_density,
                                    double neutron_density,
                                    double proton_chemical_potential,
                                    double neutron_chemical_potential,
                                    double kinectic_energy_density);

double HadronFermiMomentum(double density);

#endif /* HadronPhaseEOS_h */
