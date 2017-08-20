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

typedef struct _HadronMassAndDensitiesSolutionParams{
    double mass_guess;
    double zero_mass_chem_pot;
    double proton_density_guess;
    double neutron_density_guess;

    int max_iter;
    double abs_error;
    double rel_error;
} HadronMassAndDensitiesSolutionParams;

double HadronSolveGapEquation(double proton_density,
                              double neutron_density,
                              double proton_fermi_momentum,
                              double neutron_fermi_momentum);

double HadronProtonFraction(double proton_barionic_density,
                            double neutron_barionic_density);

double HadronScalarDensity(double mass,
                           double fermi_momentum,
                           double cutoff);

double HadronVacuumScalarDensity();

double ProtonChemicalPotentialEquation(double proton_fermi_momentum,
                                       double scalar_density,
                                       double mass,
                                       double proton_density,
                                       double neutron_density);

double NeutronChemicalPotentialEquation(double neutron_fermi_momentum,
                                        double scalar_density,
                                        double mass,
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

int HadronMassAndDensitiesSolution(double proton_chemical_potential,
                                   double neutron_chemical_potential,
                                   double * return_mass,
                                   double * return_proton_density,
                                   double * return_neutron_density);

double ProtonChemicalPotential(double barionic_chemical_potential,
                               double isovector_chemical_potential);
double NeutronChemicalPotential(double barionic_chemical_potential,
                                double isovector_chemical_potential);

double HadronPhaseAsymmetry(double proton_density, double neutron_density);

#endif /* HadronPhaseEOS_h */
