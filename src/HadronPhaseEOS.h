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
    double initial_mass_guess;
    double initial_proton_density_guess;
    double initial_neutron_density_guess;

    double zero_mass_tolerance;
    double zero_dens_tolerance;

    int max_iter;
    double abs_error;
    double rel_error;
} HadronMassAndDensitiesSolutionParams;

// TODO: Check below, may be *_renorm_chem_pot;
typedef struct _hadron_mass_and_renorm_chem_pot_input_params{

    double proton_chemical_potential;
    double neutron_chemical_potential;

} hadron_mass_and_renorm_chem_pot_input_params;

typedef struct _GridRootFinderParameters{
    int num_pts_mass;
    double min_mass;
    double max_mass;

    int num_pts_dens;
    double min_density;
    double max_density;

    double zero_tol;
} GridRootFinderParameters;
int GridRootFinder(double proton_chemical_potential,
                   double neutron_chemical_potential,
                   double * return_hadron_mass,
                   double * return_proton_density,
                   double * return_neutron_density,
                   int * return_num_solutions);

double HadronSolveGapEquation(double proton_density,
                              double neutron_density,
                              double proton_fermi_momentum,
                              double neutron_fermi_momentum);

double HadronProtonFraction(double proton_barionic_density,
                            double neutron_barionic_density);

double HadronZeroTemperatureScalarDensity(double mass,
                                          double fermi_momentum,
                                          double cutoff);

double HadronScalarDensity(double mass,
                           double renorm_chem_pot,
                           double temperature,
                           double density);

double HadronVacuumScalarDensity();

double ProtonChemicalPotentialEquation(double proton_fermi_momentum,
                                       double proton_scalar_density,
                                       double neutron_scalar_density,
                                       double mass,
                                       double proton_density,
                                       double neutron_density);

double NeutronChemicalPotentialEquation(double neutron_fermi_momentum,
                                        double proton_scalar_density,
                                        double neutron_scalar_density,
                                        double mass,
                                        double proton_density,
                                        double neutron_density);

double HadronKinecticEnergyDensity(double temperature,
                                   double mass,
                                   double proton_density,
                                   double neutron_density,
                                   double proton_renorm_chem_pot,
                                   double neutron_renorm_chem_pot);

double HadronVacuumKinecticEnergyDensity();

double HadronVacuumEnergyDensity();

double HadronEnergyDensity(double pressure,
                           double proton_chemical_potential,
                           double neutron_chemical_potential,
                           double proton_density,
                           double neutron_density);

double HadronPressure(double termodynamic_potential);

double HadronThermodynamicPotential(double mass,
                                    double proton_scalar_density,
                                    double neutron_scalar_density,
                                    double proton_barionic_density,
                                    double neutron_barionic_density,
                                    double proton_chemical_potential,
                                    double neutron_chemical_potential,
                                    double proton_renorm_chem_pot,
                                    double neutron_renorm_chem_pot,
                                    double kinectic_energy_density);

double HadronFermiMomentumFromBarionicDensity(double density);
double HadronBarionicDensityFromFermiMomentum(double fermi_momentum);

int HadronMassAndDensitiesSolution(double proton_chemical_potential,
                                   double neutron_chemical_potential,
                                   double hadron_mass_guess,
                                   double proton_density_guess,
                                   double neutron_density_guess,
                                   double * return_mass,
                                   double * return_proton_density,
                                   double * return_neutron_density);

int HadronMassAndDensitiesSolutionEquation(const gsl_vector   *x,
                                           void *params,
                                           gsl_vector *return_values);

double ProtonChemicalPotential(double barionic_chemical_potential,
                               double isovector_chemical_potential);
double NeutronChemicalPotential(double barionic_chemical_potential,
                                double isovector_chemical_potential);

double HadronPhaseAsymmetry(double proton_density, double neutron_density);

double HadronZeroedGapEquation(double mass,
                               double proton_density,
                               double neutron_density,
                               double proton_scalar_density,
                               double neutron_scalar_denstiy);


#endif /* HadronPhaseEOS_h */
