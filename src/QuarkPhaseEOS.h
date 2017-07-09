//
//  QuarkPhaseEOS.h
//  binodal
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#ifndef QuarkPhaseEOS_h
#define QuarkPhaseEOS_h

// Parameters for the determination of
// entropy by integration
typedef struct _QuarkEntropyIntegrationParameters{
    double lower_limit;
    double upper_limit;
    double abs_error;
    double rel_error;
    int max_sub_interval;
    int integration_key;
} QuarkEntropyIntegrationParameters;

// Parameterization variables
typedef struct _QuarkModelParameters{
    char * parameters_set_identifier;
    char * parameters_set_origin;   // Where the set was taken from
    double G_S;                     // scalar-isoscalar coupling (fm^2)
    double G_V;                     // vector-isoscalar coupling (fm^2)
    double cutoff;                  // \Lambda (MeV)
    double bare_mass;               // (MeV)
} QuarkModelParameters;

typedef struct _QuarkMassAndRenormChemPotSolPar{
    double up_mass_guess;
    double down_mass_guess;
    double abs_error;
    double rel_error;
    int max_iter;
} QuarkMassAndRenormChemPotSolParams;

// Parameters for the determination of
// entropy by integration
typedef struct _EntropyParameters{
    double lower_limit;
    double upper_limit;
    double abs_error;
    double rel_error;
    int max_sub_interval;
    int integration_key;
} EntropyParameters;

typedef struct _QuarkVacuumMassDeterminationParameters{
    double up_vacuum_mass_guess;
    double down_vacuum_mass_guess;

    int max_iter;

    double abs_error;
    double rel_error;
} QuarkVacuumMassDeterminationParameters;

typedef struct _QuarkThermodynamicPotMinDetermParams{
    int max_iter;
    double tolerance;
    double up_vacuum_mass_guess;
    double down_vacuum_mass_guess;
    double up_mass_step;
    double down_mass_step;
} QuarkThermodynamicPotMinDetermParams;

typedef struct _QuarkRenormChemPotSolutionParameters{
    double up_renorm_chem_pot_guess;
    double down_renorm_chem_pot_guess;

    double abs_error;
    double rel_error;
    double max_iter;
} QuarkRenormChemPotSolutionParameters;

double QuarkDensity(double mass,
                    double renormalized_chemical_potential,
                    double temperature);

double QuarkScalarDensity(double temperature,
                          double mass,
                          double renormalized_chemical_potential);

double QuarkThermodynamicPotential(double up_mass,
                                   double down_mass,
                                   double up_chemical_potential,
                                   double down_chemical_potential,
                                   double up_renormalized_chemical_potential,
                                   double down_renormalized_chemical_potential,
                                   double temperature);

double
QuarkThermodynamicPotentialFreeGasTerm(double mass,
                                       double chemical_potential,
                                       double renormalized_chemical_potential,
                                       double temperature);

double QuarkPressure(double regularized_thermodynamic_potential,
                     double temperature);

double QuarkEnergyDensity(double regularized_thermodynamic_potential,
                          double chemical_potential,
                          double barionic_density,
                          double temperature,
                          double entropy);

double QuarkEntropy(double mass,
                    double temperature,
                    double renormalized_chemical_potential);

double QuarkEntropyIntegrand(double momentum,
                             void * parameters);

typedef struct _quark_renorm_chem_pot_equation_input{
    double chemical_potential;
    double mass;
	double temperature;
} quark_renorm_chem_pot_equation_input;

int ZeroedRenormalizedQuarkChemPotEquation(const gsl_vector   *x,
                                           void *params,
                                           gsl_vector *return_values);

void QuarkSelfConsistentRenormChemPot(double up_quark_mass,
                                      double down_quark_mass,
                                      double up_chemical_potential,
                                      double down_chemical_potential,
                                      double temperature,
                                      double *return_up_renorm_chem_pot,
                                      double *return_down_renorm_chem_pot);
double QuarkZeroedGapEquation(double mass,
                              double up_scalar_density,
                              double down_scalar_density);

void QuarkThermodynamicPotentialMinimum(double * up_mass_at_minimum,
                                        double * down_mass_at_minimum);

void QuarkVacuumMassDetermination(double * up_vacuum_mass,
                                  double * down_vacuum_mass);

int QuarkVacuumMassDeterminationEquation(const gsl_vector   *x,
                                         void *params,
                                         gsl_vector *return_values);

double QuarkFermiMomentum(double mass, double renormalized_chemical_potential);

int QuarkMassAndRenormChemPotSolution(double up_chemical_potential,
                                      double down_chemical_potential,
                                      double * return_up_mass,
                                      double * return_down_mass,
                                      double * return_up_renorm_chem_pot,
                                      double * return_down_renorm_chem_pot);
#endif /* QuarkPhaseEOS_h */
