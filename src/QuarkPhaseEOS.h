//
//  EOS.h
//  NJLv
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#ifndef QuarkPhaseEOS_h
#define QuarkPhaseEOS_h

#include "SimultaneousSolution.h"

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

double QuarkSelfConsistentRenormChemPot(double chemical_potential,
                                        double up_barionic_density,
                                        double down_barionic_density,
                                        double temperature);

double QuarkBarionicDensity(double mass,
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
QuarkThermodynamicPotentialFreeGasContribution(double mass,
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

int ZeroedRenormalizedQuarkChemicalPotentialEquation(const gsl_vector   *x,
                                                        void *params,
                                                        gsl_vector *return_values);

void QuarkSelfConsistentRenormalizedChemicalPotential(RenormChemPotSolutionParameters params,
                                                      double up_quark_mass,
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

double QuarkVacuumMassEquation(double mass, void * input);

double QPFermiMomentum(double mass, double renormalized_chemical_potential);

int TestMassAndRenormChemPotSimultaneousSolution(double up_chemical_potential,
                                                 double down_chemical_potential,
                                                 double up_mass_guess,
                                                 double down_mass_guess,
                                                 double abs_error,
                                                 double rel_error,
                                                 int max_iter,
                                                 double * return_up_mass,
                                                 double * return_down_mass,
                                                 double * return_up_renorm_chem_pot,
                                                 double * return_down_renorm_chem_pot);
#endif /* QuarkPhaseEOS_h */
