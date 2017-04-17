//
//  EOS.h
//  NJLv
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#ifndef QuarkPhaseEOS_h
#define QuarkPhaseEOS_h

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


double QuarkBarionicDensity(double mass,
                            double renormalized_chemical_potential,
                            double temperature);

double QuarkScalarDensity(double temperature,
                          double mass,
                          double renormalized_chemical_potential);

double QuarkThermodynamicPotential(double mass,
                                   double chemical_potential,
                                   double renormalized_chemical_potential,
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

typedef struct _quark_gap_eq_input_params{
    double renormalized_chemical_potential;
	double temperature;
} quark_gap_eq_input_params;

typedef struct _quark_renorm_chem_pot_equation_input{
    double chemical_potential;
    double mass;
	double temperature;
} quark_renorm_chem_pot_equation_input;

double QuarkZeroedGapEquation(double mass,
                         void * params);

double
QuarkZeroedRenormalizedChemicalPotentialEquation(double  renorm_chemical_potential,
                                                 void   *params);

double QuarkVacuumMassDetermination();
double QuarkVacuumMassEquation(double mass, void * input);

double QPFermiMomentum(double mass, double renormalized_chemical_potential);
#endif /* QuarkPhaseEOS_h */
