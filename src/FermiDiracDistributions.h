//
//  FermiDiracDistributions.h
//  NJLv
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#ifndef FermiDiracDistributions_h
#define FermiDiracDistributions_h

typedef struct _fermi_dirac_distrib_integrand{
    double mass;
    double chemical_potential;
    double temperature;
} fermi_dirac_distrib_integrand;

double FermiDiracDistributionForParticles(double energy,
                                          double chemical_potential,
                                          double temperature);

double FermiDiracDistributionForAntiparticles(double energy,
                                              double chemical_potential,
                                              double temperature);

double
FermiDiracDistributionIntegralFromBarionicDensity(double temperature,
                                                  double mass,
                                                  double renorm_chem_pot,
                                                  double cutoff);

double FermiDiracDistributionIntegralFromScalarDensity(double temperature,
                                                       double mass,
                                                       double renorm_chem_pot,
                                                       double cutoff);

double FermiDiracDistributionIntegralFromHadronEnergy(double temperature,
                                                      double mass,
                                                      double renorm_chem_pot,
                                                      double cutoff);

double FermiDiracDistributionIntegralFromHadronEntropy(double temperature,
                                                       double mass,
                                                       double renorm_chem_pot,
                                                       double cutoff);

double
FermiDiracDistributionIntegralFromQuarkThermodynamicPotential(double temperature,
                                                              double mass,
                                                              double renorm_chem_pot,
                                                              double cutoff);

#endif /* FermiDiracDistributions_h */
