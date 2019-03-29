//
//  FermiDiracDistributions.c
//  binodal
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#include <math.h>

#include "libdatafun/libdatafun.h"

#include "Parameters.h"
#include "FermiDiracDistributions.h"
#include "DefiniteIntegrals.h"

double FermiDiracDistributionForParticles(double energy,
                                          double chemical_potential,
                                          double temperature)
{
    return 1.0 / (1.0 + exp((energy - chemical_potential)/temperature));
}

double FermiDiracDistributionForAntiparticles(double energy,
                                              double chemical_potential,
                                              double temperature)
{
    return FermiDiracDistributionForParticles(energy,
                                              -chemical_potential,
                                              temperature);
}

double
FermiDiracDistributionIntegralFromBarionicDensityIntegrand(double momentum,
                                                           void * params)
{
    fermi_dirac_distrib_integrand * p = (fermi_dirac_distrib_integrand *) params;

    double E = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));

    double particle_term = FermiDiracDistributionForParticles(E,
                                                              p->chemical_potential,
                                                              p->temperature);
    double antiparticle_term =
        FermiDiracDistributionForAntiparticles(E,
                                               p->chemical_potential,
                                               p->temperature);
    return (particle_term - antiparticle_term) * pow(momentum, 2.0);
}

double
FermiDiracDistributionIntegralFromBarionicDensity(double temperature,
                                                  double mass,
                                                  double renorm_chem_pot,
                                                  double cutoff)
{
    fermi_dirac_distrib_integrand p;
    p.mass = mass;
    p.chemical_potential = renorm_chem_pot;
    p.temperature = temperature;

    gsl_function F;
    F.function = &FermiDiracDistributionIntegralFromBarionicDensityIntegrand;
    F.params = &p;

    parameters.fermi_dirac_integrals.upper_limit = cutoff;

    double integral = OnedimensionalIntegrator(&F, parameters.fermi_dirac_integrals);

    return integral;
}

double FermiDiracDistributionIntegralFromScalarDensityIntegrand(double momentum,
                                                                void * params)
{
    fermi_dirac_distrib_integrand * p = (fermi_dirac_distrib_integrand *) params;

    double E = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));

    double particle_term =
    FermiDiracDistributionForParticles(E,
                                       p->chemical_potential,
                                       p->temperature);

    double antiparticle_term =
    FermiDiracDistributionForAntiparticles(E,
                                           p->chemical_potential,
                                           p->temperature);

    return (particle_term + antiparticle_term) * pow(momentum, 2.0) / E;
}

double FermiDiracDistributionIntegralFromScalarDensity(double temperature,
                                                       double mass,
                                                       double renorm_chem_pot,
                                                       double cutoff)
{
    fermi_dirac_distrib_integrand p;
    p.mass = mass;
    p.chemical_potential = renorm_chem_pot;
    p.temperature = temperature;

    gsl_function F;
    F.function = &FermiDiracDistributionIntegralFromScalarDensityIntegrand;
    F.params = &p;

    parameters.fermi_dirac_integrals.upper_limit = cutoff;

    double integral =
    OnedimensionalIntegrator(&F, parameters.fermi_dirac_integrals);

    return F0(mass, cutoff) - integral;
}

double FermiDiracDistributionIntegralFromHadronEnergyIntegrand(double momentum,
                                                               void * params)
{
    fermi_dirac_distrib_integrand * p =
    (fermi_dirac_distrib_integrand *) params;

    double E = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));

    double particle_term =
    FermiDiracDistributionForParticles(E,
                                       p->chemical_potential,
                                       p->temperature);

    double antiparticle_term =
    FermiDiracDistributionForAntiparticles(E,
                                           p->chemical_potential,
                                           p->temperature);

    return (particle_term + antiparticle_term) * pow(momentum, 4.0) / E;
}

double FermiDiracDistributionIntegralFromHadronEnergy(double temperature,
                                                      double mass,
                                                      double renorm_chem_pot,
                                                      double cutoff)
{
    fermi_dirac_distrib_integrand p;
    p.mass = mass;
    p.chemical_potential = renorm_chem_pot;
    p.temperature = temperature;

    gsl_function F;
    F.function = &FermiDiracDistributionIntegralFromHadronEnergyIntegrand;
    F.params = &p;

    parameters.fermi_dirac_integrals.upper_limit = cutoff;

    double integral =
    OnedimensionalIntegrator(&F, parameters.fermi_dirac_integrals);

    return F2(mass, cutoff) - integral;
}

double FermiDiracDistributionIntegralFromHadronEntropyIntegrand(double momentum,
                                                                void * params)
{
    fermi_dirac_distrib_integrand * p =
    (fermi_dirac_distrib_integrand *) params;

    double E = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));

    double np =
    FermiDiracDistributionForParticles(E,
                                       p->chemical_potential,
                                       p->temperature);

    double nap =
    FermiDiracDistributionForAntiparticles(E,
                                           p->chemical_potential,
                                           p->temperature);

    double ln_np = log(np);
    double ln_1mnp = log1p(-np);

    double ln_nap = log(nap);
    double ln_1mnap = log1p(-nap);

    double expr = 0;
    
    // In the cases where np goes to 1 or nap goes to 0,
    // the expressions below go to zero, but we will have
    // problems with infinity. Add the expression results
    // only when the logarithms have proper values.
    if (!isinf(ln_1mnp))
        expr += (1.0 - np) * ln_1mnp;
        
    if (!isinf(ln_nap))
        expr += nap * ln_nap;

    expr += np * ln_np + (1.0 - nap) * ln_1mnap;
    
    return expr * pow(momentum, 2.0);
}

double FermiDiracDistributionIntegralFromHadronEntropy(double temperature,
                                                       double mass,
                                                       double renorm_chem_pot,
                                                       double cutoff)
{
    fermi_dirac_distrib_integrand p;
    p.mass = mass;
    p.chemical_potential = renorm_chem_pot;
    p.temperature = temperature;

    gsl_function F;
    F.function = &FermiDiracDistributionIntegralFromHadronEntropyIntegrand;
    F.params = &p;

    parameters.fermi_dirac_integrals.upper_limit = cutoff;

    double integral =
    OnedimensionalIntegrator(&F, parameters.fermi_dirac_integrals);

    return integral;
}

double
FermiDiracDistributionIntegralFromQuarkThermodynamicPotentialIntegrand(double momentum,
                                                                       void * params)
{
    fermi_dirac_distrib_integrand * p =
    (fermi_dirac_distrib_integrand *) params;

    double E = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));

    double np =
    FermiDiracDistributionForParticles(E,
                                       p->chemical_potential,
                                       p->temperature);

    double nap =
    FermiDiracDistributionForAntiparticles(E,
                                           p->chemical_potential,
                                           p->temperature);

    double ln_1mnp = log1p(-np);

    double ln_1mnap = log1p(-nap);

    return (ln_1mnap + ln_1mnp) * pow(momentum, 2.0);
}

double
FermiDiracDistributionIntegralFromQuarkThermodynamicPotential(double temperature,
                                                              double mass,
                                                              double renorm_chem_pot,
                                                              double cutoff)
{
    fermi_dirac_distrib_integrand p;
    p.mass = mass;
    p.chemical_potential = renorm_chem_pot;
    p.temperature = temperature;

    gsl_function F;
    F.function =
    &FermiDiracDistributionIntegralFromQuarkThermodynamicPotentialIntegrand;
    F.params = &p;

    parameters.fermi_dirac_integrals.upper_limit = cutoff;

    double integral =
    OnedimensionalIntegrator(&F, parameters.fermi_dirac_integrals);

    // TODO: Check this upper limit, I think this should be the hadron cutoff
    return F_E(mass, cutoff) - temperature * integral;
}
