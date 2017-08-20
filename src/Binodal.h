//
//  Binodal.h
//  binodal
//
//  Created by Clebson Graeff on 2017-06-20.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#ifndef Binodal_h
#define Binodal_h

typedef struct _BinodalPoint{
    double hadron_mass;
    double proton_density;
    double neutron_density;
    double proton_chemical_potential;
    double neutron_chemical_potential;

    double up_quark_mass;
    double down_quark_mass;
    double up_quark_density;
    double down_quark_density;
    double up_chemical_potential;
    double down_chemical_potential;

    double pressure;
} BinodalPoint;

BinodalPoint DetermineBinodalPoint(double temperature,
                                   double proton_fraction,
                                   double hadron_vacuum_thermodynamic_potential,
                                   double quark_vacuum_thermodynamic_potential);

void DetermineHadronPressureAndDensities(double proton_chemical_potential,
                                         double neutron_chemical_potential,
                                         double hadron_vacuum_potential,
                                         double *return_hadron_mass,
                                         double *return_proton_density,
                                         double *return_neutron_density,
                                         double *return_pressure);

void DetermineQuarkPressure(double up_chemical_potential,
                            double down_chemical_potential,
                            double temperature,
                            double quark_vacuum_thermodynamic_potential,
                            double *return_up_mass,
                            double *return_down_mass,
                            double *return_pressure);

double BarionicChemicalPotential(double proton_chemical_potential,
                                 double neutron_chemical_potential);

double IsovectorChemicalPotential(double proton_chemical_potential,
                                  double neutron_chemical_potential);

double
UpChemicalPotentialFromGibbsConditions(double proton_chemical_potential,
                                       double neutron_chemical_potential);

double
DownChemicalPotentialFromGibbsConditions(double proton_chemical_potential,
                                         double neutron_chemical_potential);

#endif /* Binodal_h */
 
