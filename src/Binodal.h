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

typedef struct _BinodalBissectionParameters{
    double lower_bound;
    double upper_bound;
    double step_size;
    double max_iterations;
    double abs_error;
    double rel_error;
} BinodalBissectionParameters;

BinodalPoint DetermineBinodalPoint(double temperature,
                                   double isovector_chemical_potential,
                                   double hadron_vacuum_thermodynamic_potential,
                                   double quark_vacuum_thermodynamic_potential);

BinodalPoint
DetermineBinodalPointByBissection(double temperature,
                                  double isovector_chemical_potential,
                                  double hadron_vacuum_thermodynamic_potential,
                                  double quark_vacuum_thermodynamic_potential,
                                  double transition_bar_chem_pot_lower_bound,
                                  double transition_bar_chem_pot_upper_bound,
                                  double hadron_mass_guess,
                                  double proton_density_guess,
                                  double neutron_density_guess,
                                  double up_mass_guess,
                                  double down_mass_guess);

void DetermineHadronPressureAndDensities(double proton_chemical_potential,
                                         double neutron_chemical_potential,
                                         double hadron_vacuum_potential,
                                         double hadron_mass_guess,
                                         double proton_density_guess,
                                         double neutron_density_guess,
                                         double *return_hadron_mass,
                                         double *return_proton_density,
                                         double *return_neutron_density,
                                         double *return_pressure);

void DetermineQuarkPressure(double up_chemical_potential,
                            double down_chemical_potential,
                            double temperature,
                            double quark_vacuum_thermodynamic_potential,
                            double up_mass_guess,
                            double down_mass_guess,
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
 
