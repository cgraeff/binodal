//
//  Binodal.h
//  binodal
//
//  Created by Clebson Graeff on 2017-06-20.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#ifndef Binodal_h
#define Binodal_h

void BinodalPoint(double temperature,
                  double proton_fraction,
                  double hadron_vacuum_thermodynamic_potential,
                  double quark_vacuum_thermodynamic_potential,
                  double *return_barionic_density,
                  double *return_barionic_chemical_potential,
                  double *return_isovector_chemical_potential,
                  double *return_pressure);

void DetermineHadronPressureAndChemPots(double barionic_density,
                                        double proton_fraction,
                                        double hadron_vacuum_potential,
                                        double *return_pressure,
                                        double *return_proton_chem_pot,
                                        double *return_neutron_chem_pot);

void DetermineQuarkPressure(double up_chemical_potential,
                            double down_chemical_potential,
                            double temperature,
                            double quark_vacuum_thermodynamic_potential,
                            double *return_pressure);

#endif /* Binodal_h */
 
