
void BinodalPoint(double temperature,
                  double proton_fraction,
                  double hadron_vacuum_thermodynamic_potential,
                  double quark_vacuum_thermodynamic_potential,
                  double *return_barionic_density,
                  double *return_barionic_chemical_potential,
                  double *return_isovector_chemical_potential,
                  double *return_pressure);

void DetermineHadronPressureAndChemicalPotentials(double barionic_density,
                                                  double proton_fraction,
                                                  double hadron_vacuum_potential,
                                                  double *return_pressure,
                                                  double *return_proton_chemical_potential,
                                                  double *return_neutron_chemical_potential);

void DetermineQuarkPressure(double up_chemical_potential,
                            double down_chemical_potential,
                            double temperature,
                            double quark_vacuum_thermodynamic_potential,
                            double *return_pressure);
