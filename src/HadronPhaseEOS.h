
#ifndef HadronPhaseEOS_h
#define HadronPhaseEOS_h

typedef struct _hadron_gap_eq_input_params{
    double proton_fermi_momentum;
    double neutron_fermi_momentum;
    double proton_density;
    double neutron_density;
    double renormalized_chemical_potential;
} hadron_gap_eq_input_params;


double HadronZeroedGapEquation(double mass,
                               void * params);

double HadronSolveGapEquation(UnidimensionalRootFindingParameters rootfinding_params,
                              double proton_density,
                              double neutron_density,
                              double proton_fermi_momentum,
                              double neutron_fermi_momentum);

double HadronScalarDensity(double mass,
                           double fermi_momentum,
                           double cutoff);

double HadronVacuumScalarDensity();

double ProtonChemicalPotential(double proton_fermi_momentum,
							   double scalar_density,
							   double mass,
							   double barionic_density,
							   double proton_density,
							   double neutron_density);

double NeutronChemicalPotential(double neutron_fermi_momentum,
								double scalar_density,
								double mass,
								double barionic_density,
								double proton_density,
								double neutron_density);

double HadronKinecticEnergyDensity(double mass,
                                   double proton_fermi_momentum,
                                   double neutron_fermi_momentum);

double HadronVacuumKinecticEnergyDensity();

double HadronVacuumEnergyDensity();

double HadronEnergyDensity(double pressure,
	                	   double proton_chemical_potential,
					       double neutron_chemical_potential,
					       double proton_density,
					       double neutron_density);

double HadronPressure(double termodynamic_potential);

double HadronThermodynamicPotential(double scalar_density,
                                    double barionic_density,
                                    double proton_density,
                                    double neutron_density,
                                    double proton_chemical_potential,
                                    double neutron_chemical_potential,
                                    double kinectic_energy_density);

double HPFermiMomentum(double density);

#endif /* HadronPhaseEOS_h */
