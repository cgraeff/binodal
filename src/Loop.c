//
//  Loop.c
//  NJLv
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#include "libdatafun/libdatafun.h"

#include "CommandlineOptions.h"
#include "Parameters.h"
#include "Constants.h"
#include "QuarkPhaseEOS.h"
#include "HadronPhaseEOS.h"

int SolveFiniteTemperatureEOS(){

    // Print name of parametrization
    if (options.verbose){
        printf("Calculation performed with %s quark phase parameters set\n"
               "and %s hadron phase paarameters set.\n"
               "\ttemperature: %4.2f (MeV)\n"
               "\tproton fraction: %4.2f\n",
               parameters.quark_model.parameters_set_identifier,
               parameters.hadron_model.parameters_set_identifier,
               parameters.variables.temperature,
               parameters.variables.proton_fraction);
    }

    // Vacuum mass determination
    if (options.verbose)
        printf("Determining the vacuum mass and bag constant ...\n");

    double quark_vacuum_mass = QuarkVacuumMassDetermination();

    double quark_vacuum_thermodynamic_potential =
	    QuarkThermodynamicPotential(quark_vacuum_mass,
                                    0.0,
                                    0.0,
                                    0.0);

    // TODO: should the hadron mass be calculated?
    double hadron_vacuum_mass = parameters.hadron_model.nucleon_mass;

    // double hadron_vacuum_mass = HadronVacuumMassDetermination();
    double hadron_vacuum_thermodynamic_potential = 0.0;
/*	    HadronThermodynamicPotential(hadron_vacuum_mass,
                                     0.0,
                                     0.0,
                                     0.0);
double HadronThermodynamicPotential(double scalar_density,
                                    double barionic_density,
                                    double proton_density,
                                    double neutron_density,
                                    double proton_chemical_potential,
                                    double neutron_chemical_potential,
                                    double kinectic_energy_density)
 */
    if (options.verbose){
        printf("\tQuark vacuum mass: %f\n", quark_vacuum_mass);
        printf("\tHadron vacuum mass: %f\n", hadron_vacuum_mass);
    }

    if (options.verbose)
        printf("Solving gap equation and equations of state ...\n");

    double return_barionic_density;
    double return_hadron_mass;
    double return_quark_mass;
    double return_proton_chemical_potential;
    double return_neutron_chemical_potential;
    double return_pressure;

    SimultaneousSolution(parameters.simultaneous_solution,
                         quark_vacuum_thermodynamic_potential,
                         hadron_vacuum_thermodynamic_potential,
                         &return_barionic_density,
                         &return_hadron_mass,
                         &return_quark_mass,
                         &return_proton_chemical_potential,
                         &return_neutron_chemical_potential,
                         &return_pressure);

    if (options.verbose)
        printf("Done!\n");

    return 0;
}
