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
#include "Binodal.h"

int SolveFiniteTemperatureEOS(){

    // Print name of parametrization
    if (options.verbose){
        printf("Calculation performed with %s quark phase parameters set\n"
               "and %s hadron phase paarameters set.\n"
               "\ttemperature: %4.2f (MeV)\n",
               parameters.quark_model.parameters_set_identifier,
               parameters.hadron_model.parameters_set_identifier,
               parameters.variables.temperature);
    }

    // Vacuum mass and potential determination
    if (options.verbose)
        printf("Determining the vacuum mass and bag constant ...\n");

    double hadron_vacuum_thermodynamic_potential =
    HadronVacuumEnergyDensity();

    double up_vacuum_mass;
    double down_vacuum_mass;

    QuarkVacuumMassDetermination(&up_vacuum_mass, &down_vacuum_mass);

    double quark_vacuum_thermodynamic_potential =
    QuarkThermodynamicPotential(up_vacuum_mass,
                                down_vacuum_mass,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0);

    if (options.verbose){
        printf("\tUp quark vacuum mass: %f\n", up_vacuum_mass);
        printf("\tDown quark vacuum mass: %f\n", down_vacuum_mass);
        printf("\tQuark vacuum thermodynamic potential: %f\n",
               quark_vacuum_thermodynamic_potential);
        printf("\tHadron vacuum thermodynamic potential: %f\n",
               hadron_vacuum_thermodynamic_potential);
    }

    if (options.verbose)
        printf("Solving gap equation and equations of state ...\n");

    // LOOP STARTS

    gsl_vector * proton_fraction_vector =
    gsl_vector_alloc(parameters.variables.num_points);

    gsl_vector * barionic_density_vector =
    gsl_vector_alloc(parameters.variables.num_points);

    gsl_vector * barionic_chemical_potential_vector =
    gsl_vector_alloc(parameters.variables.num_points);

    gsl_vector * isovector_chemical_potential_vector =
    gsl_vector_alloc(parameters.variables.num_points);

    gsl_vector * pressure_vector =
    gsl_vector_alloc(parameters.variables.num_points);

    double proton_fraction = parameters.variables.min_proton_fraction;
    double proton_fraction_step = Step(parameters.variables.min_proton_fraction,
                                       parameters.variables.max_proton_fraction,
                                       parameters.variables.num_points);

    for (int i = 0; i < parameters.variables.num_points; i++){

        if (options.verbose)
            printf("\tProton fraction: %f\r", proton_fraction);

        double barionic_density = NAN;
        double barionic_chemical_potential = NAN;
        double isovector_chemical_potential = NAN;
        double pressure = NAN;

        BinodalPoint(parameters.variables.temperature,
                     proton_fraction,
                     hadron_vacuum_thermodynamic_potential,
                     quark_vacuum_thermodynamic_potential,
                     &barionic_density,
                     &barionic_chemical_potential,
                     &isovector_chemical_potential,
                     &pressure);

        gsl_vector_set(barionic_chemical_potential_vector, i, barionic_chemical_potential);
        gsl_vector_set(isovector_chemical_potential_vector, i, isovector_chemical_potential);
        gsl_vector_set(barionic_density_vector, i, barionic_density);
        gsl_vector_set(pressure_vector, i, pressure);
        gsl_vector_set(proton_fraction_vector, i, proton_fraction);

        proton_fraction += proton_fraction_step;
    }
    printf("\n");


    // Write results
    SetFilePath("data");

    WriteVectorsToFile ("pressure.dat",
                        "# proton fraction, pressure at transition (MeV/fm^3)\n",
                        2,
                        proton_fraction_vector,
                        pressure_vector);

    WriteVectorsToFile ("barionic_chemical_potential.dat",
                        "# proton fraction, barionic chemical potential at transition (MeV)\n",
                        2,
                        proton_fraction_vector,
                        barionic_chemical_potential_vector);

    WriteVectorsToFile ("isovector_chemical_potential.dat",
                        "# proton fraction, barionic chemical potential at transition (MeV)\n",
                        2,
                        proton_fraction_vector,
                        isovector_chemical_potential_vector);

    if (options.verbose)
        printf("Done!\n");

    return 0;
}
