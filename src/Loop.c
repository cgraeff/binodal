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
               "\ttemperature: %4.2f (MeV)\n",
               parameters.quark_model.parameters_set_identifier,
               parameters.hadron_model.parameters_set_identifier,
               parameters.variables.temperature);
    }

    // Vacuum mass determination
    if (options.verbose)
        printf("Determining the vacuum mass and bag constant ...\n");

    double up_vacuum_mass;
    double down_vacuum_mass;
    QuarkVacuumMassDetermination(&up_vacuum_mass,
                                 &down_vacuum_mass);

    double quark_vacuum_thermodynamic_potential =
	    QuarkThermodynamicPotential(up_vacuum_mass,
                                    down_vacuum_mass,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0);

    // TODO: should the hadron mass be calculated?
    double hadron_vacuum_mass = parameters.hadron_model.nucleon_mass;

	double hadron_vacuum_thermodynamic_potential = HadronVacuumEnergyDensity();

    if (options.verbose){
        printf("\tUp quark vacuum mass: %f\n", up_vacuum_mass);
        printf("\tDown quark vacuum mass: %f\n", down_vacuum_mass);
        printf("\tHadron vacuum mass: %f\n", hadron_vacuum_mass);
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
    gsl_vector * hadron_mass_vector =
        gsl_vector_alloc(parameters.variables.num_points);
    gsl_vector * up_quark_mass_vector =
        gsl_vector_alloc(parameters.variables.num_points);
    gsl_vector * down_quark_mass_vector =
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
            printf("\tProton fraction: %f\n", proton_fraction);

        double barionic_density;
        double hadron_mass;
        double up_quark_mass;
        double down_quark_mass;

        // FIXME: We always start with the initial guess for the first point.
        //        Make it use the last result as initial guess from the second
        //        point on.
        SimultaneousSolution(parameters.simultaneous_solution,
                             quark_vacuum_thermodynamic_potential,
                             hadron_vacuum_thermodynamic_potential,
                             proton_fraction,
                             &barionic_density,
                             &hadron_mass,
                             &up_quark_mass,
                             &down_quark_mass);

        gsl_vector_set(barionic_density_vector, i, barionic_density);
        gsl_vector_set(hadron_mass_vector, i, hadron_mass);
        gsl_vector_set(up_quark_mass_vector, i, up_quark_mass);
        gsl_vector_set(down_quark_mass_vector, i, down_quark_mass);

        // Hadrons:
        double proton_density = proton_fraction * barionic_density;
        double neutron_density = (1.0 - proton_fraction) * barionic_density;

        double proton_fermi_momentum = HPFermiMomentum(proton_density);
        double neutron_fermi_momentum = HPFermiMomentum(neutron_density);

        double total_hadron_scalar_density =
            HadronScalarDensity(hadron_mass,
                                proton_fermi_momentum,
                                parameters.hadron_model.cutoff)
            + HadronScalarDensity(hadron_mass,
                                  neutron_fermi_momentum,
                                  parameters.hadron_model.cutoff);

        double proton_chemical_potential =
            ProtonChemicalPotential(proton_fermi_momentum,
                                    total_hadron_scalar_density,
                                    hadron_mass,
                                    barionic_density,
                                    proton_density,
                                    neutron_density);

        double neutron_chemical_potential =
            NeutronChemicalPotential(neutron_fermi_momentum,
                                     total_hadron_scalar_density,
                                     hadron_mass,
                                     barionic_density,
                                     proton_density,
                                     neutron_density);

        // Barionic (isoscalar) and isovector chemical potentials
        double barionic_chemical_potential = (proton_chemical_potential
                                              + neutron_chemical_potential)
                                              / 2.0;

        double isovector_chemical_potential = proton_chemical_potential
                                              - neutron_chemical_potential;

        gsl_vector_set(barionic_chemical_potential_vector,
                       i,
                       barionic_chemical_potential);
        gsl_vector_set(isovector_chemical_potential_vector,
                       i,
                       isovector_chemical_potential);

        double hadron_kinectic_energy_density =
            HadronKinecticEnergyDensity(hadron_mass,
                                        proton_fermi_momentum,
                                        neutron_fermi_momentum);

        double hadron_thermodynamic_potential =
            HadronThermodynamicPotential(total_hadron_scalar_density,
                                         barionic_density,
                                         proton_density,
                                         neutron_density,
                                         proton_chemical_potential,
                                         neutron_chemical_potential,
                                         hadron_kinectic_energy_density);

        double hadron_pressure = HadronPressure(hadron_thermodynamic_potential);

        gsl_vector_set(pressure_vector, i, hadron_pressure);

        gsl_vector_set(proton_fraction_vector, i, proton_fraction);
        proton_fraction += proton_fraction_step;
    } // LOOP ENDS
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
