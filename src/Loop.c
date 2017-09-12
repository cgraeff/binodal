//
//  Loop.c
//  binodal
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

int SolveBinodalForBarionicChemicalPotentialRange(){

    // Print name of parametrization
    if (options.verbose){
        printf("Calculation performed with %s quark phase parameters set\n"
               "and %s hadron phase paarameters set.\n"
               "\ttemperature: %4.2f (MeV)\n",
               parameters.quark.model.parameters_set_identifier,
               parameters.hadron.model.parameters_set_identifier,
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

    gsl_vector * hadron_proton_fraction_vector =
    gsl_vector_alloc(parameters.variables.num_points);

    gsl_vector * quark_proton_fraction_vector =
    gsl_vector_alloc(parameters.variables.num_points);

    gsl_vector * hadron_asymmetry_vector =
    gsl_vector_alloc(parameters.variables.num_points);

    gsl_vector * quark_asymmetry_vector =
    gsl_vector_alloc(parameters.variables.num_points);

    gsl_vector * barionic_density_vector =
    gsl_vector_alloc(parameters.variables.num_points);

    gsl_vector * barionic_chemical_potential_vector =
    gsl_vector_alloc(parameters.variables.num_points);

    gsl_vector * isovector_chemical_potential_vector =
    gsl_vector_alloc(parameters.variables.num_points);

    gsl_vector * pressure_vector =
    gsl_vector_alloc(parameters.variables.num_points);

    double isovector_chemical_potential =
    parameters.variables.min_isovector_chemical_potential;

    double isovector_chemical_potential_step =
    Step(parameters.variables.min_isovector_chemical_potential,
         parameters.variables.max_isovector_chemical_potential,
         parameters.variables.num_points);

    for (int i = 0; i < parameters.variables.num_points; i++){

        if (options.verbose)
            printf("\tIsovector chemical potential: %f\r",
                   isovector_chemical_potential);

        BinodalPoint point =
        DetermineBinodalPoint(parameters.variables.temperature,
                              isovector_chemical_potential,
                              hadron_vacuum_thermodynamic_potential,
                              quark_vacuum_thermodynamic_potential);

        double quark_proton_fraction =
        QuarkProtonFraction(point.up_quark_density,
                            point.down_quark_density);

        double barionic_density = point.proton_density + point.neutron_density;
        double barionic_chemical_potential =
        BarionicChemicalPotential(point.proton_density,
                                  point.neutron_chemical_potential);

        gsl_vector_set(hadron_asymmetry_vector,
                       i,
                       HadronPhaseAsymmetry(point.proton_density,
                                            point.neutron_density));

        gsl_vector_set(quark_asymmetry_vector,
                       i,
                       QuarkPhaseAsymmetry(point.up_quark_density,
                                           point.down_quark_density));

        gsl_vector_set(barionic_chemical_potential_vector,
                       i,
                       barionic_chemical_potential);
        gsl_vector_set(isovector_chemical_potential_vector,
                       i,
                       isovector_chemical_potential);

        gsl_vector_set(barionic_density_vector,
                       i,
                       barionic_density);
        gsl_vector_set(pressure_vector, i, point.pressure);
        gsl_vector_set(hadron_proton_fraction_vector,
                       i,
                       point.proton_density / barionic_density);
        gsl_vector_set(quark_proton_fraction_vector,
                       i,
                       quark_proton_fraction);

        isovector_chemical_potential += isovector_chemical_potential_step;
    }
    printf("\n");


    // Write results
    SetFilePath("output/data");

    WriteVectorsToFile("pressure_vs_isovector_chemical_potential.dat",
                       "# isovector chemical potential (MeV), "
                       "pressure at transition (MeV/fm^3)\n",
                       2,
                       isovector_chemical_potential_vector,
                       pressure_vector);

    WriteVectorsToFile("pressure_vs_hadron_proton_fraction.dat",
                       "# hadron proton fraction, "
                       "pressure at transition (MeV/fm^3) \n",
                       2,
                       hadron_proton_fraction_vector,
                       pressure_vector);

    WriteVectorsToFile("pressure_vs_quark_proton_fraction.dat",
                       "# hadron proton fraction, "
                       "pressure at transition (MeV/fm^3) \n",
                       2,
                       quark_proton_fraction_vector,
                       pressure_vector);

    WriteVectorsToFile("pressure_vs_hadron_asymmetry.dat",
                       "# asymmetry, pressure at transition (MeV/fm^3) \n",
                       2,
                       hadron_asymmetry_vector,
                       pressure_vector);

    WriteVectorsToFile("pressure_vs_quark_asymmetry.dat",
                       "# asymmetry, pressure at transition (MeV/fm^3) \n",
                       2,
                       quark_asymmetry_vector,
                       pressure_vector);

    WriteVectorsToFile("barionic_density_vs_isovector_chemical_potential.dat",
                       "# isovector chemical potential (MeV), "
                       "barionic density at transition (fm^{-3})\n",
                       2,
                       isovector_chemical_potential_vector,
                       barionic_density_vector);

    WriteVectorsToFile("barionic_chemical_potential_vs_isovector_chemical_potential.dat",
                       "# isovector chemical potential (MeV), "
                       "barionic chemical potential at transition (MeV)\n",
                       2,
                       isovector_chemical_potential_vector,
                       barionic_chemical_potential_vector);

    WriteVectorsToFile("hadron_proton_fraction_vs_isovector_chemical_potential.dat",
                       "# isovector chemical potential (MeV), "
                       "hadron proton fraction at transition (MeV)\n",
                       2,
                       isovector_chemical_potential_vector,
                       hadron_proton_fraction_vector);

    WriteVectorsToFile("quark_proton_fraction_vs_isovector_chemical_potential.dat",
                       "# isovector chemical potential (MeV), "
                       "quark proton fraction at transition (MeV)\n",
                       2,
                       isovector_chemical_potential_vector,
                       quark_proton_fraction_vector);


    if (options.verbose)
        printf("Done!\n");

    return 0;
}

int SolveBinodalForBarionicAndIsovectorChemicalPotentialGrid()
{
    // Set path to write results
    SetFilePath("output/data");

    FILE * file_h = OpenFile("hadron_pressure.dat");
    FILE * file_q = OpenFile("quark_pressure.dat");
    FILE * file_t = OpenFile("transition.dat");

    // Print file headers
    fprintf(file_h,
            "# barionic_chemical_potential (MeV), "
            "isovector chemical potential (MeV), "
            "hadron pressure (MeV/fm^3)\n");
    fprintf(file_q,
            "# barionic_chemical_potential (MeV), "
            "isovector chemical potential (MeV), "
            "quark pressure (MeV/fm^3)\n");
    fprintf(file_t,
            "# barionic_chemical_potential (MeV), "
            "isovector chemical potential (MeV), "
            "transition pressure pressure (MeV/fm^3)\n");

    double hadron_vacuum_potential = HadronVacuumEnergyDensity();

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

    double barionic_chemical_potential_step =
    Step(parameters.variables.min_barionic_chemical_potential,
         parameters.variables.max_barionic_chemical_potential,
         parameters.variables.num_points);

    double isovector_chemical_potential_step =
    Step(parameters.variables.min_isovector_chemical_potential,
         parameters.variables.max_isovector_chemical_potential,
         parameters.variables.num_points);

    double isovector_chemical_potential =
    parameters.variables.min_isovector_chemical_potential;

    HadronMassAndDensitiesSolutionParams ph =
    parameters.hadron.mass_and_densities_solution;

    QuarkMassAndRenormChemPotSolParams pq =
    parameters.quark.mass_and_renorm_chem_pot_solution;

    double hadron_mass_initial_guess = ph.initial_mass_guess;
    double proton_density_initial_guess = ph.initial_proton_density_guess;
    double neutron_density_initial_guess = ph.initial_neutron_density_guess;

    double up_mass_initial_guess = pq.initial_up_mass_guess;
    double down_mass_initial_guess = pq.initial_down_mass_guess;

    for (int i = 0; i < parameters.variables.num_points; i++){

        double barionic_chemical_potential =
        parameters.variables.min_barionic_chemical_potential;

        double hadron_mass_guess = hadron_mass_initial_guess;
        double proton_density_guess = proton_density_initial_guess;
        double neutron_density_guess = neutron_density_initial_guess;

        double up_mass_guess = up_mass_initial_guess;
        double down_mass_guess = down_mass_initial_guess;

        for (int j = 0; j < parameters.variables.num_points; j++){

            double proton_chemical_potential =
            ProtonChemicalPotential(barionic_chemical_potential,
                                    isovector_chemical_potential);

            double neutron_chemical_potential =
            NeutronChemicalPotential(barionic_chemical_potential,
                                     isovector_chemical_potential);

            double hadron_mass = NAN;
            double hadron_pressure = NAN;
            double proton_density = NAN;
            double neutron_density = NAN;
            DetermineHadronPressureAndDensities(proton_chemical_potential,
                                                neutron_chemical_potential,
                                                hadron_vacuum_potential,
                                                hadron_mass_guess,
                                                proton_density_guess,
                                                neutron_density_guess,
                                                &hadron_mass,
                                                &proton_density,
                                                &neutron_density,
                                                &hadron_pressure);

            double up_chemical_potential =
            UpChemicalPotentialFromGibbsConditions(proton_chemical_potential,
                                                   neutron_chemical_potential);

            double down_chemical_potential =
            DownChemicalPotentialFromGibbsConditions(proton_chemical_potential,
                                                     neutron_chemical_potential);

            double up_quark_mass;
            double down_quark_mass;
            double quark_pressure;

            DetermineQuarkPressure(up_chemical_potential,
                                   down_chemical_potential,
                                   parameters.variables.temperature,
                                   quark_vacuum_thermodynamic_potential,
                                   up_mass_guess,
                                   down_mass_guess,
                                   &up_quark_mass,
                                   &down_quark_mass,
                                   &quark_pressure);

            fprintf(file_h,
                    "%20.15E\t%20.15E\t%20.15E\n",
                    barionic_chemical_potential,
                    isovector_chemical_potential,
                    hadron_pressure);

            fprintf(file_q,
                    "%20.15E\t%20.15E\t%20.15E\n",
                    barionic_chemical_potential,
                    isovector_chemical_potential,
                    quark_pressure);

            if (fabs(quark_pressure - hadron_pressure)
                < parameters.variables.pressure_tolerance){
                fprintf(file_t,
                        "%20.15E\t%20.15E\t%20.15E\n",
                        barionic_chemical_potential,
                        isovector_chemical_potential,
                        hadron_pressure);
            }

            // Update guesses for next point in the scanline
            hadron_mass_guess = hadron_mass;
            proton_density_guess = proton_density;
            neutron_density_guess = neutron_density;

            up_mass_guess = up_quark_mass;
            down_mass_guess = down_quark_mass;

            // Store solution of first point for use as initial guess in the
            // next scanline
            if (j == 0){
                hadron_mass_initial_guess = hadron_mass;
                proton_density_initial_guess = proton_density;
                neutron_density_initial_guess = neutron_density;

                up_mass_initial_guess = up_quark_mass;
                down_mass_initial_guess = down_quark_mass;
            }

            barionic_chemical_potential += barionic_chemical_potential_step;
        }
        fprintf(file_h, "\n");
        fprintf(file_q, "\n");

        isovector_chemical_potential += isovector_chemical_potential_step;
    }

    fclose(file_h);
    fclose(file_q);

    return 0;
}
