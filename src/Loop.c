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

static
int GridLoop(double initial_mass_guess,
             double initial_proton_density_guess,
             double initial_neutron_density_guess,
             double initial_up_mass_guess,
             double initial_down_mass_guess,
             double hadron_vacuum_thermodynamic_potential,
             double quark_vacuum_thermodynamic_potential,
             double min_barionic_chemical_potential,
             double min_isovector_chemical_potential,
             int barionic_chem_pot_num_pts,
             int isovector_chem_pot_num_pts,
             double barionic_chemical_potential_step,
             double isovector_chemical_potential_step,
             gsl_vector * hadron_pressure_grid_barionic_chem_pot,
             gsl_vector * hadron_pressure_grid_isovector_chem_pot,
             gsl_vector * hadron_pressure_grid_pressure,
             gsl_vector * quark_pressure_grid_barionic_chem_pot,
             gsl_vector * quark_pressure_grid_isovector_chem_pot,
             gsl_vector * quark_pressure_grid_pressure,
             gsl_vector * binodal_barionic_chem_pot,
             gsl_vector * binodal_isovector_chem_pot,
             gsl_vector * binodal_pressure,
             gsl_vector * binodal_hadron_asymmetry,
             gsl_vector * binodal_quark_asymmetry,
             int * count,
             int * num_successfull_pts);

int SolveBinodalForBarionicAndIsovectorChemicalPotentialsGrid()
{
    // Set path to write results
    SetFilePath("output/data");

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

    double hadron_vacuum_potential = HadronVacuumEnergyDensity();

    double up_vacuum_mass;
    double down_vacuum_mass;

    QuarkVacuumMassDetermination(&up_vacuum_mass, &down_vacuum_mass);

    double quark_vacuum_potential =
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
               quark_vacuum_potential);
        printf("\tHadron vacuum thermodynamic potential: %f\n",
               hadron_vacuum_potential);
    }

    if (options.verbose)
        printf("Solving gap equation and equations of state ...\n");


    double barionic_chemical_potential_step =
    Step(parameters.variables.min_barionic_chemical_potential,
         parameters.variables.max_barionic_chemical_potential,
         parameters.variables.num_points);

    double isovector_chemical_potential_step =
    Step(parameters.variables.min_isovector_chemical_potential,
         parameters.variables.max_isovector_chemical_potential,
         parameters.variables.num_points);

    HadronMassAndDensitiesSolutionParams ph =
    parameters.hadron.mass_and_densities_solution;

    QuarkMassAndRenormChemPotSolParams pq =
    parameters.quark.mass_and_renorm_chem_pot_solution;

    int num_grid_pts =
    parameters.variables.num_points * parameters.variables.num_points;

    int count;
    int num_successfull_pts;

    gsl_vector * hadron_pressure_grid_barionic_chem_pot =
        gsl_vector_alloc(num_grid_pts);
    gsl_vector * hadron_pressure_grid_isovector_chem_pot =
        gsl_vector_alloc(num_grid_pts);
    gsl_vector * hadron_pressure_grid_pressure =
        gsl_vector_alloc(num_grid_pts);

    gsl_vector * quark_pressure_grid_barionic_chem_pot =
        gsl_vector_alloc(num_grid_pts);
    gsl_vector * quark_pressure_grid_isovector_chem_pot =
        gsl_vector_alloc(num_grid_pts);
    gsl_vector * quark_pressure_grid_pressure =
        gsl_vector_alloc(num_grid_pts);

    gsl_vector * binodal_barionic_chem_pot = gsl_vector_alloc(num_grid_pts);
    gsl_vector * binodal_isovector_chem_pot = gsl_vector_alloc(num_grid_pts);
    gsl_vector * binodal_pressure = gsl_vector_alloc(num_grid_pts);
    gsl_vector * binodal_hadron_asymmetry = gsl_vector_alloc(num_grid_pts);
    gsl_vector * binodal_quark_asymmetry = gsl_vector_alloc(num_grid_pts);

    if (options.verbose)
            printf("Positive section:\n");

    int status =
    GridLoop(ph.initial_mass_guess,
             ph.initial_proton_density_guess,
             ph.initial_neutron_density_guess,
             pq.initial_up_mass_guess,
             pq.initial_down_mass_guess,
             hadron_vacuum_potential,
             quark_vacuum_potential,
             parameters.variables.min_barionic_chemical_potential,
             parameters.variables.min_isovector_chemical_potential,
             parameters.variables.num_points,
             parameters.variables.num_points,
             barionic_chemical_potential_step,
             isovector_chemical_potential_step,
             hadron_pressure_grid_barionic_chem_pot,
             hadron_pressure_grid_isovector_chem_pot,
             hadron_pressure_grid_pressure,
             quark_pressure_grid_barionic_chem_pot,
             quark_pressure_grid_isovector_chem_pot,
             quark_pressure_grid_pressure,
             binodal_barionic_chem_pot,
             binodal_isovector_chem_pot,
             binodal_pressure,
             binodal_hadron_asymmetry,
             binodal_quark_asymmetry,
             &count,
             &num_successfull_pts);

    if (status == 0){
        WriteVectorsToFileUpToIndex("binodal_positive_section.dat",
                                    "# isovector chemical potential (MeV), "
                                    "barionic chemical potential (MeV) \n",
                                    count,
                                    2,
                                    binodal_isovector_chem_pot,
                                    binodal_barionic_chem_pot);

        WriteVectorsToFileUpToIndex("pressure_vs_hadron_asymmetry_pos_sec.dat",
                                    "# hadron phase asymmetry, "
                                    "pressure at transition (MeV) \n",
                                    count,
                                    2,
                                    binodal_hadron_asymmetry,
                                    binodal_pressure);

        WriteVectorsToFileUpToIndex("pressure_vs_quark_asymmetry_pos_sec.dat",
                                    "# quark phase asymmetry, "
                                    "pressure at transition (MeV) \n",
                                    count,
                                    2,
                                    binodal_quark_asymmetry,
                                    binodal_pressure);

        WriteVectorsToFileUpToIndex("hadron_pressure_grid_positive_section.dat",
                                    "# barionic chemical potential (MeV), "
                                    "isovector chemical potential (MeV), "
                                    "hadron pressure (MeV/fm^3) \n",
                                    num_successfull_pts,
                                    3,
                                    hadron_pressure_grid_barionic_chem_pot,
                                    hadron_pressure_grid_isovector_chem_pot,
                                    hadron_pressure_grid_pressure);
                                   
        WriteVectorsToFileUpToIndex("quark_pressure_grid_positive_section.dat",
                                    "# barionic chemical potential (MeV), "
                                    "isovector chemical potential (MeV), "
                                    "quark pressure (MeV/fm^3) \n",
                                    num_successfull_pts,
                                    3,
                                    quark_pressure_grid_barionic_chem_pot,
                                    quark_pressure_grid_isovector_chem_pot,
                                    quark_pressure_grid_pressure);
    }

    if (options.verbose)
        printf("Negative section: \n");

    int status_g =
    GridLoop(ph.initial_mass_guess,
             ph.initial_proton_density_guess,
             ph.initial_neutron_density_guess,
             pq.initial_up_mass_guess,
             pq.initial_down_mass_guess,
             hadron_vacuum_potential,
             quark_vacuum_potential,
             parameters.variables.min_barionic_chemical_potential,
             parameters.variables.min_isovector_chemical_potential,
             parameters.variables.num_points,
             parameters.variables.num_points,
             barionic_chemical_potential_step,
             -isovector_chemical_potential_step,
             hadron_pressure_grid_barionic_chem_pot,
             hadron_pressure_grid_isovector_chem_pot,
             hadron_pressure_grid_pressure,
             quark_pressure_grid_barionic_chem_pot,
             quark_pressure_grid_isovector_chem_pot,
             quark_pressure_grid_pressure,
             binodal_barionic_chem_pot,
             binodal_isovector_chem_pot,
             binodal_pressure,
             binodal_hadron_asymmetry,
             binodal_quark_asymmetry,
             &count,
             &num_successfull_pts);

    if (status_g == 0){
        WriteVectorsToFileUpToIndex("binodal_negative_section.dat",
                                    "# isovector chemical potential (MeV), "
                                    "barionic chemical potential (MeV) \n",
                                    count,
                                    2,
                                    binodal_isovector_chem_pot,
                                    binodal_barionic_chem_pot);

        WriteVectorsToFileUpToIndex("pressure_vs_hadron_asymmetry_neg_sec.dat",
                                    "# hadron phase asymmetry, "
                                    "pressure at transition (MeV) \n",
                                    count,
                                    2,
                                    binodal_hadron_asymmetry,
                                    binodal_pressure);

        WriteVectorsToFileUpToIndex("pressure_vs_quark_asymmetry_neg_sec.dat",
                                    "# quark phase asymmetry, "
                                    "pressure at transition (MeV) \n",
                                    count,
                                    2,
                                    binodal_quark_asymmetry,
                                    binodal_pressure);

        WriteVectorsToFileUpToIndex("hadron_pressure_grid_negative_section.dat",
                                    "# barionic chemical potential (MeV), "
                                    "isovector chemical potential (MeV), "
                                    "hadron pressure (MeV/fm^3) \n",
                                    num_successfull_pts,
                                    3,
                                    hadron_pressure_grid_barionic_chem_pot,
                                    hadron_pressure_grid_isovector_chem_pot,
                                    hadron_pressure_grid_pressure);
                                   
        WriteVectorsToFileUpToIndex("quark_pressure_grid_negative_section.dat",
                                    "# barionic chemical potential (MeV), "
                                    "isovector chemical potential (MeV), "
                                    "quark pressure (MeV/fm^3) \n",
                                    num_successfull_pts,
                                    3,
                                    quark_pressure_grid_barionic_chem_pot,
                                    quark_pressure_grid_isovector_chem_pot,
                                    quark_pressure_grid_pressure);
    }

    return 0;
}

static
int GridLoop(double initial_mass_guess,
             double initial_proton_density_guess,
             double initial_neutron_density_guess,
             double initial_up_mass_guess,
             double initial_down_mass_guess,
             double hadron_vacuum_thermodynamic_potential,
             double quark_vacuum_thermodynamic_potential,
             double min_barionic_chemical_potential,
             double min_isovector_chemical_potential,
             int barionic_chem_pot_num_pts,
             int isovector_chem_pot_num_pts,
             double barionic_chemical_potential_step,
             double isovector_chemical_potential_step,
             gsl_vector * hadron_pressure_grid_barionic_chem_pot,
             gsl_vector * hadron_pressure_grid_isovector_chem_pot,
             gsl_vector * hadron_pressure_grid_pressure,
             gsl_vector * quark_pressure_grid_barionic_chem_pot,
             gsl_vector * quark_pressure_grid_isovector_chem_pot,
             gsl_vector * quark_pressure_grid_pressure,
             gsl_vector * binodal_barionic_chem_pot,
             gsl_vector * binodal_isovector_chem_pot,
             gsl_vector * binodal_pressure,
             gsl_vector * binodal_hadron_asymmetry,
             gsl_vector * binodal_quark_asymmetry,
             int * binodal_count,
             int * num_successfull_pts)
{
    double hadron_mass_scanline_initial_guess = initial_mass_guess;
    double proton_dens_scanline_initial_guess = initial_proton_density_guess;
    double neutron_dens_scanline_initial_guess = initial_neutron_density_guess;

    double up_mass_scanline_initial_guess = initial_up_mass_guess;
    double down_mass_scanline_initial_guess = initial_down_mass_guess;

    double isovector_chemical_potential = min_isovector_chemical_potential;
    
    *num_successfull_pts = 0;
    *binodal_count = 0;
    for (int i = 0; i < isovector_chem_pot_num_pts; i++){

        double barionic_chemical_potential = min_barionic_chemical_potential;

        double hadron_mass_guess = hadron_mass_scanline_initial_guess;
        double proton_density_guess = proton_dens_scanline_initial_guess;
        double neutron_density_guess = neutron_dens_scanline_initial_guess;

        double up_mass_guess = up_mass_scanline_initial_guess;
        double down_mass_guess = down_mass_scanline_initial_guess;

        for (int j = 0; j < barionic_chem_pot_num_pts; j++){

            if (options.verbose){
                printf("\r\t%3.0f%%",
                       (j + 1 + i * barionic_chem_pot_num_pts) * 100
                       / pow(barionic_chem_pot_num_pts, 2.0));
                fflush(stdout);
            }

            BinodalPoint point;

            int status =
            BinodalPointCandidate(barionic_chemical_potential,
                                  isovector_chemical_potential,
                                  parameters.variables.temperature,
                                  hadron_vacuum_thermodynamic_potential,
                                  hadron_mass_guess,
                                  proton_density_guess,
                                  neutron_density_guess,
                                  quark_vacuum_thermodynamic_potential,
                                  up_mass_guess,
                                  down_mass_guess,
                                  &point);

            // When BinodalPointCandidate returns zero, we have a successfull
            // calculation. In this case, we will test if the pressure tolerance
            // is met. If it's not, the point is skipped. If the calculation is
            // not successfull, the point is also skipped.
            if (!status){
            
                // save values of quark and hadron pressures
                gsl_vector_set(hadron_pressure_grid_barionic_chem_pot,
                               *num_successfull_pts,
                               barionic_chemical_potential);
                               
                gsl_vector_set(hadron_pressure_grid_isovector_chem_pot,
                               *num_successfull_pts,
                               isovector_chemical_potential);
                               
                gsl_vector_set(hadron_pressure_grid_pressure,
                               *num_successfull_pts,
                               point.hadron_pressure);
                               
                gsl_vector_set(quark_pressure_grid_barionic_chem_pot,
                               *num_successfull_pts,
                               barionic_chemical_potential);
                               
                gsl_vector_set(quark_pressure_grid_isovector_chem_pot,
                               *num_successfull_pts,
                               isovector_chemical_potential);
                               
                gsl_vector_set(quark_pressure_grid_pressure,
                               *num_successfull_pts,
                               point.quark_pressure);
                               
                (*num_successfull_pts)++;
                               
                               
                if (fabs(point.quark_pressure - point.hadron_pressure)
                    < parameters.variables.pressure_tolerance){

                    gsl_vector_set(binodal_barionic_chem_pot,
                                   *binodal_count,
                                   barionic_chemical_potential);
                    gsl_vector_set(binodal_isovector_chem_pot,
                                   *binodal_count,
                                   isovector_chemical_potential);
                    gsl_vector_set(binodal_pressure,
                                   *binodal_count,
                                   point.hadron_pressure);

                    gsl_vector_set(binodal_hadron_asymmetry,
                                   *binodal_count,
                                   HadronPhaseAsymmetry(point.proton_density,
                                                        point.neutron_density));

                    gsl_vector_set(binodal_quark_asymmetry,
                                   *binodal_count,
                                   QuarkPhaseAsymmetry(point.up_quark_density,
                                                       point.down_quark_density));

                    (*binodal_count)++;
                }
            }
            else{
                printf("\n%s:%d: Determination of binodal point candidate "
                       "unsuccessfull, skipping point.\n",
                       __FILE__,
                       __LINE__);
            }


            // Update guesses for next point in the scanline.
            // The masses are not supposed to have bigger values than
            // the initial guess provided to the function. If this happens,
            // reject the guess (maybe even the solution should be discarded
            // as I see no reason for that to happen). The update is performed
            // only if we have a successfull binodal point candidate calculation
            if (!status){

                hadron_mass_guess = point.hadron_mass;
                proton_density_guess = point.proton_density;
                neutron_density_guess = point.neutron_density;
                down_mass_guess = point.down_quark_mass;
                up_mass_guess = point.up_quark_mass;

                // Save first point results to use as initial guesses
                // for the neighbouring scanline
                if (j == 0){
                    hadron_mass_scanline_initial_guess = hadron_mass_guess;
                    proton_dens_scanline_initial_guess = proton_density_guess;
                    neutron_dens_scanline_initial_guess = neutron_density_guess;

                    up_mass_scanline_initial_guess = up_mass_guess;
                    down_mass_scanline_initial_guess = down_mass_guess;
                }
            }

            barionic_chemical_potential += barionic_chemical_potential_step;
        }

        isovector_chemical_potential += isovector_chemical_potential_step;
    }

    if (options.verbose)
        printf("\n");

    return 0;
}

