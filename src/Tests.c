//
//  Tests.c
//  binodal
//
//  Created by Clebson Graeff on 2016-06-08.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <stdbool.h>

#include "libdatafun/libdatafun.h"

#include "Tests.h"

#include "Parameters.h"
#include "QuarkPhaseEOS.h"
#include "HadronPhaseEOS.h"
#include "Binodal.h"

void RunTests()
{
    SetParametersSet(NULL, NULL);

    // Determine thermodynamic potential map as function of quark masses:
    //      To plot on gnuplot use
    //          splot "quark_thermodynamic_potential_mass_map.dat" w pm3d
    //      To make a top view with contour lines use
    //          set view map
    //          set contour
    //          set cntrparam levels 15
    //          splot "quark_thermodynamic_potential_mass_map.dat" w pm3d
    if(false)
    {
        SetParametersSet("Buballa_1", NULL);
        SetFilePath("tests/quark_thermodynamic_potential/data");

        double min_mass = 0.0;
        double max_mass = 1000.0;
        int n_pts = 100;

        double up_chemical_potential = 0.0;
        double down_chemical_potential = 0.0;
        double up_renormalized_chemical_potential = 0.0;
        double down_renormalized_chemical_potential = 0.0;
        double temperature = 0.0;

        FILE * file = OpenFile("quark_thermodynamic_potential_mass_map.dat");
        fprintf(file,
                "# up mass (MeV), "
                "down mass (MeV), "
                "quark thermodynamic potential (MeV/fm^3)\n");

        double mass_step = Step (min_mass, max_mass, n_pts);

        double down_mass = 0.0;
        for (int i = 0; i < n_pts; i++){

            double up_mass = 0.0;

            for (int j = 0; j < n_pts; j++){

                double omega =
                QuarkThermodynamicPotential(up_mass,
                                            down_mass,
                                            up_chemical_potential,
                                            down_chemical_potential,
                                            up_renormalized_chemical_potential,
                                            down_renormalized_chemical_potential,
                                            temperature);

                fprintf(file,
                        "%20.15E\t%20.15E\t%20.15E\t\n",
                        up_mass,
                        down_mass,
                        omega);

                up_mass += mass_step;
            }

            fprintf(file, "\n"); // This is needed to tell gnuplot that we are
                                 // going to scan anoter line of the surface

            down_mass += mass_step;
        }

        fclose(file);

        // Find minimum
        double up_mass_at_minimum;
        double down_mass_at_minimum;

        QuarkThermodynamicPotentialMinimum(&up_mass_at_minimum,
                                           &down_mass_at_minimum);

        double potential_at_minimum =
        QuarkThermodynamicPotential(up_mass_at_minimum,
                                    down_mass_at_minimum,
                                    up_chemical_potential,
                                    down_chemical_potential,
                                    up_renormalized_chemical_potential,
                                    down_renormalized_chemical_potential,
                                    temperature);

        file =
        OpenFile("quark_thermodynamic_potential_mass_map_minimum_location.dat");
        fprintf(file,
                "# up mass (MeV), "
                "down mass (MeV), "
                "quark thermodynamic potential (MeV/fm^3)\n");

        fprintf(file,
                "%20.15E\t%20.15E\t%20.15E\n",
                up_mass_at_minimum,
                down_mass_at_minimum,
                potential_at_minimum);

        fclose(file);

        SetParametersSet(NULL, NULL);
    }

    // Determine quark pressure as functions of quark chemical potentials.
    // The same functions used to determine the binodal are used here.
    //      To plot on gnuplot use
    //          splot "quark_pressure_chem_pot_map.dat" w pm3d
    //      To make a top view with contour lines use
    //          set view map
    //          set contour
    //          set cntrparam levels 15
    //          splot "quark_pressure_map.dat" w pm3d
    if(false)
    {
        SetParametersSet("Buballa_1", NULL);
        SetFilePath("tests/quark_pressure/data");

        double min_up_chemical_potential = 0.0;
        double max_up_chemical_potential = 1000.0;
        double min_down_chemical_potential = 0.0;
        double max_down_chemical_potential = 1000.0;

        int n_pts = 100;

        double temperature = 0.0;

        double initial_up_mass_guess = 313.0;
        double initial_down_mass_guess = 313.0;

        double up_vacuum_mass;
        double down_vacuum_mass;

        QuarkVacuumMassDetermination (&up_vacuum_mass, &down_vacuum_mass);

        double quark_vacuum_thermodynamic_potential =
        QuarkThermodynamicPotential(up_vacuum_mass,
                                    down_vacuum_mass,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0);

        FILE * file =
        OpenFile("quark_pressure_chem_pot_map.dat");
        fprintf (file,
                 "# up chemical potential (MeV), "
                 "down chemical potential (MeV), "
                 "quark pressure (MeV/fm^3)\n");

        double up_chem_pot_step = Step(min_up_chemical_potential,
                                       max_up_chemical_potential,
                                       n_pts);
        double down_chem_pot_step = Step(min_down_chemical_potential,
                                         max_down_chemical_potential,
                                         n_pts);

        double up_mass_guess = initial_up_mass_guess;
        double down_mass_guess = initial_down_mass_guess;

        double up_chemical_potential = 0.0;
        for (int i = 0; i < n_pts; i++){

            double down_chemical_potential = 0.0;

            for (int j = 0; j < n_pts; j++){

                double up_quark_mass;
                double down_quark_mass;
                double quark_pressure;

                DetermineQuarkPressure(up_chemical_potential,
                                       down_chemical_potential,
                                       temperature,
                                       quark_vacuum_thermodynamic_potential,
                                       up_mass_guess,
                                       down_mass_guess,
                                       &up_quark_mass,
                                       &down_quark_mass,
                                       &quark_pressure);

                fprintf(file,
                        "%20.15E\t%20.15E\t%20.15E\n",
                        up_chemical_potential,
                        down_chemical_potential,
                        quark_pressure);

                // Update guesses
                up_mass_guess = up_quark_mass;
                down_mass_guess = down_quark_mass;

                down_chemical_potential += down_chem_pot_step;
            }
            fprintf(file, "\n");

            up_chemical_potential += up_chem_pot_step;
        }

        fclose(file);

        SetParametersSet(NULL, NULL);
    }

    // Determine pressures for hadron and quark phases for a
    // particular value of the isovector chemical potential
    if(false)
    {
        SetParametersSet("PCP-0.2", "eNJL2mSigmaRho1");
        SetFilePath("tests/binodal_point_graph/data");

        FILE * file_h = OpenFile("hadron_pressure.dat");
        FILE * file_q = OpenFile("quark_pressure.dat");

        // printf file headers
        fprintf(file_h,
                "# barionic_chemical_potential (MeV), "
                "hadron pressure (MeV/fm^3)\n");
        fprintf(file_q,
                "# barionic_chemical_potential (MeV), "
                "quark pressure (MeV/fm^3)\n");

        int points_number = 3000;

        double min_barionic_chemical_potential = 950.0;
        double max_barionic_chemical_potential = 2000.0;

        double isovector_chemical_potential = 0.0;
        double temperature = 0.0;

        double initial_hadron_mass_guess = 1000.0;
        double initial_proton_density_guess = 0.05;
        double initial_neutron_density_guess = 0.05;

        double initial_up_mass_guess = 313.0;
        double initial_down_mass_guess = 313.0;

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
        Step(min_barionic_chemical_potential,
             max_barionic_chemical_potential,
             points_number);

        double hadron_mass_guess = initial_hadron_mass_guess;
        double proton_density_guess = initial_proton_density_guess;
        double neutron_density_guess = initial_neutron_density_guess;

        double up_mass_guess = initial_up_mass_guess;
        double down_mass_guess = initial_down_mass_guess;

        double barionic_chemical_potential = min_barionic_chemical_potential;
        for (int i = 0; i < points_number; i++){

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
                                   temperature,
                                   quark_vacuum_thermodynamic_potential,
                                   up_mass_guess,
                                   down_mass_guess,
                                   &up_quark_mass,
                                   &down_quark_mass,
                                   &quark_pressure);

            fprintf(file_h,
                    "%20.15E\t%20.15E\n",
                    barionic_chemical_potential,
                    hadron_pressure);

            fprintf(file_q,
                    "%20.15E\t%20.15E\n",
                    barionic_chemical_potential,
                    quark_pressure);

            // Update guesses
            hadron_mass_guess = hadron_mass;
            proton_density_guess = proton_density;
            neutron_density_guess = neutron_density;

            up_mass_guess = up_quark_mass;
            down_mass_guess = down_quark_mass;

            barionic_chemical_potential += barionic_chemical_potential_step;
        }

        fclose(file_h);
        fclose(file_q);

        SetParametersSet(NULL, NULL);
    }

    // Map zeroed gap equations for vacuum
    //      To plot on gnuplot use
    //          splot "quark_vacuum_mass_map.dat" w pm3d
    //      To make a top view with contour lines use
    //          set view map
    //          set contour
    //          set cntrparam levels 15
    //          splot "quark_vacuum_mass_map.dat" w pm3d
    if (false)
    {
        SetParametersSet("PCP-0.0", NULL);
        SetFilePath("tests/quark_vacuum_mass_maps/data");

        double min_mass = 0.0;
        double max_mass = 400.0;
        int n_pts = 100;

        FILE * file = OpenFile("quark_vacuum_mass_maps.dat");
        fprintf(file,
                "# up mass (MeV), down mass (MeV), zeroed gap equation (MeV)\n");

        // Save a side view of f(m_u, m_d) with m_u = m_d
        FILE * equal_mass_view = OpenFile("equal_mass_view.dat");
        fprintf(equal_mass_view,
                "# mass (MeV), zeroed gap equation (MeV)\n");

        FILE * intersection_file = OpenFile ("intersection.dat");

        double mass_step = Step (min_mass, max_mass, n_pts);

        gsl_vector * x = gsl_vector_alloc(2);
        gsl_vector * return_values = gsl_vector_alloc(2);

        double down_mass = 0.0;
        for (int i = 0; i < n_pts; i++){

            // the square root is to account for the
            // mapping done on QuarkVacuumMassDeterminationEquation
            // below
            gsl_vector_set(x, 0, pow(down_mass, 1.0 / 2.0));

            double up_mass = 0.0;

            for (int j = 0; j < n_pts; j++){

                gsl_vector_set(x, 1, pow(up_mass, 1.0 / 2.0)); // pow: same as above

                QuarkVacuumMassDeterminationEquation(x,
                                                     NULL,
                                                     return_values);

                if (gsl_vector_get(return_values, 0)
                    == gsl_vector_get(return_values, 1)){

                        fprintf(intersection_file,
                                "%20.15E\t%20.15E\t%20.15E\n",
                                up_mass,
                                down_mass,
                                gsl_vector_get(return_values, 0));
                }

                fprintf(file,
                        "%20.15E\t%20.15E\t%20.15E\t%20.15E\n",
                        up_mass,
                        down_mass,
                        gsl_vector_get(return_values, 0),
                        gsl_vector_get(return_values, 1));

                if (up_mass == down_mass){
                    fprintf(equal_mass_view,
                            "%20.15E\t%20.15E\n",
                            up_mass,
                            gsl_vector_get(return_values, 0));
                }

                up_mass += mass_step;
            }
            fprintf(file, "\n");

            down_mass += mass_step;
        }

        gsl_vector_free(x);
        gsl_vector_free(return_values);

        fclose(file);
        fclose(equal_mass_view);

        // Find common zero
        double up_vacuum_mass;
        double down_vacuum_mass;

        QuarkVacuumMassDetermination(&up_vacuum_mass, &down_vacuum_mass);

        file = OpenFile("quark_vacuum_mass_location.dat");

        fprintf (file,
                 "%20.15E\t%20.15E\t%20.15E\n",
                 up_vacuum_mass,
                 down_vacuum_mass,
                 0.0);

        fclose(file);

        SetParametersSet(NULL, NULL);
    }

    // Maps for determination of quark masses
    // and renormalized chemical potentials
    if (false)
    {
        SetParametersSet("PCP-0.0", NULL);
        SetFilePath("tests/quark_mass_maps/data");

        int n_pts = 100;

        double min_up_mass = 0.0;
        double max_up_mass = 600.0;
        double min_down_mass = 0.0;
        double max_down_mass = 600.0;

        double temperature = 0.0;

        double up_chemical_potential = 4.068543850286816E+02;
        double down_chemical_potential = 3.821859585511042E+02;

        FILE * output = OpenFile("maps.dat");
        FILE * line = OpenFile("line.dat");

        double up_mass_step = Step (min_up_mass, max_up_mass, n_pts);
        double down_mass_step = Step (min_down_mass, max_down_mass, n_pts);

        double down_mass = min_down_mass;

        for (int i = 0; i < n_pts; i++){

            double up_mass = min_up_mass;

            for (int j = 0; j < n_pts; j++){

                double up_renorm_chem_pot;
                double down_renorm_chem_pot;

                // for given parameters, determine the renormalized
                // chemical potentials
                QuarkSelfConsistentRenormChemPot(up_mass,
                                                 down_mass,
                                                 up_chemical_potential,
                                                 down_chemical_potential,
                                                 temperature,
                                                 &up_renorm_chem_pot,
                                                 &down_renorm_chem_pot);

                // Gap equations:
                double up_scalar_density =
                    QuarkScalarDensity(parameters.variables.temperature,
                                       up_mass,
                                       up_renorm_chem_pot);

                double down_scalar_density =
                    QuarkScalarDensity(parameters.variables.temperature,
                                       down_mass,
                                       down_renorm_chem_pot);

                double up_quark_zeroed_gap_eq =
                    QuarkZeroedGapEquation(up_mass,
                                           up_scalar_density,
                                           down_scalar_density);

                double down_quark_zeroed_gap_eq =
                    QuarkZeroedGapEquation(down_mass,
                                           up_scalar_density,
                                           down_scalar_density);

                fprintf(output,
                        "%20.15E\t%20.15E\t%20.15E\t%20.15E\n",
                        up_mass,
                        down_mass,
                        up_quark_zeroed_gap_eq,
                        down_quark_zeroed_gap_eq);

                if (up_quark_zeroed_gap_eq == down_quark_zeroed_gap_eq)
                    fprintf(line,
                            "%20.15E\t%20.15E\t%20.15E\n",
                            up_mass,
                            down_mass,
                            up_quark_zeroed_gap_eq);

                up_mass += up_mass_step;
            }
            fprintf(output, "\n"); // separate scanlines

            down_mass += down_mass_step;
        }

        fclose(output);
        fclose(line);

        SetParametersSet(NULL, NULL);
    }

    // Maps of hadron simultaneous equations
    if (true)
    {
        SetParametersSet(NULL, "eNJL2mSigmaRho1");
        SetFilePath("tests/hadron_simultaneous_equations/data");

        const int num_pts_mass = 180;
        const double min_mass = 0.0;
        const double max_mass = 1000;

        const int num_pts_dens = 200;
        const double min_density = 0.0;
        const double max_density = 1.5;

        const double barionic_chemical_potential = 1800.0;
        const double isovector_chemical_potential = 150.0;

        const double zero_tol = 10.0;

        ///

        hadron_mass_and_renorm_chem_pot_input_params p;
        p.proton_chemical_potential =
        ProtonChemicalPotential(barionic_chemical_potential,
                                isovector_chemical_potential);

        p.neutron_chemical_potential =
        NeutronChemicalPotential(barionic_chemical_potential,
                                 isovector_chemical_potential);

        const double hadron_mass_step = Step(min_mass, max_mass, num_pts_mass);
        const double density_step = Step(min_density, max_density, num_pts_dens);

        const int dimension = 3;

        FILE * output_file[3];
        FILE * zero_region_file[3];

        gsl_vector * input = gsl_vector_alloc(dimension);
        gsl_vector * output = gsl_vector_alloc(dimension);

        double hadron_mass = min_mass;
        for(int i = 0; i < num_pts_mass; i++){

            double sols_p_dens[num_pts_dens * num_pts_dens];
            double sols_n_dens[num_pts_dens * num_pts_dens];
            int solutions = 0;

            char filename[256];
            sprintf(filename, "mass_%d.dat", i);
            FILE * mass_file = OpenFile (filename);
            fprintf(mass_file, "\"m = %1.4f\"", hadron_mass);
            fclose(mass_file);

            for (int index = 0; index < dimension; index++){
                sprintf(filename, "equation_%d_%d.dat", index, i);
                output_file[index] = OpenFile(filename);

                sprintf(filename, "zero_reg_%d_%d.dat", index, i);
                zero_region_file[index] = OpenFile(filename);
            }

            double proton_density = min_density;
            for(int j = 0; j < num_pts_dens; j++){

                double neutron_density = min_density;
                for (int k = 0; k < num_pts_dens; k++){

                    gsl_vector_set(input,
                                   0,
                                   sqrt(hadron_mass));

                    gsl_vector_set(input,
                                   1,
                                   sqrt(proton_density));

                    gsl_vector_set(input,
                                   2,
                                   sqrt(neutron_density));

                    HadronMassAndDensitiesSolutionEquation(input,
                                                           (void *)&p,
                                                           output);

                    bool is_solution = true;
                    for (int index = 0; index < dimension; index++){

                        double val = gsl_vector_get(output, index);

                        fprintf(output_file[index],
                                "%20.15E\t%20.15E\t%20.15E\n",
                                proton_density,
                                neutron_density,
                                val);

                        if (fabs(val) > zero_tol){

                            is_solution = false;
                            //break;
                        }

                        if (fabs(val) < zero_tol){
                            fprintf(zero_region_file[index],
                                    "%20.15E\t%20.15E\t%20.15E\n",
                                    proton_density,
                                    neutron_density,
                                    val);
                        }
                    }
                    if (is_solution){
                        sols_p_dens[solutions] = proton_density;
                        sols_n_dens[solutions] = neutron_density;
                        solutions++;
                    }

                    neutron_density += density_step;
                }

                proton_density += density_step;
            }

            if (solutions > 0){

                double sum_p_dens = 0;
                double sum_n_dens = 0;
                for (int n = 0; n < solutions; n++){
                    sum_p_dens += sols_p_dens[n];
                    sum_n_dens += sols_n_dens[n];
                }

                double mean_p_dens = sum_p_dens / (double)solutions;
                double mean_n_dens = sum_n_dens / (double)solutions;

                FILE * intersection_file;
                sprintf(filename, "intersection_%d.dat", i);
                intersection_file = OpenFile(filename);

                fprintf(intersection_file,
                        "%20.15E\t%20.15E\n",
                        mean_p_dens,
                        mean_n_dens);
            }

            hadron_mass += hadron_mass_step;

            for (int index = 0; index < dimension; index++){
                fclose(output_file[index]);
                fclose(zero_region_file[index]);
            }
        }

        SetParametersSet(NULL, NULL);
    }
    return;
}


