//
//  Tests.c
//  binodal
//
//  Created by Clebson Graeff on 2016-06-08.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include <stdbool.h>

#include "libdatafun/libdatafun.h"

#include "CommandlineOptions.h"

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
    if(true)
    {
        SetParametersSet("Buballa_1", NULL);
        SetFilePath("tests/quark_thermodynamic_potential/data");

        printf("- Quark Thermodynamic Potential\n");

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
    if(true)
    {
        SetParametersSet("Buballa_1", NULL);
        SetFilePath("tests/quark_pressure/data");

        printf("- Quark pressures as functions of quark chemical potentials\n");

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

                DetermineQuarkPressureAndMasses(up_chemical_potential,
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

    // Determine quark and hadron pressures as functions of chemical potential
    // for a given proton fraction.
    //
    // Note that this will be the same graphic as the binodal point, but here
    // we want it to cover a bigger range of values of barionic chemical
    // potential.
    if(true)
    {
        SetParametersSet("BuballaR_2", "eNJL2mSigmaRho1");
        SetFilePath("tests/quark_and_hadron_pressures/data");

        printf("- Quark and hadron pressures as functions of "
               "chemical potential\n");

        int n_pts = 100;

        double min_barionic_chemical_potential = 1050.0; // Try to widen this range
        double max_barionic_chemical_potential = 2000.0;
        double isovector_chemical_potential = 0.0; // write in terms of
                                                   // proton_fraction?
        double initial_proton_density_guess = 0.1;
        double initial_neutron_density_guess = 0.1;
        double initial_hadron_mass_guess = 900.0;

        double initial_up_mass_guess = 313.0;
        double initial_down_mass_guess = 313.0;

        double temperature = 0.0;


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
            printf("\tSolving gap equation and equations of state ...\n");

        FILE * file =
        OpenFile("hadron_and_quark_pressure_vs_barionic_chem_pot.dat");
        fprintf (file,
                 "# isovector_chemical_potential: %f\n",
                 isovector_chemical_potential);
        fprintf (file,
                 "# barionic chemical potential (MeV/fm^3), "
                 "hadron pressure (MeV/fm^3), "
                 "hadron mass (MeV), "
                 "proton density (fm^{-3}), "
                 "neutron density (fm^{-3}), "
                 "quark pressure (MeV/fm^3), "
                 "up quark mass (MeV), "
                 "down quark mass (MeV)\n");

        double barionic_chemical_potential_step =
        Step(min_barionic_chemical_potential,
             max_barionic_chemical_potential,
             n_pts);


        double hadron_mass_guess = initial_hadron_mass_guess;
        double proton_density_guess = initial_proton_density_guess;
        double neutron_density_guess = initial_neutron_density_guess;

        double up_mass_guess = initial_up_mass_guess;
        double down_mass_guess = initial_down_mass_guess;

        double barionic_chemical_potential = min_barionic_chemical_potential;
        for (int i = 0; i < n_pts; i++){

            // Determine hadron pressure for given barionic chemical potential

            double proton_chemical_potential =
            ProtonChemicalPotential(barionic_chemical_potential,
                                    isovector_chemical_potential);

            double neutron_chemical_potential =
            NeutronChemicalPotential(barionic_chemical_potential,
                                     isovector_chemical_potential);

            double hadron_mass;
            double hadron_pressure;
            double proton_density;
            double neutron_density;

            int status_h =
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

            if (status_h){
                printf("%s:%d: Problems with determination of hadron pressure.\n",
                       __FILE__,
                       __LINE__);
                abort();
            }

            double up_chemical_potential =
            UpChemicalPotentialFromGibbsConditions(proton_chemical_potential,
                                                   neutron_chemical_potential);

            double down_chemical_potential =
            DownChemicalPotentialFromGibbsConditions(proton_chemical_potential,
                                                     neutron_chemical_potential);


            // Determine quark pressure for given barionic chemical potential

            double up_quark_mass;
            double down_quark_mass;
            double quark_pressure;

            int status_q =
            DetermineQuarkPressureAndMasses(up_chemical_potential,
                                            down_chemical_potential,
                                            temperature,
                                            quark_vacuum_potential,
                                            up_mass_guess,
                                            down_mass_guess,
                                            &up_quark_mass,
                                            &down_quark_mass,
                                            &quark_pressure);

            if (status_q){
                printf("%s:%d: Problems with determination of quark pressure.\n",
                       __FILE__,
                       __LINE__);
                abort();
            }

            fprintf(file,
                    "%20.15E\t%20.15E\t%20.15E\t%20.15E"
                    "\t%20.15E\t%20.15E\t%20.15E\t%20.15E\n",
                    barionic_chemical_potential,
                    hadron_pressure,
                    hadron_mass,
                    proton_density,
                    neutron_density,
                    quark_pressure,
                    up_quark_mass,
                    down_quark_mass);

            // Use solutions as guesses for next point
            hadron_mass_guess = hadron_mass;
            proton_density_guess = proton_density;
            neutron_density_guess = neutron_density;

            up_mass_guess = up_quark_mass;
            down_mass_guess = down_quark_mass;

            barionic_chemical_potential += barionic_chemical_potential_step;
        }

        fclose(file);

        SetParametersSet(NULL, NULL);
    }

    // Determine pressures for hadron and quark phases for a
    // particular value of the isovector chemical potential
    if(true)
    {
        printf("- Binodal point\n");

        int points_number = 1000;
        int n_isovec_chem_pot_pts = 1;

        double min_isovector_chemical_potential = 0.0;
        double max_isovector_chemical_potential = 0.0;

        double min_barionic_chemical_potential = 1050.0;
        double max_barionic_chemical_potential = 2000.0;

        double temperature = 0.0;

        double starting_hadron_mass_guess = 900.0;
        double starting_proton_density_guess = 0.3;
        double starting_neutron_density_guess = 0.3;

        double starting_up_mass_guess = 313.0;
        double starting_down_mass_guess = 313.0;

        const int n_h_sets = 2;
        const int n_q_sets = 6;

        char quark_sets[6][256] =
        {
            "PCP-0.0",
            "PCP-0.1",
            "PCP-0.2",
            "BuballaR_2",
            "Buballa_1",
            "Buballa_2"
        };

        char hadron_sets[2][256] =
        {
            "eNJL2mSigmaRho1",
            "eNJL3SigmaRho1"
        };

        double isovec_chem_pot_step =
        Step(min_isovector_chemical_potential,
             max_isovector_chemical_potential,
             n_isovec_chem_pot_pts);

        double barionic_chemical_potential_step =
        Step(min_barionic_chemical_potential,
             max_barionic_chemical_potential,
             points_number);

        double initial_hadron_mass_guess = starting_hadron_mass_guess;
        double initial_proton_density_guess = starting_proton_density_guess;
        double initial_neutron_density_guess = starting_neutron_density_guess;

        double initial_up_mass_guess = starting_up_mass_guess;
        double initial_down_mass_guess = starting_down_mass_guess;

        for (int h_set = 0; h_set < n_h_sets; h_set++){
            for (int q_set = 0; q_set < n_q_sets; q_set++){

                printf("\t%s-%s\n", quark_sets[q_set], hadron_sets[h_set]);

                SetParametersSet(quark_sets[q_set], hadron_sets[h_set]);

                char dir_name[256];
                sprintf(dir_name,
                        "tests/binodal_point_graph/data/%s-%s",
                        quark_sets[q_set],
                        hadron_sets[h_set]);

                SetFilePath(dir_name);

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

                printf("\t\tUp quark vacuum mass: %f\n", up_vacuum_mass);
                printf("\t\tDown quark vacuum mass: %f\n", down_vacuum_mass);
                printf("\t\tQuark vacuum thermodynamic potential: %f\n",
                       quark_vacuum_potential);
                printf("\t\tHadron vacuum thermodynamic potential: %f\n",
                       hadron_vacuum_potential);

                double isovector_chemical_potential =
                min_isovector_chemical_potential;
                for (int k = 0; k < n_isovec_chem_pot_pts; k++){

                    char filename[256];
                    sprintf(filename,"hadron_pressure_%d.dat", k);
                    FILE * file_h = OpenFile(filename);

                    sprintf(filename, "quark_pressure_%d.dat", k);
                    FILE * file_q = OpenFile(filename);

                    // printf file headers
                    fprintf(file_h,
                            "# barionic_chemical_potential (MeV), "
                            "hadron pressure (MeV/fm^3)\n");
                    fprintf(file_q,
                            "# barionic_chemical_potential (MeV), "
                            "quark pressure (MeV/fm^3)\n");

                    double hadron_mass_guess = initial_hadron_mass_guess;
                    double proton_density_guess = initial_proton_density_guess;
                    double neutron_density_guess = initial_neutron_density_guess;

                    double up_mass_guess = initial_up_mass_guess;
                    double down_mass_guess = initial_down_mass_guess;

                    double barionic_chemical_potential =
                    min_barionic_chemical_potential;

                    bool quark_pressure_is_bigger = false;

                    for (int i = 0; i < points_number; i++){

                        BinodalPoint point;

                        int status =
                        BinodalPointCandidate(barionic_chemical_potential,
                                              isovector_chemical_potential,
                                              temperature,
                                              hadron_vacuum_potential,
                                              hadron_mass_guess,
                                              proton_density_guess,
                                              neutron_density_guess,
                                              quark_vacuum_potential,
                                              up_mass_guess,
                                              down_mass_guess,
                                              &point);

                        if (options.debug)
                            if (status){
                                printf("%s:%d: Problems with determination of the "
                                       "binodal point candidate.\n",
                                       __FILE__,
                                       __LINE__);
                                printf("h_set: %d\nq_set: %d\nk: %d\ni: %d\n",
                                       h_set,
                                       q_set,
                                       k,
                                       i);

                                printf("bar_chm_pot: %f\n", barionic_chemical_potential);
                                printf("iso_chm_pot: %f\n", isovector_chemical_potential);

                                barionic_chemical_potential +=
                                barionic_chemical_potential_step;
                            }

                        if (quark_pressure_is_bigger == false)
                            if (point.quark_pressure > point.hadron_pressure){

                                printf("\t\tPressure at transition: %f\n",
                                       point.hadron_pressure);
                                printf("\t\tChemical potential at transition: %f\n",
                                       BarionicChemicalPotential(point.proton_chemical_potential,
                                                                 point.neutron_chemical_potential));
                                printf("\t\tBarionic density at transition: %f\n",
                                       point.proton_density + point.neutron_density);

                                quark_pressure_is_bigger = true;
                            }

                        fprintf(file_h,
                                "%20.15E\t%20.15E\t%20.15E\t%20.15E\t%20.15E\t%20.15E\n",
                                barionic_chemical_potential,
                                point.hadron_pressure,
                                isovector_chemical_potential,
                                point.hadron_mass,
                                point.proton_density,
                                point.neutron_density);

                        fprintf(file_q,
                                "%20.15E\t%20.15E\t%20.15E\t%20.15E\t%20.15E\t%20.15E\t%20.15E\n",
                                barionic_chemical_potential,
                                point.quark_pressure,
                                isovector_chemical_potential,
                                point.up_quark_mass,
                                point.down_quark_mass,
                                point.up_quark_density,
                                point.down_quark_density);

                        // Update guesses
                        hadron_mass_guess = point.hadron_mass;

                        proton_density_guess = point.proton_density;

                        neutron_density_guess = point.neutron_density;

                        up_mass_guess = point.up_quark_mass;

                        down_mass_guess = point.down_quark_mass;

                        barionic_chemical_potential +=
                        barionic_chemical_potential_step;

                        if (k == 0 && i == 0){
                            initial_hadron_mass_guess = hadron_mass_guess;

                            initial_proton_density_guess = proton_density_guess;

                            initial_neutron_density_guess = neutron_density_guess;

                            initial_up_mass_guess = up_mass_guess;

                            initial_down_mass_guess = down_mass_guess;
                        }
                    }

                    fclose(file_h);
                    fclose(file_q);
                    isovector_chemical_potential += isovec_chem_pot_step;
                }
            }
        }

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
    if (true)
    {
        SetParametersSet("BuballaR_2", NULL);
        SetFilePath("tests/quark_vacuum_mass_maps/data");

        printf("- Quark gap equations\n");

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
    if (true)
    {
        SetParametersSet("BuballaR_2", NULL);
        SetFilePath("tests/quark_mass_maps/data");

        printf("- Quark mass and renormalized chemical potentials maps\n");
        int n_pts = 1000;

        double min_up_mass = 0.0;
        double max_up_mass = 600.0;
        double min_down_mass = 0.0;
        double max_down_mass = 600.0;

        double zero_tol = 1.0;

        double temperature = 0.0;

//        double up_chemical_potential = 366.404035;
//        double down_chemical_potential = 433.062365;

        double up_chemical_potential = 400.0;
        double down_chemical_potential = 400.0;

        FILE * output = OpenFile("maps.dat");
        FILE * line = OpenFile("line.dat");
        FILE * near_zero1 = OpenFile("near_zero1.dat");
        FILE * near_zero2 = OpenFile("near_zero2.dat");

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

                if (fabs(up_quark_zeroed_gap_eq) < zero_tol)
                    fprintf(near_zero1,
                            "%20.15E\t%20.15E\n",
                            up_mass,
                            down_mass);

                if (fabs(down_quark_zeroed_gap_eq) < zero_tol)
                    fprintf(near_zero2,
                            "%20.15E\t%20.15E\n",
                            up_mass,
                            down_mass);

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

        printf("- Hadron gap equations map\n");

        const int num_pts_mass = 180;
        const double min_mass = 0.0;
        const double max_mass = 1000;

        const int num_pts_dens = 200;
        const double min_density = 0.0;
        const double max_density = 1.5;

        const double barionic_chemical_potential = 1052.8509503167727;
        const double isovector_chemical_potential = 168.42105263157902;

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

    // Figure of hadron mass gap equation
    if (true)
    {
        SetParametersSet(NULL,"eNJL2mSigmaRho1");
        SetFilePath("tests/hadron_gap_equation/data");

        printf("- Hadron gap equation\n");

        int n_pts = 1000;
        int n_dens_pts = 30;

        double min_mass = 0.0;
        double max_mass = 1200.0;

        double min_density = 0.0;
        double max_density = 1.0;

        double proton_fraction = 0.2;

        double dens_step = Step(min_density, max_density, n_dens_pts);
        double mass_step = Step(min_mass, max_mass, n_pts);

        char filename[256];

        double dens = min_density;
        for (int n = 0; n < n_dens_pts; n++){

            sprintf(filename, "hadron_zeroed_gap_eq_%d.dat", n);
            FILE * output = OpenFile(filename);

            double mass = min_mass;
            for (int i = 0; i < n_pts; i++){

                double proton_density = proton_fraction * dens;
                double neutron_density = (1.0 - proton_fraction) * dens;

                double proton_fermi_momentum =
                HadronFermiMomentum(proton_density);

                double neutron_fermi_momentum =
                HadronFermiMomentum(neutron_density);

                hadron_gap_eq_input_params  gap_params;
                gap_params.proton_fermi_momentum = proton_fermi_momentum;
                gap_params.neutron_fermi_momentum = neutron_fermi_momentum;
                gap_params.proton_density = proton_density;
                gap_params.neutron_density = neutron_density;

                double gap_equation =
                HadronZeroedGapEquation(mass,
                                        (void *)&gap_params);

                fprintf(output,
                        "%20.15E\t%20.15E\n",
                        mass,
                        gap_equation);

                mass += mass_step;
            }

            dens += dens_step;
            fclose(output);
        }
        SetFilePath(NULL);
        SetParametersSet(NULL, NULL);
    }

    // Make graph of quark zeroed functions as functions
    // of mass. In this test we should see if there are
    // multiple solutions for given chemical potentials.
    if (true)
    {
        SetParametersSet("BuballaR_2", NULL);

        SetFilePath("tests/quark_zeroed_functions_for_bissection/data/");

        printf("- Quark gap equation\n");

        int num_points = 1000;

        double temperature = 0.0;

        double min_mass = 0.1;
        double max_mass = 600.0;

//        double up_chemical_potential = 366.404035;
//        double down_chemical_potential = 433.062365;


        double up_chemical_potential = 350.0;
        double down_chemical_potential = 350.0;

        FILE * test_file = OpenFile("quark_zeroed_functions.dat");
        FILE * solution_file = OpenFile("solution.dat");

        // Set up parameters to be passed to helper function
        quark_mass_and_renorm_chem_pot_input_params p;
        p.up_chemical_potential = up_chemical_potential;
        p.down_chemical_potential = down_chemical_potential;
        p.up_renorm_chem_pot = NAN;
        p.down_renorm_chem_pot = NAN;

        double mass_step = Step(min_mass, max_mass, num_points);

        double mass = min_mass;
        for (int i = 0; i < num_points; i++){

            fprintf(test_file,
                    "%20.15E\t%20.15E\n",
                    mass,
                    MyAdapterFunction(mass, (void *)&p));

            mass += mass_step;
        }

        double up_mass;
        double down_mass;
        double up_renorm_chem_pot;
        double down_renorm_chem_pot;

        int status =
        QuarkMassAndRenormChemPotSolutionBissection(up_chemical_potential,
                                                    down_chemical_potential,
                                                    temperature,
                                                    &up_mass,
                                                    &down_mass,
                                                    &up_renorm_chem_pot,
                                                    &down_renorm_chem_pot);

        if (!status)
            fprintf(solution_file,
                    "%20.15E\t0.0\n",
                    up_mass);
    }

    return;
}


