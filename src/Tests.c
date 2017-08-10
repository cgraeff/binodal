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
    if(true)
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
    if(true)
    {
        SetParametersSet("Buballa_1", NULL);
        SetFilePath("tests/quark_pressure/data");

        double min_up_chemical_potential = 0.0;
        double max_up_chemical_potential = 1000.0;
        double min_down_chemical_potential = 0.0;
        double max_down_chemical_potential = 1000.0;

        int n_pts = 100;

        double temperature = 0.0;

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
                                       &up_quark_mass,
                                       &down_quark_mass,
                                       &quark_pressure);

                fprintf(file,
                        "%20.15E\t%20.15E\t%20.15E\n",
                        up_chemical_potential,
                        down_chemical_potential,
                        quark_pressure);

                down_chemical_potential += down_chem_pot_step;
            }
            fprintf(file, "\n");

            up_chemical_potential += up_chem_pot_step;
        }

        fclose(file);

        SetParametersSet(NULL, NULL);
    }

    // Determine pressures for hadron and quark phases for a
    // particular value of proton fraction
    if(true)
    {
        SetParametersSet("PCP-0.0", "eNJL1");
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

        double min_barionic_density = 0.0;
        double max_barionic_density = 0.8;

        double proton_fraction = 0.5;
        double temperature = 0.0;

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

        double bar_dens_step = Step(min_barionic_density,
                                    max_barionic_density,
                                    points_number);

        double barionic_density = min_barionic_density;
        for (int i = 0; i < points_number; i++){

            double hadron_mass;
            double hadron_pressure;
            double proton_chemical_potential;
            double neutron_chemical_potential;
            DetermineHadronPressureAndChemPots(barionic_density,
                                               proton_fraction,
                                               hadron_vacuum_potential,
                                               &hadron_mass,
                                               &proton_chemical_potential,
                                               &neutron_chemical_potential,
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
                                   &up_quark_mass,
                                   &down_quark_mass,
                                   &quark_pressure);

            double barionic_chemical_potential =
            (proton_chemical_potential + neutron_chemical_potential) / 2.0;

            fprintf(file_h,
                    "%20.15E\t%20.15E\n",
                    barionic_chemical_potential,
                    hadron_pressure);

            fprintf(file_q,
                    "%20.15E\t%20.15E\n",
                    barionic_chemical_potential,
                    quark_pressure);

            barionic_density += bar_dens_step;
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
    if (true)
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

    return;
}


