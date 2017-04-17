//
//  Tests.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-06-08.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include "Tests.h"

#include "libdatafun/libdatafun.h"

#include "Parameters.h"
#include "SimultaneousSolution.h"
#include "QuarkPhaseEOS.h"
#include "HadronPhaseEOS.h"

void RunTests()
{
    // Write functions to be zeroed in some way
    if (true)
    {
        int n_pts = 1000;
        double temperature = 0.0;
        double proton_fraction = 0.5;
        double barionic_density = 0.15;

        double h_mass_min = 0.1;
        double h_mass_max = 1000.0;

        double u_q_mass_min = 0.1;
        double u_q_mass_max = 600.0;

        double d_q_mass_min = 0.1;
        double d_q_mass_max = 600.0;

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

        multi_dim_root_params params;
        params.temperature = temperature;
        params.proton_fraction = proton_fraction;
        params.quark_vacuum_thermodynamic_potential =
        quark_vacuum_thermodynamic_potential;
        params.hadron_vacuum_thermodynamic_potential =
                            hadron_vacuum_thermodynamic_potential;

        double h_mass_step = Step(h_mass_min, h_mass_max, n_pts);
        double u_q_mass_step = Step(u_q_mass_min, u_q_mass_max, n_pts);
        double d_q_mass_step = Step(d_q_mass_min, d_q_mass_max, n_pts);

        gsl_vector * input = gsl_vector_alloc(4);
        gsl_vector * output = gsl_vector_alloc(4);

        FILE * files[4];
        SetFilePath ("tests/map");
        files[0] = OpenFile("first.dat");
        files[1] = OpenFile("second.dat");
        files[2] = OpenFile("third.dat");
        files[3] = OpenFile("fourth.dat");

        double h_mass = h_mass_min;
        double u_q_mass = u_q_mass_min;
        double d_q_mass = d_q_mass_min;

        for (int j = 0; j < n_pts; j++){
            for (int k = 0; k < n_pts; k++){
                for (int l = 0; l < n_pts; l++){

                    gsl_vector_set(input, 0, barionic_density);
                    gsl_vector_set(input, 1, h_mass);
                    gsl_vector_set(input, 2, u_q_mass);
                    gsl_vector_set(input, 3, d_q_mass);

                    MultiDimensionalRootFinderHelperFunction(input,
                                                             (void *)&params,
                                                             output);

                    for (int m = 0; m < 4; m++){
                        fprintf(files[m],
                                "%20.15E\t%20.15E\t%20.15E\t%20.15E\n",
                                h_mass,
                                u_q_mass,
                                d_q_mass,
                                gsl_vector_get(output, m));
                    }

                    d_q_mass += u_q_mass_step;
                }
                u_q_mass += d_q_mass_step;
            }
            h_mass += h_mass_step;
        }

        SetFilePath (NULL);
    }

    if (true)
    {
        int n_pts = 1000;

        double mu_r_min = 0.0;
        double mu_r_max = 1000.0;

        double mu_r_step = Step(mu_r_min, mu_r_max, n_pts);

        renorm_chem_pot_equation_input params;
        params.chemical_potential = 45.595386;
        params.mass = 0.01;
	    params.temperature = 0.0;

        SetFilePath ("tests");
        FILE * file = OpenFile("zeroed_q_renorm_chem_pot.dat");

        double mu_r = 0;
        for (int i = 0; i < n_pts; i++){
            double result =
                ZeroedRenormalizedChemicalPotentialEquation(mu_r, (void *)&params);

            fprintf(file, "%20.15E\t%20.15E\n", mu_r, result);

            mu_r += mu_r_step;
        }
    }

    return;
}
