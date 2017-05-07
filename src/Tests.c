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
    SetParametersSet(NULL, NULL);

    // Write functions to be zeroed in some way
    if (false)
    {
        SetParametersSet(NULL, NULL);

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
        SetParametersSet(NULL, NULL);

        int n_pts = 1000;

        double mu_r_min = 0.0;
        double mu_r_max = 1000.0;

        double mu_r_step = Step(mu_r_min, mu_r_max, n_pts);

        renorm_chem_pot_equation_input params;
        params.chemical_potential = 315.0;
        params.mass = 400.0;
	    params.temperature = 0.0;

        SetFilePath ("tests/renor_chem_pot");
        FILE * file = OpenFile("zeroed_q_renorm_chem_pot.dat");

        double mu_r = 0;
        for (int i = 0; i < n_pts; i++){
            double result =
                ZeroedRenormalizedQuarkChemicalPotentialEquation(mu_r, (void *)&params);

            fprintf(file, "%20.15E\t%20.15E\n", mu_r, result);

            mu_r += mu_r_step;
        }
    }

    // Test pressures
    if(true){

        SetParametersSet(NULL, NULL);

        const double barionic_density_min = 0.01;
        const double barionic_density_max = 0.3;
        const double hadron_mass = 939.0;
        const double up_quark_mass = 400.0;
        const double down_quark_mass = 400.0;
        const double proton_fraction = 0.5;
        const double temperature = 0;
        const int n_pts = 1000;

        SetFilePath ("tests/Pressures");

        FILE * hadron_pressure_file = OpenFile("hadron_pressure_file.dat");
        FILE * quark_pressure_file = OpenFile ("quark_pressure_file.dat");

        double quark_vacuum_mass = QuarkVacuumMassDetermination();
        double quark_vacuum_thermodynamic_potential =
            QuarkThermodynamicPotential(quark_vacuum_mass,
                                        0.0,
                                        0.0,
                                        0.0);

        double barionic_density_step = Step(barionic_density_min,
                                            barionic_density_max,
                                            n_pts);

        double barionic_density = barionic_density_min;

        for (int i = 0; i < n_pts; i++){

            // Hadrons:
            double proton_density = proton_fraction * barionic_density;
            double neutron_density = (1.0 - proton_fraction) * barionic_density;

            double proton_fermi_momentum = HPFermiMomentum(proton_density);
            double neutron_fermi_momentum = HPFermiMomentum(neutron_density);

            double total_hadron_scalar_density = HadronScalarDensity(hadron_mass,
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

            // Gibbs Conditions:

            // Use chemical potential equality to determine
            // quark chemical potentials
            double up_chemical_potential = (2.0 * proton_chemical_potential
                                            - neutron_chemical_potential) / 3.0;

            double down_chemical_potential = (-proton_chemical_potential
                                              + 2.0 * neutron_chemical_potential) / 3.0;

            // up quark
            double up_renormalized_chemical_potential =
                QuarkSelfConsistentRenormChemPot(up_quark_mass, up_chemical_potential, temperature);

            // down quark
            double down_renormalized_chemical_potential =
                QuarkSelfConsistentRenormChemPot(down_quark_mass, down_chemical_potential, temperature);

            // Gibbs' conditions:
            // Determination of hadron pressure

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

            // Determination o quark pressure
            double up_quark_thermodynamic_potential =
                    QuarkThermodynamicPotential(up_quark_mass,
                                                up_chemical_potential,
                                                up_renormalized_chemical_potential,
                                                parameters.variables.temperature);
            double down_quark_thermodynamic_potential =
                    QuarkThermodynamicPotential(down_quark_mass,
                                                down_chemical_potential,
                                                down_renormalized_chemical_potential,
                                                parameters.variables.temperature);

            double regularized_thermodynamic_potential = up_quark_thermodynamic_potential
                                                         + down_quark_thermodynamic_potential
                                                         - 2.0 * quark_vacuum_thermodynamic_potential;

            double quark_pressure = QuarkPressure(regularized_thermodynamic_potential,
                                                    parameters.variables.temperature);

            fprintf(hadron_pressure_file,
                    "%20.15E\t%20.15E\n",
                    barionic_density,
                    hadron_pressure);

            fprintf(quark_pressure_file,
                    "%20.15E\t%20.15E\n",
                    barionic_density,
                    quark_pressure);

            barionic_density += barionic_density_step;
        }

        fclose(hadron_pressure_file);
        fclose(quark_pressure_file);

        SetFilePath(NULL);
    }

    return;
}
