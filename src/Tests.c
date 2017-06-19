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

int TestMassAndRenormChemPotSimultaneousSolution(double up_chemical_potential,
                                                 double down_chemical_potential,
                                                 double up_mass_guess,
                                                 double down_mass_guess,
                                                 double abs_error,
                                                 double rel_error,
                                                 int max_iter,
                                                 double * return_up_mass,
                                                 double * return_down_mass,
                                                 double * return_up_renorm_chem_pot,
                                                 double * return_down_renorm_chem_pot);

void BinodalPoint(double temperature,
                  double proton_fraction,
                  double *return_barionic_chemical_potential,
                  double *return_isovector_chemical_potential,
                  double *return_pressure);

void RunTests()
{
    SetParametersSet(NULL, NULL);

    if (true){
        SetParametersSet("Buballa_1", "eNJL1");

        double return_barionic_chemical_potential = NAN;
        double return_isovector_chemical_potential = NAN;
        double return_pressure = NAN;

        BinodalPoint(0.0,
                     0.5,
                     &return_barionic_chemical_potential,
                     &return_isovector_chemical_potential,
                     &return_pressure);

    }
    // Determine thermodynamic potential map:
    //      To plot on gnuplot use
    //          set view map
    //          splot "quark_vacuum_mass_map.dat" w pm3d
    //      To make contour lines and adjust the number of lines
    //          set contour
    //          set cntrparam levels 35
    if(true)
    {
        SetParametersSet("PCP-0.0", NULL);
        SetFilePath ("tests/thermodynamic_potential_map");

        double min_mass = 0.0;
        double max_mass = 1000.0;
        int n_pts = 100;

        double up_chemical_potential = 0.0;
        double down_chemical_potential = 0.0;
        double up_renormalized_chemical_potential = 0.0;
        double down_renormalized_chemical_potential = 0.0;
        double temperature = 0.0;

        FILE * file = OpenFile("thermodynamic_potential_map.dat");

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

        file = OpenFile("thermodynamic_potential_minimum_location.dat");

        fprintf (file,
                 "%20.15E\t%20.15E\t%20.15E\n",
                 up_mass_at_minimum,
                 down_mass_at_minimum,
                 potential_at_minimum);

        fclose(file);
    }

    // Determine quark thermodynamic potential and pressure as functions of mass
    // for m_u = m_d
    if(true)
    {
        SetParametersSet("PCP-0.0", NULL);
        SetFilePath ("tests/quark_thermodynamic_potential_equal_masses");

        double min_mass = 0.0;
        double max_mass = 1000.0;
        int n_pts = 100;

        double up_chemical_potential = 0.0;
        double down_chemical_potential = 0.0;
        double up_renormalized_chemical_potential = 0.0;
        double down_renormalized_chemical_potential = 0.0;
        double temperature = 0.0;

        double up_vacuum_mass;
        double down_vacuum_mass;

        QuarkVacuumMassDetermination (&up_vacuum_mass, &down_vacuum_mass);

        double vacuum_potential =
            QuarkThermodynamicPotential(up_vacuum_mass,
                                        down_vacuum_mass,
                                        0.0,
                                        0.0,
                                        0.0,
                                        0.0,
                                        0.0);

        FILE * file = OpenFile("quark_thermodynamic_potential_vs_mass_for_equal_masses.dat");
        FILE * filep = OpenFile("quark_pressure_vs_mass_for_equal_masses.dat");

        double mass_step = Step (min_mass, max_mass, n_pts);

        double mass = 0.0;
        for (int i = 0; i < n_pts; i++){

            // FIXME: the chemical_potential and the renormalized_chemical_potential,
            // maybe even the barionic density, should depend on the value of mass
            //
            double omega =
                QuarkThermodynamicPotential(mass,
                                            mass,
                                            up_chemical_potential,
                                            down_chemical_potential,
                                            up_renormalized_chemical_potential,
                                            down_renormalized_chemical_potential,
                                            temperature);

                fprintf(file,
                        "%20.15E\t%20.15E\t\n",
                        mass,
                        omega - vacuum_potential);

                fprintf(filep,
                        "%20.15E\t%20.15E\n",
                        mass,
                        HadronPressure (omega - vacuum_potential));

                mass += mass_step;
            }

        fclose(file);

        SetParametersSet(NULL, NULL);
    }

    // Determine quark thermodynamic potential and pressure as functions of
    // chemical potential for m_u = m_d
    if(true)
    {
        SetParametersSet("Buballa_1", NULL);
        SetFilePath ("tests/quark_thermodynamic_potential_equal_masses");

        double min_up_chemical_potential = 0.0;
        double max_up_chemical_potential = 1000.0;
        double min_down_chemical_potential = 0.0;
        double max_down_chemical_potential = 1000.0;

        double up_mass_guess = 300.0;
        double down_mass_guess = 300.0;
        double abs_error = 1.0E-5;
        double rel_error = 1.0E-5;
        int max_iter = 1000;

        int n_pts = 100;

        double temperature = 0.0;

        double up_vacuum_mass;
        double down_vacuum_mass;

        QuarkVacuumMassDetermination (&up_vacuum_mass, &down_vacuum_mass);

        printf("vacuum masses: %f, %f\n", up_vacuum_mass, down_vacuum_mass);

        double vacuum_potential =
            QuarkThermodynamicPotential(up_vacuum_mass,
                                        down_vacuum_mass,
                                        0.0,
                                        0.0,
                                        0.0,
                                        0.0,
                                        0.0);

        printf("Vacuum potential: %f\n", vacuum_potential);

        FILE * file = OpenFile("quark_thermodynamic_potential_"
                               "vs_chem_pot_for_equal_masses.dat");
        FILE * file_equal = OpenFile("quark_thermodynamic_potential_"
                                     "vs_chem_pot_for_equal_masses_and_chem_pot.dat");
        FILE * filep_equal = OpenFile("pressure_"
                                     "vs_chem_pot_for_equal_masses_and_chem_pot.dat");

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

                double return_up_mass = NAN;
                double return_down_mass = NAN;
                double return_up_renorm_chem_pot = NAN;
                double return_down_renorm_chem_pot = NAN;

                TestMassAndRenormChemPotSimultaneousSolution(up_chemical_potential,
                                                             down_chemical_potential,
                                                             up_mass_guess,
                                                             down_mass_guess,
                                                             abs_error,
                                                             rel_error,
                                                             max_iter,
                                                             &return_up_mass,
                                                             &return_down_mass,
                                                             &return_up_renorm_chem_pot,
                                                             &return_down_renorm_chem_pot);

                double omega =
                    QuarkThermodynamicPotential(return_up_mass,
                                                return_down_mass,
                                                up_chemical_potential,
                                                down_chemical_potential,
                                                return_up_renorm_chem_pot,
                                                return_down_renorm_chem_pot,
                                                temperature);

                fprintf(file,
                        "%20.15E\t%20.15E\t%20.15E\t\n",
                        up_chemical_potential,
                        down_chemical_potential,
                        omega - vacuum_potential);


                if (up_chemical_potential == down_chemical_potential){

                /*    if (return_up_mass != return_down_mass){
                        printf("%s:%d: Uh, I expected to see equal masses here.\n",
                        __FILE__,
                        __LINE__);

                        printf("up_mass: %f\ndown_mass: %f\n", return_up_mass, return_down_mass);

                        abort();
                    }
                 */

                    fprintf(file_equal,
                            "%20.15E\t%20.15E\n",
                            up_chemical_potential * 3.0,
                            omega - vacuum_potential);

                    fprintf(filep_equal,
                            "%20.15E\t%20.15E\n",
                            up_chemical_potential * 3.0,
                            QuarkPressure(omega - vacuum_potential, 0.0));
                }

                down_chemical_potential += down_chem_pot_step;
            }
            fprintf(file, "\n");

            up_chemical_potential += up_chem_pot_step;
        }

        fclose(file);

        SetParametersSet(NULL, NULL);
    }

    // Determine hadron thermodynamic potential as function of mass
    if(false)
    {
        SetParametersSet(NULL, "eNJL1");
        SetFilePath ("tests/hadron_thermodynamic_potential");

        FILE * file = OpenFile ("hadron_thermodynamic_potential_mass.dat");

        int points_number = 1000;

        double mass_min = 0.0;      // (MeV)
        double mass_max = 3.0E3;    // (MeV)

        double proton_fraction = 0.5;
        double barionic_density = 0.0; // (MeV) // Should be zero

        double mass_step = Step(mass_min, mass_max, points_number);

        double proton_density = proton_fraction * barionic_density;
        double neutron_density = (1.0 - proton_fraction) * barionic_density;

        double proton_fermi_momentum = HPFermiMomentum(proton_density);
        double neutron_fermi_momentum = HPFermiMomentum(neutron_density);

        double mass = mass_min;
        for (int i = 0; i < points_number; i++){

            double kinectic_energy_density =
                HadronKinecticEnergyDensity(mass,
                                            proton_fermi_momentum,
                                            neutron_fermi_momentum);

            double scalar_density =
                HadronScalarDensity(mass, proton_fermi_momentum, parameters.hadron_model.cutoff)
                + HadronScalarDensity(mass, neutron_fermi_momentum, parameters.hadron_model.cutoff);

            // If barionic density == 0, the chemical
            // potential functions must return zero.
            double proton_chemical_potential =
                ProtonChemicalPotential(proton_fermi_momentum,
                                        scalar_density,
                                        mass,
                                        barionic_density,
                                        proton_density,
                                        neutron_density);

            double neutron_chemical_potential =
                NeutronChemicalPotential(neutron_fermi_momentum,
                                         scalar_density,
                                         mass,
                                         barionic_density,
                                         proton_density,
                                         neutron_density);

            double potential =
                HadronThermodynamicPotential(scalar_density,
                                             barionic_density,
                                             proton_density,
                                             neutron_density,
                                             proton_chemical_potential,
                                             neutron_chemical_potential,
                                             kinectic_energy_density);

                fprintf(file,
                        "%20.15E\t%20.15E\n",
                        mass,
                        potential);

                mass += mass_step;
            }


        fclose(file);

        SetParametersSet(NULL, NULL);
    }

    // hadron gap equation as funcion of mass
    if(true)
    {
        SetParametersSet(NULL, "eNJL1");
        SetFilePath ("tests/hadron_thermodynamic_potential");

        FILE * file = OpenFile ("hadron_gap_test.dat");

        int points_number = 1000;

        double barionic_denisty = 0.15;

        double proton_fraction = 0.5;

        double proton_density = proton_fraction * barionic_denisty;
        double neutron_density = (1.0 - proton_fraction) * barionic_denisty;

        double proton_fermi_momentum = HPFermiMomentum (proton_density);
        double neutron_fermi_momentum = HPFermiMomentum (neutron_density);

        double min_mass = 0.0;
        double max_mass = 3000.0;

        double mass_step = Step(min_mass,
                                max_mass,
                                points_number);

        double mass = min_mass;
        for (int i = 0; i < points_number; i++){

            hadron_gap_eq_input_params p;
            p.proton_density = proton_density;
            p.neutron_density = neutron_density;
            p.proton_fermi_momentum = proton_fermi_momentum;
            p.neutron_fermi_momentum = neutron_fermi_momentum;

            double gap = HadronZeroedGapEquation (mass, (void *)&p);

            fprintf(file,
                    "%20.15E\t%20.15E\n",
                    mass,
                    gap);

                mass += mass_step;
            }

        fclose(file);

        SetParametersSet(NULL, NULL);
    }


    // Determine hadron thermodynamic potential as function of barionic_density
    if(true)
    {
        SetParametersSet(NULL, "eNJL1");
        SetFilePath ("tests/hadron_thermodynamic_potential");

        FILE * file = OpenFile ("hadron_thermodynamic_potential_bar_density.dat");
        FILE * filep = OpenFile ("pressure_bar_density.dat");

        int points_number = 1000;

        double min_barionic_density = 0.0;
        double max_barionic_density = 0.8;

        double proton_fraction = 0.5;

        double vacuum_potential = HadronVacuumEnergyDensity();

        UnidimensionalRootFindingParameters rootfinding_params;
        rootfinding_params.max_iterations = 2000;
        rootfinding_params.lower_bound = 0.01;
        rootfinding_params.upper_bound = 4000;
        rootfinding_params.abs_error = 1E-5;
        rootfinding_params.rel_error = 1E-5;

        double bar_dens_step = Step(min_barionic_density,
                                    max_barionic_density,
                                    points_number);

        double barionic_density = min_barionic_density;
        for (int i = 0; i < points_number; i++){

            double proton_density = proton_fraction * barionic_density;
            double neutron_density = (1.0 - proton_fraction) * barionic_density;

            double proton_fermi_momentum = HPFermiMomentum(proton_density);
            double neutron_fermi_momentum = HPFermiMomentum(neutron_density);

            double mass = HadronSolveGapEquation(rootfinding_params,
                                                 proton_density,
                                                 neutron_density,
                                                 proton_fermi_momentum,
                                                 neutron_fermi_momentum);

            double kinectic_energy_density =
                HadronKinecticEnergyDensity(mass,
                                            proton_fermi_momentum,
                                            neutron_fermi_momentum);

            double scalar_density =
                HadronScalarDensity(mass,
                                    proton_fermi_momentum,
                                    parameters.hadron_model.cutoff)
                + HadronScalarDensity(mass,
                                      neutron_fermi_momentum,
                                      parameters.hadron_model.cutoff);

            // If barionic density == 0, the chemical
            // potential functions must return zero.
            double proton_chemical_potential =
                ProtonChemicalPotential(proton_fermi_momentum,
                                        scalar_density,
                                        mass,
                                        barionic_density,
                                        proton_density,
                                        neutron_density);

            double neutron_chemical_potential =
                NeutronChemicalPotential(neutron_fermi_momentum,
                                         scalar_density,
                                         mass,
                                         barionic_density,
                                         proton_density,
                                         neutron_density);

            double potential =
                HadronThermodynamicPotential(scalar_density,
                                             barionic_density,
                                             proton_density,
                                             neutron_density,
                                             proton_chemical_potential,
                                             neutron_chemical_potential,
                                             kinectic_energy_density)
                - vacuum_potential;

            fprintf(file,
                    "%20.15E\t%20.15E\n",
                    barionic_density,
                    potential);

            fprintf(filep,
                    "%20.15E\t%20.15E\n",
                    barionic_density,
                    HadronPressure (potential));

                barionic_density += bar_dens_step;
            }


        fclose(file);

        SetParametersSet(NULL, NULL);
    }

    // Determine hadron thermodynamic potential as function of chemical_potential
    // (this test and the one running on the barionic density could be fused,
    // I think.)
    if(true)
    {
        SetParametersSet(NULL, "eNJL1");
        SetFilePath ("tests/hadron_thermodynamic_potential");

        FILE * file = OpenFile ("hadron_thermodynamic_potential_chem_pot.dat");
        FILE * filep = OpenFile ("pressure_chem_pot.dat");

        int points_number = 1000;

        double min_barionic_density = 0.0;
        double max_barionic_density = 0.8;

        double proton_fraction = 0.5;

        UnidimensionalRootFindingParameters rootfinding_params;
        rootfinding_params.max_iterations = 2000;
        rootfinding_params.lower_bound = 0.01;
        rootfinding_params.upper_bound = 4000;
        rootfinding_params.abs_error = 1E-5;
        rootfinding_params.rel_error = 1E-5;

        double vacuum_potential = HadronVacuumEnergyDensity();

        double bar_dens_step = Step(min_barionic_density,
                                    max_barionic_density,
                                    points_number);

        double barionic_density = min_barionic_density;
        for (int i = 0; i < points_number; i++){

            double proton_density = proton_fraction * barionic_density;
            double neutron_density = (1.0 - proton_fraction) * barionic_density;

            double proton_fermi_momentum = HPFermiMomentum(proton_density);
            double neutron_fermi_momentum = HPFermiMomentum(neutron_density);

            double mass = HadronSolveGapEquation(rootfinding_params,
                                                 proton_density,
                                                 neutron_density,
                                                 proton_fermi_momentum,
                                                 neutron_fermi_momentum);

            double kinectic_energy_density =
                HadronKinecticEnergyDensity(mass,
                                            proton_fermi_momentum,
                                            neutron_fermi_momentum);

            double scalar_density =
                HadronScalarDensity(mass,
                                    proton_fermi_momentum,
                                    parameters.hadron_model.cutoff)
                + HadronScalarDensity(mass,
                                      neutron_fermi_momentum,
                                      parameters.hadron_model.cutoff);

            // If barionic density == 0, the chemical
            // potential functions must return zero.
            double proton_chemical_potential =
                ProtonChemicalPotential(proton_fermi_momentum,
                                        scalar_density,
                                        mass,
                                        barionic_density,
                                        proton_density,
                                        neutron_density);

            double neutron_chemical_potential =
                NeutronChemicalPotential(neutron_fermi_momentum,
                                         scalar_density,
                                         mass,
                                         barionic_density,
                                         proton_density,
                                         neutron_density);

            double potential =
                HadronThermodynamicPotential(scalar_density,
                                             barionic_density,
                                             proton_density,
                                             neutron_density,
                                             proton_chemical_potential,
                                             neutron_chemical_potential,
                                             kinectic_energy_density)
                - vacuum_potential;

            fprintf(file,
                    "%20.15E\t%20.15E\n",
                    (proton_chemical_potential + neutron_chemical_potential) / 2.0,
                    potential);

            fprintf(filep,
                    "%20.15E\t%20.15E\n",
                    (proton_chemical_potential + neutron_chemical_potential) / 2.0,
                    HadronPressure(potential));

                barionic_density += bar_dens_step;
            }


        fclose(file);

        SetParametersSet(NULL, NULL);
    }

    // Determine scalar density
    if(true)
    {
        SetParametersSet("PCP-0.0", NULL);
        SetFilePath ("tests/scalar_density");

        FILE * file = OpenFile ("scalar_density.dat");

        int points_number = 1000;

        double temperature = 0.0;
        double renormalized_chemical_potential = 0.0;

        double mass_min = 0.0;      // (MeV)
        double mass_max = 3.0E3;    // (MeV)

        double mass_step = Step(mass_min, mass_max, points_number);

        double mass = mass_min;
        for (int i = 0; i < points_number; i++){

            double scalar_density =
                QuarkScalarDensity(temperature,
                                   mass,
                                   renormalized_chemical_potential);

                fprintf(file,
                        "%20.15E\t%20.15E\n",
                        mass,
                        scalar_density);

                mass += mass_step;
            }

        fclose(file);

        SetParametersSet(NULL, NULL);
    }

    // Map zeroed gap equations
    //      To plot on gnuplot use
    //          set view map
    //          splot "quark_vacuum_mass_map.dat" w pm3d
    //      To make contour lines and adjust the number of lines
    //          set contour
    //          set cntrparam levels 35
    if (true)
    {
        SetParametersSet("PCP-0.0", NULL);
        SetFilePath ("tests/quark_vacuum_mass_maps");

        double min_mass = 0.0;
        double max_mass = 400.0;
        int n_pts = 100;

        FILE * file = OpenFile("quark_vacuum_mass_maps.dat");

        // Save a side view of f(m_u, m_d) with m_u = m_d
        FILE * equal_mass_view = OpenFile("equal_mass_view.dat");

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

                gsl_vector_set(x, 1, pow(up_mass, 1.0 / 2.0)); // same as above

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

        fclose(file);
        fclose(equal_mass_view);

        // Find common zero
        double up_vacuum_mass;
        double down_vacuum_mass;

        QuarkVacuumMassDetermination(&up_vacuum_mass, &down_vacuum_mass);

        file = OpenFile("quark_vacuum_mass_location.dat");

        fprintf (file,
                 "%20.15E\t%20.15E\n",
                 up_vacuum_mass,
                 down_vacuum_mass);

        fclose(file);

        SetParametersSet(NULL, NULL);
    }

    // Trying to see the dependence of each equation to be zeroed in simultaneous_solution
    // with each of the variables
    if (false)
    {
        SetParametersSet("PCP-0.0", "eNJL1");

        int n_pts = 1000;

        double min_dens = 0.0;
        double max_dens = 10.0;

        multi_dim_root_params params;
        params.temperature = 0.0;
        params.proton_fraction = 0.5;
        params.quark_vacuum_thermodynamic_potential = 0.0;
        params.hadron_vacuum_thermodynamic_potential = 0.0;

        double dens_step = Step (min_dens, max_dens, n_pts);

        gsl_vector * x = gsl_vector_alloc(4);
        gsl_vector * return_values = gsl_vector_alloc(4);

        gsl_vector_set(x, 1, 939.0); // hadron_mass
        gsl_vector_set(x, 2, 300.0); // up mass
        gsl_vector_set(x, 3, 300.0); // down mass

        SetFilePath ("tests/simultaneous_solution");
        FILE * file = OpenFile ("data.dat");

        // running on density
        double dens = min_dens;
        for (int i = 0; i < n_pts; i++){
            gsl_vector_set(x, 0, dens);
            MultiDimensionalRootFinderHelperFunction(x,
                                                     &params,
                                                     return_values);
            fprintf(file,
                    "%20.15E\t%20.15E\n",
                    dens,
                    gsl_vector_get(return_values, 0)); // we should verify the four
                                                       // equations for the four
                                                       // variables
            dens += dens_step;
        }

        SetParametersSet(NULL, NULL);
    }

/*    // Write functions to be zeroed in some way
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
                                        quark_vacuum_mass,
                                        0.0,
                                        0.0,
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
*/
/*    // Test pressures
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
                                        quark_vacuum_mass,
                                        0.0,
                                        0.0,
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
            double quark_thermodynamic_potential =
                QuarkThermodynamicPotential(up_quark_mass,
                                            down_quark_mass,
                                            up_chemical_potential,
                                            down_chemical_potential,
                                            up_renormalized_chemical_potential,
                                            down_renormalized_chemical_potential,
                                            temperature);

            double regularized_thermodynamic_potential = quark_thermodynamic_potential
                                                         - quark_vacuum_thermodynamic_potential;

            double quark_pressure = QuarkPressure(regularized_thermodynamic_potential,
                                                    parameters.variables.temperature);

            fprintf(hadron_pressure_file,
                    "%20.15E\t%20.15E\n",
                    barionic_density,
                    hadron_pressure);

            fprintf(quark_pressure_file,
                    "%20.15E\t%20.15E\n",
                    barionic_density,
                    quark_pressure);8

            barionic_density += barionic_density_step;
        }

        fclose(hadron_pressure_file);
        fclose(quark_pressure_file);

        SetFilePath(NULL);
    }
*/
    return;
}

typedef struct _TestMassAndRenormChemPot{
    double up_chemical_potential;
    double down_chemical_potential;

    double up_renorm_chem_pot;
    double down_renorm_chem_pot;
} TestMassAndRenormChemPot;

int TestMassAndRenormChemPotSimultaneousSolutionEquation(const gsl_vector   *x,
                                                         void *params,
                                                         gsl_vector *return_values);

int TestMassAndRenormChemPotSimultaneousSolution(double up_chemical_potential,
                                                 double down_chemical_potential,
                                                 double up_mass_guess,
                                                 double down_mass_guess,
                                                 double abs_error,
                                                 double rel_error,
                                                 int max_iter,
                                                 double * return_up_mass,
                                                 double * return_down_mass,
                                                 double * return_up_renorm_chem_pot,
                                                 double * return_down_renorm_chem_pot)
{
    // Set up parameters to be passed to helper function
    TestMassAndRenormChemPot params;
    params.up_chemical_potential = up_chemical_potential;
    params.down_chemical_potential = down_chemical_potential;

    // Set dimension (number of equations|variables to solve|find)
    const int dimension = 2;

    gsl_multiroot_function f;
    f.f = &TestMassAndRenormChemPotSimultaneousSolutionEquation;
    f.n = dimension;
    f.params = (void *)&params;

    gsl_vector * initial_guess = gsl_vector_alloc(dimension);
    gsl_vector * return_results = gsl_vector_alloc(dimension);

    gsl_vector_set(initial_guess,
                   0,
                   sqrt(up_mass_guess));
    gsl_vector_set(initial_guess,
                   1,
                   sqrt(down_mass_guess));

    int status =
        MultidimensionalRootFinder(dimension,
                                   &f,
                                   initial_guess,
                                   abs_error,
                                   rel_error,
                                   max_iter,
                                   return_results);

    if (status != 0){
        printf("%s:%d: Something is wrong with the rootfinding.\n",
               __FILE__,
               __LINE__);
        abort();
    }

    // Save results in return variables,
    // taking care of the mappinps
    *return_up_mass = pow(gsl_vector_get(return_results, 0), 2.0);
    *return_down_mass = pow(gsl_vector_get(return_results, 1), 2.0);

    *return_up_renorm_chem_pot = params.up_renorm_chem_pot;
    *return_down_renorm_chem_pot = params.down_renorm_chem_pot;

    // Free vectors
    gsl_vector_free(initial_guess);
    gsl_vector_free(return_results);

    return status;
}

int TestMassAndRenormChemPotSimultaneousSolutionEquation(const gsl_vector   *x,
                                                         void *params,
                                                         gsl_vector *return_values)
{
    const double up_mass = pow(gsl_vector_get(x, 0), 2.0);
    const double down_mass = pow(gsl_vector_get(x, 1), 2.0);

    TestMassAndRenormChemPot * p = (TestMassAndRenormChemPot *)params;

    double up_renorm_chem_pot;
    double down_renorm_chem_pot;

    QuarkSelfConsistentRenormalizedChemicalPotential(parameters.simultaneous_solution.renorm_chem_pot_solution,
                                                     up_mass,
                                                     down_mass,
                                                     p->up_chemical_potential,
                                                     p->down_chemical_potential,
                                                     parameters.variables.temperature,
                                                     &up_renorm_chem_pot,
                                                     &down_renorm_chem_pot);
    // save renormalized chemical potentials
    p->up_renorm_chem_pot = up_renorm_chem_pot;
    p->down_renorm_chem_pot = down_renorm_chem_pot;

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

    gsl_vector_set(return_values, 0, up_quark_zeroed_gap_eq);
    gsl_vector_set(return_values, 1, down_quark_zeroed_gap_eq);

    return GSL_SUCCESS;
}

typedef struct _BinodalParameters {
    double proton_fraction;
    double temperature;
} BinodalParameters;

double BinodalPointEquation(double  barionic_density,
                            void   *params);

void BinodalPoint(double temperature,
                  double proton_fraction,
                  double *return_barionic_chemical_potential,
                  double *return_isovector_chemical_potential,
                  double *return_pressure)
{

    BinodalParameters params;
    params.temperature = temperature;
    params.proton_fraction = proton_fraction;

    gsl_function F;
    F.function = &BinodalPointEquation;
    F.params = &params;

    double return_result;

    int status = UnidimensionalRootFinder(&F,
                                          parameters.binodal_rootfinding_params,
                                          &return_result);

    if (status != 0){
        printf("%s:%d: Problem with rootfinding.\n", __FILE__, __LINE__);
        abort();
    }

    // Determine return things

    double barionic_density = return_result;

    double vacuum_potential = HadronVacuumEnergyDensity();

    double proton_density = proton_fraction * barionic_density;
    double neutron_density = (1.0 - proton_fraction) * barionic_density;

    double proton_fermi_momentum = HPFermiMomentum(proton_density);
    double neutron_fermi_momentum = HPFermiMomentum(neutron_density);

    double mass = HadronSolveGapEquation(parameters.rootfinding_params,
                                         proton_density,
                                         neutron_density,
                                         proton_fermi_momentum,
                                         neutron_fermi_momentum);

    double kinectic_energy_density =
    HadronKinecticEnergyDensity(mass,
                                proton_fermi_momentum,
                                neutron_fermi_momentum);

    double scalar_density =
    HadronScalarDensity(mass,
                        proton_fermi_momentum,
                        parameters.hadron_model.cutoff)
    + HadronScalarDensity(mass,
                          neutron_fermi_momentum,
                          parameters.hadron_model.cutoff);

    // If barionic density == 0, the chemical
    // potential functions must return zero.
    double proton_chemical_potential =
    ProtonChemicalPotential(proton_fermi_momentum,
                            scalar_density,
                            mass,
                            barionic_density,
                            proton_density,
                            neutron_density);

    double neutron_chemical_potential =
    NeutronChemicalPotential(neutron_fermi_momentum,
                             scalar_density,
                             mass,
                             barionic_density,
                             proton_density,
                             neutron_density);

    double potential =
    HadronThermodynamicPotential(scalar_density,
                                 barionic_density,
                                 proton_density,
                                 neutron_density,
                                 proton_chemical_potential,
                                 neutron_chemical_potential,
                                 kinectic_energy_density)
    - vacuum_potential;

    double hadron_pressure = HadronPressure(potential);

    *return_barionic_chemical_potential = (proton_chemical_potential
                                           + neutron_chemical_potential) / 2.0;
    *return_isovector_chemical_potential = (proton_chemical_potential
                                            - neutron_chemical_potential);
    *return_pressure = hadron_pressure;

    printf("density: %f\n", barionic_density);
    printf("pressure: %f\n", hadron_pressure);
    printf("bar chemical potential: %f\n", *return_barionic_chemical_potential);
    printf("isovector chemical potential: %f\n", *return_isovector_chemical_potential);

    return;
}

double BinodalPointEquation(double  barionic_density,
                            void   *params)
{

    BinodalParameters * p = (BinodalParameters *)params;

    double proton_fraction = p->proton_fraction;
    double temperature = p->temperature;

    double hadron_vacuum_potential = HadronVacuumEnergyDensity();

    double proton_density = proton_fraction * barionic_density;
    double neutron_density = (1.0 - proton_fraction) * barionic_density;

    double proton_fermi_momentum = HPFermiMomentum(proton_density);
    double neutron_fermi_momentum = HPFermiMomentum(neutron_density);

    double mass = HadronSolveGapEquation(parameters.rootfinding_params,
                                         proton_density,
                                         neutron_density,
                                         proton_fermi_momentum,
                                         neutron_fermi_momentum);

    double kinectic_energy_density =
    HadronKinecticEnergyDensity(mass,
                                proton_fermi_momentum,
                                neutron_fermi_momentum);

    double scalar_density =
    HadronScalarDensity(mass,
                        proton_fermi_momentum,
                        parameters.hadron_model.cutoff)
    + HadronScalarDensity(mass,
                          neutron_fermi_momentum,
                          parameters.hadron_model.cutoff);

    // If barionic density == 0, the chemical
    // potential functions must return zero.
    double proton_chemical_potential =
    ProtonChemicalPotential(proton_fermi_momentum,
                            scalar_density,
                            mass,
                            barionic_density,
                            proton_density,
                            neutron_density);

    double neutron_chemical_potential =
    NeutronChemicalPotential(neutron_fermi_momentum,
                             scalar_density,
                             mass,
                             barionic_density,
                             proton_density,
                             neutron_density);

    double potential =
    HadronThermodynamicPotential(scalar_density,
                                 barionic_density,
                                 proton_density,
                                 neutron_density,
                                 proton_chemical_potential,
                                 neutron_chemical_potential,
                                 kinectic_energy_density)
    - hadron_vacuum_potential;

    double hadron_pressure = HadronPressure(potential);

    // From Gibbs conditions
    double up_chemical_potential = (2.0 * proton_chemical_potential
                                    - neutron_chemical_potential) / 3.0;

    double down_chemical_potential = (-proton_chemical_potential
                                      + 2.0 * neutron_chemical_potential) / 3.0;

    double return_up_mass = NAN;
    double return_down_mass = NAN;
    double return_up_renorm_chem_pot = NAN;
    double return_down_renorm_chem_pot = NAN;

    TestMassAndRenormChemPotSimultaneousSolution(up_chemical_potential,
                                                 down_chemical_potential,
                                                 parameters.other_rootfinding.up_mass_guess,
                                                 parameters.other_rootfinding.down_mass_guess,
                                                 parameters.other_rootfinding.abs_error,
                                                 parameters.other_rootfinding.rel_error,
                                                 parameters.other_rootfinding.max_iter,
                                                 &return_up_mass,
                                                 &return_down_mass,
                                                 &return_up_renorm_chem_pot,
                                                 &return_down_renorm_chem_pot);

    double up_vacuum_mass;
    double down_vacuum_mass;

    QuarkVacuumMassDetermination (&up_vacuum_mass, &down_vacuum_mass);

    double quark_vacuum_potential =
    QuarkThermodynamicPotential(up_vacuum_mass,
                                down_vacuum_mass,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0);

    double omega =
    QuarkThermodynamicPotential(return_up_mass,
                                return_down_mass,
                                up_chemical_potential,
                                down_chemical_potential,
                                return_up_renorm_chem_pot,
                                return_down_renorm_chem_pot,
                                temperature);

    double quark_pressure = QuarkPressure(omega - quark_vacuum_potential,
                                          temperature);

    return hadron_pressure - quark_pressure;
}
