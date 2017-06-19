
#include "libdatafun/libdatafun.h"

#include "SimultaneousSolution.h"

#include "Constants.h"

#include "Parameters.h"
#include "QuarkPhaseEOS.h"
#include "HadronPhaseEOS.h"
#include "Functions.h"

int MultiDimensionalRootFinderHelperFunction(const gsl_vector   *x,
                                             void               *params,
                                             gsl_vector         *return_values);

void SimultaneousSolution(SimultaneousSolutionParameters params,
                          double quark_vacuum_thermodynamic_potential,
                          double hadron_vacuum_thermodynamic_potential,
                          double proton_fraction,
                          double *return_barionic_density,
                          double *return_hadron_mass,
                          double *return_up_quark_mass,
                          double *return_down_quark_mass)
{
    // Set up parameters to be passed to helper function
    multi_dim_root_params p;
    p.temperature = parameters.variables.temperature;
    p.proton_fraction = proton_fraction;
    p.quark_vacuum_thermodynamic_potential = quark_vacuum_thermodynamic_potential;
    p.hadron_vacuum_thermodynamic_potential = hadron_vacuum_thermodynamic_potential;

    // Set dimension (number of equations|variables to solve|find)
    const int dimension = 4;

    gsl_multiroot_function f;
    f.f = &MultiDimensionalRootFinderHelperFunction;
    f.n = dimension;
    f.params = (void *)&p;

    gsl_vector * initial_guess = gsl_vector_alloc(dimension);
    gsl_vector * return_results = gsl_vector_alloc(dimension);

    gsl_vector_set(initial_guess, 0, sqrt(params.barionic_density_guess));
    gsl_vector_set(initial_guess, 1, sqrt(params.hadron_mass_guess));
    gsl_vector_set(initial_guess, 2, sqrt(params.up_quark_mass_guess));
    gsl_vector_set(initial_guess, 3, sqrt(params.down_quark_mass_guess));

    int status =
        MultidimensionalRootFinder(dimension,
                                   &f,
                                   initial_guess,
                                   parameters.simultaneous_solution.abs_error,
                                   parameters.simultaneous_solution.rel_error,
                                   parameters.simultaneous_solution.max_iter,
                                   return_results);

    if (status != 0){
        printf("%s:%d: Something is wrong with the rootfinding.\n",
               __FILE__,
               __LINE__);
        abort();
    }

    // Save results in return variables,
    // taking care of the mappinps
    *return_barionic_density = pow(gsl_vector_get(return_results, 0), 2.0);
    *return_hadron_mass = pow(gsl_vector_get(return_results, 1), 2.0);
    *return_up_quark_mass = pow(gsl_vector_get(return_results, 2), 2.0);
    *return_down_quark_mass = pow(gsl_vector_get(return_results, 3), 2.0);

    // Free vectors
    gsl_vector_free(initial_guess);
    gsl_vector_free(return_results);

    return;
}

int MultiDimensionalRootFinderHelperFunction(const gsl_vector   *x,
                                             void               *params,
                                             gsl_vector         *return_values)
{
    multi_dim_root_params *p = (multi_dim_root_params *)params;

    double proton_fraction = p->proton_fraction;

    // The parameters will be passed to the subfunctions Zeroed*

    // Mappings:
    //      The variables for the root finding are assumed to cover the range
    //      (-\infty, +\infty), but that is not the case for the variables
    //      that we are trying to solve.
    //      To solve that, we use the mappings:
    //          $\rho_B = x^2$
    //          $m_H = x^2$
    //          $m_Q = x^2$
    //      The initial guesses must be transformed by inverting the relations
    //      above

    const double barionic_density = pow(gsl_vector_get(x, 0), 2.0);
    const double hadron_mass = pow(gsl_vector_get(x, 1), 2.0);
    const double up_quark_mass = pow(gsl_vector_get(x, 2), 2.0);
    const double down_quark_mass = pow(gsl_vector_get(x, 3), 2.0);
/*
    printf("%s:%d\n"
           "\tRHO: %f\n"
           "\tM_H: %f\n"
           "\tUP_M: %f\n"
           "\tDN_M: %f\n",
           __FILE__,
           __LINE__,
           barionic_density,
           hadron_mass,
           up_quark_mass,
           down_quark_mass);
*/
    /*
     *  Hadrons
     */

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

    // Set up parameters for gap equation
	hadron_gap_eq_input_params h_gap_input;
    h_gap_input.proton_fermi_momentum = proton_fermi_momentum;
    h_gap_input.neutron_fermi_momentum = neutron_fermi_momentum;
    h_gap_input.proton_density = proton_density;
    h_gap_input.neutron_density = neutron_density;

    double zeroed_gap_eq = HadronZeroedGapEquation(hadron_mass, &h_gap_input);

    // Prepare return vector
   	gsl_vector_set(return_values, 0, zeroed_gap_eq);

    /*
     * Gibbs Conditions
     */

    // Use chemical potential equality to determine
    // quark chemical potentials
    // TODO: Verify that this is correct
    double up_chemical_potential = (2.0 * proton_chemical_potential
                                    - neutron_chemical_potential) / 3.0;

    double down_chemical_potential = (-proton_chemical_potential
                                      + 2.0 * neutron_chemical_potential) / 3.0;

    /*
     *  Quarks
     */

    // Self consistent determination of renormalized
    // chemical potentials

    double up_renorm_chem_pot;
    double down_renorm_chem_pot;

    QuarkSelfConsistentRenormalizedChemicalPotential(parameters.simultaneous_solution.renorm_chem_pot_solution,
                                                      up_quark_mass,
                                                      down_quark_mass,
                                                      up_chemical_potential,
                                                      down_chemical_potential,
                                                      parameters.variables.temperature,
                                                      &up_renorm_chem_pot,
                                                      &down_renorm_chem_pot);

    // Gap equations:
    double up_scalar_density =
        QuarkScalarDensity(parameters.variables.temperature,
                           up_quark_mass,
                           up_renorm_chem_pot);

    double down_scalar_density =
        QuarkScalarDensity(parameters.variables.temperature,
                           down_quark_mass,
                           down_renorm_chem_pot);

    double up_quark_zeroed_gap_eq =
        QuarkZeroedGapEquation(up_quark_mass,
                               up_scalar_density,
                               down_scalar_density);

    double down_quark_zeroed_gap_eq =
        QuarkZeroedGapEquation(down_quark_mass,
                               up_scalar_density,
                               down_scalar_density);

    gsl_vector_set(return_values, 1, up_quark_zeroed_gap_eq);
    gsl_vector_set(return_values, 2, down_quark_zeroed_gap_eq);

    /*
     * Gibbs' conditions for pressure
     */

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

    double hadron_pressure = HadronPressure(hadron_thermodynamic_potential
                                            - p->hadron_vacuum_thermodynamic_potential);

    // Determination o quark pressure
    double quark_thermodynamic_potential =
        QuarkThermodynamicPotential(up_quark_mass,
                                    down_quark_mass,
                                    up_chemical_potential,
                                    down_chemical_potential,
                                    up_renorm_chem_pot,
                                    down_renorm_chem_pot,
                                    parameters.variables.temperature);

    double regularized_thermodynamic_potential = quark_thermodynamic_potential
                                                 - p->quark_vacuum_thermodynamic_potential;

    double quark_pressure = QuarkPressure(regularized_thermodynamic_potential,
                                            parameters.variables.temperature);

    double zeroed_pressure = hadron_pressure - quark_pressure;

//    printf("\tPRESSURE: %f\n", hadron_pressure);
//    printf("\tMU_B: %f\n", (proton_chemical_potential + neutron_chemical_potential) * 0.5);

    gsl_vector_set(return_values, 3, zeroed_pressure);

    return GSL_SUCCESS;
}

