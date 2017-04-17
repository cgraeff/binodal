
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

double QuarkSelfConsistentRenormChemPot(double quark_mass,
                                        double chemical_potential,
                                        double temperature);



void SimultaneousSolution(SimultaneousSolutionParameters params,
                          double quark_vacuum_thermodynamic_potential,
                          double vacuum_thermodynamic_potential,
                          double *return_barionic_density,
                          double *return_hadron_mass,
                          double *return_quark_mass,
                          double *return_proton_chemical_potential,
                          double *return_neutron_chemical_potential,
                          double *return_pressure)
{
    // Set up parameters to be passed to helper function
    multi_dim_root_params p;
    p.temperature = parameters.variables.temperature;
    p.proton_fraction = parameters.variables.proton_fraction;

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
    *return_barionic_density = pow(gsl_vector_get(return_results, 0), 2.0); // TODO return this?
    *return_hadron_mass = pow(gsl_vector_get(return_results, 1), 2.0);
    *return_quark_mass = pow(gsl_vector_get(return_results, 2), 2.0);

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

    printf("bar_dens: %f\n"
           "hadron_mass: %f\n"
           "up_q_mass: %f\n"
           "dn_q_mass: %f\n",
           barionic_density,
           hadron_mass,
           up_quark_mass,
           down_quark_mass);

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

    // Set up parameters for gap equation
	hadron_gap_eq_input_params h_gap_input;
    h_gap_input.proton_fermi_momentum = proton_fermi_momentum;
    h_gap_input.neutron_fermi_momentum = neutron_fermi_momentum;
    h_gap_input.proton_density = proton_density;
    h_gap_input.neutron_density = neutron_density;

    double zeroed_gap_eq = HadronGapEquation(hadron_mass, &h_gap_input);

    // Prepare return vector
   	gsl_vector_set(return_values, 0, zeroed_gap_eq);

   	// Gibbs Conditions: // TODO verify: are those up and down barionic chemical potentials? should they be?
    // Use chemical potential equality to determine
    // quark chemical potentials
    double up_chemical_potential = (2.0 * proton_chemical_potential
                                    - neutron_chemical_potential) / 3.0;

    double down_chemical_potential = (-proton_chemical_potential
                                      + 2.0 * neutron_chemical_potential) / 3.0;

   	// Quarks:
    quark_gap_eq_input_params q_gap_input;

    // up quark
    double up_renormalized_chemical_potential =
        QuarkSelfConsistentRenormChemPot(up_quark_mass, up_chemical_potential, p->temperature);

    q_gap_input.renormalized_chemical_potential = up_renormalized_chemical_potential;
    q_gap_input.temperature = p->temperature;

    double up_quark_zeroed_gap_eq = QuarkZeroedGapEquation(up_quark_mass, &q_gap_input);

    gsl_vector_set(return_values, 1, up_quark_zeroed_gap_eq);

    // down quark
    double down_renormalized_chemical_potential =
        QuarkSelfConsistentRenormChemPot(down_quark_mass, down_chemical_potential, p->temperature);

    q_gap_input.renormalized_chemical_potential = down_renormalized_chemical_potential;
    q_gap_input.temperature = p->temperature;

    double down_quark_zeroed_gap_eq = QuarkZeroedGapEquation(down_quark_mass, &q_gap_input);

    gsl_vector_set(return_values, 2, down_quark_zeroed_gap_eq);

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
            QuarkThermodynamicPotential(up_quark_mass,
                                        up_chemical_potential,
                                        up_renormalized_chemical_potential,
                                        parameters.variables.temperature);

    double regularized_thermodynamic_potential = up_quark_thermodynamic_potential
                                                 + down_quark_thermodynamic_potential
                                                 - p->quark_vacuum_thermodynamic_potential;

    double quark_pressure = QuarkPressure(regularized_thermodynamic_potential,
                                            parameters.variables.temperature);

    double zeroed_pressure = hadron_pressure - quark_pressure;

    gsl_vector_set(return_values, 3, zeroed_pressure);

    printf("eq 1: %f\n"
           "eq 2: %f\n"
           "eq 3: %f\n"
           "eq 4: %f\n",
           gsl_vector_get(return_values, 0),
           gsl_vector_get(return_values, 1),
           gsl_vector_get(return_values, 2),
           gsl_vector_get(return_values, 3));

    return GSL_SUCCESS;
}

double QuarkZeroedGapEquation(double mass,
                         void * params)
{
    quark_gap_eq_input_params * p = (quark_gap_eq_input_params *)params;

    double term = 2.0 * parameters.quark_model.G_S * CONST_HBAR_C
                  * QuarkScalarDensity(parameters.variables.temperature,
                                  mass,
                                  p->renormalized_chemical_potential);

    return mass - parameters.quark_model.bare_mass + term;
}

double HadronZeroedGapEquation(double mass,
                               void * params)
{
    hadron_gap_eq_input_params * p = (hadron_gap_eq_input_params *)params;

    double barionic_density = p->proton_density + p->neutron_density;
    double rho_3 = p->proton_density - p->neutron_density;

	double scalar_density = HadronScalarDensity(mass,
                                                p->proton_fermi_momentum,
                                                parameters.hadron_model.cutoff)
                            + HadronScalarDensity(mass,
                                                  p->neutron_fermi_momentum,
                                                  parameters.hadron_model.cutoff);

	double gap_1st_term = parameters.hadron_model.G_S * scalar_density;
	double gap_2nd_term = - parameters.hadron_model.G_SV
                            * scalar_density
                            * pow(barionic_density, 2.0);
    double gap_3rd_term = - parameters.hadron_model.G_SRHO
                            * scalar_density
                            * pow(rho_3, 2.0);

	return mass
           + 2.0 * CONST_HBAR_C
             * (gap_1st_term + gap_2nd_term + gap_3rd_term)
           - parameters.hadron_model.bare_mass;
}

double QuarkSelfConsistentRenormChemPot(double quark_mass,
                                        double chemical_potential,
                                        double temperature)
{

    renorm_chem_pot_equation_input p;
    p.chemical_potential = chemical_potential;
    p.mass = quark_mass;
    p.temperature = temperature;

    gsl_function F;
    F.params = (void *)&p;
    F.function = &ZeroedRenormalizedChemicalPotentialEquation;

    double result = 0;

    int status = UnidimensionalRootFinder(&F,
                                          parameters.q_renorm_chem_pot_finding,
                                          &result);

    if (status != 0){
        printf("Problems in rootfinding in %s, line %d\n", __FILE__, __LINE__);
        abort();
    }

    return result;
}

double ZeroedRenormalizedChemicalPotentialEquation(double renormalized_chemical_potential,
                                                   void * params)
{
    renorm_chem_pot_equation_input * p = (renorm_chem_pot_equation_input *)params;

    double term = 2.0 * parameters.quark_model.G_V * NUM_Q_COLORS * CONST_HBAR_C
                  * QuarkBarionicDensity(p->mass,
                                         renormalized_chemical_potential,
                                         p->temperature);

    return renormalized_chemical_potential - p->chemical_potential + term;
}
