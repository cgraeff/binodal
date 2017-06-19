
#ifndef SimultaneousSolution_h
#define SimultaneousSolution_h

typedef struct _RenormChemPotSolutionParameters{
    double up_renorm_chem_pot_guess;
    double down_renorm_chem_pot_guess;

    double abs_error;
    double rel_error;
    double max_iter;
} RenormChemPotSolutionParameters;

typedef struct _SimultaneousSolutionParameters{
    double temperature;
    double barionic_density_guess;
    double hadron_mass_guess;
    double up_quark_mass_guess;
    double down_quark_mass_guess;

    RenormChemPotSolutionParameters renorm_chem_pot_solution;

    int max_iter;

    double abs_error;
    double rel_error;
} SimultaneousSolutionParameters;

typedef struct _multi_dim_root_params{
    double temperature;
    double proton_fraction;
    double quark_vacuum_thermodynamic_potential;
    double hadron_vacuum_thermodynamic_potential;
} multi_dim_root_params;

typedef struct _renorm_chem_pot_equation_input{
    double up_quark_mass;
    double down_quark_mass;
    double up_chemical_potential;
    double down_chemical_potential;
	double temperature;
} renorm_chem_pot_equation_input;

void SimultaneousSolution(SimultaneousSolutionParameters params,
                          double quark_vacuum_thermodynamic_potential,
                          double vacuum_thermodynamic_potential,
                          double proton_fraction,
                          double *return_barionic_density,
                          double *return_hadron_mass,
                          double *return_up_quark_mass,
                          double *return_down_quark_mass);

int MultiDimensionalRootFinderHelperFunction(const gsl_vector   *x,
                                             void               *params,
                                             gsl_vector         *return_values);

#endif /* SimultaneousSolution_h */

