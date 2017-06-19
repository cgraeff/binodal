//
//  Parameters.h
//  NJLv
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#ifndef Parameters_h
#define Parameters_h

#include <stdbool.h>
#include <stdio.h>

#include "libdatafun/libdatafun.h"

#include "HadronPhaseEOS.h"
#include "QuarkPhaseEOS.h"
#include "SimultaneousSolution.h"
#include "FermiDiracDistributions.h"
#include "Loop.h"

// Parameterization variables
typedef struct _QuarkModelParameters{
    char * parameters_set_identifier;
    char * parameters_set_origin;   // Where the set was taken from
    double G_S;                     // scalar-isoscalar coupling (fm^2)
    double G_V;                     // vector-isoscalar coupling (fm^2)
    double cutoff;                  // \Lambda (MeV)
    double bare_mass;               // (MeV)
} QuarkModelParameters;

// Hadron phase parameters
typedef struct _HadronModelParameters{
    char * parameters_set_identifier;
    char * parameters_set_origin;

    double G_S;				// scalar-isoscalar coupling (fm^2)
	double G_V;				// vector-isoscalar coupling (fm^2)
    double G_RHO;			// vector-isovector coupling (fm^2)
    double G_SV;			// scalar-isovector coupling (fm^8)
    double G_VRHO;			// (fm^8)
    double G_SRHO;			// (fm^8)

    double cutoff;			// (MeV)
    double nucleon_mass;	// (MeV)
    double bare_mass;		// (MeV)
} HadronModelParameters;

// Parameters for the variable for which we are
// performing the calculation of the EOS
typedef struct _VariableParameters{
    double temperature; // (MeV)
    double min_proton_fraction;		// (no dimension)
    double max_proton_fraction;
    int num_points;
} VariableParameters;

typedef struct _QuarkVacuumMassDeterminationParameters{
    double up_vacuum_mass_guess;
    double down_vacuum_mass_guess;

    int max_iter;

    double abs_error;
    double rel_error;
} QuarkVacuumMassDeterminationParameters;

typedef struct _OtherRootfindingParameters{
    double up_mass_guess;
    double down_mass_guess;
    double abs_error;
    double rel_error;
    int max_iter;
} OtherRootFindingParameters;

typedef struct _parameters
{
    VariableParameters variables;
    HadronModelParameters hadron_model;
    QuarkModelParameters quark_model;

    QuarkVacuumMassDeterminationParameters vacuum_mass_determination;
    UnidimensionalRootFindingParameters q_renorm_chem_pot_finding;
    SimultaneousSolutionParameters simultaneous_solution;

    IntegratorParameters fermi_dirac_integrals;
    IntegratorParameters therm_pot_free_gas_integral;

    //my tests
    UnidimensionalRootFindingParameters rootfinding_params;
    OtherRootFindingParameters other_rootfinding;
    UnidimensionalRootFindingParameters binodal_rootfinding_params;

} Parameters;

extern Parameters parameters;

void ParametersSetsSetup(void);
void SetParametersSet(char *quark_set_identifier,
                      char *hadron_set_identifier);

void PrintParametersToFile(FILE * file);

#endif /* Parameters_h */
