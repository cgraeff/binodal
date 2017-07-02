//
//  Parameters.h
//  binodal
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#ifndef Parameters_h
#define Parameters_h

#include "libdatafun/libdatafun.h"

#include "HadronPhaseEOS.h"
#include "QuarkPhaseEOS.h"

// Parameters for the variables for which we are
// performing the calculation of the EOS
typedef struct _VariableParameters{
    double temperature;             // (MeV)
    double min_proton_fraction;     // (no dimension)
    double max_proton_fraction;
    int num_points;
} VariableParameters;

typedef struct _parameters
{
    VariableParameters variables;

    struct _hadron {
        HadronModelParameters model;
        UnidimensionalRootFindingParameters gap_eq_solution_params;
    } hadron;

    struct _quark {
        QuarkModelParameters model;

        QuarkVacuumMassDeterminationParameters vacuum_mass_determination;
        QuarkRenormChemPotSolutionParameters renorm_chem_pot_solution;

        QuarkMassAndRenormChemPotSolParams mass_and_renorm_chem_pot_solution;

        QuarkThermodynamicPotMinDetermParams therm_pot_minimum;

        QuarkEntropyIntegrationParameters entropy_integration_params;
    } quark;

    IntegratorParameters fermi_dirac_integrals;
    IntegratorParameters therm_pot_free_gas_integral;

    UnidimensionalRootFindingParameters binodal_rootfinding_params;

} Parameters;

extern Parameters parameters;

void ParseModelParametersSets(void);
void SetParametersSet(char *quark_set_identifier,
                      char *hadron_set_identifier);

#endif /* Parameters_h */
