//
//  main.c
//  binodal
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_errno.h>

#include "CommandlineOptions.h"
#include "Parameters.h"
#include "Tests.h"
#include "Loop.h"

int main(int argc, char * argv[])
{
    // TODO: Check why we need that:
    // initialize random number generator
    srand(time(NULL));

    // Turn off abort-on-error for GSL
    gsl_set_error_handler_off();

    CommandlineOptionsParse(argc, argv);
    ParseModelParametersSets();

    if (options.tests){
        RunTests();
        exit(EXIT_SUCCESS);
    }

    // If options -q and/or -h are invoked when the program is run,
    // set quark and/or hadron parameter-set accordingly.
    // If those options are not used, 'options.quark_parameterization'
    // and/or 'options.hadron_parameterization' will be NULL, and that
    // corresponds to the default parameter-set, which is the first to
    // be declared in 'Parameters.c'
    SetParametersSet(options.quark_parameterization,
                     options.hadron_parameterization);

    // If the temperature was chosen using
    // commandline options, use it.
    // (-1.0 is a place holder value)
    if (options.temp != -1.0){
        if (options.temp >= 0.0){
            parameters.variables.temperature = options.temp;
        }
        else{
            printf("Values of temperature must be non-negative.\n");
            printf("(%f was provided).\n",
                   options.temp);
            exit(EXIT_FAILURE);
        }
    }

    // The solution may be TODO: verify how exactly this solution works
    if (options.solution_for_barionic_chemical_potential_range){

        SolveBinodalForBarionicChemicalPotentialRange();

        return 0;
    }

    SolveBinodalForBarionicAndIsovectorChemicalPotentialsGrid();

    return 0;
}

