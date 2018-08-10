//
//  CommandlineOptions.h
//  binodal
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#ifndef CommandlineOptions_h
#define CommandlineOptions_h

#include <stdbool.h>

typedef struct _options{
    // List options and flags that will be acessible
    // during the execution
    bool abort_on_error;
    bool verbose;
    bool debug;
    bool dirs;
    bool tests;
    bool list_available_parameterizations;
    bool solution_for_barionic_chemical_potential_range;
    char * quark_parameterization;
    char * hadron_parameterization;
    double temp;
    double proton_fraction;
} Options;

extern Options options;

int CommandlineOptionsParse(int argc, char * argv[]);
void CommandlineOptionsPrintUsage();

#endif /* CommandlineOptions.h */
