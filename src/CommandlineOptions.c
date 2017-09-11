//
//  CommandlineOptions.c
//  binodal
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#include <stddef.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <getopt.h>

#include "CommandlineOptions.h"

// Default values for options and flags that will be acessible
// during the execution (specified in order of declaration).
Options options = {true, false, false, false, false, NULL, NULL, -1.0, -1.0};

int CommandlineOptionsParse(int argc, char * argv[])
{
    // Short options must be declared in the next variable.
    // In short_options, each short option is declared in a string by the character that
    // will invoke the option. If the option takes an argument, a colon (:) must be
    // follow the character. When a character contained in this string is found in the
    // options, the number that corresponds to this caracter is returned by getopt()
    // (that is "a" -> 'a', ...).
    // EXAMPLE:
    // char * short_options = "a:bvuh";
    //
    // BEWARE: if an option takes many arguments, the spaces must be escaped, otherwise
    // the arguments after the first will be misinterpreted as unknown, or unclaimed.
    // This particular implementation will stop if there are any unprocessed arguments.

    char * short_options = "rq:h:t:lsdau";

    int opt;
    while ((opt = getopt(argc, argv, short_options)) != -1){

        // If an option have an argument, it is accessed through 'optarg'
        switch (opt){
            case 'r':
                options.solution_for_barionic_chemical_potential_range = true;
                break;
            case 'q':
                options.quark_parameterization = optarg;
                break;
            case 'h':
                options.hadron_parameterization = optarg;
                break;
            case 't':
                options.temp = atof(optarg);
                break;
            case 'y':
                options.proton_fraction = atof(optarg);
                break;
            case 'l':
                options.list_available_parameterizations = true;
                break;
            case 's':
                options.verbose = false;
                break;
            case 'd':
                options.dirs = true;
                break;
            case 'a':
                options.dirs = true;
                options.tests = true;
                break;
            case 'u':
                CommandlineOptionsPrintUsage();
                exit(EXIT_SUCCESS);
                break;
            case '?':
                if (optopt == 'p'){
                    fprintf (stderr,
                             "Option -%c requires an argument. Use -u for help.\n",
                             optopt);
                }else if (isprint (optopt)){
                    fprintf (stderr,
                             "Unknown option `-%c'.  Use -u for help.\n",
                             optopt);
                }else{
                    fprintf (stderr,
                             "Unknown option character `\\x%x'.  Use -u for help.\n",
                             optopt);
                }
                exit(EXIT_FAILURE);
                break;
            default:
                printf("Use -u to print a list of accepted options.\n");
                exit(EXIT_FAILURE);
        }
    }


    // Print any remaining command line arguments (not options)
    // and let the user know that they are invalid
    if (optind < argc){
        printf ("%s: invalid arguments -- ", argv[0]);
        while (optind < argc)
            printf ("%s ", argv[optind++]);
        putchar ('\n');
        printf("Use -u to print a list of accepted options.\n");
        exit(EXIT_FAILURE);
    }

    return 0;
}

void CommandlineOptionsPrintUsage()
{
    printf("This code calculates the binodal "
           "using NJL and eNJL models\n\n");
    printf("Usage: eos [options]\n");
    printf("Options:\n"
           "\t-r: Solve for range of barionic chemical potential;\n"
           "\t-q par: Chooses a builtin quark parameterization;\n"
           "\t-h par: Chooses a builtin hadron parameterization;\n"
           "\t-t temp: Chooses a temperature;\n"
           "\t-y frac: Chooses a proton fraction;\n"
           "\t-l: Lists available builtin parameterizations;\n"
           "\t-s: silent (supress information written to standard out);"
           "\t-d: write results using a dir structure;\n"
           "\t-a: run tests. Automatically sets -d;\n"
           "\t-u: Prints this message;\n\n");
    printf("The source code is available at github.com/cgraeff/binodal\n");

    return;
}

