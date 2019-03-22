#!/bin/bash

### Job name
#PBS -N binodal

### Email on abort (a), begins (b), end (e)
#PBS -M cgraeff@gmail.com
#PBS -m abe

### Set max wallclock time
#PBS -l walltime=100:00:00

### Save stderr and stdout
#PBS -o output.txt
#PBS -e error.txt

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

cd $PBS_O_WORKDIR

./binodal