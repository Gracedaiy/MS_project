#!/bin/sh

# see where the job is being run
hostname

# print date and time
date

module load R

Rscript try.R ${SGE_TASK_ID}
