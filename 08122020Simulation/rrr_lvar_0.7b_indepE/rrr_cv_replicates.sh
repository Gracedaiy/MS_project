#!/bin/sh

# see where the job is being run
hostname

# print date and time
date

module load R

Rscript rrr_cv_replicates.R ${SGE_TASK_ID}
