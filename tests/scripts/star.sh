#!/bin/bash
#
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
# https://github.com/nesro/sparse-matrices
#
# This is a helper file for running jobs on our school server STAR.
# Please do NOT to try hack this in any way.

cat >> _$$.sh << __EOF__
#!/bin/sh

#  ===========================================================================
# |                                                                           |
# |             COMMAND FILE FOR SUBMITTING SGE JOBS                          |
# |                                                                           |
# |                                                                           |
# | SGE keyword statements begin with #$                                      |
# |                                                                           |
# | Comments begin with #                                                     |
# | Any line whose first non-blank character is a pound sign (#)              |
# | and is not a SGE keyword statement is regarded as a comment.              |
#  ===========================================================================

# Request Bourne shell as shell for job
#$ -S /bin/sh

# Execute the job from the current working directory.
#$ -cwd

# Defines  or  redefines  the  path used for the standard error stream of the job.
#$ -e /mnt/data/nesrotom/sparse-matrices

# The path used for the standard output stream of the job.
#$ -o /mnt/data/nesrotom/sparse-matrices

# Do not change.
#$ -pe ompi 1

# Do not change.
#$ -q 12c_1slots_per_host.q

time ./tests/scripts/star_job.sh

__EOF__

chmod +x _$$.sh
/opt/bin/qrun.sh 12c 1 1slots_per_host _$$.sh