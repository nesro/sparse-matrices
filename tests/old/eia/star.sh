#!/bin/bash

if [[ "$1" == "mmm" ]]; then
	c=m
elif [[ "$1" == "mvm" ]]; then
	c=v
else
	echo Invalid arument!
	exit 1
fi

cat >> $c$$.star.sh << __EOF__
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
#$ -e /mnt/data/nesrotom/ 

# The path used for the standard output stream of the job.
#$ -o /mnt/data/nesrotom/ 

# Do not change.
#$ -pe ompi 1

# Do not change.
#$ -q 12c_1slots_per_host.q

time ./tests/eia/eia.sh $1 $2 $3 $4 $5
echo "^TIME of ./tests/eia/eia.sh $1 $2 $3 $4 $5"

__EOF__

chmod +x $c$$.star.sh
/opt/bin/qrun.sh 12c 1 1slots_per_host $c$$.star.sh

