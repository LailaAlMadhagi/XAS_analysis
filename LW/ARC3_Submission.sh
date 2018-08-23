#Submission script to run serial job
#name the job
#$ -N Liquid_samples 
#Email address for job notification
# -M fy11lham@leeds.ac.uk
#Send email when jobs begins, ends and aborts
# -m bae
#Use current working directory and export variables
#$ -cwd -V
#Specify memory
# -l h_vmem=1.0G
#Specify Wallclock time
#$ -l h_rt=05:00:00
/nobackup/fy11lham/miniconda3/bin/python /nobackup/fy11lham/XAS_Analysis/LW/LW.py /nobackup/fy11lham/XAS_Analysis/LW/L_args_Linux.txt
