#! /bin/bash
#$ -q long
#$ -pe smp 1
#$ -t 0-26:1

#run permutation tests for random forest network inference
#run in chunks - 9 hr run time, split into reasonable number of nodes
cd ~/ferdig_lab
module load python

#want to say... for ii in (${SGE_TASK_ID})10 to (${SGE_TASK_ID})10+9...
#for ((i=(${SGE_TASK_ID})10));i<=(${SGE_TASK_ID})10+9;i++));
run=${SGE_TASK_ID}
start=$(echo "$((run * 20))")
end=$(echo "$((run * 20 + 9))")
for (( file_num=$start; file_num<=$end; file_num++ ))
do
python rf_aracne_perms.py -i ~/ferdig_lab/perm_nets/hu_perm_"$file_num".txt -o ~/ferdig_lab/perm_nets/sampling_reults_perm_"$file_num".txt
done
