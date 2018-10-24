#loop over all clustering results files and run validation/GO enrichment code on each
#do per network inference method
for file in "$@"
do
num=$(echo $file | sed 's/consensus_inflation//')
inflation=$(echo $num | sed -E 's/([0-9]{1})([0-9]{1})/\1\.\2/')
python main_FDR_V1.py -in consensus_edges_weights.txt -o cons -ig PID_6_10_NEW.csv -n 3 -a 0.05 -ic ../CRC_MCL_clustering/consensus_clustering_results/$file -d consensus_results/inflation$num
done
