#!/bin/bash

treecmp=$HOME/Downloads/TreeCmp/build/libs/TreeCmp-2.0-b103.jar
cwd=$(pwd)
newick_data=$cwd/trees
out_dir=$cwd/treecmp_out
log=$out_dir/log
tmp=$out_dir/tmp
meta=$cwd/reports/sampled_cds.csv

mkdir -p $log
mkdir -p $tmp

# Generate a csv file to store the results
echo "Threshold,Sample,RefTree,Tree,RefTree_taxa,Tree_taxa,Common_taxa,RF(0.5)"	> $out_dir/treecmp_results.csv
echo "Threshold,Method,Avg,Std,Min,Max,N_Subsets,N_Genes" > $out_dir/treecmp_summary.csv

for i in {1..190}; do
	j=$(printf	"%03d" $i)
	j=$(echo "scale=3; $j/1000" | bc)
	j=$(echo $j | sed 's/0\{1,\}$//')
	if [ -f "$newick_data/0$j.newick" ]; then
		echo "Running TreeCMP on 0$j.newick" > $log/0$j.log
		cat	"$newick_data/0$j"_{1..6}".newick" > $tmp/catsamples_0$j.newick
		java -jar $treecmp -r $newick_data/0$j.newick -i $tmp/catsamples_0$j.newick -d rf -o $out_dir/0$j.out -P -I >> $log/0$j.log

		#	Extract the results
		sed -i"" -e 's/\t/,/g' $out_dir/0$j.out
		awk '/^[0-9]/{print "0'$j'," $0}' $out_dir/0$j.out >> $out_dir/treecmp_results.csv
		awk '/^[RF]/{print "0'$j',"	$0}' $out_dir/0$j.out > $tmp/tmp.csv

		# Count the number of genes that fall within the threshold
		count=$(awk -F, -v thresh=0$j 'NR==1 {for (i=1; i<=NF; i++) {if ($i=="snps_per_base") col=i}} NR>1 {if (col && sprintf("%.2f",$col) >= thresh) count++} END {if (count) printf count; else printf "0"}' $meta)
		echo "$(cat $tmp/tmp.csv),$count" >> $out_dir/treecmp_summary.csv

	else
		continue
		fi
done

# Cleanup
rm -rf $tmp
rm -rf $out_dir/0*.out
rm -rf $out_dir/0*.out-e
