

perl find_consistent_TADs_TALL.pl test_data/TALL-ALL116_domains.k\=001.bed test_data/T_cell-donor1_domains.k\=001.bed 3 40000 test_results/TALL-ALL16-vs-T_cell1_TALL-TADs.tsv test_results/TALL-ALL16-vs-T_cell1_Tcell-TADs.tsv

./run_comparison.sh test_data/TALL-ALL116/ test_data/T_cell-donor1/ test_results/TALL-ALL16-vs-T_cell1_TALL-TADs.tsv test_results/TALL-ALL16-vs-T_cell1_Tcell-TADs.tsv test_results/results

R --no-save test_results/results/final_results.tsv gene-name.tsv FALSE FALSE \
	"Name_Sample_1" "Name_Sample_2" 40000 400000 test_results/results/final_results < differential_tad_activity_expression.r

