Some bash commands that help with pipeline management

Find pipeline steps that don't have a 'results' directory

`analysis_dir/pipeline$ find . -maxdepth 1 -type d ! -exec test -e "{}/results" \; -print`

Find and delete all 'errors' directories (be careful with this!)

`analysis_dir/pipeline$ find . -mindepth 2  -type d -name "errors" -exec rm -rf {} \;`

Find and delete all 'results' directories except for the alignments (be carefule with this!)

`analysis_dir/pipeline$ find . -type d -name "results" ! -path "*align*" -exec rm -rf {} \;`

Find files from specific branch of a results directory and run a custom script on them

`find "$path_to_results_dir" -path "*/peaks.by_sample.macs_broad/*" -name "venn.txt" -print0 | xargs -0 $custom_script`
