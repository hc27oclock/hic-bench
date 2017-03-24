Some bash commands that help with pipeline management

Find pipeline steps that don't have a 'results' directory

`analysis_dir/pipeline$ find . -maxdepth 1 -type d ! -exec test -e "{}/results" \; -print`

Find and delete all 'errors' directories (be careful with this!)

`analysis_dir/pipeline$ find . -mindepth 2  -type d -name "errors" -exec rm -rf {} \;`

Find and delete all 'results' directories except for the alignments (be careful with this!)

`analysis_dir/pipeline$ find . -type d -name "results" ! -path "*align*" -exec rm -rf {} \;`

Find files from specific branch of a results directory and run a custom script on them

`find "$path_to_results_dir" -path "*/peaks.by_sample.macs_broad/*" -name "venn.txt" -print0 | xargs -0 $custom_script`


Backup all `results` dirs to run the pipeline fresh, except for the alignments.

```bash
backup_results () {
    local results_dir="$1"
    (
    cd "$(dirname "$results_dir")"
    code/file_backup.sh "$results_dir" old
    )
}

# run from the 'pipeline' dir
pipeline$ find -maxdepth 2 -type d -name "results" ! -path "*align/*" | while read item; do
    backup_results "$item"
done


```



Get rid of all the `_Sxx` entries in the inputs sample directories

```bash
cd inputs/fastq
bad_dirs="$(find . -maxdepth 1 -type d -name "*_S*")"
for oldname in $bad_dirs; do 
newname="$(echo "$oldname" | sed -e 's|\(_S[[:digit:]]*\)$||g')"; 
mkdir -p "$newname"
mv ${oldname}/* "${newname}/"
rm -rf "${oldname}"
done


```

