Some bash commands that help with pipeline management

Find pipeline steps that don't have a 'results' directory

`analysis_dir/pipeline$ find . -maxdepth 1 -type d ! -exec test -e "{}/results" \; -print`

Find and delete all 'errors' directories (be careful with this!)

`analysis_dir/pipeline$ find . -mindepth 2  -type d -name "errors" -exec rm -rf {} \;`

Find and delete all 'results' directories except for the alignments (be careful with this!)

`analysis_dir/pipeline$ find . -type d -name "results" ! -path "*align*" -exec rm -rf {} \;`

Find files from specific branch of a results directory and run a custom script on them

`find "$path_to_results_dir" -path "*/peaks.by_sample.macs_broad/*" -name "venn.txt" -print0 | xargs -0 $custom_script`


Move all the pipeline results to a new directory called 'results_old'

```
# switch to the 'pipeline' dir
cd ~/projects/ChIPSeq_project_2016-06-06/pipeline

# get a list of results directories
FILES="$(find . -maxdepth 2 -mindepth 2 -type d -name "results" -exec readlink -f {} \; | tr '\n' ' ')"

# iterate over the directories
for i in $FILES; do
# the old results dir
tmp_dir="$i"
echo "$tmp_dir"

# get the basename of the dir
tmp_basename_old="$(basename $tmp_dir)"
echo "$tmp_basename_old"

# append the timestamp
tmp_basename_new="${tmp_basename_old}_$(date -u +%Y%m%dt%H%M)"
echo "$tmp_basename_new"

# set the path for the results_old dir
old_dir="${tmp_dir}_old"
echo "$old_dir"

# make the results_old dir
mkdir -p "$old_dir"

# move the results dir to the results_old dir
mv "$tmp_dir" "${old_dir}/${tmp_basename_new}"

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

