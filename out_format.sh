#parse network validation results files
#network validation prints to screen, > to file when run
#parse here to get single comma delimited row output/file, append all together
for file in "$@"
do
var=$(cat $file | grep -Ev \{*\} | grep \: | cut -d \: -f 2 | sed 's/ //g' | awk 'BEGIN { ORS = "," } { print }')
echo $file,$var	>> all_output.csv
done
