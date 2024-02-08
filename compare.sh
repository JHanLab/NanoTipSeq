#!/bin/bash

# Check if correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <file1> <file2>"
    exit 1
fi

file1="$1"
file2="$2"

# Extract filenames without path
filename1=$(basename "$file1")
filename2=$(basename "$file2")

# Check if both files exist
if [ ! -f "$file1" ]; then
    echo "Error: $filename1 not found!"
    exit 1
fi

if [ ! -f "$file2" ]; then
    echo "Error: $filename2 not found!"
    exit 1
fi

# Perform the diff
echo "Detailed differences between $filename1 and $filename2:"
diff_output=$(diff -u "$file1" "$file2")

# Output the differences to a text file
report_file="diff_report.txt"

# Prepend each line with the corresponding file name
echo "$diff_output" | awk -v f1="$filename1" -v f2="$filename2" '{
    if ($0 ~ /^---/) {
        print "\033[1m" f1 ": (source)\033[0m"
    } else if ($0 ~ /^\+\+\+/) {
        print "\033[1m" f2 ": (source)\033[0m"
    } else if ($0 ~ /^[-+]/) {
        if ($0 ~ /^-/) {
            print "\033[31m" f1 ": " $0 "\033[0m"
        } else if ($0 ~ /^\+/) {
            print "\033[32m" f2 ": " $0 "\033[0m"
        }
    }
}' > "$report_file"

echo "Detailed differences written to $report_file"
