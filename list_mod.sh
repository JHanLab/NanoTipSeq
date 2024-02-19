#!/bin/bash

input_file="$1"  # Replace with your input file
output_file="$2"          # Replace with your desired output file
# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "Input file '$input_file' not found."
  exit 1
fi

# Initialize color counter
color_counter=127  # Decimal equivalent of ff7f00 in base 10

# Create the header line
header="Type\tShape\tChr\tStart\tEnd\tcolor"

# Initialize the output file with the header
echo -e "$header" > "$output_file"

# Loop through each line in the input file
while IFS=$'\t' read -r line; do
  # Split the line into separate columns
  read -a columns <<< "$line"

  # Extract values from the columns
  chr="${columns[0]}"
  start="${columns[1]}"
  end="${columns[2]}"
  
  # Append "insertion" for the Type column
  type="insertion"
  
  # Append "circle" for the Shape column
  shape="circle"
  
  # Convert color_counter to hexadecimal
  hex_color=$(printf "%02x" "$color_counter")
  color="e86051"
  
  # Increment color_counter
  ((color_counter++))
  
  # Write the formatted line to the output file
  echo -e "$type\t$shape\t$chr\t$start\t$end\t$color" >> "$output_file"
done < "$input_file"

echo "Output written to '$output_file'"
