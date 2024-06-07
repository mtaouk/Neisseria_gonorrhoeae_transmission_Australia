
#!/bin/bash

main_file="dates.txt"

for subset_file in $(cat clustersover5.txt); do
  output_file="dates_${subset_file}"

  mkdir "/home/taouk/NGtransmission/Snippy_subset_2/LSD/$subset_file"
  
  awk -F'\t' 'FNR==NR{a[$1]; next} $1 in a{print $1,$2}' "/home/taouk/NGtransmission/Snippy_subset_2/inputs/$subset_file.filtered.txt" "$main_file" | awk -v count=$(wc -l < "/home/taouk/NGtransmission/Snippy_subset_2/inputs/$subset_file.filtered.txt") 'BEGIN{print count} {print}' > "/home/taouk/NGtransmission/Snippy_subset_2/LSD/$subset_file/$output_file"
done 