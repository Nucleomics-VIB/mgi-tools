# initialize results
results="mgi_classification.txt"
cat /dev/null > ${results}

# loop through 4 lanes
for lane in L*_DemulX; do

# loop through read2 fastq
for fq in $(find ${lane} -name "V*_2.fq.gz" | sort); do

# extract info
fname=$(basename ${fq})
#echo ${fname}

pfx=${fname%_2.fq.gz}
#echo ${pfx}

label=$(echo ${pfx} | awk 'BEGIN{FS="_";OFS=""}{print $3}')
# echo ${label}

bioawk -c fastx -v class="${label}" 'BEGIN{OFS="\t"}{print $name, class}' ${fq} >> ${results}
done

done
