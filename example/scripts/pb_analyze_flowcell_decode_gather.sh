#!/bin/bash

# Specialized script to gather the results of scattering the phase1 read decoding.

decodeScatter="$1"
outDir="$2"
scatterDir="${outDir}/scatter"

for tag in bounds detail1 detail2 stats; do
    scatterFiles="$(ls ${scatterDir}/*.*@*.${tag}.txt | sort -V)"
    if [ ! -z "${scatterFiles}" ]; then
        outFile=${outDir}/"$(echo ${scatterFiles} | awk -v N=${decodeScatter} '{ n = split($1,p,"/"); f = p[n]; sub("[.]1@"N"[.]",".",f); print f }')"
        if [ "${tag}" == "stats" ]; then
            # The stats file is summed, not concatenated
            cat ${scatterFiles} | body grep -v SAMPLE \
                | awk -v FS="\t" -v OFS="\t" \
                  '{ if (NR == 1) { print; next; }
                     else { S[1] = $1; for (i = 2; i <= NF; i++) { S[i] += $i } } }
                   END { for (i = 1; i <= length(S); i++) { $i = S[i] }; print }' \
                > ${outFile} || exit 1
        else 
            cat ${scatterFiles} | body grep -v ^SAMPLE > ${outFile} || exit 1
        fi
    fi
done
