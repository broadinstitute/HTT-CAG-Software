#!/bin/bash

rcThreshold=10
mdThreshold=4
outRightThreshold=0.20

# Scattered decoding is not implemented in this example analysis
#decodeScatter=10
decodeScatter=0

flowcell=m64298e_220829_033656
#bamPath=data/${flowcell}.ds3.bam
bamPath=data/${flowcell}.bam
outDir=out

workspace=".."
scriptDir="${workspace}/dist/R"
export SV_DIR="../svtoolkit"
export SV_CLASSPATH="${workspace}/dist/HTTCAGSoftware.jar:${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar"

# To run this example, you need to download and supply here a suitable hg38 reference fasta file.
# Ideally, this should be the GRCh38 "no alt" fasta file which was used to align the pacbio data.
# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
# In practice, any hg38 reference with the same chr1 sequence will work, including the standard GRCh38 reference fasta.
# Note that although the example data used here is Q111 mouse, this mouse model has a human HTT exon 1 and  mouse data was aligned to a human hg38 reference.
# The HTT analysis code currently assumes a GRCh38 compatible reference, with respect to the sequence at HTT.
referenceFile=/humgen/cnp04/svtoolkit/references/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta

# Utility function
body() {
    IFS= read -r header
    printf '%s\n' "$header"
    "$@"
}

mkdir -p ${outDir} || exit 1

echo $(date) "Analyzing flowcell ${flowcell} bam file ${bamPath} ..."

prefix="$(basename ${bamPath} | sed 's/.bam$//')"
if [ ! -e ${outDir}/${prefix}.stats.txt ]; then
    if [ -z "${decodeScatter}" -o "${decodeScatter}" == 0 ]; then
        echo $(date) "Running PBSCHTTAnalyzer ..."
        java -Xmx2g -cp ${SV_CLASSPATH} org.broadinstitute.sv.custom.htt.PBSCHTTAnalyzer \
            -R ${referenceFile} \
            -I ${bamPath} \
            -statsOutputFile ${outDir}/${prefix}.stats.txt \
            -detailOutputFile1 ${outDir}/${prefix}.detail1.txt \
            -detailOutputFile2 ${outDir}/${prefix}.detail2.txt \
            -boundaryOutputFile ${outDir}/${prefix}.bounds.txt \
            || exit 1
    else
        echo $(date) "Scattering PBSCHTTAnalyzer ..."
        scatterDir=${outDir}/scatter
        scatterArgs="$(seq 1 ${decodeScatter} | awk -v N=${decodeScatter} '{ print $1 "@" N }')"
        scatterArgs="$(echo ${scatterArgs} | sed 's/ /,/g')"
        rm -rf ${scatterDir}/logs
        mkdir -p ${scatterDir}/logs || exit 1
        /stanley/genome_strip/tools/scripts/run_queue_scatter.sh \
            -jobLogDir ${scatterDir}/logs \
            scripts/scatter_by_substitution.sh \
                ${scatterArgs} \
                "java -Xmx2g -cp ${SV_CLASSPATH} org.broadinstitute.sv.custom.htt.PBSCHTTAnalyzer \
                     -R ${referenceFile} \
                     -I ${bamPath} \
                     -scatter @@ \
                     -statsOutputFile ${scatterDir}/${prefix}.@@.stats.txt \
                     -detailOutputFile1 ${scatterDir}/${prefix}.@@.detail1.txt \
                     -detailOutputFile2 ${scatterDir}/${prefix}.@@.detail2.txt \
                     -boundaryOutputFile ${scatterDir}/${prefix}.@@.bounds.txt" \
            || exit 1
        echo $(date) "Gathering PBSCHTTAnalyzer ..."
        scripts/pb_analyze_flowcell_decode_gather.sh ${decodeScatter} ${outDir} || exit 1
        rm -rf ${scatterDir}
    fi
fi

if [ ! -e ${outDir}/${prefix}.out1.txt ]; then
    echo $(date) "Running PBSCHTTAnalyzer2 ..."
    # Currently we do not use the CBC or STAMP files in this step
    # Copy the index file to the output directory to match our standard directory layout
    indexFile=${outDir}/${prefix}.indexes.txt
    cp data/${flowcell}.indexes.txt ${indexFile} || exit 1
    java -Xmx2g -cp ${SV_CLASSPATH} org.broadinstitute.sv.custom.htt.PBSCHTTAnalyzer2 \
        -R ${referenceFile} \
        -indexFile ${indexFile} \
        -detailFile1 ${outDir}/${prefix}.detail1.txt \
        -detailFile2 ${outDir}/${prefix}.detail2.txt \
        -outputFile1 ${outDir}/${prefix}.out1.txt \
        -outputFile2 ${outDir}/${prefix}.out2.txt \
        || exit 1
fi

if [ ! -e ${outDir}/${prefix}.readinfo.txt ]; then
    # Filter out reads with non-unique alignments
    echo $(date) "Merging to create readinfo.txt file ..."
    cat ${outDir}/${prefix}.out2.txt | cut -f 1 | body sort | uniq -c | awk '$1 == 1 { print $2 }' > ${outDir}/${prefix}.unique.list || exit 1
    for tag in out1 out2; do
        cat ${outDir}/${prefix}.${tag}.txt | body sort | join -t $'\t' ${outDir}/${prefix}.unique.list - > ${outDir}/${prefix}.${tag}.tmp || exit 1
    done
    cut -f 1,2,3,8- ${outDir}/${prefix}.out1.tmp | join -t $'\t' - ${outDir}/${prefix}.out2.tmp > ${outDir}/${prefix}.readinfo.txt || exit 1
    rm -f ${outDir}/${prefix}.out1.tmp ${outDir}/${prefix}.out2.tmp
fi

stampMapFile=${outDir}/${prefix}.stampfiles.txt
if [ ! -e ${stampMapFile} ]; then
    # Copy the stamp map file to the output directory to match our standard directory layout
    cp data/${flowcell}.stampfiles.txt ${stampMapFile} || exit 1
fi

if [ ! -e ${outDir}/qc/${prefix}.t=${rcThreshold}.umi_summary.txt ]; then
    echo $(date) "Running analyze_umis t=${rcThreshold} for QC (no stamps) ..."
    mkdir -p ${outDir}/qc || exit 1
    Rscript ${scriptDir}/analyze_umis.R \
        --minReadPurity 0.90 \
        --minReadCount ${rcThreshold} \
        --stampMinReadCount ${rcThreshold} \
        --outputReadInfo ${outDir}/qc/${prefix}.t=${rcThreshold}.readinfo.qc.txt \
        ${outDir}/${prefix}.readinfo.txt \
        ${outDir}/qc/${prefix}.t=${rcThreshold}.umi_summary.txt \
        |& tee ${outDir}/qc/${prefix}.t=${rcThreshold}.analyze_umis.log \
        || exit 1
fi

if [ ! -e ${outDir}/filtered/${prefix}.t=${rcThreshold}.umi_summary.txt -a -e ${stampMapFile} ]; then
    echo $(date) "Running analyze_umis t=${rcThreshold} d=${mdThreshold} ..."
    mkdir -p ${outDir}/filtered || exit 1
    Rscript ${scriptDir}/analyze_umis.R \
        --minReadPurity 0.90 \
        --maxUMILength 28 \
        --minReadCount ${rcThreshold} \
        --minDistanceThreshold ${mdThreshold} \
        --outRightThreshold ${outRightThreshold} \
        --stampMapFile ${stampMapFile} \
        --stampMinReadCount ${rcThreshold} \
        --outputReadInfo ${outDir}/filtered/${prefix}.t=${rcThreshold}.readinfo.qc.txt \
        ${outDir}/${prefix}.readinfo.txt \
        ${outDir}/filtered/${prefix}.t=${rcThreshold}.umi_summary.txt \
        |& tee ${outDir}/filtered/${prefix}.t=${rcThreshold}.analyze_umis.log \
        || exit 1
fi

if [ ! -e ${outDir}/filtered2/${prefix}.t=${rcThreshold}.umi_summary.txt -a -e ${stampMapFile} ]; then
    echo $(date) "Running analyze_umis t=${rcThreshold} d=${mdThreshold} --assignReactionsFromStamps true ..."
    mkdir -p ${outDir}/filtered2 || exit 1
    Rscript ${scriptDir}/analyze_umis.R \
        --assignReactionsFromStamps true \
        --minReadPurity 0.90 \
        --maxUMILength 28 \
        --minReadCount ${rcThreshold} \
        --minDistanceThreshold ${mdThreshold} \
        --outRightThreshold ${outRightThreshold} \
        --stampMapFile ${stampMapFile} \
        --stampMinReadCount ${rcThreshold} \
        --outputReadInfo ${outDir}/filtered2/${prefix}.t=${rcThreshold}.readinfo.qc.txt \
        ${outDir}/${prefix}.readinfo.txt \
        ${outDir}/filtered2/${prefix}.t=${rcThreshold}.umi_summary.txt \
        |& tee ${outDir}/filtered2/${prefix}.t=${rcThreshold}.analyze_umis.log \
        || exit 1
fi

echo $(date) "Analysis done."
