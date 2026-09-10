#!/bin/bash

declare -A SAMPLE_ROLE
SAMPLE_ROLE["HG00406"]="child"
SAMPLE_ROLE["HG00407"]="father"
SAMPLE_ROLE["HG00408"]="mother"

REF="/home/BCG2026_exam/chr20"
FA="/home/BCG2026_exam/chr20.fa"
BED="/home/BCG2026_exam/chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed"
VEP_CACHE="/data/vep_cache"

for trio in trio_{1..5};  do
  echo "=== Processing $trio ==="
  TRIO_DIR="$trio"

  for sample in HG00406 HG00407 HG00408; do
    role=${SAMPLE_ROLE[$sample]}
    bam="${TRIO_DIR}/${trio}_${role}.bam"

    echo "--- Aligning $role ($sample) ---"
    bowtie2 -1 ${TRIO_DIR}/${sample}.targets_R1.fq.gz -2 ${TRIO_DIR}/${sample}.targets_R2.fq.gz -x ${REF} --rg-id "${role}" --rg "SM:${role}" | samtools view -Sb | samtools sort -o ${bam}

    echo "--- Indexing $role ---"
    samtools index ${bam}

    echo "--- FastQC on $role ---"
    fastqc ${bam}

    echo "--- Qualimap on $role ---"
    qualimap bamqc -bam ${bam} --feature-file ${BED} -outdir ${TRIO_DIR}/${trio}_${role}
  done

  echo "--- Variant calling ---"
  freebayes -f ${FA} -m 20 -C 5 -Q 10 --min-coverage 10 ${TRIO_DIR}/${trio}_child.bam ${TRIO_DIR}/${trio}_father.bam ${TRIO_DIR}/${trio}_mother.bam > ${TRIO_DIR}/${trio}.vcf

  bgzip ${TRIO_DIR}/${trio}.vcf
  bcftools index ${TRIO_DIR}/${trio}.vcf.gz

  if [[ $trio == "trio_1" || $trio == "trio_4" || $trio == "trio_5" ]]; then
    FILTER='GT[0]="AA" && GT[1]="RA" && GT[2]="RA"'
  elif [[ $trio == "trio_2" ]]; then
    FILTER='GT[0]="RA" && GT[1]="RR" && GT[2]="RA"'
  elif [[ $trio == "trio_3" ]]; then
    FILTER='GT[0]="RA" && GT[1]="RR" && GT[2]="RR"'
  fi

  echo "--- Filtering variants ---"
  bcftools view -R ${BED} ${TRIO_DIR}/${trio}.vcf.gz | bcftools view -S samples.txt | bcftools view -i "${FILTER}" | bcftools filter -i 'QUAL>20' -Ov -o ${TRIO_DIR}/${trio}.cand.vcf 

  echo "--- Annotating with VEP ---"
  vep -i ${TRIO_DIR}/${trio}.cand.vcf -o ${TRIO_DIR}/${trio}.vep_annotated.vcf --vcf --cache --offline --assembly GRCh38 --dir_cache ${VEP_CACHE} --use_given_ref --mane --pick_allele --af --af_1kg --af_gnomade --max_af --sift b --polyphen b --no_fasta

  echo "--- Filtering VEP output ---"
  filter_vep -i ${TRIO_DIR}/${trio}.vep_annotated.vcf -o ${TRIO_DIR}/${trio}.vep_filtered.vcf --filter "(IMPACT is HIGH or IMPACT is MODERATE) and (not MAX_AF or MAX_AF < 0.0001)"

  echo "--- Computing coverage tracks ---"
  for sample in HG00406 HG00407 HG00408; do
    role=${SAMPLE_ROLE[$sample]}
    bam="${TRIO_DIR}/${trio}_${role}.bam"
    bedtools genomecov -ibam ${bam} -bg -trackline -trackopts "name=\"${role}\"" -max 100 > ${TRIO_DIR}/${trio}_${role}Cov.bg
  done

done

echo "--- Running MultiQC ---"
multiqc trio_*/
echo "====== DONE! ======="
