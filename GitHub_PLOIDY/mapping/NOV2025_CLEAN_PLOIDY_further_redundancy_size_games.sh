###########################################################################################################
###########################################################################################################
######################### Control species 
###########################################################################################################
###########################################################################################################


# download all SRA files (.sra)
cat accessions.txt | xargs prefetch

# convert to FASTQ (gzipped) quickly
for acc in $(cat accessions.txt); do fasterq-dump ./$acc/*.sra --split-files ; done

# renaming


## Tried all (RF, FR, UNSTR) - because according to https://docs.google.com/presentation/d/1glkGK1dI8bkFvmVxTtvFp43CkpzgH926aMGrUisRpF0/edit?slide=id.g2e387da3bdc_0_8#slide=id.g2e387da3bdc_0_8 
# strandedness setting maters a lot
 
for i in $(cat vocabulary_ploid_controls.txt | cut -f 2) ; do mv $i*_1.fastq "$(grep $i vocabulary_ploid_controls.txt | cut -f 1)"_1.fastq ; done
for i in $(cat vocabulary_ploid_controls.txt | cut -f 2) ; do mv $i*_2.fastq "$(grep $i vocabulary_ploid_controls.txt | cut -f 1)"_2.fastq ; done
for i in $(cat vocabulary_ploid_controls.txt | cut -f 2) ; do mv $i.fastq "$(grep $i vocabulary_ploid_controls.txt | cut -f 1)".fastq ; done
 
## paired-end 
 
# RF 
for i in $(cat vocabulary_ploid_controls.txt | cut -f 1) ; do rnaspades -1 "$i"_1.fastq -2 "$i"_2.fastq --ss-rf -o ../control_ass/RF_STR_RNAspades_control/"$i"/ ; done

# FR 
for i in $(cat vocabulary_ploid_controls.txt | cut -f 1) ; do rnaspades -1 "$i"_1.fastq -2 "$i"_2.fastq --ss-fr -o ../control_ass/FR_STR_RNAspades_control/"$i"/ ; done

#UNSTR 
for i in $(cat vocabulary_ploid_controls.txt | cut -f 1) ; do rnaspades -1 "$i"_1.fastq -2 "$i"_2.fastq -o ../control_ass/UNSTR_RNAspades_control/"$i"/ ; done


## single-end 
 
# RF 
for i in $(cat vocabulary_ploid_controls.txt | cut -f 1) ; do rnaspades -s "$i".fastq --ss-rf -o ../control_ass/RF_STR_RNAspades_control/"$i"/ ; done

# FR 
for i in $(cat vocabulary_ploid_controls.txt | cut -f 1) ; do rnaspades -s "$i".fastq --ss-fr -o ../control_ass/FR_STR_RNAspades_control/"$i"/ ; done

#UNSTR 
for i in $(cat vocabulary_ploid_controls.txt | cut -f 1) ; do rnaspades -s "$i".fastq -o ../control_ass/UNSTR_RNAspades_control/"$i"/ ; done


##continue assemblies that gor interrupted due to lack of RAM
cd /scratch/data2/dzavadska/QUANTITATIVE_DATASET/PLOIDY/control_ass/RF_STR_RNAspades_control/

for i in $(echo "Cornospumella_fuschlensis_DIPL" "Leishmania_infantum_DIPL" "Emiliania_huxleyi_DIPL") ; do rnaspades --continue -o "$i" ; done
rnaspades -1 ../../control_reads/Poteriospumella_lacustris_TRIPL_1.fastq -2 ../../control_reads/Poteriospumella_lacustris_TRIPL_2.fastq --ss-rf -o Poteriospumella_lacustris_TRIPL

cd /scratch/data2/dzavadska/QUANTITATIVE_DATASET/PLOIDY/control_ass/FR_STR_RNAspades_control
rnaspades -1 ../../control_reads/Poteriospumella_lacustris_TRIPL_1.fastq -2 ../../control_reads/Poteriospumella_lacustris_TRIPL_2.fastq --ss-fr -o Poteriospumella_lacustris_TRIPL



#copying and renaming all

for i in ./* ; do cp $i/transcripts.fasta ../"$i"_FR_contigs.fasta ; done
for i in ./* ; do cp $i/transcripts.fasta ../"$i"_RF_contigs.fasta ; done
for i in ./* ; do cp $i/transcripts.fasta ../../"$i"_UNSTR_contigs.fasta ; done
# ./Cornospumella_fuschlensis_DIPL RF runs out of memory

#!!!!!!!!!!!!!!!!!!!!!!!!!
## Decontamination skipped 
#!!!!!!!!!!!!!!!!!!!!!!!!!

cd /scratch/data2/dzavadska/QUANTITATIVE_DATASET/PLOIDY/transdecoder_BUSCO_control/
 cp ../control_ass/*FR_contigs* ./
  cp ../control_ass/*RF_contigs* ./
    cp ../control_ass/*UNSTR_contigs* ./
 
for file in ./*FR_contigs*.fasta ; do sudo TransDecoder.LongOrfs -t $file ; done
for file in ./*FR_contigs*.fasta ; do sudo TransDecoder.Predict -t $file ; done

for file in ./*RF_contigs*.fasta ; do TransDecoder.LongOrfs -t $file ; done
for file in ./*RF_contigs*.fasta ; do TransDecoder.Predict -t $file ; done

for file in ./*UNSTR_contigs*.fasta ; do TransDecoder.LongOrfs -t $file ; done
for file in ./*UNSTR_contigs*.fasta ; do TransDecoder.Predict -t $file ; done

for file in *FR_contigs*.pep ; do sudo busco -i $file -l eukaryota_odb10 -o "$file"_FR_contigs -m proteins --cpu 8 ; done

for file in *RF_contigs*.pep ; do sudo busco -i $file -l eukaryota_odb10 -o "$file"_RF_contigs -m proteins --cpu 8 ; done

for file in *UNSTR_contigs*.pep ; do sudo busco -i $file -l eukaryota_odb10 -o "$file"_UNSTR_contigs -m proteins --cpu 8 ; done


## manually overview BUSCO values, and copy the assenblies with highhest BUSCO to the next folder
cp *RF_contigs*.pep ../control_cdhit_mapping/
cd /scratch/data2/dzavadska/QUANTITATIVE_DATASET/PLOIDY/control_cdhit_mapping
rm Leishmania_infantum_DIPL* Leishmania_major_DIPL* Monocercomonoides_exilis_HAPL*

cp ../transdecoder_BUSCO_control/Leishmania_infantum_DIPL*UNSTR_contigs*.pep ./
cp ../transdecoder_BUSCO_control/Leishmania_major_DIPL*UNSTR_contigs*.pep ./
cp ../transdecoder_BUSCO_control/Monocercomonoides_exilis_HAPL*UNSTR_contigs*.pep ./

cp ../transdecoder_BUSCO_control/Cornospumella_fuschlensis_DIPL*FR_contigs*.pep ./




cp *RF_contigs*.cds ../control_cdhit_mapping/
cd /scratch/data2/dzavadska/QUANTITATIVE_DATASET/PLOIDY/control_cdhit_mapping
rm Leishmania_infantum_DIPL*RF_contigs* Leishmania_major_DIPL*RF_contigs* Monocercomonoides_exilis_HAPL*RF_contigs*

cp ../transdecoder_BUSCO_control/Leishmania_infantum_DIPL*UNSTR_contigs*.cds ./
cp ../transdecoder_BUSCO_control/Leishmania_major_DIPL*UNSTR_contigs*.cds ./
cp ../transdecoder_BUSCO_control/Monocercomonoides_exilis_HAPL*UNSTR_contigs*.cds ./

cp ../transdecoder_BUSCO_control/Cornospumella_fuschlensis_DIPL*FR_contigs*.cds ./

cp ../control_reads/vocabulary_ploid_controls.txt  ./



# cd-hit
cd /scratch/data2/dzavadska/QUANTITATIVE_DATASET/PLOIDY/control_cdhit_mapping

for i in ./*.pep ; do cd-hit -i "$i" -o ./"$i"_NR90.txt -c 0.90 -n 5 -T 4 ; done


######################################################
####################### mapping our species #######################
######################################################


cd /scratch/data2/dzavadska/QUANTITATIVE_DATASET/PLOIDY/mapping

#copying CDS from TransDecoder 

for i in $(echo "RhMon" "G65133" "D44" "16Ckin" "13G" "10D" "Gemkin") ; do cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/ALL_SPECIES_1STSELECTION_CLEAN_DATA/CDS/"$i"*.cds ./"$i"_ORIGINAL_CDS.fasta ; done
for i in $(echo "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/NICELY_CLEANED_ASS_batch2/"$i"*.cds ./"$i"_ORIGINAL_CDS.fasta ; done


!!!!!!!!!!!!!!!!!!!!!!
#remove incorrect-stranded assembly of GemKin
rm Gemkin_doubleclean_contigs.fasta.transdecoder.cds


#####redundancy thing
#throwing away predicted proteins that were removed during redundancy reduction by aa pid =90 from initial CDS files 

#for i in $(echo "RhMon" "G65133" "D44" "16Ckin" "13G" "10D" "Gemkin") ; do grep ">" /scratch/data2/dzavadska/QUANTITATIVE_DATASET/NR90_to_annotate/"$i".fasta >> "$i"_NR90_tokeep.txt ; done

#for i in $(echo "RhMon" "G65133" "D44" "16Ckin" "13G" "10D" "Gemkin") ; do seqkit grep -n -f ./"$i"_NR90_tokeep.txt ./"$i"*.cds -o "$i"_NR90_CDS.fasta ; done

#####redundancy thing


#################################################################################
###################### MAPPING KINETOS ##########################################
#################################################################################
for i in $(echo "RhMon" "G65133" "D44" "16Ckin" "13G" "10D" "Gemkin") ; do bowtie2-build ./"$i"_ORIGINAL_CDS.fasta "$i"_ORIGINAL_CDS_ref ; done
for i in $(echo "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do bowtie2-build ./"$i"_ORIGINAL_CDS.fasta "$i"_ORIGINAL_CDS_ref ; done


for i in $(echo "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do
bowtie2 -p 32 -X 2000 -x "$i"_ORIGINAL_CDS_ref -1 /scratch/data2/dzavadska/QUANTITATIVE_DATASET/reads/"$i"trimmed_1.fastq.gz -2 /scratch/data2/dzavadska/QUANTITATIVE_DATASET/reads/"$i"trimmed_2.fastq.gz | samtools sort -m 10G -@ 10 -T "$i"_tmp1 -o "$i"_ORIGINAL_CDS_ref.bam ; done

#specially for Gemkin
bowtie2 -p 32 -X 2000 -x Gemkin_ORIGINAL_CDS_ref -1 /scratch/data2/dzavadska/QUANTITATIVE_DATASET/reads/GK1_R1.fastq.gz -2 /scratch/data2/dzavadska/QUANTITATIVE_DATASET/reads/GK1_R2.fastq.gz | samtools sort -m 10G -@ 10 -T "$i"_tmp -o Gemkin_ORIGINAL_CDS_ref.bam


for i in $(echo "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do
bowtie2 -p 32 -X 2000 -x "$i"_ORIGINAL_CDS_ref -1 /scratch/data2/dzavadska/COMPARATIVE_DATASET/reads_and_co/"$i"_R1_trimmed.fastq.gz -2 /scratch/data2/dzavadska/COMPARATIVE_DATASET/reads_and_co/"$i"_R2_trimmed.fastq.gz | samtools sort -m 10G -@ 10 -T "$i"_tmp1 -o "$i"_ORIGINAL_CDS_ref.bam ; done


#indexing .bam files (for visualising in IGV)
for i in $(echo "RhMon" "G65133" "D44" "16Ckin" "13G" "10D" "Gemkin") ; do samtools index "$i"_ORIGINAL_CDS_ref.bam ; done
for i in $(echo "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do samtools index "$i"_ORIGINAL_CDS_ref.bam ; done


##############################################
# subsetting reads (not exactly to the uniquely mapped reads, but by mapping quality)

bowtie2 -x index_prefix -1 reads_1.fastq -2 reads_2.fastq -S output.sam --very-sensitive -k 1 --no-unal

for i in $(echo "RhMon" "G65133" "D44" "16Ckin" "13G" "10D") ; do echo "$i" ; samtools stats "$i"_ORIGINAL_CDS_ref.bam | grep "average quality" ; done
for i in $(echo "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do echo "$i" ; samtools stats "$i"_ORIGINAL_CDS_ref.bam | grep "average quality" ; done
#RhMon
#SN	average quality:	36.0
#G65133
#SN	average quality:	36.1
#D44
#SN	average quality:	36.0
#16Ckin
#SN	average quality:	36.0
#13G
#SN	average quality:	35.9
#10D
#SN	average quality:	35.9
#BABKIN
#SN	average quality:	39.3
#INUKIN
#SN	average quality:	39.3
#KUCKIN
#SN	average quality:	39.3
#OB23
#SN	average quality:	35.1


for i in $(echo "RhMon" "G65133" "D44" "16Ckin" "13G" "10D" "Gemkin") ; do samtools view -b -q 25 "$i"_ORIGINAL_CDS_ref.bam > "$i"_ORIGINAL_CDS_Q25.bam ; done

for i in $(echo "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do samtools view -b -q 25 "$i"_ORIGINAL_CDS_ref.bam > "$i"_ORIGINAL_CDS_Q25.bam ; done


#estimating average coverage for each gene

for i in $(echo "RhMon" "G65133" "D44" "16Ckin" "13G" "10D" "Gemkin") ; do samtools depth -a "$i"_ORIGINAL_CDS_Q25.bam | awk '{cov[$1]+=$3; len[$1]++} END {for (t in cov) print t, cov[t]/len[t]}' >  "$i"_ORIGINAL_CDS_cov_perCDS_Q25.txt ; done


for i in $(echo "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do samtools depth -a "$i"_ORIGINAL_CDS_Q25.bam | awk '{cov[$1]+=$3; len[$1]++} END {for (t in cov) print t, cov[t]/len[t]}' >  "$i"_ORIGINAL_CDS_cov_perCDS_Q25.txt ; done



#### SNP calling itself 
#### reads filtered by Q25!!!

for i in $(echo"D44" "16Ckin" "13G" "10D" "Gemkin") ; do bcftools mpileup -a AD,DP -Ou -f "$i"_ORIGINAL_CDS.fasta "$i"_ORIGINAL_CDS_Q25.bam | bcftools view -V indels -Ov -o "$i"_ORIGINAL_CDS_Q25_BCFtools_view.vcf ; done
for i in $(echo "RhMon" "G65133" ) ; do bcftools mpileup -A -a AD,DP -Ou -f "$i"_ORIGINAL_CDS.fasta "$i"_ORIGINAL_CDS_Q25.bam | bcftools view -V indels -Ov -o "$i"_ORIGINAL_CDS_Q25_BCFtools_view.vcf ; done

for i in $(echo  "BABKIN" "INUKIN" "KUCKIN" "OB23" ) ; do bcftools mpileup -A -a AD,DP -Ou -f "$i"_ORIGINAL_CDS.fasta "$i"_ORIGINAL_CDS_Q25.bam | bcftools view -V indels -Ov -o "$i"_ORIGINAL_CDS_Q25_BCFtools_view.vcf ; done



# variants were called while ignoring indels,ﬁltering for a minimum read depth of 5, a minimum of four alternate bases before calling a SNP, accepting only bi-allelicSNPs, and only accepting SNPs with a MAF(=minor allel frequency) of 2%, and maximum missing data of 50%

--min-af 0.02
for i in $(echo "RhMon" "G65133" "D44" "16Ckin" "13G" "10D" "Gemkin") ; do bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\t%DP\t%QUAL\n' "$i"_ORIGINAL_CDS_Q25_BCFtools_view.vcf >  "$i"_ORIGINAL_CDS_Q25_BCFtools_view.tsv ; done
for i in $(echo "RhMon" "G65133" "D44" "16Ckin" "13G" "10D" "Gemkin") ; do sed -i 's/,/\t/g' "$i"_ORIGINAL_CDS_Q25_BCFtools_view.tsv ; done
for i in $(echo "RhMon" "G65133" "D44" "16Ckin" "13G" "10D" "Gemkin") ; do awk  -F'\t' 'NF == 10' "$i"_ORIGINAL_CDS_Q25_BCFtools_view.tsv > "$i"_ORIGINAL_CDS_Q25_BCFtools_viewfiltered.tsv ; done

for i in $(echo  "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\t%DP\t%QUAL\n' "$i"_ORIGINAL_CDS_Q25_BCFtools_view.vcf >  "$i"_ORIGINAL_CDS_Q25_BCFtools_view.tsv ; done
for i in $(echo  "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do sed -i 's/,/\t/g' "$i"_ORIGINAL_CDS_Q25_BCFtools_view.tsv ; done
for i in $(echo  "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do awk  -F'\t' 'NF == 10' "$i"_ORIGINAL_CDS_Q25_BCFtools_view.tsv > "$i"_ORIGINAL_CDS_Q25_BCFtools_viewfiltered.tsv ; done





#estimating the length of every transcript
for i in $(echo "RhMon" "G65133" "D44" "16Ckin" "13G" "10D" "Gemkin") ; do grep ">" "$i"_ORIGINAL_CDS.fasta > "$i"_ORIGINAL_CDS_length.txt ; done
#estimating the length of every transcript
for i in $(echo "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do grep ">" "$i"_ORIGINAL_CDS.fasta > "$i"_ORIGINAL_CDS_length.txt ; done



for i in $(echo  "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do scp nisaba:/scratch/data2/dzavadska/QUANTITATIVE_DATASET/PLOIDY/mapping/"$i"*viewfiltered.tsv ./ ; done
for i in $(echo  "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do scp nisaba:/scratch/data2/dzavadska/QUANTITATIVE_DATASET/PLOIDY/mapping/"$i"*cov_perCDS_Q25.txt ./ ; done
for i in $(echo  "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do scp nisaba:/scratch/data2/dzavadska/QUANTITATIVE_DATASET/PLOIDY/mapping/"$i"*Q25_BCFtools_viewfiltered.tsv ./ ; done
for i in $(echo  "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do scp nisaba:/scratch/data2/dzavadska/QUANTITATIVE_DATASET/PLOIDY/mapping/"$i"*_ORIGINAL_CDS_length.txt ./ ; done







######For PCRs

for i in $(echo  "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do scp nisaba:/scratch/data2/dzavadska/COMPARATIVE_DATASET/NICELY_CLEANED_ASS_batch2/"$i"*_doubleclean_contigs.fasta.transdecoder.pepEUGLENOZOApep/run_euglenozoa_odb10/full_table.tsv ./"$i"_BUSCO_EUGLoriginal.tsv ; done

for i in $(echo  "BABKIN" "INUKIN" "KUCKIN" "OB23") ; do scp nisaba:/scratch/data2/dzavadska/COMPARATIVE_DATASET/NICELY_CLEANED_ASS_batch2/"$i"*_NR90_EUGLENOZOApep/run_euglenozoa_odb10/full_table.tsv ./"$i"_BUSCO_EUGL.tsv ; done





