## download all reads which were not previously used 

cd /scratch/data2/dzavadska/QUANTITATIVE_DATASET/BUSCO_TITRATION_CLEAN_NOV2025/new_assemblies


#manually create accessions.txt from google spreadsheet table

# download all SRA files (.sra)
cat accessions.txt | xargs prefetch

# convert to FASTQ (gzipped) quickly
for acc in $(cat accessions.txt); do fasterq-dump ./$acc/*.sra --split-files ; done

# renaming


## Tried all (RF, FR, UNSTR) - because according to https://docs.google.com/presentation/d/1glkGK1dI8bkFvmVxTtvFp43CkpzgH926aMGrUisRpF0/edit?slide=id.g2e387da3bdc_0_8#slide=id.g2e387da3bdc_0_8 
# strandedness setting maters a lot
 
 
for i in $(cat vocabulary_busco_controls.txt | cut -f 2) ; do mv $i*_1.fastq "$(grep $i vocabulary_busco_controls.txt | cut -f 1)"_1.fastq ; done
for i in $(cat vocabulary_busco_controls.txt | cut -f 2) ; do mv $i*_2.fastq "$(grep $i vocabulary_busco_controls.txt | cut -f 1)"_2.fastq ; done
for i in $(cat vocabulary_busco_controls.txt | cut -f 2) ; do mv $i.fastq "$(grep $i vocabulary_busco_controls.txt | cut -f 1)".fastq ; done
 
## paired-end 
# RF 
for i in $(cat vocabulary_busco_controls.txt | cut -f 1) ; do rnaspades -1 "$i"_1.fastq -2 "$i"_2.fastq --ss-rf -o ./RF_STR_RNAspades_new_assemblies/"$i"/ ; done

# FR 
for i in $(cat vocabulary_busco_controls.txt | cut -f 1) ; do rnaspades -1 "$i"_1.fastq -2 "$i"_2.fastq --ss-fr -o ./FR_STR_RNAspades_new_assemblies/"$i"/ ; done

#UNSTR 
for i in $(cat vocabulary_busco_controls.txt | cut -f 1) ; do rnaspades -1 "$i"_1.fastq -2 "$i"_2.fastq -o ./UNSTR_RNAspades_new_assemblies/"$i"/ ; done


#copy from each directory the corresponding Spades assembly

for i in * ; do cp ./"$i"/transcripts.fasta ./"$i"_UNSTR_contigs.fasta ; done
for i in * ; do cp ./"$i"/transcripts.fasta ./"$i"_FR_contigs.fasta ; done
for i in * ; do cp ./"$i"/transcripts.fasta ./"$i"_RF_contigs.fasta ; done


for file in ./*_contigs.fasta ; do sudo TransDecoder.LongOrfs -t $file -S ; done
for file in ./*_contigs.fasta ; do sudo TransDecoder.Predict -t $file ; done


for i in *transdecoder.pep ; do sudo busco -i $i -l eukaryota_odb10 -o "$i"pep -m proteins --cpu 8 ; done

# for all assemblies which were completed, those with RF orientation all have highest BUSCO

#cd-hit
# cd /scratch/data2/dzavadska/QUANTITATIVE_DATASET/BUSCO_TITRATION_CLEAN_NOV2025/new_assemblies/RF_STR_RNAspades_new_assemblies

cat ../vocabulary_busco_controls.txt | cut -f 1 | while read i ; do cd-hit -i ./"$i"*transdecoder.pep -o "$i"_NR90_denovo.fasta; done


#### copy to corresponding dir

cp ./*_NR90_denovo.fasta /scratch/data2/dzavadska/QUANTITATIVE_DATASET/BUSCO_TITRATION_CLEAN_NOV2025/ALL_NR90
#rename RF de nov new assemblies and copy them into final destination 
cat ../vocabulary_busco_controls.txt | cut -f 1 | while read i ; do cp ./"$i"*fasta.transdecoder.pep /scratch/data2/dzavadska/QUANTITATIVE_DATASET/BUSCO_TITRATION_CLEAN_NOV2025/ALL_ORIGINAL/"$i"_ORIGINAL_denovo.fasta ; done










###########
# Downloading reference genomes

#cd /scratch/data2/dzavadska/QUANTITATIVE_DATASET/BUSCO_TITRATION_CLEAN_NOV2025/new_genomes

cat new_genomes.txt | while read i ; do datasets download genome accession $i --include rna,cds,protein,genome ; unzip -o ./ncbi_dataset.zip ; cp ./ncbi_dataset/data/GC*/protein.faa ./"$i"_protein.faa ; cp ./ncbi_dataset/data/GC*/*cds_from_genomic.fna ./"$i"_cds_from_genomic.fna ; mv ./ncbi_dataset  "$i" ; done 

for i in $(cat vocabulary_genomes_controls.txt | cut -f 2) ; do mv $i*_protein.faa "$(grep $i vocabulary_genomes_controls.txt | cut -f 1)"_ORIGINAL_REF.fasta ; done


#CD-hit references

for i in $(cat vocabulary_genomes_controls.txt | cut -f 1) ; do cd-hit -i ./"$i"*_ORIGINAL_REF.fasta -o "$i"_NR90_REF.fasta ; done






########
## copying assemblies from phylogenomics dataset
# cd /scratch/data2/dzavadska/QUANTITATIVE_DATASET/BUSCO_TITRATION_CLEAN_NOV2025/ALL_NR90
cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/NR_proteins_90_49taxa/*fasta ./


#copying contigs for transcriptome mode

for i in $(cat vocabulary_genomes_controls.txt | cut -f 2) ; do mv $i/data/GC*/rna.fna "$(grep $i vocabulary_genomes_controls.txt | cut -f 1)"_ORIGINAL_REF_contigs.fasta ; done
for i in $(echo "GCA_000188675.2" "GCA_000146045.2" "GCA_041937235.1" ) ; do mv $i/data/GC*/cds_from_genomic.fna "$(grep $i vocabulary_genomes_controls.txt | cut -f 1)"_ORIGINAL_REF_contigs.fasta ; done


# Trypanosoma_cruzi	GCA_000188675.2
# Saccharomyces_cerevisiae	GCA_000146045.2
# Willaertia_magna	GCA_041937235.1

# cd /scratch/data2/dzavadska/QUANTITATIVE_DATASET/BUSCO_TITRATION_CLEAN_NOV2025/new_assemblies/ 
for i in $(cat vocabulary_busco_controls.txt | cut -f 1) ; do
cp /scratch/data2/dzavadska/QUANTITATIVE_DATASET/BUSCO_TITRATION_CLEAN_NOV2025/new_assemblies/RF_STR_RNAspades_new_assemblies/"$i"*contigs.fasta ../NUCLEOTIDE_ORIGINAL/"$i"_ORIGINAL_denovo_contigs.fasta ; done


/scratch/data2/dzavadska/QUANTITATIVE_DATASET/BUSCO_TITRATION_CLEAN_NOV2025/NUCLEOTIDE_ORIGINAL/


#copy de novo assembled contigs from ploidy 
# cd /scratch/data2/dzavadska/QUANTITATIVE_DATASET/BUSCO_TITRATION_CLEAN_NOV2025/NUCLEOTIDE_ORIGINAL
for i in $(cat de_novo_ploidy.txt ) ; do cp /scratch/data2/dzavadska/QUANTITATIVE_DATASET/PLOIDY/control_cdhit_mapping/"$i"*ORIGINAL_CDS.fasta "$i"_ORIGINAL_denovo_contigs.fasta ; done


#and from batch2 set
for i in $(echo "OB23" "BABKIN" "KUCKIN" "INUKIN" ) ; do cp /scratch/data2/dzavadska/COMPARATIVE_DATASET/NICELY_CLEANED_ASS_batch2/"$i"_doubleclean_contigs.fasta ./"$i"_ORIGINAL_REF_contigs.fasta ; done











#renaming 
cat names_ref_comp.txt | while read i ; do mv $i*.fasta "$i"_NR90_REF_comp.fasta ; done
cat names_NONref_comp.txt | while read i ; do mv $i*.fasta "$i"_NR90_NON_comp.fasta ; done

##coping assemblies from ploidy dataset


cat de_novo_ploidy.txt | while read i ; do cp /scratch/data2/dzavadska/QUANTITATIVE_DATASET/PLOIDY/control_cdhit_mapping/"$i"*transdecoder.pep_NR90.txt ./"$i"_NR90_denovo.fasta ; done

cat de_novo_ploidy.txt | while read i ; do cp /scratch/data2/dzavadska/QUANTITATIVE_DATASET/PLOIDY/control_cdhit_mapping/"$i"*transdecoder.pep ./"$i"_ORIGINAL_denovo.fasta ; done







######
#######RUNNING BUSCO ON TITRATION DATASETS

for i in *fasta ; do sudo busco -i $i -l eukaryota_odb10 -o "$i"_EUK -m proteins --cpu 8 ; done

for i in *fasta ; do sudo busco -i $i -l euglenozoa_odb10 -o "$i"_EUGL -m proteins --cpu 8 ; done

#FOR CONTIGS; run on dazhbog because in nisaba refuses to work
for i in *fasta ; do busco -i $i -f -l eukaryota_odb10 -o "$i"_EUK -m tran --cpu 4 ; done



for i in *_EUK ; do cp "$i"/short_summary*txt ./*_EUK_busco_titr_summary ; done
for i in *_EUGL ; do cp "$i"/short_summary*txt ./*_EUGL_busco_titr_summary ; done



echo -e "file\tC\tS\tD\tF\tM\tn"
for f in *txt ; do
    vals=$(grep -E 'Complete BUSCOs \(C\)|single-copy|duplicated|Fragmented|Missing|Total BUSCO groups' "$f" \
           | awk '{print $1}' | paste -sd' ' -)
    echo -e "$f\t$vals" >> NR90_EUK_busco_titr_summary.tsv
done

echo -e "file\tC\tS\tD\tF\tM\tn"
for f in *txt ; do
    vals=$(grep -E 'Complete BUSCOs \(C\)|single-copy|duplicated|Fragmented|Missing|Total BUSCO groups' "$f" \
           | awk '{print $1}' | paste -sd' ' -)
    echo -e "$f\t$vals"  >> NUCLEOTIDE_EUGL_busco_titr_summary.tsv
done

scp nisaba:/scratch/data2/dzavadska/QUANTITATIVE_DATASET/BUSCO_TITRATION_CLEAN_NOV2025/ALL_NR90/NR90_EUK_busco_titr_summary/NR90_EUK_busco_titr_summary.tsv ./


#manually edit the list in google table

#Complete BUSCOs (C)			   
#Complete and single-copy BUSCOs (S)	   
#Complete and duplicated BUSCOs (D)	   
#Fragmented BUSCOs (F)			   
#Missing BUSCOs (M)			   
#Total BUSCO groups searched	



######## estimate coverage of de novo assemblies, based on spades-approximated estimate in fasta header ;)

for i in *denovo.fasta ; do 
awk '/^>/{ 
        # extract length and coverage from NODE_..._length_XXX_cov_YYY
        match($0, /_length_([0-9]+)_cov_([0-9.]+)/, arr)
        if(arr[1] != "" && arr[2] != ""){
            len = arr[1]
            cov = arr[2]
            total_len += len
            weighted_cov += len * cov
        }
     }
     END{
        if(total_len>0)
            print "Weighted average coverage:", weighted_cov/total_len
        else
            print "No coverage info found."
     }' "$i" ; echo "$i" ; done


#NR90 assemblies are here: /scratch/data2/dzavadska/QUANTITATIVE_DATASET/BUSCO_TITRATION_CLEAN_NOV2025/ALL_NR90

# lengths estimation 
seqkit stats *fasta




####################################################################
####################################################################
####################################################################
# Orthology inference and verification
####################################################################
####################################################################
####################################################################


# OrthoFinder 49 taxa, CD-hit 90%, diamond_ultra_sens , msa

#here /scratch/data2/dzavadska/COMPARATIVE_DATASET/SINGLE_GENE_TREES_2ND_TRY/NR_proteins_90_49taxa
/home/agalvez/programs/OrthoFinder/orthofinder -f Orthofinder -S diamond_ultra_sens -M msa -n "BUSCO_Titration_Results_diamond_ultra_sens_FASTTREE" -a 16 -t 16

#scp nisaba:/scratch/data2/dzavadska/QUANTITATIVE_DATASET/BUSCO_TITRATION_CLEAN_NOV2025/Orthofinder/OrthoFinder/Results_BUSCO_Titration_Results_diamond_ultra_sens_FASTTREE/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv ./

#scp nisaba:/scratch/data2/dzavadska/QUANTITATIVE_DATASET/BUSCO_TITRATION_CLEAN_NOV2025/Orthofinder/OrthoFinder/Results_BUSCO_Titration_Results_diamond_ultra_sens_FASTTREE/Orthogroups/Orthogroups.GeneCount.tsv ./

# Lengths of proteins


for i in *fasta ; do seqkit fx2tab -n -l $i | awk -v f="$i" '{print f"\t"$0}'; done > Ortho_NR90_length.tsv 

for i in *fasta ; do seqkit fx2tab -n -l $i | awk -v f="$i" '{print f"\t"$0}'; done > BUSCO_NR90_length.tsv 


####################################################################
####################################################################
# removing singletons and unassigned short proteins for the final length estimation
####################################################################
####################################################################












