fastq1=$1
output=$2
#picardJars="/gpfs/software/genomics/Picard/1.117/"
#javaParams="/usr/bin/java -Djava.io.tmpdir=/gpfs/ngsdata/scratch/tmp/"
module load MetaPhlAn2/2.6.0
module load Bowtie2/2.3.5.1

##outputdir
prefix=$3
outputdir=${output}
fastqFile=${outputdir}/${prefix}
log="${outputdir}/logs"

#Create output directory.  If it does not exist then it creates it 
mkdir ${outputdir}
mkdir ${log}

#align fastq files
metaphlan2.py $fastq1 --input_type fastq --nproc 8 -t rel_ab_w_read_stats --tax_lev a --no_map --ignore_viruses --ignore_eukaryotes --ignore_archaea -o ${fastqFile}_a.txt

#metaphlan2.py $fastq1 --input_type fastq --nproc 8 -t rel_ab_w_read_stats --tax_lev k --no_map --ignore_viruses --ignore_eukaryotes --ignore_archaea -o ${fastqFile}_k.txt

#metaphlan2.py $fastq1 --input_type fastq --nproc 8 -t rel_ab_w_read_stats --tax_lev p --no_map --ignore_viruses --ignore_eukaryotes --ignore_archaea -o ${fastqFile}_p.txt

#metaphlan2.py $fastq1 --input_type fastq --nproc 8 -t rel_ab_w_read_stats --tax_lev f --no_map --ignore_viruses --ignore_eukaryotes --ignore_archaea -o ${fastqFile}_f.txt

#metaphlan2.py $fastq1 --input_type fastq --nproc 8 -t rel_ab_w_read_stats --tax_lev g --no_map --ignore_viruses --ignore_eukaryotes --ignore_archaea -o ${fastqFile}_g.txt

#metaphlan2.py $fastq1 --input_type fastq --nproc 8 -t rel_ab_w_read_stats --tax_lev s --no_map --ignore_viruses --ignore_eukaryotes --ignore_archaea -o ${fastqFile}_s.txt