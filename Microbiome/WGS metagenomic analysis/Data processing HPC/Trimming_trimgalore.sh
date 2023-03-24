fastq1=$1
output=$2
picardJars="/gpfs/software/genomics/Picard/1.117/"
javaParams="/usr/bin/java -Djava.io.tmpdir=/gpfs/ngsdata/scratch/tmp/"
module load TrimGalore/v0.6.5

##outputdir
prefix=$3
outputdir=${output}
fastqFile=${outputdir}/${prefix}
log="${outputdir}/logs"

#Create output directory.  If it does not exist then it creates it 
mkdir ${outputdir}
mkdir ${log}

#Trim fastq files

trim_galore $fastq1 -j 8 --length 35 --quality 25 -o /gpfs/projects/mharis_Lab/Data/microbiome/trimmed