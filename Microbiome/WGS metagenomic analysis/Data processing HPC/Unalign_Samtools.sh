bam1=$1
output=$2
picardJars="/gpfs/software/genomics/Picard/1.117/"
javaParams="/usr/bin/java -Djava.io.tmpdir=/gpfs/ngsdata/scratch/tmp/"
module load SAMtools/1.14

##outputdir
prefix=$3
outputdir=${output}
bamFile=${outputdir}/${prefix}
log="${outputdir}/logs"

#Create output directory.  If it does not exist then it creates it 
mkdir ${outputdir}
mkdir ${log}

#Create bam file with unmapped read whose mate is mapped.

samtools view -u  -f 4 -F 264 $bam1 -o ${bamFile}.tmp1.bam

#Create bam file with umapped read whos mate is unmapped.
samtools view -u  -f 8 -F 260 $bam1 -o ${bamFile}.tmp2.bam

#Create bam file with both reads of the pair unmapped.
samtools view -u  -f 12 -F 256 $bam1 -o ${bamFile}.tmp3.bam

#Create merged unmapped bam.
samtools merge ${bamFile}.unaligned.bam ${bamFile}.tmp1.bam ${bamFile}.tmp2.bam ${bamFile}.tmp3.bam

#Create sorted merged unmapped bam.
samtools sort -@ 32 -m 1G -o ${bamFile}.unaligned.sorted.bam ${bamFile}.unaligned.bam

#Index the sorted unmapped bam.
samtools index ${bamFile}.unaligned.sorted.bam

#Create fastq files from unaligned bam.
samtools fastq  -@ 16 -o unaligned_fastq/${prefix}_unaligned.fastq ${bamFile}.unaligned.sorted.bam