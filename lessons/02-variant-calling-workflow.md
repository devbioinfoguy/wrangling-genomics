# Lesson

Automating a workflow
===================

Learning Objectives:
-------------------
#### What's the goal for this lesson?

* Use a series of command line tools to perform a variant calling workflow
* Use a For loop from the previous lesson to help automate repetitive tasks
* Group a series of sequential commands into a script to automate a workflow

## Running a Workflow

To get started with this lesson, ensure you are logged into the cluster and are working
in an interactive session ona compute node. Next, we will need to grab some data from an outside
server using `wget` on the command line.

Make sure you are in the dc_workshop drectory first

    cd /n/regal/datac/$USER/
    wget http://devbioinfoguy.github.io/wrangling-genomics-HPC/data/variant_calling.tar.gz

The file 'variant_calling.tar.gz' is what is commonly called a "tarball", which is
a compressed archive similar to the .zip files we have seen before.  We can decompress
this archive using the command below.

    tar xvf variant_calling.tar.gz

This will create a directory tree that contains some input data (reference genome and fastq files)
and a shell script that details the series of commands used to run the variant calling workflow.

<!--  need to fix the directory structure -->
<pre>
variant_calling
├── ref_genome
│   └── ecoli_rel606.fasta
├── run_variant_calling.sh
└── trimmed_fastq
    ├── SRR097977.fastq
    ├── SRR098026.fastq
    ├── SRR098027.fastq
    ├── SRR098028.fastq
    ├── SRR098281.fastq
    └── SRR098283.fastq
</pre>

Without getting into the details yet, the variant calling workflow will do the following steps

1. Index the reference genome for use by bwa and samtools
2. Align reads to reference genome
3. Convert the format of the alignment to sorted BAM, with some intermediate steps.
4. Calculate the read coverage of positions in the genome
5. Detect the single nucleotide polymorphisms (SNPs)
6. Filter and report the SNP variants in VCF (variant calling format)

We'll first perform the commands in the workflow. We'll next create a script for the 
commands and test this. Finally, we'll modify the script to run on the cluster.
So let's get started.

The first command is to change to our working directory
so the script can find all the files it expects

     cd ./variant_calling

Before we start using software, we have to load the environments for each software
package. On clusters, this is typically done using a module system. For what we need, 
you can execute the following commands

     source new-modules.sh
     module load bwa
     module load samtools
     module load XXX
     
We need to index the refernce genome for bwa and samtools. bwa
and samtools are programs that are pre-installed on our server.

     bwa index data/ref_genome/ecoli_rel606.fasta
     samtools data/ref_genome/ecoli_rel606.fasta

Create output paths for various intermediate and result files
The -p option means mkdir will create the whole path if it
does not exist and refrain from complaining if it does exist

     mkdir -p results/sai
     mkdir -p results/sam
     mkdir -p results/bam
     mkdir -p results/bcf
     mkdir -p results/vcf

In the script, we will eventually loop over all of our files and have the cluster work
on each one in parallel. For now, we're going to work on just one to set up our workflow:

     ls -al data/trimmed_fastq/SRR097977.fastq

<!-- Insert more details on each of the steps!! -->
Align the reads to the reference genome 

    bwa aln genome fastq > SAIfile
    
So

    bwa aln data/ref_genome/ecoli_rel606.fasta data/trimmed_fastq/SRR097977.fastq \
      > results/sai/SRR097977_aligned.sai

Convert the output to the SAM format

    bwa samse genome SAIfile fastq > SAMfile

So

    bwa samse data/ref_genome/ecoli_rel606.fasta results/sai/SRR097977_aligned.sai \
      data/trimmed_fastq/SRR097977.fastq > results/sam/SRR097977_aligned.sam

Convert the SAM file to BAM format

    samtools view -S -b results/sam/$SRR097977_aligned.sam \
      > results/bam/SRR097977_aligned.bam

Sort the BAM file

    samtools sort -f results/bam/SRR097977_aligned.bam \
      results/bam/SRR097977_aligned_sorted.bam

Index the BAM file for display purposes

    samtools index results/bam/SRR097977_aligned_sorted.bam

Do the first pass on variant calling by counting
read coverage

    samtools mpileup -g -f data/ref_genome/ecoli_rel606.fasta \
      results/bam/SRR097977_aligned_sorted.bam > results/bcf/SRR097977_raw.bcf

Do the SNP calling with bcftools

    bcftools view -bvcg results/bcf/SRR097977_raw.bcf > results/bcf/SRR097977_variants.bcf

Filter the SNPs for the final output

    bcftools view results/bcf/SRR097977_variants.bcf | vcfutils.pl \
      varFilter - > results/vcf/SRR097977_final_variants.vcf

That's a lot of work, yes? But you have five more FASTQ files to go...

Exercises:
- Try running this workflow on a different FASTQ file. What did you have to do differently
in order to get this workflow to work?
- Remembering what commands *and* what parameters to type can be pretty daunting. What can
you do to help yourself out in this regard?
- If you were to automate this process, what additional bits of information might you need?


## Automating this Workflow with a Bash Script

The easiest way for you to be able to repeat this process is to capture the steps that
you've performed in a bash script. And you've already learned how to do this in previous
lessons. So...

Exercises:
- Using your command history, create a script file that will repeat these commands
for you. Name your script *run_variant_call_on_file.sh*. Delete your results 
directories, and run your script. Do you get all the proper output files?

One additional command we can put in the top of the script to allow you to see what
is going on is the `set -x` bash command. This debugging tool will print out every
step before it is executed.

- Insert the debugging command in your script and re-run it. How is the output different?
If you're comfortable with how this looks and runs, then comment out this line.

- In order run this workflow on another file, you'll need to make changes. Copy this file,
giving the file a similar name, and make appropriate changes to run on another input
FASTQ file. What did you have to do differently in order to get this workflow to work?

- Knowing techniques that you've learned in previous lessons, what can we do to make this
workflow more friendly to different sets of input files?

- Again, reviewing your two scripts, are there additional commonalities across scripts
or within scripts that we could optimize?


## Granting our Workflow More Flexibility

A couple of changes need to be made to make this script more friendly to both changes
in the workflow and changes in files. 

The first is major change is allowing a change in the filename. Thus at the start of 
the script let's capture an input parameter that must be supplied with the script name.
This input parameter will be the name of the file we want to work on:

     fq="$1"

And we'll add a shortcut to store the location to the genome reference FASTA file:

     # location to genome reference FASTA file
     genome=data/ref_genome/ecoli_rel606.fasta

Now, walking thru the code, we can make some substitutions. To index with bwa and samtools
we can run those commands on the genome variable ($genome) so these value aren't 
static & hardcoded:

     bwa index $genome
     samtools faidx $genome

We'll keep the output paths creation, as it looks fine. (Though really, we could
put results/ in a variable and declare that at the top, so we can change where the
results will be as well. We'll leave that for an optional exercise)

     # make all of our output directories
     mkdir -p results/sai
     mkdir -p results/sam
     mkdir -p results/bam
     mkdir -p results/bcf
     mkdir -p results/vcf

In the script, it is a good idea to use echo for debugging/reporting to the screen

    echo "Processing file $fq ..."

We also need to use one special trick, to extract the base name of the file
(without the path and .fastq extension). We'll assign it
to the $base variable

    # grab base of filename for future naming
    base=$(basename $fq .fastq)
    echo "basename is $base"

Since we've already created our output directories, we can now specify all of our
output files in their proper locations. We will assign various file names to
 variables both for convenience but also to make it easier to see what 
is going on in the sommand below.

    # set up output filenames and locations
    fq=data/trimmed_fastq/$base\.fastq
    sai=results/sai/$base\_aligned.sai
    sam=results/sam/$base\_aligned.sam
    bam=results/bam/$base\_aligned.bam
    sorted_bam=results/bam/$base\_aligned_sorted.bam
    raw_bcf=results/bcf/$base\_raw.bcf
    variants=results/bcf/$base\_variants.bcf
    final_variants=results/vcf/$base\_final_variants.vcf

Our data are now staged.  We now need to change the series of command below
to use our variables so that it will run with more flexibility the steps of the 
analytical workflow>

    # Align the reads to the reference genome
    bwa aln $genome $fq > $sai

    # Convert the output to the SAM format
    bwa samse $genome $sai $fq > $sam

    # Convert the SAM file to BAM format
    samtools view -S -b $sam > $bam

    # Sort the BAM file
    samtools sort -f $bam $sorted_bam

    # Index the BAM file for display purposes
    samtools index $sorted_bam

    # Do the first pass on variant calling by counting read coverage
    samtools mpileup -g -f $genome $sorted_bam > $raw_bcf

    # Do the SNP calling with bcftools
    bcftools view -bvcg $raw_bcf > $variants

    # And finally, filter the SNPs for the final output
    bcftools view $variants | /usr/share/samtools/vcfutils.pl varFilter - > $final_variants

And finally, we need to add our SLURM directives in the *top* of the file so that 
the scheduler knows what resources we need to use in order to run our job on the
compute nodes. And we have to say what software we want to load. 
So the top of the file should look like:

    #!/bin/bash
    #
    #SBATCH -p serial_requeue   # Partition to submit to (comma separated)
    #SBATCH -n 1                # Number of cores
    #SBATCH -N 1                # Ensure that all cores are on one machine
    #SBATCH -t 0-1:00           # Runtime in D-HH:MM (or use minutes)
    #SBATCH --mem 100           # Memory in MB
    #SBATCH -J frog_fastqc      # Job name
    #SBATCH -o fastqc.out       # File to which standard out will be written
    #SBATCH -e fastqc.err       # File to which standard err will be written
    #SBATCH --mail-type=ALL     # Type of email notification: BEGIN,END,FAIL,ALL
    #SBATCH --mail-user=rmf@123.com # Email to which notifications will be sent 

    # set up our software environment...
    source new-modules.sh
    module load bwa
    module load samtools
    module load bcftools

    # now do work...

Test running on one file... blah

XXX fill in details here.
     
Now we need to actually run this to compute the values for all our FASTQ files. And
now is where we'll use the for loop with the power of the cluster. What we'd like to do
is run this script on a compute node for every FASTQ -- pleasantly parallelizing our
workflow.

    for fq in data/trimmed_fastq/*.fastq
    do
      sbatch run_variant_call_on_file.sh $fq
      sleep 1
    done

What you should see on the output of your screen would be the jobIDs that are returned
from the scheduler for each of the jobs that you submitted.

You can see their progress by using the squeue command (though there is a lag of
about 60 seconds between what is happening and what is reported).

Don't forget about the scancel command, should something go wrong and you need to
cancel your jobs.


Exercises:
- Change the script so that one can include an additional variable to point to
 a results directory.
    
<!--
****
**Exercise**
Run the script run_variant_calling.sh
****
-->


