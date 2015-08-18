# Lesson QC of Sequence Read Data

Quality Control of NGS Data
===================

Learning Objectives:
-------------------
#### What's the goal for this lesson?

* Understand how the FastQ format encodes quality
* Be able to evaluate a FastQC report
* Use Trimmommatic to clean FastQ reads
* Use a For loop to automate operations on multiple files


## Details on the FASTQ format

Although it looks complicated  (and maybe it is), its easy to understand the [fastq](https://en.wikipedia.org/wiki/FASTQ_format) format with a little decoding. Some rules about the format include...

|Line|Description|
|----|-----------|
|1|Always begins with '@' and then information about the read|
|2|The actual DNA sequence|
|3|Always begins with a '+' and sometimes the same info in line 1|
|4|Has a string of characters which represent the quality scores; must have same number of characters as line 2|

so for example in our data set, one complete read is:
```
$ head -n4 SRR098281.fastq 
@SRR098281.1 HWUSI-EAS1599_1:2:1:0:318 length=35
CNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+SRR098281.1 HWUSI-EAS1599_1:2:1:0:318 length=35
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
```
This is a pretty bad read. 

Notice that line 4 is:
```
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
```
As mentioned above, line 4 is a encoding of the quality. In this case, the code is the [ASCII](https://en.wikipedia.org/wiki/ASCII#ASCII_printable_code_chart) character table. According to the chart a '#' has the value 35 and '!' has the value 33. If only it were that simple. There are actually several historical differences in how Illumina and other players have encoded the scores. Heres the chart from wikipedia:

```
  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  |                         |    |        |                              |                     |
 33                        59   64       73                            104                   126
  0........................26...31.......40                                
                           -5....0........9.............................40 
                                 0........9.............................40 
                                    3.....9.............................40 
  0.2......................26...31........41                              

 S - Sanger        Phred+33,  raw reads typically (0, 40)
 X - Solexa        Solexa+64, raw reads typically (-5, 40)
 I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
 J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
     with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
     (Note: See discussion above).
 L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
 ```
 So using the Illumina 1.8 encouding, which is what you will mostly see from now on, our first c is called with a Phred score of 0 and our Ns are called with a score of 2. Read quality is assessed using the Phred Quality Score.  This score is logarithmically based and the score values can be interpreted as follows:

|Phred Quality Score |Probability of incorrect base call |Base call accuracy|
|:-------------------|:---------------------------------:|-----------------:|
|10	|1 in 10 |	90%|
|20	|1 in 100|	99%|
|30	|1 in 1000|	99.9%|
|40	|1 in 10,000|	99.99%|
|50	|1 in 100,000|	99.999%|
|60	|1 in 1,000,000|	99.9999%|

## FastQC
FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

The main functions of FastQC are
* Import of data from BAM, SAM or FastQ files (any variant)
* Providing a quick overview to tell you in which areas there may be problems
* Summary graphs and tables to quickly assess your data
* Export of results to an HTML based permanent report
* Offline operation to allow automated generation of reports without running the interactive application


## Running FASTQC
###A. Stage your data

Create a working directory for your analysis
    
    cd
    mkdir dc_workshop

Give it three three subdirectories

    mkdir dc_workshop/data
    mkdir dc_workshop/docs
    mkdir dc_workshop/results

The sample data we will be working with is in a hidden directory, we will
move it to our working directory

    mv ~/.dc_sampledata_lite/untrimmed_fastq/ ~/dc_workshop/data/

###B. Run FastQC  

Before we run FastQC, let's start an interactive session on the cluster

	srun -p interact --pty --mem 500 -t 0-06:00 /bin/bash

*An interactive session is a very useful to test tools, workflows, run jobs that open new interactive windows (X11-forwarding) and so on.*

Once your interactive job starts, notice that the command prompt no longer says rclogin; this is because we are not working on the login node any more.

    cd ~/dc_workshop/data/untrimmed_fastq/  

To run the FastQC program, we first need to load the appropriate module.

	module load centos6/fastqc-0.10.1

Once a module for a tool is loaded, you have essentially made it directly available to you like any other basic UNIX command.

FastQC will accept multiple file names as input, so we can use the *.fastq wildcard.

	fastqc *.fastq

*Did you notice how each file was processed serially? How do we speed this up?*

Exit the interactive session and start a new one with 3 cores, and use the multi-threading funcionality of FastQC to run 3 jobs at once.

	exit      #exit the current interactive session
	
	srun -p interact -n 3 --pty --mem 500 -t 0-06:00 /bin/bash      #start a new one with 3 cpus (-n 3)
	
	module load centos6/fastqc-0.10.1      #you'll have to reload the module for the new session
	
	fastqc -t 3 *.fastq       #note the extra parameter we specified for 3 threads


Now, let's create a home for our results

    mkdir ~/dc_workshop/results/fastqc_untrimmed_reads

...and move them there (recall, we are still in '~/dc_workshop/data/untrimmed_fastq/')

    mv *.zip ~/dc_workshop/results/fastqc_untrimmed_reads/
    mv *.html ~/dc_workshop/results/fastqc_untrimmed_reads/

###C. Results
   
	~/dc_workshop/results/fastqc_untrimmed_reads/
	ls 
   
The zip files need to be unpacked with the 'unzip' program.  If we try to do them all
at once.

    unzip *.zip

Did it work? No, because 'unzip' expects to get only one zip file.  Welcome to the real world.
We *could* do each file, one by one, but what if we have 500 files?  There is a smarter way.
We can save time by using a simple shell 'for loop' to iterate through the list of files in *.zip.
After you type the first line, you will get a special '>' prompt to type next next lines.  
You start with 'do', then enter your commands, then end with 'done' to execute the loop.

    $ for zip in *.zip
    > do
    > unzip $zip
    > done

Note that, in the first line, we create a variable named 'zip'.  After that, we call that variable
with the syntax $zip.  $zip is assigned the value of each item (file) in the list *.zip, once for each
iteration of the loop.

This loop is basically a simple program.  When it runs, it will run unzip 
once for each file (whose name is stored in the $zip variable). The contents of each file will 
be unpacked into a separate directory by the unzip program.

The for loop is intepreted as a multipart command.  If you press the up arrow
on your keyboard to recall the command, it will be shown like so:

    for zip in *.zip; do echo File $zip; unzip $zip; done

When you check your history later, it will help your remember what you did!

#### Document your work

To save a record, let's cat all fastqc summary.txts into one full_report.txt and move this to ~/dc_workshop/docs. 
You can use wildcards in paths as well as file names.  Do you remember how we said 'cat' is
really meant for concatenating text files?
    
    cat */summary.txt > ~/dc_workshop/docs/fastqc_summaries.txt


##How to clean reads using *Trimmomatic*
###A detailed explanation of features

Once we have an idea of the quality of our raw data, it is time to trim away adapters and filter out poor quality score reads. To accomplish this task we will use [*Trimmomatic*](http://www.usadellab.org/cms/?page=trimmomatic).

*Trimmomatic* is a java based program that can remove sequencer specific reads and nucleotides that fall below a certain threshold. *Trimmomatic* can be multithreaded to run quickly. 

Because *Trimmomatic* is java based, it is run using the command:

**_java jar trimmomatic-0.30.jar_**

What follows this are the specific commands that tells the program exactly how you want it to operate. *Trimmomatic* has a variety of options and parameters:

* **_-threads_** How many processors do you want *Trimmomatic* to run with?
* **_SE_** or **_PE_** Single End or Paired End reads?
* **_-phred33_** or **_-phred64_** Which quality score do your reads have?
* **_SLIDINGWINDOW_** Perform sliding window trimming, cutting once the average quality within the window falls below a threshold.
* **_LEADING_** Cut bases off the start of a read, if below a threshold quality.
* **_TRAILING_** Cut bases off the end of a read, if below a threshold quality.
* **_CROP_** Cut the read to a specified length.
* **_HEADCROP_** Cut the specified number of bases from the start of the read.
* **_MINLEN_** Drop an entire read if it is below a specified length.
* **_TOPHRED33_** Convert quality scores to Phred-33.
* **_TOPHRED64_** Convert quality scores to Phred-64.

A generic command for *Trimmomatic* looks like this:

**_java jar trimmomatic-0.30.jar SE -thr _**

A complete command for *Trimmomatic* will look something like this:

**_java jar trimmomatic-0.30.jar SE -threads 4 -phred64 SRR_0156.fastq SRR_1056_trimmed.fastq ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20 _**

This command tells *Trimmomatic* to run on a Single End file (``SRR_0156.fastq``, in this case), the output file will be called ``SRR_0156_trimmed.fastq``,  there is a file with Illumina adapters called ``SRR_adapters.fa``, and we are using a sliding window of size 4 that will remove those bases if their phred score is below 20.


## Exercise - Running Trimmomatic

Go to the untrimmed fastq data location:

     cd ~/dc_workshop/data/untrimmed_fastq

The command line incantation for trimmomatic is more complicated.  This is where what you have been learning about accessing your command line history will start to become important.

Let's load the trimmomatic module:

	module load centos6/Trimmomatic-0.30

The general form of the command on this cluster is:

    java -jar $TRIMMOMATIC/trimmomatic-0.30.jar SE -phred64 inputfile outputfile OPTION:VALUE... # DO NOT RUN THIS

'java -jar' calls the Java program, which is needed to run trimmotic, which lived in a 'jar' file (trimmomatic-0.30.jar), a special kind of java archive that is often used for programs written in the Java programing language.  If you see a new program that ends in '.jar', you will know it is a java program that is executed 'java -jar program name'.  The 'SE' argument is a keyword that specifies we are working with single-end reads.

NOTE: The "$TRIMMOMATIC" variable denoting the path to the .jar file is created on this cluster for ease of use, and is specific to this set up.

The next two arguments are input file and output file names. These are then followed by a series of options. The specifics of how options are passed to a program are different depending on the program. You will always have to read the manual of a new program to learn which way it expects its command-line aruments to be composed.


So, for the single fastq input file 'SRR098283.fastq', the command would be:

    java -jar $TRIMMOMATIC/trimmomatic-0.30.jar SE -phred33 SRR098283.fastq \
    SRR098283.fastq_trim.fastq SLIDINGWINDOW:4:20 MINLEN:20

	TrimmomaticSE: Started with arguments: -phred33 SRR098283.fastq SRR098283.fastq_trim.fastq SLIDINGWINDOW:4:20 MINLEN:20
	Input Reads: 250000 Surviving: 193900 (77.56%) Dropped: 56100 (22.44%)
	TrimmomaticSE: Completed successfully

So that worked and we have a new fastq file.

    ls SRR098283*
    SRR098283.fastq  SRR098283.fastq_trim.fastq

Now we know how to run trimmomatic but there is some good news and bad news.  
One should always ask for the bad news first.  Trimmomatic only operates on 
one input file at a time and we have more than one input file.  The good news?
We already know how to use a for loop to deal with this situation.

    for infile in *.fastq
    >do
    >outfile=$infile\_trim.fastq
    >java -jar $TRIMMOMATIC/trimmomatic-0.30.jar SE $infile $outfile SLIDINGWINDOW:4:20 MINLEN:20
    >done

Do you remember how the first specifies a variable that is assigned the value of each item in the list in turn?  We can call it whatever we like.  This time it is called infile.  Note that the third line of this for loop is creating a second variable called outfile.  We assign it the value of $infile
with '_trim.fastq' appended to it.  The '\' escape character is used so the shell knows that whatever
follows \ is not part of the variable name $infile.  There are no spaces before or after the '='.







