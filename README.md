### MiniProject user guide
Author: Zihan Zheng
Date 3/4/2021


### How to install this to your GitHub repository
git clone with your GitHub repository link

### Usage
#After navigating to the directory of pipeline.py
python pipeline.py /accession file path/GenBank file path/sample information csv file path/blast database name

#e.g.

#To full run for HCMV samples with Betaherpesvirinae subfamily as blast database
python pipeline.py HCMV/accessions HCMV/sequence.gb HCMV/sample_info.csv Betaherpesvirinae HCMV/sequence.fasta
#To run test samples 
python pipeline.py test/accessions test/sequence.gb test/sample_info.csv Betaherpesvirinae test/sequence.fasta

#Test data included in Test.zip in this repo#

### Data preparation

#Save accession numbers to a file for loop usage to have:
#e.g.
cat accessions
#Which output:
#SRR5660030
#SRR5660033
#SRR5660044
#SRR5660045

#Make a csv file for sample information
#The first column has to be the information of sample
#The second column has to be the information of condition
#The third column has to be the information of donor. If all samples are from the same donor, set them as Donor 1.
#The first row has to be column names.
#Columns have to be separated by ,

#Get your virus' FASTA file
#e.g. Get FASTA file of HCMV
wget https://www.ebi.ac.uk/ena/browser/api/fasta/EF999921.1
#Or simply download on the website (NCBI) and put into main directory
#Check if the file is in fasta suffix or change it by:
#e.g.
mv EF999921.1 EF999921.1.fa

#Download the Genbank format of your virus (.gb file)
#Search your virus GenBank id in nucleotide database and download full GenBank information 
#e.g.
#To download the GenBank format of EF999921.1 (i.e. HCMV)
#Go to https://www.ncbi.nlm.nih.gov/nuccore/EF999921
#Click Send to -> Choose Complete Record -> Choose File in Choose Destination ->  Choose GenBank(full) in Format -> Click Create File -> Move the downloaded .gb file into current directory

#Make a local database of just sequences from your target (sub)family
#e.g. make a nr database from Betaherpesvirinae subfamily. 
#When searching sequences of Betaherpesvirinae subfamily in NCBI nucleotide database (by search with: txid10357[Organism:exp]), to have a nr nucleotide database, no source database should be excluded. They are downloaded by 
“click Send to -> Select Complete Record -> Select File in Choose Destination -> Select 
FASTA as format, and Sequence Length in Sort by -> Click Create File”



### Packages and softwares installation

#Many softwares will be installed from conda. Make sure you have had conda on your machine.
#Check if you have conda by showing its current version
#This step is really important
conda -V

#To convert SRA files to paired-end fastq files, fastq-dump is used
#See download page of fastq-dump: https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
#Download sratoolkit (which contains fastq-dump)
cd ~/bin
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.9/sratoolkit.2.10.9-centos_linux64.tar.gz 
#Unzip the file
tar -xzvf sratoolkit.2.10.9-centos_linux64.tar.gz
#Navigate to the bin directory inside the unzipped directory
cd ~/bin/sratoolkit.2.10.9-centos_linux64/bin
pwd
#Copy the path generated from pwd 
#Add it to the environment variable
export PATH=<paste the pwd path here>:$PATH 
source ~/.bashrc
#Check fastq-dump is successfully installed
fastq-dump -V 
#Use fastq-dump to split to paired-end fastq files. Fastq-dump supports to do this without having local SRA but using accession. This saves tones of time from downloading SRA files which is in hundreds of MB large.

#Install biopython for Python code
pip install biopython

#To quantify the TPM of each CDS in each transcriptome, Kallisto is used, and is installed from condo:
conda install -c bioconda kallisto

#Install rstudio in conda to run r script
conda install -c r rstudio

#Install packages used in R script from conda
conda install --channel bioconda r-sleuth
conda install -c r r-tidyverse
conda install -c conda-forge r-devtools
conda install -c bioconda r-annotables
conda install -c r r-data.table

#Install bowtie2 from conda. 
conda install -c bioconda bowtie2

#To using the output from bowtie2, which is SAM, in SPAdes, the format should be transferred to BAM. And it is done by using samtools
#Install samtools from conda
conda install samtools

#Spades installation
conda install SPAdes
#Check if run properly
spades.py --test

#Install blast+ 
#Go to the BLAST download page at: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ and download the compiled program for your machine

#Unzip the zip file

tar -xvzf ncbi-blast-2.11.0+-x64-linux.tar.gz

#Navigate to its bin directory

pwd

#Copy the output path and save into PATH

export PATH=<the copied path>:$PATH
 
#Refresh bashrc file
source ~/.bashrc

### Make sure you installed the tools required, inappropriate installation may cause errors.






 








