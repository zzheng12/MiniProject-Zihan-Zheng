from Bio import SeqIO
import sys
import os


if len(sys.argv) != 6:
    raise Exception('Missing input argument(s)\npython pipeline.py <accession file path> <GenBank file path> <sample information csv file path> <blast database name>')
#Raise error if number of input arguments is not 4
  
accession = sys.argv[1]
input1 = sys.argv[2]
sample_info = sys.argv[3]
db_name = sys.argv[4]
seqdb = sys.argv[5]
#Read input file paths from stdin 
     
    

if not os.path.isfile('./'+accession.strip()):
    raise Exception('Given accession file does not exist.')
if not os.path.isfile('./'+input1.strip()):
    raise Exception('Given GenBank file does not exist.')
if not os.path.isfile('./'+sample_info.strip()):
    raise Exception('Given csv file does not exist.')
#Check the existence of input

os.system('cat '+accession+' | while read line; do fastq-dump ${line} --split-files --gzip; done')
#Run fastq-dump to split to paired-end fastq files.
#Make new directory if not exist
if not os.path.isdir('./miniProject_Zihan_Zheng'):
    
    os.mkdir('./miniProject_Zihan_Zheng')
    #Make output directory
#Make new directory if not exist   
   
seqfa = './miniProject_Zihan_Zheng/sequence.fa'
seqfa_dna = './miniProject_Zihan_Zheng/sequence_dna.fa'
#Generate output file path 

miniProject = './miniProject.log'
#Generate log file path  
   
with open(seqfa, 'w') as op, \
    open(seqfa_dna, 'w') as op2:
    ct = 0
    for j in SeqIO.parse(open(input1.strip(), 'r'), 'genbank'):
        print("Dealing with GenBank record %s"%(j.id))
        for feature in j.features:
            if feature.type == 'CDS':
                ct += 1
        print("Current number of CDS: %i"%(ct))
         
        op.write('>%s\n%s\n'%(j.id, j.seq.transcribe()))
        op2.write('>%s\n%s\n'%(j.id, j.seq))
        #Write transcriptome into the output
    #Read GenBank input and count the number of CDS
        
os.rename(seqfa, './miniProject_Zihan_Zheng/'+j.id+'.fa')
#Change the output file name to the id       
    
with open(miniProject, 'w') as log:
    log.write('The HCMV genome (EF999921) has %i CDS\n\n'%(ct))
#Write number of CDS to log file 
    

os.system("kallisto index -i miniProject_Zihan_Zheng/transcripts.idx miniProject_Zihan_Zheng/" + j.id + ".fa")
# Building an index by kallisto
# where:
# -i indicates index file path


os.mkdir('./miniProject_Zihan_Zheng/3_kallisto')
#Make kallisto directory


#Quantify the TPM of each CDS in each transcriptome by kallisto
os.system("cat "+accession+" | while read line; do kallisto quant -i miniProject_Zihan_Zheng/transcripts.idx -o miniProject_Zihan_Zheng/3_kallisto/${line} -b 100 ${line}_1.fastq.gz ${line}_2.fastq.gz; done")
#abundaces are reported in ‘est_counts’ and ‘TPM’ in abundance.tsv 

#Initialize header for sleuth output
with open(miniProject, 'a') as op:
    op.write('\ntarget_id\ttest_stat\tpval\tqval\n')    
#Run r script to run sleuth to get details for each significant transcript (FDR < 0.05)
os.system("Rscript diff_expressed_gene.R")



os.mkdir('./miniProject_Zihan_Zheng/4_bowtie')
#Make directory for bowtie2 data
os.system("bowtie2-build ./miniProject_Zihan_Zheng/sequence_dna.fa miniProject_Zihan_Zheng/4_bowtie/bowtie2_index")
#Do indexing to run bowtie2

os.system("cat "+accession+" | while read line; do bowtie2 -x miniProject_Zihan_Zheng/4_bowtie/bowtie2_index -1 ${line}_1.fastq.gz -2 ${line}_2.fastq.gz -S miniProject_Zihan_Zheng/4_bowtie/reads_aligned_${line}.sam >> miniProject_Zihan_Zheng/4_bowtie/bowtie2.log; done")
#Do mapping 
# -x indicates the name of and path to the index
# -1 -2 indicates the name of and path to aligned reads
# -S indicates that we want the mapping results in a SAM format file with given output name. 
# The STDOUT output of the command is stored into .log file.
# The output SAM file is unsorted


os.system("cat "+accession+" | while read line; do samtools view -bSF12 -o miniProject_Zihan_Zheng/4_bowtie/reads_aligned_${line}.bam  miniProject_Zihan_Zheng/4_bowtie/reads_aligned_${line}.sam; done")
#Convert SAM file to BAM as well as filter excluding 4 flag and 8 flag which are the mapped reads

# -b indicates the output should be BAM
# -S indicates the input is in SAM
# -F12 indicates that the 4 flag and 8 flag will not be included


os.system("cat "+accession+" | while read line; do samtools fastq -1 miniProject_Zihan_Zheng/4_bowtie/${line}_mapped_1.fastq -2 miniProject_Zihan_Zheng/4_bowtie/${line}_mapped_2.fastq miniProject_Zihan_Zheng/4_bowtie/reads_aligned_${line}.bam; done")
#Covert BAM with only mapped reads to fastq

#Write the number of reads in each transcriptome before and after the Bowtie2 mapping
#Get the file path into a list
bef_paths = {}
aft_paths = {}
csv_infor = {}
col_name = [] 
# Note that by rule for this script:
#    col_name[0]: sample accession
#    col_name[1]: condition
#    col_name[2]: donor (if all from same donor, all of them are donor1)


with open(sample_info, 'r') as si:
    header = si.readline().strip().split(',')
    for col in header:
        csv_infor[col] = []
        col_name.append(col)
    for line in si:
        element = line.strip().split(',')
        count = 0
        for ele in element:
            csv_infor[col_name[count]].append(ele)
            count += 1   
#Get all infor from sample_infor and save into a dictionary

for sample_name in csv_infor[col_name[0]]:
    bp = './'+sample_name+'_1.fastq.gz'
    ap = './miniProject_Zihan_Zheng/4_bowtie/'+sample_name+'_mapped_1.fastq'
    os.system('gunzip -c '+bp+' > '+bp[:len(bp)-3])
    
    num_read_bp = 0
    for bp_line in SeqIO.parse(open(bp[:len(bp)-3].strip(), 'r'), 'fastq'):
        num_read_bp += 1
    #Count the number of reads for file before mapping
    num_read_ap = 0
    for ap_line in SeqIO.parse(open(ap.strip(), 'r'), 'fastq'):
        num_read_ap += 1 
    #Count the number of reads for file after mapping
    bef_paths[sample_name] = num_read_bp
    aft_paths[sample_name] = num_read_ap
    #Store counts into corresponding dictionaries' value
#Read samples from dictionary, get the file paths and number of reads for them

with open(miniProject, 'a') as log:
    log.write('\n')
    ind = 0
    while ind < len(bef_paths):
        log.write(csv_infor[col_name[2]][ind]+' ('+col_name[1][ind]+') had '+str(bef_paths[csv_infor[col_name[0]][ind]])+' read pairs before Bowtie2 filtering and '+str(aft_paths[csv_infor[col_name[0]][ind]])+' read pairs after.\n\n\n')
        ind += 1
#Append information with given format

os.mkdir('./miniProject_Zihan_Zheng/5_assembly')
#Make directory for assembly
os.system("samtools merge miniProject_Zihan_Zheng/5_assembly/merge.bam miniProject_Zihan_Zheng/4_bowtie/reads_aligned_*.bam")
print('All mapped reads have been merged into one')
#Merge all mapped reads from transcriptomes 

os.system("samtools fastq -1 miniProject_Zihan_Zheng/5_assembly/merge_mapped_1.fastq -2 miniProject_Zihan_Zheng/5_assembly/merge_mapped_2.fastq miniProject_Zihan_Zheng/5_assembly/merge.bam")
print('Bam filehas been transferred into splited fastq file')
#Transfer merged bam file into splited fastq file


os.system("spades.py -k 127 -t 16 --only-assembler -1 miniProject_Zihan_Zheng/5_assembly/merge_mapped_1.fastq -2 miniProject_Zihan_Zheng/5_assembly/merge_mapped_2.fastq -o miniProject_Zihan_Zheng/5_assembly/Spades")
#Do RNA SPAdes 
#os.system("spades.py -k 21,33,55,77 -t 16 --only-assembler -1 miniProject_Zihan_Zheng/5_assembly/merge_mapped_1.fastq -2 miniProject_Zihan_Zheng/5_assembly/merge_mapped_2.fastq -o miniProject_Zihan_Zheng/5_assembly/Spades")
# -k indicates parameters to seed the algorithm for de novo assembly. 
# They are k-mer sizes SPAdes uses to build de Bruijn graphs and find the best assembly.
# -t indicates the number of processors
# -- only-assembler indicates only do assembly
# -1, -2 indicates the input forward and reverse fastq files
# -o indicates the output directory



os.system("echo 'spades.py -k 127 -t 16 --only-assembler -1 miniProject_Zihan_Zheng/5_assembly/merge_mapped_1.fastq -2 miniProject_Zihan_Zheng/5_assembly/merge_mapped_2.fastq -o miniProject_Zihan_Zheng/5_assembly/Spades' >> miniProject.log")
#Save SPAdes command to log file
#os.system("echo 'spades.py -k 21,33,55,77 -t 16 --only-assembler -1 miniProject_Zihan_Zheng/5_assembly/merge_mapped_1.fastq -2 miniProject_Zihan_Zheng/5_assembly/merge_mapped_2.fastq -o miniProject_Zihan_Zheng/5_assembly/Spades' >> miniProject.log")


ip_path = './miniProject_Zihan_Zheng/5_assembly/Spades/scaffolds.fasta'
op2_path = './miniProject_Zihan_Zheng/longest_contig.fasta' #the output path for the longest contig
#Set input path
 
count_1000 = 0
bp = 0
max_len = 0
max_contig = ''
max_seq = ''
current_contig = ''
contigs = {}
with open(ip_path, 'r') as ip, \
    open(miniProject, 'a') as op, \
        open(op2_path, 'w') as op2:
    for line in ip:
        if line.startswith('>'): #Store current contig name
            current_contig = line[1:].strip().strip('\n')
        else:
            if current_contig not in contigs: #indicate current contig has not been added into contigs dict
                contigs[current_contig] = line.strip().strip('\n')
            else: #indicate current contig have had line(s) in contigs dict
                contigs[current_contig] += line.strip().strip('\n')
    for k, v in contigs.items():
        if len(v) > 1000:
            count_1000 += 1 #count if the length of line is > 1000
            bp += len(v)
            if len(v) > max_len:
                max_len = len(v)
                max_contig = k
                max_seq = v
    op.write('\n\nThere are %i contigs > 1000 bp in the assembly.\n'%(count_1000))
    op.write('\nThere are %i bp in the assembly.\n\n'%(bp))
    op2.write('>%s\n%s'%(max_contig, max_seq))
    print('There are %i contigs > 1000 bp in the assembly.\n'%(count_1000))
    print('There are %i bp in the assembly.\n'%(bp))
    print('The longest contig is %s in %i length. Its sequence is %s'
          %(max_contig, max_len, max_seq))
#Calculate the number of contigs with a length > 1000


print('Making blast database')
makeblast_cmd = 'makeblastdb -in '+seqdb+' -out miniProject_Zihan_Zheng/blastdb/' + db_name.strip() + ' -title '+ db_name.strip() + ' -dbtype nucl'
# -in indicates input file path
# -out indicates output file path 
# -title indicates database name
# -dbtype indicates the type of database


os.system(makeblast_cmd)
#Make blast database

db = './miniProject_Zihan_Zheng/blastdb/'+db_name.strip()
ip_path = './miniProject_Zihan_Zheng/longest_contig.fasta'
op_path = './miniProject_Zihan_Zheng/blastResult_for_longest_fasta.csv'
blast_cmd = 'blastn -query ' + ip_path + ' -db ' + db + ' -out ' + op_path + ' -outfmt "10 sacc pident length qstart qend sstart send bitscore evalue stitle"'
##Define file paths
#-query indicates input fasta file path
#-db indicates database path 
#-out indicates output csv file path
#-outfmt indicates output format. 10 indicates csv format. 
#sacc : Subject accession, pident: Percent identity, length = Alignment length, 
#qstart : Start of alignment in query, qend : End of alignment in query, 
#sstart : Start of alignment in subject, send : End of alignment in subject,
#bitscore : Bit score, evalue : E-value, stitle : Subject Title.


os.system(blast_cmd)
print('Blastn has been successfully run. See the output in miniProject_Zihan_Zheng/blastResult_for_longest_fasta.csv')
#Run blastn

with open(op_path, 'r') as csv, \
    open(miniProject, 'a') as log:
        log.write('\nsacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n')
        count = 0
        for line in csv:
            if count < 10:
                log.write('\t'.join(line.strip().split(','))+'\n')
                count += 1
            else:
                break
#Save top 10 hits into log file as final result
print('Top 10 hits have been written into miniProject.log')
#Print final result



    

    

