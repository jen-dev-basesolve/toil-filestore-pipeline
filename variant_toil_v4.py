from toil.common import Toil
from toil.job import Job
import os
import tempfile 
import requests
import shlex
import subprocess

def main(job,memory="2G", cores=2, disk="3G"):
    ref_job = job.addChildJobFn(ref_genome)
    ref_job.addChildJobFn(index_ref,path=ref_job.rv(0),inputFileID=ref_job.rv(1))            #,tempdir_path=ref_job.rv(1))
    ref_job.addFollowOnJobFn(fastq_files,path=ref_job.rv(0),tempdir_path=ref_job.rv(1))


def ref_genome(job):
    
    container_path='~/singularity_containers'
    
    # parent_dir = '/home/bioinfo/Desktop/Basesolve_jm/Python_R_Toil_bash/Toil_dev/output'
    child_dir = 'dc_workshop_1/data/ref_genome'
    # path = os.path.join(parent_dir,child_dir)
    
    tempdir_path = job.fileStore.getLocalTempDir()
    
    path=tempdir_path + child_dir
    

    # os.makedirs(path,exist_ok=True)
    # Downloading fasta file
    url = "http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz"
    local_filename = 'ecoli_rel606.fasta.gz'
    # the stream=True parameter below

    with requests.get(url, stream=True) as r:
        r.raise_for_status()                            # raises exception if the link not found, or any other error. 
        with open(f"{job.tempDir}/{local_filename}", 'wb') as fi:
            for chunk in r.iter_content(chunk_size=2000): #iter_content is used to break down the content into smaller chunks. Concurrency is also achieved by this. And not the whole of response is loaded into memory at once, hence avoiding memory overflowing.
                fi.write(chunk)
    

    # writing a file into filestore
    fasta = job.fileStore.writeGlobalFile(localFileName=f'{path}/{local_filename}',cleanup=True)
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',inputFileID) 
    
    # extracting fasta file
    comms = ("gunzip %s/%s"%(path,local_filename))      # decomressing the fasta files
    os.system(comms)              # why is read() used over here?
    return parent_dir, inputFileID


# Indexing reference fasta file
def index_ref(job,path,memory="2G", cores=2, disk="3G"):
    comms = "singularity exec %s/bwa_latest.sif bwa index %s/dc_workshop_1/data/ref_genome/ecoli_rel606.fasta"%(path,path)
    process = os.popen(comms)
    process.read()
    job.log("".join(["Following is the path >>>>",str(path)]))

# Downloading trimmed FastQ files for faster operations
def fastq_files(job,tempdir_path,path,memory="4G", cores=4, disk="5G"):
    child_dir = 'dc_workshop_1/data/trimmed_fastq_small'            
    path = os.path.join(tempdir_path,child_dir)
    os.makedirs(path,exist_ok=True)
    
    comms = ("curl -L -o %s/sub.tar.gz https://ndownloader.figshare.com/files/14418248"%path)       # tested, runs
    os.system(comms) #.read()

    comms = ("tar xvf %s/sub.tar.gz -C %s"%(path,path))     # tested, runs
    os.system(comms)

    comms = "cd %s/dc_workshop_1 ; mkdir -p results/sam results/bam results/bcf results/vcf"%tempdir_path     # to be tested in this function
    os.system(comms) #.read()     
    
    job.addChildJobFn(align_fastq,tempdir_path=tempdir_path,path=path)


# Aligning sample sequences to reference genome
def align_fastq(job,path,tempdir_path,memory="2G", cores=2, disk="3G"):
    comms = f"singularity exec {path}/bwa_latest.sif bwa mem {tempdir_path}/dc_workshop_1/data/ref_genome/ecoli_rel606.fasta {tempdir_path}/dc_workshop_1/data/trimmed_fastq_small/sub/SRR2584866_1.trim.sub.fastq {tempdir_path}/dc_workshop_1/data/trimmed_fastq_small/sub/SRR2584866_2.trim.sub.fastq > {tempdir_path}/dc_workshop_1/results/sam/SRR2584866.aligned.sam"
    os.system(comms)
    # process.read()
    job.addChildJobFn(convrt_sort,path=path,tempdir_path=tempdir_path)

# Converting SAM file to BAM file using view option in samtools
# Sorting bam files
def convrt_sort(job,path,tempdir_path,memory="2G", cores=2, disk="3G"):
    # convert sam file to bam
    comms = f"singularity exec {path}/samtools_latest.sif samtools view -S -b {tempdir_path}/dc_workshop_1/results/sam/SRR2584866.aligned.sam > {tempdir_path}/dc_workshop_1/results/bam/SRR2584866.aligned.bam"
    os.system(comms) #.read()

    # sort bam file
    comms = f"singularity exec {path}/samtools_latest.sif samtools sort -o {tempdir_path}/dc_workshop_1/results/bam/SRR2584866.aligned.sorted.bam {tempdir_path}/dc_workshop_1/results/bam/SRR2584866.aligned.bam"
    os.system(comms) #.read()

    job.addChildJobFn(variant_call,path=path,tempdir_path=tempdir_path)

# Variant calling
def variant_call(job,path,tempdir_path, memory="2G", cores=2, disk="3G"):
    comms = f"singularity exec {path}/bcftools_latest.sif bcftools mpileup -O b -o {tempdir_path}/dc_workshop_1/results/bcf/SRR2584866_raw.bcf -f {tempdir_path}/dc_workshop_1/data/ref_genome/ecoli_rel606.fasta {tempdir_path}/dc_workshop_1/results/bam/SRR2584866.aligned.sorted.bam"
    os.system(comms)  #.read()      # works, tested.

    # Detecting Single Nucleotide Variants(SNV) from VCF file.
    # ploidy - number of chromosome sets in nucleus
    comms = f"singularity exec {path}/bcftools_latest.sif bcftools call --ploidy 1 -m -v -o {tempdir_path}/dc_workshop_1/results/vcf/SRR2584866_variants.vcf {tempdir_path}/dc_workshop_1/results/bcf/SRR2584866_raw.bcf"
    os.system(comms)     #.read()

    # Filter and report the SNV variants in variant calling format (VCF)
    comms = f"singularity exec {path}/bcftools_latest.sif vcfutils.pl varFilter {tempdir_path}/dc_workshop_1/results/vcf/SRR2584866_variants.vcf  > {tempdir_path}/dc_workshop_1/results/vcf/SRR2584866_final_variants.vcf"
    os.system(comms)     #.read()


    # saving the local copy of file
    # inputFileID = job.fileStore.writeGlobalFile(localFileName=f'{path}/dc_workshop_1/results/vcf/SRR2584866_variants.vcf',cleanup=True)
    # print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',inputFileID)


    inputFileID=f'{tempdir_path}/dc_workshop_1/results/vcf/SRR2584866_variants.vcf'
    
    with job.fileStore.readGlobalFileStream(inputFileID, encoding='utf-8') as fi:
            with job.fileStore.writeGlobalFileStream(encoding='utf-8') as (fo, outputFileID):
                fo.write(fi.read())                  # readGlobalTempFile can read only temp files and not local files
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',outputFileID)

    job.addChildJobFn(index_fastq,path=path,tempdir_path=tempdir_path)
    return outputFileID



# Indexing for visualization
def index_fastq(job,path,tempdir_path,memory="2G", cores=2, disk="3G"):
    comms = f"singularity exec {path}/samtools_latest.sif samtools index {tempdir_path}/dc_workshop_1/results/bam/SRR2584866.aligned.sorted.bam"
    os.system(comms)     #.read()
    

# Visualising SNV using t_view
# def t_view():
#   path = '/home/bioinfo/singularity/variant_analysis'
#   comm = f"singularity exec {path}/samtools_latest.sif samtools tview {path}/dc_workshop_1/results/bam/SRR2584866.aligned.sorted.bam {path}/dc_workshop_1/data/ref_genome/ecoli_rel606.fasta"
#   os.popen(comm).read()


if __name__ == "__main__":
    jobstore: str = tempfile.mkdtemp("tutorial_staging")
    os.rmdir(jobstore)
    options = Job.Runner.getDefaultOptions(jobstore)
    options.logLevel = "Debug"
    options.clean = "always"
    job = Job.wrapJobFn(main)

    with Toil(options) as toil:
        toil.start(job)