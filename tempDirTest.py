from toil.common import Toil
from toil.job import Job
import os
import tempfile 
import requests
import shlex
import subprocess, shlex




container_path='/home/bioinfo/singularity_containers'

def path_generator(job,file_id,fileName):
    file_path = job.fileStore.readGlobalFile(fileStoreID=file_id,userPath=fileName)
    return file_path

def makeFileId(job,filePath):
    fileId = job.fileStore.writeGlobalFile(filePath)
    return fileId


def download_url(url, work_dir='.', name=None):
    """
    Downloads URL, can pass in file://, http://, s3://, or ftp://
    If downloading S3 URLs, the S3AM binary must be on the PATH

    :param str url: URL to download from
    :param str work_dir: Directory to download file to
    :param str name: Name of output file, if None, basename of URL is used
    :param str s3_key_path: Path to 32-byte encryption key if url points to S3 file that uses SSE-C
    :return: Path to the downloaded file
    :rtype: str
    """
    file_path = os.path.join(work_dir, name) if name else os.path.join(work_dir, os.path.basename(url))
    subprocess.check_call(['curl', '-fs', '--retry', '5', '--create-dir', url, '-o', file_path])    # check_call is similar to popen.wait()
    assert os.path.exists(file_path)        # checks if the file path exists.
    return file_path


def main(job,memory="2G", cores=2, disk="3G"):
    ref_job = job.addChildJobFn(ref_genome)
    indexedFa_job = ref_job.addChildJobFn(index_ref,fasta_id=ref_job.rv())
    trimmedFastqJob = ref_job.addFollowOnJobFn(trimmed_fastq_files)
    alignedFastqJob = trimmedFastqJob.addChildJobFn(align_fastq,indexedFilesId=indexedFa_job.rv(),trimmedFastqPath=trimmedFastqJob.rv())
    convrt_sort_job = alignedFastqJob.addChildJobFn(convrt_sort,samFileId=alignedFastqJob.rv())
    variant_call_job = convrt_sort_job.addChildJobFn(variant_call,indexedFilesDict=indexedFa_job.rv(),sortedBamId=convrt_sort_job.rv())
    return convrt_sort_job.rv()

def ref_genome(job):    
    """
    Downloads fasta files and gunzips them and returns fileID

    Returns:
        _type_: _description_
    """
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz"
    local_filename = 'ecoli_rel606.fna.gz'
    "Downloading fna.gz files and storing them in global filestore"
    work_dir = job.fileStore.getLocalTempDir()
    fpath = download_url(url,work_dir=work_dir,name=local_filename)     # returns filepath of the downloaded files
    # fastq_id = job.fileStore.writeGlobalFile(fpath)                     # save local file to filestore and return its filestore ID
    subprocess.run(['gunzip',fpath])    
    # writing a file into filestore
    fpath = os.path.join(work_dir,'ecoli_rel606.fna')
    fasta_id = job.fileStore.writeGlobalFile(localFileName=fpath)
    return fasta_id


# Indexing reference fasta file
def index_ref(job,fasta_id,memory="2G", cores=2, disk="3G"):
    """
    Indexes fasta file, returns fileID of the indexed file

    Args:
        job (_type_): _description_
        fasta_id (_type_): _description_
        memory (str, optional): _description_. Defaults to "2G".
        cores (int, optional): _description_. Defaults to 2.
        disk (str, optional): _description_. Defaults to "3G".

    Returns:
        _type_: _description_
    """

    indexedFilesPath = os.path.join(job.fileStore.getLocalTempDir(), 'indexedFilesPath')
    os.makedirs(indexedFilesPath,exist_ok=True)
    fastaFile = job.fileStore.readGlobalFile(fasta_id, os.path.join(indexedFilesPath, 'ecoli_rel606.fasta'))       
    comms = f"singularity exec {container_path}/bwa_latest.sif bwa index {fastaFile}"
    args= shlex.split(comms)
    subprocess.run(args)
    indexFilesNameList = subprocess.run(['ls',indexedFilesPath],capture_output=True,text=True).stdout.splitlines()

    indexedFilesIdDict = {}
    for fi in indexFilesNameList:
        fileid = makeFileId(job,os.path.join(indexedFilesPath,fi))
        indexedFilesIdDict[fi]=fileid
    return indexedFilesIdDict


# Downloading trimmed FastQ files
def trimmed_fastq_files(job,memory="4G", cores=4, disk="5G"):
    """
    _summary_ : Downloads trimmed fastq files, adds align fastq files as the child job

    Args:
        job (_type_): _description_
        memory (str, optional): _description_. Defaults to "4G".
        cores (int, optional): _description_. Defaults to 4.
        disk (str, optional): _description_. Defaults to "5G".
    """

    # Download trimmed fastq files into a temporary directory and return the path to that temporary directory     
    work_dir = os.path.join(job.fileStore.getLocalTempDir(),'trimmed_fastq_files')
    os.makedirs(work_dir,exist_ok=True)
    assert os.path.exists(work_dir)        # checks if the file path exists.
    comms = ("curl -L -o %s/sub.tar.gz https://ndownloader.figshare.com/files/14418248"%work_dir)       # tested, runs
    os.system(comms)
    #unzip files
    comms = ("tar xvf %s/sub.tar.gz -C %s"%(work_dir,work_dir))     # extract the downloaded files
    os.system(comms)
    trimmedFastqPath = os.path.join(work_dir,'sub')
    subprocess.run(['ls',trimmedFastqPath])
    return trimmedFastqPath

    # comms = "cd %s/dc_workshop_1 ; mkdir -p results/sam results/bam results/bcf results/vcf"%container_path     # to be tested in this function
    # os.system(comms) #.read()     
    # job.addChildJobFn(align_fastq)


# Aligning sample sequences to reference genome
def align_fastq(job,trimmedFastqPath,indexedFilesId,memory="2G", cores=2, disk="3G"):
    """_summary_: Created sam file using bwa tool. Returns fileID of sam File

    Args:
        job (_type_): _description_
        memory (str, optional): _description_. Defaults to "2G".
        cores (int, optional): _description_. Defaults to 2.
        disk (str, optional): _description_. Defaults to "3G".

    Returns:
        _type_: _description_
    """
    
    # generate index files directory
    work_dir = os.path.join(job.fileStore.getLocalTempDir(),'indexed_files')
    os.makedirs(work_dir,exist_ok=True)
    for filename in indexedFilesId:
        path_generator(job,file_id=indexedFilesId[filename],fileName=os.path.join(work_dir,filename))
   
    samFilePath = os.path.join(job.fileStore.getLocalTempDir(),'results_sam')
    os.makedirs(samFilePath,exist_ok=True)
    comms = f"singularity exec {container_path}/bwa_latest.sif bwa mem {work_dir}/'ecoli_rel606.fasta' {trimmedFastqPath}/SRR2584866_1.trim.sub.fastq {trimmedFastqPath}/SRR2584866_2.trim.sub.fastq > {samFilePath}/SRR2584866.aligned.sam"
    os.system(comms)
    samFileId = makeFileId(job,os.path.join(samFilePath,'SRR2584866.aligned.sam'))
    return samFileId

# Converting SAM file to BAM file using view option in samtools
# Sorting bam files
def convrt_sort(job,samFileId,memory="2G", cores=2, disk="3G"):
    
    samFilePath = path_generator(job,file_id=samFileId,fileName='SRR2584866.aligned.sam')
    subprocess.run(['echo',"##########################################"])
    subprocess.run(['cat',samFilePath])
    
    # creating temp directory to store bam file
    bamFilePath = os.path.join(job.fileStore.getLocalTempDir(),'results_bam')
    os.makedirs(bamFilePath,exist_ok=True)
    # convert sam file to bam
    comms = f"singularity exec {container_path}/samtools_latest.sif samtools view -S -b {samFilePath} > {bamFilePath}/SRR2584866.aligned.bam"
    os.system(comms)
    assert os.path.exists(bamFilePath+'/SRR2584866.aligned.bam')

    # sort bam file
    comms = f"singularity exec {container_path}/samtools_latest.sif samtools sort -o {bamFilePath}/SRR2584866.aligned.sorted.bam {bamFilePath}/SRR2584866.aligned.bam"
    os.system(comms)

    subprocess.run(['ls',bamFilePath])
    sortedBamId = makeFileId(job,filePath=os.path.join(bamFilePath,'SRR2584866.aligned.sorted.bam'))
    
    # save bam file to host storage.
    # bamPath=  path_generator(job,file_id=sortedBamId,fileName='/home/bioinfo/Desktop/Basesolve_jm/toil_filestore/global/SRR2584866.aligned.sorted.bam')
    # print('>>>>>>>>>>>>>>>>>>>>>>>BamPath',bamPath)
    # job.fileStore.deleteLocalFile(sortedBamId)
    return sortedBamId


# Variant calling
def variant_call(job,indexedFilesDict,sortedBamId,memory="2G", cores=2, disk="3G"):
    # generate index files directory
    indexFiles_dir = os.path.join(job.fileStore.getLocalTempDir(),'indexed_files')
    os.makedirs(indexFiles_dir,exist_ok=True)
    for filename in indexedFilesDict:
        path_generator(job,file_id=indexedFilesDict[filename],fileName=os.path.join(indexFiles_dir,filename))

    bamPath = path_generator(job,file_id=sortedBamId,fileName='SRR2584866.aligned.sorted.bam')
    baseDir = '/home/bioinfo/Desktop/Basesolve_jm/toil_filestore/global'
    comms = f"singularity exec {container_path}/bcftools_latest.sif bcftools mpileup -O b -o  {baseDir}/SRR2584866_raw.bcf -f {indexFiles_dir}/ecoli_rel606.fasta {bamPath}"
    os.system(comms)  #.read()      # works, tested.

    # Detecting Single Nucleotide Variants(SNV) from VCF file.
    # ploidy - number of chromosome sets in nucleus
    comms = f"singularity exec {container_path}/bcftools_latest.sif bcftools call --ploidy 1 -m -v -o {baseDir}/SRR2584866_variants.vcf {baseDir}/SRR2584866_raw.bcf"
    os.system(comms)     #.read()

    # Filter and report the SNV variants in variant calling format (VCF)
    comms = f"singularity exec {container_path}/bcftools_latest.sif vcfutils.pl varFilter {baseDir}/SRR2584866_variants.vcf  > {baseDir}/SRR2584866_final_variants.vcf"
    os.system(comms)     #.read()

    #generate final vcf fileId
    finalVcfId = makeFileId(job,filePath=baseDir+'/SRR2584866_final_variants.vcf')



if __name__ == "__main__":
    jobstore: str = tempfile.mkdtemp("tutorial_staging")
    os.rmdir(jobstore)
    options = Job.Runner.getDefaultOptions(jobstore)
    options.logLevel = "Debug"
    options.clean = "always"
    job = Job.wrapJobFn(main)

    with Toil(options) as toil:
        file_id = toil.start(job)
        toil.exportFile(file_id, "file://" + '/home/bioinfo/Desktop/Basesolve_jm/toil_filestore/global/test.sam')








