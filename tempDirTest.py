from toil.common import Toil
from toil.job import Job
import os
import tempfile 
import requests
import shlex
import subprocess, shlex




container_path='/home/bioinfo/singularity_containers'

def path_generator(job,file_id):
    file_path = job.fileStore.readGlobalFile(fileStoreID=file_id,userPath='test.txt')
    return file_path


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
    alignedFastqJob = trimmedFastqJob.addChildJobFn(align_fastq)

    return alignedFastqJob.rv()

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
    test_path = job.wrapJobFn(path_generator,fasta_id).rv()
    print('>>>>>>>>>>>>>>>>>>test_path',test_path)

    indexedFa_path = job.fileStore.readGlobalFile(fasta_id, os.path.join(job.tempDir, 'indexedFasta.fasta'))        # creating a local copy from filestore using fileID
    comms = f"singularity exec {container_path}/bwa_latest.sif bwa index {indexedFa_path}"
    args= shlex.split(comms)
    subprocess.run(args)

    # writing a file into filestore
    indexed_id = job.fileStore.writeGlobalFile(localFileName=indexedFa_path)
    return indexed_id


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
    
    # child_dir = 'dc_workshop_1/data/trimmed_fastq_small'            
    
    # path = os.path.join(tempdir_path,child_dir)
    
    work_dir = os.path.join(job.tempDir,'trimmed_fastq_files')

    os.makedirs(work_dir,exist_ok=True)
    assert os.path.exists(work_dir)        # checks if the file path exists.

    comms = ("curl -L -o %s/sub.tar.gz https://ndownloader.figshare.com/files/14418248"%work_dir)       # tested, runs
    os.system(comms) #.read()

    comms = ("tar xvf %s/sub.tar.gz -C %s"%(work_dir,work_dir))     # extract the downloaded files
    os.system(comms)

    # comms = "cd %s/dc_workshop_1 ; mkdir -p results/sam results/bam results/bcf results/vcf"%container_path     # to be tested in this function
    # os.system(comms) #.read()     
    # job.addChildJobFn(align_fastq)


# Aligning sample sequences to reference genome
def align_fastq(job,memory="2G", cores=2, disk="3G"):
    """_summary_: Created sam file using bwa tool. Returns fileID of sam File

    Args:
        job (_type_): _description_
        memory (str, optional): _description_. Defaults to "2G".
        cores (int, optional): _description_. Defaults to 2.
        disk (str, optional): _description_. Defaults to "3G".

    Returns:
        _type_: _description_
    """
    
    work_dir = job.tempDir
    sam_file_dir = os.path.join(work_dir,'results','sam')
    os.makedirs(sam_file_dir,exist_ok=True)

    comms = f"singularity exec {container_path}/bwa_latest.sif bwa mem {work_dir}/ecoli_rel606.fna {work_dir}/trimmed_fastq_files/sub/SRR2584866_1.trim.sub.fastq {work_dir}/trimmed_fastq_files/sub/SRR2584866_2.trim.sub.fastq > {sam_file_dir}/SRR2584866.aligned.sam"
    os.system(comms)

    # assert os.path.exists(f'{work_dir}/ecoli_rel606.fna')
    # assert os.path.exists(f'{work_dir}/trimmed_fastq_files/sub/SRR2584866_1.trim.sub.fastq')
    # assert os.path.exists(f'{sam_file_dir}/SRR2584866.aligned.sam')

    # writing a file into filestore
    fpath = os.path.join(sam_file_dir,'SRR2584866.aligned.sam')
    samFile_id = job.fileStore.writeGlobalFile(localFileName=fpath)
    return samFile_id

    # process.read()
    # job.addChildJobFn(convrt_sort,path=path,tempdir_path=tempdir_path)



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








