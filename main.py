import os
import tempfile

from toil.common import Toil
from toil.job import Job

# create a tem file
# write into it
# send its copy locally
# delete its temp copy


def run(job):
    
    #create a new directory with getLocalTempDir > returns path
    dir_path = job.fileStore.getLocalTempDir()
    child_dir = 'dc_workshop_1/data/ref_genome'
    path = dir_path+child_dir
    os.makedirs(path,exist_ok=True)

    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>job.tempDir',job.tempDir)

    #create a new temporary file with getLocalTempFile() and store its path
    temp_path = job.fileStore.getLocalTempFile(prefix='mytext',suffix='.txt')
    
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',temp_path)
    #write to this file
    with open(temp_path,'wb') as fi:
        fi.write(b'Yaayyy')
    
    # using writeGlobalFile() to understand its features
    inputFileID = job.fileStore.writeGlobalFile(localFileName='/home/bioinfo/Desktop/Basesolve_jm/toil_filestore/in.txt',cleanup=True)
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',inputFileID) # assuming that this function creates a temp copy of a local file.
    # inputFileID is FileID object, but it looks like a path to the file.

    #understanding function of readGlobalFile() 
    read_return = job.fileStore.readGlobalFile(fileStoreID=inputFileID,userPath='test.txt',
                                               cache=False,mutable='False',symlink=False) # so basically readGlobalFile, takes in file stored in filestore and stores its copy in temporary directory.
    # so what purpose is this function serving?   

    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',read_return)

    # lets try deleting files
    # job.fileStore.deleteLocalFile(inputFileID)
    # job.fileStore.deleteGlobalFile(inputFileID)

    with job.fileStore.readGlobalFileStream(inputFileID, encoding='utf-8') as fi:
            with job.fileStore.writeGlobalFileStream(encoding='utf-8') as (fo, outputFileID):
                fo.write(fi.read() + 'World!')                  # readGlobalTempFile can read only temp files and not local files
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',outputFileID)
    return outputFileID




if __name__ == "__main__":
    jobstore: str = tempfile.mkdtemp("tutorial_staging")
    os.rmdir(jobstore)
    options = Job.Runner.getDefaultOptions(jobstore)
    options.logLevel = "Debug"
    options.clean = "always"
    job = Job.wrapJobFn(run)

    with Toil(options) as toil:
        inputFileID_export=toil.start(job)
        toil.exportFile(inputFileID_export, "file://" + os.path.abspath(os.path.join('/home/bioinfo/Desktop/Basesolve_jm/toil_filestore/global', "out.txt")))
        
        
        
        # if not toil.options.restart:
        #     ioFileDirectory = os.path.join(os.path.dirname(os.path.abspath(__file__)), "stagingExampleFiles")
        #     inputFileID = toil.importFile("file://" + os.path.abspath(os.path.join(ioFileDirectory, "in.txt")))
        #     outputFileID = toil.start(HelloWorld(inputFileID))
        # else:
        #     outputFileID = toil.restart()

        # toil.exportFile(outputFileID, "file://" + os.path.abspath(os.path.join(ioFileDirectory, "out.txt")))