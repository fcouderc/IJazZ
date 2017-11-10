
import os
import shutil
import subprocess


dirAFS = os.getcwd()
dirScriptBatch  = '%s/tmp/' % dirAFS


bsuboptions    = "-q 2nd"


def prepareGridCert():
    if not os.path.exists(dirScriptBatch):
        print " mkdir ",dirScriptBatch
        os.mkdir( dirScriptBatch )
    
#    u = commands.getoutput('voms-proxy-init')
    userid = os.environ['KRB5CCNAME'].split('_')[1]
    pathToGridCert='/tmp/x509up_u%s' % userid
    gridCert = os.path.basename(pathToGridCert)
    shutil.copy(pathToGridCert,dirScriptBatch)
    pathToGridCert = '%s/%s' % (dirScriptBatch,gridCert)
    return pathToGridCert



def config( commandexe, scriptName ):
    script = open(scriptName,'w')
    script.write('#!/bin/bash\n')
    script.write('cd %s\n'% dirAFS) 
    script.write('source etc/scripts/setup.sh\n')
    pathToGridCert = prepareGridCert()
    script.write('export X509_USER_PROXY=%s\n'%pathToGridCert)
    script.write('%s\n'%commandexe)
    script.close()
    st=os.stat(scriptName)
    os.chmod(scriptName, st.st_mode | 0o111)
 #   print os.stat(scriptName)


def bsubJob( script, filenameJobList ):
    fileJobList = open( filenameJobList, 'a' )
    os.chdir(dirScriptBatch)
    proc = subprocess.Popen( 'bsub %s %s' % ( bsuboptions,script) , shell = True, stdout=subprocess.PIPE)
    tmp = proc.stdout.read().split(' ')
    jobid = tmp[1].split('<')[1].split('>')[0]
    fileJobList.write( '%s  %s\n' % (jobid,script) )
    fileJobList.close()
    os.chdir( dirAFS )








