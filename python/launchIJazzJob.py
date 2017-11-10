#!/usr/bin/python

import argparse
import sys
import os
import subprocess

import python.batch as batch

if __name__ == "__main__":

    configDef = 'etc/config/ijazzEtaScale_RunII.conf'
    parser = argparse.ArgumentParser(description='ijazz fit')
    parser.add_argument('--ecalElf'     , type = str, default = None     , help = 'ecalElf file')
    parser.add_argument('-p', '--period', type = int, default = -1        , help = 'period')
    parser.add_argument('--config'      , type = str, default = configDef , help = 'config file' )
    parser.add_argument('--isMC'        , action='store_true', help = 'MC run' )
    parser.add_argument('--batch'       , action='store_true', help = 'batch' )
    parser.add_argument('-d', '--dry'   , action='store_true', help = 'batch' )

    args = parser.parse_args()

    if args.ecalElf is None :
       print ' IJazZ launcher requires an ecalElf file'
       sys.exit(1)



    fEcalElf = open( args.ecalElf, 'r')
    ecalElfLines = fEcalElf.readlines()


    periodMatch = {
        1 : { 'name': 'RunB'      , 'sel': ['s1 ','d1 ', 'd2 ']       },
        2 : { 'name': 'RunC'      , 'sel': ['s1 ','d3 ', 'd4 ', 'd5 ']  },
        3 : { 'name': 'RunD'      , 'sel': ['s1 ','d6 ']  },
        4 : { 'name': 'RunE'      , 'sel': ['s1 ','d7 ']       },
        5 : { 'name': 'RunF'      , 'sel': ['s1 ','d8 ']  },
        6 : { 'name': 'RunBCDEF'  , 'sel': ['s1 ','d1 ','d2 ', 'd3 ', 'd4 ', 'd5 ', 'd6 ', 'd7 ', 'd8 ']  },
    }
    
    if args.period > 0 :
        tmppath = 'etc/data/2017/tmp/'
        try: os.makedirs( tmppath )
        except: print 'Directory tmp already exists'
        fEcalElfName = tmppath + os.path.basename(args.ecalElf) + '.' + periodMatch[args.period]['name'] + '.tmp'
        fEcalElfOut = open( fEcalElfName, 'w' )
        for period in periodMatch[args.period]['sel']:
            for line in ecalElfLines:
                if 'selected' in line and period in line:
                    line = line.replace('d1 ','d').replace('d2 ','d').replace('d3 ','d').replace('d4 ','d')
                    line = line.replace('d5 ','d').replace('d6 ','d').replace('d7 ','d').replace('d8 ','d')
                    line = line.replace('s1 ','s').replace('s2 ','s').replace('s3 ','s').replace('s4 ','s')
                    fEcalElfOut.write(line)
                    fEcalElfOut.write('\n')

        print ' - Created file : ', fEcalElfName
        version = os.path.basename(args.ecalElf).replace('.dat','') + '_' +  periodMatch[args.period]['name']

        command  = 'source etc/scripts/setup.sh; '
        command += 'IJazZexe %s --ecalElf %s -v %s ' % (args.config, fEcalElfName, version )
        if args.isMC: command += ' --isMC'

        
        print '  --> I will submit/execute' 
        print command

        if args.dry: sys.exit(0)

        if args.batch :
            strMC = '.data'
            if args.isMC: strMC = '.mc'
            strEcalElf = batch.dirScriptBatch + os.path.basename(args.ecalElf).replace('.dat','')

            fscript = strEcalElf + '_' +  periodMatch[args.period]['name'] + strMC + '.sh'
            flog    = strEcalElf + '_' +  periodMatch[args.period]['name'] + strMC + '.log'
            print ' ** script: ', fscript
            print ' ** log   : ', flog
            batch.prepareGridCert()            
            batch.config(  command, fscript )
            print ' --> batch submit: ', fscript
            batch.bsubJob( fscript, flog )
            
        else:
            fileout = os.path.basename(args.ecalElf).replace('.dat','') + '_' +  periodMatch[args.period]['name'] + '.out'
            fileerr = os.path.basename(args.ecalElf).replace('.dat','') + '_' +  periodMatch[args.period]['name'] + '.err'
            print ' --> log   file: ', fileout
            print ' --> error file: ', fileerr
            subprocess.Popen(command, shell=True,
                                stdout=open(fileout,'w'),
                                stderr=open(fileerr,'w')) 
