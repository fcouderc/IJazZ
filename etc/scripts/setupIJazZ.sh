cname=`uname -n | head -c6`
if [ "$cname" == "lxplus" ]; then
    echo " --- IJazZ for LXPLUS --- "
    DIRSLC6=/afs/cern.ch/cms/slc6_amd64_gcc491/
    LCG=/afs/cern.ch/sw/lcg/
    
    ROOTSYS=$LCG/app/releases/ROOT/6.02.10/x86_64-slc6-gcc49-opt/root/
    XRDCP=$LCG/external/xrootd/3.1.0p2/x86_64-slc6-gcc49-opt/
    . $LCG/external/gcc/4.9.1/x86_64-slc6/setup.sh
    . $ROOTSYS/bin/thisroot.sh
    #LD_LIBRARY_PATH=`pwd`/lib:$XRDCP/lib64:$LD_LIBRARY_PATH:.
    LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH:.
    LD_LIBRARY_PATH=$DIRSLC6/external/boost/1.47.0/lib/:$LD_LIBRARY_PATH
else
    echo " --- IJazZ local --- "
fi

touch ./python/__init__.py

IJAZZDIR=$PWD
export LD_LIBRARY_PATH=$IJAZZDIR/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DYLD_LIBRARY_PATH
export PATH=$IJAZZDIR/bin:$PATH

