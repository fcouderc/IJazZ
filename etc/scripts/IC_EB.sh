#!/bin/bash


config=etc/config/ijazzIC10x10_Xtal_FitSample.conf
userDir=ICEB5x5
etaScaleCorr=undefined
version=RegrV5_EB
runRangesFile=etc/data/runRangesAll.txt
noFit=0

optionsub="-q 2nw "
dirafs=`pwd`
dirscript=tmpBatchOut/

mkdir -p $dirscript

jobConfig() {
    script=$1
    mc=$2
    vers=$3
    conf=$4
    runrange=$5
    EEpart=$6
    optMC=""
    if [ $mc -eq 1 ]; then
	optMC="--isMC"
    fi
    exe="./bin/IJazZexe -a 0 $optMC -v $vers -u $userDir $conf --runRange $runrange  --noEB=$6 --subSampleOnly=1 -d 2"
    echo $exe
    cat > $script<<EOF
#!/bin/bash
cd $dirafs 
source etc/scripts/setupIJazZ.sh
echo "ROOTSYS: \$ROOTSYS"
gcccom="\`which gcc\`"
echo "gcc:" \$gcccom
echo "where am I:\`pwd\`"

echo  $exe
$exe
EOF
    chmod +x $script
}

configTmp=$dirscript/ijazzEE_IC_v0.conf
cp -f $config $configTmp
if [ -f $etaScaleCorr ]; then
    echo "ecalRespCorr = "$etaScaleCorr >> $configTmp
fi
    
### send the jobs in batch for each run ranges
if [ $noFit -eq 0 ]; then
for runrange in $( cat $runRangesFile | awk '{print $1}' ); do 
    rm -f IJazZ_output/FitResults/${userDir}/$runrange/IJazZ_E[BE]_*_${version}.*
    config=$configTmp
    
    cd $dirscript
    script=script$$
    script=${script}_EBIC_Run${runrange}_data.sh
    jobConfig $script 0 $version $config $runrange -1    
    bsub $optionsub $script


    
#    script=script$$
#    script=${script}_EscalesHgg_Run${runrange}_mc.sh
#    jobConfig $script 1 $version $config $runrange
#    bsub $optionsub $script    
    cd -

#    sleep 30
done
fi

