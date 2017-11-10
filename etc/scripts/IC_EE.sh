#!/bin/bash


#config=etc/config/ijazzIC_Xtal_FitSample.conf
config=etc/config/ijazzIC_Xtal_MC.conf
userDir=IC_perXtal
version=RegrV5
runRangesFile=etc/data/runRangesAll.txt
noFit=0
isMC=1

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
    
### send the jobs in batch for each run ranges
if [ $noFit -eq 0 ]; then
for runrange in $( cat $runRangesFile | awk '{print $1}' ); do 
    rm -f IJazZ_output/FitResults/${userDir}/$runrange/IJazZ_E[BE]_*_${version}.*
    config=$configTmp
    
    cd $dirscript
    script=script$$
    script=${script}_EEmIC_Run${runrange}_data.sh
    jobConfig $script $isMC $version $config $runrange 2    
    bsub $optionsub $script

    script=script$$
    script=${script}_EEpIC_Run${runrange}_data.sh
    jobConfig $script $isMC $version $config $runrange 3    
    bsub $optionsub $script

    
#    script=script$$
#    script=${script}_EscalesHgg_Run${runrange}_mc.sh
#    jobConfig $script 1 $version $config $runrange
#    bsub $optionsub $script    
    cd -

#    sleep 30
done
fi

