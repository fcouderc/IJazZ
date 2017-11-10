#!/bin/bash


config=etc/config/ijazzResolutionHgg.conf
userDir=FitEscalesHgg-phoRegr
runRangesFile=etc/data/runranges_Moriond2013.dat
etaScaleCorr=undefined
version=EscaleHggCoarse

noFit=1

optionsub="-q 1nd "
dirafs=`pwd`
dirscript=tmpBatchOut/

mkdir -p $dirscript

jobConfig() {
    script=$1
    mc=$2
    vers=$3
    conf=$4
    runrange=$5
    optMC=""
    if [ $mc -eq 1 ]; then
	optMC="--isMC"
    fi
    exe="./bin/IJazZexe -a 0 $optMC -v $vers -u $userDir $conf --runRange $runrange --eleEnCorr 2"
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

configTmp=$dirscript/ijazzEscalesHgg_v0.conf
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
    script=${script}_EscalesHgg_Run${runrange}_data.sh
    jobConfig $script 0 $version $config $runrange
    
    bsub $optionsub $script
    
    script=script$$
    script=${script}_EscalesHgg_Run${runrange}_mc.sh
    jobConfig $script 1 $version $config $runrange
    bsub $optionsub $script    
    cd -

    sleep 30
done
fi

### compute the eta scale for each run ranges
for runrange in $( cat $runRangesFile | awk '{print $1}' ); do 
    respDataEB=IJazZ_output/FitResults/${userDir}/$runrange/IJazZ_EB_Data_${version}.root.fittedResp
    respMCEB=IJazZ_output/FitResults/${userDir}/$runrange/IJazZ_EB_MC_${version}.root.fittedResp
    respDataEE=IJazZ_output/FitResults/${userDir}/$runrange/IJazZ_EE_Data_${version}.root.fittedResp
    respMCEE=IJazZ_output/FitResults/${userDir}/$runrange/IJazZ_EE_MC_${version}.root.fittedResp
    
    while [ ! -f $respDataEB ]; do 
	echo "Waiting file: $respDataEB to be done..."
	sleep 200; 
    done
    while [ ! -f $respMCEB ]; do 
	echo "Waiting file: $respMCEB to be done..."
	sleep 200; 
    done
    while [ ! -f $respDataEE ]; do 
	echo "Waiting file: $respDataEE to be done..."
	sleep 200; 
    done
    while [ ! -f $respMCEE ]; do 
	echo "Waiting file: $respMCEE to be done..."
	sleep 200; 
    done

    ./bin/IJazZetaScale -v $version -u $userDir -r $runrange
done
    
echo ""
echo "######################################################"
echo "########## For E scale correction, add to config file:"
find ../IJazZ/IJazZ_output/FitResults/${userDir} -name "*${version}*.etaScaleDataOverMC" | awk '{print "ecalRespCorr = "$1}'
