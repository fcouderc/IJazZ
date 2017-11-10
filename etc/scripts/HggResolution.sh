#!/bin/bash



config=etc/config/ijazzResolutionHgg.conf
userDir=FitResolutionHgg

etaScaleCorr=IJazZ_output/FitResults/ecalEtaScale/0-999999/IJazZ_EB_Data_ICcombTCcomb_v3.root.fittedResp.etaScaleDataOverMC

optionsub="-q 1nd "
dirafs=`pwd`
dirscript=tmpBatchOut/

mkdir -p $dirscript

jobConfig() {
    script=$1
    mc=$2
    vers=$3
    conf=$4
    optMC=""
    if [ $mc -eq 1 ]; then
	optMC="--isMC"
    fi
    exe="./bin/IJazZexe -a 0 $optMC -v $vers -u $userDir $conf"
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

configTmp=$dirscript/ijassResolutionHgg_Iter0.conf
cp -f $config $configTmp
if [ -f $etaScaleCorr ]; then
    echo "ecalRespCorr = "$etaScaleCorr >> $configTmp
fi
    

for it in $( seq 1 1 3 ); do 
    version=ResoHgg_Iter${it}
    config=$configTmp
    configTmp=$dirscript/ijassResolutionHgg_Iter${it}.conf
    rm -f $configTmp
    cp $config $configTmp

    if [ $it -eq 1 ]; then
	resoData=IJazZ_output/FitResults/${userDir}/0-999999/IJazZ_EB_Data_${version}.root.fittedReso
	resoDataEE=IJazZ_output/FitResults/${userDir}/0-999999/IJazZ_EE_Data_${version}.root.fittedReso
	rm -f IJazZ_output/FitResults/${userDir}/0-999999/IJazZ_E[BE]_Data_${version}.root.fittedReso 
    fi
    
    resoMC=IJazZ_output/FitResults/${userDir}/0-999999/IJazZ_EB_MC_${version}.root.fittedReso
    rm -f IJazZ_output/FitResults/${userDir}/0-999999/IJazZ_E[BE]_MC_${version}.root.fittedReso
    
    cd $dirscript
    script=script$$
    script=${script}_ResoHggIter${it}_mc.sh
    jobConfig $script 1 ResoHgg_Iter${it} $configTmp
    bsub $optionsub $script    
    if [ $it -eq 1 ]; then
	script=script$$
	script=${script}_ResoHggIter${it}_data.sh
	jobConfig $script 0 $version $configTmp
	bsub $optionsub $script 
    fi
    cd -    

    resoMCEE=IJazZ_output/FitResults/${userDir}/0-999999/IJazZ_EE_MC_${version}.root.fittedReso
 
    ### wait for batch jobs to be done before next iteration
    while [ ! -f $resoData ]; do 
	echo "Waiting file: $resoData to be done..."
	sleep 200; 
    done
    while [ ! -f $resoDataEE ]; do 
	echo "Waiting file: $resoDataEE to be done..."
	sleep 200; 
    done
    while [ ! -f $resoMC ]; do 
	echo "Waiting file: $resoMC to be done..."
	sleep 200; 
    done
    while [ ! -f $resoMCEE ]; do 
	echo "Waiting file: $resoMCEE to be done..."
	sleep 200; 
    done

    if [ $it -eq 1 ]; then
	echo "resoData = "$resoData >> $configTmp
    fi
    echo "resoMC   = "$resoMC   >> $configTmp

done

