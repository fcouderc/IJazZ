#!/bin/bash



config=etc/config/ijazzResolutionHgg.conf
userDir=FitResolutionHggTail

optionsub="-q 1nw "
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
    exe="./bin/IJazZexe -a 0 $optMC -v $vers -u $userDir $conf --numCPU=12"
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

version=ResoHgg


cd $dirscript

script=script$$
script=${script}_ResoHgg_mc.sh
jobConfig $script 1 $version  $config
bsub $optionsub $script    

script=script$$
script=${script}_ResoHgg_data.sh
jobConfig $script 0 $version  $config
bsub $optionsub $script    

cd -