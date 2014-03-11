#! /bin/bash

export SOURCEADD="srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat"
export DIR="/store/user/hbakhshi/TauPlusX/V03-09-01_TauPlusX-Run2012A/312837b6d79fb3f6c6ab5915d0138724/"


#    export MainDirectoryName="/MultiJet/V03-09-01_MultiJet-Run2012D-PromptReco-v1/2c909eaf42bba0fead9d3c5a78eda5f3/"

#    export SOURCEADD="srm://t3se01.psi.ch:8443/srm/managerv2?SFN="
#    export DIR="pnfs/psi.ch/cms/trivcat/store/user/$USER/$MainDirectoryName"


if [ ! -f DSList ]; then
    lcg-ls "$SOURCEADD/$DIR" > DSList
    echo "DSList created!"
else
    echo "DSList exists!"
#    exit 1
fi

a=0
for DSName in $( cat DSList ); do
    export DSName=`basename $DSName`
#    let a='a%2'
#    if [ $a -eq 1 ]; then
    echo $DSName
#    python DBSInvalidateFile.py --DBSURL=https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet --lfn=$DIR$DSName
    lcg-del -l $SOURCEADD$DIR/$DSName
#    fi
    (( a++ ))
done

lcg-del -ld $SOURCEADD$DIR





