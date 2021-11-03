#!/bin/bash
set -e

# # SET 1
# M=50
# N=40
# L=2000
# numTrials=20
# methodName=FW
# code=SNR

# SET 2
M=50
N=40
SNR=10
numTrials=1
methodName=FW
code=memory_lambda_0


main() {
    expCodeName="$code"_"$L"
    pathDir=./results/$expCodeName/$methodName/
    mkdir -p $pathDir

    # {
    #     /usr/bin/time -v matlab -nodisplay -nodesktop \
    #     -r "methodName='$methodName'; expCodeName='$expCodeName'; M=$M; N=$N; L=$L; SNR=$SNR; numTrials=$numTrials; isCalibration=1; expSetup__script; quit" &
    # } 2>&1 | tee  "$pathDir"/logging_calibration.txt &

    {
        /usr/bin/time -v matlab -nodisplay -nodesktop \
        -r "methodName='$methodName'; expCodeName='$expCodeName'; M=$M; N=$N; L=$L; SNR=$SNR; numTrials=$numTrials; isCalibration=0; expSetup__script; quit" 
    } 2>&1 | tee  "$pathDir"/logging.txt 
}

for L in 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 
do
    main
done


