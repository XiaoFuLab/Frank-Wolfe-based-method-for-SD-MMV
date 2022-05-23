#!/bin/bash
set -e
methodName=MERIT
main() {
    {
        /usr/bin/time -v matlab -nodisplay -nodesktop -r "method_name='$methodName'; N=$N; run_it; quit" &
    } 2>&1 | tee  ./logging_"$methodName"_"$N".txt &
}

for N in 200 1000
do
    main
done

