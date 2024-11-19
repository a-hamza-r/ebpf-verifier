#!/bin/bash

### Update this path to the root of your prevail repository.
PREVAIL_ROOT=$(pwd)


EBPF_BENCHMARKS=${PREVAIL_ROOT}/ebpf-samples
PREVAIL_CHECK=${PREVAIL_ROOT}/check
DOMAINS="inference"
PREFIX=prevail_$(date +"%m%d%y%H%M")

for dom in $DOMAINS
do
    rm -f log_${dom}.txt
    echo -n "Running Prevail with $dom ... "
    echo "File,Result,Cpu,Mem"  1>> ${PREFIX}_${dom}.csv
    for f in ${EBPF_BENCHMARKS}/*/*.o
    do
      ${PREVAIL_CHECK} $f -l 2> /dev/null | while read -r s; do
      if [ -n "$s" ]; then
        values=$(echo "$s" | awk -F' |=' '{print $2, $4}')
        read section function <<< "$values"
        echo "${PREVAIL_CHECK} ${f} ${section} ${function} --domain=${dom}" >> log_${dom}.txt
        echo -n $f:$section:$function 1>> ${PREFIX}_${dom}.csv
        o=$(${PREVAIL_CHECK} ${f} ${section} ${function} --domain=${dom} 2>>log_${dom}.txt)
        echo -n ",$o" 1>> ${PREFIX}_${dom}.csv
        echo 1>> ${PREFIX}_${dom}.csv
      fi
		done
    done
    echo "DONE"
done
