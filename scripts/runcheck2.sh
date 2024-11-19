#!/bin/bash

### Update this path to the root of your prevail repository.
PREVAIL_ROOT=$(pwd)


EBPF_BENCHMARKS=${PREVAIL_ROOT}/ebpf-samples
PREVAIL_CHECK=${PREVAIL_ROOT}/check
DOMAINS=$1
PREFIX=prevail_$(date +"%m%d%y%H%M")

for dom in $DOMAINS
do
    rm -f ${PREFIX}_${dom}.txt
    echo -n "Running Prevail with $dom ... "
    for f in ${EBPF_BENCHMARKS}/*/*.o
    do
		${PREVAIL_CHECK} $f -l 2> /dev/null | while read -r s; do
    	  if [ -n "$s" ]; then
    	  values=$(echo "$s" | awk -F' |=' '{print $2, $4}')
    	  read section function <<< "$values"
    	  echo "${PREVAIL_CHECK} ${f} ${section} ${function} --domain=${dom} -v" >> ${PREFIX}_${dom}.txt
		  # echo -n $f:$section:$function 1>> ${PREFIX}_${dom}.txt
	      echo -n ""
		  o=$(${PREVAIL_CHECK} ${f}  ${section} ${function} --domain=${dom} -v >>${PREFIX}_${dom}.txt 2>&1)
		  echo 1>> ${PREFIX}_${dom}.txt
		  fi
		done
	done
    echo "DONE"
done
