#!/bin/bash

# Usage: $0 folder1/ folder2 folder3 ...
set -e

for C in $*
do
	pushd ${C} > /dev/null 2>&1
	for i in *pstats_ts.*
	do
		B=$(echo ${i}|sed 's/^.*\.b\([0-9]\+\)\..*$/\1/')
		S=$(echo ${i}|sed -e 's/^.*\.\(set[0-9a-f][0-9a-f]\).txt$/\1/')
		if $(echo ${S}|grep pstats > /dev/null)
		then
				S='set00'
		fi
		CT=$(echo ${C}|sed -e 's:^.*/\([^/]\+$\):\1:' -e 's:/$::')

		TZERO=$(head -n 2 ${i}|tail -n1)
		TFINAL=$(tail -n1 ${i})
		echo -e "${CT}\t${B}\t${S}\t${TZERO}"
		echo -e "${CT}\t${B}\t${S}\t${TFINAL}"
	done
	popd > /dev/null 2>&1
done
