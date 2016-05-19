#!/bin/bash
if [ $# -ne 4 ]
then
	echo "Usage: $0 <Tool> <Case File> <Comparison Directory> <Comparison Case Number>"
	exit 1
fi

TOOL=$1
CASE=$2
BASEDIR=$3
CASENUM=$4


TESTOUTPATH=${BASEDIR}/Case\ ${CASENUM}
BASENAME=$(grep -E "^basename" ${CASE}|sed -e 's/^[^=]*[= ]*//')
BASENAME=$(basename "${BASENAME}")
CASEOUTPATH=$(grep -E "^output" ${CASE}|sed -e 's/^[^=]*[= ]*//')
PARTICLEFILE=$(grep -E "^particles" ${CASE}|sed -e 's/^[^=]*[= ]*//')

echo "Running ${TOOL} for case ${CASE}..."

mkdir -p "${TESTOUTPATH}"
./${TOOL} ${CASE}

echo "Copying files to ${TESTOUTPATH}..."
cp "$CASEOUTPATH/${BASENAME}.tracks.txt" "${TESTOUTPATH}"
cp "$CASEOUTPATH/${BASENAME}".*.vtu "${TESTOUTPATH}"
cp "${CASE}" "${TESTOUTPATH}"
cp "${PARTICLEFILE}" "${TESTOUTPATH}"

echo Plotting to "${TESTOUTPATH}/${BASENAME}.png..."
./plot_tools/plot_combined.py -nm -t ${CASE} -o "${TESTOUTPATH}/${BASENAME}.png"

echo "Test run complete. Output files in ${TESTOUTPATH}"
