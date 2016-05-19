#!/bin/bash
set -e
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

if [ $# -ne 2 ]; then
	cat <<-EOF
	Usage: $(basename $0) <searchbase> <VTU>
	       searchbase  = Base folder for PVD search
	       VTU = VTU file to use as background
	EOF
	exit 1
fi

TMPNAME=$(mktemp)
exec 3<>${TMPNAME}
DISPLAY=:0 Xephyr -screen 1200x800 -glamor -no-host-grab -title Xephyr -displayfd 3 &
jobs
sleep 5
read XEPHDISPLAY < ${TMPNAME}
echo ${XEPHDISPLAY}
export DISPLAY=:${XEPHDISPLAY}
rm ${TMPNAME}
find $1 -name '*-tracks.pvd'|xargs ${DIR}/../python_tools/render_tracks_pv.py $2
kill %1
