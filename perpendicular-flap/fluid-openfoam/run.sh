#!/usr/bin/env bash
set -e -u

touch probedLocations.dat
rm probedLocations.dat

touch preciceForceWrite.dat
rm preciceForceWrite.dat

touch preciceDisplacementRead.dat
rm preciceDisplacementRead.dat

. ../../tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

blockMesh

../../tools/run-openfoam.sh "$@"
. ../../tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
