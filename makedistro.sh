#!/bin/bash
# This script tries to create a new Penlab distribution by deleting
# unnecessary files from the default checkout. It still should be
# checked that there is nothing strange left.
# Usage:
#   ./makedistro.sh [dirname]
# Should be run from the checkout directory.
# If dirname is missing, default directory is DISTRIBUTION.

# default directory where to create a distribution
distrodir=PENLABv101a

# choose directory
if [ -n "$1" ]; then
  distrodir="$1"
fi

echo Disribution directory is \""$distrodir"\"
if [ -d "$distrodir" ]; then
  echo ERROR: distribution directory already exists, will not overwrite!
  exit 1
fi

mkdir "$distrodir"
echo " "

# copy (almost) all

echo Copying all items to the new location
find . -maxdepth 1 -mindepth 1 ! -name "$distrodir" -a ! -name attic -a ! -name "DISTR*" -print
find . -maxdepth 1 -mindepth 1 ! -name "$distrodir" -a ! -name attic -a ! -name "DISTR*" -exec cp -r {} "$distrodir" \;
echo " "

# deleting all unnecessary stuff

echo Deleting all .svn directories
find "$distrodir" -type d -name .svn -print
find "$distrodir" -type d -name .svn -exec rm -rf {} \;
echo " "

echo Deleting all forgotten log files
find "$distrodir" \( -name 'penm_log*.txt*' -o -name asl_error.log -o -name 'logtest*.txt' -o -name 'diary' \) -exec rm -v {} \;
echo " "

echo Deleting all temp files \(vim: backup and swap, matlab: asv\)
find "$distrodir" \( -name '*.swp' -o -name '*~' -o -name '*.asv' \) -exec rm -v {} \;
echo " "

echo Deleting all DCF files
find "$distrodir/datafiles" -name '*.dcf' -exec rm -v {} \;
echo " "

echo Deleting all bigger NL files
find "$distrodir/datafiles" -name '*.nl' -size +34k -exec rm -v {} \;
echo " "

echo Deleting all bigger GEO files
find "$distrodir/applications" \( -name '*.geo' -size +8k -a ! -name 'tb1.geo' -o -name 'tenbar.geo' -o -name 'wheel.geo' \) -exec rm -v {} \;
echo " "

echo Deleting all tex files except the PDF
find "$distrodir/tex" -type f ! -name '*.pdf' -exec rm -v {} \;
find "$distrodir/tex" -mindepth 2 -type d -exec rm -v -rf {} \;
echo " "

echo Deleting files explicitly listed
rm -v "$distrodir/makedistro.sh"
rm -v "$distrodir/doc/internal"
rm -v "$distrodir/doc/manual.txt"
rm -v "$distrodir/doc/todo.txt"     # perhaps leave it?
rm -v "$distrodir/penlab_flyer_arial.docx"
echo " "

# anything to remove from source?
# anything to remove from utilities?

# check if the version in the manual is the same as in the solver
vmanual=`grep 'PenLab, version' "$distrodir/README.txt" | sed 's/, version//'`
vsrc=`grep 'solvername *=' "$distrodir/source/@penlab/penlab.m" | sed "s/[^']*'\(.*\)'.*/\1/"`
if [ "$vmanual" == "$vsrc" ]; then
  echo Current version is: "$vsrc"
else
  echo !!!!!!!
  echo WARNING, version in manual is "$vmanual" and in source "$vsrc" ! 
  echo !!!!!!!
fi
echo " "

echo Make sure that the html doc generated from manual.m is up to date!
echo " "
echo Done.

