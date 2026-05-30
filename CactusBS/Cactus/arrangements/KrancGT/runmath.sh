#! /bin/bash

# Abort on errors
set -e
# null-expand non-matching globs
set -u

MATHEMATICA="math"

script=$1

if test -z "$script"; then
    echo "Usage:"
    echo "$0 <script.m>"
    exit 2
fi

KRANCPATH="${KRANCPATH:+$KRANCPATH}"
if [ -z "$KRANCPATH" ] ; then
    KRANCPATHS=(. ../../../../arrangements/../repos ../../../arrangements/../repos $HOME)
    # look for Kranc in a number of likely locations
    # choice 2 and 3 are $CCTK_HOME if the thorn resides in arrangements or
    # git-repos
    for i in ${KRANCPATHS[@]} ; do
        if [ -d $i/kranc/Tools ] ; then
            KRANCPATH=$i/kranc
            break
        elif [ -d $i/Kranc/Tools ] ; then
            KRANCPATH=$i/Kranc
            break
        fi
    done
fi
if [ -z "$KRANCPATH" ] ; then # try some more clever things...
    KRANCBIN="$(which kranc || true)"
    if [ -h "$KRANCBIN" ] ; then
        KRANCBIN=$(readlink "$KRANCBIN")
    fi
    if [ -d "$(dirname "$KRANCBIN")/../Tools" ] ; then
        KRANCPATH="$(dirname "$KRANCBIN")/.."
    fi
fi
if [ -n "$KRANCPATH" ] ; then
    LOCALKRANCPATH="$KRANCPATH/Tools/CodeGen:$KRANCPATH/Tools/MathematicaMisc"
    
    if [ -z "${MATHPATH:+}" ]
    then
        MATHPATH=$LOCALKRANCPATH
    else
        MATHPATH=$LOCALKRANCPATH:$MATHPATH
    fi
    
    export MATHPATH
else
    echo -e "Could not find Kranc. Expect trouble unless you have made sure Mathematica\nwill find Kranc. You might also consider setting \$KRANCPATH."
fi

echo "MATHPATH = $MATHPATH"

error=$(basename $script .m).err
output=$(basename $script .m).out

rm -f $output

# Run Mathematica to regenerate the code
< $script "$MATHEMATICA" | tee $error

if grep 'KrancError' $error; then
    echo
    echo "There was an error when running Kranc on $script."
    echo "The file $error contains details."
    echo
    echo "*** The Cactus thorns have NOT been updated. ***"
    echo
    exit 1
fi

mv $error $output
