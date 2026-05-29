#! /bin/bash

################################################################################
# Prepare
################################################################################

# Set up shell
if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    set -x                      # Output commands
fi
set -e                          # Abort on errors



################################################################################
# Configure Cactus
################################################################################

if [ -z "${OPENCL_DIR}" ]; then
    echo "BEGIN MESSAGE"
    echo "OpenCL selected, but OPENCL_DIR not set. Checking some places..."
    echo "END MESSAGE"

    DIRS="/usr /usr/local /opt/local /usr/local/packages usr/local/apps /opt/local ${HOME} c:/packages/"
    # look into each directory
    for dir in $DIRS; do
        # libraries might have different file extensions
        for libext in a dll dll.a dylib so; do
            # libraries can be in /lib or /lib64
            for libdir in lib64 lib; do
                FILES="$libdir/libOpenCL.${libext}"
                # assume this is the one and check all needed files
                OPENCL_DIR="$dir"
                for file in $FILES; do
                    # discard this directory if one file was not found
                    if [ ! -r "$dir/$file" ]; then
                        unset OPENCL_DIR
                        break
                    fi
                done
                # don't look further if all files have been found
                if [ -n "$OPENCL_DIR" ]; then
                    break
                fi
            done
           # don't look further if all files have been found
            if [ -n "$OPENCL_DIR" ]; then
                break
            fi
        done
        # don't look further if all files have been found
        if [ -n "$OPENCL_DIR" ]; then
            break
        fi
    done
    
    # Prefer the system OpenCL on Mac OSX
    if [ -r /System/Library/Frameworks/OpenCL.framework ]; then
        OPENCL_DIR=/System/Library/Frameworks/OpenCL.framework
        OPENCL_INC_DIRS=/System/Library/Frameworks/OpenCL.framework/Headers
        OPENCL_LIB_DIRS=/System/Library/Frameworks/OpenCL.framework/Libraries
        OPENCL_LIBS="-Wl,-framework -Wl,OpenCL"
    fi
fi



################################################################################
# Configure Cactus
################################################################################

# Set options
: ${OPENCL_INC_DIRS:=${OPENCL_DIR}/include}
: ${OPENCL_LIB_DIRS:=${OPENCL_DIR}/lib}
: ${OPENCL_LIBS:=OpenCL}

OPENCL_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh ${OPENCL_INC_DIRS})"
OPENCL_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh ${OPENCL_LIB_DIRS})"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "OPENCL_DIR      = ${OPENCL_DIR}"
echo "OPENCL_INC_DIRS = ${OPENCL_INC_DIRS}"
echo "OPENCL_LIB_DIRS = ${OPENCL_LIB_DIRS}"
echo "OPENCL_LIBS     = ${OPENCL_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(OPENCL_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(OPENCL_LIB_DIRS)'
echo 'LIBRARY           $(OPENCL_LIBS)'
