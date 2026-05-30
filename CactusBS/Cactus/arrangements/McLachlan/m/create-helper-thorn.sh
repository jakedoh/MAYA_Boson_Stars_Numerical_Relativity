#! /bin/bash

# Abort on errors
set -e
set -o pipefail
set -u

orig="$1"
name="$2"

# Create a thorn $name_Helper

rm -rf "${name}_Helper"
cp -R "prototype/${orig}_Helper" "${name}_Helper"
find "${name}_Helper" -name '*~' | xargs rm -f
find "${name}_Helper" -type f | xargs perl -pi -e "s/${orig}/${name}/g"

true
