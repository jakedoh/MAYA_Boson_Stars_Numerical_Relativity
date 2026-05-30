#! /bin/bash

# Abort on errors
set -e
set -o pipefail
set -u

script="${1:-}"

if [[ -z $script ]]; then
    echo "Usage:"
    echo "$0 <script.m>"
    exit 2
fi

outfile=${2:-}

if [[ -z $outfile ]]; then
    outfile=$(basename "$script" .m)
fi

error="$outfile.err"
output="$outfile.out"

rm -f "$output"

# Run Kranc to regenerate the code
../../../repos/Kranc/Bin/kranc "$script" | tee "$error"
(( PIPESTATUS )) && exit "$PIPESTATUS"

mv "$error" "$output"

true
