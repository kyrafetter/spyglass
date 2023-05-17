#!/bin/bash

# This script contains command line tests for mypileup

die()
{
    BASE=$(basename "$0")
    echo "$BASE error: $1" >&2
    exit 1
}

runcmd_pass()
{
    echo "[runcmd_pass]: $1"
    sh -c "$1" >/dev/null 2>&1 || die "Error running: $1"
}

runcmd_fail()
{
    echo "[runcmd_fail]: $1"
    sh -c "$1" >/dev/null 2>&1 && die "Command should have failed: $1"
}

if [ $# -eq 0 ]; then
    # use default example location
    EXDATADIR="example-files"
elif (( $# != 1 )) ; then
    echo "usage: cmdline_tests.sh {example_dir}" 2>&1
    echo "Expected 1 arguments but recieved $#" 2>&1
    exit 1
else
    EXDATADIR=$1
fi

TMPDIR=$(mktemp -d -t tmp-XXXXXXXXXX)

echo "Saving tmp files in ${TMPDIR}"

runcmd_pass "mypileup -f ${EXDATADIR}/test.fa ${EXDATADIR}/test.bam"

runcmd_fail "mypileup -f ${EXDATADIR}/testXYZ.fa ${EXDATADIR}/test.bam"
