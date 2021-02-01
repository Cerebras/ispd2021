#!/bin/bash

EXEC=$1; shift;
ISPD_DIR=$1; shift;
TEST_SUITE=$1; shift;

if [ "all" == ${TEST_SUITE} ]
then
    for i in ${ISPD_DIR}/tests/*/; do
        for j in ${i}*.in; do
            echo ${j}
            ${EXEC} "${j}" $@ < "${j%.*}.out"
            echo ""
        done
    done
else
    TEST_CASE=$1; shift;
    if [[ "all" == ${TEST_CASE} || "" == ${TEST_CASE} ]]
    then
        for i in ${ISPD_DIR}/tests/${TEST_SUITE}/*.in; do
            echo ${i}
            ${EXEC} "${i}" $@ < "${i%.*}.out"
            echo ""
        done
    else
        echo ${ISPD_DIR}/tests/${TEST_SUITE}/${TEST_CASE}.in
        ${EXEC} "${ISPD_DIR}/tests/${TEST_SUITE}/${TEST_CASE}.in" $@ < "${ISPD_DIR}/tests/${TEST_SUITE}/${TEST_CASE}.out"
        echo ""
    fi
fi
