#!/bin/bash

EXEC=$1; shift;
ISPD_DIR=$1; shift;
TEST_SUITE=$1; shift;

if [ "all" == ${TEST_SUITE} ]
then
    for i in ${ISPD_DIR}/tests/*/; do
        for j in ${i}*/; do
            for k in ${j}*.in; do
                if [ ! -e ${k} ]
                then
                    continue
                fi
                echo ${k}
                ${EXEC} "${k}" $@ < "${k%.*}.out"
                echo ""
            done
        done
    done
else
    TEST_CASE=$1; shift;
    if [[ "all" == ${TEST_CASE} || "" == ${TEST_CASE} ]]
    then
        for i in ${ISPD_DIR}/tests/*/; do
            for j in ${i}${TEST_SUITE}/*.in; do
                if [ ! -e ${j} ]
                then
                    continue
                fi
                echo ${j}
                ${EXEC} "${j}" $@ < "${j%.*}.out"
                echo ""
            done
        done
    else
        if [ -e ${ISPD_DIR}/tests/2d/${TEST_SUITE}/${TEST_CASE}.in ]
        then
            echo ${ISPD_DIR}/tests/2d/${TEST_SUITE}/${TEST_CASE}.in
            ${EXEC} "${ISPD_DIR}/tests/2d/${TEST_SUITE}/${TEST_CASE}.in" $@ < "${ISPD_DIR}/tests/2d/${TEST_SUITE}/${TEST_CASE}.out"
            echo ""
        fi

        if [ -e ${ISPD_DIR}/tests/3d/${TEST_SUITE}/${TEST_CASE}.in ]
        then
            echo ${ISPD_DIR}/tests/3d/${TEST_SUITE}/${TEST_CASE}.in
            ${EXEC} "${ISPD_DIR}/tests/3d/${TEST_SUITE}/${TEST_CASE}.in" $@ < "${ISPD_DIR}/tests/3d/${TEST_SUITE}/${TEST_CASE}.out"
            echo ""
        fi
    fi
fi
