#!/bin/bash
error=0
for file in $(find $TARGETS -type f); do
    filename=${file##*/}
    if [[ ! $(grep -F "MIT License" $file) ]]; then
        if [[ "${filename##*.}" == "h" ]]; then
            echo $file does not contains MIT License
            error=1
        fi
        if [[ "${filename##*.}" == "cpp" ]]; then
            echo $file does not contains MIT License
            error=1
        fi
        if [[ "${filename##*.}" == "hpp" ]]; then
            echo $file does not contains MIT License
            error=1
        fi
    fi
done

if [[ $error -eq 1 ]]; then
    exit 1 # exit with an error
fi
