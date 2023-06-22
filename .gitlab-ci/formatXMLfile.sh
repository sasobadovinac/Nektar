#!/bin/bash

# Indentation level
if [[ $XML_INDENT ]]; then
    indent=$XML_INDENT
else
    indent="    "
fi

# Filename
if [[ $1 ]]; then
    filename=$1
else
    echo "Please specify xml filename"
    exit 1 # exit with error
fi


# Exit if not xml file
[[ "$filename" =~ ^.*\.xml$ ]] || exit 1

# Format file
rm -f tmp
tab=""
dos2unix $filename 2> /dev/null
while read -r line || [[ -n "$line" ]]; do
    # Print empty line
    if [[ $line == "" ]]; then
        echo "" >> tmp
        continue
    fi

    # Just print line if comment
    if [[ $line == "<!--"*"-->" ]]; then
        echo "${tab}${line}" >> tmp
        continue
    fi

    # Remove tab
    if [[ $line == "</"*">"* ]] || [[ $line == "/>"* ]] || 
       [[ $line == "-->" ]]; then
        tab=${tab::-4}
    fi

    # Print line
    echo "${tab}${line}" >> tmp

    # Add tab
    if [[ $line == "<!--" ]] || [[ $line == "<!--"* ]]; then
        # Comment tag
        tab=${tab}"    "
    elif [[ $line == "<"* ]] && [[ $line == !("</"*) ]]; then
        # Opening tag
        if [[ $line == !("<?"*"?>") ]]; then
            tab=${tab}"    "
        fi
    fi

    # Remove tab
    if  [[ $line == "-->" ]]; then
        # Comment tag
        continue
    elif  [[ $line == *"-->" ]] && [[ $line == !(*">"*"-->") ]] && 
          [[ $line == !(*"</"*"-->") ]]; then
        tab=${tab::-4}
    elif [[ $line == "<"*"</"*">"* ]] || [[ $line == "<"*"/>"* ]] || 
     ( [[ $line == !("<"*) ]] && [[ $line == !("/>"*) ]] && [[ $line == *"/>"* ]] ); then
        # Closing tag
        tab=${tab::-4}
    fi

done < $filename

# Rename file
if [[ $2 == !(--no-output) ]]; then
    if [[ $2 ]]; then
        mv tmp $2
    else 
        cat tmp
    fi
else
    if cmp -s tmp $filename; then 
       rm -f tmp
    else
        echo "$filename not formatted properly"
        rm -f tmp
        exit 2 # exit with error
    fi
fi
rm -f tmp

# Exit script
exit 0
