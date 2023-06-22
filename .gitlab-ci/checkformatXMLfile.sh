#!/bin/bash
XML_INDEX="    " find -type f -name "*.xml" -print0 | xargs -0 -I % $(dirname "$0")/formatXMLfile.sh % --no-output
exit $?
