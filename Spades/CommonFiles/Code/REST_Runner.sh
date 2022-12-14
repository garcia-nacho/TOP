#!/bin/sh

# 1-Sequence
# 2-URL
# 3-OUTPUT

(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${1}; echo '"}') | \
curl -s -H "Content-Type: application/json" -X POST ${2} -d @- > ${3}
