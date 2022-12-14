#!/bin/sh

url=$(Rscript /home/docker/CommonFiles/Code/MLSTurlParser.R)

(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${1}; echo '"}') | \
curl -s -H "Content-Type: application/json" -X POST ${url} -d @- > \
    testsh_rmlst.json
