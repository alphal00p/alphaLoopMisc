#!/bin/bash

function json2yaml {
    python3 -c 'import sys, yaml, json; print(yaml.dump(json.loads(sys.stdin.read())))'
}



for file in *.json; do cat $file | json2yaml > "${file%%.json}.yaml"; done
