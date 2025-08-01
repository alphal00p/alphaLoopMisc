#/usr/bin/env python3
import sys
import json
data = None
with open(f"./{sys.argv[1]}",'r') as f:
    data = eval(f.read())
print(json.dumps(data[0]))
