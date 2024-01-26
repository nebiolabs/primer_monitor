import json
import sys

if len(sys.argv) < 2:
    sys.stderr.write("usage: python extract_dataset_version.py <dataset JSON>\n")
    sys.exit(1)

with open(sys.argv[1]) as f:
    data = json.load(f)

print(data['version']['tag'])