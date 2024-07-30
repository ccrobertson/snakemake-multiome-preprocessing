#!/usr/bin/env python

import os
import sys
import json
import re

FASTQ_FILES = sys.argv[1:]

def parse_fastq_name(f):
    RE = '^(.*)-GEX.*S\d+.*([IR]\d+)_001.fastq.gz'
    if not re.match(RE, os.path.basename(f)):
        raise ValueError(f'Could not parse filename {f}')
    library, read = re.match(RE, os.path.basename(f)).groups()
    return {'library': library, 'readgroup': 'rg1', 'read': read}


libraries = set([parse_fastq_name(f)['library'] for f in FASTQ_FILES])
libraries = {library: {'genome': ['hg38'], 'readgroups': {}} for library in libraries}

for f in FASTQ_FILES:
    info = parse_fastq_name(f)
    library = info['library']
    readgroup = '{library}_{readgroup}'.format(**info)
    if readgroup not in libraries[library]['readgroups']:
        libraries[library]['readgroups'][readgroup] = dict()
    read = info['read']
    if read in ['I1', 'I2']:
        continue
    read = read.replace('R', '')
    libraries[library]['readgroups'][readgroup][read] = f

CONFIG = {
    'libraries': libraries
}
print(json.dumps(CONFIG, sort_keys = True, indent = 4))
