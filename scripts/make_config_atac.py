#!/usr/bin/env python

import os
import sys
import json
import re

FASTQ_FILES = sys.argv[1:]

def parse_fastq_name(f):
    RE = '^(.*-VD-\d+)-ATAC_(S\d+)_([IR]\d+)_001.fastq.gz'
    if not re.match(RE, os.path.basename(f)):
        raise ValueError(f'Could not parse filename {f}')
    library, readgroup, read = re.match(RE, os.path.basename(f)).groups()
    return {'library': library, 'readgroup': readgroup, 'read': read}


libraries = set([parse_fastq_name(f)['library'] for f in FASTQ_FILES])
libraries = {library: {'genome': ['hg38'], 'readgroups': {}} for library in libraries}

for f in FASTQ_FILES:
    info = parse_fastq_name(f)
    library = info['library']
    readgroup = '{library}_{readgroup}'.format(**info)
    if readgroup not in libraries[library]['readgroups']:
        libraries[library]['readgroups'][readgroup] = dict()
    read = info['read']
    if read == 'R1':
        read = '1'
    elif read == 'R2':
        read = '2'
    elif read == 'I2':
        read = 'index'
    else:
        sys.stderr.write(f'Skipping file {f}\n')
        continue
    libraries[library]['readgroups'][readgroup][read] = f

CONFIG = {'libraries': libraries}
print(json.dumps(CONFIG, sort_keys = True, indent = 4))
