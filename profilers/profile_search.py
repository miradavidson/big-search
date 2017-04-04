#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Profile the Matcher Module
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals


def profile_genome_search():
    """Profile the sort bucket step"""
    genome_file  = 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
    pattern = b'TGGATGTGAAATGAGTCAAG'
    maxmismatch = 3
    # TODO Your implementation

def main():
    profile_genome_search()

if __name__ == '__main__':
    main()