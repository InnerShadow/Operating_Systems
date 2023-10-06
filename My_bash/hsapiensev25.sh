#!/usr/bin/env bash

grep -P "^[\w\.]+\s+\S+\sgene\s" gencode.v25.primary_assembly.annotation.gtf | cut -f 9 | awk '{split($0,a,"; "); print a[2]}' | sort | uniq | wc -l

