#!/usr/bin/env python

from binSortCountMapReduce import combiner
import sys

if __name__ == '__main__':
	combiner(sys.stdin, sys.stdout)