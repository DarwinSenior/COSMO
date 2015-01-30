#!/usr/bin/env python

from binSortCountMapReduce import mapper
import sys

if __name__ == '__main__':
	mapper(sys.stdin, sys.stdout)