#!/usr/bin/env python

from binSortCountMapReduce import reducer
import sys

if __name__ == '__main__':
	reducer(sys.stdin, sys.stdout)