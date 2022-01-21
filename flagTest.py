#!/usr/bin/env python


import os
import sys

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--noFilter', dest="filterTrigger",default=True,action="store_false",required=False, help='Input files after running detector.py.')
args = parser.parse_args()

print("no Filter",args.filterTrigger)