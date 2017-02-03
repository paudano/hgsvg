#!/usr/bin/env python

import sys

import argparse

ap = argparse.ArgumentParser(description="Given a vcf merged by Zev Kronenberg, output in bed5 format, with the name of the event the callers supporting the event")
ap.add_argument("--vcf", help="Input vcf file",default="/dev/stdin")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()
