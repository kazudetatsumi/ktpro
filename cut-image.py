#!/usr/bin/env python
import Image
import sys

argvs = sys.argv
im = Image.open(argvs[1])
values=argvs[1].split(" ")
w,h = im.size
im.crop((3,360,w-3,h-310)).save(values[2])
#im.crop((3,80,w-3,h-33)).save(values[2])
#im.crop((3,125,w-3,h-75)).save(values[2])
