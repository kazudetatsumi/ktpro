#!/usr/bin/env python
import Image
import sys

argvs = sys.argv
im = Image.open(argvs[1])
values=argvs[1].split(".")
values2=values[0].split("test")
w,h = im.size
if int(values2[1])%2==1: 
        print "kisuu"
	box = (80,160, w-50, h-110)
else:
	box = (40,160, w-70, h-110)

print w,h,values[0]
#im.crop((50,50,w-100,h-100)).save(values[2])
im.crop(box).save(values[0]+"rs.png")
#im.crop((3,125,w-3,h-75)).save(values[2])
