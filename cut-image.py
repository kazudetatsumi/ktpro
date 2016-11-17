#!/usr/bin/env python
import Image
import sys

argvs = sys.argv
im = Image.open(argvs[1])
values=argvs[1].split(".")
values2=values[0].split("test")
w,h = im.size
if int(values2[1])%2==0: 
	box = (70,150, w-60, h-110)
	box = (80,150, w-50, h-110)
else:
	box = (30,150, w-100, h-110)

print w,h,values[0]
#im.crop((50,50,w-100,h-100)).save(values[2])
im.crop(box).save(values[0]+"rs.jpg")
#im.crop((3,125,w-3,h-75)).save(values[2])
