#!/usr/bin/python
import sys

ifile = file( sys.argv[1], 'r' )
for line in ifile:
#	words = line.split(',')
#	if words[0] == 't':
#		print 'ZONE I=2, J=2'
#	else:
#		print line,

     words= line.split(",")
     if words[0] == "ZONE":
     	print "t, 0"
     else:
        print line,
ifile.close()
