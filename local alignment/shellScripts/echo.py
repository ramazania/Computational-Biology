"""
    This script will take in a variable number of command line arguments
    and will print the the total number of arguments along with each
    one on a different line.

    NOTE: This script is written for python 2.7.

    Here is an example run of the program and the output it will produce
    
    $ python echo.py 1 2 3
    Total number of arguments:  3
    Argument  0  :  1
    Argument  1  :  2
    Argument  2  :  3
 """

import sys

num_params = len(sys.argv)-1

print "Total number of arguments:",num_params

for v in xrange(num_params):
    print "Argument",v,":",sys.argv[v+1]
