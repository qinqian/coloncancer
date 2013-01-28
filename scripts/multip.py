#!/usr/bin/env python
"""
multiprocess or threads
threads: for a single CPU
process: for multiple CPU
Global Interpreter Lock(GIL) in a language effectively limits the amount of parallelism reachable through concurrency of a single interpreter process with multiple threads
* small collections of data, using subprocess.Popen
* compute bound, using multiprocessing module
* I/O bound, using threading module

"""
# QinQian

##
import time
import threading
from subprocess import Popen
from multiprocessing import Pool

def q(i):
    print i
    time.sleep(0.01)

def main():
    poll = Pool(processes=8)
    for i in range(100000, 101000):
        result = poll.apply_async(q, (i,))
        #result = poll.map(q, [i, i+1])
    poll.close()
    poll.join()
    #print len(result)  ## no successful()
    if result.successful():
        print "successful"

if __name__ == "__main__":
    main()
