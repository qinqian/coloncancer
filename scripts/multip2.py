#!/usr/bin/env python
"""
"""
# QinQian

##

from multiprocessing import Pool, Manager

list = [1,2,3,10]
d = {}
def test(x, d):
    for xx in range(100):
        for xxx in range(100):
            d[x]=xx*xxx



if __name__ == '__main__':
    pool = Pool(processes=4)
    mgr = Manager()
    d = mgr.dict()
    for N in list:
        pool.apply_async(test, (N, d))

    # Mark pool as closed -- no more tasks can be added.
    pool.close()

    # Wait for tasks to exit
    pool.join()

    # Output results
    print d
