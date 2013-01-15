## explict queue
from multiprocessing import Process, Queue
from collections import defaultdict


def computeCopyNum(queue, val):
    queue.put(val) # can also put a tuple of thread-id and value if we would like to

procs=list()

queue = Queue() ## with Queue
for i in range(1,3):
    p = Process(target=computeCopyNum, args=(queue, i))
    procs.append(p)
    p.start()

for _ in procs:
    val = queue.get()
    # do whatever with val
    print val

for p in procs:
    p.join()



## for multiple output value
#def slave(queue):
#    for i in range(128): # just for example
#        val = #some calculated result
#        queue.put(val)
#
#    queue.put(None) # add a sentinel value to tell the master we're done
#
#queue = Queue()
#
## spawn 32 slave processes
#num_procs = 32
#procs = [Process(target=slave, args=(queue, )) for _ in range(num_procs)]
#for proc in procs: 
#    proc.start()
#
#finished = 0
#while finished < num_procs:
#    item = queue.get()
#    if item is None: 
#        finished += 1
#    else: 
#        # do something with item
#
#for proc in procs:
#    proc.join()
#
### implicit memory copying between process using Manager, not recommend
#manager = Manager()
#a = manager.dict({1: 0, 2: -1})
#b = manager.list((1, 2, 3))
