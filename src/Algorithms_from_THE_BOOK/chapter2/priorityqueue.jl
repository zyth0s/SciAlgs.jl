# Author: Kenneth Lange @ University of California, Los Angeles 

using DataStructures

pq = PriorityQueue() # empty queue
pq['a'] = 10 # enqueue or push
pq['b'] = 5 
pq['c'] = 15
peek(pq)
dequeue!(pq) # dequeue or pop
