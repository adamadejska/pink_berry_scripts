# Coed from: https://www.geeksforgeeks.org/bellman-ford-algorithm-dp-23/
# Python3 program for Bellman-Ford's single source
# shortest path algorithm.
 
# Class to represent a graph
class Graph:
 
    def __init__(self, vertices):
        self.V = vertices # No. of vertices
        self.graph = []
 
    # function to add an edge to graph
    def addEdge(self, u, v, w):
        self.graph.append([u, v, w])
         
    # utility function used to print the solution
    def printArr(self, dist):
        print("Vertex Distance from Source")
        for i in range(self.V):
            print("{0}\t\t{1}".format(i, dist[i]))
     
    # The main function that finds shortest distances from src to
    # all other vertices using Bellman-Ford algorithm. The function
    # also detects negative weight cycle
    def BellmanFord(self, src):
 
        # Step 1: Initialize distances from src to all other vertices
        # as INFINITE
        dist = [float("Inf")] * self.V
        dist[src] = 0
        p = [-1] * self.V      # Solve what nodes create the shortest path
 
 
        # Step 2: Relax all edges |V| - 1 times. A simple shortest
        # path from src to any other vertex can have at-most |V| - 1
        # edges
        for _ in range(self.V - 1):
            # Update dist value and parent index of the adjacent vertices of
            # the picked vertex. Consider only those vertices which are still in
            # queue
            for u, v, w in self.graph:
                if dist[u] != float("Inf") and dist[u] + w < dist[v]:
                        dist[v] = dist[u] + w
                        p[v] = u
                         
        # Print shortest distance to the last node
        print('shortest distance to last node: %d' %dist[self.V-1])
        
        
        # Print path from start to the last node
        path = []
        v = self.V-1
        while v != -1:
            path.append(v)
            v = p[v]

        path.reverse()
        print(path)
 
g = Graph(14)
# (start, end, weight)
g.addEdge(0, 1, 0)
g.addEdge(0, 2, 23)
g.addEdge(0, 3, 27)
g.addEdge(1, 4, 0)
g.addEdge(1, 5, 10)
g.addEdge(2, 4, 10)
g.addEdge(2, 5, 24)
g.addEdge(2, 6, 10)
g.addEdge(3, 5, 10)
g.addEdge(3, 6, 29)
g.addEdge(4, 7, 1)
g.addEdge(4, 8, 10)
g.addEdge(5, 7, 10)
g.addEdge(5, 8, 1)
g.addEdge(5, 9, 10)
g.addEdge(6, 8, 10)
g.addEdge(6, 9, 5)
g.addEdge(7, 10, 6)
g.addEdge(7, 11, 10)
g.addEdge(8, 10, 10)
g.addEdge(8, 11, 0)
g.addEdge(8, 12, 10)
g.addEdge(9, 11, 10)
g.addEdge(9, 12, 3)
g.addEdge(10, 13, 10)
g.addEdge(11, 13, 0)
g.addEdge(12, 13, 1)

# Print the solution
g.BellmanFord(0)

