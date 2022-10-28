"""
Assign 03 - <Zack Khalidov>

Directions:
    * Complete the graph algorithm functions given below. Note that it may be
      helpful to define auxiliary/helper functions that are called from the
      functions below.  Refer to the README.md file for additional info.

    * NOTE: As with other assignments, please feel free to share ideas with
      others and to reference sources from textbooks or online. However, be sure
      to **cite your resources in your code. Also, do your best to attain a
      reasonable grasp of the algorithm that you are implementing as there will
      very likely be questions related to it on quizzes/exams.

    * NOTE: Remember to add a docstring for each function, and that a reasonable
      coding style is followed (e.g. blank lines between functions).
      Your program will not pass the tests if this is not done!
"""
# could be useful for Dijkstra
from queue import PriorityQueue

# for timing checks
import time

# use a very large number as placeholder for infinity
import sys
INF = sys.maxsize


def adjMatFromFile(filename):
    """ Create an adj/weight matrix from a file with verts, neighbors, and weights. """
    f = open(filename, "r")
    n_verts = int(f.readline())
    print(f" n_verts = {n_verts}")
    adjmat = [[INF] * n_verts for i in range(n_verts)]
    for i in range(n_verts):
        adjmat[i][i] = 0
    for line in f:
        int_list = [int(i) for i in line.split()]
        vert = int_list.pop(0)
        assert len(int_list) % 2 == 0
        n_neighbors = len(int_list) // 2
        neighbors = [int_list[n] for n in range(0, len(int_list), 2)]
        distances = [int_list[d] for d in range(1, len(int_list), 2)]
        for i in range(n_neighbors):
            adjmat[vert][neighbors[i]] = distances[i]
    f.close()
    return adjmat


def printAdjMat(mat, width=3):
    """ Print an adj/weight matrix padded with spaces and with vertex names. """
    res_str = '     ' + ' '.join([str(v).rjust(width, ' ') for v in range(len(mat))]) + '\n'
    res_str += '    ' + '-' * ((width + 1) * len(mat)) + '\n'
    for i, row in enumerate(mat):
        row_str = [str(elem).rjust(width, ' ') if elem < INF else ' 999' for elem in row]
        res_str += ' ' + str(i).rjust(2, ' ') + ' |' + ' '.join(row_str) + "\n"
    print(res_str)


def floyd(W):
    """ Carry out Floyd's algorithm using W as a weight/adj matrix. """
    vert_count = len(W)

    # initialize matrix for floyd's algo
    matrix = [[INF] * vert_count for i in range(vert_count)]
    for i in range(vert_count):
        matrix[i][i] = 0
    for i in range(vert_count):
        for j in range(vert_count):
            matrix[i][j] = W[i][j]

    for temp in range(vert_count):
        for i in range(vert_count):
            for j in range(vert_count):
                matrix[i][j] = min(matrix[i][j], matrix[i][temp] + matrix[temp][j])

    return matrix


def dijkstra_w_pri_queue(W, sv):
    """ Dijkstra's using a priority queue w/ adj matrix W and sv as starting vert. """
    vert_count = len(W)
    visited = [False] * vert_count
    distance = [INF] * vert_count
    distance[sv] = 0
    pq = PriorityQueue()

    # initialize matrix for dijkstra_w_pri_queue
    matrix_d_2 = [[INF] * vert_count for i in range(vert_count)]
    for i in range(vert_count):
        matrix_d_2[i][i] = 0
    for i in range(vert_count):
        for j in range(vert_count):
            matrix_d_2[i][j] = W[i][j]

    # put weight and vertex into priority queue
    for i in range(vert_count):
        pq.put((matrix_d_2[sv][i], i))

    # run until pq is empty to find optimal path
    while (not pq.empty()):
        node = pq.get()
        vertex = int(node[1])
        visited[vertex] = True

        for i in range(len(matrix_d_2[sv])):
            if(matrix_d_2[vertex][i] > 0):
                if(not visited[i]):
                    if(distance[i] > distance[vertex] + matrix_d_2[vertex][i]):
                        distance[i] = distance[vertex] + matrix_d_2[vertex][i]
                        pq.put((distance[i], i))

    return distance


def min_distance(matrix_d_1, distance, visited):
    """ This method is used to find the minimum distance and return the index. """
    min = INF

    for i in range(len(matrix_d_1)):
        if distance[i] < min and not visited[i]:
            min = distance[i]
            index = i

    return index


def dijkstra_w_array(W, sv):
    """ Dijkstra's using an array w/ adj matrix W and sv as starting vert. """
    vert_count = len(W)
    visited = [False] * vert_count
    distance = [INF] * vert_count
    distance[sv] = 0

    # initialize matrix for dijkstra_w_array
    matrix_d_1 = [[INF] * vert_count for i in range(vert_count)]
    for i in range(vert_count):
        matrix_d_1[i][i] = 0
    for i in range(vert_count):
        for j in range(vert_count):
            matrix_d_1[i][j] = W[i][j]

    # loop and find shortest distance to all verts from sv
    for count in range(len(matrix_d_1[sv])):
        min_index = min_distance(matrix_d_1[sv], distance, visited)
        visited[min_index] = True

        for i in range(len(matrix_d_1[sv])):
            if(matrix_d_1[min_index][i] > 0):
                if(not visited[i]):
                    if(distance[i] > distance[min_index] + matrix_d_1[min_index][i]):
                        distance[i] = distance[min_index] + matrix_d_1[min_index][i]

    return distance


def assign03_main():
    """ Demonstrate the functions, starting with creating the graph. """
    g = adjMatFromFile("graph_20verts.txt")
    if len(g) <= 20:
        printAdjMat(g, width=4)

    # Run Floyd's algorithm
    start_time = time.time()
    res_floyd = floyd(g)
    elapsed_time_floyd = time.time() - start_time
    if len(g) <= 20 and res_floyd:
        printAdjMat(res_floyd, width=4)

    # Run Dijkstra's for a single starting vertex, v2
    start_vert = 2
    res_dijkstra_single1 = dijkstra_w_pri_queue(g, start_vert)
    print(f" dijkstra_w_pri_queue (for starting vert {start_vert})= : {res_dijkstra_single1}")
    res_dijkstra_single2 = dijkstra_w_array(g, start_vert)
    print(f" dijkstra_w_array (for starting vert {start_vert})= : {res_dijkstra_single2}")
    if len(g) <= 20 and res_floyd:
        print(f" floyd's algorithm    (for starting vert {start_vert}): {res_floyd[start_vert]}")
        print(f" dijkstra_w_pri_queue (for starting vert {start_vert}): {res_dijkstra_single1}")
        print(f" dijkstra_w_array     (for starting vert {start_vert}): {res_dijkstra_single2}")
    elif res_floyd:
        print(f" adjacency matrix for starting vert {start_vert}      : {g[start_vert][:20]}")
        print(f" floyd's algorithm    (for starting vert {start_vert}): {res_floyd[start_vert][:20]}")
        print(f" dijkstra_w_pri_queue (for starting vert {start_vert}): {res_dijkstra_single1[:20]}")
        print(f" dijkstra_w_array     (for starting vert {start_vert}): {res_dijkstra_single2[:20]}")

    # Check that the two produce same results by comparing result from
    # dijkstra with the corresponding row from floyd's output matrix
    assert res_floyd[start_vert] == res_dijkstra_single1
    assert res_floyd[start_vert] == res_dijkstra_single2

    # Run Dijkstra's overall starting points (note this is not the intended way
    # to utilize this algorithm, however we are using it as point of comparion).
    res_dijkstra1 = [[None] * len(g) for i in range(len(g))]
    start_time = time.time()
    for sv in range(len(g)):
        res_dijkstra1[sv] = dijkstra_w_pri_queue(g, sv)
    elapsed_time_dijkstra1 = time.time() - start_time

    res_dijkstra2 = [[None] * len(g) for i in range(len(g))]
    start_time = time.time()
    for sv in range(len(g)):
        res_dijkstra2[sv] = dijkstra_w_array(g, sv)
    elapsed_time_dijkstra2 = time.time() - start_time

    # Compare times, Dijkstra's should be slower (for non-trivial sized graphs)
    print(f"  Dijkstra's w/ pri queue elapsed time (for all starts): {elapsed_time_dijkstra1:.2f}")
    print(f"  Dijkstra's w/ array elapsed time (for all starts): {elapsed_time_dijkstra2:.2f}")
    print(f"  Floyd's elapsed time: {elapsed_time_floyd:.2f}")

    # Double check again that the results are all the same
    assert res_floyd == res_dijkstra1
    assert res_floyd == res_dijkstra2


# Check if the program is being run directly (i.e. not being imported)
if __name__ == '__main__':
    assign03_main()
