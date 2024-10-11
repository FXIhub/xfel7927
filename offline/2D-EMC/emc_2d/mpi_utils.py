from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()



def allgather(array, axis):
    for r in range(size):
        indices = list(range(r, array.shape[axis], size))
        s = len(array.shape) * [None,]
        s[axis] = indices
        s = tuple(s)
        array[s] = comm.bcast(array[s], root=r)
