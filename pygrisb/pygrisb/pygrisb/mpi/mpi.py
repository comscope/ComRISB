from mpi4py import MPI


def get_k_range_list(numk):
    comm = MPI.COMM_WORLD
    ncpu = comm.Get_size()
    nk_per_cpu = numk//ncpu
    nk_left = numk%ncpu
    nk_list = [nk_per_cpu if i >= nk_left else nk_per_cpu+1
            for i in range(ncpu)]
    k_range_list = [[sum(nk_list[:i]), sum(nk_list[:i+1])] \
            for i in range(ncpu)]
    return k_range_list


def get_myk_range(numk):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ncpu = comm.Get_size()
    nk_per_cpu = numk//ncpu
    nk_left = numk%ncpu
    if rank < nk_left:
        k_start = (nk_per_cpu+1)*rank
        k_end = (nk_per_cpu+1)*(rank+1)
    else:
        k_start = nk_per_cpu*rank + nk_left
        k_end = k_start + nk_per_cpu
    return [k_start, k_end]


if __name__ == "__main__":
    numk = 111
    k_range_list = get_k_range_list(numk)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        print(k_range_list)
    k_range = get_myk_range(numk)
    print(f"rank = {rank}: k_range = {k_range}")
