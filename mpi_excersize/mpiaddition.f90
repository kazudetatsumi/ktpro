program mpiaddition
use mpi
implicit none
integer :: mpi_rank,mpi_size,ierr,root=0
integer :: nvals, i, tot_nvals
real    :: sum_local_vals, sum_vals_reduce, sum_vals_gather
real, dimension(:), allocatable :: vals, all_vals
integer, dimension(:), allocatable :: all_nvals, disp_nvals

call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world,mpi_rank,ierr)
call mpi_comm_size(mpi_comm_world,mpi_size,ierr)

! Every process generates this number of data
nvals = mod(1000 - mpi_rank + 1, mpi_size+5) + 1
print *, "mpi_size:",mpi_size
print *, "nvals from my proc #", mpi_rank, ":", nvals

! Allocate and fill the array on each process
allocate( vals (nvals) )
do i=1,nvals
    vals(i) = i
end do

print*,mpi_rank,mpi_size,nvals
call mpi_barrier(mpi_comm_world, ierr)

! allocate recvcounts and displs arrays
if (mpi_rank .eq. root) then
    allocate(all_nvals(mpi_size))
    allocate(disp_nvals(mpi_size))
end if

! Gather the number of data to receive from each process
call mpi_gather(nvals, 1, MPI_INT, all_nvals, 1, MPI_INT, root, mpi_comm_world, ierr)

! Calculate displacements
if (mpi_rank .eq. root) then
    tot_nvals = 0
    do i = 1, mpi_size
        disp_nvals(i) = tot_nvals
        tot_nvals = tot_nvals + all_nvals(i)
    end do
    ! allocate receive buffer
    allocate(all_vals(tot_nvals))

    print*, "[root] nvals", all_nvals
    print*, "[root] disps", disp_nvals
    print*, "[root] size(all_vals)", size(all_vals)
end if

! Gather the final data
call mpi_gatherv(vals, nvals, MPI_REAL, all_vals, all_nvals, disp_nvals, MPI_REAL, root, mpi_comm_world, ierr)

! Perform a final sum on process 0
if (mpi_rank .eq. root) then
    sum_vals_gather = sum(all_vals)
end if

! Compare the sum of gathered values with the sum obtained using mpi_reduce
do i = 1, nvals
    sum_local_vals = sum(vals)
end do

call mpi_reduce(sum_local_vals, sum_vals_reduce, 1, MPI_REAL, MPI_SUM, root, mpi_comm_world, ierr)

if (mpi_rank .eq. root) then
    print*, "[root] sum by gather: ",  sum_vals_gather
    print*, "[root] sum by reduce: ",  sum_vals_reduce
end if

call mpi_finalize(ierr)
deallocate (vals)
if (mpi_rank .eq. root) then
    deallocate (all_nvals)
    deallocate (disp_nvals)
    deallocate (all_vals)
end if
end program
