program main
 use mpi 
 implicit none
 integer :: PETOT, my_rank, ierr, i
 double precision, dimension(5) :: VECp, VECs
 double precision sumAR, sum0, sumR

 call MPI_INIT(ierr)
 call MPI_COMM_SIZE(MPI_COMM_WORLD, PETOT, ierr)
 call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

 sumAR = 0
 sumR = 0
 do i=1,5
    VECp(i) = 2.d0
    VECs(i) = 3.d0
 enddo

 sum0 = 0
 do i=1,5
    sum0=sum0+VECp(i)*VECs(i)
 enddo

 call MPI_REDUCE(sum0, sumR, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
 call MPI_allREDUCE(sum0, sumAR, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

 print *, 'before BCAST', my_rank, sumAR, sumR, sum0

 call MPI_BCAST (sumR, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

 print *, 'after BCAST', my_rank, sumAR, sumR, sum0




 call MPI_FINALIZE(ierr)
 stop
 end
