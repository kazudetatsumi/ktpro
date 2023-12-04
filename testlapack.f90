module linearalgebra
    implicit none
    interface eigen
        module procedure eigen_dble_sym_eigenvalues
        module procedure eigen_dble_sym_both
        module procedure eigen_complex_her_eigenvalues
        module procedure eigen_complex_her_both
    end interface eigen

    contains

    subroutine eigen_dble_sym_eigenvalues(H,e)
        implicit none
        real(8),intent(in)::H(:,:)
        real(8),intent(out)::e(:)

        real(8),allocatable::V(:,:),work(:)
        integer::n,lwork,info

        n = ubound(H,1)
        lwork = 3*n-1
        write(*,*) n
        allocate(V(1:n,1:n))
        allocate(work(lwork))
        V = H
        call dsyev('N', 'U', n, V, n, e, work, lwork, info)
        if (info .ne. 0) then
            write(*,*) "error! in eigen. info = ",info
        end if
        deallocate(V,work)

        return

    end subroutine eigen_dble_sym_eigenvalues


    subroutine eigen_dble_sym_both(H,e,V)
        implicit none
        real(8),intent(in)::H(:,:)
        real(8),intent(out)::e(:)
        real(8),intent(out)::V(:,:)

        real(8),allocatable::work(:)
        integer::n,lwork,info


        n = ubound(H,1)
        lwork = 3*n-1
        allocate(work(lwork))
        V = H
        call dsyev('V', 'U', n, V, n, e, work, lwork, info)
        if (info .ne. 0) then
            write(*,*) "error! in eigen. info = ",info
        end if
        deallocate(work)

        return

    end subroutine eigen_dble_sym_both  

    subroutine eigen_complex_her_eigenvalues(H,e)
        implicit none
        complex(8),intent(in)::H(:,:)
        real(8),intent(out)::e(:)

        complex(8),allocatable::V(:,:),work(:)
        integer::n,lwork,info
        integer,allocatable::rwork(:)

        n = ubound(H,1)
        lwork = 3*n-1
        allocate(V(1:n,1:n))
        allocate(work(lwork))
        allocate(rwork(2*n-2))
        V = H
        call zheev('N', 'U', n, V, n, e, work, lwork, rwork, info)
        if (info .ne. 0) then
            write(*,*) "error! in eigen. info = ",info
        end if

        deallocate(work,rwork,V)

        return

    end subroutine eigen_complex_her_eigenvalues  

    subroutine eigen_complex_her_both(H,e,V)
        implicit none
        complex(8),intent(in)::H(:,:)
        real(8),intent(out)::e(:)

        complex(8),intent(out)::V(:,:)
        complex(8),allocatable::work(:)
        integer::n,lwork,info
        integer,allocatable::rwork(:)

        n = ubound(H,1)
        lwork = 3*n-1
        allocate(work(lwork))
        allocate(rwork(2*n-2))
        V = H
        call zheev('V', 'U', n, V, n, e, work, lwork, rwork, info)
        if (info .ne. 0) then
            write(*,*) "error! in eigen. info = ",info
        end if

        deallocate(work,rwork)

        return

    end subroutine eigen_complex_her_both  

end module linearalgebra


subroutine test()
    use linearalgebra
    implicit none
    integer::n
    real(8),allocatable::Hs(:,:)
    real(8),allocatable::evec(:)
    real(8),allocatable::Vs(:,:)
    complex(8),allocatable::V(:,:)
    complex(8),allocatable::H(:,:)
    real(8),allocatable::Hr(:,:)
    real(8),allocatable::Hi(:,:)
    complex,parameter::ci = (0d0,1d0)
    n = 40000
    !allocate(Hs(1:n,1:n))

    allocate(V(1:n,1:n))
    !allocate(Vs(1:n,1:n))
    allocate(evec(1:n))

    call random_number(Hs)

    !Hs = (Hs + transpose(Hs))/2
    !!write(*,*) Hs

    !call eigen(Hs,evec)
    !write(*,*) evec

    !call eigen(Hs,evec,Vs)
    !write(*,*) evec

    allocate(H(1:n,1:n))
    allocate(Hi(1:n,1:n))
    allocate(Hr(1:n,1:n))
    call random_number(Hr)
    call random_number(Hi)
    H = Hr + ci*Hi

    H = (H + conjg(transpose(H)))/2

    call eigen(H,evec)
    write(*,*) evec

    !call eigen(H,evec,V)
    !write(*,*) evec



    return
end subroutine


program main
    use linearalgebra
    implicit none
    call test()
end program
