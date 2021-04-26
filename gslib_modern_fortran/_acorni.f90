
program test_random_number

    use OMP_LIB
    integer, parameter :: l = 3
    real :: rr(l)
    integer, allocatable :: seed(:)
    integer :: n
    
    call random_seed(size = n)
    allocate(seed(n*3))
    seed = 878768768
    call random_seed(put=seed)

    !$OMP PARALLEL private(rr)

        call random_number(rr)

        print *, 'Process', OMP_GET_THREAD_NUM(), rr

    !$OMP END PARALLEL
    
    ! this parallel code will produce this code in two runs:
    ! 
    !    PS C:\Users\AMartinez\Desktop\fortran_test> ./random.exe
    !     Process           2  0.162345409      0.840058029       9.94555354E-02
    !     Process           0  0.334234238      0.642744005      0.283221900
    !     Process           4   7.85911679E-02  0.241505861      0.400358200
    !     Process           3  0.766133487      0.591100812      0.382381201
    !     Process           6   1.73936486E-02  0.396930397      0.280517757
    !     Process           5  0.186948657      0.335481524      0.356720865
    !     Process           7  0.252353966      0.576505065       2.22165585E-02
    !     Process           1  0.844142258      0.983546615      0.513327420
    !    PS C:\Users\AMartinez\Desktop\fortran_test> ./random.exe
    !     Process           6  0.162345409      0.840058029       9.94555354E-02    
    !     Process           0  0.334234238      0.642744005      0.283221900
    !     Process           3   7.85911679E-02  0.241505861      0.400358200
    !     Process           2   1.73936486E-02  0.396930397      0.280517757
    !     Process           1  0.766133487      0.591100812      0.382381201
    !     Process           5  0.186948657      0.335481524      0.356720865
    !     Process           7  0.252353966      0.576505065       2.22165585E-02
    !     Process           4  0.844142258      0.983546615      0.513327420
    !
    ! the sequences of random numbers are the same but obtained in a different order 
    ! 
    ! for conditional simulations, you can generate the entire array in one call
    
end program