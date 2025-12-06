!---------------------------------< Percolation >----------------------------------!
!                                                                                  !    
! Contains: Computation of the probability of percolation for different            !   
!           probabilities of occupation and number of Penrose substitutions        !
!                                                                                  !
! Last revision:    06/12/2025                                                     !                  
!                                                                                  !
!----------------------------------------------------------------------------------!


MODULE Percolation
    USE Penrose
    USE omp_lib
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: master_percolations

    INTEGER, PARAMETER :: nSims = 10000


    CONTAINS

    ! =======================================================================
    ! Subroutine that occupies grid and computes whether it is spanned
    ! =======================================================================
    SUBROUTINE percolate(n_nodes, neighbors, num_neighbors, boundary_flags, p_occup, is_spanned)
        IMPLICIT NONE
        INTEGER, INTENT(IN)                 :: n_nodes
        INTEGER, DIMENSION(:,:), INTENT(IN) :: neighbors
        INTEGER, DIMENSION(:), INTENT(IN)   :: num_neighbors    
        LOGICAL, DIMENSION(:,:), INTENT(IN) :: boundary_flags
        REAL(8), INTENT(IN)                 :: p_occup
        LOGICAL, INTENT(OUT)                :: is_spanned

        LOGICAL, ALLOCATABLE :: occupied(:)
        INTEGER, ALLOCATABLE :: clusters(:)
        LOGICAL, ALLOCATABLE :: cluster_flags(:,:)

        INTEGER :: i, j, k
        REAL(8) :: rand_num

        ALLOCATE(occupied(n_nodes))
        ALLOCATE(clusters(n_nodes))
        ALLOCATE(cluster_flags(2, n_nodes))

        is_spanned = .FALSE.
        occupied = .FALSE.
        cluster_flags = .FALSE.
        clusters = 0

        DO i = 1, n_nodes
            CALL RANDOM_NUMBER(rand_num)
            IF (rand_num < p_occup) THEN
                occupied(i) = .TRUE.
                clusters(i) = i

                cluster_flags(1, i) = boundary_flags(i, 1)
                cluster_flags(2, i) = boundary_flags(i, 2)
            END IF
        END DO


        DO i = 1, n_nodes

            IF (.NOT. occupied(i)) CYCLE

            DO j = 1, num_neighbors(i)
                k = neighbors(i, j)
                IF (occupied(k)) THEN
                    CALL unite_clusters(i,k)
                END IF
            END DO

        END DO

        
        DO i = 1, n_nodes
            IF (occupied(i)) THEN
                IF (clusters(i) == i) THEN
                    IF (cluster_flags(1, i) .AND. cluster_flags(2, i)) THEN
                        is_spanned = .TRUE.
                        EXIT
                    END IF
                END IF
            END IF
        END DO

        DEALLOCATE(occupied, clusters, cluster_flags)





        CONTAINS

            SUBROUTINE unite_clusters(node1, node2)
                INTEGER, INTENT(IN) :: node1, node2
                INTEGER :: root1, root2

                root1 = find_root(node1)
                root2 = find_root(node2)

                IF (root1 /= root2) THEN
                    clusters(root2) = root1
                    cluster_flags(1, root1) = cluster_flags(1, root1) .OR. cluster_flags(1, root2)
                    cluster_flags(2, root1) = cluster_flags(2, root1) .OR. cluster_flags(2, root2)
                END IF
            END SUBROUTINE unite_clusters


            RECURSIVE FUNCTION find_root(node) RESULT(root)
                INTEGER, INTENT(IN) :: node
                INTEGER :: root

                IF (clusters(node) /= node) THEN
                    clusters(node) = find_root(clusters(node))
                END IF
                root = clusters(node)
            END FUNCTION find_root


    END SUBROUTINE percolate




    ! =======================================================================
    ! Runs nSims simulations and returns probability of percolation 
    ! =======================================================================
    SUBROUTINE simulate_percolations(n_nodes, neighbors, num_neighbors, boundary_flags, p_occup, p_perc)
        IMPLICIT NONE
        INTEGER, INTENT(IN)                 :: n_nodes
        INTEGER, DIMENSION(:,:), INTENT(IN) :: neighbors
        INTEGER, DIMENSION(:), INTENT(IN)   :: num_neighbors    
        LOGICAL, DIMENSION(:,:), INTENT(IN) :: boundary_flags
        REAL(8), INTENT(IN)                 :: p_occup
        REAL(8), INTENT(OUT)                :: p_perc

        INTEGER :: i, n_spanned
        LOGICAL :: is_spanned

        INTEGER :: n_threads, max_seed_size, th_id
        INTEGER, ALLOCATABLE :: thread_seeds(:,:), temp_seed(:)

        n_spanned = 0


        !$OMP MASTER
        n_threads = OMP_GET_MAX_THREADS()
        CALL RANDOM_SEED(size = max_seed_size)
        ALLOCATE(thread_seeds(n_threads, max_seed_size))
        ALLOCATE(temp_seed(max_seed_size))
        CALL RANDOM_SEED(GET=temp_seed) 
        DO i = 1, n_threads
            temp_seed(1) = temp_seed(1) + (i * 1000) + (i * i * 37)
            CALL RANDOM_SEED(PUT=temp_seed)
            CALL RANDOM_SEED(GET=thread_seeds(i, :))
        END DO
        DEALLOCATE(temp_seed)
        !$OMP END MASTER
        !$OMP BARRIER


        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, th_id, is_spanned) REDUCTION(+:n_spanned)
            th_id = OMP_GET_THREAD_NUM() + 1
            CALL RANDOM_SEED(PUT=thread_seeds(th_id,:))

            !$OMP DO SCHEDULE(DYNAMIC)
            DO i = 1, nSims
                CALL percolate(n_nodes, neighbors, num_neighbors, boundary_flags, p_occup, is_spanned)
                IF (is_spanned) THEN
                    n_spanned = n_spanned + 1
                END IF
            END DO
            !$OMP END DO

        !$OMP END PARALLEL

        DEALLOCATE(thread_seeds)

        p_perc = REAL(n_spanned, 8) / REAL(nSims, 8)

    END SUBROUTINE simulate_percolations



    ! =======================================================================
    ! Runs simulations for given grid size & occupation probability 
    ! =======================================================================
    SUBROUTINE master_percolations(n_iter_vec, p_occup_vec, xtip, ytip)
        IMPLICIT NONE
        INTEGER, DIMENSION(:), INTENT(IN) :: n_iter_vec
        REAL(8), DIMENSION(:), INTENT(IN) :: p_occup_vec
        REAL(8), INTENT(IN)               :: xtip, ytip

        INTEGER :: i, k, n_iter, n_grids, n_probs, iounit
        REAL(8) :: p_occup
        REAL(8), ALLOCATABLE :: p_perc_mtx(:,:)

        CHARACTER(LEN=64), PARAMETER :: filename = 'percolation_data.txt'
        CHARACTER(LEN=256) :: message

        REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: A, B, C, D
        REAL(8), ALLOCATABLE :: nodes(:,:)
        INTEGER, ALLOCATABLE :: neighbours(:,:), num_neighbours(:)
        LOGICAL, ALLOCATABLE :: boundary_flags(:,:)

        n_grids = SIZE(n_iter_vec)
        n_probs = SIZE(p_occup_vec)
        ALLOCATE(p_perc_mtx(n_grids,n_probs))
        p_perc_mtx = -1.0d0

        CALL CONSOLE("Starting percolation simulations...")
        WRITE(message, '(A,I0,A)') "Total grid configurations: ", n_grids, "."

        DO k = 1, n_grids

            n_iter = n_iter_vec(k)
            WRITE(message, '(A,I0,A)') "Processing grid with: n_iterations = ", n_iter, "..."
            CALL CONSOLE(message)

            ALLOCATE(C(1,3,2))
            ALLOCATE(D(1,3,2))

            C(1,1,1) = 0.0d0
            C(1,1,2) = ytip
            C(1,2,1) = xtip
            C(1,2,2) = 0.0d0
            C(1,3,1) = -xtip
            C(1,3,2) = 0.0d0

            D(1,1,1) = 0.0d0
            D(1,1,2) = -ytip
            D(1,2,1) = -xtip
            D(1,2,2) = 0.0d0
            D(1,3,1) = xtip
            D(1,3,2) = 0.0d0

            DO i = 1, n_iter
                CALL penrose_substitution(A, B, C, D)
            END DO

            CALL construct_grid(A, B, C, D, n_iter, xtip, ytip, nodes, neighbours, num_neighbours)

            IF (ALLOCATED(A)) DEALLOCATE(A)
            IF (ALLOCATED(B)) DEALLOCATE(B)
            IF (ALLOCATED(D)) DEALLOCATE(D)
            IF (ALLOCATED(C)) DEALLOCATE(C)
            
            CALL generate_boundary(nodes, n_iter, xtip, ytip, boundary_flags)

            WRITE(message, '(A,I0,A)') "Grid generated with #nodes = ", SIZE(nodes, 1), "."
            CALL CONSOLE(message)
            CALL CONSOLE("Running Monte Carlo simulation...")

            DO i = 1, n_probs
                p_occup = p_occup_vec(i)

                CALL simulate_percolations(SIZE(nodes,1), neighbours, num_neighbours, boundary_flags, p_occup, p_perc_mtx(k,i))
                WRITE(message, '(A,F6.4,A,F6.4)') "      p_occup = ", p_occup, " => p_perc=", p_perc_mtx(k,i)
                CALL CONSOLE(message)


            END DO

            CALL CONSOLE("Monte Carlo finished for this grid size.")

            DEALLOCATE(nodes, neighbours, num_neighbours, boundary_flags)


        END DO

        CALL CONSOLE("All simulations finished. Exporting data...")

        OPEN(NEWUNIT=iounit, FILE=TRIM(filename), STATUS='REPLACE', ACTION='WRITE')

        WRITE(iounit, '(A10)', ADVANCE='NO') '# p_occup'
        DO k = 1, n_grids
            WRITE(iounit, '(2X, A9, I2, 1X)', ADVANCE='NO') '       n=', n_iter_vec(k)
        END DO
        WRITE(iounit, *)

        WRITE(iounit, '(F10.6)', ADVANCE='NO') 0.0d0
        DO k = 1, n_grids
            WRITE(iounit, '(2X, F12.10)', ADVANCE='NO') 0.0d0
        END DO
        WRITE(iounit, *)

        DO i = 1, n_probs
            WRITE(iounit, '(F10.6)', ADVANCE='NO') p_occup_vec(i)
            DO k = 1, n_grids
                WRITE(iounit, '(2X, F12.10)', ADVANCE='NO') p_perc_mtx(k,i)
            END DO
            WRITE(iounit, *)
        END DO

        WRITE(iounit, '(F10.6)', ADVANCE='NO') 1.0d0
        DO k = 1, n_grids
            WRITE(iounit, '(2X, F12.10)', ADVANCE='NO') 1.0d0
        END DO
        WRITE(iounit, *)

        CLOSE(iounit)
        CALL CONSOLE("Data succesfully exported.")

        DEALLOCATE(p_perc_mtx)


    END SUBROUTINE master_percolations




    SUBROUTINE CONSOLE(text)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: text
        CHARACTER(LEN=10) :: timestr
        INTEGER :: h, m, s
        
        CALL DATE_AND_TIME(TIME=timestr)

        READ(timestr(1:2),*) h
        READ(timestr(3:4),*) m
        READ(timestr(5:6),*) s

        WRITE(*,'(A,I2.2,A,I2.2,A,I2.2,A,1X,A)') '[', h, ':', m, ':', s, ']:', TRIM(text)
    END SUBROUTINE CONSOLE



END MODULE Percolation