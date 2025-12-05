!--------------------------------< Random  Walk >----------------------------------!
!                                                                                  !    
! Contains: Random walks and their statisitcs                                      !
!                                                                                  !
! Last revision:    17/11/2025                                                     !                  
!                                                                                  !
!----------------------------------------------------------------------------------!



MODULE RandomWalk
    USE OMP_LIB
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: simulate_random_walks

    INTEGER, PARAMETER :: nWalks = 10000

    CONTAINS

    ! =======================================================================
    ! Random walk subroutines
    ! =======================================================================


    SUBROUTINE random_walk_simple(nsteps, edges, num_edges, start_ID, final_ID)
        IMPLICIT NONE
        INTEGER, INTENT(IN)                         :: nsteps
        INTEGER, DIMENSION(:,:), INTENT(IN)         :: edges
        INTEGER, DIMENSION(:), INTENT(IN)           :: num_edges
        INTEGER, INTENT(IN)                         :: start_ID
        INTEGER, INTENT(OUT)                        :: final_ID

        INTEGER :: current_ID, i, nedge, rand_i
        REAL(8) :: rand_val

        current_ID = start_ID

        DO i = 1, nsteps
            nedge = num_edges(current_ID)
            CALL RANDOM_NUMBER(rand_val)
            rand_i = INT(rand_val * REAL(nedge)) + 1
            current_ID = edges(current_ID, rand_i)
        END DO

        final_ID = current_ID

    END SUBROUTINE random_walk_simple



    SUBROUTINE random_walk_noback(nsteps, edges, num_edges, start_ID, final_ID, nsteps_done)
        IMPLICIT NONE
        INTEGER, INTENT(IN)                         :: nsteps
        INTEGER, DIMENSION(:,:), INTENT(IN)         :: edges
        INTEGER, DIMENSION(:), INTENT(IN)           :: num_edges
        INTEGER, INTENT(IN)                         :: start_ID
        INTEGER, INTENT(OUT)                        :: final_ID
        INTEGER, INTENT(OUT)                        :: nsteps_done

        INTEGER :: current_ID, prev_ID, nedge, rand_i, num_candidates
        INTEGER :: i, j
        INTEGER :: candidates(10)
        REAL(8) :: rand_val

        nsteps_done = 0
        current_ID = start_ID
        prev_ID = -1

        DO i = 1, nsteps
            nedge = num_edges(current_ID)
            num_candidates = 0

            DO j = 1, nedge
                IF (edges(current_ID, j) /= prev_ID) THEN
                    num_candidates = num_candidates + 1
                    candidates(num_candidates) = edges(current_ID, j)
                END IF
            END DO

            IF (num_candidates == 0) THEN
                nsteps_done = i - 1
                final_ID = current_ID
                RETURN
            ELSE 
                CALL RANDOM_NUMBER(rand_val)
                rand_i = INT(rand_val * REAL(num_candidates)) + 1
                prev_ID = current_ID
                current_ID = candidates(rand_i)
            END IF
        END DO

        nsteps_done = nsteps
        final_ID = current_ID

    END SUBROUTINE random_walk_noback


    SUBROUTINE random_walk_onlyonce(nsteps, edges, num_edges, start_ID, final_ID, nsteps_done)
        IMPLICIT NONE
        INTEGER, INTENT(IN)                         :: nsteps
        INTEGER, DIMENSION(:,:), INTENT(IN)         :: edges
        INTEGER, DIMENSION(:), INTENT(IN)           :: num_edges
        INTEGER, INTENT(IN)                         :: start_ID
        INTEGER, INTENT(OUT)                        :: final_ID
        INTEGER, INTENT(OUT)                        :: nsteps_done

        INTEGER :: current_ID, nedge, rand_i, num_candidates
        INTEGER :: i, j
        LOGICAL, ALLOCATABLE :: visited(:)
        INTEGER :: candidates(10)
        REAL(8) :: rand_val

        nsteps_done = 0
        current_ID = start_ID
        ALLOCATE(visited(SIZE(num_edges)))
        visited = .FALSE.
        visited(start_ID) = .TRUE.

        DO i = 1, nsteps
            nedge = num_edges(current_ID)
            num_candidates = 0

            DO j = 1, nedge
                IF (.NOT. visited(edges(current_ID, j))) THEN
                    num_candidates = num_candidates + 1
                    candidates(num_candidates) = edges(current_ID, j)
                END IF
            END DO

            IF (num_candidates == 0) THEN
                nsteps_done = i - 1
                final_ID = current_ID
                RETURN
            ELSE 
                CALL RANDOM_NUMBER(rand_val)
                rand_i = INT(rand_val * REAL(num_candidates)) + 1
                current_ID = candidates(rand_i)
                visited(current_ID) = .TRUE.
            END IF
        END DO

        DEALLOCATE(visited)

        nsteps_done = nsteps
        final_ID = current_ID

    END SUBROUTINE random_walk_onlyonce




    ! =======================================================================
    ! Random walk simulation subroutines for 'nWalks' walks
    ! =======================================================================

    SUBROUTINE simulate_random_walks_simple(nsteps, grid_points, nedges, num_edges, start_ID, R_values)
        IMPLICIT NONE
        INTEGER, INTENT(IN)                                  :: nsteps
        REAL(8), DIMENSION(:,:), INTENT(IN)                  :: grid_points
        INTEGER, DIMENSION(:,:), INTENT(IN)                  :: nedges
        INTEGER, DIMENSION(:), INTENT(IN)                    :: num_edges
        INTEGER, INTENT(IN)                                  :: start_ID
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT)      :: R_values

        INTEGER :: i, th_end_ID
        REAL(8) :: start_x, start_y, end_x, end_y, R

        INTEGER :: n_threads, max_seed_size, th_id
        INTEGER, ALLOCATABLE :: thread_seeds(:,:), temp_seed(:)

        start_x = grid_points(start_ID, 1)
        start_y = grid_points(start_ID, 2)

        ALLOCATE(R_values(nWalks))

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

        !$OMP PARALLEL PRIVATE(th_id, th_end_ID, end_x, end_y, R) DEFAULT(SHARED) 
        th_id = OMP_GET_THREAD_NUM() + 1
        CALL RANDOM_SEED(PUT=thread_seeds(th_id, :))

        !$OMP DO SCHEDULE(DYNAMIC)
        DO i = 1, nWalks
            CALL random_walk_simple(nsteps, nedges, num_edges, start_ID, th_end_ID)
            end_x = grid_points(th_end_ID, 1)
            end_y = grid_points(th_end_ID, 2)
            R = SQRT( (end_x - start_x)**2 + (end_y - start_y)**2 )
            R_values(i) = R
        END DO
        !$OMP END DO

        !$OMP END PARALLEL

        DEALLOCATE(thread_seeds)
        
    END SUBROUTINE simulate_random_walks_simple



    SUBROUTINE simulate_random_walks_noback(nsteps, grid_points, nedges, num_edges, start_ID, R_values, tot_attempts)
        IMPLICIT NONE
        INTEGER, INTENT(IN)                                  :: nsteps
        REAL(8), DIMENSION(:,:), INTENT(IN)                  :: grid_points
        INTEGER, DIMENSION(:,:), INTENT(IN)                  :: nedges
        INTEGER, DIMENSION(:), INTENT(IN)                    :: num_edges
        INTEGER, INTENT(IN)                                  :: start_ID
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT)      :: R_values
        INTEGER, INTENT(OUT)                                 :: tot_attempts

        INTEGER :: i, th_end_ID, th_nsteps_done
        REAL(8) :: start_x, start_y, end_x, end_y, R
        LOGICAL :: th_success

        INTEGER :: n_threads, max_seed_size, th_id
        INTEGER, ALLOCATABLE :: thread_seeds(:,:), temp_seed(:)

        start_x = grid_points(start_ID, 1)
        start_y = grid_points(start_ID, 2)
        tot_attempts = 0

        ALLOCATE(R_values(nWalks))

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

        !$OMP PARALLEL PRIVATE(th_id, th_end_ID, th_nsteps_done, th_success, end_x, end_y, R) DEFAULT(SHARED) 
        th_id = OMP_GET_THREAD_NUM() + 1
        CALL RANDOM_SEED(PUT=thread_seeds(th_id, :))

        !$OMP DO SCHEDULE(DYNAMIC) REDUCTION(+:tot_attempts) 
        DO i = 1, nWalks
            th_success = .FALSE.

            DO WHILE (.NOT. th_success)
                tot_attempts = tot_attempts + 1

                CALL random_walk_noback(nsteps, nedges, num_edges, start_ID, th_end_ID, th_nsteps_done)
                IF (th_nsteps_done == nsteps) THEN
                    th_success = .TRUE.
                END IF
            END DO

            end_x = grid_points(th_end_ID, 1)
            end_y = grid_points(th_end_ID, 2)
            R = SQRT( (end_x - start_x)**2 + (end_y - start_y)**2 )
            R_values(i) = R

        END DO
        !$OMP END DO

        !$OMP END PARALLEL

        DEALLOCATE(thread_seeds)

    END SUBROUTINE simulate_random_walks_noback



    SUBROUTINE simulate_random_walks_onlyonce(nsteps, grid_points, nedges, num_edges, start_ID, R_values, tot_attempts)
        IMPLICIT NONE
        INTEGER, INTENT(IN)                                  :: nsteps
        REAL(8), DIMENSION(:,:), INTENT(IN)                  :: grid_points
        INTEGER, DIMENSION(:,:), INTENT(IN)                  :: nedges
        INTEGER, DIMENSION(:), INTENT(IN)                    :: num_edges
        INTEGER, INTENT(IN)                                  :: start_ID
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT)      :: R_values
        INTEGER, INTENT(OUT)                                 :: tot_attempts

        INTEGER :: i, th_end_ID, th_nsteps_done
        REAL(8) :: start_x, start_y, end_x, end_y, R
        LOGICAL :: th_success

        INTEGER :: n_threads, max_seed_size, th_id
        INTEGER, ALLOCATABLE :: thread_seeds(:,:), temp_seed(:)

        start_x = grid_points(start_ID, 1)
        start_y = grid_points(start_ID, 2)
        tot_attempts = 0

        ALLOCATE(R_values(nWalks))

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

        !$OMP PARALLEL PRIVATE(th_id, th_end_ID, th_nsteps_done, th_success, end_x, end_y, R) DEFAULT(SHARED) 
        th_id = OMP_GET_THREAD_NUM() + 1
        CALL RANDOM_SEED(PUT=thread_seeds(th_id, :))

        !$OMP DO SCHEDULE(DYNAMIC) REDUCTION(+:tot_attempts) 
        DO i = 1, nWalks
            th_success = .FALSE. 

            DO WHILE (.NOT. th_success)
                tot_attempts = tot_attempts + 1

                CALL random_walk_onlyonce(nsteps, nedges, num_edges, start_ID, th_end_ID, th_nsteps_done)
                IF (th_nsteps_done == nsteps) THEN
                    th_success = .TRUE.
                END IF
            END DO

            end_x = grid_points(th_end_ID, 1)
            end_y = grid_points(th_end_ID, 2)
            R = SQRT( (end_x - start_x)**2 + (end_y - start_y)**2 )
            R_values(i) = R 

        END DO
        !$OMP END DO

        !$OMP END PARALLEL

        DEALLOCATE(thread_seeds)

    END SUBROUTINE simulate_random_walks_onlyonce



    ! =======================================================================
    ! Master simulation subroutine
    ! =======================================================================

    SUBROUTINE simulate_random_walks(nsteps_array, start_ID, grid_points, edges, num_edges)
        IMPLICIT NONE
        INTEGER, DIMENSION(:), INTENT(IN)                :: nsteps_array
        INTEGER, INTENT(IN)                              :: start_ID
        REAL(8), DIMENSION(:,:), INTENT(IN)              :: grid_points
        INTEGER, DIMENSION(:,:), INTENT(IN)              :: edges
        INTEGER, DIMENSION(:), INTENT(IN)                :: num_edges

        INTEGER :: i, j, k, nstep, iounit
        INTEGER :: total_att_nb, total_att_oo
        INTEGER :: nstepsN
        CHARACTER(LEN=80) :: filename_r, filename_stats
        CHARACTER(LEN=256) :: message

        REAL*8, ALLOCATABLE :: R_simple(:), R_noback(:), R_onlyonce(:)
        REAL*8, ALLOCATABLE :: mean_R(:,:), std_R(:,:)
        REAL(8) :: meanR, stdR

        CALL SYSTEM("mkdir -p Rvalues")

        nstepsN = SIZE(nsteps_array)
        ALLOCATE(mean_R(nstepsN,3), std_R(nstepsN,3))
        filename_stats = "Rstatistics.txt"

        DO k = 1, nstepsN
            nstep = nsteps_array(k)

            WRITE(message, '(A,I0,A)') "Simulating random walks for nsteps = ", nstep, " ..."
            CALL CONSOLE(message)        
            
            CALL CONSOLE('Computing simple random walks...')
            CALL simulate_random_walks_simple(nstep, grid_points, edges, num_edges, start_ID, R_simple)
            CALL CONSOLE('Computing simple random walks finished.')

            CALL CONSOLE('Computing no-backtracking random walks...')
            CALL simulate_random_walks_noback(nstep, grid_points, edges, num_edges, start_ID, R_noback, total_att_nb)
            CALL CONSOLE('Computing no-backtracking random walks finished.')
            WRITE(message, '(A,I0,A,I0)') 'Total attempts for no-backtracking walks: ', total_att_nb, ' out of ', nWalks
            CALL CONSOLE(message)

            CALL CONSOLE('Computing only-once random walks...')
            CALL simulate_random_walks_onlyonce(nstep, grid_points, edges, num_edges, start_ID, R_onlyonce, total_att_oo)
            CALL CONSOLE('Computing only-once random walks finished.')
            WRITE(message, '(A,I0,A,I0)') 'Total attempts for only-once walks: ', total_att_oo, ' out of ', nWalks
            CALL CONSOLE(message)

            CALL CONSOLE('Writing results to files...')
            WRITE(filename_r, '(A,I0.5,A)') 'Rvalues/R', nstep, '.txt'
            OPEN(NEWUNIT=iounit, FILE=filename_r, STATUS='REPLACE', ACTION='WRITE')
            WRITE(iounit, '(A)') '# Index, R_simple, R_noback, R_onlyonce'
            DO j = 1, nWalks
                WRITE(iounit, '(I10, 3(2X, E22.14))') j, R_simple(j), R_noback(j), R_onlyonce(j)
            END DO
            CLOSE(iounit)
            CALL CONSOLE('Writing results to files finished.')

            CALL CONSOLE('Computing R statistics...')

            meanR = SUM(R_simple) / REAL(nWalks)
            stdR = SQRT( ABS( SUM(R_simple**2)/REAL(nWalks) - meanR**2) )
            mean_R(k,1) = meanR
            std_R(k,1) = stdR

            meanR = SUM(R_noback) / REAL(nWalks)
            stdR = SQRT( ABS( SUM(R_noback**2)/REAL(nWalks) - meanR**2) )
            mean_R(k,2) = meanR
            std_R(k,2) = stdR

            meanR = SUM(R_onlyonce) / REAL(nWalks)
            stdR = SQRT( ABS( SUM(R_onlyonce**2)/REAL(nWalks) - meanR**2) )
            mean_R(k,3) = meanR
            std_R(k,3) = stdR

            DEALLOCATE(R_simple, R_noback, R_onlyonce)

             CALL CONSOLE('Computing R statistics finished.')
         
        END DO 

        CALL CONSOLE('Writing R statistics to file...')
        OPEN(NEWUNIT=iounit, FILE=filename_stats, STATUS='REPLACE', ACTION='WRITE')
        WRITE(iounit, '(A)') '# nsteps, mean_R_simple, std_R_simple, mean_R_noback, std_R_noback, mean_R_onlyonce, std_R_onlyonce'
        DO i = 1, nstepsN
            WRITE(iounit, '(I10, 6(2X, E22.14))') nsteps_array(i), mean_R(i,1), std_R(i,1), mean_R(i,2), std_R(i,2), mean_R(i,3), std_R(i,3)
        END DO
        CLOSE(iounit)
        CALL CONSOLE('Writing R statistics to file finished.')
        DEALLOCATE(mean_R, std_R)
        CALL CONSOLE('All random walk simulations finished.')

    END SUBROUTINE simulate_random_walks






    SUBROUTINE CONSOLE(message)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: message
        CHARACTER(LEN=10) :: timestr
        INTEGER :: h, m, s
        
        CALL DATE_AND_TIME(TIME=timestr)

        READ(timestr(1:2),*) h
        READ(timestr(3:4),*) m
        READ(timestr(5:6),*) s

        WRITE(*,'(A,I2.2,A,I2.2,A,I2.2,A,1X,A)') '[', h, ':', m, ':', s, ']:', TRIM(message)
    END SUBROUTINE CONSOLE

END MODULE RandomWalk