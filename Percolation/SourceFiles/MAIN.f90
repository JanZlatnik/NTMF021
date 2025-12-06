!-------------------------------------< MAIN >-------------------------------------!
!                                                                                  !    
! Contains: Runns all modules, exports data                                        !
!                                                                                  !
! Last revision:    06/12/2025                                                     !                  
!                                                                                  !
!----------------------------------------------------------------------------------!



PROGRAM MAIN
    USE Penrose
    USE Percolation
    IMPLICIT NONE


    ! =======================================================================
    ! Computation settings & parameters
    ! =======================================================================
    LOGICAL, PARAMETER              :: test_grid_plotting = .TRUE.
    LOGICAL, PARAMETER              :: percolation_computation = .TRUE.

    INTEGER, PARAMETER              :: n_plotting_iterations = 7

    INTEGER, PARAMETER              :: n_iter_array(*) = (/ 7, 9, 11, 13 /)
    INTEGER, PARAMETER              :: n_prob_points = 200
    REAL(8), PARAMETER              :: p_crit_expected = 0.638d0
    REAL(8), PARAMETER              :: a_steep = 8.0d0
    REAL(8), PARAMETER              :: range_lim = 0.995d0

    ! =======================================================================
    ! Grid generation parameters
    ! =======================================================================
    REAL(8), PARAMETER :: xtip = 1.0d0
    REAL(8), PARAMETER :: ytip = SQRT( ((2.0d0*xtip)**2 / (phi * phi)) - 1.0d0 )

    ! =======================================================================
    
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: A, B, C, D
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: grid_points
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: boundary_flags
    INTEGER, ALLOCATABLE :: edges(:,:), num_edges(:)
    INTEGER :: i, j, id1, id2
    INTEGER :: iounit
    REAL*8, ALLOCATABLE :: p_occup_array(:)
    REAL*8 :: x, prob



    ! =======================================================================
    ! Generate Penrose grid and export to files for plotting and testing
    ! =======================================================================
    IF (test_grid_plotting) THEN

        CALL CONSOLE("======================================================================")
        CALL CONSOLE("                   STARTING PENROSE GRID TEST                         ")
        CALL CONSOLE("======================================================================")

        CALL CONSOLE("Generating Penrose grid for plotting...")

        ! Initialize Penrose tiling with two large triangle
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

        DO i = 1, n_plotting_iterations
            CALL penrose_substitution(A, B, C, D)
        END DO

        CALL CONSOLE("Penrose substitutuon finished, constructing grid...")

        CALL construct_grid(A, B, C, D, n_plotting_iterations, xtip, ytip, grid_points, edges, num_edges)

        CALL CONSOLE("Grid construction finished, generationg boundary...")

        CALL generate_boundary(grid_points, n_plotting_iterations, xtip, ytip, boundary_flags)

        CALL CONSOLE("Exporting data to files...")

        OPEN(NEWUNIT=iounit, FILE='nodes.txt', STATUS='REPLACE', ACTION='WRITE')
        DO i = 1, SIZE(grid_points, 1)
            WRITE(iounit, '(2(F20.10, 2X))') grid_points(i, 1), grid_points(i, 2)
        END DO
        CLOSE(iounit)

        OPEN(NEWUNIT=iounit, FILE='edges.txt', STATUS='REPLACE', ACTION='WRITE')
        DO id1 = 1, SIZE(edges, 1)
            DO j = 1, num_edges(id1)
                id2 = edges(id1, j)
                IF (id1 < id2) THEN
                    WRITE(iounit, '(2(I10, 2X))') id1, id2
                END IF
            END DO
        END DO
        CLOSE(iounit)

        OPEN(NEWUNIT=iounit, FILE='boundary_SW.txt', STATUS='REPLACE', ACTION='WRITE')
        DO i = 1, SIZE(grid_points, 1)
            IF (boundary_flags(i, 1)) THEN
                WRITE(iounit, '(2(F20.10, 2X))') grid_points(i, 1), grid_points(i, 2)
            END IF
        END DO
        CLOSE(iounit)

        OPEN(NEWUNIT=iounit, FILE='boundary_NE.txt', STATUS='REPLACE', ACTION='WRITE')
        DO i = 1, SIZE(grid_points, 1)
            IF (boundary_flags(i, 2)) THEN
                WRITE(iounit, '(2(F20.10, 2X))') grid_points(i, 1), grid_points(i, 2)
            END IF
        END DO
        CLOSE(iounit)
        
        DEALLOCATE(A, B, C, D, grid_points, edges, num_edges, boundary_flags)

        CALL CONSOLE("Export complete: nodes.txt, edges.txt, boundary_SW.txt, boundary_NE.txt")


        CALL CONSOLE("======================================================================")
        CALL CONSOLE("                     PENROSE GRID TEST COMPLETED                      ")
        CALL CONSOLE("======================================================================")

    END IF


    ! =======================================================================
    ! Percolation computation
    ! =======================================================================

    IF (percolation_computation) THEN

        CALL CONSOLE("======================================================================")
        CALL CONSOLE("                   STARTING PERCOLATION SIMULATION                    ")
        CALL CONSOLE("======================================================================")

        ALLOCATE(p_occup_array(n_prob_points))

        DO i = 1, n_prob_points
            x = -range_lim + REAL((i-1),8) / REAL(n_prob_points-1,8) * 2.0d0 * range_lim
            prob = p_crit_expected + ATANH(x) / a_steep
            p_occup_array(i) = MAX(0.001d0, MIN(0.999d0, prob))
        END DO

        CALL CONSOLE("Initializing computation modules...")
        CALL master_percolations(n_iter_array, p_occup_array, xtip, ytip)

        CALL CONSOLE("======================================================================")
        CALL CONSOLE("                       SIMULATION COMPLETED                           ")
        CALL CONSOLE("======================================================================")

    END IF



    ! =======================================================================
    CONTAINS

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

    
END PROGRAM MAIN