!-------------------------------------< MAIN >-------------------------------------!
!                                                                                  !    
! Contains: Runns all modules, exports data                                        !
!                                                                                  !
! Last revision:    05/12/2025                                                     !                  
!                                                                                  !
!----------------------------------------------------------------------------------!



PROGRAM MAIN
    USE Penrose
    IMPLICIT NONE


    ! =======================================================================
    ! Computation settings & parameters
    ! =======================================================================
    LOGICAL, PARAMETER              :: test_grid_plotting = .TRUE.
    LOGICAL, PARAMETER              :: percolation_computation = .FALSE.

    INTEGER, PARAMETER              :: n_plotting_iterations = 7

    INTEGER, PARAMETER              :: nsteps_array(*) = (/ 5, 10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250 /)

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
    CHARACTER(LEN=256) :: message



    ! =======================================================================
    ! Generate Penrose grid and export to files for plotting and testing
    ! =======================================================================
    IF (test_grid_plotting) THEN

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

    END IF


    ! =======================================================================
    ! Percolation computation
    ! =======================================================================

    IF (percolation_computation) THEN

        CALL CONSOLE("xxx...")


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