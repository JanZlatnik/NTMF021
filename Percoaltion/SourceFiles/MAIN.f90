!-------------------------------------< MAIN >-------------------------------------!
!                                                                                  !    
! Contains: Runns all modules, exports data                                        !
!                                                                                  !
! Last revision:    17/11/2025                                                     !                  
!                                                                                  !
!----------------------------------------------------------------------------------!



PROGRAM MAIN
    USE Penrose
    USE RandomWalk
    IMPLICIT NONE


    ! =======================================================================
    ! Computation settings & parameters
    ! =======================================================================
    LOGICAL, PARAMETER              :: test_grid_plotting = .TRUE.
    LOGICAL, PARAMETER              :: regenerate_penrose_grid = .TRUE. 
    LOGICAL, PARAMETER              :: random_walk_computation = .TRUE.

    INTEGER, PARAMETER              :: n_plotting_iterations = 9
    INTEGER, PARAMETER              :: n_penrose_iterations = 14

    INTEGER, PARAMETER              :: nsteps_array(*) = (/ 5, 10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250 /)

    ! =======================================================================
    ! Grid generation parameters
    ! =======================================================================
    REAL(8), PARAMETER :: xtip = 1.0d0
    REAL(8), PARAMETER :: ytip = SQRT( ((2.0d0*xtip)**2 / (phi * phi)) - 1.0d0 )

    ! =======================================================================
    
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: A, B, C, D
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: grid_points
    REAL(8) :: scaling_factor
    INTEGER, ALLOCATABLE :: edges(:,:), num_edges(:)
    INTEGER :: i, j, id1, id2
    INTEGER :: iounit
    INTEGER :: dim1, dim2
    INTEGER :: start_id
    REAL(8) :: min_dist, dist
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

        CALL CONSOLE("Grid construction finished, exporting grid data to files...")


        scaling_factor = (phi ** n_plotting_iterations) / SQRT(xtip**2 + ytip**2)
        grid_points = grid_points * scaling_factor


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
        
        DEALLOCATE(A, B, C, D, grid_points, edges, num_edges)

        CALL CONSOLE("Grid data exported to files: penrose.txt, nodes.txt, edges.txt")

    END IF


    ! =======================================================================
    ! Generation of Penrose grid for computations
    ! =======================================================================
    IF (regenerate_penrose_grid) THEN

        CALL CONSOLE("Regenerating Penrose grid for computations...")

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

        DO i = 1, n_penrose_iterations
            CALL penrose_substitution(A, B, C, D)
        END DO

        CALL CONSOLE("Penrose substitutuon finished, constructing grid...")

        CALL construct_grid(A, B, C, D, n_penrose_iterations, xtip, ytip, grid_points, edges, num_edges)

        CALL CONSOLE("Grid construction finished, exporting grid data to files...")

        dim1 = SIZE(grid_points, 1)
        dim2 = SIZE(grid_points, 2)
        OPEN(NEWUNIT=iounit, FILE='grid_points.bin', FORM='UNFORMATTED', STATUS='REPLACE')
        WRITE(iounit) dim1, dim2  
        WRITE(iounit) grid_points 
        CLOSE(iounit)

        dim1 = SIZE(edges, 1)
        dim2 = SIZE(edges, 2)
        OPEN(NEWUNIT=iounit, FILE='edges.bin', FORM='UNFORMATTED', STATUS='REPLACE')
        WRITE(iounit) dim1, dim2  
        WRITE(iounit) edges 
        CLOSE(iounit)

        dim1 = SIZE(num_edges)
        OPEN(NEWUNIT=iounit, FILE='num_edges.bin', FORM='UNFORMATTED', STATUS='REPLACE')
        WRITE(iounit) dim1  
        WRITE(iounit) num_edges 
        CLOSE(iounit)

        CALL CONSOLE("Grid data exported to binary files: grid_points.bin, edges.bin, num_edges.bin")

        DEALLOCATE(A, B, C, D, grid_points, edges, num_edges)

    END IF


    ! =======================================================================
    ! Generation of Penrose grid for computations
    ! =======================================================================

    IF (random_walk_computation) THEN

        CALL CONSOLE("Loading Penrose grid...")

        OPEN(NEWUNIT=iounit, FILE='grid_points.bin', FORM='UNFORMATTED', STATUS='OLD')
        READ(iounit) dim1, dim2
        ALLOCATE(grid_points(dim1, dim2))
        READ(iounit) grid_points
        CLOSE(iounit)

        OPEN(NEWUNIT=iounit, FILE='edges.bin', FORM='UNFORMATTED', STATUS='OLD')
        READ(iounit) dim1, dim2
        ALLOCATE(edges(dim1, dim2))
        READ(iounit) edges
        CLOSE(iounit)

        OPEN(NEWUNIT=iounit, FILE='num_edges.bin', FORM='UNFORMATTED', STATUS='OLD')
        READ(iounit) dim1
        ALLOCATE(num_edges(dim1))
        READ(iounit) num_edges
        CLOSE(iounit)

        CALL CONSOLE("Penrose grid succesfully loaded.")

        scaling_factor = (phi ** n_penrose_iterations) / SQRT(xtip**2 + ytip**2)
        grid_points = grid_points * scaling_factor
        CALL CONSOLE("Grid rescaled to d=1.")

        CALL CONSOLE("Searching for starting point...")
        min_dist = 1.0d100
        start_id = -1
        DO i = 1, SIZE(grid_points, 1)
            dist = (grid_points(i, 1))**2 + (grid_points(i, 2))**2 
            IF (dist < min_dist) THEN
                min_dist = dist
                start_id = i
            END IF
        END DO
        WRITE(message, '(A,I0,A,F0.6,A,F0.6,A)') "Starting point ID found: ", start_id, " at (", grid_points(start_id, 1), ", ", grid_points(start_id, 2), ")"
        CALL CONSOLE(message)

        CALL CONSOLE("Simulating random walks...")
        CALL simulate_random_walks(nsteps_array, start_id, grid_points, edges, num_edges)
        CALL CONSOLE("Random walk simulation finished.")

        DEALLOCATE(grid_points, edges, num_edges)


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