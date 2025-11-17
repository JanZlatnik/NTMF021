!-----------------------------------< Penrose >------------------------------------!
!                                                                                  !    
! Contains: Create Penrose tiling by substitution                                  !
!                                                                                  !
! Last revision:    17/11/2025                                                     !                  
!                                                                                  !
!----------------------------------------------------------------------------------!


MODULE Penrose
    IMPLICIT NONE
    PRIVATE 
    PUBLIC :: penrose_substitution, construct_grid, phi, iphi, ophi

    REAL(8), PARAMETER :: phi = (1.0d0 + SQRT(5.0d0)) / 2.0d0 
    REAL(8), PARAMETER :: iphi = 1.0d0 / phi
    REAL(8), PARAMETER :: ophi = 1.0d0 - iphi


    CONTAINS

    SUBROUTINE penrose_substitution(A,B,C,D)
        REAL*8, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: A, B, C, D

        REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: A_new, B_new, C_new, D_new
        INTEGER :: nA, nB, nC, nD
        INTEGER :: nA_new, nB_new, nC_new, nD_new
        INTEGER :: i
        REAL(8) :: x1, y1, x2, y2

        IF (ALLOCATED(A)) THEN
            nA = SIZE(A, 1)
        ELSE
            nA = 0
        END IF

        IF (ALLOCATED(B)) THEN
            nB = SIZE(B, 1)
        ELSE
            nB = 0
        END IF

        IF (ALLOCATED(C)) THEN
            nC = SIZE(C, 1)
        ELSE
            nC = 0
        END IF

        IF (ALLOCATED(D)) THEN
            nD = SIZE(D, 1)
        ELSE
            nD = 0
        END IF

        nA_new = nA + nC
        nB_new = nB + nD
        nC_new = nB + nC + nD
        nD_new = nA + nC + nD

        ALLOCATE(A_new(nA_new, 3, 2))
        ALLOCATE(B_new(nB_new, 3, 2))
        ALLOCATE(C_new(nC_new, 3, 2))
        ALLOCATE(D_new(nD_new, 3, 2))

        DO i = 1, nA
            x1 = ophi * A(i,1,1) + iphi * A(i,2,1)
            y1 = ophi * A(i,1,2) + iphi * A(i,2,2)
            A_new(i,1,:) = A(i,3,:)
            A_new(i,2,1) = x1
            A_new(i,2,2) = y1
            A_new(i,3,:) = A(i,2,:)
            D_new(i,1,1) = x1
            D_new(i,1,2) = y1
            D_new(i,2,:) = A(i,3,:)
            D_new(i,3,:) = A(i,1,:)
        END DO

        DO i = 1, nB
            x1 = ophi * B(i,1,1) + iphi * B(i,3,1)
            y1 = ophi * B(i,1,2) + iphi * B(i,3,2)
            B_new(i,1,:) = B(i,2,:)
            B_new(i,2,:) = B(i,3,:)
            B_new(i,3,1) = x1
            B_new(i,3,2) = y1
            C_new(i,1,1) = x1
            C_new(i,1,2) = y1
            C_new(i,2,:) = B(i,1,:)
            C_new(i,3,:) = B(i,2,:) 
        END DO

        DO i = 1, nC
            x1 = ophi * C(i,3,1) + iphi * C(i,2,1)
            y1 = ophi * C(i,3,2) + iphi * C(i,2,2)
            x2 = ophi * C(i,3,1) + iphi * C(i,1,1)
            y2 = ophi * C(i,3,2) + iphi * C(i,1,2)
            A_new(nA + i,1,1) = x1
            A_new(nA + i,1,2) = y1
            A_new(nA + i,2,1) = x2
            A_new(nA + i,2,2) = y2
            A_new(nA + i,3,:) = C(i,1,:)
            C_new(nB + i,1,1) = x1
            C_new(nB + i,1,2) = y1
            C_new(nB + i,2,:) = C(i,1,:)
            C_new(nB + i,3,:) = C(i,2,:)
            D_new(nA + i,1,1) = x2
            D_new(nA + i,1,2) = y2
            D_new(nA + i,2,1) = x1
            D_new(nA + i,2,2) = y1
            D_new(nA + i,3,:) = C(i,3,:)
        END DO

        DO i = 1, nD
            x1 = ophi * D(i,2,1) + iphi * D(i,3,1)
            y1 = ophi * D(i,2,2) + iphi * D(i,3,2)
            x2 = ophi * D(i,2,1) + iphi * D(i,1,1)
            y2 = ophi * D(i,2,2) + iphi * D(i,1,2)
            B_new(nB + i,1,1) = x1
            B_new(nB + i,1,2) = y1
            B_new(nB + i,2,:) = D(i,1,:)
            B_new(nB + i,3,1) = x2
            B_new(nB + i,3,2) = y2
            C_new(nB + nC + i,1,1) = x2
            C_new(nB + nC + i,1,2) = y2
            C_new(nB + nC + i,2,:) = D(i,2,:)
            C_new(nB + nC + i,3,1) = x1
            C_new(nB + nC + i,3,2) = y1
            D_new(nA + nC + i,1,1) = x1
            D_new(nA + nC + i,1,2) = y1
            D_new(nA + nC + i,2,:) = D(i,3,:)
            D_new(nA + nC + i,3,:) = D(i,1,:)
        END DO


        IF(ALLOCATED(A)) DEALLOCATE(A)
        IF(ALLOCATED(B)) DEALLOCATE(B)
        IF(ALLOCATED(C)) DEALLOCATE(C)
        IF(ALLOCATED(D)) DEALLOCATE(D)
        
        ALLOCATE(A(nA_new, 3, 2))
        A = A_new
        DEALLOCATE(A_new)
        ALLOCATE(B(nB_new, 3, 2))
        B = B_new
        DEALLOCATE(B_new)
        ALLOCATE(C(nC_new, 3, 2))
        C = C_new
        DEALLOCATE(C_new)
        ALLOCATE(D(nD_new, 3, 2))
        D = D_new
        DEALLOCATE(D_new)

    END SUBROUTINE penrose_substitution






    SUBROUTINE construct_grid(A,B,C,D,n,xtip,ytip, nodes,neighbors,num_neighbors)
        REAL*8, DIMENSION(:,:,:), INTENT(IN)    :: A, B, C, D
        INTEGER, INTENT(IN)                     :: n
        REAL(8), INTENT(IN)                     :: xtip, ytip

        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: nodes 
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: neighbors
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)   :: num_neighbors   

        INTEGER, PARAMETER :: max_neighbors = 10

        REAL(8) :: astand, amin, grid_cell_size, eps, x, y
        REAL(8) :: x_min, x_max, y_min, y_max
        INTEGER :: xN, yN
        INTEGER, ALLOCATABLE :: hashgrid(:,:)

        REAL*8, ALLOCATABLE :: nodes_temp(:,:)
        INTEGER :: num_nodes

        INTEGER, DIMENSION(:,:), ALLOCATABLE :: Anodes(:,:), Bnodes(:,:), Cnodes(:,:), Dnodes(:,:)
        INTEGER :: nA, nB, nC, nD
        INTEGER :: i, j, ID

        astand = SQRT(xtip**2 + ytip**2) / phi**n
        amin = astand * iphi
        grid_cell_size = amin / 5.0d0
        eps = amin / 100.0d0

        x_max = xtip + amin
        x_min = -x_max
        y_max = ytip + amin
        y_min = -y_max

        xN = CEILING((x_max - x_min) / grid_cell_size) + 1
        yN = CEILING((y_max - y_min) / grid_cell_size) + 1

        nA = SIZE(A, 1)
        nB = SIZE(B, 1)
        nC = SIZE(C, 1)
        nD = SIZE(D, 1) 

        ALLOCATE(hashgrid(xN, yN))
        hashgrid = -1

        ALLOCATE(nodes_temp(3*(nA + nB + nC + nD), 2))
        num_nodes = 0


        ALLOCATE(Anodes(nA, 3), Bnodes(nB, 3), Cnodes(nC, 3), Dnodes(nD, 3))


        DO i = 1, nA
            DO j = 1, 3
                x = A(i,j,1)
                y = A(i,j,2)
                ID = find_or_add_node(x, y)
                Anodes(i,j) = ID
            END DO
        END DO

        DO i = 1, nB
            DO j = 1, 3
                x = B(i,j,1)
                y = B(i,j,2)
                ID = find_or_add_node(x, y)
                Bnodes(i,j) = ID
            END DO
        END DO

        DO i = 1, nC
            DO j = 1, 3
                x = C(i,j,1)
                y = C(i,j,2)
                ID = find_or_add_node(x, y)
                Cnodes(i,j) = ID
            END DO
        END DO

        DO i = 1, nD
            DO j = 1, 3
                x = D(i,j,1)
                y = D(i,j,2)
                ID = find_or_add_node(x, y)
                Dnodes(i,j) = ID
            END DO
        END DO

        ALLOCATE(nodes(num_nodes, 2))
        nodes = nodes_temp(1:num_nodes, :)
        DEALLOCATE(nodes_temp)

        ALLOCATE(neighbors(num_nodes, max_neighbors))
        ALLOCATE(num_neighbors(num_nodes))
        neighbors = -1
        num_neighbors = 0

        DO i = 1, nA
            CALL add_neighbour(Anodes(i,1), Anodes(i,2))
            CALL add_neighbour(Anodes(i,1), Anodes(i,3))
        END DO

        DO i = 1, nB
            CALL add_neighbour(Bnodes(i,1), Bnodes(i,2))
            CALL add_neighbour(Bnodes(i,1), Bnodes(i,3))
        END DO

        DO i = 1, nC
            CALL add_neighbour(Cnodes(i,1), Cnodes(i,2))
            CALL add_neighbour(Cnodes(i,1), Cnodes(i,3))
        END DO

        DO i = 1, nD
            CALL add_neighbour(Dnodes(i,1), Dnodes(i,2))
            CALL add_neighbour(Dnodes(i,1), Dnodes(i,3))
        END DO

        DEALLOCATE(Anodes, Bnodes, Cnodes, Dnodes)
        DEALLOCATE(hashgrid)

    CONTAINS

        SUBROUTINE get_grid_indices(xx,yy,ixx,iyy)
            REAL*8, INTENT(IN) :: xx, yy
            INTEGER, INTENT(OUT) :: ixx, iyy

            ixx = INT(FLOOR((xx - x_min) / grid_cell_size)) + 1
            iyy = INT(FLOOR((yy - y_min) / grid_cell_size)) + 1

            IF (ixx < 1) ixx = 1
            IF (iyy < 1) iyy = 1
            IF (ixx > xN) ixx = xN
            IF (iyy > yN) iyy = yN

        END SUBROUTINE get_grid_indices



        FUNCTION find_or_add_node(xx,yy) RESULT(node_id)
            REAL*8, INTENT(IN) :: xx, yy
            INTEGER :: node_id, candidate_id
            INTEGER :: ixx, iyy
            INTEGER :: in, jn
            REAL(8) :: dist

            CALL get_grid_indices(xx, yy, ixx, iyy)

            DO in = MAX(1, ixx-1), MIN(xN, ixx+1)
                DO jn = MAX(1, iyy-1), MIN(yN, iyy+1)
                    candidate_id = hashgrid(in, jn)
                    IF (candidate_id > 0) THEN
                        dist = SQRT((nodes_temp(candidate_id,1) - xx)**2 + (nodes_temp(candidate_id,2) - yy)**2)
                        IF (dist < eps) THEN
                            node_id = candidate_id
                            RETURN
                        ELSE
                            PRINT*, "Distance check failed: ", dist, ' compared to amin = ', amin
                        END IF
                    END IF
                END DO
            END DO

            num_nodes = num_nodes + 1
            node_id = num_nodes
            nodes_temp(node_id,1) = xx
            nodes_temp(node_id,2) = yy
            hashgrid(ixx, iyy) = node_id

        END FUNCTION find_or_add_node


        SUBROUTINE add_neighbour(id1, id2)
            INTEGER, INTENT(IN) :: id1, id2
            INTEGER :: k, kn
            LOGICAL :: exists

            exists = .FALSE.
            kn = num_neighbors(id1)
            DO k = 1, kn
                IF (neighbors(id1, k) == id2) THEN
                    exists = .TRUE.
                    EXIT
                END IF
            END DO

            IF (.NOT. exists) THEN
                num_neighbors(id1) = num_neighbors(id1) + 1
                neighbors(id1, num_neighbors(id1)) = id2

                num_neighbors(id2) = num_neighbors(id2) + 1
                neighbors(id2, num_neighbors(id2)) = id1
            END IF

        END SUBROUTINE add_neighbour





    END SUBROUTINE construct_grid




END MODULE Penrose