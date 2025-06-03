module mod_schemas

    use mod_param
    use mod_solution

    implicit none

    contains

        subroutine solution_numerique_2D(Nx, Ny, tmax, Lx, Ly, f0, u0, frames)
            integer, intent(in)                     :: Nx, Ny, frames
            real(PR), intent(in)                    :: tmax, Lx, Ly, f0, u0
            real(PR), dimension(0:Nx+1,0:Ny+1)      :: Pn, Pnp1, Pnm1, G
            real(PR)                                :: dx, dy, dt, tn, P0, y
            integer                                 :: ct, i, j, Ny1, Ny2
            character(20)                           :: chtn, chNx, chNy, chLx, chLy
            character(len=200)                      :: fichier

            write(chNx,'(1I4)') Nx
            write(chNy,'(1I4)') Ny
            write(chLx,'(1F5.2)') Lx*100
            write(chLy,'(1F5.2)') Ly*100 

            
            !Paramètres du schéma
            dx = Lx/(Nx+1)
            dy = Ly/(Ny+1) 
            tn = 0._PR
            ct = 0

            Ny1 = int(Ny/2 - Ny*(d/Ly)/2)  
            Ny2 = int(Ny/2 + Ny*(d/Ly)/2) 
            dt = 0.5_PR * min(dx,dy)/(c_eau*sqrt(2._PR))

            print *, "dt :", dt
            print *, tmax
            Pnp1 = 0._PR
            Pn = 0._PR
            Pnm1 = 0._PR

            P0 = u0 * (2*PI*f0) * rho_eau * c_eau
            do while (tn < tmax)

                G = calcul_laplacien(Pn, dx, dy, Nx, Ny)

                
                !Condition de Dirichlet à gauche
                do j = 1, Ny
                    y = j*dy - Ly/2._PR
                    Pnp1(0,j) = P0 * cos(2._PR*PI*f0*tn) * exp(-(y/(Ly/10._PR))**2)
                end do

                !Schéma numérique
                do j = 1, Ny
                    do i = 1, Nx
                        Pnp1(i,j) = 2._PR*Pn(i,j) - Pnm1(i,j) &
                            + c_eau**2 * dt**2 * G(i,j)
                    end do
                end do
                !CL à gauche
                Pnp1(0,Ny2:Ny+1) = Pn(1,Ny2:Ny+1)
                Pnp1(0,0:Ny1) = Pn(1,0:Ny1)

                !CL en bas/haut 
                Pnp1(:,Ny+1) = Pnp1(:,Ny)
                Pnp1(:,0) = Pnp1(:,1) 

                !CL à droite
                Pnp1(Nx+1,Ny2:Ny+1) = Pnp1(Nx,Ny2:Ny+1)
                Pnp1(Nx+1,0:Ny1) = Pnp1(Nx,0:Ny1)
                Pnp1(Nx+1,Ny1:Ny2) = 0._PR
                
                Pnm1 = Pn
                Pn = Pnp1

                if (mod(ct, frames) == 0) then
                    write (chtn,'(1F10.7)') tn

                    fichier = 'data/P_Lx_'// trim(adjustl(chLx)) &
                    & // '_Ly_'// trim(adjustl(chLy)) &
                    & // '_Nx_' // trim(adjustl(chNX)) &
                    & // '_Ny_' // trim(adjustl(chNy)) &
                    & // '_tn_' // trim(adjustl(chtn)) //'.dat'
                
                    open(unit=100, file=fichier)

                    do i = 0, Nx+1
                        do j = 0, Ny+1
                            write(100, *) i*dx, j*dy, Pn(i,j)
                        end do
                    end do
                end if
                    

                ct = ct + 1
                tn = tn + dt
            end do

            close(100)
        end subroutine solution_numerique_2D

        function calcul_F(dx, dy, x0, y0, Nx, Ny) result(F) 

            real(PR) :: dx, dy, x0, y0, xi, yj
            real(PR), dimension(Nx,Ny) :: F
            integer :: Nx, Ny, i, j
            
            F = 0._PR

            do i = 1, Nx
                xi = x0 + (i-1)*dx
                do j = 1, Ny
                    yj = y0 + (j-1)*dy
                    F(i,j) = exp(-7*(xi**2 + yj**2))
                end do
            end do
            
        end function calcul_F

        function calcul_laplacien(U, dx, dy, Nx, Ny) result(G)
            implicit none
            integer :: i, j
            real(PR), dimension(0:Nx+1,0:Ny+1) :: U, G
            real(PR) :: dx, dy
            integer :: Nx, Ny

            G = 0._PR

            do i = 1, Nx
                do j = 1, Ny
                    G(i,j) = (U(i+1,j) - 2._PR*U(i,j) + U(i-1,j))/(dx**2) + &
                             (U(i,j+1) - 2._PR*U(i,j) + U(i,j-1))/(dy**2)
                end do
            end do

        end function calcul_laplacien


        subroutine solution_numerique_2D_alu(Nx, Ny, tmax, Lx, Ly, d, h, b, La, frames)
            integer, intent(in)                     :: Nx, Ny, frames
            real(PR), intent(in)                    :: tmax, Lx, Ly, d, h, b, La
            real(PR), dimension(0:Nx+1,0:Ny+1)      :: Uxn, Uxnp1, Uxnm1, Uyn, Uynp1, Uynm1, Gx, Gy 
            real(PR)                                :: dx, dy, dt, tn
            integer                                 :: ct, i, j
            character(20)                           :: chtn, chNx, chNy, chLx, chLy, chd, chb, ch
            character(len=200)                      :: fichier

            write(chNx,'(1I4)') Nx
            write(chNy,'(1I4)') Ny
            write(chLx,'(1F5.2)') Lx
            write(chLy,'(1F5.2)') Ly 
            write(chLy,'(1F5.2)') h 
            write(chLy,'(1F5.2)') Ly 


            !Paramètres du schéma
            dx = Lx/(Nx+1)
            dy = Ly/(Ny+1) 
            tn = 0._PR
            ct = 0

            dt = min(dx,dy)/(c_eau*sqrt(2._PR))
            Uxnm1 = 0._PR
            Uxn = 0._PR
            Uxnp1 = 0._PR
            Uynm1 = 0._PR
            Uyn = 0._PR
            Uynp1 = 0._PR

            do while (tn < tmax)

                Uxnp1(1:Nx,1:Ny) = 2._PR*Uxn(1:Nx,1:Ny) - Uxnm1(1:Nx,1:Ny) &
                    + c_eau**2 * dt**2 * Gx(1:Nx,1:Ny)

                Uynp1(1:Nx,1:Ny) = 2._PR*Uyn(1:Nx,1:Ny) - Uynm1(1:Nx,1:Ny) &
                    + c_eau**2 * dt**2 * Gy(1:Nx,1:Ny)

                Uxnm1(1:Nx,1:Ny) = Uxn(1:Nx,1:Ny)
                Uxn(1:Nx,1:Ny) = Uxnp1(1:Nx,1:Ny)

                Uynm1(1:Nx,1:Ny) = Uyn(1:Nx,1:Ny)
                Uyn(1:Nx,1:Ny) = Uynp1(1:Nx,1:Ny)

                if (mod(ct, frames) == 0) then
                    write (chtn,'(1F5.2)') tn
                end if
                    

                ct = ct + 1
                tn = tn + dt
            end do

            close(100)
        end subroutine solution_numerique_2D_alu

end module mod_schemas