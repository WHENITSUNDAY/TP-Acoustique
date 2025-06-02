program ondes

    use mod_param
    use mod_schemas
    use mod_solution

    implicit none

    real(PR) :: Lx, Ly, d, h, b, tmax, La
    integer :: Nx, Ny, frames

    ! Paramètres physiques:
    Lx = 0.1_PR !Longueur horizontale (m)
    Ly = 0.06_PR !Longueur verticale (m)
    d = 0.02_PR !Diamètre du capteur piezométrique (m)
    h = 0.04_PR !Largeur de la plaque d'alu (m)
    b = 0.005_PR !Epaisseur de la plaque d'alu (m)
    La = Lx/2-b/2

    !Paramètres numériques :
    Nx = 100 !Taille du maillage selon x
    Ny = 60  !Taille du maillage selon y
    tmax = 1e-5 !Temps d'étude maximal
    frames = 100 !Nombre d'images qu'on va afficher

    call solution_numerique_2D(Nx, Ny, tmax, Lx, Ly, d, h, b, La, frames)

end program ondes