program ondes

    use mod_param
    use mod_schemas
    use mod_solution

    implicit none

    real(PR) :: Lx, Ly, h, b, tmax, La, f0, u0
    integer :: Nx, Ny

    ! Paramètres physiques:
    Lx = 0.15_PR !Longueur horizontale du domaine (m)
    Ly = 0.1_PR !Longueur verticale (m)
    h = 0.04_PR !Largeur de la plaque d'alu (m)
    b = 0.005_PR !Epaisseur de la plaque d'alu (m)
    La = Lx/2-b/2 !Distance entre l'émetteur et la plaque (m)
    f0 = 2e6_PR !Fréquence de l'onde (Hz)
    u0 = 1e-9_PR !Amplitude de l'onde (m)

    !Paramètres numériques :
    Nx = 150 !Taille du maillage selon x
    Ny = 100  !Taille du maillage selon y
    tmax = 4e-4 !Temps d'étude maximal

    call solution_numerique_2D(Nx, Ny, tmax, Lx, Ly, f0, u0, cfl = 0.5_PR, frames=40)

end program ondes