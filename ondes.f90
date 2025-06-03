program ondes

    use mod_param
    use mod_schemas
    use mod_solution

    implicit none

    real(PR) :: Lx, Ly, h, b, tmax, La, f0, u0
    integer :: Nx, Ny

    ! Paramètres physiques:
    Lx = 0.1_PR !Longueur horizontale du domaine (m)
    Ly = 0.06_PR !Longueur verticale (m)
    h = 0.04_PR !Largeur de la plaque d'alu (m)
    b = 0.005_PR !Epaisseur de la plaque d'alu (m)
    La = Lx/2-b/2 !Distance entre l'émetteur et la plaque (m)
    f0 = 2e6_PR !Fréquence de l'onde (Hz)
    u0 = 1e-9_PR !Amplitude de l'onde (m)

    !Paramètres numériques :
    Nx = 100 !Taille du maillage selon x
    Ny = 60  !Taille du maillage selon y
    tmax = 1e-3 !Temps d'étude maximal

    call solution_numerique_2D(Nx, Ny, tmax, Lx, Ly, f0, u0, frames=20)

end program ondes