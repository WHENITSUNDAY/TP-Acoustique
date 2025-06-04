module mod_param

    implicit none

    integer, parameter  :: PR = 8
    real(PR), parameter :: PI = 4*atan(1._PR)
    real(PR), parameter :: E = 2e8_PR !Module d'Young de l'alu (Pa)
    real(PR), parameter :: nu = 0.33_PR !Coefficient de Poisson de l'alu
    real(PR), parameter :: rho_eau = 1000._PR !Masse volumique de l'eau (kg/m3)
    real(PR), parameter :: d = 0.025_PR !Diamètre du cylindre (m)
    real(PR), parameter :: rho_a = 2700._PR !Masse volumique de l'alu (kg/m3)
    real(PR), parameter :: c_eau = 1480._PR !Célérité de l'eau (m/s)
    real(PR), parameter :: c_a = 6000._PR !Célérité de l'aluminium (m/s)
    
end module mod_param