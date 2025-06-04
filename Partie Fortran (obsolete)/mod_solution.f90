module mod_solution

    use mod_param

    implicit none

    contains

        function F(x, y) 

            real(PR) :: x, y, F
            
            F = exp(-7*(x**2 + y**2))
        end function F


        function calcul_h(t, f0) result(h)

            real(PR) :: h, t, f0

            h = f0 * (2._PR * pi**2 * (f0*t-1._PR)**2 - 1) *&
             exp(-pi**2 * (f0*t-1._PR)**2)

        end function calcul_h
end module mod_solution