module work_function
    implicit none
        integer, parameter :: mp = 8
    contains

    function read_matrix(N, file_name) result(input_matrix)
        integer :: N, i
        character(12) :: file_name
        real(mp), dimension(N, N) :: input_matrix

        open(12, file=file_name)

        do i = 1, N
            read(12, *) input_matrix(i, :)
        end do

        close(12)

    end function read_matrix

    subroutine write_matrix(num_flux, input_matrix)
        integer :: N, i, num_flux
        real(mp), dimension(0:, 0:), intent(in) :: input_matrix

        N = size(input_matrix, dim=1) - 2

        do i = 1, N
            write(num_flux, *) input_matrix(i, 1:N)
        end do

        write(num_flux, *)

    end subroutine write_matrix


    function internal_energy(T_matrix, R, M) result(omega_matrix)
        real(mp), dimension(:, :) :: T_matrix
        real(mp) :: R, M
        real(mp), dimension(1:size(T_matrix, dim=1), 1:size(T_matrix, dim=1)) :: omega_matrix

        omega_matrix = 3.0_mp / 2.0_mp * R / M * T_matrix

    end function internal_energy

    function full_energy(omega_matrix, u_matrix, v_matrix) result(eps_matrix)
        real(mp), dimension(1:, 1:) :: omega_matrix
        real(mp), dimension(0:, 0:) :: u_matrix, v_matrix
        real(mp), dimension(1:size(u_matrix, dim=1) - 2, 1:size(u_matrix, dim=1) - 2) :: eps_matrix
        integer :: N

        N = size(u_matrix, dim=1) - 2
        eps_matrix = omega_matrix + (u_matrix(1:N, 1:N) ** 2.0_mp + v_matrix(1:N, 1:N) ** 2.0_mp) / 2.0_mp

    end function full_energy

    function pressure(gamma, rho_matrix, omega_matrix) result(p_matrix)
        real(mp) :: gamma
        real(mp), dimension(0:, 0:) :: rho_matrix
        real(mp), dimension(1:, 1:) :: omega_matrix
        real(mp), dimension(1:size(rho_matrix, dim=1) - 2, 1:size(rho_matrix, dim=1) - 2) :: p_matrix
        integer :: N

        N = size(rho_matrix, dim=1) - 2

        p_matrix = (gamma - 1.0_mp) * rho_matrix(1:N, 1:N) * omega_matrix(1:N, 1:N)

    end function pressure

    subroutine boundary_conditions(u_matrix, v_matrix, eps_matrix, rho_matrix, p_matrix)
        real(mp), dimension(0:, 0:), intent(inout) :: u_matrix, v_matrix, eps_matrix, rho_matrix, p_matrix

        integer :: N, i

        N = size(u_matrix, dim=1) - 2

        do i = 1, N
            u_matrix(i, N + 1) = u_matrix(i, N)
            v_matrix(i, N + 1) = -v_matrix(i, N)
            eps_matrix(i, N + 1) = eps_matrix(i, N)
            rho_matrix(i, N + 1) = rho_matrix(i, N)
            p_matrix(i, N + 1) = p_matrix(i, N)
        end do

        do i = 1, N
            u_matrix(N + 1, i) = -u_matrix(N, i)
            v_matrix(N + 1, i) = v_matrix(N, i)
            eps_matrix(N + 1, i) = eps_matrix(N, i)
            rho_matrix(N + 1, i) = rho_matrix(N, i)
            p_matrix(N + 1, i) = p_matrix(N, i)
        end do

        do i = 1, N
            u_matrix(i, 0) = u_matrix(i, 1)
            v_matrix(i, 0) = -v_matrix(i, 1)
            eps_matrix(i, 0) = eps_matrix(i, 1)
            rho_matrix(i, 0) = rho_matrix(i, 1)
            p_matrix(i, 0) = p_matrix(i, 1)
        end do

        do i = 1, N
            u_matrix(0, i) = -u_matrix(1, i)
            v_matrix(0, i) = v_matrix(1, i)
            eps_matrix(0, i) = eps_matrix(1, i)
            rho_matrix(0, i) = rho_matrix(1, i)
            p_matrix(0, i) = p_matrix(1, i)
        end do

    end subroutine boundary_conditions

    function sgn(x)
        real(mp) :: x
        integer :: sgn

        if (x > 0.0_mp) then
            sgn = 1
        end if

        if (x < 0.0_mp) then
            sgn = -1
        end if

        if (x == 0.0_mp) then
            sgn = 0
        end if

    end function sgn

end module work_function