module work_program
        use :: work_function
        implicit none
    contains

    subroutine method_Davydov_Belotserkovsky(gamma, g, R, M, h, tau, u_matrix, v_matrix, T_matrix, rho_matrix)
        real(mp), intent(in) :: gamma, g, R, M, h, tau
        real(mp), intent(inout), dimension(0:, 0:) :: u_matrix, v_matrix, rho_matrix
        real(mp), intent(inout), dimension(1:, 1:) :: T_matrix

        real(mp), dimension(0:size(u_matrix, dim=1) - 1, 0:size(u_matrix, dim=1) - 1) :: &
        & u_tilde_matrix, v_tilde_matrix, eps_tilde_matrix

        real(mp), dimension(1:size(rho_matrix, dim=1) - 2, 1:size(rho_matrix, dim=1) - 2) :: omega_matrix

        integer, dimension(1:size(rho_matrix, dim=1) - 2, 1:size(rho_matrix, dim=1) - 2) :: s_r, s_u
        real(mp), dimension(0:size(rho_matrix, dim=1) - 2, 0:size(rho_matrix, dim=1) - 2) :: Pi_r, Pi_u
        real(mp), dimension(0:size(rho_matrix, dim=1) - 1, 0:size(rho_matrix, dim=1) - 1) :: rho_matrix_new

        real(mp), dimension(0:size(rho_matrix, dim=1) - 1, 0:size(rho_matrix, dim=1) - 1) :: u_matrix_new
        real(mp), dimension(0:size(rho_matrix, dim=1) - 1, 0:size(rho_matrix, dim=1) - 1) :: v_matrix_new
        real(mp), dimension(0:size(rho_matrix, dim=1) - 1, 0:size(rho_matrix, dim=1) - 1) :: eps_matrix_new

        real(mp) :: t1, t2

        integer :: N, k

        N = size(u_matrix, dim=1) - 2

        open(12, file='rho_result_matrix.dat', status='replace')
        open(13, file='u_result_matrix.dat', status='replace')
        open(14, file='v_result_matrix.dat', status='replace')
        open(15, file='internal_energy_result_matrix.dat', status='replace')

        call cpu_time(t1)

        do k = 0, 1000000
            if (mod(k, 1000) == 0) then
                call cpu_time(t2)
                omega_matrix = internal_energy(T_matrix, R, M)
                write(*, *) t2 - t1, k

                call write_matrix(12, rho_matrix)

                call write_matrix(13, u_matrix)

                call write_matrix(14, v_matrix)

                call write_matrix(15, omega_matrix)
                
                call cpu_time(t1)
            end if

            call Eulerian_stage(gamma, g, R, M, h, tau, u_matrix, v_matrix, T_matrix, rho_matrix, &
            & u_tilde_matrix, v_tilde_matrix, eps_tilde_matrix)

            call mass_transport(h, tau, rho_matrix, Pi_r, Pi_u, s_r, s_u, u_tilde_matrix, v_tilde_matrix, rho_matrix_new)

            call gas_param_transport(h, tau, rho_matrix, Pi_r, Pi_u, s_r, s_u, rho_matrix_new, u_tilde_matrix, &
            & u_matrix_new)
            call gas_param_transport(h, tau, rho_matrix, Pi_r, Pi_u, s_r, s_u, rho_matrix_new, v_tilde_matrix, &
            & v_matrix_new)
            call gas_param_transport(h, tau, rho_matrix, Pi_r, Pi_u, s_r, s_u, rho_matrix_new, eps_tilde_matrix, &
            & eps_matrix_new)

            rho_matrix = rho_matrix_new
            u_matrix = u_matrix_new
            v_matrix = v_matrix_new
            T_matrix = 2.0_mp / 3.0_mp * M / R * (eps_matrix_new(1:N, 1:N) - &
            & (u_matrix(1:N, 1:N) ** 2.0_mp + v_matrix(1:N, 1:N) ** 2.0_mp) / 2.0_mp)

        end do

        close(12)
        close(13)
        close(14)
        close(15)

    end subroutine method_Davydov_Belotserkovsky

    subroutine Eulerian_stage(gamma, g, R, M, h, tau, u_matrix, v_matrix, T_matrix, rho_matrix, &
        & u_tilde_matrix, v_tilde_matrix, eps_tilde_matrix)
        real(mp), intent(in) :: gamma, g, h, tau, R, M
        real(mp), intent(inout), dimension(0:, 0:) :: u_matrix, v_matrix, rho_matrix
        real(mp), intent(in), dimension(1:, 1:) :: T_matrix

        real(mp), intent(out), dimension(0:size(u_matrix, dim=1) - 1, 0:size(u_matrix, dim=1) - 1) :: &
        & u_tilde_matrix, v_tilde_matrix, eps_tilde_matrix

        real(mp), dimension(0:size(u_matrix, dim=1) - 1, 0:size(u_matrix, dim=1) - 1) :: p_matrix, eps_matrix
        real(mp), dimension(1:size(u_matrix, dim=1) - 2, 1:size(u_matrix, dim=1) - 2) :: omega_matrix
        real(mp) :: p_bound_l, p_bound_r, p_bound_u, p_bound_d
        real(mp) :: u_bound_l, u_bound_r, v_bound_u, v_bound_d

        integer :: N, i, j

        N = size(u_matrix, dim=1) - 2

        omega_matrix = internal_energy(T_matrix, R, M)
        eps_matrix(1:N, 1:N) = full_energy(omega_matrix, u_matrix, v_matrix)
        p_matrix(1:N, 1:N) = pressure(gamma, rho_matrix, omega_matrix)

        call boundary_conditions(u_matrix, v_matrix, eps_matrix, rho_matrix, p_matrix)

        do i = 1, N
            do j = 1, N
                p_bound_r = (p_matrix(i, j) + p_matrix(i + 1, j)) / 2.0_mp
                p_bound_l = (p_matrix(i, j) + p_matrix(i - 1, j)) / 2.0_mp
                p_bound_u = (p_matrix(i, j) + p_matrix(i, j + 1)) / 2.0_mp
                p_bound_d = (p_matrix(i, j) + p_matrix(i, j - 1)) / 2.0_mp

                u_bound_r = (u_matrix(i, j) + u_matrix(i + 1, j)) / 2.0_mp
                u_bound_l = (u_matrix(i, j) + u_matrix(i - 1, j)) / 2.0_mp

                v_bound_u = (v_matrix(i, j) + v_matrix(i, j + 1)) / 2.0_mp
                v_bound_d = (v_matrix(i, j) + v_matrix(i, j - 1)) / 2.0_mp

                u_tilde_matrix(i, j) = u_matrix(i, j) - (p_bound_r - p_bound_l) / rho_matrix(i, j) * tau / h
                v_tilde_matrix(i, j) = v_matrix(i, j) - (p_bound_u - p_bound_d) / rho_matrix(i, j) * tau / h &
                & - g * tau
                eps_tilde_matrix(i, j) = eps_matrix(i, j) - ((p_bound_r * u_bound_r - p_bound_l * u_bound_l) / rho_matrix(i, j) + &
                & (p_bound_u * v_bound_u - p_bound_d * v_bound_d) / rho_matrix(i, j)) * tau / h

            end do
        end do

        call boundary_conditions(u_tilde_matrix, v_tilde_matrix, eps_tilde_matrix, rho_matrix, p_matrix)

    end subroutine Eulerian_stage

    subroutine mass_transport(h, tau, rho_matrix, Pi_r, Pi_u, s_r, s_u, u_tilde_matrix, v_tilde_matrix, rho_matrix_new)
        real(mp), intent(in) :: h, tau
        real(mp), intent(in), dimension(0:, 0:) :: rho_matrix
        real(mp), intent(out), dimension(0:size(rho_matrix, dim=1) - 2, 0:size(rho_matrix, dim=1) - 2) :: Pi_r, Pi_u
        real(mp), intent(out), dimension(0:size(rho_matrix, dim=1) - 1, 0:size(rho_matrix, dim=1) - 1) :: u_tilde_matrix, &
        & v_tilde_matrix, rho_matrix_new

        real(mp), dimension(1:size(rho_matrix, dim=1) - 2, 1:size(rho_matrix, dim=1) - 2) :: u_tilde_matrix_r, v_tilde_matrix_u
        real(mp), dimension(1:size(rho_matrix, dim=1) - 2, 1:size(rho_matrix, dim=1) - 2) :: rho_tilde_matrix_r, rho_tilde_matrix_u
        integer, intent(out), dimension(1:size(rho_matrix, dim=1) - 2, 1:size(rho_matrix, dim=1) - 2) :: s_r, s_u


        integer :: N, i, j

        N = size(rho_matrix, dim=1) - 2

        do i = 1, N
            do j = 1, N

                s_r(i, j) = sgn(u_tilde_matrix(i, j) + u_tilde_matrix(i + 1, j))
                s_u(i, j) = sgn(v_tilde_matrix(i, j) + v_tilde_matrix(i, j + 1))
            
            end do
        end do

        Pi_r = 0.0_mp
        Pi_u = 0.0_mp

        u_tilde_matrix_r = 0.0_mp
        v_tilde_matrix_u = 0.0_mp

        do i = 1, N
            do j = 1, N

                if (s_r(i, j) == 1) then
                    u_tilde_matrix_r(i, j) = u_tilde_matrix(i, j) & 
                    & + (u_tilde_matrix(i + 1, j) - u_tilde_matrix(i - 1, j)) / 4.0_mp
                    rho_tilde_matrix_r(i, j) = rho_matrix(i, j) & 
                    & + (rho_matrix(i + 1, j) - rho_matrix(i - 1, j)) / 4.0_mp
                end if

                if (s_r(i, j) == -1) then
                    u_tilde_matrix_r(i, j) = u_tilde_matrix(i + 1, j) & 
                    & - (u_tilde_matrix(i + 2, j) - u_tilde_matrix(i, j)) / 4.0_mp
                    rho_tilde_matrix_r(i, j) = rho_matrix(i + 1, j) & 
                    & - (rho_matrix(i + 2, j) - rho_matrix(i, j)) / 4.0_mp
                end if

                if (s_u(i, j) == 1) then
                    v_tilde_matrix_u(i, j) = v_tilde_matrix(i, j) & 
                    & + (v_tilde_matrix(i, j + 1) - v_tilde_matrix(i, j - 1)) / 4.0_mp
                    rho_tilde_matrix_u(i, j) = rho_matrix(i, j) & 
                    & + (rho_matrix(i, j + 1) - rho_matrix(i, j - 1)) / 4.0_mp
                end if

                if (s_u(i, j) == -1) then
                    v_tilde_matrix_u(i, j) = v_tilde_matrix(i, j + 1) & 
                    & - (v_tilde_matrix(i, j + 2) - v_tilde_matrix(i, j)) / 4.0_mp
                    rho_tilde_matrix_u(i, j) = rho_matrix(i, j + 1) & 
                    & - (rho_matrix(i, j + 2) - rho_matrix(i, j)) / 4.0_mp
                end if

            end do
        end do

        do i = 1, N
            do j = 1, N

                if (s_r(i, j) == sgn(u_tilde_matrix_r(i, j))) then
                    Pi_r(i, j) = h * tau * u_tilde_matrix_r(i, j) * rho_tilde_matrix_r(i, j)
                    Pi_u(i, j) = h * tau * v_tilde_matrix_u(i, j) * rho_tilde_matrix_u(i, j)
                end if

            end do
        end do

        do i = 1, N
            do j = 1, N

                rho_matrix_new(i, j) = rho_matrix(i, j) - (Pi_r(i, j) - Pi_r(i - 1, j)) / h ** 2.0_mp - &
                & (Pi_u(i, j) - Pi_u(i, j - 1)) / h ** 2.0_mp

            end do
        end do

    end subroutine mass_transport

    subroutine gas_param_transport(h, tau, rho_matrix, Pi_r, Pi_u, s_r, s_u, rho_matrix_new, f_matrix, f_matrix_new)
        real(mp), intent(in) :: h, tau
        real(mp), intent(in), dimension(0:, 0:) :: rho_matrix, rho_matrix_new, Pi_r, Pi_u, f_matrix
        integer, intent(in), dimension(:, :) :: s_r, s_u
        real(mp), intent(out), dimension(0:size(rho_matrix, dim=1) - 1, 0:size(rho_matrix, dim=1) - 1) :: f_matrix_new
        
        integer :: N, i, j
        real(mp) :: transpote_F

        real(mp), dimension(0:size(rho_matrix, dim=1) - 1, 0:size(rho_matrix, dim=1) - 1) :: f_matrix_r, f_matrix_u

        N = size(rho_matrix, dim=1) - 2

        f_matrix_u = 0.0_mp
        f_matrix_r = 0.0_mp

        do i = 1, N
            do j = 1, N

                if (s_r(i, j) == 1) then
                    f_matrix_r(i, j) = f_matrix(i, j)
                end if

                if (s_r(i, j) == -1) then
                    f_matrix_r(i, j) = f_matrix(i + 1, j)
                end if

                if (s_u(i, j) == 1) then
                    f_matrix_u(i, j) = f_matrix(i, j)
                end if

                if (s_u(i, j) == -1) then
                    f_matrix_u(i, j) = f_matrix(i, j + 1)
                end if

            end do
        end do

        do i = 1, N
            do j = 1, N

                transpote_F = (Pi_r(i, j) * f_matrix_r(i, j) - Pi_r(i - 1, j) * f_matrix_r(i - 1, j)) + &
                            & (Pi_u(i, j) * f_matrix_u(i, j) - Pi_u(i, j - 1) * f_matrix_u(i, j - 1))

                f_matrix_new(i, j) = rho_matrix(i, j) / rho_matrix_new(i, j) * f_matrix(i, j) &
                & - transpote_F / h ** 2.0_mp / rho_matrix_new(i, j)

            end do
        end do

    end subroutine gas_param_transport

end module work_program