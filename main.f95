program main
    use :: work_function
    use :: work_program
    implicit none

    integer :: N, i
    real(mp), dimension(:, :), allocatable :: u_matrix, v_matrix, rho_matrix, T_matrix
    real(mp) :: gamma, g, a, h, tau, R, M

    open(11, file='file_parameters.dat')

    read(11, *) N
    read(11, *) gamma
    read(11, *) g
    read(11, *) a
    read(11, *) h
    read(11, *) tau
    read(11, *) R
    read(11, *) M

    close(11)

    allocate(T_matrix(0:N + 1, 0:N + 1), u_matrix(0:N+1, 0:N+1), v_matrix(0:N+1, 0:N+1), rho_matrix(0:N+1, 0:N+1))

    T_matrix = 0.0_mp
    u_matrix = 0.0_mp
    v_matrix = 0.0_mp
    rho_matrix = 0.0_mp

    T_matrix(1:N, 1:N) = read_matrix(N, 't_matrix.dat')
    u_matrix(1:N, 1:N) = read_matrix(N, 'u_matrix.dat')
    v_matrix(1:N, 1:N) = read_matrix(N, 'v_matrix.dat')
    rho_matrix(1:N, 1:N) = read_matrix(N, 'r_matrix.dat')

    call method_Davydov_Belotserkovsky(gamma, g, R, M, h, tau, u_matrix, v_matrix, T_matrix, rho_matrix)

    deallocate(T_matrix, u_matrix, v_matrix, rho_matrix)

end program main