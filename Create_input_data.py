import numpy as np


def fun_T_distribution(y, a, b):
    return a * np.power(y, 2.0) + b


if __name__ == '__main__':
    k_B = 1.380649E-16
    a_0 = 0.529E-8

    N = 200
    gamma = 5.0 / 3.0
    rho = 0.17847E-3
    g = 1.0

    R = 83144626.1815324
    M = 4.0E+3

    h = 10.0
    tau = 1.0E-7

    T_f = -25.0 + 273.0
    T_0 = -15.0 + 273.0

    size = h * N

    kappa = 3.0 * k_B / 3.0 / np.power(np.pi, 3.0 / 2.0) / np.power(a_0, 2.0) * np.sqrt(R * (T_0 + T_f) / 2.0 / M)
    c_p = 14.304 * 1.0E+3 * 1.0E+7 / 1.0E+3

    xi = kappa / c_p / rho
    beta = 2.0 / (T_0 + T_f)
    nu = 20.0 * 1.0E-6 * 1.0E+1 / rho

    vazkost = beta / nu / xi

    T = (T_0 + T_f) / 2.0

    print(np.power(1000.0 / vazkost / g / T, 1.0 / 3.0))

    rho_matrix = np.zeros((N, N))
    omega_matrix = np.zeros((N, N))

    v_matrix = np.zeros((N, N))
    u_matrix = np.zeros((N, N))

    epsilon_matrix = np.zeros((N, N))
    pressure_matrix = np.zeros((N, N))

    cond_matrix = np.zeros((N, N))
    cond_matrix_2 = np.zeros((N, N))
    cond_matrix_3 = np.zeros((N, N))

    rho_matrix.fill(rho)

    print(T_f - T_0)

    T_matrix = np.zeros((N, N))

    for j in range(N):
        T_matrix[j, :].fill(fun_T_distribution(h * (j + 1), (T_f - T_0) / np.power(size, 2.0), T_0))

    omega_matrix = 3.0 / 2.0 * R / M * T_matrix

    for j in range(N):
        for i in range(N):
            v_matrix[j, i] = np.random.normal(0.0, 1.0E+2)
            u_matrix[j, i] = np.random.normal(0.0, 1.0E+2)

    epsilon_matrix = omega_matrix + (np.power(v_matrix, 2.0) + np.power(u_matrix, 2.0)) / 2.0
    pressure_matrix = (gamma - 1.0) * rho_matrix * omega_matrix

    for j in range(N - 1):
        for i in range(N - 1):
            cond_matrix[j, i] = abs(T_matrix[j + 1, i] - T_matrix[j, i]) / h >= abs(
                (gamma - 1.0) / gamma * T_matrix[j, i] / pressure_matrix[j, i] * rho_matrix[j, i] * g)
            cond_matrix_2[j, i] = abs(T_matrix[j + 1, i] - T_matrix[j, i]) / h >= g / (c_p)

    cond_matrix_3 = abs(u_matrix) * rho_matrix * tau / h + abs(v_matrix) * rho_matrix * tau / h

    print(cond_matrix)

    print(cond_matrix_2)

    print(cond_matrix_3)

    file_main_parameter = open('file_parameters.dat', 'w')
    file_velocity_x = open('u_matrix.dat', 'w')
    file_velocity_y = open('v_matrix.dat', 'w')
    file_density = open('r_matrix.dat', 'w')
    file_temperature = open('t_matrix.dat', 'w')

    file_main_parameter.write(f'{N}\n')
    file_main_parameter.write(f'{gamma}\n')
    file_main_parameter.write(f'{g}\n')
    file_main_parameter.write(f'{size}\n')
    file_main_parameter.write(f'{h}\n')
    file_main_parameter.write(f'{tau}\n')
    file_main_parameter.write(f'{R}\n')
    file_main_parameter.write(f'{M}\n')

    file_main_parameter.close()

    for j in range(N):
        for i in range(N):
            file_velocity_x.write(f'{u_matrix[j, i]}  ')
        file_velocity_x.write('\n')

    file_velocity_x.close()

    for j in range(N):
        for i in range(N):
            file_velocity_y.write(f'{v_matrix[j, i]}  ')
        file_velocity_y.write('\n')

    file_velocity_y.close()

    for j in range(N):
        for i in range(N):
            file_density.write(f'{rho_matrix[j, i]}  ')
        file_density.write('\n')

    file_density.close()

    for j in range(N):
        for i in range(N):
            file_temperature.write(f'{T_matrix[j, i]}  ')
        file_temperature.write('\n')

    file_temperature.close()
