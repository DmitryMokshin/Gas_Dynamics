import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


def update_graphics(i, img, data_file):
    img.set_data(data_file[i])
    return img


def animation_parameter(file_input, parameter):
    file_for_analys = open(file=file_input)

    data_analys = file_for_analys.readlines()

    matrix_param = []

    matrix_time = []

    for line in data_analys:
        if len(line) != 1:
            line_matrix = []
            for symbol in line.split():
                line_matrix.append(float(symbol))
            matrix_param.append(line_matrix)
        if len(line) == 1:
            matrix_time.append(np.array(matrix_param))
            matrix_param = []

    print('animation build')

    print('Num of steps: ', len(matrix_time))

    fig, ax = plt.subplots()
    img = ax.imshow(matrix_time[0])

    ani = animation.FuncAnimation(fig, update_graphics, fargs=(img, matrix_time), frames=len(matrix_time), interval=1,
                                  repeat=False, save_count=len(matrix_time))
    fig.colorbar(img)

    writergif = animation.PillowWriter(fps=30)
    ani.save(parameter + '.gif', writer=writergif)

    plt.close()


if __name__ == '__main__':
    animation_parameter('rho_result_matrix.dat', 'density')
