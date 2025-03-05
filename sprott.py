import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d
from matplotlib.widgets import Button, TextBox
import subprocess

def load_and_plot_data(event=None):
    # Загрузка данных RK4 (4 колонки: t, x, y, z)
    rk4_data = np.loadtxt('sprott_rk4.txt', usecols=(0, 1, 2, 3))
    t_rk4 = rk4_data[:, 0]
    x_rk4 = rk4_data[:, 1]
    y_rk4 = rk4_data[:, 2]
    z_rk4 = rk4_data[:, 3]

    # Загрузка данных dopri5 с адаптивным шагом (6 колонок: t, x, y, z, локальная ошибка, dt)
    dopri5_data = np.loadtxt('sprott_dopri5_adaptive.txt', usecols=(0, 1, 2, 3, 4, 5))
    t_dopri5_a = dopri5_data[:, 0]
    x_dopri5_a = dopri5_data[:, 1]
    y_dopri5_a = dopri5_data[:, 2]
    z_dopri5_a = dopri5_data[:, 3]
    local_err = dopri5_data[:, 4]  
    dt_dopri5 = dopri5_data[:, 5]    

    # Очистка предыдущих графиков
    ax.cla()
    ax2.cla()
    ax3.cla()
    ax4.cla()

    # Построение 3D-графика
    ax.plot(x_rk4, y_rk4, z_rk4, label='RK4', lw=0.5, color='blue')
    ax.plot(x_dopri5_a, y_dopri5_a, z_dopri5_a, label='dopri5 adaptive', lw=0.5, color='red')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.view_init(elev=30, azim=45)
    ax.set_title('Аттрактор Sprott: RK4 vs dopri5')
    ax.legend()

    # Интерполируем данные dopri5 на временную сетку RK4 для вычисления ошибки
    interp_x = interp1d(t_dopri5_a, x_dopri5_a, kind='linear', fill_value="extrapolate")
    interp_y = interp1d(t_dopri5_a, y_dopri5_a, kind='linear', fill_value="extrapolate")
    interp_z = interp1d(t_dopri5_a, z_dopri5_a, kind='linear', fill_value="extrapolate")

    x_dopri5_interp = interp_x(t_rk4)
    y_dopri5_interp = interp_y(t_rk4)
    z_dopri5_interp = interp_z(t_rk4)

    # Вычисление глобальной ошибки (евклидова норма разности состояний)
    error_global = np.sqrt((x_rk4 - x_dopri5_interp)**2 +
                           (y_rk4 - y_dopri5_interp)**2 +
                           (z_rk4 - z_dopri5_interp)**2)

    # Построение графика глобальной ошибки
    ax2.plot(t_rk4, error_global, label='Global Error', color='magenta')
    ax2.set_xlabel('Время')
    ax2.set_ylabel('Ошибка (евклидова норма)')
    ax2.set_title('Глобальная разность между RK4 и dopri5')
    ax2.legend()
    
    # Построение графика локальной ошибки dopri5
    ax3.plot(t_dopri5_a, local_err, label='Local Error', color='green')
    ax3.set_xlabel('Время')
    ax3.set_ylabel('Локальная ошибка')
    ax3.set_title('Локальная ошибка dopri5')
    ax3.legend()
    
    # Построение графика шага интегрирования dopri5
    ax4.plot(t_dopri5_a, dt_dopri5, label='dt', color='orange')
    ax4.set_xlabel('Время')
    ax4.set_ylabel('Шаг интегрирования')
    ax4.set_title('Шаг интегрирования dopri5')
    ax4.legend()

    plt.draw()

def aaa(event=None):
    try:
        arg = int(text_box.text)
    except ValueError:
        return
    
    if 0 < arg < 5000:
        run_cpp_sprott(event=event)
    
    return

def run_cpp_sprott(event=None, exe='./a.out'):
    arg = text_box.text
    subprocess.run([exe, arg])
    load_and_plot_data()

# Создание фигуры и осей
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(1, 4, (1,2), projection='3d')
ax2 = fig.add_subplot(2, 4, (3, 4))
ax3 = fig.add_subplot(247)
ax4 = fig.add_subplot(248)

# Текстовое поле для ввода аргумента
ax_text_box = plt.axes([0.25, 0.01, 0.15, 0.05])
text_box = TextBox(ax_text_box, 'Время')
text_box.on_text_change(aaa)
text_box.on_submit(run_cpp_sprott)

try:
    load_and_plot_data()
except FileNotFoundError:
    subprocess.run(['./a.out', '1000'])
    load_and_plot_data()

plt.tight_layout()
plt.show()
