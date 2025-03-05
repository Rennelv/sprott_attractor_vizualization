#include <algorithm>  // для std::max
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

// Структура для хранения состояния системы
struct State {
    double x, y, z;

    State operator*(double k) const {
        return {x * k, y * k, z * k};
    }

    State operator+(const State &s) const {
        return {x + s.x, y + s.y, z + s.z};
    }

    State operator-(const State &s) const {
        return {x - s.x, y - s.y, z - s.z};
    }
};

#define A 2.07
#define B 1.79

State f(const State &s) {
    // x = y + axy + xz
    // y = 1 - bx^2 + yz
    // z = x - x^2 - y^2

    return {
        s.y + A * s.x * s.y + s.x * s.z,
        1.0 - B * s.x * s.x + s.y * s.z,
        s.x - s.x * s.x - s.y * s.y,
    };
}

void rk4(double T, double dt, State s, double t_start = 0.0, std::string save_file = "sprott_rk4.txt") {
    double t = 0.0;  // начальное время

    std::ofstream fout(save_file);
    if (!fout) {
        std::cerr << "Ошибка открытия файла для записи\n";
        throw std::runtime_error("Ошибка открытия файла для записи");
    }

    // Используем строковый поток для кэширования данных
    std::ostringstream buffer;

    // Записываем начальное состояние
    if (t >= t_start) buffer << t << " " << s.x << " " << s.y << " " << s.z << "\n";

    // Интегрирование методом RK4
    while (t + dt < T) {
        State k1 = f(s) * dt;
        State k2 = f(s + k1 * 0.5) * dt;
        State k3 = f(s + k2 * 0.5) * dt;
        State k4 = f(s + k3) * dt;

        s.x += (k1.x + 2 * k2.x + 2 * k3.x + k4.x) / 6.0;
        s.y += (k1.y + 2 * k2.y + 2 * k3.y + k4.y) / 6.0;
        s.z += (k1.z + 2 * k2.z + 2 * k3.z + k4.z) / 6.0;

        t += dt;
        if (t >= t_start) buffer << t << " " << s.x << " " << s.y << " " << s.z << "\n";
    }

    fout << buffer.str();
    fout.close();
}

void dopri5_adaptive(double T, double dt, State s, double t_start = 0.0, std::string save_file = "sprott_dopri5_adaptive.txt") {
    double t = 0.0;       // начальное время
    double tol = 1e-6;    // допустимый локальный шаговой допуск
    double safety = 0.9;  // коэффициент безопасности для изменения шага
    double p = 5.0;       // порядок метода для оценки (экспонента = 1/5)

    std::ofstream fout(save_file);
    if (!fout) {
        std::cerr << "Ошибка открытия файла для записи\n";
        throw std::runtime_error("Ошибка открытия файла для записи");
    }

    // Используем строковый поток для кэширования данных
    std::ostringstream buffer;

    // Записываем начальное состояние
    if (t >= t_start) buffer << t << " " << s.x << " " << s.y << " " << s.z << " " << 0.0 << " " << dt << "\n";

    // Основной цикл интегрирования
    while (t < T) {
        // Если следующий шаг выходит за T корректируем dt
        if (t + dt > T) {
            dt = T - t;
        }

        bool stepAccepted = false;
        double dt_new = dt;  // новый шаг, вычисляемый на основе ошибки

        // Повторяем попытки шага, пока ошибка не удовлетворит допуск
        while (!stepAccepted) {
            // Вычисляем 7 стадий метода Dormand–Prince:
            State k1 = f(s) * dt;
            State k2 = f(s + k1 * (1.0 / 5.0)) * dt;
            State k3 = f(s + k1 * (3.0 / 40.0) + k2 * (9.0 / 40.0)) * dt;
            State k4 = f(s + k1 * (44.0 / 45.0) - k2 * (56.0 / 15.0) + k3 * (32.0 / 9.0)) * dt;
            State k5 = f(s + k1 * (19372.0 / 6561.0) - k2 * (25360.0 / 2187.0) + k3 * (64448.0 / 6561.0) - k4 * (212.0 / 729.0)) * dt;
            State k6 = f(s + k1 * (9017.0 / 3168.0) - k2 * (355.0 / 33.0) + k3 * (46732.0 / 5247.0) + k4 * (49.0 / 176.0) - k5 * (5103.0 / 18656.0)) * dt;
            // Седьмая стадия используется для оценки 4‑порядкового решения
            State k7 = f(s + k1 * (35.0 / 384.0) + k3 * (500.0 / 1113.0) + k4 * (125.0 / 192.0) - k5 * (2187.0 / 6784.0) + k6 * (11.0 / 84.0)) * dt;

            // Вычисляем 5‑порядковое решение (принимаемое)
            State s5 = s + k1 * (35.0 / 384.0) + k3 * (500.0 / 1113.0) + k4 * (125.0 / 192.0) - k5 * (2187.0 / 6784.0) + k6 * (11.0 / 84.0);

            // Вычисляем 4‑порядковое решение для оценки ошибки
            State s4 = s + k1 * (5179.0 / 57600.0) + k3 * (7571.0 / 16695.0) + k4 * (393.0 / 640.0) - k5 * (92097.0 / 339200.0) + k6 * (187.0 / 2100.0) +
                       k7 * (1.0 / 40.0);

            // Оценка локальной ошибки по макс ошибке координат
            State error = s5 - s4;
            double err = std::max({std::fabs(error.x), std::fabs(error.y), std::fabs(error.z)});

            // Расчёт нового шага dt_new на основе ошибки
            dt_new = dt * safety * std::pow(tol / (err + 1e-10), 1.0 / p);

            // Если ошибка удовлетворяет допуск, шаг принимается
            if (err <= tol) {
                t += dt;
                s = s5;
                stepAccepted = true;
                if (t >= t_start) buffer << t << " " << s.x << " " << s.y << " " << s.z << " " << err << " " << dt << "\n";
            } else {
                // Если ошибка слишком велика, шаг отклоняется и уменьшается
                dt = dt_new;
            }
        }
        // Обновляем шаг для следующей итерации
        dt = dt_new;
        // Дополнительные ограничения на диапазон dt (при необходимости)
        if (dt < 1e-10) dt = 1e-10;
        if (dt > 0.1) dt = 0.1;
    }

    fout << buffer.str();
    fout.close();
}

int main(int argc, char *argv[]) {
    double T = 1000.0;  // конечное время
    double dt = 0.01;   // начальный шаг интегрирования
    double t_start = 0.0;
    std::string save_file_rk4 = "sprott_rk4.txt";
    std::string save_file_dopri5 = "sprott_dopri5_adaptive.txt";

    // Начальные условия
    State s{0.63, 0.47, -0.54};

    // T
    if (argc > 1) {
        T = std::atof(argv[1]);
        if (T <= 0) {
            return 1;
        }
    }
    // dt
    if (argc > 2) {
        dt = std::atof(argv[2]);
        if (dt <= 0) {
            return 1;
        }
    }
    // t start
    if (argc > 3) {
        t_start = std::atof(argv[3]);
        if (t_start < 0) {
            return 1;
        }
    }

    std::cout << "Время моделирования: " << T << "\n";

    auto start = std::clock();
    rk4(T, dt, s, t_start, save_file_rk4);
    auto end = std::clock();
    std::cout << "Время выполнения метода RK4: " << double(end - start) / CLOCKS_PER_SEC << " секунд\n";

    start = std::clock();
    dopri5_adaptive(T, dt, s, t_start, save_file_dopri5);
    end = std::clock();
    std::cout << "Время выполнения метода Dopri5: " << double(end - start) / CLOCKS_PER_SEC << " секунд\n";
    return 0;
}