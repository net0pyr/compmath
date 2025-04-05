import numpy as np
import sympy as sp

# Определение переменных и функций
x, y = sp.symbols('x y')
f1 = sp.tan(x * y) - x**2
f2 = 0.8 * x**2 + 2 * y**2 - 1

# Вычисление якобиана
J = sp.Matrix([[f1.diff(x), f1.diff(y)], [f2.diff(x), f2.diff(y)]])

# Преобразование символьных функций в численные
f1_num = sp.lambdify((x, y), f1, 'numpy')
f2_num = sp.lambdify((x, y), f2, 'numpy')
J_num = sp.lambdify((x, y), J, 'numpy')

# Метод Ньютона
def newton_method(f1, f2, J, x0, y0, tol=0.01, max_iter=100):
    x, y = x0, y0
    iterations = []
    for i in range(max_iter):
        F = np.array([f1(x, y), f2(x, y)])
        J_inv = np.linalg.inv(J(x, y))
        delta = -J_inv @ F
        x_new = x + delta[0]
        y_new = y + delta[1]
        iterations.append((x, y, F, J(x, y), delta, x_new, y_new))
        if np.linalg.norm(delta) < tol:
            break
        x, y = x_new, y_new
    return x, y, iterations

# Начальные приближения
x0, y0 = 0.5, 0.5

# Решение системы
x_sol, y_sol, iterations = newton_method(f1_num, f2_num, J_num, x0, y0)

# Генерация отчета
latex_output = r"""
\documentclass{article}
\usepackage{amsmath}
\begin{document}

\section*{Решение системы нелинейных уравнений методом Ньютона}

Дана система уравнений:
\[
\begin{cases}
\tan(xy) = x^2 \\
0.8x^2 + 2y^2 = 1
\end{cases}
\]

\subsection*{Метод Ньютона}

Метод Ньютона заключается в итерационном уточнении решения с использованием якобиана системы. Якобиан системы:
\[
J(x, y) = \begin{pmatrix}
\frac{\partial f_1}{\partial x} & \frac{\partial f_1}{\partial y} \\
\frac{\partial f_2}{\partial x} & \frac{\partial f_2}{\partial y}
\end{pmatrix}
\]
где
\[
\frac{\partial f_1}{\partial x} = y \cdot \sec^2(xy) - 2x, \quad
\frac{\partial f_1}{\partial y} = x \cdot \sec^2(xy)
\]
\[
\frac{\partial f_2}{\partial x} = 1.6x, \quad
\frac{\partial f_2}{\partial y} = 4y
\]

\subsection*{Итерации}
"""

for i, (x_n, y_n, F, J_n, delta, x_new, y_new) in enumerate(iterations):
    latex_output += f"""
\subsubsection*{{Итерация {i + 1}}}

Текущее приближение:
\[
x_n = {x_n:.6f}, \quad y_n = {y_n:.6f}
\]

Значения функций:
\[
f_1(x_n, y_n) = \\tan({x_n:.6f} \\cdot {y_n:.6f}) - ({x_n:.6f})^2 = {F[0]:.6f}
\]
\[
f_2(x_n, y_n) = 0.8 \\cdot ({x_n:.6f})^2 + 2 \\cdot ({y_n:.6f})^2 - 1 = {F[1]:.6f}
\]

Якобиан на текущей итерации:
\[
J(x_n, y_n) = \\begin{{pmatrix}}
{J_n[0, 0]:.6f} & {J_n[0, 1]:.6f} \\\\
{J_n[1, 0]:.6f} & {J_n[1, 1]:.6f}
\\end{{pmatrix}}
\]

Решаем систему линейных уравнений для нахождения приращений:
\[
\\begin{{pmatrix}}
{J_n[0, 0]:.6f} & {J_n[0, 1]:.6f} \\\\
{J_n[1, 0]:.6f} & {J_n[1, 1]:.6f}
\\end{{pmatrix}}
\\begin{{pmatrix}}
\\Delta x_n \\\\
\\Delta y_n
\\end{{pmatrix}}
= -\\begin{{pmatrix}}
{F[0]:.6f} \\\\
{F[1]:.6f}
\\end{{pmatrix}}
\]

Приращения:
\[
\\Delta x_n = {delta[0]:.6f}, \quad \\Delta y_n = {delta[1]:.6f}
\]

Новое приближение:
\[
x_{{n+1}} = x_n + \\Delta x_n = {x_n:.6f} + {delta[0]:.6f} = {x_new:.6f}
\]
\[
y_{{n+1}} = y_n + \\Delta y_n = {y_n:.6f} + {delta[1]:.6f} = {y_new:.6f}
\]
"""

latex_output += r"""
\subsection*{Результат}

После нескольких итераций метод Ньютона сходится к решению:
\[
x \approx """ + f"{x_sol:.6f}" + r", \quad y \approx """ + f"{y_sol:.6f}" + r"""
\]

\end{document}
"""

# Сохранение LaTeX-кода в файл
with open("newton_solution_report_detailed.tex", "w") as f:
    f.write(latex_output)

print("LaTeX-код сгенерирован и сохранен в файл newton_solution_report_detailed.tex")
