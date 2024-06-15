from typing import List, Union
from sympy import lambdify
import numpy as np
from scipy.optimize import fsolve


class SoODESolver:
    def __init__(self, init_vals: List[Union[float, int]],
                 right_part: List[str],
                 left_b: Union[float, int],
                 right_b: Union[float, int],
                 n_points: int):
        self.dim = len(init_vals)

        assert self.dim == len(right_part), ("number of init vals should be equal to the number of"
                                             "equations in the right part")
        self.vars = ['t']
        if self.dim == 0:
            raise ValueError('system of equations should consist at least of one equation!')
        elif self.dim == 1:
            self.vars.append('y')
        else:
            self.vars.extend(['y' + str(i) for i in range(1, self.dim + 1)])
        self.right_part_eqs = [lambdify(self.vars, f_right) for f_right in right_part]
        self.domain = np.linspace(left_b, right_b, n_points)
        self.n_points = n_points
        self.init_vals = np.array(init_vals)
        self.init_vals.resize(1, self.dim)
        self.dt = (right_b - left_b) / (n_points - 1)

    def init_uniform_array(self) -> np.ndarray:
        answer = np.zeros((self.n_points, self.dim), np.float32)
        answer[0] = self.init_vals
        return answer

    def explicit_euler(self) -> np.ndarray:
        answer = self.init_uniform_array()
        for iter_i in range(self.n_points - 1):
            f_term = self.dt * np.array([right_p(self.domain[iter_i], *answer[iter_i])
                                         for right_p in self.right_part_eqs])
            answer[iter_i + 1] = answer[iter_i] + f_term
        return answer

    def implicit_euler(self) -> np.ndarray:
        answer = self.init_uniform_array()

        def equation(y_vals, *data):
            t, y_prev = data
            n = len(self.right_part_eqs)
            return [y_vals[i] - y_prev[i] - self.dt*self.right_part_eqs[i](t, *y_vals)
                    for i in range(n)]

        for iter_i in range(self.n_points - 1):
            answer[iter_i + 1] = fsolve(equation, answer[iter_i], args=(self.domain[iter_i + 1], answer[iter_i]))
        return answer

    def explicit_rungekutta4(self) -> np.ndarray:
        answer = self.init_uniform_array()
        for iter_i in range(self.n_points - 1):
            k1 = np.array([right_p(self.domain[iter_i], *answer[iter_i])
                                         for right_p in self.right_part_eqs])
            k2 = np.array([right_p(self.domain[iter_i] + self.dt/2, *answer[iter_i] + self.dt*k1/2)
                                         for right_p in self.right_part_eqs])
            k3 = np.array([right_p(self.domain[iter_i] + self.dt/2, *answer[iter_i] + self.dt*k2/2)
                                         for right_p in self.right_part_eqs])
            k4 = np.array([right_p(self.domain[iter_i] + self.dt, *answer[iter_i] + self.dt * k3)
                           for right_p in self.right_part_eqs])
            f_term = self.dt/6 * (k1+2*k2+2*k3+k4)
            answer[iter_i + 1] = answer[iter_i] + f_term
        return answer

