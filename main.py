from SourceCode.SoODESolver import SoODESolver
import numpy as np

right_part = ["y"]
init_vals = [1]
true_sol = lambda t: np.exp(t)
left = 0
right = 1
n_points = 1001
model = SoODESolver(init_vals, right_part, left, right, n_points)
ans = model.explicit_euler()
#ans = model.implicit_euler()

norm_Cab_error = lambda exact_sol, approx: np.max(np.abs(exact_sol - approx))
true_vals = true_sol(model.domain)
true_vals.resize(n_points, 1)
error = norm_Cab_error(true_vals, ans)
print(error)

