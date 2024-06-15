import numpy as np

get_error = lambda exact_sol, approx: np.sqrt(np.power((exact_sol - approx), 2).sum(axis=1))
norm_Cab_error = lambda exact_sol, approx: np.max(get_error(exact_sol, approx))