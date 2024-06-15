from SourceCode import *

def test1_expl_euler():
    right_part = ["y"]
    a = 1.00
    init_vals = [a]
    true_sol = lambda t: a*np.exp(t)
    left = 0
    right = 1
    n_points = 101
    model = SoODESolver(init_vals, right_part, left, right, n_points)
    ans = model.explicit_euler()
    print(ans[-1])

    true_vals = true_sol(model.domain)
    true_vals.resize(n_points, 1)
    error = norm_Cab_error(true_vals, ans)
    print(error)


def test2_impl_euler():
    right_part = ["y"]
    init_vals = [1]
    true_sol = lambda t: np.exp(t)
    left = 0
    right = 1
    n_points = 101
    model = SoODESolver(init_vals, right_part, left, right, n_points)
    ans = model.implicit_euler()

    true_vals = true_sol(model.domain)
    true_vals.resize(n_points, 1)
    error = norm_Cab_error(true_vals, ans)
    print(error)


def test3_system_expl_euler():
    right_part = ['-y2', 'y1+cos(t)']
    init_vals = [0, 0]
    true_sol = [lambda t: -0.5 * t * np.sin(t), lambda t: 0.5 * (t * np.cos(t) + np.sin(t))]
    left = 0
    right = 5
    n_points = 10001
    model = SoODESolver(init_vals, right_part, left, right, n_points)
    ans = model.explicit_euler()
    true_vals = np.stack([f(model.domain) for f in true_sol], axis=1)
    error = norm_Cab_error(true_vals, ans)
    print(error)


def test4_system_expl_runge_kutta4():
    right_part = ['-y2', 'y1+cos(t)']
    init_vals = [0, 0]
    true_sol = [lambda t: -0.5 * t * np.sin(t), lambda t: 0.5 * (t * np.cos(t) + np.sin(t))]
    left = 0
    right = 5
    n_points = 101
    model = SoODESolver(init_vals, right_part, left, right, n_points)
    ans = model.explicit_rungekutta4()
    true_vals = np.stack([f(model.domain) for f in true_sol], axis=1)
    error = norm_Cab_error(true_vals, ans)
    print(error)


def test5_system_expl_runge_kutta4():
    right_part = ['-y2', 'y1+cos(t)']
    init_vals = [0, 0]
    true_sol = [lambda t: -0.5 * t * np.sin(t), lambda t: 0.5 * (t * np.cos(t) + np.sin(t))]
    left = 0
    right = 1
    n_points = 101
    model = SoODESolver(init_vals, right_part, left, right, n_points)
    ans = model.explicit_rungekutta4()
    true_vals = np.stack([f(model.domain) for f in true_sol], axis=1)
    error = norm_Cab_error(true_vals, ans)
    print(error)


def test6_expl_runge_kutta4():
    right_part = ['exp(-t/5)*cos(t)-0.2*y']
    init_vals = [0]
    true_sol = lambda t: np.exp(-t / 5) * np.sin(t)
    left = 0
    right = 1
    n_points = 10
    model = SoODESolver(init_vals, right_part, left, right, n_points)
    ans = model.explicit_rungekutta4()
    true_vals = true_sol(model.domain)
    true_vals.resize(n_points, 1)
    error = norm_Cab_error(true_vals, ans)
    print(error)


def test7_expl_euler():
    right_part = ['exp(-t/5)*cos(t)-0.2*y']
    init_vals = [0]
    true_sol = lambda t: np.exp(-t / 5) * np.sin(t)
    left = 0
    right = 5
    n_points = 100
    model = SoODESolver(init_vals, right_part, left, right, n_points)
    ans = model.explicit_euler()
    true_vals = true_sol(model.domain)
    true_vals.resize(n_points, 1)
    error = norm_Cab_error(true_vals, ans)
    print(error)


# test1_expl_euler()
# test2_impl_euler()
# test3_system_expl_euler()
test6_expl_runge_kutta4()
