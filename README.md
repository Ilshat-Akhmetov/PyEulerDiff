# PyEulerDiff

In this repository you may find implementation of various numerical 
method for solving Cauchy problem for ODEs and systems of ODEs. 
Currently, repository includes implemented implicit 
and explicit classic Euler methods, classic explicit 
Runge-Kutta method AKA 'RK4'. In the future I have plans to extend 
current list of implemented numerical methods with 
others like modified Euler, Adams–Bashforth methods and 
add a possibility to select dynamic step size to make sure 
error does not exceed a threshold given by the user and so on.
I have plans to implement a nice GUI for this repo to make 
this program available to people not familiar with python.

## References

More about Cauchy problem

https://en.wikipedia.org/wiki/Cauchy_problem

Adam-Bashforth method

https://en.wikipedia.org/wiki/Linear_multistep_method

Runge-Kutta method

https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

Classic Euler method

https://en.wikipedia.org/wiki/Euler_method

## How to use?

I assume you are in PyEulerDiff directory. First of all, import everything from SourceCode 
with `from SourceCode import *`. Below you may find 2 examples how to use this repo. 
API is simple and intuitive. Program assumes that you would 
provide it with an equation in the form 

$$ y'=f(x,y) $$

You only have to provide the right part f(x,y).

Let's start with the first example.
Say, you want to solve a simple problem:

$$ y'- y = 0 $$

$$ y(0) = 1 $$

Exact solution:

$$ y = exp(x) $$

Then you would have to write the following code.

```python
from SourceCode import *

right_part = ["y"]
a = 1.00
init_vals = [a]
true_sol = lambda t: a*np.exp(t)
left = 0
right = 1
n_points = 101
model = SoODESolver(init_vals, right_part, left, right, n_points)
ans = model.explicit_euler()

true_vals = true_sol(model.domain)
true_vals.resize(n_points, 1)
error = norm_Cab_error(true_vals, ans)
print(error)
```

comparing with the exact solution is optional, because normally exact solution 
is not know and cannot be obtained.

second example

$$ y1'=-y2 $$

$$ y2'=y1+cos(t) $$

$$ y1(0)=0, y2(0)=0 $$

That is how you can solve this equation with this repo:
```python
from SourceCode import *

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
```

Keep in mind that in case if you need to solve a system of equation 
the program expects you to denote unknown functiona as y1, y2 .. yn where n - number of equations 
but if you need to solve just one equation then you should denote an unknown function simply as y. 
In both cases you should denote variable as t. You should not to 
explicitly specify that y is a function of t as y(t) in the right part, program assumes it 
by default. 

Current implementation supports only constant uniform step, eventually I will 
add some methods supporting an adaptive step.





