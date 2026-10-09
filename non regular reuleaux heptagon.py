import numpy as np
from scipy import integrate

# =================================================================
#                             Newton loop
# =================================================================

def integrate_piecewise(integrand, breakpoints, limit=200):
    """integrate `integrand` over [breakpoints[0], breakpoints[-1]] taking into account the seven intervals"""
    total = 0.0
    for i in range(len(breakpoints) - 1):
        lo, hi = breakpoints[i], breakpoints[i + 1]
        val, _ = integrate.quad(integrand, lo, hi, limit=limit)
        total += val
    return total


def F_and_J(a, b, h_func, breakpoints):
    """Computation of F(a,b)=(f1,f2) and of the Jacobian J(a,b)"""

    def D(t):
        return h_func(t) - a * np.cos(t) - b * np.sin(t)

    f1 = integrate_piecewise(lambda t: np.cos(t) / D(t) ** 3, breakpoints)
    f2 = integrate_piecewise(lambda t: np.sin(t) / D(t) ** 3, breakpoints)

    j11 = integrate_piecewise(lambda t: 3 * np.cos(t) ** 2 / D(t) ** 4, breakpoints)
    j12 = integrate_piecewise(lambda t: 3 * np.cos(t) * np.sin(t) / D(t) ** 4, breakpoints)
    j22 = integrate_piecewise(lambda t: 3 * np.sin(t) ** 2 / D(t) ** 4, breakpoints)

    F = np.array([f1, f2])
    J = np.array([[j11, j12],
                  [j12, j22]])
    return F, J


def santalo_newton(h_func, breakpoints, s0=(0.0, 0.0),
                    tol=1e-13, maxiter=50, verbose=True):
    """Resolution of F(a,b)=0 via Newton"""
    a, b = s0
    history = []

    for k in range(maxiter):
        F, J = F_and_J(a, b, h_func, breakpoints)
        history.append((k, a, b, np.linalg.norm(F)))
        if verbose:
            print(f"iter {k:2d} : s=({a: .14f},{b: .14f})  ||F||={np.linalg.norm(F):.3e}")
        if np.linalg.norm(F) < tol:
            break
        delta = np.linalg.solve(J, F)
        a, b = a - delta[0], b - delta[1]

    def D(t):
        return h_func(t) - a * np.cos(t) - b * np.sin(t)
    polar_area = 0.5 * integrate_piecewise(lambda t: 1.0 / D(t) ** 2, breakpoints)

    return (a, b), polar_area, history


def body_area(h_func, dh_func, breakpoints):
    """Area of K :  A = 1/2 ∫ (h(t)^2 - h'(t)^2) dt"""
    def integrand(t):
        return h_func(t) ** 2 - dh_func(t) ** 2
    return 0.5 * integrate_piecewise(integrand, breakpoints)


# =================================================================
# Definition of the support function h(t) (seven pieces over [0,pi])
#        (maple export of switch times and support function)
# =================================================================

t1 = 0.2499109628
t2 = 0.7339481128
t3 = 1.205338126
t4 = 1.408149043
t5 = 1.875328710
t6 = 2.642468712
pi = np.pi

breakpoints_half = [0, t1, t2, t3, t4, t5, t6, pi]

# coefficients (a_i, b_i, c_i)  h_i(t) = a_i*cos(t) + b_i*sin(t) + c_i
a1, b1, c1 = 0.5,            0.0,            0.0
a2, b2, c2 = -0.468934446,  -0.247317689,    1.0
a3, b3, c3 = 0.2736012789,   0.4224887735,   0.0
a4, b4, c4 = -0.0837759943, -0.5114713352,   1.0
a5, b5, c5 = 0.0781551217,   0.4753307289,   0.0
a6, b6, c6 = 0.3780022292,  -0.4786565415,   1.0
a7, b7, c7 = -0.500000001,  -4e-10,          0.0

pieces_half = [(a1, b1, c1), (a2, b2, c2), (a3, b3, c3), (a4, b4, c4),
               (a5, b5, c5), (a6, b6, c6), (a7, b7, c7)]

# extension to [pi, 2pi] via h(t+pi) = 1 - h(t)  ->  same (a,b), c -> 1-c
pieces_full = pieces_half + [(a, b, 1.0 - c) for (a, b, c) in pieces_half]
breakpoints_full = breakpoints_half + [pi + x for x in breakpoints_half[1:]]


def which_piece(t):
    tt = t % (2 * pi)
    for i in range(len(breakpoints_full) - 1):
        if breakpoints_full[i] - 1e-9 <= tt <= breakpoints_full[i + 1] + 1e-9:
            return i
    raise ValueError(f"t={t} hors domaine")


def h(t):
    i = which_piece(t)
    a, b, c = pieces_full[i]
    return a * np.cos(t) + b * np.sin(t) + c


def dh(t):
    i = which_piece(t)
    a, b, c = pieces_full[i]
    return -a * np.sin(t) + b * np.cos(t)


# =================================================================
#             Santalo point computation
# =================================================================

A = body_area(h, dh, breakpoints_full)
print("Area of K :", A)

s0 = (0.0, 0.0)
(a_star, b_star), polar_area, hist = santalo_newton(
    h, breakpoints_full, s0=s0, tol=1e-13, verbose=True
)

print("Santalo point  s* =", (a_star, b_star))
print("Area of K° at s* :", polar_area)
print("Area product A(K) * A(K°) :", A * polar_area)
