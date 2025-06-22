import numpy as np
import plotext as plt
import caltechschool as cs
import argparse


def euler(
    x: np.ndarray, u: np.ndarray, t: float, dt: float, E: cs.Field, B: cs.Field
) -> tuple[np.ndarray, np.ndarray]:
    xlast = x[-1]
    ulast = u[-1]
    Eval = E(xlast, t)
    Bval = B(xlast, t)
    gammalast = np.sqrt(1.0 + np.linalg.norm(ulast) ** 2)
    xnew = xlast + dt * ulast / gammalast
    unew = ulast + dt * (Eval + np.cross(ulast, Bval) / gammalast)
    return xnew, unew


def rk4(
    x: np.ndarray, u: np.ndarray, t: float, dt: float, E: cs.Field, B: cs.Field
) -> tuple[np.ndarray, np.ndarray]:
    xlast = x[-1]
    ulast = u[-1]
    Eval = E(xlast, t)
    Bval = B(xlast, t)
    gammalast = np.sqrt(1.0 + np.linalg.norm(ulast) ** 2)

    def rhs_x(u: np.ndarray) -> np.ndarray:
        return u / np.sqrt(1.0 + np.linalg.norm(u) ** 2)

    def rhs_u(u: np.ndarray) -> np.ndarray:
        return Eval + np.cross(u, Bval) / np.sqrt(1.0 + np.linalg.norm(u) ** 2)

    k1x = rhs_x(ulast)
    k1u = rhs_u(ulast)

    k2x = rhs_x(ulast + 0.5 * dt * k1u)
    k2u = rhs_u(ulast + 0.5 * dt * k1u)

    k3x = rhs_x(ulast + 0.5 * dt * k2u)
    k3u = rhs_u(ulast + 0.5 * dt * k2u)

    k4x = rhs_x(ulast + dt * k3u)
    k4u = rhs_u(ulast + dt * k3u)

    xnew = xlast + dt * (k1x + 2 * k2x + 2 * k3x + k4x) / 6
    unew = ulast + dt * (k1u + 2 * k2u + 2 * k3u + k4u) / 6

    return xnew, unew


def boris(
    x: np.ndarray, u: np.ndarray, t: float, dt: float, E: cs.Field, B: cs.Field
) -> tuple[np.ndarray, np.ndarray]:
    xlast = x[-1]
    ulast = u[-1]
    Eval = E(xlast, t)
    Bval = B(xlast, t)

    um = ulast + 0.5 * dt * Eval
    gammam = np.sqrt(1.0 + np.linalg.norm(um) ** 2)

    kappa = 0.5 * dt / gammam

    uprime = um + kappa * np.cross(um, Bval)
    uplus = um + 2 * kappa / (1 + kappa**2 * np.linalg.norm(Bval) ** 2) * np.cross(
        uprime, Bval
    )
    unew = uplus + 0.5 * dt * Eval
    xnew = xlast + dt * unew / np.sqrt(1.0 + np.linalg.norm(unew) ** 2)
    return xnew, unew


def implicit(
    x: np.ndarray, u: np.ndarray, t: float, dt: float, E: cs.Field, B: cs.Field
) -> tuple[np.ndarray, np.ndarray]:
    xlast = x[-1]
    ulast = u[-1]
    Eval = E(xlast, t + dt)
    Bval = B(xlast, t + dt)
    xprime = xlast
    uprime = ulast
    for _ in range(10):
        uint = 0.5 * (uprime + ulast)
        gammaint = np.sqrt(1.0 + np.linalg.norm(uint) ** 2)
        xprime = xlast + dt * uint / gammaint
        uprime = ulast + dt * (Eval + np.cross(uint, Bval) / gammaint)
    return xprime, uprime


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ExB tests")
    parser.add_argument(
        "--method",
        choices=["euler", "implicit", "rk4", "boris"],
        help="Integration method to use",
    )
    parser.add_argument(
        "--dt", type=float, default=0.01, help="Time step for integration"
    )
    parser.add_argument(
        "--tmax", type=float, default=50.0, help="Maximum time for integration"
    )
    parser.add_argument(
        "--x0",
        type=float,
        nargs=3,
        default=(0.0, 0.0, 0.0),
        help="Initial position: x, y, z",
    )
    parser.add_argument(
        "--u0",
        type=float,
        nargs=3,
        default=(0.05, 0.0, 0.0),
        help="Initial 4-velocity: ux, uy, uz",
    )
    parser.add_argument(
        "--emag", type=float, default=0.01, help="Magnitude of the electric field"
    )
    parser.add_argument(
        "--bmag", type=float, default=1.0, help="Magnitude of the magnetic field"
    )
    parser.add_argument(
        "--size", type=int, nargs=2, default=(80, 30), help="Dimensions of the plot"
    )
    parser.add_argument(
        "--xlim", type=float, nargs=2, default=(-np.inf, np.inf), help="X-axis limits"
    )
    parser.add_argument(
        "--ylim", type=float, nargs=2, default=(-np.inf, np.inf), help="Y-axis limits"
    )
    parser.add_argument("--xaxis", type=str, default="x", help="Quantity for x-axis")
    parser.add_argument("--yaxis", type=str, default="y", help="Quantity for y-axis")
    args = parser.parse_args()

    plt.theme("pro")
    plt.plot_size(*args.size)

    def Efunc(x: np.ndarray, t: float) -> np.ndarray:
        return np.array([0.0, args.emag, 0.0])

    def Bfunc(x: np.ndarray, t: float) -> np.ndarray:
        return np.array([0.0, 0.0, args.bmag])

    integrator = None
    name = None
    x0, u0 = args.x0, args.u0
    if args.method == "euler":
        integrator = euler
        name = "Explicit Euler"
    elif args.method == "implicit":
        integrator = implicit
        name = "Implicit"
    elif args.method == "rk4":
        integrator = rk4
        name = "Runge-Kutta 4th Order"
    elif args.method == "boris":
        integrator = boris
        name = "Boris"

    # roll back half step for leapfrog
    if name in ["Boris"]:
        _, u0 = implicit(
            np.array([args.x0]), np.array([args.u0]), 0.0, -args.dt / 2, Efunc, Bfunc
        )

    if integrator is None:
        raise ValueError("Invalid integration method specified.")

    ts, xs, us = cs.integrate(
        x0=np.array(x0),
        u0=np.array(u0),
        dt=args.dt,
        tmax=args.tmax,
        E=Efunc,
        B=Bfunc,
        integrator=integrator,
    )

    def get_quantity(q: str) -> np.ndarray:
        if q == "x":
            return xs[:, 0]
        elif q == "y":
            return xs[:, 1]
        elif q == "z":
            return xs[:, 2]
        elif q == "ux":
            return us[:, 0]
        elif q == "uy":
            return us[:, 1]
        elif q == "uz":
            return us[:, 2]
        else:
            raise ValueError(f"Unknown quantity: {q}")

    def get_extent(values: np.ndarray) -> tuple[float, float]:
        dv = np.max(values) - np.min(values)
        return np.min(values) - 0.1 * dv, np.max(values) + 0.1 * dv

    xaxis = get_quantity(args.xaxis)
    yaxis = get_quantity(args.yaxis)
    plt.plot(xaxis, yaxis)
    if args.xlim[0] == -np.inf and args.xlim[1] == np.inf:
        plt.xlim(*get_extent(xaxis))
    else:
        plt.xlim(*args.xlim)

    if args.ylim[0] == -np.inf and args.ylim[1] == np.inf:
        plt.ylim(*get_extent(yaxis))
    else:
        plt.ylim(*args.ylim)
    plt.xlabel(args.xaxis)
    plt.ylabel(args.yaxis)
    plt.title(f"{name} dt={args.dt}")
    plt.show()
    plt.cld()
