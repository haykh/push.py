import numpy as np
import plotext as plt
import caltechschool as cs
import argparse


#
# Methods
#
def euler(
    x: np.ndarray,
    u: np.ndarray,
    t: float,
    dt: float,
    E: cs.Field,
    B: cs.Field,
    drag: cs.DragTerm,
    gammarad: float,
) -> tuple[np.ndarray, np.ndarray]:
    xlast = x[-1]
    ulast = u[-1]
    Eval = E(xlast, t)
    Bval = B(xlast, t)
    gammalast = np.sqrt(1.0 + np.linalg.norm(ulast) ** 2)
    xnew = xlast + dt * ulast / gammalast
    unew = ulast + dt * (Eval + np.cross(ulast, Bval) / gammalast)
    if drag is not None:
        unew += dt * drag(ulast, Eval, Bval, gammarad)
    return xnew, unew


def rk4(
    x: np.ndarray,
    u: np.ndarray,
    t: float,
    dt: float,
    E: cs.Field,
    B: cs.Field,
    drag: cs.DragTerm,
    gammarad: float,
) -> tuple[np.ndarray, np.ndarray]:
    xlast = x[-1]
    ulast = u[-1]
    Eval = E(xlast, t)
    Bval = B(xlast, t)

    def rhs_x(u: np.ndarray) -> np.ndarray:
        return u / np.sqrt(1.0 + np.linalg.norm(u) ** 2)

    def rhs_u(u: np.ndarray) -> np.ndarray:
        lorentz = Eval + np.cross(u, Bval) / np.sqrt(1.0 + np.linalg.norm(u) ** 2)
        if drag is None:
            return lorentz
        else:
            return lorentz + drag(u, Eval, Bval, gammarad)

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
    x: np.ndarray,
    u: np.ndarray,
    t: float,
    dt: float,
    E: cs.Field,
    B: cs.Field,
    drag: cs.DragTerm,
    gammarad: float,
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

    if drag is not None:
        unew += dt * drag(0.5 * (unew + ulast), Eval, Bval, gammarad)

    xnew = xlast + dt * unew / np.sqrt(1.0 + np.linalg.norm(unew) ** 2)
    return xnew, unew


def implicit(
    x: np.ndarray,
    u: np.ndarray,
    t: float,
    dt: float,
    E: cs.Field,
    B: cs.Field,
    drag: cs.DragTerm,
    gammarad: float,
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
        if drag is not None:
            uprime += dt * drag(uint, Eval, Bval, gammarad)
    return xprime, uprime


#
# Presets
#
def mirror_field(x: np.ndarray, t: float) -> tuple[np.ndarray, np.ndarray]:
    e = np.zeros(3)
    b = np.zeros(3)
    L = 1.0
    D = 0.1
    Bin = 0.5
    Bout = 2.0
    if np.abs(x[0]) >= L + D:
        b[0] = Bout
        b[1] = 0.0
    elif np.abs(x[0]) <= L:
        b[0] = Bin
        b[1] = 0.0
    else:
        b[0] = (
            Bin
            + (Bout - Bin)
            * (np.sin(np.pi * (np.abs(x[0]) - L + 0.5 * D) / D) + 1)
            * 0.5
        )
        phi = np.arctan2(x[2], x[1])
        br = (
            -0.5
            * np.sqrt(x[2] ** 2 + x[1] ** 2)
            * (
                0.5
                * (Bout - Bin)
                * np.cos(np.pi * (np.abs(x[0]) - L + 0.5 * D) / D)
                * (np.pi / D)
                * np.sign(x[0])
            )
        )
        b[1] = br * np.cos(phi)
        b[2] = br * np.sin(phi)
    return e, b


def betatron_field(x: np.ndarray, t: float) -> tuple[np.ndarray, np.ndarray]:
    e = np.array([0.1, 0.0, 0.0])
    b = np.array([0.0, 0.0, np.sign(x[1])])
    return e, b


def ExB_field(x: np.ndarray, t: float) -> tuple[np.ndarray, np.ndarray]:
    e = np.array([0.0, 0.01, 0.0])
    b = np.array([0.0, 0.0, 1.0])
    return e, b


def gradB_field(x: np.ndarray, t: float) -> tuple[np.ndarray, np.ndarray]:
    e = np.array([0.0, 0.0, 0.0])
    b = np.array([0.0, 0.0, 2.0 + np.max([x[0], 0.0])])
    return e, b


#
# Drag terms
#
def synchrotron(
    u: np.ndarray, e: np.ndarray, b: np.ndarray, gammarad: float
) -> np.ndarray:
    gamma = np.sqrt(1 + np.linalg.norm(u) ** 2)
    beta = u / gamma
    kappaR = np.cross(e + np.cross(beta, b), b) + np.dot(beta, e) * e
    chiR2 = np.linalg.norm(e + np.cross(beta, b)) ** 2 - np.dot(beta, e) ** 2
    return (kappaR - gamma**2 * chiR2 * beta) / gammarad**2


def inverseCompton(
    u: np.ndarray, e: np.ndarray, b: np.ndarray, gammarad: float
) -> np.ndarray:
    gamma = np.sqrt(1 + np.linalg.norm(u) ** 2)
    beta = u / gamma
    return -(gamma**2) * beta**2 / gammarad**2


PRESETS = {
    "mirror": {
        "fields": mirror_field,
    },
    "ExB": {
        "fields": ExB_field,
    },
    "gradB": {
        "fields": gradB_field,
    },
    "betatron": {
        "fields": betatron_field,
    },
}

METHODS = {
    "euler": {
        "name": "Explicit Euler",
        "integrator": euler,
    },
    "implicit": {
        "name": "Implicit",
        "integrator": implicit,
    },
    "rk4": {
        "name": "Runge-Kutta 4th Order",
        "integrator": rk4,
    },
    "boris": {
        "name": "Boris",
        "integrator": boris,
    },
}

DRAGS = {
    "sync": synchrotron,
    "ic": inverseCompton,
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Push particles in electromagnetic fields"
    )
    parser.add_argument(
        "--method",
        choices=METHODS.keys(),
        help="Integration method to use",
    )
    parser.add_argument(
        "--drag",
        type=str,
        choices=list(DRAGS.keys()) + [None],
        default=None,
        help="Radiative drag to impose",
    )
    parser.add_argument(
        "--gammarad", type=float, default=np.inf, help="Radiative drag gamma factor"
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
        default=(0.0, 0.0, 0.0),
        help="Initial 4-velocity: ux, uy, uz",
    )
    parser.add_argument(
        "--e",
        type=str,
        default="0,0,0",
        help="Electric field [format: ex,ey,ez; can use expressions with x, y, z, t]",
    )
    parser.add_argument(
        "--b",
        type=str,
        default="0,0,0",
        help="Magnetic field [format: bx,by,bz; can use expressions with x, y, z, t]",
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
    parser.add_argument(
        "--xscale",
        type=str,
        choices=["log", "linear"],
        default="linear",
        help="X-axis scale",
    )
    parser.add_argument(
        "--yscale",
        type=str,
        choices=["log", "linear"],
        default="linear",
        help="Y-axis scale",
    )
    parser.add_argument(
        "--preset",
        type=str,
        choices=list(PRESETS.keys()) + [None],
        default=None,
        help="preset E, B configuration",
    )
    parser.add_argument("--xaxis", type=str, default="x", help="Quantity for x-axis")
    parser.add_argument("--yaxis", type=str, default="y", help="Quantity for y-axis")
    parser.add_argument(
        "--movie",
        type=float,
        default=None,
        help="Play the simulation as a movie with a given framerate",
    )
    args = parser.parse_args()

    def Efunc(x: np.ndarray, t: float) -> np.ndarray:
        if args.preset is not None:
            e, _ = PRESETS[args.preset]["fields"](x, t)
            return e
        e = np.zeros(3)
        for i, comp in enumerate(args.e.split(",")):
            e[i] = eval(comp, {"x": x[0], "y": x[1], "z": x[2], "t": t, "np": np})
        return e

    def Bfunc(x: np.ndarray, t: float) -> np.ndarray:
        if args.preset is not None:
            _, b = PRESETS[args.preset]["fields"](x, t)
            return b
        b = np.zeros(3)
        for i, comp in enumerate(args.b.split(",")):
            b[i] = eval(comp, {"x": x[0], "y": x[1], "z": x[2], "t": t, "np": np})
        return b

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
            np.array([args.x0]),
            np.array([args.u0]),
            0.0,
            -args.dt / 2,
            Efunc,
            Bfunc,
            None,
            np.inf,
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
        drag=DRAGS.get(args.drag, None),
        gammarad=args.gammarad,
    )

    def get_quantity(q: str) -> np.ndarray:
        return eval(
            q,
            {
                "x": xs[:, 0],
                "y": xs[:, 1],
                "z": xs[:, 2],
                "ux": us[:, 0],
                "uy": us[:, 1],
                "uz": us[:, 2],
                "g": np.sqrt(1 + np.linalg.norm(us, axis=1) ** 2),
                "t": ts,
                "np": np,
            },
        )

    def get_extent(values: np.ndarray) -> tuple[float, float]:
        dv = np.max(values) - np.min(values)
        return np.min(values) - 0.1 * dv, np.max(values) + 0.1 * dv

    def print_vector(vec: np.ndarray) -> str:
        return f"[{', '.join(f'{v:.2f}' for v in vec)}]"

    def print_field(field: str) -> str:
        return f"[{field.replace(',', ', ')}]"

    def frame(ti: int):
        print(f"dt={args.dt}, tmax={args.tmax}")
        print(f"x0={print_vector(args.x0)}, u0={print_vector(args.u0)}")
        if args.preset is None:
            print(f"E={print_field(args.e)}, B={print_field(args.b)}")
        else:
            print(f"using preset {args.preset}")

        xaxis = get_quantity(args.xaxis)
        yaxis = get_quantity(args.yaxis)
        if args.xscale == "log":
            xaxis = np.log10(xaxis)
        if args.yscale == "log":
            yaxis = np.log10(yaxis)
        plt.plot(xaxis[:ti], yaxis[:ti])
        if args.xlim[0] == -np.inf and args.xlim[1] == np.inf:
            plt.xlim(*get_extent(xaxis))
        else:
            plt.xlim(*args.xlim)

        if args.ylim[0] == -np.inf and args.ylim[1] == np.inf:
            plt.ylim(*get_extent(yaxis))
        else:
            plt.ylim(*args.ylim)
        plt.xlabel(("log " if args.xscale == "log" else "") + args.xaxis)
        plt.ylabel(("log " if args.yscale == "log" else "") + args.yaxis)
        plt.title(f"{name} dt={args.dt}")

    plt.theme("pro")
    plt.plot_size(*args.size)

    if args.movie is not None:
        delta = 1.0 / args.movie
        stopit = False
        while not stopit:
            for ti in range(len(ts)):
                plt.cld()
                plt.clt()
                frame(ti)
                plt.show()
                plt.sleep(delta)
    else:
        frame(len(ts) - 1)
        plt.show()
