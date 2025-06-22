import numpy as np
from typing import TypeAlias, Callable

Field: TypeAlias = Callable[[np.ndarray, float], np.ndarray]
Integrator: TypeAlias = Callable[
    [np.ndarray, np.ndarray, float, float, Field, Field], tuple[np.ndarray, np.ndarray]
]


def integrate(
    x0: np.ndarray,
    u0: np.ndarray,
    dt: float,
    tmax: float,
    E: Field,
    B: Field,
    integrator: Integrator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xarr = np.array([x0])
    uarr = np.array([u0])
    nsteps = int(tmax / dt)
    tarr = np.linspace(0, tmax, nsteps + 1)
    for t in tarr[1:]:
        xnew, unew = integrator(xarr, uarr, t, dt, E, B)
        xarr = np.vstack((xarr, xnew))
        uarr = np.vstack((uarr, unew))

    return tarr, xarr, uarr
