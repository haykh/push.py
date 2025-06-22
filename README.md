first deploy the virtual python environment and install all the dependencies:
```sh
# either use the built-in script:
./deploy.sh

# or do it manually
python3 -m venv .venv
source .venv/bin/activate
pip install numpy plotext
```

use the script to run the particle pusher with given params and plot:
```sh
python3 exb_test.py [arguments]
```

example:
```sh
python exb-tests.py --method implicit --dt 0.1 --tmax 20.0
```

all the possible arguments are listed below:
```sh
usage: exb-tests.py [-h] [--method {euler,implicit,rk4,boris}] [--dt DT] [--tmax TMAX] [--x0 X0 X0 X0] [--u0 U0 U0 U0] [--emag EMAG] [--bmag BMAG] [--size SIZE SIZE] [--xlim XLIM XLIM] [--ylim YLIM YLIM] [--xaxis XAXIS] [--yaxis YAXIS]

Run ExB tests

options:
  -h, --help            show this help message and exit
  --method {euler,implicit,rk4,boris}
                        Integration method to use
  --dt DT               Time step for integration
  --tmax TMAX           Maximum time for integration
  --x0 X0 X0 X0         Initial position: x, y, z
  --u0 U0 U0 U0         Initial 4-velocity: ux, uy, uz
  --emag EMAG           Magnitude of the electric field
  --bmag BMAG           Magnitude of the magnetic field
  --size SIZE SIZE      Dimensions of the plot
  --xlim XLIM XLIM      X-axis limits
  --ylim YLIM YLIM      Y-axis limits
  --xaxis XAXIS         Quantity for x-axis
  --yaxis YAXIS         Quantity for y-axis
```
