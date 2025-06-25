# `push.py`

python command-line implementation of popular pusher algorithms for charges in electromagnetic fields. 

### usage

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
python3 push.py [ARGS]
```

### things to try

1. play with the `--method [METHOD]` flag trying different options, and vary the `--dt [TIMESTEP]`
2. use predefined fields with `--preset [PRESET]`, e.g., `--preset mirror`
3. try custom fields with `--e [EFIELD]` and `--b [BFIELD]`, e.g., `--b 0.0,0.0,2.0+np.tanh(y)` (go crazy)
4. play it as a movie(!!!) by adding `--movie [FRAMERATE]`
5. add synchrotron drag to the problem by providing `--drag sync` and `--gammarad [NUMBER]`

### examples

<details>
<summary><pre>python push.py --method boris --dt 5.0 --size 80 40 --u0 20.0 0.0 0.0 --tmax 500 --e 0.0,0.5,0.0 --b 0.0,0.0,1.0</pre></summary>

![demo-1](demos/demo-1.png)

</details>

<details>
<summary><pre>python push.py --method implicit --u0 0.05 0.1 0.0 --tmax 150 --preset mirror --dt 0.05</pre></summary>

![demo-2](demos/demo-2.png)

</details>

<details>
<summary><pre>python push.py --method boris --u0 0.0 1.0 0.0 --x0 0.0 0.0 0.0 --tmax 40 --preset betatron --dt 0.02</pre></summary>

![demo-3](demos/demo-3.png)

</details>

<details>
<summary><pre>python push.py --method rk4 --dt 0.2 --size 80 40 --u0 1.0 0.0 0.0 --tmax 100 --e "0.0,0.0,0.0" --b "0.0,0.0,1/(1+t/20)"</pre></summary>

![demo-4](demos/demo-4.png)

</details>


### all options

all the possible arguments are listed below:
```sh
usage: push.py [-h] [--method {euler,implicit,rk4,boris}] [--drag {sync,ic,None}] [--gammarad GAMMARAD] [--dt DT] [--tmax TMAX] [--x0 X0 X0 X0] [--u0 U0 U0 U0] [--e E] [--b B]
               [--size SIZE SIZE] [--xlim XLIM XLIM] [--ylim YLIM YLIM] [--xscale {log,linear}] [--yscale {log,linear}] [--preset {mirror,ExB,gradB,betatron,None}] [--xaxis XAXIS]
               [--yaxis YAXIS] [--movie MOVIE]

Push particles in electromagnetic fields

options:
  -h, --help            show this help message and exit
  --method {euler,implicit,rk4,boris}
                        Integration method to use
  --drag {sync,ic,None}
                        Radiative drag to impose
  --gammarad GAMMARAD   Radiative drag gamma factor
  --dt DT               Time step for integration
  --tmax TMAX           Maximum time for integration
  --x0 X0 X0 X0         Initial position: x, y, z
  --u0 U0 U0 U0         Initial 4-velocity: ux, uy, uz
  --e E                 Electric field [format: ex,ey,ez; can use expressions with x, y, z, t]
  --b B                 Magnetic field [format: bx,by,bz; can use expressions with x, y, z, t]
  --size SIZE SIZE      Dimensions of the plot
  --xlim XLIM XLIM      X-axis limits
  --ylim YLIM YLIM      Y-axis limits
  --xscale {log,linear}
                        X-axis scale
  --yscale {log,linear}
                        Y-axis scale
  --preset {mirror,ExB,gradB,betatron,None}
                        preset E, B configuration
  --xaxis XAXIS         Quantity for x-axis
  --yaxis YAXIS         Quantity for y-axis
  --movie MOVIE         Play the simulation as a movie with a given framerate
```
