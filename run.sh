#!/usr/bin/env bash

if [ "$1" == "1" ]; then
  python push.py --method euler --preset ExB --u0 0.5 0.0 0.0 --x0 0.0 0.0 0.0 --dt 0.01
elif [ "$1" == "2" ]; then
  python push.py --method euler --preset ExB --u0 0.5 0.0 0.0 --x0 0.0 0.0 0.0 --dt 0.01 --xaxis t --yaxis g
elif [ "$1" == "3" ]; then
  python push.py --method euler --preset ExB --u0 0.5 0.0 0.0 --x0 0.0 0.0 0.0 --dt 0.25 --movie 60
elif [ "$1" == "4" ]; then
  python push.py --method implicit --preset ExB --u0 0.5 0.0 0.0 --x0 0.0 0.0 0.0 --dt 0.25
elif [ "$1" == "5" ]; then
  python push.py --method implicit --preset ExB --u0 0.5 0.0 0.0 --x0 0.0 0.0 0.0 --dt 1.25
elif [ "$1" == "6" ]; then
  python push.py --method rk4 --preset ExB --u0 0.5 0.0 0.0 --x0 0.0 0.0 0.0 --dt 0.25
elif [ "$1" == "7" ]; then
  python push.py --method rk4 --preset ExB --u0 0.5 0.0 0.0 --x0 0.0 0.0 0.0 --dt 1.25
elif [ "$1" == "8" ]; then
  python push.py --method boris --preset ExB --u0 0.5 0.0 0.0 --x0 0.0 0.0 0.0 --dt 0.25
elif [ "$1" == "9" ]; then
  python push.py --method boris --dt 0.2 --u0 50.0 0.0 50.0 --tmax 100 --e 0.0,0.0,0.0 --b 0.0,0.0,1.0 --drag sync --gammarad 10 --xaxis t --yaxis g --yscale log
elif [ "$1" == "10" ]; then
  python push.py --method boris --dt 0.2 --u0 50.0 0.0 50.0 --tmax 100 --e 0.0,0.0,0.0 --b 0.0,0.0,1.0 --drag sync --gammarad 10 --movie 60
elif [ "$1" == "11" ]; then
  python push.py --method boris --dt 0.2 --u0 50.0 0.0 50.0 --tmax 100 --e 0.0,0.0,0.0 --b 0.0,0.0,1.0 --drag sync --gammarad 10 --xaxis t --yaxis uz/g --ylim 0 1
fi
