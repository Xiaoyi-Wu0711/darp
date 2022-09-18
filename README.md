---

## Compilation

If you want to run in parallel mode, set `#define PARALLEL_MODE 1` in `definitions.h`. Otherwise, set `#define PARALLEL_MODE 0`.
For the parallelization, we use `OpenMP` framework
To compile. On linux:
```
g++ --std=c++17 -g -fopenmp Project08-07_reconstruct.cpp -o idarp
```

If you are running on Mac or Windows, check the equivalent compilation command (in particular, check how to use OpenMP).
---

## Run

Example of run command line(the last argument is the maximum tolerating coefficient)
```
./idarp.exe input/one_hundred_cust.txt input/ListeOfLocations_10000_Cust.txt input/Discretisation_10000_Cust.txt input/load_file.txt 1.5

```
run a series of experienments by python:

```
python python_script.py

```


## Results
./result/Graphic1.txt: this file contains the following information the number of vehicle of global best particle per iteration

./result/Graphic2.txt: Average Travel Time: the Average Travel Time that for each rider of global best particle per iteration

./result/Graphic3.txt: contains the number of vehicles in each line, per each PSO iteration, of the best particle per iteration.

./result/Graphic4.txt: contains the waiting time in each line upon iteration, of the best particle per iteration 

./result/Graphic5.txt: contains the number of active stations in each line upon iterations, of the global best 

./result/Graphic6.txt: contains the number of riders in each rider types upon each iteration, with riders types
	R_W_PT_RS: riders that will use an integrated trasport composed of Walk -> Public Transport -> Ride Sharing
	R_RS_PT_W: riders that will use an integrated trasport composed of Ride Sharing -> Public Transport -> Walk
	R_PT: riders that will use only Public Transport -> Ride Sharing
	R_RS: riders that will use only Ride Sharing
	walk: only walkers 

./result/Graphic7.txt: this file visualizes the average traveling time in PT and RS from a zone to another through iterations

./result/Graphic8.txt: visualizes the number of users of each line (See functions Visualize_Number_of_users_of_each_line_1 and Visualize_Number_of_users_of_each_line_2))

./result/Graphic12.txt: record each particle's performance in each iteration

./result/Graphic15.txt: record number of active station number in each line of G_best


./result/Graphic19.txt: contains the number of rider in each line from one zone to another zone

## visualization
./plot.xlsx visualize the results above

./javascript_visualization/index.html visualize the location of each inactive/active stations
