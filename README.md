---

## Compilation

If you want to run in parallel mode, set `#define PARALLEL_MODE 1` in `definitions.h`. Otherwise, set `#define PARALLEL_MODE 0`.
For the parallelization, we use `OpenMP` framework
To compile. On linux:
```
g++ --std=c++17 -g -fopenmp Project08-07_reconstruct.cpp -o idarp
g++ GenerateUserFinal.cpp
```

If you are running on Mac or Windows, check the equivalent compilation command (in particular, check how to use OpenMP).
---

## Run
Example of run command line
```
./idarp.exe Ten_thousand_customers.txt ListeOfLocations_10000_Cust.txt Discretisation_10000_Cust.txt load_file.txt

./idarp.exe one_thousand_cust.txt ListeOfLocations_10000_Cust.txt Discretisation_10000_Cust.txt load_file.txt

./idarp.exe one_hundred_cust.txt ListeOfLocations_10000_Cust.txt Discretisation_10000_Cust.txt load_file.txt

./a.exe
```

### 
ListOfCustomers:	"username;cus_orig_lat;cus_orig_lon;cus_orig_x;cus_orig_y;cus_orig_id;cus_dest_lat;cus_dest_lon;cus_dest_x;cus_dest_y;cus_dest_id;EAT;LAT;location_orig;location_dest";
ListOfCustomers: "username;cus_orig_lat;cus_orig_lon;cus_orig_x;cus_orig_y;cus_orig_id;cus_dest_lat;cus_dest_lon;cus_dest_x;cus_dest_y;cus_dest_id;EAT;LAT;location_orig;location_dest";

Stations_100:	"Station_Number   Station_lat    Station_lon    Station_x    Station_y    Station_x_discr    Station_y_discr    Station_Id";

## Results
Graphic1.txt: this file contains the following information the number of vehicle of global best particle per iteration

Graphic2.txt: Average Travel Time: the Average Travel Time that for each rider of global best particle per iteration

Graphic3.txt: contains the number of vehicles in each line, per each PSO iteration, of the best particle per iteration.

Graphic4.txt: contains the waiting time in each line upon iteration, of the best particle per iteration 

Graphic5.txt: contains the number of active stations in each line upon iterations, of the global best 

Graphic6.txt: contains the number of riders in each rider types upon each iteration, with riders types
	R_W_PT_RS: riders that will use an integrated trasport composed of Walk -> Public Transport -> Ride Sharing
	R_RS_PT_W: riders that will use an integrated trasport composed of Ride Sharing -> Public Transport -> Walk
	R_PT: riders that will use only Public Transport -> Ride Sharing
	R_RS: riders that will use only Ride Sharing
	walk: only walkers 

Graphic7.txt: this file visualizes the average traveling time in PT and RS from a zone to another through iterations

Graphic8.txt: visualizes the number of users of each line (See functions Visualize_Number_of_users_of_each_line_1 and Visualize_Number_of_users_of_each_line_2))

Figure 1 :
1st plot: we get the number of users of each line from Graphic8.txt and the number of ride sharing users from Graphic6 column 4. The torelated extratime 
          can be modified in Maximum_Trip_Coefficient
		  
Figure 2: 1st plot, Graphic6
Figure 3: Graphic16


