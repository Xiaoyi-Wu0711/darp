/* IDENTIFICATION OF THE PHYSICAL LOCATIONS
--------------------------------------------------------------------------------------------------------------
A. User locations:			0 -> n_cust_locations - 1
B. Depot location:			n_cust_locations
C. Transfer stations:		n_cust_locations + 1 -> n_cust_locations + n_stations
-------------------------------------------------------------------------------------------------------------- */

/* IDENTIFICATION OF THE NODES
--------------------------------------------------------------------------------------------------------------
A. User origins:			0 -> n_ride_sharing_requests - 1
B. User destinations:		n_ride_sharing_requests -> 2 * n_ride_sharing_requests - 1
C. Artificial depot:		2 * n_ride_sharing_requests
D. Transfer nodes:			2 * n_ride_sharing_requests + 1 -> 3 * n_ride_sharing_requests
-------------------------------------------------------------------------------------------------------------- */

#include "definitions.h"




// Structs related to the network:
class Location
{ public:
		double x_coord;
		double y_coord;
		double angle;
		std::string type;

		Location(){
			x_coord=0; y_coord=0; angle=0; type="not-a-location";
		}
};

//include station/line/customer configuration, and time matrix from one stop to another stop
class Network{
	public:	
		int n_cust_locations;
		int n_stations;
		int n_used_locations;
		vector<int> map_key_to_location;
		vector<double> time_matrix;
		vector<Location> locations;
		//constructor
		Network(){
			n_cust_locations=-1; n_stations=-1; n_used_locations=-1;}

		void compute_time_matrix(double circuity_factor, double delta_x, double delta_y);
		void read_data_network(std::string locations_file, std::string stations_file_2);
};

// Structs related to the problem:

/**
 * This is any location that can be visited by riders or drivers.
 * It can be an origin, destination or transfer node.
 * See the occurrences of the word "node" in the paper
 * 	Molenbruch, Y., Braekers, K., Hirsch, P., & Oberscheider, M. (2020). Analyzing the benefits of an integrated mobility system using a matheuristic routing algorithm. European Journal of Operational Research.
 */
class node {
	public:
	int loc; double lower_tw; double upper_tw; int load;
	node(){loc=-1; lower_tw=-1; upper_tw=-1; load=-1;}
};




class request
{
	public:
	int original_id=-1; std::vector<int> pickup_nodes; std::vector<int> delivery_nodes; std::string type="invalid_type";
	void check_consistency() const
	{
			if (Consts::debug_mode)
			{
				if (original_id==-1 || type.compare("invalid_type")==0 )
				{
						std::stringstream err_msg;
						err_msg<<"\nLine "<<__LINE__<<": Error: the request is invalid since "
							<<"original_id="<<original_id
							<<" and type="<<type<<endl;
						cout<<err_msg.str();
						throw err_msg;
				}
			}
	}
};

class Global_solution;


class Problem
{
	public:
		int n_ride_sharing_requests=-1;
		int n_nodes=-1;
		vector<double> dissimilarity;
		vector<node> nodes;
		vector<request> requests;

		// We use const to be sure we do not modify the network
		void compute_user_dissimilarity(const Network &n, int n_ride_sharing_requests);

		void construct_data_problem(const Network &n, const Global_solution& glob_sol, int n_ride_sharing_requests,
							const vector<stop>& set_of_stations
							);
}; // End of class Problem





// Structs related to the solution:
struct route { std::vector<int> route_structure; };


class Solution {
	public:
		double obj_value=-1;
		int vehicles_used=-1;
		int max_vehicles=-1;
		vector<int> selected_pickup;
		vector<int> selected_delivery;
		vector<int> vehicle;
		vector<int> load;
		vector<double> schedule;
		std::vector<int> request_bank;
		vector<route> routes;

		// Solution management:
		void initialize_solution(const Network &n, const Problem &p, const int n_ride_sharing_requests);

		// Copy the solution
		void copy_solution_from(const Network &n, const Problem &p, const Solution &source_solution);
	};

// Initialization of the network:
void initialize_network(std::string locations, std::string stations);
void read_data_network(Network &n, std::string locations, std::string stations);
void compute_time_matrix(Network &n);

// Initialization of the problem:
int run_lns_algorithm(int particle_index,const Network &n, std::string results, int ind_sol,
	int n_ride_sharing_requests, vector<customer>& cus, const vector<stop> set_of_stations,
	const vector<Global_solution>& Pop);


// Removal operators:
void random_request_removal(const Network &n, const Problem &p, Solution &s);
void related_request_removal(const Network &n, const Problem &p, Solution &s);
void worst_request_removal(const Network &n, const Problem &p, Solution &s);
void random_route_removal(const Network &n, const Problem &p, Solution &s);
void remove_request(const Network &n, const Problem &p, Solution &s, int r);

// Insertion operators:
void random_order_best_insertion(Network &n, Problem &p, Solution &s);
void perform_best_insertion(Network &n, Problem &p, Solution &s, int r);

// Feasibility checks:
bool check_load_and_schedule(const Network &n, const Problem &p, Solution &s, int v, int start_pos = 1);
void update_load_and_schedule(const Network &n, const Problem &p, Solution &s, int v, int start_pos = 1);

// Other functions:
void print_solution(Problem &p, Solution &s);



/**
	This solution includes the network structure (lines and frequencies)
	and ride sharing car routing
*/
class Global_solution
{public:
	int nume=-1;    // index of the solution (0 for the base solution and PT+1 for the population solutions
	vector<line> L; 
	float G[number_of_lines][Number_stations][Number_stations];
	float headway_time[number_of_lines];          // we store the AWTs according to the frequencies
	Evl obj;
	int nk=-1;
	chang CH[number_of_lines]; // dynamic array 
	int number_of_stations;
	float travel_time_by_category[Number_categories];
	
	// Number of passengers using ride sharing, either in the first mile, last mile 
	// or for their entire trips. People only walking or only using transit are not
	// counted here.
	int RideSharing_Users=0; 
	int R_w_pt_rs_counter=0; 
	int R_rs_pt_w_counter=0; 	
	int R_pt_counter=0; 	
	int R_rs_counter=0; 


	std::vector<float> PT_Estimated_EDT_LAT_To_Station; // map a customer position to its PT_Estimated_EDT_LAT_To_Station (Earliest Departure Time, Latest Arrival Time)

	std::vector<rider_category> Category; 	


	// Number of passengers using ride sharing, either in the first mile, last mile
	// or for their entire trips. People only walking or only using transit are not
	// counted here.

	int ftn=-1; //fitness within the PSO framework
	vector<vector<float> > DIStance; 	// Matrix of size Number_stations x Number_stations
																		// It is the distance between any pair of stations
																		// The unit is kilometer

	double ATT=-1;          //Compute the average travel time spent at PT+RS
	double CTT_PT=0;      //Compute travel time spent at PT
	double CTT_RS=0;       //Compute travel time spent at RS

	// Below is the route, i.e., a sequence of stations
	list <int> list_of_stations[Number_stations][Number_stations];
	list <int> list_of_lines[Number_stations][Number_stations];
	set <int> set_of_used_lines[Number_stations][Number_stations];

	// The set below contains all the lines used to go from the origin to the destination	double ATT;
	set <int>  v[Number_stations][Number_stations];
	
	time_transp Avg_TT_0[9];    //Avg_TT[0]=P->P;Avg_TT[1]=P->PC;Avg_TT[2]=P->GC;Avg_TT[3]=PC->P;Avg_TT[4]=PC->PC;Avg_TT[5]=PC->GC;Avg_TT[6]=GC->P;Avg_TT[7]=GC->PC;Avg_TT[8]=GC->GC;
	int category_zone[9][Number_categories];
	vector <double> velocity;
	
	
	vector <customer> cus; 
	// cus is an attribute of GlobalSolution. By doing so,
	// every GlobalSolution can be optimized in parallel, each with 
	// its own customer statistics.


	Global_solution(int number_of_stations_, const vector <customer> cus_,
				const vector<line> L_
	);
	
	void check_consistency() const;
	
	void initialisation_avg_TT_and_category_zone();
	void Initialise_velocity();
	
	double compute_travel_times_PT();
	
	double compute_travel_times_RS(const Problem &p, const Solution &s);
	
	/**
	 * Returns 1 if station stn is activated, 0 otherwise
	 */
	int station_active_on_any_line(int stn) const;
	
	int station_active_on_specific_line(int l, int stn) const;
	
	/**
	 * Cette fonnction permet d'identifier le type de customer à partir des deux stations actives les plus proches
	 */
	void Customer_category(int indice_cust);
	
	void Define_customer_type(int indice_customer, const vector<stop>& set_of_stations);
	void compteur_2(int iter);
	
	float end_to_end_time(int index_line, const vector<stop>& set_of_stations) const;
	
	int nbr_Pt_V(int index, int waiting_time, const vector<stop>& set_of_stations) const;
	
	float calculate_frequency(int index_line, const vector<stop>& set_of_stations);


	int Generate_NV(int indice_lin, int initial_headway, const vector<stop>& set_of_stations) const;
	
		
	/**
	 * Updates the average waiting time for the lines within the global solution
	 */
	void update_headway_time(const vector<stop>& set_of_stations);
	
	void dijkstra(const vector<stop>& set_of_stations,const vector<line> LIN_,
			const vector< vector<float> >& travel_time_between_stations);
	
	void DPSO(int nm, const Global_solution& G_best, const vector<Global_solution>& P_best_,
		const vector<stop>& set_of_stations, float CR1, double CR2, double CR3);
		
	void PSO(int indice_ind, const Global_solution& G_best,
		vector<Global_solution>& P_best_, const vector<stop>& set_of_stations);
		
	double Eval(const Network &n, const vector<stop>& set_of_stations, int particle_index);
		
	//double moyenne()const;
		
	int write_ride_sharing_results(const vector<stop>& set_of_stations) const;
	
	int run_lns_algorithm(int particle_index, const Network &n,
			const vector<stop> set_of_stations);
			
	void proba(int Cmp, int nl, int bit, Global_solution& G_best, 
			const vector<Global_solution>& P_best_,int H);
			
	void Compute_Customer_PT_trajectory();

	void Update_nbr_usrs();
		
};//class Global_solution;


void line::check_consistency() const
{
	if (Consts::debug_mode)
	{
		// For the moment, nothing we have to check
	}
}

