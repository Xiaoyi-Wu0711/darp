#define PARALLEL_MODE 0


#define walk_speed 3     //     1.4 m/s     ==>   5.04 km/H
//using namespace std;
#define MYINFINITY 9999
#define MYINFINITY_2 999999 
#define t_ingress 2          //time for the passenger to enter to the platform (min)
#define t_egress 1           //time for the passenger to exit to the platform (min)
#define TS 10        //time spent by the bus to a stop (min)
#define PTV_speed 60   //Public transit speed
#define RSV_speed 30   //Ride sharing vehicle speed
#define Number_stations 63      // Number of stations
#define Number_categories 5 //5 category customers

#define number_of_lines 7      //Number of line
#define ST 2					//the time required to access the platform center to the docks
#define PI 3.14159265358979323846
#define maximum_walk 2.48       // 30 min of walk    htPopulation_sizes://d2wwhj0amomscw.cloudfront.net/temps-parcours/WEB-A4-distances-temps-paris-marche-pied-FOND-IMAGE-V2.pdf
#define Maximum_Trip_Coefficient  1  // tolerate time 1; 1.25; 1.5; 1.75; 2; 3
#define Planification_Horizon 2    // Number of hours in the horizon of planification
#define Circuity 1.255
#define Circuity_walker 1.391
#define number_of_possible_locations 32000

/**************    PSO Parameters   *************/
#define PSO_Numbe_Iteration 10
#define errorCon 0.01           // Termination Condition By Error
#define W 0.9                 // Inertia Coefficient
#define C1 0.6                 // Acceleration Coefficient
#define C2 1.4           // Acceleration Coefficient
#define C_1 2
#define C_2 2

#include <map>
#include<stdio.h>
#include <stdlib.h>

class Consts
{public:

	// Input related to the network:
	static const string stations_file;
	
	//Number of different movement types between zones P-P, P-PC, P-GC..  where P=Paris GC=Grande coronne PC=Petite coronne
	static const int num_of_zonal_movement_types;

	// Input related to the problem:
	static const string result_file;
	static const int max_line_stations;//Maximum of stations on a line

	// Algorithm-related parameters:
	static const int LNS_iterations;
	static const double max_removal_pct;
	static const double max_deterioration_factor;

	// Fleet-related parameters:
	static const int vehicle_cap;
	static const float vehicle_cost;
	static const float circuity_factor;
	static const int vehicle_speed;
	static const int initial_PT_vehicle_number;
	static const int initial_headway;

	static const int max_dissimilarity;
	static const bool debug_mode;
	static const int seed;

	// PSO parameters
	static const int population_size;
};



enum rider_category
{
	INVALID=-1, // Invalid
	R_W_PT_RS=0, // Walk + public transit + ride sharing
	R_RS_PT_W=1, // Ride sharing + public transit + walk
	R_W=2, // Only walking
	R_PT=3, // Walking + public transit
	R_RS=4 // Only ride sharing
};

struct scores
{
	int TNV;
	int VRS;
	int VPT;
	int N;
	double TT;
	float TA[number_of_lines]; //headway time of each line
	double Avg_TT[9];
};


vector<scores>scor;
class line{ //TODO aa: COuld you please add some explanation to
	// these attributes
	public:
		int NV; // # of buses allocated to this line
		vector<int> s;  // stop id 
		vector<int> active; // active status of stops on line 
		// int nbr_stat; // number of stations on each line
		int nbr_bus;
		vector<double> V_1;
		vector<double> V_0;
		float freq;
		int waitingtime;
		int index;
		int Max_NV;
		int Min_NV;
		int Nbr_users;
		bool is_it_equal_to(const line& other) const;
		line(int NV_, vector<int> s_, vector<int> active_);
		void check_consistency() const;
		
		vector<int> get_active_station_indices() const {
			vector<int> res;
			for (int i=0; i<s.size(); i++) {
				if (s.at(i)>=0 && active.at(i)>0)
					res.push_back(i);
			}
			return res;
		};

		vector<int> get_active_station_ids() const {
			vector<int> res;
			for (int i=0; i<s.size(); i++) {
				if (s.at(i)>=0 && active.at(i)>0)
					res.push_back(this->s.at(i));
			}
			return res;
		};

		int get_nbr_stat() const {
			int count = 0;
			for (int i=0; i<this->s.size(); i++) {
				if (this->s.at(i)>=0 && this->active.at(i)>0) 
					count++;
			}
			return count;
		}
};

struct Point{
	float lat;
	float lon;
	int indx; //index of the point
	int x; // It is the discretized latitude
	int y; // It is the discretized longitude
	string name;
	int Identif; //identifier within the list of discretized point
							// It can happen that two different Point with different
							// lat and lon, have the same discretized version, thus
							// having the same Identif
	int location; //It identifier a macro-zone in whcih this point lies.
								//Indeed, we divide the entire surface into multiple
								//macro-zones, i.e., Paris Grande Couronne, Petite Couronne,
								// centre-ville
};

struct PointDistance{
    Point pt;
    float dist;
};

struct stop{
	int station_id; //station id: from 0 to 63
	Point p;
	int EL[number_of_lines];//Set of station's line, if it is -1 means not on the line
};

struct arc{
	int line;
	int nbr_stat;
};

struct chang{
	int vlr;
	int pos;
};

struct Evl{
	double fit_1;
	int fit_2;
};

struct time_transp{
	double RS;
	double PT;
};


struct Heure{
	int Hr=-1;
	int Minut=-1;
};

struct Trans{
	int dx=-1;  // station point index
	float wd=-1;
};

struct Sol_hor{
	int T_o=-1;
	int T_sd=-1;
	int T_d=-1;
	int T_sa=-1;
};

class customer{
	public:
	int indxe = -1;
	Point orig;
	Point dest;
	rider_category RT;        //rider type: R_W_PT_RS:0; R_RS_PT_W:1; R_W:2; R_PT:3; R_RS:4 
	int EDT; // earliest departure time
	int LAT;//latest arrive time
	float dstnce; //distance bewteen orig to dest
	Trans DS;         //Departure station
	Trans AS;         //Aiival Station
	vector<PointDistance> FM; //this vector contains the all possible first miles distances that means the distance between the origin of a customer and all the stations (because not always we take the closest station because it can be desactivated)
	vector<PointDistance> LM; //this vector contains the all possible last miles distances that means the distance between the destination of a customer and all the stations (because not always we take the closest station because it can be desactivated)
	Sol_hor SH;
	float fm; //we initialise the real first mile distance by -1
	float lm; //we initialise the real last mile distance by -1
	vector<string> traj_stat;      // La solution pour chaque population ainsi que le réseau original
	vector<string> traj_line;
	float WT;          // waiting time
	time_transp TTT; //total travel time
	//int trans_point;

	customer()
	{
		RT=INVALID;
		EDT=-1;
		LAT=-1;
		dstnce=-1;
		fm=-1; lm=-1; WT=-1;
		for (int i=0; i<Consts::population_size+1; i++)
		{
			traj_stat.push_back("void");
			traj_line.push_back("void");
		}
	}
};

/**
 * This represents the serving of a station within a certain line
 */
class StopPoint
{
public:
	int station_id=-1; // station id 
	int line_id; // lines' list include this station

	StopPoint(int line_id_, int station_id_)
	{
		this->line_id = line_id_;
		this->station_id = station_id_;
	}
};

struct Path
{
int solut=-1;
int dep=-1;			//station de départ
int ariv=-1;			//station d'arrivée
int cost=-1;			//trajectory t
};

class DiscretizationHandler
{
	public:
		double delta_x=-1;
		double delta_y=-1;
		DiscretizationHandler(const std::string discretization_info_file);
};







