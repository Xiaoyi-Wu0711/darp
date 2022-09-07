#include<iostream>
#include<stdio.h>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <sstream>
#include <chrono>
#include <set>
#include <list>
using namespace std; 
using namespace chrono;
#include <algorithm>
#include <fstream>
#include <iterator>
#include <numeric>
#include <random>
#include <cstddef>
#include <valarray>
#include <time.h>
#include <unordered_map>
#include "YvesCalculator.h"

#if PARALLEL_MODE
#include <omp.h> 	//Parallelization is inspired by this great tutorial:
									// https://youtu.be/_1QNzaWPYOE
#endif

#define minHT  2                //minimum headway time
#define maxHT  15              //maximal headway time


#include <assert.h>     /* assert */

using namespace std;


//Avg_TT[0]=P->P;Avg_TT[1]=P->PC;Avg_TT[2]=P->GC;Avg_TT[3]=PC->P;Avg_TT[4]=PC->PC;Avg_TT[5]=PC->GC;Avg_TT[6]=GC->P;Avg_TT[7]=GC->PC;Avg_TT[8]=GC->GC;




#define min(X,Y) ((X)<(Y)?(X):(Y))
#define max(X,Y) ((X)>(Y)?(X):(Y))



/////////////////////////////////////////////////////////////
///////////// CONSTANTS /////////////////////////////////////
/////////////////////////////////////////////////////////////

//All these things are defined in class Consts in definitions.h
//int max_num_customers=10000; aa: I removed it as global variable

// Input related to the network:
const string Consts::stations_file = "./input/Stations_10000.txt"; // there are 63 stations in total
// Input related to the problem:
const string Consts::result_file = "./result/Results_1.txt";

// Algorithm-related parameters:
const int Consts::LNS_iterations = 100;       //10000;
const double Consts::max_removal_pct = 0.10;
const double Consts::max_deterioration_factor = 1.05;

// Fleet-related parameters:
//int max_vehicles = 3000; aa: I deactivated this global variable
const int Consts::vehicle_cap = 5;
const float Consts::vehicle_cost = 1000.0;
const int Consts::initial_PT_vehicle_number = 50;  // initial bus number on each line
const float Consts::circuity_factor = 1.0; // the time used by ridersharing

const int Consts::max_dissimilarity = 10000.0;
const bool Consts::debug_mode = true;
const int Consts::seed = 3;
const int Consts::vehicle_speed = 30;
const int Consts::initial_headway = 2; // initial headway 
const int Consts::num_of_zonal_movement_types = 9;
const int Consts::max_line_stations = 18;

//PSO parameters, particle number
const int Consts::population_size=10;
float HeadwayTime[PSO_Numbe_Iteration][number_of_lines];
int ActivStat[PSO_Numbe_Iteration][number_of_lines]; //change for xiaoyi: change to G_best.L[i].get_nbr_stat()
int TypeRider[PSO_Numbe_Iteration][4];        // We have four types of riders
int TypeRider2[PSO_Numbe_Iteration][4];   
vector<line> LIN;


void load_state(string path, vector<line>& state) {
    string l;
    ifstream readFile(path);
    while (getline(readFile, l)) {

		size_t pos = l.find(';');
		if (pos == -1) break;
		// cout << pos;
		int vehicle_number = std::stoi( l.substr(0, pos) );
		l = l.substr(pos+1);

		// station ids
		pos = l.find(';');
		size_t tmp_pos;
		vector<int> station_ids;
		string tmp_l = l.substr(0, pos)+",";
		while ((tmp_pos = tmp_l.find(',')) != string::npos) {
			station_ids.push_back(std::stoi(tmp_l.substr(0, tmp_pos)));
			tmp_l = tmp_l.substr(tmp_pos+1);
		}
		l = l.substr(pos+1);

		// station states
		pos = l.find(';');
		vector<int> station_states;
		tmp_l = l.substr(0, pos)+",";
		while ((tmp_pos = tmp_l.find(',')) != string::npos) {
			station_states.push_back(std::stoi(tmp_l.substr(0, tmp_pos)));
			tmp_l = tmp_l.substr(tmp_pos+1);
		}
		l = l.substr(pos+1);

		// number of active states
		pos = l.find(';');
		int active_states = std::stoi( l.substr(0, pos) );

		line nl = line(vehicle_number, station_ids, station_states);
		state.push_back(nl);
	}
}

//update set_of_station and generate 63 stations with their station id 
void generation_stops(const vector<line> Line, vector<stop>& set_of_stations){
	if (Consts::debug_mode){
		if (Line.size() != number_of_lines)
		{std::stringstream err_msg;
			err_msg<<"\nLine "<<__LINE__<<": Error: number_of_lines="<<number_of_lines<<
				", Line.size()="<<Line.size()<<endl;
			cout<<err_msg.str();
			throw err_msg.str();
		}
	}

	int station_index, line_index, c, d;

	for(int station_index=0; station_index<Number_stations; station_index++){
		if (Consts::debug_mode){
			if (station_index >= set_of_stations.size()){
				cout<<"\n Line:"<<__LINE__<<": Error: trying to access a="<<station_index<<
					", while set_of_stations.size()="<<set_of_stations.size()<<endl;
				exit(0);
			}
		}
		set_of_stations.at(station_index).station_id=station_index;
		for(line_index=0; line_index<number_of_lines; line_index++)
			set_of_stations.at(station_index).EL[line_index]=-1;
	}

	for(line_index=0; line_index<number_of_lines; line_index++){
		vector<int> active_staion_ids = Line[line_index].get_active_station_ids();
		for(int active_staion_id:Line[line_index].get_active_station_ids()){
			d=0;
			while(set_of_stations.at(active_staion_id).EL[d]!=-1  && d!=number_of_lines)
				d++;
			set_of_stations.at(active_staion_id).EL[d]=line_index;
		}
	}	
	// for(int i=0; i<set_of_stations.size(); i++){
	// 	cout<<"\n i \t"<<i<<endl;
	// 	for(line_index=0; line_index<number_of_lines; line_index++){
	// 		cout<<"--------------"<<endl;
	// 		cout<<"\n line_index \t"<<line_index<<endl;
	// 		cout<<"\n set_of_stations.at(station_index).EL[line_index] \t"<<set_of_stations.at(i).EL[line_index]<<endl;
	// 		cout<<"--------------"<<endl;
	// 	}
	// }
	// exit(0);
}

// Read information on all station locations:
void read_data_stations(int index_network, vector<stop>& set_of_stations, std::string stations_file){

	int station_index=0;
	ifstream infile;
	infile.open(stations_file);
	if(infile.fail()){
		cout << "\n failed to read the stations_file"<<stations_file<<endl;
		exit(EXIT_FAILURE);
	}
	int id; double lon; double lat;; double x_coord; double y_coord; int x_discr; int y_discr; int key_loc_discr;
	while (infile >> id >> lat >> lon >> x_coord >> y_coord >> x_discr >> y_discr >> key_loc_discr){
		set_of_stations[station_index].p.indx=id;
		set_of_stations[station_index].p.lat=lat;
		set_of_stations[station_index].p.lon=lon;
		set_of_stations[station_index].p.x=x_coord;
		set_of_stations[station_index].p.y=y_coord;
		set_of_stations[station_index].p.Identif=key_loc_discr;
		station_index++;
	}
	if (Consts::debug_mode){
		if (set_of_stations.size() != Number_stations){
			cout<<"\nLine "<<__LINE__<<": Error: set_of_stations.size() ="<<set_of_stations.size()<<endl;
			exit(0);
		}
	}
	cout<<"\n station file size is:\t"<<set_of_stations.size();
	infile.close();
}

//read customer loc, demands
vector <customer> read_customers_from_file(const string Liste_Of_Customers_file){

	vector <customer> cus;
	ifstream infile;
	infile.open (Liste_Of_Customers_file);

	if(infile.fail()){
		cout << "\n failed to read the Liste_Of_Customers_file"<<Liste_Of_Customers_file<<endl;
		exit(EXIT_FAILURE);
	}

	float latitude_orig, longitude_orig, latitude_dest, longitude_dest;
	double x_coord_orig, y_coord_orig, Identif_orig, x_coord_dest, y_coord_dest, Identif_dest;
	int Id, EDT, LAT, orig_loc, dest_loc;

	while(infile >> Id >> latitude_orig >> longitude_orig >> x_coord_orig >> y_coord_orig >> Identif_orig >> latitude_dest >> longitude_dest >> x_coord_dest >> y_coord_dest >> Identif_dest >> EDT >> LAT >> orig_loc >> dest_loc ){
		customer Cs;
		Cs.indxe=Id;
		Cs.orig.lat=latitude_orig;
		Cs.orig.lon=longitude_orig;
		Cs.orig.x=x_coord_orig;
		Cs.orig.y=y_coord_orig;
		Cs.orig.Identif=Identif_orig;
		Cs.dest.lat=latitude_dest;
		Cs.dest.lon=longitude_dest;
		Cs.dest.x=x_coord_dest;
		Cs.dest.y=y_coord_dest;
		Cs.dest.Identif=Identif_dest;
		Cs.EDT=EDT;
		Cs.LAT=LAT;
		Cs.orig.location=orig_loc;
		Cs.dest.location=dest_loc;
		cus.push_back(Cs);
		}

	int Number_customers = cus.size();
    cout<<"\n customer size is:\t"<<Number_customers<<endl;
	infile.close();
	return cus;
}

//distance = sin^2((point_1.lat-point_2.lat)*pi/180)/2)+cos(dlat)
//check here: https://www.movable-type.co.uk/scripts/latlong.html
float calculateDistance(const Point point_1, const Point point_2){
	double dlat1=point_1.lat*(PI/180);
	double dlong1=point_1.lon*(PI/180);
	double dlat2=point_2.lat*(PI/180);
	double dlong2=point_2.lon*(PI/180);
	double dLong=(point_1.lon-point_2.lon)*(PI/180);
	double dLat=(point_1.lat-point_2.lat)*(PI/180);
	double aHarv= pow(sin(dLat/2.0),2.0)+cos(dlat1)*cos(dlat2)*pow(sin(dLong/2),2);
	double cHarv=2*atan2(sqrt(aHarv),sqrt(1.0-aHarv));
	float distance=6372*cHarv;        //earth_radius = 6372.795477598
	return distance;
}

//This function compare two distances and sort in ascending distance order
bool DistanceMore(const PointDistance& lhs, const PointDistance& rhs){
	return lhs.dist < rhs.dist;
}

//calculate the distance between each customer's origin/destionation and each station 
void Customer_LastMile_FirstMile (vector <customer> &cus, const vector<stop>& set_of_stations){
	int customer_index;
	for(customer_index=0; customer_index<cus.size(); customer_index++){
		// calculate the first mile distance
		for (int station_index = 0; station_index < Number_stations; station_index++){
			PointDistance pd;
			pd.pt = set_of_stations[station_index].p;
		    pd.dist = calculateDistance(cus[customer_index].orig, pd.pt)*Circuity_walker;
		    cus[customer_index].FM.push_back(pd);
		}
		std::sort(cus[customer_index].FM.begin(), cus[customer_index].FM.end(), DistanceMore);

        // cout<<cus[customer_index].FM[0].pt.indx<<endl;
		// cout<<cus[customer_index].FM[62].pt.indx<<endl;
		// exit(0);

        //calculate the last mile distance
		for (int station_index = 0; station_index < Number_stations; station_index++){
			PointDistance pd2;
			pd2.pt = set_of_stations[station_index].p;
			pd2.dist = calculateDistance(cus[customer_index].dest, pd2.pt)*Circuity_walker;
			cus[customer_index].LM.push_back(pd2);
		}
		std::sort(cus[customer_index].LM.begin(), cus[customer_index].LM.end(), DistanceMore);
		
		if (Consts::debug_mode){
			if (cus[customer_index].LM.size() != Number_stations ||cus[customer_index].FM.size() != Number_stations  ){
				cout<<"\nLine "<<__LINE__<<": Error: cus[customer_index].FM="<<cus[customer_index].FM.size()<<
					" or cus[customer_index].LM="<<cus[customer_index].LM.size()<<endl;
				exit(0);
			}
		}

		// cout<<cus[customer_index].LM[0].pt.indx<<endl;
		// cout<<cus[customer_index].LM[62].pt.indx<<endl;
		// exit(0);

	}
}


// initialize the network with 100 customer and 63 stops+1 depot locations
// customer_test.txt ListeOfLocations_10000_Cust.txt
void Network::read_data_network(std::string locations_file, std::string stations_file){

	ifstream infile;
	infile.open(locations_file);
	if(infile.fail()){
		cout << "\n failed to read the locations_file"<<locations_file<<endl;
		exit(EXIT_FAILURE);
	}
	
	//ListeofLocatios_10000_Cust with 1714 customer locations and 63 station locations 
	infile >> n_cust_locations >> n_stations;

	cout<<"n_cust_locations "<<n_cust_locations<<endl;
	cout<<"n_stations "<< n_stations<<endl;

	n_used_locations = n_cust_locations + n_stations + 1;

	for (int i=0; i<n_used_locations+1; i++){
		Location loc;
		locations.push_back(loc);
	}

	// This array maps the keys of the locations to a location id in this code:
	for (int i=0; i< number_of_possible_locations; i++)
		map_key_to_location.push_back(-1);

	int key; double x_coord; double y_coord;
	int counter = 0;

	// Read information on all potential customer locations:
	while (infile >> key >> x_coord >> y_coord){
		map_key_to_location.at(key) = counter;
		Location loc;
		loc.x_coord = x_coord;
		loc.y_coord = y_coord;
		loc.type = "cust_location";
		locations.at(counter)=loc;
		counter += 1;
		// cout<<"\n counter_2 \t"<<counter<<endl;
	}
	infile.close();



	// Depot location:
	int depot_index = n_cust_locations;

	Location depot_loc;
	depot_loc.x_coord = 0.0;
	depot_loc.y_coord = 0.0;
	depot_loc.type = "depot";
	locations.at(counter)=depot_loc;
	assert(locations.at(n_cust_locations).type.compare("depot")==0 );//verification

	//stations_file with loc of 63 stations
	infile.open(stations_file);
	if(infile.fail()){
		cout << "failed to read the stations_file!";
		exit(EXIT_FAILURE);
	}
	int id; double lon; double lat; int x_discr; int y_discr; int key_loc_discr;
	while (infile >> id >> lat >> lon >> x_coord >> y_coord >> x_discr >> y_discr >> key_loc_discr){
		int stat_idx = id + n_cust_locations + 1;
		Location station_loc;
		station_loc.x_coord = x_coord;
		station_loc.y_coord = y_coord;
		station_loc.type = "station";
		locations.at(stat_idx)=station_loc;
		counter += 1;
		map_key_to_location.at(key_loc_discr) = counter;
	}
	infile.close();
	cout<<"\n  map_key_to_location.at(0+100+1) \t"<< map_key_to_location.at(0+1715+1)<<endl;

	// cout<<"\n number_of_possible_locations= "<<number_of_possible_locations<<endl;
	// cout<<"\n map_key_to_location.at(key_loc_discr)= "<<map_key_to_location.at(16427)<<endl;
}

// calculate the travel time between n_used_location in the network(including 63 stations, customer, i depot) by private car
void Network::compute_time_matrix(double circuity_factor, double delta_x, double delta_y){
	//todo for xiaoyi
	for (int loc1 = 0; loc1 < n_used_locations; loc1++)
		for (int loc2 = 0; loc2 <= n_used_locations; loc2++)
			time_matrix.push_back(0);

	// Compute all travel times based on the coordinates:
	for (int loc1 = 0; loc1 < n_used_locations; loc1++){
		for (int loc2 = 0; loc2 <= n_used_locations; loc2++){

			if (Consts::debug_mode){
				if (loc1 >= locations.size() || loc2>=locations.size()){
					cout<<"\n Line "<<__LINE__<<": Try to access loc1="<<loc1<<
						" and loc2="<<loc2<<", but location.size()="<<locations.size()
						<<endl;
					exit(0);
				}
			}

			double travel_time =
					sqrt(pow(delta_x*(locations.at(loc1).x_coord - locations.at(loc2).x_coord),2) + pow(delta_y*(locations.at(loc1).y_coord - locations.at(loc2).y_coord),2));
			travel_time = travel_time * circuity_factor*60*1./Consts::vehicle_speed;		// Note: the assumed speed is now 60 km/u
			// cout<<"\n Line "<<__LINE__<<": travel_time1="<<travel_time<<endl;
			time_matrix.at(loc1 * n_used_locations + loc2) = travel_time;
			time_matrix.at(loc2 * n_used_locations + loc1) = travel_time;
			// cout<<"time_matrx\t"<<time_matrix.size()<<endl;
			// exit(0);
		
		}
	}

	int depot = n_cust_locations;

	// Zero travel time on all links to and from the artificial depot:
	for (int loc = 0; loc < n_used_locations; loc++){
		time_matrix.at(loc * n_used_locations + depot) = 0.0;
		time_matrix.at(depot * n_used_locations + loc) = 0.0;
	}

}; //Network::compute_time_matrix

////include station/line/customer configuration, and time matrix from one stop to another stop
void initialize_network(Network &n, const std::string locations_file, std::string stations_file, const DiscretizationHandler& discretizationHandler){	
	n.read_data_network(locations_file, stations_file);
	n.compute_time_matrix(Consts::circuity_factor, discretizationHandler.delta_x, discretizationHandler.delta_y);
	}

//travel time between 2 active stations(next to each other) on the same line, according to middle part of the equation 1
//if 2 stations are not on the same line or one of them are deactive, it would be infinity
vector< vector<float> > initialisation_of_travel_time_on_arc(const vector<stop>& set_of_stations){
	
	assert(Number_stations == set_of_stations.size() );
	vector< vector<float> > Travel_time_on_arc;

	for(int s_1=0; s_1<Number_stations; s_1++){
		vector<float> vec;
		for(int s_2=0; s_2<Number_stations; s_2++)
			vec.push_back(MYINFINITY); //aa: This replaces J[s][cn]=MYINFINITY;
		Travel_time_on_arc.push_back(vec);
	}

	for(int l=0; l<number_of_lines; l++){
		vector<int> active_station_ids;
		for (int i:LIN[l].get_active_station_ids()){
			active_station_ids.push_back(i);
		}

		for (int active_station_index=0; active_station_index<active_station_ids.size()-1;active_station_index++) {

			int active_station_id = active_station_ids.at(active_station_index);
			int next_active_station_id = active_station_ids.at(active_station_index+1);

			const Point p1 = set_of_stations.at(active_station_id).p;
			const Point p2 = set_of_stations.at(next_active_station_id).p;
			
			Travel_time_on_arc.at(active_station_id).at(next_active_station_id) =
					(calculateDistance(p1, p2 )/PTV_speed)*60;
			Travel_time_on_arc.at(next_active_station_id).at(active_station_id)=
					Travel_time_on_arc.at(active_station_id).at(next_active_station_id);
			// } else break; 
			// cout<<Travel_time_on_arc[1][2]<<endl;
			// cout<<Travel_time_on_arc[1][3]<<endl;

			// exit(0);

		}
	}
	return Travel_time_on_arc;
}//initialisation_of_travel_time_on_arc


void Global_solution::check_consistency() const
{
		if (Consts::debug_mode)
		{
			if (number_of_stations < 2)
			{
				std::stringstream err_msg;
				err_msg<<"\nLine "<<__LINE__<<": Error: number_of_stations="<<number_of_stations
					<<". They seem to be too few. Are you sure it is not an error?"<<endl;
				cout<<err_msg.str();
				throw err_msg.str();
			}
					
			
			if (L.size() != number_of_lines)
			{
				std::stringstream err_msg;
				err_msg<<"\nLine "<<__LINE__<<": Error: lines.size()="<<L.size()
					<<" and number_of_liness="<<number_of_lines<<": they should instead be "
					<<"the same"<<endl;
				cout<<err_msg.str();
				throw err_msg;
			}

			if (number_of_stations != Number_stations)
			{
				std::stringstream err_msg;
				err_msg<<"\nLine "<<__LINE__<<": Error: Number_stations="<<Number_stations
					<<" and number_of_stations="<<number_of_stations<<": they should instead be "
					<<"the same"<<endl;
				cout<<err_msg.str();
				throw err_msg;
			}

			if (DIStance.size()!=number_of_stations)
			{
				std::stringstream err_msg;
				err_msg<<"\nLine "<<__LINE__<<": Error: DIStance.size()="<<DIStance.size()
					<<" and number_of_stations="<<number_of_stations<<": they should instead be "
					<<"the same"<<endl;
				cout<<err_msg.str();
				throw err_msg;
			}

			for (int i=0; i<DIStance.size(); i++)
			{
				vector<float> row = DIStance.at(i);
				if (row.size()!=number_of_stations)
				{
					std::stringstream err_msg;
					err_msg<<"\nLine "<<__LINE__<<": Error: When we take the row corresponding to all the "
						<<"distances from station "<<i<<", we see that row.size()="<<row.size()
						<<" and number_of_stations="<<number_of_stations<<": they should instead be "
						<<"the same"<<endl;
					cout<<err_msg.str();
					throw err_msg;
				}
			}

			/*		for (int i=0; i<cus.size(); i++)
					{   cout<<" category("<<i<<")="<<this->Category.at(i)<<"    RT="<<this->cus.at(i).RT;
						if (this->cus.at(i).RT != this->Category.at(i) )
						{
							std::stringstream err_msg;
							err_msg<<"\n Line "<<__LINE__<<": cus.at("<<i<<").RT="<<this->cus.at(i).RT
								<<" Category(i)="<< this->Category.at(i)<<". They should instead be the same"<<endl;
							cout<<err_msg.str();
							throw err_msg;
						}

					}

					}*/

					//TODO aa: pour moi meme if the error above is never raised, then we can remove Category
		}
} //Global_solution::check_consistency()


//initial velocity for each line is 0
void Global_solution::Initialise_velocity(){
	this->velocity.resize(number_of_lines, 0);
	this->check_consistency();
}

// This function return the travel time from the first active station to last active station on line
float Global_solution::end_to_end_time(int line_index, const vector<stop>& set_of_stations) const{

    float cum=0;
	float cum_2=0;

	vector<int> station_ids = this->L.at(line_index).s;

	vector<int> active_station_ids;
	
	for(int i:this->L.at(line_index).get_active_station_ids()){
		active_station_ids.push_back(i);
		cout<<"\n i \t"<<i<<endl;
	}

	int active_nrb=L.at(line_index).get_nbr_stat();

	cout<<"\n active_nrb \t"<<active_nrb<<endl;
	if(active_nrb ==0 || active_nrb ==1 ){
		if(active_nrb ==1){
			cout<<"\n active_nrb==1 line_index\t"<<line_index;
		}
		cum = MYINFINITY; 
	}else{

		int first_active_station_id = active_station_ids[0];
		int first_active_station_index = 0;
		cout<<" \n first_active_station_id \t"<<first_active_station_id<<endl;
		for(int station_index = 0; station_index < station_ids.size(); station_index++){
			if (station_ids[station_index] == first_active_station_id){
				first_active_station_index = station_index; //find the index of first active station in LIN[i].s
				break;
			}
		}
		cout<<" \n first_active_station_index \t"<<first_active_station_index<<endl;

		int last_active_station_id =  active_station_ids.back();
		int last_active_station_index = 0;
		cout<<" \n last_station_id \t"<<last_active_station_id<<endl;
		for(int station_index = 0; station_index < station_ids.size(); station_index++){
			if (station_ids[station_index] == last_active_station_id){
				last_active_station_index = station_index; //find the index of last station in LIN[i].s
				break;
			}
		}
		cout<<" \n last_active_station_index \t"<<last_active_station_index<<endl;

		vector<int> end_to_end_vec;


		end_to_end_vec = {station_ids.begin()+first_active_station_index, station_ids.begin()+last_active_station_index+1}; // get the subvector started from the frist active station to the end, it is a mixture of active and deacitve stations
		// cout<<" \n line \t"<<__LINE__<<endl;
		// cout<<"\n end_to_end_vec.size()"<<end_to_end_vec.size()<<endl;;
		// cout<<"\n end_to_end_vec.size() line_index \t"<<line_index<<endl;

		// for(int tmp =0; tmp<end_to_end_vec.size(); tmp++){
		// 	cout<<"\n end_to_end_vec[tmp] \t"<<end_to_end_vec[tmp]<<endl;
		// }
		for(int station_index = 0; station_index < end_to_end_vec.size()-1; station_index++){
			// cout<<"---------------"<<endl;
			int station_u = end_to_end_vec[station_index];
			int station_v = end_to_end_vec[station_index+1];
			const Point p1 = set_of_stations.at(station_u).p;
			const Point p2 = set_of_stations.at(station_v).p;
			cum+= (calculateDistance(p1,p2)/PTV_speed)*60; 
			// cout<<"\n (calculateDistance(p1,p2)/PTV_speed)*60 \t"<<calculateDistance(p1,p2)/PTV_speed*60<<endl;
			if(std::count(active_station_ids.begin(), active_station_ids.end(), station_u)) //if the depart station(station u on Arc) is active, add one as accelerating time 
				cum+=1;
			if(std::count(active_station_ids.begin(), active_station_ids.end(), station_v))//if the following station(station v on Arc) is active, add one as decelerating time 
				cum+=1;
			// cout<<"\n cum \t"<<cum<<endl;
			// cout<<"---------------"<<endl;

		}
	}
	
    this->check_consistency();
	cout<<"\n cum \t"<<cum<<endl;
	// exit(0);

    return cum;
}


/****    Update freq of the lines within the global solution according to equation 1, freq = N_l/(2t_{end to end time}) ****/
float Global_solution::calculate_frequency(int line_index, const vector<stop>& set_of_stations){	
	cout<<"\n line_index\t"<<line_index<<endl;
	cout<<"\n line \t"<<__LINE__<<endl;
	cout<<"\n L[line_index].NV \t"<<L[line_index].NV<<endl;
	if(this->end_to_end_time(line_index, set_of_stations)>9998){
		this->L[line_index].freq=1/MYINFINITY;
	}else{
		this->L[line_index].freq=this->L[line_index].NV/(2*this->end_to_end_time(line_index, set_of_stations));
		cout<<"\n line \t"<<__LINE__<<endl;
	}
	if(this->L[line_index].freq==0)
		this->L[line_index].freq=1/MYINFINITY;
	this->check_consistency();
	cout<<"\n line \t"<<__LINE__<<endl;
	cout<<"\n L[line_index].freq \t"<<L[line_index].freq<<endl;
	return this->L[line_index].freq;
}

/****  calculate the vehicle number on each line by taking the desired headyway time as input, 
 		according to equation 8: NV= 2*end_to_end_time/ headway_time ****/  
int Global_solution::Generate_NV(int line_index, int headway_time_for_each_line, const vector<stop>& set_of_stations) const{
	this->check_consistency();
	float bus_number;
	bus_number = 2* this->end_to_end_time(line_index, set_of_stations)/headway_time_for_each_line;
	int  bus_number_on_line  = floor(bus_number)+1;
	cout<<"  bus_number_on_line   "<<bus_number_on_line<<endl;
	// exit(0);
	// //check for xiaoyi
	// if(this->headway_time[line_index] > 9998 || end_to_end_time(line_index, set_of_stations) > 9998){
	// 		cout<<"\n this->L.at(line_index).get_nbr_stat() \t"<<this->L.at(line_index).get_nbr_stat()<<endl;
	// 		cout<<"\n line index \t"<<line_index<<endl;
	// 		cout<<"\n this->end_to_end_time(line_index, set_of_stations) \t"<<this->end_to_end_time(line_index, set_of_stations) <<endl;
	// 		cout<<"\n L[line_index].freq \t"<<this->L[line_index].freq<<endl;
	// 		cout<<"\n headway_time[line_index] \t"<<this->headway_time[line_index]<<endl;
	// 		cout<<"\n bus_number_on_each_line \t"<<bus_number_on_line<<endl;
	// 		exit(0);
	// }
	return bus_number_on_line;
}




//update headway time 
void Global_solution::update_headway_time(const vector<stop>& set_of_stations){
	for(int line_index=0; line_index<number_of_lines; line_index++){
		cout<<"\n line_index \t"<<line_index<<endl;
		float freq =this->calculate_frequency(line_index, set_of_stations);
		if (freq ==0){
			headway_time[line_index]=MYINFINITY;
			// cout<<"\n this->L.at(line_index).get_nbr_stat() \t"<<this->L.at(line_index).get_nbr_stat()<<endl;
			// cout<<"\n line index \t"<<line_index<<endl;
			// cout<<"\n this->end_to_end_time(line_index, set_of_stations) \t"<<this->end_to_end_time(line_index, set_of_stations) <<endl;
			// this->calculate_frequency(line_index,set_of_stations);
			// cout<<"\n L[line_index].freq \t"<<L[line_index].freq <<endl;
			// cout<<"\n headway_time[line_index] \t"<<headway_time[line_index]<<endl;
			// exit(0);
		}else{
			headway_time[line_index]=1/freq;
		}
		cout<<"\n headway_time[line_index] \t"<<headway_time[line_index]<<endl;
	}
	this->check_consistency();
}

/**** This function checks if a station is active on line l ****/
int Global_solution::station_active_on_specific_line(int line_index, int station_index) const{
	int active_status;
	vector<int> active_station_ids = L[line_index].get_active_station_ids();
	if(std::count(active_station_ids.begin(), active_station_ids.end(), station_index))
		return active_status=1;
	this->check_consistency();
	return active_status=0;
}

bool contains(vector<int> vec, const int & elem)
{
    bool result = false;
    if( find(vec.begin(), vec.end(), elem) != vec.end() )
    {
        result = true;
    }
    return result;
}

/*
	This function computes the shortest time path between any pair of nodes. The results are written into the following vectors:
	In this->DIStance.at(s).at(t) we write the average time along the shortest time path to go from station s to station t
	In this->list_of_stations[s][t] we write the sequence of stations visited along the shortest time path from s to t
	In this->list_of_lines[s][t]  we write the sequence of lines visited along the shortest time path from s to t
*/
//xiaoyi: remove LIN_ in dijkstra
void Global_solution::dijkstra(const vector<stop>& set_of_stations, const vector<line> LIN_, const vector< vector<float> >& travel_time_on_arc){
	this->check_consistency();
	int i, j, line_index, stop_index_1, stop_index_2;
	vector<StopPoint> V; //a vector of active stops
	vector< list <int> > succ;	//succ: all possible following stations
	vector < vector < float > > disPopulation_sizeT; //the travel time from one stop to get to another stop
	vector<list <int> > path;
	path.clear();

	for(line_index=0; line_index<number_of_lines; line_index++){// xiaoyi: change here find all active and add to V
		line current_line = this->L.at(line_index);
		cout<<"\n current_line.get_nbr_stat()\t"<<current_line.get_nbr_stat()<<endl;
		if (current_line.get_nbr_stat() ==0 || current_line.get_nbr_stat() ==1){
			cout<<"\n get_nbr_stat() ==0 current_line \t"<< line_index ;
			// exit(0);
		}
		if (current_line.get_nbr_stat() >1){
			for (int active_station_id: current_line.get_active_station_ids()){
				StopPoint stop_point(line_index, active_station_id);
				V.push_back(stop_point);
			}
		}
	}	



	//find the possible following stations for each node
	for(stop_index_1=0; stop_index_1<V.size(); stop_index_1++){
		list <int> suc;
		for(stop_index_2=0; stop_index_2<V.size(); stop_index_2++)
			if(	
				(travel_time_on_arc.at(V.at(stop_index_1).station_id).at(V.at(stop_index_2).station_id) != MYINFINITY  &&  
					V.at(stop_index_1).line_id==V.at(stop_index_2).line_id) ||  // stop_index_2 is the node on the same Arc as stop_index_1
				(V.at(stop_index_1).station_id == V.at(stop_index_2).station_id &&  
					V.at(stop_index_1).line_id!=V.at(stop_index_2).line_id) // stop_index_2 is stop_index_1's interchange node
			){
				suc.push_back(stop_index_2);
			}
		succ.push_back(suc);
	}

	this->update_headway_time(set_of_stations);

	for(i=0; i<V.size(); i++){
		vector<float>Ve; //travel time beween any 2 stations from i to j
		for(j=0; j<V.size(); j++){
			if(V[i].station_id==V[j].station_id && V[i].line_id!=V[j].line_id){ //if the 2 stations have same station id but on different line, they are the interchange stop
				Ve.push_back(this->headway_time[V[j].line_id]); //the travel time is the headway time
			}else{
				if(V.at(i).line_id == V.at(j).line_id){//if stations are on same line
					Ve.push_back(travel_time_on_arc.at(V.at(i).station_id).at(V.at(j).station_id)); //the travel time is the travel time on Arc
				}else
					Ve.push_back(MYINFINITY);
			}
		}
	disPopulation_sizeT.push_back(Ve);
	}


	// s:departing station
	for(int s=0; s<Number_stations; s++){
		cout<<"\n s \t"<<s<<endl;
		vector<float> d; // label for each node of the graph, the shortest distance between source to any pending node 
		vector<int> pred; // predecessor for each node of the road network graph (in the shortest paths)
		set<int> pending; // set of pending nodes, one hop candidate set 
		int u;
		float best,travel_time=0;
		d.resize(V.size());
		pred.resize(V.size());
		
		for (int i = 0; i < V.size(); i++)
			d[i] = MYINFINITY;  					//Initiatization to the MYINFINITY

		for(int cn=0; cn<Number_stations;cn++)     // arrive at cn
			this->DIStance.at(s).at(cn)=MYINFINITY; // Matrix of size Number_stations x Number_stations, It is the distance between any pair of stations(km)

		for (int i = 0; i < V.size(); i++)         //Not yet visited node
		{	
			if (V[i].station_id == s){ //the starting station id 
				d[i] = t_ingress + this->headway_time[V[i].line_id]; //the time to go to transit and headway waiting time
				if(d[i]>MYINFINITY){
					cout<<"\n i \t"<<i<<endl;
					cout<<"\n V[i].station_id \t"<<V[i].station_id<<endl;
					cout<<"\n V[i].line_id \t"<<V[i].line_id<<endl;
					cout<<"\n d[i] \t"<<d[i]<<endl;
					cout<<"\n Line \t"<<__LINE__<<endl;
					cout<<"\n warning: d[i]>MYINFINITY!!!"<<endl;
					exit(1);
					// cout<<"end_to_end_time(V[i].line_id, set_of_stations)"<<this->end_to_end_time(V[i].line_id, set_of_stations)<<endl;
				}
				pending.insert(i);
			}
		}

		int k=0,dep;
		while (!pending.empty()){
			// search the node u with minimal value d[u] in the pending set
			best = MYINFINITY;
			u = -1;
			for (set<int>::iterator iter=pending.begin(); iter!=pending.end(); iter++){
				if (d[*iter] < best){
					best = d[*iter];
					u = *iter;
				}
			}
			// if( u ==-1){
			// 	for (set<int>::iterator iter=pending.begin(); iter!=pending.end(); iter++){
			// 		cout<<"\n iter \t"<<*iter<<endl;
			// 		cout<<"\n pending.size() \t"<<pending.size()<<endl;
			// 		cout<<"\n d[*iter] \t "<<d[*iter]<<endl;
			// 		cout<<"\n u==-1 line_index"<<V[*iter].line_id<<endl;

			// 		cout<<"\n u==-1 LIN_.at(line_index).get_nbr_stat() \t"<< LIN_.at(V[*iter].line_id).get_nbr_stat()<<endl;
			// 		cout<<"this->L.at(V[*iter].line_id).get_nbr_stat()\t"<<this->L.at(V[*iter].line_id).get_nbr_stat()<<endl;

			// 		cout<<"\n u==-1 LIN_.at(V[i].line_id).freq \t"<< LIN_.at(V[*iter].line_id).freq<<endl;
			// 		cout<<"this->LIN.at(V[*iter].line_id).freq \t"<<this->L.at(V[*iter].line_id).freq<<endl;

			// 		cout<<"end_to_end_time(V[i].line_id, set_of_stations)"<<this->end_to_end_time(V[i].line_id, set_of_stations)<<endl;

			// 	}
			// 	exit(1);
			// }

			// retrieve u from the pending set and update the value of the sucessors of u
			// when the value is updated, the node is added to the set (no impact if it was already in the set)
			if(k==0) 
				dep=u;

			pending.erase(u);
			k++;

			for (list<int>::iterator iter = succ[u].begin(); iter != succ[u].end(); iter++){
				// cout<<"\n line \t"<< __LINE__<<endl;
				// cout<<"\n u \t"<< u<<"\n line\t"<<__LINE__<<endl;
				// cout<<"\n *iter \t"<<*iter<<"\n line\t"<<__LINE__<<endl;

				if(V[u].line_id == V[*iter].line_id){


					float travel_time = d[u];
					cout<<"\n line \t"<< __LINE__<<endl;

					if ( this->station_active_on_specific_line( V[*iter].line_id, V[*iter].station_id) == 1)//todo for xiaoyi: this constraint will always true,since all station in V is active
						travel_time = travel_time + TS;  //TS: time spent by the bus to a stop (min)
									// cout<<"\n line \t"<< __LINE__<<endl;

					travel_time = travel_time + disPopulation_sizeT[u][*iter];
									// cout<<"\n line \t"<< __LINE__<<endl;

					if (d[*iter] > travel_time){
										// cout<<"\n line \t"<< __LINE__<<endl;

						d[*iter] = travel_time;
						// pred[*iter] = u;
					
					if(*iter<pred.size()){
						pred[*iter] = u;
							}else{
								for(int tmp=0; tmp<pred.size(); tmp++){
									cout<<"\n line\t"<<__LINE__<<endl;
									cout<<"\n pred[tmp]_1 \t"<<pred[tmp]<<endl;
								}
								cout<<"\n pred.size() \t"<<pred.size()<<endl;
								cout<<"\n i \t"<<i<<endl;
								cout<<"\n iter \t"<<*iter<<endl;
								cout<<"\n v.size() \t"<<V.size()<<endl;
								cout<<"\n u_1 \t"<<u<<endl;
								exit(1);
							}
						pending.insert(*iter); ///no risk of duplication: set does not generate duplicates
						// cout<<"\n line \t"<< __LINE__<<endl;

					}
				}else{ //case transfer from one line to other line 
					if (	
						V[u].station_id==V[*iter].station_id
						&& this->station_active_on_specific_line(V[u].line_id, V[u].station_id) == 1
						&& station_active_on_specific_line(V[*iter].line_id, V[*iter].station_id) == 1
					) //both lines stop at this station
					{
						// cout<<"\n line \t"<< __LINE__<<endl;

						if (d[*iter] > d[u] + disPopulation_sizeT[u][*iter]){
											cout<<"\n line \t"<< __LINE__<<endl;

							d[*iter] = d[u] + disPopulation_sizeT[u][*iter];
							pred[*iter] = u;
							if(*iter<pred.size()){
								pred[*iter] = u;
							}else{
								for(int tmp=0; tmp<pred.size(); tmp++){
									cout<<"\n line\t"<<__LINE__<<endl;
									cout<<"\n pred[tmp]_2 \t"<<pred[tmp]<<endl;
								}
								cout<<"\n pred.size() \t"<<pred.size()<<endl;
								cout<<"\n i \t"<<i<<endl;
								cout<<"\n iter \t"<<*iter<<endl;
								cout<<"\n v.size() \t"<<V.size()<<endl;
								cout<<"\n u_2 \t"<<u<<endl;
								exit(1);
							}
						// try{
						// 	pred[*iter] = u;
						// 	// throw 100;
						// }
						// catch(std::exception& e){
						// 	for(int i=0; i<pred.size(); i++){
						// 		cout<<"\n pred[i]_2 \t"<<pred[i]<<endl;
						// 	}
						// 	cout<<"\n u_2 \t"<<u<<endl;
						// }
							pending.insert(*iter); //no risk of duplication: set does not generate duplicates
						}
					}
				}
			}
		}

					cout<<"\n line \t"<< __LINE__<<endl;

		//cout<<"\n d[12]="<<d[12];end
		// V:set of public transit nodes, i.e., pair (station,line)
		for (int i = 0; i < V.size(); i++)	//for (int i = 0; i < V.size(); i++)
		{	


				if(V[i].station_id!=s){   
					list<int> P;
					int t = V[i].station_id;
					//cout<<"\n"<<DIStance[s][t];
					if (d[i] + t_egress < this->DIStance.at(s).at(t) ) // distance from s to t (the result of algo ; ‡ initialization ‡ infinity)
					{   
						this->DIStance.at(s).at(t) = d[i]+t_egress;
						// compute and save the associated path
						//path[s][t].clear();
						P.clear();

						// With the two lines above we reset the public transit routes from the origin to
						// the destination
						this->list_of_stations[s][t].clear(); //s:departing_station; t:arrival_station
						this->set_of_used_lines[s][t].clear();
						this->list_of_lines[s][t].clear();
						u = i;                       
																		// cout<<"\n u_3 \t"<< u<<endl;

						while (V[u].station_id != s){
																		// cout<<"\n line \t"<< __LINE__<<endl;

							P.push_front(u);
							this->list_of_stations[s][t].push_front(V[u].station_id);
																		// cout<<"\n line \t"<< __LINE__<<endl;

							this->list_of_lines[s][t].push_front(V[u].line_id);
																		// cout<<"\n line \t"<< __LINE__<<endl;

							this->set_of_used_lines[s][t].insert(V[u].line_id);
																		// cout<<"\n line \t"<< __LINE__<<endl;
																		// cout<<"\n pred[u] \t"<< pred[u]<<endl;
															cout<<"\n line \t"<< __LINE__<<endl;

														
							if(u<pred.size()){
								u = pred[u];	
																cout<<"\n line \t"<< __LINE__<<endl;

							}else{
																cout<<"\n line \t"<< __LINE__<<endl;

								for(int tmp=0; tmp<pred.size(); tmp++){
									cout<<"\n line\t"<<__LINE__<<endl;
									cout<<"\n pred[tmp]_3 \t"<<pred[tmp]<<endl;
								}
								cout<<"\n pred.size() \t"<<pred.size()<<endl;
								cout<<"\n i \t"<<i<<endl;
								cout<<"\n v.size() \t"<<V.size()<<endl;
								cout<<"\n u_3 \t"<<u<<endl;
								cout<<"warning in line "<<__LINE__<<endl;
								exit(1);
							}
						}
																		// cout<<"\n line \t"<< __LINE__<<endl;

						P.push_front(dep);
																		// cout<<"\n line \t"<< __LINE__<<endl;
																		// cout<<"\n dep \t"<<dep<<endl;
																		// cout<<"\n s \t"<<s<<endl;

						// cout<<"\n t \t"<<t<<endl;
						this->list_of_stations[s][t].push_front(V[dep].station_id);
					}
					path.push_back(P);
				}
			}	

	}
	this->check_consistency();
} // end Dijikstra


// It resets the number of users using each line to 0.
void Global_solution::Update_nbr_usrs(){
	for(int i=0; i<number_of_lines;i++)
			this->L[i].Nbr_users=0;
	this->check_consistency();
}



int Global_solution::station_active_on_any_line(int station_index) const {
	int active_status=0;
	for(int l=0; l<number_of_lines; l++){
		vector<int> active_station_ids = L[l].get_active_station_ids();
		if(std::count(active_station_ids.begin(),active_station_ids.end(),station_index)){
			active_status=1;
			break;
		}
	}
	return active_status;	
}
//To get which Dij matrix I have to use
//customer index not to change the type for everyone and the nearest station
void Global_solution::Customer_category (int customer_index){
	//todo for xiaoyi: no MTTC_C
	
	this->check_consistency();
	int tst=0, i=0;
	//R_W_PT_RS=0, R_RS_PT_W=1, R_W=2,R_PT=3,R_RS=4
	cus.at(customer_index).dstnce=calculateDistance(cus.at(customer_index).orig, cus.at(customer_index).dest)*Circuity_walker;
	cus.at(customer_index).RT=INVALID;

	// cout<<"cus.at(customer_index).FM.size()"<<cus.at(customer_index).FM.size()<<endl;
	// cout<<"cus.at(customer_index).LM.size()"<<cus.at(customer_index).LM.size()<<endl;


	// if d(i, i+n)<= d_{max}, c belongs to R_W
	if(cus.at(customer_index).dstnce < maximum_walk)
		cus.at(customer_index).RT=R_W;
	else{
		while(tst ==0 || i== Number_stations){
			if (station_active_on_any_line(cus.at(customer_index).FM[i].pt.indx)==1){
				//find the closest first mile station
				cus.at(customer_index).DS.dx=cus.at(customer_index).FM[i].pt.indx;
				cus.at(customer_index).fm=cus.at(customer_index).FM[i].dist;
				tst=1;
			}else
				i++;
		}
		
		tst=0,i=0;
		while(tst==0 || i== Number_stations){
			if (station_active_on_any_line(cus.at(customer_index).LM[i].pt.indx) == 1){
				//find the closest last mile station
				cus.at(customer_index).AS.dx=cus.at(customer_index).LM[i].pt.indx;
				cus.at(customer_index).lm=cus.at(customer_index).LM[i].dist;
				tst=1;
			}else
				i++;
		}


		if (cus.at(customer_index).fm+cus.at(customer_index).lm > maximum_walk){
			if(cus.at(customer_index).fm> maximum_walk && cus.at(customer_index).lm > maximum_walk)
				cus.at(customer_index).RT=R_RS;
			else{
				if(cus.at(customer_index).fm > maximum_walk)
					cus.at(customer_index).RT=R_RS_PT_W;
				else
					if(cus.at(customer_index).lm > maximum_walk)
						cus.at(customer_index).RT=R_W_PT_RS;
					else
						if (cus.at(customer_index).fm > cus.at(customer_index).lm )
							cus.at(customer_index).RT = R_RS_PT_W;
						else cus.at(customer_index).RT = R_W_PT_RS;
			}
		}
		else //if first_mile+last_mile < maximum walk, custom belongs to R_PT
			cus.at(customer_index).RT = R_PT;

		if(cus.at(customer_index).DS.dx==cus.at(customer_index).AS.dx)
			cus.at(customer_index).RT = R_RS;
	}
	this->check_consistency();
}//Global_solution::Customer_category


void Global_solution::Define_customer_type(int customer_index, const vector<stop>& set_of_stations){
	this->check_consistency();
	float MTT_c = Maximum_Trip_Coefficient*((this->cus.at(customer_index).dstnce*Circuity)/RSV_speed)*60+20;
	float travel_time = 0;
	// cout<<"\n customer_index \t"<<customer_index<<endl;
	// cout<<"\n MTT_c \t"<<MTT_c<<endl;
    if(this->cus.at(customer_index).DS.dx==this->cus.at(customer_index).AS.dx   && this->cus.at(customer_index).RT !=2){
    	this->cus.at(customer_index).RT = R_RS;
    	this->Category.at(customer_index) = R_RS;
    }
	//todo for xiaoyi: use cus.at(customer_index).RT replace Category.at(customer_index)
    if(this->cus.at(customer_index).RT == R_W_PT_RS){ //if last mile > maximum_walk distance
		// change for xiaoyi: calculate the travel time of ride sharing according to the distance between arrival station and destination
		// cout<<"\n calculateDistance(set_of_stations[cus[customer_index].AS.dx].p,cus[customer_index].dest)/RSV_speed*60 \t"<<calculateDistance(set_of_stations[cus[customer_index].AS.dx].p,cus[customer_index].dest)/RSV_speed*60;
		// cout<<"\n cus[customer_index].AS.dx \t"<<cus[customer_index].AS.dx<<"\t set_of_stations[cus[customer_index].AS.dx].p \t"<<set_of_stations[cus[customer_index].AS.dx].p.Identif<<"\t cus[customer_index].dest \t"<<cus[customer_index].dest.Identif;
    	
		travel_time=this->DIStance.at(this->cus.at(customer_index).DS.dx).at(this->cus.at(customer_index).AS.dx)+
				(this->cus.at(customer_index).DS.wd*60)/walk_speed+
				calculateDistance(set_of_stations[cus[customer_index].AS.dx].p,cus[customer_index].dest)/RSV_speed*60;

    	this->PT_Estimated_EDT_LAT_To_Station.at(customer_index)=floor(// earliest RS departure time at station s 
			this->cus.at(customer_index).EDT+
			this->DIStance.at(this->cus.at(customer_index).DS.dx).at(this->cus.at(customer_index).AS.dx)//average time from DS to AS
			+(this->cus.at(customer_index).DS.wd*60)/walk_speed 
  	  	);

    	if(this->PT_Estimated_EDT_LAT_To_Station.at(customer_index)>1440)
    		this->PT_Estimated_EDT_LAT_To_Station.at(customer_index)=
				this->PT_Estimated_EDT_LAT_To_Station.at(customer_index)-1440;
    	
		this->Category.at(customer_index) = R_W_PT_RS;
		
		// cout<<"\n travel_time_1 \t"<<travel_time<<endl;
    }
    	if(this->cus.at(customer_index).RT == R_RS_PT_W){
			// change for xiaoyi: calculate the travel time of ride sharing according to the distance between origin and departure destination 
    		travel_time=this->DIStance.at(this->cus.at(customer_index).DS.dx).at(this->cus.at(customer_index).AS.dx)+
					   (this->cus.at(customer_index).DS.wd*60)/walk_speed+
					   	calculateDistance(cus[customer_index].orig, set_of_stations[cus[customer_index].DS.dx].p)/RSV_speed*60;

    		this->PT_Estimated_EDT_LAT_To_Station.at(customer_index)=floor( //RS latest arrival time to station s 
																this->cus.at(customer_index).LAT-
																this->DIStance.at(this->cus.at(customer_index).DS.dx).at(this->cus.at(customer_index).AS.dx)-
																(this->cus.at(customer_index).AS.wd*60)/walk_speed 
																);

    		if(this->PT_Estimated_EDT_LAT_To_Station.at(customer_index)<0)
    			this->PT_Estimated_EDT_LAT_To_Station.at(customer_index)=
					1440-this->PT_Estimated_EDT_LAT_To_Station.at(customer_index);
			this->Category.at(customer_index) = R_RS_PT_W;
			
			// cout<<"\n travel_time_2 \t"<<travel_time<<endl;

    	}

    	if(this->cus.at(customer_index).RT == R_PT){
    		travel_time=this->DIStance.at(this->cus.at(customer_index).DS.dx).at(this->cus.at(customer_index).AS.dx);
    		this->PT_Estimated_EDT_LAT_To_Station.at(customer_index) =
				this->cus.at(customer_index).EDT+
				this->DIStance.at(this->cus.at(customer_index).DS.dx).at(this->cus.at(customer_index).AS.dx);
    		this->Category.at(customer_index) = R_PT;
			// cout<<"\n travel_time_3 \t"<<travel_time<<endl;
    	}

    	if(this->cus.at(customer_index).RT == R_RS){
    		this->PT_Estimated_EDT_LAT_To_Station.at(customer_index)=0;
    		this->Category.at(customer_index) = R_RS;
    	}

    	if(this->cus.at(customer_index).RT == R_W){
			this->PT_Estimated_EDT_LAT_To_Station.at(customer_index)=0;
			this->Category.at(customer_index) = R_W;
		}
    	
		if(travel_time > MTT_c   && this->cus.at(customer_index).RT != R_W){
    		this->cus.at(customer_index).RT = R_RS; 
    		this->Category.at(customer_index) = R_RS;
    	}

    	if(this->cus.at(customer_index).RT == R_W_PT_RS ||
    			this->cus.at(customer_index).RT==R_RS_PT_W ||
    			this->cus.at(customer_index).RT == R_PT
    	){
    		set_of_used_lines[cus[customer_index].DS.dx][cus[customer_index].AS.dx];
    		for (set<int>::iterator p = 
    				this->set_of_used_lines[this->cus.at(customer_index).DS.dx][this->cus.at(customer_index).AS.dx].begin(); p !=
    				this->set_of_used_lines[this->cus.at(customer_index).DS.dx][this->cus.at(customer_index).AS.dx].end(); p++
    			)
    		    this->L[*p].Nbr_users++; // We update the number of users that use every line
    	}
	// cout<<"\n >Category.at(customer_index) \t"<<Category.at(customer_index)<<endl;

    this->check_consistency();
 }//Global_solution::Define_customer_type()

//this function compute cust trajectory by the dijkstra algorithm 
void Global_solution::Compute_Customer_PT_trajectory(){
	this->check_consistency();
	int customer_index;
	int ca;
	for(customer_index=0; customer_index<this->cus.size(); customer_index++)
		if(this->cus.at(customer_index).RT!=2 && this->cus.at(customer_index).RT!=4){ //does not including only walking or ride sharing modes
		list<int>::const_iterator p = this->list_of_stations[this->cus.at(customer_index).DS.dx][this->cus.at(customer_index).AS.dx].begin();
		ca=0;
        while(p != this->list_of_stations[this->cus.at(customer_index).DS.dx][this->cus.at(customer_index).AS.dx].end()){
        	this->cus.at(customer_index).traj_stat[0]+= std::to_string(*p); // does not include the departure sation
        	ca++;
        	if(ca != this->list_of_stations[this->cus.at(customer_index).DS.dx][this->cus.at(customer_index).AS.dx].size())
        		this->cus.at(customer_index).traj_stat[0] += ",";
        	p++;
        }
	}
	this->check_consistency();
}

//write ride sharing O-D and earliest departure and latest arrive time, not include PT and Walker					
int Global_solution::write_ride_sharing_results(const vector<stop>& set_of_stations) const{
	string const nomFichier(Consts::result_file);
	ofstream monFlux(nomFichier.c_str());
	int i,flt,aaa,bbb,A;
	Heure HH; // time
	if(monFlux)
	{
		monFlux << "write ride sharing O-D and earliest departure and latest arrive time";
		monFlux <<"\n customer index \t orgin id \t destination id \t earliest depart time \t lastest arrive time \t customer arrive/depart station id \t RS latest/earliest arrive time"<<endl;
	    for ( i = 0; i < this->cus.size(); i++){
	    	if(this->cus.at(i).RT==0) // R_W_PT_RS
	    		monFlux<<"\n"<<i<<"\t"<<this->cus.at(i).orig.Identif<<"\t"<<this->cus.at(i).dest.Identif<< "\t"
				<< this->cus.at(i).EDT<<"\t"<<this->cus.at(i).LAT<< "\t R_W_PT_RS \t"<<set_of_stations.at(this->cus.at(i).AS.dx).p.Identif<<"\t"
				<<this->PT_Estimated_EDT_LAT_To_Station.at(i)<<"\t 0";
	    	if(this->cus.at(i).RT==1) // R_RS_PT_w
				monFlux<<"\n"<<i<<"\t"<<this->cus.at(i).orig.Identif<<"\t"<<this->cus.at(i).dest.Identif<< "\t"
				<< this->cus.at(i).EDT<<"\t"<<this->cus.at(i).LAT<< "\t R_W_PT_RS \t"<<set_of_stations.at(this->cus.at(i).DS.dx).p.Identif<<"\t"
				<<this->PT_Estimated_EDT_LAT_To_Station.at(i)<<"\t 0";
	    	if(this->cus.at(i).RT==4) // R_RS
	    		monFlux<<"\n"<<i<<"\t"<<this->cus.at(i).orig.Identif<<"\t"<<this->cus.at(i).dest.Identif<< "\t"
				<< this->cus.at(i).EDT<<"\t"<<this->cus.at(i).LAT<< "\t R_W_PT_RS \t"<<set_of_stations.at(this->cus.at(i).AS.dx).p.Identif<<"\t"
				<<this->PT_Estimated_EDT_LAT_To_Station.at(i)<<"\t 0";
			}
	}else{
		cout << "ERROR: Unable to open file. line="<<__LINE__ << endl;
				exit(0);
	}
	//cout<<"\n done"<<endl;
	this->check_consistency();
	return 0;

}

//initialize different customer type in 9 zones are 0
void Global_solution::initialisation_avg_TT_and_category_zone(){
	for(int zone_index=0; zone_index<Consts::num_of_zonal_movement_types; zone_index++){
		this->Avg_TT_0[zone_index].RS=0;
		this->Avg_TT_0[zone_index].PT=0;
		for(int b=0;b<Number_categories;b++)
			this->category_zone[zone_index][b]=0;
	}
		this->check_consistency();
}


void Problem::construct_data_problem(const Network &n, const Global_solution& glob_sol, int n_ride_sharing_requests_,
							const vector<stop>& set_of_stations
){
			// We assume that every user has one possible transfer node. Therefore, for each request we
			// have 1 origin, 1 transfer and 1 destination.
			n_ride_sharing_requests = n_ride_sharing_requests_;
			this->n_nodes = 3 * n_ride_sharing_requests + 1; //aa: why "+1"? Please, explain        // Nis:I don't know //aa: I'll ask Yves
			// cout<<"\n ^^^^^^^^Global solution particle "<< glob_sol.nume <<"; n_ride_sharing_requests="<<n_ride_sharing_requests;
			if (Consts::debug_mode)
			{
				if (n_ride_sharing_requests < 0 || n_ride_sharing_requests < 0)
				{
					cout<<"\nLine "<<__LINE__<<": Error: n_ride_sharing_requests="<<n_ride_sharing_requests<<
						" and n_nodes="<<n_nodes<<endl;
					exit(0);
				}
			}

			for (int i=0; i<n_nodes; i++)
			{
				node dummy_node;
				nodes.push_back(dummy_node);
			}

			for (int i=0; i<n_ride_sharing_requests; i++)
			{
				request dummy_req;
				this->requests.push_back(dummy_req);
			}

			//ride_sharing_request_index is the index of the list of customers who use ride_sharing
			// it's different from cus.at(i), as i instead accounts for any customer (ither using ride sharing or not).
			// A customer may have the index 6788 in the list of customer
			// and the index 3456 as ride_sharing_request_index or user of ridesharing
			// (he may also have no ride_sharing_request_index if he doesn't use RS)
			int ride_sharing_request_index = 0;


			const vector<customer>& cus = glob_sol.cus;

			// Read data on the customer requests:
			for(int i=0;i<cus.size();i++)
			{
				const customer& single_customer = cus.at(i);
				//cout<<"\nLine "<<__LINE__<<": reading the "<<i<<"-th customer"<<
				//	"whose type is "<< cus.at(i).RT <<endl;
				if( // single_customer uses a ride sharing
					single_customer.RT == R_W_PT_RS || single_customer.RT == R_RS_PT_W || single_customer.RT == R_RS
				){
					int origin_node = ride_sharing_request_index;
					int destination_node = n_ride_sharing_requests + ride_sharing_request_index;
					int transfer_node = 2 * n_ride_sharing_requests + 1 + ride_sharing_request_index;
					//Note that ride_sharing_request_index is different from i. The objective is to create a list of customers that will use
					//totally or partially ride sharing, so a customer origin_node is different from cus[i].

					//cout<<"\n n_ride_sharing_requests="<<n_ride_sharing_requests<< " counter="<<ride_sharing_request_index;
					//cout<<"\n cus[i].orig.Identif="<<cus.at(i).orig.Identif<<"   cus[i].dest.Identif="<<cus.at(i).dest.Identif<<endl;
					this->nodes.at(origin_node).loc = n.map_key_to_location.at(single_customer.orig.Identif);

					this->nodes.at(destination_node).loc = n.map_key_to_location.at(single_customer.dest.Identif);

					//cout<<"set_of_stations.at(cus.at(i).AS.dx).p.Identif="<<set_of_stations.at(cus.at(i).AS.dx).p.Identif<<endl;

					if(single_customer.RT==R_W_PT_RS)
					{
							this->nodes.at(transfer_node).loc =
							n.map_key_to_location.at(set_of_stations.at(single_customer.AS.dx).p.Identif);           // Note: the discretized location is used here, because the station id number is not provided in the file 'results'
					}
					if(single_customer.RT==R_RS_PT_W)
					{

							this->nodes.at(transfer_node).loc =
							n.map_key_to_location.at(set_of_stations.at(single_customer.DS.dx).p.Identif);           // Note: the discretized location is used here, because the station id number is not provided in the file 'results'
					//cout<<"\n 2nd if"<<endl;
					}

					//cout<<"\n here cus["<<i<<"].RT="<<cus[i].RT<<endl;
							// Request data:
					this->requests.at(ride_sharing_request_index).original_id = single_customer.indxe;
									//cout<<"\n count="<<counter<<endl;

					this->requests.at(ride_sharing_request_index).pickup_nodes.push_back(origin_node);
					this->requests.at(ride_sharing_request_index).delivery_nodes.push_back(destination_node);
							//cout<<"\n in"<<endl;
							// Origin data:
					this->nodes.at(origin_node).lower_tw = single_customer.EDT;


					int matrix_index= this->nodes.at(origin_node).loc * n.n_used_locations +
								 this->nodes.at(destination_node).loc;

					if (Consts::debug_mode)
					{

						if(origin_node >= this->n_nodes || destination_node >= this->n_nodes)
						{
										cout<<"\n origin_node ="<<origin_node <<"   p.n_nodes="<<this->n_nodes
											<<"    destination_node"<<destination_node<<endl;
										cout<<"\n"<<__LINE__<<endl;
										exit(1);
						}
									 // cout<<"\n p.nodes[origin_node].loc="<<p.nodes[origin_node].loc<<"     n.n_used_locations="<<n.n_used_locations<<"    p.nodes[destination_node].loc= "<<p.nodes[destination_node].loc;
						if(this->nodes.at(origin_node).loc < 0 || n.n_used_locations<0 || this->nodes.at(destination_node).loc<0)
						{

											cout<<"\n p.this->nodes[origin_node].loc="<<this->nodes.at(origin_node).loc<<"     n.n_used_locations="<<n.n_used_locations<<"    p.this->nodes[destination_node].loc= "<<this->nodes.at(destination_node).loc;
											exit(1);
						}

									 // cout<<" \n matrix_index="<<matrix_index<<endl;
						if(origin_node >= matrix_index)
						{
							std::stringstream err_msg;
							err_msg<<"\nLine "<<__LINE__<<" WARNING: \n origin_node ="<<origin_node <<"  matrix_index="<<matrix_index
								<<". Is it normal that origin_node is >= matrix_index???????????"<< endl;
							cout<<err_msg.str();
							// aa: Is the condition I described in this WARNING an error or not?
							// Nis: I think it's a warning, I will let Yves confirm
						}
					}//debug_mode



					this->nodes.at(origin_node).upper_tw = single_customer.LAT - (n.time_matrix.at(matrix_index))/2;    // Note: time window tightening
					this->nodes.at(origin_node).load = +1;

					this->nodes.at(origin_node).upper_tw = single_customer.LAT - (n.time_matrix.at(matrix_index))/2;    // Note: time window tightening
					this->nodes.at(origin_node).load = +1;
					//cout<<"\n Arrive ici 2 "<<endl;
					// Destination data:
					this->nodes.at(destination_node).lower_tw =
						single_customer.EDT + (n.time_matrix.at(this->nodes.at(origin_node).loc * n.n_used_locations +
						this->nodes.at(destination_node).loc) )/2;
						this->nodes.at(destination_node).upper_tw = single_customer.LAT;
					//cout<<"\n Arrive ici 3"<<endl;
					this->nodes.at(destination_node).load = -1;
					//cout<<"\n Arrive ici 2"<<endl;
					//cout<<"\n cus.at(i).RT= "<<cus.at(i).RT<<endl;
					if (single_customer.RT == R_RS)
					{

								this->requests.at(ride_sharing_request_index).type = "R_RS";
					}
					// Transfer data (if applicable):
					else
					{
								if (single_customer.RT == R_RS_PT_W)
								{

									this->nodes.at(transfer_node).lower_tw = this->nodes.at(origin_node).lower_tw +
											n.time_matrix.at(this->nodes.at(origin_node).loc * n.n_used_locations +
											this->nodes.at(transfer_node).loc);
									this->nodes.at(transfer_node).upper_tw = glob_sol.PT_Estimated_EDT_LAT_To_Station.at(i);
									this->nodes.at(transfer_node).load = -1;

									this->requests.at(ride_sharing_request_index).delivery_nodes.push_back(transfer_node);
									this->requests.at(ride_sharing_request_index).type = "R_RS_PT_W";
									//cout<<"\n p.nodes["<<transfer_node<<"].lower_tw="<<p.nodes[transfer_node].lower_tw;
									//cout<<"\n p.nodes["<<transfer_node<<"].upper_tw="<<p.nodes[transfer_node].upper_tw;
								}

								else if (single_customer.RT == R_W_PT_RS)
								{
									//cout<<" ligne="<<__LINE__<<endl;
									this->nodes.at(transfer_node).lower_tw = glob_sol.PT_Estimated_EDT_LAT_To_Station.at(i);
									//cout<<" ligne="<<__LINE__<<endl;
									//cout<<"\n transfer_node="<<transfer_node<<endl;
									//cout<<"\n p.nodes["<<destination_node<<"].upper_tw="<<nodes.at(destination_node).upper_tw<<endl;

									//cout<<"\n p.nodes["<<destination_node<<"].loc="<<nodes.at(destination_node).loc<<endl;
									//cout<<"\n nodes.at(transfer_node).loc="<<nodes.at(transfer_node).loc<<endl;
									//cout<<"\n Line "<<__LINE__<<": i="<<i<<"-th customer"<<endl;
									//cout<<"\n n.time_matrix.at(nodes.at(transfer_node).loc * n.n_used_locations +nodes.at(destination_node).loc)="<<nodes.at(transfer_node).loc * n.n_used_locations +nodes.at(destination_node).loc<<endl;
									this->nodes.at(transfer_node).upper_tw = this->nodes.at(destination_node).upper_tw -n.time_matrix.at(this->nodes.at(transfer_node).loc * n.n_used_locations +this->nodes.at(destination_node).loc);
									//cout<<" ligne="<<__LINE__<<endl;
									this->nodes.at(transfer_node).load = +1;
									//cout<<" ligne="<<__LINE__<<endl;
									this->requests.at(ride_sharing_request_index).pickup_nodes.push_back(transfer_node);
									//cout<<" ligne="<<__LINE__<<endl;
									this->requests.at(ride_sharing_request_index).type = "R_W_PT_RS";
									//cout<<" ligne="<<__LINE__<<endl;


								}
					}//else, i.e., if customer != 4

					if (Consts::debug_mode && false)
					{
						cout<<"\nLine "<<__LINE__<<": The i="<<i<<"-th customer"<<
							" has been read. Her type is "<< single_customer.RT <<
							". Her ride_sharing_request_index="<<ride_sharing_request_index<<
							". We will now check if everything is correct"<<endl;
						this->requests.at(ride_sharing_request_index).check_consistency();
						cout<<"\nLine "<<__LINE__<<": The i="<<i<<"-th customer"<<
							" has been correctly read. Her type is "<< single_customer.RT <<endl;
					}

					ride_sharing_request_index += 1;
				}//if the user uses ride_sharing
			
			}//for-loop on customers
			//cout<<"\n Out "<<endl;


			// Complete the data for the artificial depot node:
			int depot_node_id = 2 * n_ride_sharing_requests;
			node depot_node;
			depot_node.loc = n.n_cust_locations;
			depot_node.lower_tw = 0.0;
			depot_node.upper_tw = 1440.0;
			depot_node.load = 0;
			this->nodes.at(depot_node_id)=depot_node;
}//Problem::construct_data_problem



void Problem::compute_user_dissimilarity(const Network &n, int n_ride_sharing_requests){

			for (int i=0; i<n_ride_sharing_requests * n_ride_sharing_requests; i++)
				dissimilarity.push_back(Consts::max_dissimilarity);

			vector<double> time_dissimilarity(n_ride_sharing_requests * n_ride_sharing_requests);
			vector<double> space_dissimilarity(n_ride_sharing_requests * n_ride_sharing_requests);

			// Compute the difference between the earliest pickup times and between the latest delivery times for all pairs of requests:
			double largest_time_dissimilarity = 0.0;
			for (int request1 = 0; request1 < n_ride_sharing_requests; request1++)
			{
				int pickup1 = request1;
				int delivery1 = request1 + n_ride_sharing_requests;
				for (int request2 = request1 + 1; request2 < n_ride_sharing_requests; request2++) {
					int pickup2 = request2;
					int delivery2 = request2 + n_ride_sharing_requests;
					//cout<<"\n"<<__LINE__<<": nodes.at(pickup2).lower_tw="<<nodes.at(pickup2).lower_tw<<endl;
					//cout<<"\n"<<__LINE__<<": nodes.at(pickup1).lower_tw="<<nodes.at(pickup1).lower_tw<<endl;
					double t = 0.5 * (abs(nodes.at(pickup2).lower_tw - nodes.at(pickup1).lower_tw) + abs(nodes.at(pickup2).upper_tw - nodes.at(pickup1).upper_tw));
					time_dissimilarity.at(request1 * n_ride_sharing_requests + request2) = t;
					time_dissimilarity.at(request2 * n_ride_sharing_requests + request1) = t;
					if (t > largest_time_dissimilarity)
					{
						largest_time_dissimilarity = t;
					}
				}
			}

			// Compute the travel time between both origins and between both destinations for all pairs of requests:
			double largest_distance = 0.0;
			for (int request1 = 0; request1 < n_ride_sharing_requests; request1++)
			{
				int pickup1 = request1;
				int delivery1 = request1 + n_ride_sharing_requests;
				for (int request2 = request1 + 1; request2 < n_ride_sharing_requests; request2++)
				{
					int pickup2 = request2;
					int delivery2 = request2 + n_ride_sharing_requests;
					double d = 0.5 * (n.time_matrix.at(nodes.at(pickup1).loc * n.n_used_locations + nodes.at(pickup2).loc) + n.time_matrix.at(nodes.at(delivery1).loc * n.n_used_locations + nodes.at(delivery2).loc) );
					space_dissimilarity.at(request1 * n_ride_sharing_requests + request2) = d;
					space_dissimilarity.at(request2 * n_ride_sharing_requests + request1) = d;
					if (d > largest_distance)
					{
						largest_distance = d;
					}
				}
			}

			// Convert the dissimilarity of two requests to a fractional number between 0 and 1:
			for (int request1 = 0; request1 < n_ride_sharing_requests; request1++)
			{
				dissimilarity.at(request1 * n_ride_sharing_requests + request1) = 0.0;
				for (int request2 = request1 + 1; request2 < n_ride_sharing_requests; request2++)
				{
					double dissimilarity_score = 0.5 * (time_dissimilarity.at(request1 * n_ride_sharing_requests + request2) / largest_time_dissimilarity + space_dissimilarity.at(request1 * n_ride_sharing_requests + request2) / largest_distance);
					dissimilarity.at(request1 * n_ride_sharing_requests + request2) = dissimilarity_score;
					dissimilarity.at(request2 * n_ride_sharing_requests + request1) = dissimilarity_score;
				}
			}

}//Problem::compute_user_dissimilarity



void Solution::initialize_solution(const Network &n, const Problem &p, const int n_ride_sharing_requests)
{
	max_vehicles = n_ride_sharing_requests; //aa: I added this line (in the worst case we have one car per request

	for (int v=0; v<max_vehicles; v++)
	{
		route r;

		// aa: We initialize the route of vehicles by starting and terminating at the depot      
		r.route_structure.push_back(2 * n_ride_sharing_requests);
		r.route_structure.push_back(2 * n_ride_sharing_requests);
		routes.push_back(r);
	}

	obj_value = 0.0;
	vehicles_used = 0;


	//cout<<"\n "<<__LINE__<<": After for  p.n_nodes="<<p.n_nodes<<"  p.n_ride_sharing_requests="<<n_ride_sharing_requests<<endl;

	for (int i=0; i<n_ride_sharing_requests; i++) selected_pickup.push_back(-1);
	for (int i=0; i<n_ride_sharing_requests; i++) selected_delivery.push_back(-1);
	for (int i=0; i<n_ride_sharing_requests; i++) vehicle.push_back(-1);
	for (int i=0; i<p.n_nodes; i++) load.push_back(0);
	for (int i=0; i<p.n_nodes; i++) schedule.push_back(0.0);

	request_bank.clear();
	for (int r = 0; r < n_ride_sharing_requests; r++) {request_bank.push_back(r);}
	//cout<<"\n request_bank.push_back(r)="<<s.request_bank.size();

}//Solution::initialize

bool check_load_and_schedule(const Network &n, const Problem &p, Solution &s, int v, int start_pos) {

	int current_node = s.routes.at(v).route_structure.at(start_pos - 1);

	int load = s.load.at(current_node);
	double time = s.schedule.at(current_node);

	for (int i = start_pos; i < s.routes.at(v).route_structure.size() - 1; i++) {

		int preceding_node = current_node;
		current_node = s.routes.at(v).route_structure.at(i);

		load += p.nodes.at(current_node).load;
		if (load > Consts::vehicle_cap) {
			return false;
		}

		time = max(
			p.nodes.at(current_node).lower_tw,
			time
			+ n.time_matrix.at(p.nodes.at(preceding_node).loc * n.n_used_locations + p.nodes.at(current_node).loc)
		);

		if (time > p.nodes.at(current_node).upper_tw) {
			return false;
		}
	}

	return true;
}


void update_load_and_schedule(const Network &n, const Problem &p, Solution &s, int v, int start_pos) {

	int current_node = s.routes.at(v).route_structure.at(start_pos - 1);

	int load = s.load.at(current_node);
	double time = s.schedule.at(current_node);

	for (int i = start_pos; i < s.routes.at(v).route_structure.size() - 1; i++) {

		int preceding_node = current_node;
		current_node = s.routes.at(v).route_structure.at(i);

		load += p.nodes.at(current_node).load;
		s.load.at(current_node) = load;

		time = max(
			p.nodes.at(current_node).lower_tw,
			time + n.time_matrix.at(
				p.nodes.at(preceding_node).loc * n.n_used_locations + p.nodes.at(current_node).loc
				)
			);
		s.schedule.at(current_node) = time;
	}
}


void perform_best_insertion(const Network &n, const Problem &p, Solution &s, int r)
{
	if (Consts::debug_mode)
	{
		if (r >= p.nodes.size() || r+p.n_ride_sharing_requests >= p.nodes.size())
		{
			cout<<"\n Line "<<__LINE__<<": Error: in the next line we are trying to access p.nodes.at(r) and "
				<<" p.nodes.at(r+p.n_ride_sharing_requests), but p.nodes.size()="<<p.nodes.size()<<", while r="<< r
				<<" and r+p.n_ride_sharing_requests="<<r+p.n_ride_sharing_requests<<endl;
			exit(0);
		}
		p.requests.at(r).check_consistency();
	}

	int idx = p.nodes.at(r).loc * n.n_used_locations + p.nodes.at(r+p.n_ride_sharing_requests).loc;
	//what is the meaning of this?  // Nis: I don't know //aa: I'll ask Yves
	
	assert(idx < n.time_matrix.size() );

	// Assume that the cost of inserting in an empty vehicle (with an artificial depot) is an upper bound on the insertion cost:
	double upper_bound = Consts::vehicle_cost + n.time_matrix.at(idx) + 0.001;

	double best_insertion_cost = upper_bound;

	int best_vehicle = -1;
	int best_sel_pickup = -1;
	int best_sel_delivery = -1;
	int best_pos_pickup = -1;
	int best_pos_delivery = -1;

	// Loop over all vehicles:
	for (int v = 0; v < s.max_vehicles; v++)
	{

		// Loop over all possible pickup nodes:
		for (int i = 0; i < p.requests.at(r).pickup_nodes.size(); i++)
		{
			int sel_pickup = p.requests.at(r).pickup_nodes.at(i);

			// Loop over all possible delivery nodes:
			for (int j = 0; j < p.requests.at(r).delivery_nodes.size(); j++) {
				int sel_delivery = p.requests.at(r).delivery_nodes.at(j);

				// Consider all possible positions for the selected pickup node:
				for (int pos_pickup = 1; pos_pickup < s.routes.at(v).route_structure.size(); pos_pickup++) {

					int predecessor_pickup = s.routes.at(v).route_structure.at(pos_pickup - 1);
					int successor_pickup = s.routes.at(v).route_structure.at(pos_pickup);

					// Break based on time windows:
					if (s.schedule.at(predecessor_pickup)
						+ n.time_matrix.at(p.nodes.at(predecessor_pickup).loc * n.n_used_locations + p.nodes.at(sel_pickup).loc)
							> p.nodes.at(sel_pickup).upper_tw)
					{
						break;
					}

					// Preliminary load and time window check:
					if (
						s.load.at(predecessor_pickup) + p.nodes.at(sel_pickup).load <= Consts::vehicle_cap &&
						p.nodes.at(sel_pickup).lower_tw +
						n.time_matrix.at(
							p.nodes.at(sel_pickup).loc * n.n_used_locations + p.nodes.at(successor_pickup).loc
							)
							<= p.nodes.at(successor_pickup).upper_tw
					){

						// Compute the cost of inserting the selected pickup node:
						double insertion_cost_pickup = 0.0
							+ n.time_matrix.at(p.nodes.at(predecessor_pickup).loc * n.n_used_locations + p.nodes.at(sel_pickup).loc)
							+ n.time_matrix.at(p.nodes.at(sel_pickup).loc * n.n_used_locations + p.nodes.at(successor_pickup).loc)
							- n.time_matrix.at(p.nodes.at(predecessor_pickup).loc * n.n_used_locations + p.nodes.at(successor_pickup).loc);

						// Take into account the fixed vehicle activation cost if needed:
						bool empty_vehicle = false;
						if (s.routes.at(v).route_structure.size() == 2) {
							empty_vehicle = true;
							insertion_cost_pickup += Consts::vehicle_cost;
						}

						if (insertion_cost_pickup < best_insertion_cost) {

							// Insert the selected pickup node:
							s.routes.at(v).route_structure.insert(s.routes.at(v).route_structure.begin() + pos_pickup, sel_pickup);

							// Consider all possible positions for the delivery node:
							for (int pos_delivery = pos_pickup + 1; pos_delivery < s.routes.at(v).route_structure.size(); pos_delivery++)
							{
								int predecessor_delivery = s.routes.at(v).route_structure.at(pos_delivery - 1);
								int successor_delivery = s.routes.at(v).route_structure.at(pos_delivery);

								// Break based on time windows:
								if (p.nodes.at(predecessor_delivery).lower_tw +
									n.time_matrix.at(p.nodes.at(predecessor_delivery).loc * n.n_used_locations + p.nodes.at(sel_delivery).loc)
									> p.nodes.at(sel_delivery).upper_tw)
								{
									break;
								}

								// Preliminary time window check:
								if (
									p.nodes.at(sel_delivery).lower_tw
									+ n.time_matrix.at(p.nodes.at(sel_delivery).loc * n.n_used_locations + p.nodes.at(successor_delivery).loc)
									<= p.nodes.at(successor_delivery).upper_tw)
								{

									// Compute the additional cost of inserting the selected delivery node:
									double insertion_cost = insertion_cost_pickup
										+ n.time_matrix.at(p.nodes.at(predecessor_delivery).loc * n.n_used_locations + p.nodes.at(sel_delivery).loc)
										+ n.time_matrix.at(p.nodes.at(sel_delivery).loc * n.n_used_locations + p.nodes.at(successor_delivery).loc)
										- n.time_matrix.at(p.nodes.at(predecessor_delivery).loc * n.n_used_locations + p.nodes.at(successor_delivery).loc);

									if (insertion_cost < best_insertion_cost)
									{

										// Insert the selected delivery node:
										s.routes.at(v).route_structure.insert(s.routes.at(v).route_structure.begin() + pos_delivery, sel_delivery);

										if (empty_vehicle == true || check_load_and_schedule(n, p, s, v, pos_pickup) == true) {

											best_insertion_cost = insertion_cost;
											best_vehicle = v;
											best_sel_pickup = sel_pickup;
											best_sel_delivery = sel_delivery;
											best_pos_pickup = pos_pickup;
											best_pos_delivery = pos_delivery;
										}

										// Undo the insertion of the delivery node:
										s.routes.at(v).route_structure.erase(s.routes.at(v).route_structure.begin() + pos_delivery);
									}
								}
							}

							// Undo the insertion of the pickup node:
							s.routes.at(v).route_structure.erase(s.routes.at(v).route_structure.begin() + pos_pickup);
						}
					}
				}
			}
		}
	}

	// Insert the request in the optimal manner and remove it from the request bank:
	if (best_insertion_cost < upper_bound) {

		s.obj_value += best_insertion_cost;
		if (s.routes.at(best_vehicle).route_structure.size() == 2) {
			s.vehicles_used += 1;
		}

		s.vehicle.at(r) = best_vehicle;
		s.selected_pickup.at(r) = best_sel_pickup;
		s.selected_delivery.at(r) = best_sel_delivery;

		s.routes.at(best_vehicle).route_structure.insert(s.routes.at(best_vehicle).route_structure.begin() + best_pos_pickup, best_sel_pickup);
		s.routes.at(best_vehicle).route_structure.insert(s.routes.at(best_vehicle).route_structure.begin() + best_pos_delivery, best_sel_delivery);

		s.request_bank.erase(s.request_bank.begin() + distance(s.request_bank.begin(), find(s.request_bank.begin(), s.request_bank.end(), r)));

		// Update the load and schedule of the adapted route:
		update_load_and_schedule(n, p, s, best_vehicle, best_pos_pickup);
	}

}//perform_best_insertion


void random_order_best_insertion(const Network &n, const Problem &p, Solution &s)
{

	// Shuffle the request bank:
	//cout<<"\n line="<<__LINE__;
	random_shuffle(s.request_bank.begin(), s.request_bank.end());
  //  cout<<"\n line="<<__LINE__;
	// Select the next request to be inserted:
	while (s.request_bank.size() != 0){
		int r = s.request_bank.at(0);
		perform_best_insertion(n, p, s, r);
		// Break if the request could not be inserted in a feasible manner:
		if (s.vehicle.at(r) == -1) {
			cout<<" \n line"<<__LINE__<<"We did not manage to serve this request r="<<r<<"   p.request="<<p.requests[r].original_id<<"  type="<<p.requests.at(r).type;
			cout<<"\n line="<<__LINE__;
			break;
		}
	}
}

void Solution::copy_solution_from(const Network &n, const Problem &p, const Solution &s)
{

	obj_value = s.obj_value;
	vehicles_used = s.vehicles_used;
	max_vehicles = s.max_vehicles;

	assert(routes.size()==s.routes.size());
	routes.clear();
	for (int v = 0; v < max_vehicles; v++)
	{
		route rt;
		//routes.at(v).route_structure.clear();
		for (int i = 0; i < s.routes.at(v).route_structure.size(); i++)
		{
			rt.route_structure.push_back(s.routes.at(v).route_structure.at(i));
		};
		routes.push_back(rt); // We are filling routes.at(v)
	}

	for (int n = 0; n < p.n_nodes; n++)
	{
		load.at(n) = s.load.at(n);
		schedule.at(n) = s.schedule.at(n);
	}

	for (int r = 0; r < p.n_ride_sharing_requests; r++) 
	{
		vehicle.at(r) = s.vehicle.at(r);
		selected_pickup.at(r) = s.selected_pickup.at(r);
		selected_delivery.at(r) = s.selected_delivery.at(r);
	}

	request_bank.clear();
	for (int i = 0; i < s.request_bank.size(); i++)
	{
		request_bank.push_back(s.request_bank.at(i));
	};
}


void remove_request(const Network &n, const Problem &p, Solution &s, int r)
{

	int v = s.vehicle.at(r);

	int pickup = s.selected_pickup.at(r);
	int delivery = s.selected_delivery.at(r);
	int pos_pickup = distance(s.routes.at(v).route_structure.begin(), find(s.routes.at(v).route_structure.begin(), s.routes.at(v).route_structure.end(), pickup));
	int pos_delivery = distance(s.routes.at(v).route_structure.begin(), find(s.routes.at(v).route_structure.begin(), s.routes.at(v).route_structure.end(), delivery));

	int predecessor_delivery = s.routes.at(v).route_structure.at(pos_delivery - 1);
	int successor_delivery = s.routes.at(v).route_structure.at(pos_delivery + 1);

	s.obj_value = s.obj_value
		- n.time_matrix.at(p.nodes.at(predecessor_delivery).loc * n.n_used_locations + p.nodes.at(delivery).loc)
		- n.time_matrix.at(p.nodes.at(delivery).loc * n.n_used_locations + p.nodes.at(successor_delivery).loc)
		+ n.time_matrix.at(p.nodes.at(predecessor_delivery).loc * n.n_used_locations + p.nodes.at(successor_delivery).loc);

	s.routes.at(v).route_structure.erase(s.routes.at(v).route_structure.begin() + pos_delivery);

	int predecessor_pickup = s.routes.at(v).route_structure.at(pos_pickup - 1);
	int successor_pickup = s.routes.at(v).route_structure.at(pos_pickup + 1);

	s.obj_value = s.obj_value
		- n.time_matrix.at(p.nodes.at(predecessor_pickup).loc * n.n_used_locations + p.nodes.at(pickup).loc)
		- n.time_matrix.at(p.nodes.at(pickup).loc * n.n_used_locations + p.nodes.at(successor_pickup).loc)
		+ n.time_matrix.at(p.nodes.at(predecessor_pickup).loc * n.n_used_locations + p.nodes.at(successor_pickup).loc);

	s.routes.at(v).route_structure.erase(s.routes.at(v).route_structure.begin() + pos_pickup);

	s.vehicle.at(r) = -1;
	s.selected_pickup.at(r) = -1;
	s.selected_delivery.at(r) = -1;

	s.request_bank.push_back(r);

	if (s.routes.at(v).route_structure.size() == 2) {
		s.vehicles_used -= 1;
		s.obj_value -= Consts::vehicle_cost;
	}
}

void random_request_removal(const Network &n, const Problem &p, Solution &s) {

	double removal_pct = (Consts::max_removal_pct - 0.05) * ((double)rand() / (double)RAND_MAX) + 0.05;

	while (s.request_bank.size() < p.n_ride_sharing_requests * removal_pct) {
		// cout<<"p.n_ride_sharing_requests "<<p.n_ride_sharing_requests<<endl;	

		// Remove a random request:
		int request = rand() % p.n_ride_sharing_requests;
		if (s.vehicle.at(request) != -1) { remove_request(n, p, s, request); }
	}

	// Update the load and schedule of each route:
	for (int v = 0; v < s.max_vehicles; v++) {
		update_load_and_schedule(n, p, s, v);
	}
}


void worst_request_removal(const Network &n, const Problem &p, Solution &s) {

	double removal_pct = (Consts::max_removal_pct - 0.05) * ((double)rand() / (double)RAND_MAX) + 0.05;

	double *removal_benefit = new double[p.n_ride_sharing_requests];
	fill(removal_benefit, removal_benefit + p.n_ride_sharing_requests, 0.0);
	// cout<<"p.n_ride_sharing_requests "<<p.n_ride_sharing_requests<<endl;	

	// Compute the benefit of removing every request:
	for (int r = 0; r < p.n_ride_sharing_requests; r++) {

		int v = s.vehicle.at(r);

		int pickup = s.selected_pickup.at(r);
		int delivery = s.selected_delivery.at(r);
		int pos_pickup = distance(s.routes.at(v).route_structure.begin(), find(s.routes.at(v).route_structure.begin(), s.routes.at(v).route_structure.end(), pickup));
		int pos_delivery = distance(s.routes.at(v).route_structure.begin(), find(s.routes.at(v).route_structure.begin(), s.routes.at(v).route_structure.end(), delivery));

		int predecessor_pickup = s.routes.at(v).route_structure.at(pos_pickup - 1);
		int successor_pickup = s.routes.at(v).route_structure.at(pos_pickup + 1);
		int predecessor_delivery = s.routes.at(v).route_structure.at(pos_delivery - 1);
		int successor_delivery = s.routes.at(v).route_structure.at(pos_delivery + 1);

		if (s.routes.at(v).route_structure.size() == 4) {
			removal_benefit[r] = Consts::vehicle_cost + n.time_matrix.at(p.nodes.at(pickup).loc * n.n_used_locations
						+ p.nodes.at(delivery).loc);
		}

		else if (pos_delivery == pos_pickup + 1) {
			removal_benefit[r] =
				n.time_matrix.at(p.nodes.at(predecessor_pickup).loc * n.n_used_locations + p.nodes.at(pickup).loc)
				+ n.time_matrix.at(p.nodes.at(pickup).loc * n.n_used_locations + p.nodes.at(delivery).loc)
				+ n.time_matrix.at(p.nodes.at(delivery).loc * n.n_used_locations + p.nodes.at(successor_delivery).loc)
				- n.time_matrix.at(p.nodes.at(predecessor_pickup).loc * n.n_used_locations + p.nodes.at(successor_delivery).loc);
		}

		else {
			removal_benefit[r] =
				n.time_matrix.at(p.nodes.at(predecessor_pickup).loc * n.n_used_locations + p.nodes.at(pickup).loc)
				+ n.time_matrix.at(p.nodes.at(pickup).loc * n.n_used_locations + p.nodes.at(successor_pickup).loc)
				- n.time_matrix.at(p.nodes.at(predecessor_pickup).loc * n.n_used_locations + p.nodes.at(successor_pickup).loc);
				+ n.time_matrix.at(p.nodes.at(predecessor_delivery).loc * n.n_used_locations + p.nodes.at(delivery).loc)
				+ n.time_matrix.at(p.nodes.at(delivery).loc * n.n_used_locations + p.nodes.at(successor_delivery).loc)
				- n.time_matrix.at(p.nodes.at(predecessor_delivery).loc * n.n_used_locations + p.nodes.at(successor_delivery).loc);
		}
	}

	// Iteratively remove the request for which the benefit is largest:
	while (s.request_bank.size() < p.n_ride_sharing_requests * removal_pct) {
		int selected_request = -1;
		double largest_benefit = 0.0;
		for (int r = 0; r < p.n_ride_sharing_requests; r++) {
			if (removal_benefit[r] > largest_benefit) {
				largest_benefit = removal_benefit[r];
				selected_request = r;
			}
		}
		remove_request(n, p, s, selected_request);
		removal_benefit[selected_request] = 0.0;
	}

	// Update the load and schedule of each route:
	for (int v = 0; v < s.max_vehicles; v++) {
		update_load_and_schedule(n, p, s, v);
	}

	delete[] removal_benefit;
}



void related_request_removal(const Network &n, const Problem &p, Solution &s) {

	double removal_pct = (Consts::max_removal_pct - 0.05) * ((double)rand() / (double)RAND_MAX) + 0.05;


	double *dissimilarity = new double[p.n_ride_sharing_requests];
	fill(dissimilarity, dissimilarity + p.n_ride_sharing_requests, 10000.0);
	// cout<<"TEST "<<rand() % p.n_ride_sharing_requests<<endl;	

	// Remove one random request:
	int r1 = rand() % p.n_ride_sharing_requests;
	// cout<<"DEBUG "<<__LINE__<<endl;	

	remove_request(n, p, s, r1);
	dissimilarity[r1] = 10000.0;
	// cout<<"DEBUG "<<__LINE__<<endl;	
	// Look up the dissimilarity to every other request:
	for (int r2 = 0; r2 < p.n_ride_sharing_requests; r2++) {
		if (r2 != r1) {
			dissimilarity[r2] = p.dissimilarity[r1 * p.n_ride_sharing_requests + r2];
		}
	}

	// Iteratively remove the most related requests:
	while (s.request_bank.size() < p.n_ride_sharing_requests * removal_pct) {
		int selected_request = -1;
		double smallest_dissimilarity = 9999.9;
		for (int request = 0; request < p.n_ride_sharing_requests; request++) {
			if (dissimilarity[request] < smallest_dissimilarity) {
				smallest_dissimilarity = dissimilarity[request];
				selected_request = request;
			}
		}
		remove_request(n, p, s, selected_request);
		dissimilarity[selected_request] = 10000.0;
		// cout<<"DEBUG "<<__LINE__<<endl;	

	}

	// Update the load and schedule of each route:
	for (int v = 0; v < s.max_vehicles; v++) {
		update_load_and_schedule(n, p, s, v);
		// cout<<"DEBUG "<<__LINE__<<endl;	

	}

	delete[] dissimilarity;
	
}


double Global_solution::compute_travel_times_RS( const Problem &p, const Solution &s)
{
	this->check_consistency();
	double sum_travel_times = 0.0;

	for (int r = 0; r < p.n_ride_sharing_requests; r++)
	{
		int selected_pickup = s.selected_pickup.at(r);
		int selected_delivery = s.selected_delivery.at(r);

		if(Category.at( p.requests.at(r).original_id )==0 ||
			Category.at(p.requests.at(r).original_id)==1 ||
			Category.at(p.requests.at(r).original_id)==4)
		{
			cus[p.requests.at(r).original_id].TTT.RS = (s.schedule.at(selected_delivery) - s.schedule.at(selected_pickup));
			sum_travel_times += cus[p.requests.at(r).original_id].TTT.RS;
		}
		else
			cus[p.requests.at(r).original_id].TTT.RS =0;
	}
	this->check_consistency();
	return sum_travel_times;
}

void random_route_removal(const Network &n, const Problem &p, Solution &s) {

	double removal_pct = (Consts::max_removal_pct - 0.05) * ((double)rand() / (double)RAND_MAX) + 0.05;
	// cout<<"p.n_ride_sharing_requests "<<p.n_ride_sharing_requests<<endl;	

	while (s.request_bank.size() < p.n_ride_sharing_requests * removal_pct) {

		// Select a random route:
		int route = rand() % s.max_vehicles;

		// Iteratively remove the first request in the route until the route is empty:
		while (s.routes[route].route_structure.size() != 2) {

			int first_request = s.routes.at(route).route_structure.at(1);
			if (first_request > 2 * p.n_ride_sharing_requests) {
				first_request = first_request - 2 * p.n_ride_sharing_requests - 1;     // Note: this assumes that every user has a single potential transfer station (to be adapted afterwards)
			}
			remove_request(n, p, s, first_request);
		}

		// Update the load and schedule of the route:
		update_load_and_schedule(n, p, s, route);
	}
}


// Run the LNS algorithm:
// We use const Network& to ensure that we do not modify the network
int Global_solution::run_lns_algorithm(int particle_index, const Network &n,
	const vector<stop> set_of_stations){
	this->check_consistency();

	Problem p;
	int n_req = this->RideSharing_Users;
	p.construct_data_problem(n, *this, n_req, set_of_stations);
	p.compute_user_dissimilarity(n, n_req);

	// Create and initialize solutions:
	Solution s_curr;
	Solution s_best;
	Solution s_temp;
	//cout<<"\n Arrive là -1    n_req="<<n_req<<endl;
	s_curr.initialize_solution(n, p, n_req);
	//cout<<"\n ici A"<<endl;
	s_best.initialize_solution(n, p, n_req);
	//cout<<"\n ici B"<<endl;
	s_temp.initialize_solution(n, p, n_req);
	//cout<<"\n Apres Init"<<endl;
	// Construct an initial solution:
	//cout<<"\n Arrive là 0"<<endl;
	random_order_best_insertion(n, p, s_curr);
	// cout<<"\n Après random"<<endl;
	// cout<<"\n Arrive là "<<endl;
	int ii=0;
	int debug_counter = 0;
	if (Consts::debug_mode)
	{	
		// cout<<"debug_counter:	"<<debug_counter<<endl;
		if (s_curr.request_bank.size() != 0)
		{
		cout<< "\n error in line "<<__LINE__<<": request bank should be empty, since at the beginning we should always be able to satisfy all ride sharing users, because in the worst case a dedicated car will be given to each user. The reasons for this error could be: (i) one of the users has an infeasible time window constraint or (ii) teh numbe of available cars is less than the number of users"<<endl;
		exit(0);
		}
	}
	/*	while (s_curr.request_bank.size() != 0) {
			debug_counter++;
			if (debug_counter > max_vehicles)
				cout<<"\n line"<<__LINE__<<": This while loop has done "<<debug_counter<<iterations;
			cout<<"\n 1. s_curr.request_bank.size()"<<s_curr.request_bank.size()<<endl;
			random_request_removal(n, p, s_curr);
			//cout<<"\n Milieu while"<<endl;
			//cout<<"\n 2. s_curr.request_bank.size()"<<s_curr.request_bank.size()<<endl;
			random_order_best_insertion(n, p, s_curr);
			ii++;
		}*/

	// Update the best and overall best solution with the initial solution:
	s_best.copy_solution_from(n, p, s_curr);

	printf("\n line %d: Global solution particle %d New best LNS solution found in initialisation with %d vehicles and total travel time %2f \n",
			__LINE__, particle_index, s_best.vehicles_used, s_best.obj_value);
	// The fleet size used in the initial solution is the maximum fleet size during the rest of the procedure:
	//max_vehicles = s_curr.vehicles_used;

	// Perform a prefixed number of iterations of the algorithm:
	for (int i = 0; i < Consts::LNS_iterations; i++)
	{	
		// Copy the current solution into s_temp:
		s_temp.copy_solution_from(n, p, s_curr);
		// cout<<"debug	"<<__LINE__<<endl;

		// Perform an LNS removal operator on the current solution:
		int r = rand() % 4;
		switch (r) {
			case 0: 
			// cout<<"case 0 "<<__LINE__<<endl; 
				random_request_removal(n, p, s_curr);break;
			case 1:
				// cout<<"case 1 "<<__LINE__<<endl; 
				related_request_removal(n, p, s_curr);break;
			case 2: 
				// cout<<"case 2 "<<__LINE__<<endl; 
				worst_request_removal(n, p, s_curr);break;
			case 3:
				// cout<<"case 3 "<<__LINE__<<endl; 
				random_route_removal(n, p, s_curr); break;
		}
		// cout<<"debug	"<<__LINE__<<endl;

		// Perform an LNS insertion operator on the current solution:
		random_order_best_insertion(n, p, s_curr);

		// cout<<"debug	"<<__LINE__<<endl;

		// Set the solution back to s_temp if the new solution is more than x percent worse than the previous solution:
		if (s_curr.request_bank.size() != 0 || s_curr.obj_value > Consts::max_deterioration_factor * s_best.obj_value) {
			s_curr.copy_solution_from(n, p, s_temp);
		}
		// cout<<"debug	"<<__LINE__<<endl;

		// Otherwise, update the best solution if needed:
	//	if (s_curr.request_bank.size() == 0 && s_curr.obj_value < s_best.obj_value - 0.001) {
	//		s_best.copy_solution_from(n, p, s_curr);
		//	printf("\nLine %d: Global solution particle %d: New best solution found in LNS iteration %d with %d vehicles and total travel time %2f \n",
			//	__LINE__, this->nume, it + 1, s_best.vehicles_used, s_best.obj_value);
		//}
	}

	this->CTT_RS = this->compute_travel_times_RS(p, s_best);

	//Pop.at(0).ATT=0;//aa: I needed to move this line outside of this function. I moved it right before
										// every call to the Eval(..) function. Is it ok what I did? What is the purpose of this
										// line?
	                 //Nis: I try to remember ! but i couldn't .. the strange thing is I didn't find that line in previous versions !!
	
	this->check_consistency();
	// cout<<"debug	"<<__LINE__<<endl;

	return s_best.vehicles_used;
}//Global_solution::run_lns_algorithm



double Global_solution::Eval(const Network &n, const vector<stop>& set_of_stations, int particle_index){	
	int customer_type_zone[9]={0,0,0,0,0,0,0,0,0};
	for(int l=0; l<number_of_lines; l++)
		this->obj.fit_1+=this->L[l].NV;
		this->obj.fit_2 =
			this->run_lns_algorithm(particle_index,n,set_of_stations); //calculate ride sharing vehicle number 
	std::stringstream solution_identifier_string_;
	solution_identifier_string_<<"Global solution particle "<<this->nume; // nume: index of the solution (0 for the base solution and PT+1 for the population solutions
	string solution_identifier_string = solution_identifier_string_.str();
	this->initialisation_avg_TT_and_category_zone(); //initialize the category_zone[zone_index][b]=0

	for(int customer_index=0; customer_index<cus.size(); customer_index++){
		if(cus[customer_index].orig.location==1 && cus[customer_index].dest.location==1)
		{
			this->Avg_TT_0[0].PT+=cus[customer_index].TTT.PT;
			this->Avg_TT_0[0].RS+=cus[customer_index].TTT.RS;
            this ->category_zone[0][cus[customer_index].RT]++;
			customer_type_zone[0]++;
		}
		if(cus[customer_index].orig.location==1 && cus[customer_index].dest.location==2)
		{
			this->Avg_TT_0[1].PT+=cus[customer_index].TTT.PT;
			this->Avg_TT_0[1].RS+=cus[customer_index].TTT.RS;
			this ->category_zone[1][cus[customer_index].RT]++;
			customer_type_zone[1]++;
		}
		if(cus[customer_index].orig.location==1 && cus[customer_index].dest.location==3)
		{
			this->Avg_TT_0[2].PT+=cus[customer_index].TTT.PT;
			this->Avg_TT_0[2].RS+=cus[customer_index].TTT.RS;
			this ->category_zone[2][cus[customer_index].RT]++;
			customer_type_zone[2]++;
		}
		if(cus[customer_index].orig.location==2 && cus[customer_index].dest.location==1)
		{
			this->Avg_TT_0[3].PT+=cus[customer_index].TTT.PT;
			this->Avg_TT_0[3].RS+=cus[customer_index].TTT.RS;
			this ->category_zone[3][cus[customer_index].RT]++;
			customer_type_zone[3]++;
		}
		if(cus[customer_index].orig.location==2 && cus[customer_index].dest.location==2)
		{
			this->Avg_TT_0[4].PT+=cus[customer_index].TTT.PT;
			this->Avg_TT_0[4].RS+=cus[customer_index].TTT.RS;
			this ->category_zone[4][cus[customer_index].RT]++;
			customer_type_zone[4]++;
		}
		if(cus[customer_index].orig.location==2 && cus[customer_index].dest.location==3)
		{
			this->Avg_TT_0[5].PT+=cus[customer_index].TTT.PT;
			this->Avg_TT_0[5].RS+=cus[customer_index].TTT.RS;
			this ->category_zone[5][cus[customer_index].RT]++;
			customer_type_zone[5]++;
		}
		if(cus[customer_index].orig.location==3 && cus[customer_index].dest.location==1)
		{
			this->Avg_TT_0[6].PT+=cus[customer_index].TTT.PT;
			this->Avg_TT_0[6].RS+=cus[customer_index].TTT.RS;
			this ->category_zone[6][cus[customer_index].RT]++;
			customer_type_zone[6]++;
		}
		if(cus[customer_index].orig.location==3 && cus[customer_index].dest.location==2)
		{
			this->Avg_TT_0[7].PT+=cus[customer_index].TTT.PT;
			this->Avg_TT_0[7].RS+=cus[customer_index].TTT.RS;
			this ->category_zone[7][cus[customer_index].RT]++;
			customer_type_zone[7]++;
		}
		if(cus[customer_index].orig.location==3 && cus[customer_index].dest.location==3)
		{
			this->Avg_TT_0[8].PT+=cus[customer_index].TTT.PT;
			this->Avg_TT_0[8].RS+=cus[customer_index].TTT.RS;
			this ->category_zone[8][cus[customer_index].RT]++;
			customer_type_zone[8]++;
		}
	}

	for(int i=0; i<9; i++){
		this->Avg_TT_0[i].PT=this->Avg_TT_0[i].PT*1/customer_type_zone[i];
		this->Avg_TT_0[i].RS=this->Avg_TT_0[i].RS*1/customer_type_zone[i];
		//cout <<"\n Pop["<<index_network<<"].Avg_TT_0["<<i<<"].PT="<<Pop.at(index_network).Avg_TT_0[i].PT<<"  RS="<<Pop.at(index_network).Avg_TT_0[i].RS;
	}
	this->check_consistency();
	return this->ftn;

}

void generate_initial_population(Global_solution& G_best, vector<Global_solution>& Pop,
	vector<Global_solution>& P_best_, const vector<stop>& set_of_stations){
	int rnd;
	for(int particle_index=1; particle_index<Consts::population_size+1; particle_index++){
	  Pop.at(particle_index).nume = particle_index;
				cout<<"\n line \t"<<__LINE__<<endl;

	  for(int line_index=0; line_index<number_of_lines; line_index++){
		Pop.at(particle_index).L.at(line_index).NV= LIN[line_index].NV;
		int nbr_stat = LIN[line_index].get_nbr_stat(); // Pop.at(particle_index).L.at(line_index).nbr_stat=LIN[line_index].nbr_stat; change for xiaoyi
		rnd=rand()%nbr_stat; // generate a random number between(0, nbr_stat-1)
		cout<<"\n line \t"<<__LINE__<<endl;
		for(int station_index=0; station_index< LIN[line_index].s.size(); station_index++) // change for xiaoyi{
		{
			Pop.at(particle_index).L.at(line_index).active[station_index]=LIN[line_index].active[station_index];
			Pop.at(particle_index).L.at(line_index).s[station_index]=LIN[line_index].s[station_index];
		}
	/*
		for(k=0; k< LIN[j].nbr_stat; k++)
		{
			Pop.at(i).L.at(j).active[k]=LIN[j].active[k];
			Pop.at(i).L.at(j).s[k]=LIN[j].s[k];
		}
	*/
	
		Pop.at(particle_index).L.at(line_index).active[rnd]=rand()%2;
		Pop.at(particle_index).L.at(line_index).Min_NV= Pop.at(particle_index).Generate_NV(line_index, minHT, set_of_stations);
		Pop.at(particle_index).L.at(line_index).Max_NV= Pop.at(particle_index).Generate_NV(line_index, maxHT, set_of_stations);
				cout<<"\n line \t"<<__LINE__<<endl;
				cout<<"\n Pop.at(particle_index).L.at(line_index).Min_NV \t"<<Pop.at(particle_index).L.at(line_index).Min_NV<<endl;
				cout<<"\n Pop.at(particle_index).L.at(line_index).Max_NV \t"<<Pop.at(particle_index).L.at(line_index).Max_NV<<endl;
		int difference = (Pop.at(particle_index).L.at(line_index).Min_NV-Pop.at(particle_index).L.at(line_index).Max_NV);
		double rd;
		if(difference ==0){
			rd= rand()%Pop.at(particle_index).L.at(line_index).Max_NV;
		}
		else{
			rd=rand()%difference+Pop.at(particle_index).L.at(line_index).Max_NV;
		}

		Pop.at(particle_index).CH[line_index].pos=rnd;
				cout<<"\n line \t"<<__LINE__<<endl;

		Pop.at(particle_index).CH[line_index].vlr=Pop.at(particle_index).L.at(line_index).active[rnd];
		Pop.at(particle_index).calculate_frequency(line_index, set_of_stations);
	  } //end of line index iteration 

		P_best_.push_back(Pop.at(particle_index));  
	} //end of particle index iteration

	G_best=Pop.at(0);
	G_best.ftn=Pop.at(0).ftn;
	G_best.ATT=Pop.at(0).ATT;
	G_best.obj.fit_1=Pop.at(0).obj.fit_1;
	G_best.obj.fit_2=Pop.at(0).obj.fit_2;
	G_best.nume=0;
}

void compute_customer_amount(vector<Global_solution>& Pop, const vector<stop>& set_of_stations){
	for(int particle_index=1; particle_index< Consts::population_size+1; particle_index++){
		Pop.at(particle_index).RideSharing_Users=0;
		Pop.at(particle_index).R_w_pt_rs_counter=0;
		Pop.at(particle_index).R_rs_pt_w_counter=0;
		Pop.at(particle_index).R_pt_counter=0;
		Pop.at(particle_index).R_rs_counter=0;
		Pop.at(particle_index).Update_nbr_usrs();
		
		for(int customer_index=0; customer_index<Pop.at(particle_index).cus.size(); customer_index++){
			Pop.at(particle_index).Customer_category(customer_index);
			Pop.at(particle_index).Define_customer_type(customer_index,set_of_stations);
			if(Pop.at(particle_index).cus[customer_index].RT==0 || Pop.at(particle_index).cus[customer_index].RT==1 || Pop.at(particle_index).cus[customer_index].RT==4)
				Pop.at(particle_index).RideSharing_Users++;
			switch (Pop.at(0).Category.at(customer_index)){
				case R_W_PT_RS:
					Pop.at(particle_index).R_w_pt_rs_counter++;
					break;
				case R_RS_PT_W:
					Pop.at(particle_index).R_rs_pt_w_counter++;
					break;
				case R_PT:
					Pop.at(particle_index).R_pt_counter++;
					break;
				case R_RS:
					Pop.at(particle_index).R_rs_counter++;
					break;
			}
		}
	} 
}


int Print_Number_of_Vehicles(int iter, const vector<Global_solution>& Pop,const Global_solution& G_best){
	string const nomFichier(".\\result\\Graphic1.txt");
	if(iter==-1){
		ofstream monFlux(nomFichier.c_str());
		if(monFlux){	
			monFlux<<"\n"<<"vehicle cost for the global best particle"<<endl;
			monFlux <<"\n"<<"iter \t Indice_best \t TN_Cost \t PT_Veh \t RS_Veh \t headway" << endl;
			monFlux<<"\n"<<-1<<"\t"<<0<<"\t"<<Pop.at(0).ftn<<"\t"<<Pop.at(0).obj.fit_1<<"\t"<<Pop.at(0).obj.fit_2<<endl;
			monFlux<<"\n"<<-1<<"\t"<<G_best.nume<<"\t"<<G_best.ftn<<"\t"<<G_best.obj.fit_1<<"\t"<<G_best.obj.fit_2<<endl;
		}else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}else{
		ofstream monFlux(nomFichier.c_str(), ios::app);
		if(monFlux)
			monFlux<<"\n"<<iter<<"\t"<<G_best.nume<<"\t"<<G_best.ftn<<"\t"<<G_best.obj.fit_1<<"\t"<<G_best.obj.fit_2<<endl;
		else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}
	return 0;
}

int Print_Average_Travel_Time(int iter, const vector<Global_solution>& Pop,const Global_solution& G_best){
	string const nomFichier(".\\result\\Graphic2.txt");
	if(iter==-1){
		ofstream monFlux(nomFichier.c_str());
		if(monFlux){	
			monFlux << "\n iter \t index \t Average Travel Time" << endl;
			monFlux<<"\n"<<-1<<"\t"<<0<<"\t"<< Pop.at(0).ATT<<endl;
			monFlux<<"\n"<<-1<< "\t"<<G_best.nume<< G_best.ATT<<endl;
		}else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}else{
		ofstream monFlux(nomFichier.c_str(), ios::app);
		if(monFlux)
			monFlux<<"\n"<<iter<< "\t"<<G_best.nume<< G_best.ATT<<endl;
		else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}
	return 0;
}

int Print_Number_of_vehicle_by_Line(int iter, const vector<Global_solution>& Pop,const Global_solution& G_best){
	string const nomFichier(".\\result\\Graphic3.txt");
	if(iter==-1){
		ofstream monFlux(nomFichier.c_str());
		if(monFlux){	
			monFlux<<" \n \n \n Number of vehicles in each line upon iterations \n ";
        	monFlux<<" \n iteration \t index\t line 1 \t line 2 \t line 3 \t line 4 \t line 5 \t line 6 \t line 7 ";		
			monFlux<<"\n"<<-1<<"\t"<<0<<"\t"<<Pop.at(0).L[0].NV<<"\t"<<Pop.at(0).L[1].NV<<"\t"<<Pop.at(0).L[2].NV<<"\t"<<Pop.at(0).L[3].NV<<"\t"<<Pop.at(0).L[4].NV<<"\t"<<Pop.at(0).L[5].NV<<"\t"<<Pop.at(0).L[6].NV<<endl;
			monFlux<<"\n"<<-1<<"\t"<<G_best.nume<<"\t"<<G_best.L[0].NV<<"\t"<<G_best.L[1].NV<<"\t"<<G_best.L[2].NV<<"\t"<<G_best.L[3].NV<<"\t"<<G_best.L[4].NV<<"\t"<<G_best.L[5].NV<<"\t"<<G_best.L[6].NV<<endl;
		}else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}else{
		ofstream monFlux(nomFichier.c_str(), ios::app);
		if(monFlux)
			monFlux<<"\n"<<iter<<"\t"<<G_best.nume<<"\t"<<G_best.L[0].NV<<"\t"<<G_best.L[1].NV<<"\t"<<G_best.L[2].NV<<"\t"<<G_best.L[3].NV<<"\t"<<G_best.L[4].NV<<"\t"<<G_best.L[5].NV<<"\t"<<G_best.L[6].NV<<endl;
		else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}
	return 0;
}

int print_particle_performance(int iter, const vector<Global_solution>& Pop){
	string const nomFichier(".\\result\\Graphic12.txt");
	if(iter==-1){
		ofstream monFlux(nomFichier.c_str());
		if(monFlux){
			monFlux << "iter \t particle_index \t TN_Cost \t PT_Veh \t RS_Veh" << endl;
			for(int particle_index=0; particle_index< Consts::population_size+1; particle_index++)
				monFlux<<"\n"<<-1<< "\t"<<particle_index<< "\t"<< Pop.at(particle_index).ftn<<"\t"<<Pop.at(particle_index).obj.fit_1<< "\t"<<Pop.at(particle_index).obj.fit_2;
		}else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}else{
		ofstream monFlux(nomFichier.c_str(), ios::app);
		if(monFlux){
			for(int particle_index=1; particle_index< Consts::population_size+1; particle_index++)
				monFlux<<"\n"<<iter<< "\t"<<particle_index<< "\t"<< Pop.at(particle_index).ftn<<"\t"<<Pop.at(particle_index).obj.fit_1<< "\t"<<Pop.at(particle_index).obj.fit_2;
		}else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}
	return 0;
}


int Visualize_Number_of_users_of_each_line(int iter, const vector<Global_solution>& Pop,const Global_solution& G_best){
	string const nomFichier(".\\result\\Graphic8.txt");
	if(iter==-1){
		ofstream monFlux(nomFichier.c_str());
		if(monFlux){
			monFlux<<" \n \n \n Nbr of users of each line \n \n \n";
			monFlux<<" \n iteration \t particle \t line 1 \t line 2 \t line 3 \t line 4 \t line 5 \t line 6 \t line 7 "<<endl;
			monFlux<<"\n"<<-1<<"\t"<<Pop.at(0).L.at(0).Nbr_users<<"\t"<< Pop.at(0).L.at(1).Nbr_users<<"\t"<<Pop.at(0).L.at(2).Nbr_users<<" \t"<< Pop.at(0).L.at(3).Nbr_users<<"\t"<<Pop.at(0).L.at(4).Nbr_users<<"\t"<< Pop.at(0).L.at(5).Nbr_users<<"\t"<< Pop.at(0).L.at(6).Nbr_users<<endl;
			monFlux<<"\n"<<-1<<"\t"<<G_best.nume<<"\t"<<G_best.L.at(0).Nbr_users<<"\t"<< G_best.L.at(1).Nbr_users<<"\t"<<G_best.L.at(2).Nbr_users<<" \t"<< G_best.L.at(3).Nbr_users<<"\t"<<G_best.L.at(4).Nbr_users<<"\t"<< G_best.L.at(5).Nbr_users<<"\t"<< G_best.L.at(6).Nbr_users<<endl;
		}else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}else{
		ofstream monFlux(nomFichier.c_str(), ios::app);
		if(monFlux){
			monFlux<<"\n";
	   		monFlux<<"\n"<<iter<<"\t"<<G_best.nume<<"\t"<<G_best.L.at(0).Nbr_users<<"\t"<< G_best.L.at(1).Nbr_users<<"\t"<<G_best.L.at(2).Nbr_users<<" \t"<< G_best.L.at(3).Nbr_users<<"\t"<<G_best.L.at(4).Nbr_users<<"\t"<< G_best.L.at(5).Nbr_users<<"\t"<< G_best.L.at(6).Nbr_users<<endl;
		}else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}
	return 0;
}

int Visualize_number_passenger_category_zone(int iter, const vector<Global_solution>& Pop,const Global_solution& G_best){
	string const nomFichier(".\\result\\Graphic19.txt");
	if(iter==-1){
		ofstream monFlux(nomFichier.c_str());
		if(monFlux){
			monFlux<<"\n iteration "<<-1<<"particle"<<0<<endl;
			for (int k = 0; k < Number_categories; k++){
				monFlux<<"\n Category "<<k<<":\t";
				for (int j = 0; j < 9; j++)    // 9  because we have 3 (P, PC,GC) zone so the different possible movement between zones are 9
					monFlux<< Pop.at(0).category_zone[j][k]<<"\t";
			}
			monFlux<<"\n iteration "<<-1<<"particle"<<G_best.nume<<endl;
			for ( int k = 0; k < Number_categories; k++){
				monFlux<<"\n Category "<<k<<":\t";
				for (int j = 0; j < 9; j++)    // 9  because we have 3 (P, PC,GC) zone so the different possible movement between zones are 9
					monFlux<< G_best.category_zone[j][k]<<"\t";
			}
		}else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}else{
		ofstream monFlux(nomFichier.c_str(), ios::app);
		if(monFlux){
			monFlux<<"\n iteration "<<iter<<"particle"<<G_best.nume<<endl;
			for (int k = 0; k < Number_categories; k++){
				monFlux<<"\n Category "<<k<<":\t";
				for (int j = 0; j < 9; j++)    // 9  because we have 3 (P, PC,GC) zone so the different possible movement between zones are 9
					monFlux<< G_best.category_zone[j][k]<<"\t";
			}
		}else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}
	return 0;
}


int Print_Active_Stations_in_Gbest_by_Line(int iter, const vector<Global_solution>& Pop, const Global_solution& G_best){
	string const nomFichier(".\\result\\Graphic15.txt");
	if(iter==-1){
		ofstream monFlux(nomFichier.c_str());
		if(monFlux){
			monFlux<<" \n \n \n Number of active stations in each line upon iterations \n ";
	   		monFlux<<" \n iteration \t line 1 \t line 2 \t line 3 \t line 4 \t line 5 \t line 6 \t line 7 \n  ";
	   		monFlux<<"\n"<<-1<<"\t"<<Pop.at(0).L.at(0).get_nbr_stat()<<"\t"<<Pop.at(0).L.at(1).get_nbr_stat()<<"\t"<<Pop.at(0).L.at(2).get_nbr_stat()<<" \t"<< Pop.at(0).L.at(3).get_nbr_stat()<<"\t"<<Pop.at(0).L.at(4).get_nbr_stat()<<"\t"<< Pop.at(0).L.at(5).get_nbr_stat()<<"\t"<< Pop.at(0).L.at(6).get_nbr_stat()<<endl;
	   		monFlux<<"\n"<<-1<<"\t"<<G_best.L.at(0).get_nbr_stat()<<"\t"<<G_best.L.at(1).get_nbr_stat()<<"\t"<<G_best.L.at(2).get_nbr_stat()<<" \t"<< G_best.L.at(3).get_nbr_stat()<<"\t"<<G_best.L.at(4).get_nbr_stat()<<"\t"<< G_best.L.at(5).get_nbr_stat()<<"\t"<< G_best.L.at(6).get_nbr_stat()<<endl;
		}
		else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}else{
		ofstream monFlux(nomFichier.c_str(), ios::app);
		if(monFlux){
	   		monFlux<<"\n"<<iter<<"\t"<<G_best.L.at(0).get_nbr_stat()<<"\t"<<G_best.L.at(1).get_nbr_stat()<<"\t"<<G_best.L.at(2).get_nbr_stat()<<" \t"<< G_best.L.at(3).get_nbr_stat()<<"\t"<<G_best.L.at(4).get_nbr_stat()<<"\t"<< G_best.L.at(5).get_nbr_stat()<<"\t"<< G_best.L.at(6).get_nbr_stat()<<endl;
		}else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}
	return 0;
}

int Print_Number_Riders_by_Ride_Type(int iter, const vector<Global_solution>& Pop, const Global_solution& G_best){
	string const nomFichier(".\\result\\Graphic6.txt");
	if(iter==-1){
		ofstream monFlux(nomFichier.c_str());
		if(monFlux){
			monFlux<<" \n \n \n Number of riders in each type upon each iteration \n \n \n";
	   		monFlux<<" \n  iteration\t index \t R_W_PT_RS_\t R_RS_PT_W \t R_PT \t R_RS\n";
	   		monFlux<<"\n"<<-1 <<"\t"<< 0<<"\t"<<Pop.at(0).R_w_pt_rs_counter<<"\t"<<Pop.at(0).R_rs_pt_w_counter<<"\t"<<Pop.at(0).R_pt_counter<<"\t"<<Pop.at(0).R_rs_counter<<"\t";
	   		monFlux<<"\n"<<-1 <<"\t"<< G_best.nume<<"\t"<<G_best.R_w_pt_rs_counter<<"\t"<<G_best.R_rs_pt_w_counter<<"\t"<<G_best.R_pt_counter<<"\t"<<G_best.R_rs_counter<<"\t";
		}
		else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}else{
		ofstream monFlux(nomFichier.c_str(), ios::app);
		if(monFlux){
	   		monFlux<<"\n";
			monFlux<<"\n"<<iter <<"\t"<< G_best.nume<<"\t"<<G_best.R_w_pt_rs_counter<<"\t"<<G_best.R_rs_pt_w_counter<<"\t"<<G_best.R_pt_counter<<"\t"<<G_best.R_rs_counter<<"\t";
		}else{
			cout << "ERREUR: Impossible d'ouvrir le fichier. ligne="<<__LINE__ << endl;
			exit(0);
		}
	}
	return 0;
}

void print_results(int iter, const vector<Global_solution>& Pop, const Global_solution& G_best){
	Print_Number_of_Vehicles(iter, Pop, G_best);
	Print_Average_Travel_Time(iter, Pop, G_best);
	Print_Number_of_vehicle_by_Line(iter, Pop, G_best);
	print_particle_performance(iter, Pop);
	Visualize_Number_of_users_of_each_line(iter, Pop, G_best);
	Visualize_number_passenger_category_zone(iter, Pop, G_best);
	Print_Active_Stations_in_Gbest_by_Line(iter, Pop,G_best);
	Print_Number_Riders_by_Ride_Type(iter, Pop,G_best);
}


float segm(float velo)
{
	return 1./(1+exp(-velo));
}


/*************************algorithm 2, change the activation/deativation of stops*************************/
void Global_solution::proba(int particle, int line_index, int stop_index, Global_solution& G_best, 
	const vector<Global_solution>& P_best_,int iter){
	cout<<"\n line \t"<<__LINE__<<endl;

	double r_1,r_2,d1_1,d1_0,d2_1,d2_0,inertia,V,r;
	r_1=(rand()%100)*1./100;
    r_2=(rand()%100)*1./100;
	//check for xiaoyi, check the P_best_[particle].L[line_index].s[stop_index] >=0
    if(P_best_[particle].L[line_index].active[stop_index]==1){
    	d1_1= C_1*r_1;
        d1_0= -C_1*r_1;
    }
    else{
    	d1_0= C_1*r_1;
        d1_1= -C_1*r_1;
    }
    if(G_best.L[line_index].active[stop_index]==1){
    	d2_1= C_2*r_2;
    	d2_0= -C_2*r_2;
    }
    else{
    	d2_0= C_2*r_2;
    	d2_1= -C_2*r_2;
    }
    inertia=((rand()%200)-100)*1./100;
    if(iter==0){
    	this->L[line_index].V_1[stop_index] = (rand()%100)*1./100;
    	this->L[line_index].V_0[stop_index] = (rand()%100)*1./100;
    }
    this->L[line_index].V_1[stop_index] = inertia*this->L[line_index].V_1[stop_index]+d1_1+d2_1;
    this->L[line_index].V_0[stop_index] = inertia*this->L[line_index].V_0[stop_index]+d1_0+d2_0;
    if(this->L[line_index].active[stop_index]==1) //change for xiaoyi: if(this->L[line_index].active[stop_index]==0)
    	V=this->L[line_index].V_1[stop_index];
    else
    	V=this->L[line_index].V_0[stop_index];
    r=(rand()%100)*1./100;
    if(r<segm(V))
    {
    	if (this->L[line_index].active[stop_index]==0 )
    		this->L[line_index].active[stop_index]=1;
    	else
    		this->L[line_index].active[stop_index]=0;
    }
    this->check_consistency();
 }//proba

 //Algorithm 3: change the number of buses per line
void Global_solution::DPSO(int particle_index, const Global_solution& G_best,
	const vector<Global_solution>& P_best_,
	const vector<stop>& set_of_stations, float CR1, double CR2, double CR3)         // indice de Pbest et indice de Gbest
{
	this->check_consistency();
	int auxiliary_1[number_of_lines]; 
	int auxiliary_2[number_of_lines];
	int auxiliary_3[number_of_lines];
	int line_index;
	double ac,ac1;
	for(line_index=0; line_index<number_of_lines; line_index++){
		ac=rand()%100*1./100;
		ac1=round(CR1*ac);
						cout<<"\n line \t"<<__LINE__<<endl;

		// Spread of partial best solution within the swarm:
		if (ac1==0)
			auxiliary_1[line_index]=this->L[line_index].NV; //auxiliary_1 is equal to the current particle p
		else
			auxiliary_1[line_index]=P_best_[particle_index].L[line_index].NV; //auxiliary_1 is equal to the local best version
								
								cout<<"\n line \t"<<__LINE__<<endl;

		//Spread of global solution within the swarm:
		ac=rand()%100*1./100;
		ac1=round(CR2*ac);
		if(ac1==0)
			auxiliary_2[line_index]=auxiliary_1[line_index];
		else
			auxiliary_2[line_index]=G_best.L[line_index].NV;
								
								cout<<"\n line \t"<<__LINE__<<endl;

		//Spread of randomness within the swarm:
		ac=rand()%100*1./100;
		ac1=round(CR3*ac);
		if(ac1==0)
				auxiliary_3[line_index]=auxiliary_2[line_index]; //todo for xiaoyi:mismatch with the pseudocode, N_p^{p}= =auxiliary_2[line_index] for all lines
		else{
			this->L[line_index].Min_NV = this->Generate_NV(line_index,minHT, set_of_stations);
			this->L[line_index].Max_NV= this->Generate_NV(line_index,maxHT, set_of_stations);
			int difference = (this->L[line_index].Min_NV-this->L[line_index].Max_NV);
			double rd;
			if(difference ==0){
				rd= rand()%this->L[line_index].Max_NV;
			}
			else{
				rd=rand()%difference+this->L[line_index].Max_NV;
			}
			auxiliary_3[line_index]=rd;
			}
		this->L[line_index].NV=auxiliary_3[line_index];
								cout<<"\n line \t"<<__LINE__<<endl;

		this->check_consistency();
	}
} //DPSO

void compteur(const Global_solution& sl, int iter){
	int a,type_1=0,type_2=0,type_3=0,type_4=0;
	for(a=0; a<sl.cus.size(); a++){
		if(sl.Category.at(a)==0)
			type_1++;
		else
			if(sl.Category.at(a)==1)
				type_2++;
			else
				if(sl.Category.at(a)==3)
					type_3++;
				else
					if(sl.Category.at(a)==4)
						type_4++;
	}

	TypeRider[iter][0]=type_1;
	TypeRider[iter][1]=type_2;
	TypeRider[iter][2]=type_3;
	TypeRider[iter][3]=type_4;
	// TypeRider[iter][4]=type_4;
}


void Global_solution:: compteur_2(int iter)
{
	int a,type_1=0,type_2=0,type_3=0,type_4=0;
	for(a=0; a<cus.size(); a++)
		if(this->cus.at(a).RT ==0)
			type_1++;
		else
			if(this->cus.at(a).RT ==1)
			type_2++;
			else
				if(this->cus.at(a).RT==3)
					type_3++;
				else
					if(this->cus.at(a).RT==4)
						type_4++;

	TypeRider2[iter][0]=type_1++;
	TypeRider2[iter][1]=type_2++;
	TypeRider2[iter][2]=type_3++;
	TypeRider2[iter][3]=type_4++;

}


/**************************** Update of g_best ******************************/
tuple<int, double, int> Update_G_best (int iter, Global_solution& G_best, const vector<Global_solution>& Pop){
	int i,j;
	double sum = 0, mean =0;
	int LastValueGbest =0;
	int Act_best=Pop[1].ftn;
	int L=G_best.ftn;
	for(i=1; i<Consts::population_size+1; i++){
		if(Pop.at(i).ftn < G_best.ftn )//here G_best is the G_best in 10 particles in iteration -1 
			G_best=Pop.at(i);
		if(Pop.at(i).ftn < Act_best)
			Act_best=Pop.at(i).ftn;
		sum+=Pop.at(i).ftn;
	}
    mean=sum/Consts::population_size;

	if (G_best.ftn == L)
		LastValueGbest++;
	else
		LastValueGbest=0;

	for(i=0; i<number_of_lines; i++){
		HeadwayTime[iter][i]=G_best.headway_time[i];
	}
	compteur(G_best,iter);
	G_best.compteur_2(iter);
	cout<<"\n check for xiaoyi 1:iter \t"<<iter<<" mean \t"<<mean<<"\t LastValueGbest \t"<<LastValueGbest<<"\t Act_best \t "<<Act_best<<endl;
	// cout<<"\n check for xiaoyi TypeRider[iter][0] \t"<< TypeRider[iter][0]<<"\t TypeRider[iter][1] \t"<< TypeRider[iter][1]<<"\t TypeRider[iter][2] \t "<< TypeRider[iter][2]<<"\t TypeRider[iter][3] \t"<< TypeRider[iter][3]<<endl;
	// cout<<"\n check for xiaoyi TypeRider2[iter][0] \t"<< TypeRider2[iter][0]<<"\t TypeRider2[iter][1] \t"<< TypeRider2[iter][1]<<"\t TypeRider2[iter][2] \t "<< TypeRider2[iter][2]<<"\t TypeRider2[iter][3] \t"<< TypeRider2[iter][3]<<endl;
	return {mean, LastValueGbest, Act_best};
}

/******* change CR1, CR2, CR3 ********/
void change(int iter, double mean, int LastValueGbest, int Act_best, 
			float CR1, double CR2, double CR3 ){
	cout<<"\n check for xiaoyi 2:iter \t"<<iter<<" mean \t"<<mean<<"\t LastValueGbest \t"<<LastValueGbest<<"\t Act_best \t "<<Act_best<<endl;
	if(LastValueGbest>=10)
		CR1+=CR1*36/100;
	if(CR1 > 0.55+0.55*37/100)
		CR1=0.55;
	if(mean-Act_best< (20*Act_best)/100 && LastValueGbest==15)
		CR2-=CR2*31/100;
	if(CR2 < 0.65-0.65*7.7/100)
		CR2=0.65;
	if(LastValueGbest>=20)
		CR3+=CR3*3.7/100;
	if(CR3 > 0.53+0.53*22.6/100)
		CR3=0.53;
}

void BPSO(const Network &n, Global_solution& G_best, vector<Global_solution>& Pop, vector<Global_solution>& P_best_,
	const vector<stop>& set_of_stations, const vector<line> LIN_, float CR1, double CR2, double CR3,
	const vector< vector<float> >& travel_time_on_arc){
		cout<<"\n line \t"<<__LINE__<<endl;

		generate_initial_population(G_best, Pop, P_best_, set_of_stations);
		cout<<"\n line \t"<<__LINE__<<endl;

		for(int particle_index=1; particle_index<Consts::population_size+1; particle_index++){
			// cout<<"\n particle_index \t"<<particle_index<<endl;
			for(int line_index=0; line_index<number_of_lines; line_index++){
				// cout<<"\n line_index \t"<<line_index<<endl;
				// cout<<"\n line \t"<<__LINE__<<endl;

				// cout<<"\n  Pop.at(0).L.at(line_index \t"<< Pop.at(0).L.at(line_index).active[1]<<endl;

				P_best_.at(particle_index).L.at(line_index)= Pop.at(0).L.at(line_index);
				// cout<<"\n  Pop.at(0).L.at(line_index \t"<<P_best_.at(particle_index).L.at(line_index).active[1]<<endl;
				P_best_.at(particle_index).check_consistency();
			}
			cout<<"\n line \t"<<__LINE__<<endl;
			if (Consts::debug_mode){
					int last_element_added = P_best_.size()-1;
					if ( &(P_best_.at(last_element_added) ) == &(Pop.at(particle_index)) )
						{
							std::stringstream err_msg;
							err_msg<<"\nLine "<<__LINE__<<": Error: the addrees &(P_best_.at(last_elemtn_added)) is "<<&(P_best_.at(last_element_added))
										<<", which is the same of &(Pop.at(i))="<<&(Pop.at(particle_index))<<
										". last_element_added="<<last_element_added<<"; i="<<particle_index<<
										". Therefore, P_best_.at(last_element_added) and Pop.at(i) are "<<
										"the same object (what you wanted, instead, is to have two SEPARATE objects that look"
										<<" the same"<<endl;
							cout<<err_msg.str();
							throw err_msg.str();
						}
				}

			if (Consts::debug_mode)
			{
				if ( P_best_.size() != Pop.size() )
				{
							std::stringstream err_msg;
							err_msg<<"\nLine "<<__LINE__<<"Warning: P_best_.size()="<<P_best_.size()
										<<", Pop.size()="<<Pop.size()<<". Should they be the same instead?????????????????"
										<<" The fact that they are not the same will produce an error later."<<endl;
							cout<<err_msg.str();
							throw err_msg.str();
				}
			}

		}


	#if PARALLEL_MODE
	#pragma omp parallel for if (PARALLEL_MODE)
	#endif
		for(int particle_index=1; particle_index<Consts::population_size+1; particle_index++)
			Pop.at(particle_index).dijkstra(set_of_stations,Pop.at(particle_index).L, travel_time_on_arc);
			cout<<"\n line \t"<<__LINE__<<endl;

	// initialize solution of particle 1-10 before PSO 

	compute_customer_amount(Pop, set_of_stations);
			cout<<"\n line \t"<<__LINE__<<endl;

	#if PARALLEL_MODE
	#pragma omp parallel for if (PARALLEL_MODE)
	#endif
		for(int particle_index=1; particle_index< Consts::population_size+1; particle_index++){
			// cout<<"\n particle_index \t"<<particle_index<<endl;

			Pop.at(particle_index).Eval(n, set_of_stations,particle_index);
			Pop.at(particle_index).Initialise_velocity();
			for(int line_index=0; line_index< number_of_lines;line_index++){
				Pop.at(particle_index).calculate_frequency(line_index, set_of_stations);
				// cout<<"\n line_index\t"<<line_index<<endl;
				}

		}

	for(int i=1;i<Consts::population_size+1;i++){
		if(Pop.at(i).ftn < G_best.ftn)
			G_best=Pop.at(i);
	}

	cout<< "\n G_best before PSO process: "<< G_best.ftn <<" corresponding to particle : "<< G_best.nume;
	
	print_results(-1, Pop, G_best);

	// beginning of the epoches！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！1 
	int iter = 0;
	while(iter < PSO_Numbe_Iteration){
		cout<<"\n iteration \t"<<iter<<endl;

		auto t1 = high_resolution_clock::now();

		#if PARALLEL_MODE
		#pragma omp parallel for if (PARALLEL_MODE)
		#endif
			for(int particle_index=1; particle_index<Consts::population_size+1; particle_index++){
				// cout<<"\n particle_index \t"<<particle_index<<endl;
				for(int line_index=0; line_index< number_of_lines; line_index++){
					cout<<"\n line_index\t"<<line_index<<endl;
					cout<<"\n Pop[particle_index].L[line_index].s.size() \t"<<Pop[particle_index].L[line_index].s.size()<<endl;
					//todo for xiaoyi: for each stops of line l, should we remove -1 in s.at(i)?
					for(int station_index=0; station_index< Pop[particle_index].L[line_index].s.size(); station_index++){ //change for xiaoyi: Pop[particle_index].L[line_index].nbr_stat
						// cout<<"\n line \t"<<__LINE__<<endl;
						Pop.at(particle_index).proba(particle_index,line_index, station_index, G_best, P_best_,iter);//algorithm2
						// cout<<"\n station_index \t"<<station_index<<endl;
					}	
					cout<<"\n line \t"<<__LINE__<<endl;

					Pop.at(particle_index).calculate_frequency(line_index, set_of_stations);
					cout<<"\n line \t"<<__LINE__<<endl;
					Pop.at(particle_index).DPSO(particle_index, G_best, P_best_, set_of_stations, CR1, CR2, CR3);// algorithm3
					cout<<"\n line \t"<<__LINE__<<endl;

				}
			}

		#if PARALLEL_MODE
		#pragma omp parallel for if (PARALLEL_MODE)
		#endif
			for(int particle_index=1; particle_index<Consts::population_size+1; particle_index++){
				cout<<"\n particle_index \t"<<particle_index<<endl;
				Pop.at(particle_index).dijkstra(set_of_stations, LIN_, travel_time_on_arc);
			}
			//todo for xiaoyi: does this use LIN or Pop.at(i).L?
	
	compute_customer_amount(Pop, set_of_stations);

	#if PARALLEL_MODE
	#pragma omp parallel for if (PARALLEL_MODE)
	#endif
		for(int particle_index=1; particle_index< Consts::population_size+1; particle_index++){
			Pop.at(particle_index).write_ride_sharing_results(set_of_stations);
			// cout<<"\n particle_index \t"<<particle_index<<endl;

		}
	
	#if PARALLEL_MODE
	#pragma omp parallel for if (PARALLEL_MODE)
	#endif
		for(int particle_index=1; particle_index< Consts::population_size+1; particle_index++){
			Pop.at(particle_index).Eval(n, set_of_stations,  particle_index);
			if( Pop.at(particle_index).ftn < P_best_.at(particle_index).ftn)
				P_best_.at(particle_index) = Pop.at(particle_index);
		}
    
	cout << "\n iter \t"<<iter<<"\t before Update_G_best G_best= " << G_best.ftn;
	auto [Act_best, mean, LastValueGbest]= Update_G_best(iter, G_best, Pop);
	cout << "\n iter \t"<<iter<<"\t after Update_G_best G_best= " << G_best.ftn;
	print_results(iter, Pop, G_best);
	iter++;
	change(iter, mean, LastValueGbest, Act_best, CR1, CR2, CR3);
	
	auto t2 = high_resolution_clock::now();
	auto diff = duration_cast< duration < double > > (t2 - t1);
	cout<<"\n time use for  "<<iter<<"th  "<<"iteration is\t"<<diff.count()<<endl;
	}
}//BPSO


DiscretizationHandler::DiscretizationHandler(const std::string discretization_info_file){
	ifstream infile;
	infile.open (discretization_info_file);
	float Delta_x, Delta_y;
	while (infile >> Delta_x >> Delta_y)
	{
		delta_x=Delta_x;
		delta_y=Delta_x;
	}
	infile.close();
}

Global_solution::Global_solution(int number_of_stations_, const vector <customer> cus_,
	const vector<line> L_){
	ftn=-1;
	ATT=-1;
	nume=0;
	CTT_PT=0; 
	CTT_RS=0; 
	number_of_stations = number_of_stations_;
	for (int i=0; i<number_of_stations; i++){
		vector<float> row;
		for (int j=0; j<number_of_stations; j++)
			row.push_back(-1);
		DIStance.push_back(row);
	}

	cus = cus_;
	for (int i=0; i<cus.size(); i++){
		Category.push_back(cus_.at(i).RT);
		PT_Estimated_EDT_LAT_To_Station.push_back(-1);
	}
	this -> L = L_;
		
		// cout<<"line  \t"<<L.at(0).s.size()<<endl;
		// cout<<"line \t"<<L.at(0).s[0]<<endl;
			
	if (Consts::debug_mode) check_consistency();
}


int main(int argc, char *argv[]){
	
	if (argc != 5){
		cout<<"The number of parameters is not correct. You have instead to call the program as"<<endl;
		cout<<"      <name_of_the_executable> <list_of_customers_file> <locations_file> <discretization_file><solution_file>"<<endl;
		exit(1);
	}

    string Liste_Of_Customers_file = argv[1];
	std::string locations_file = argv[2];
	std::string discretization_info_file = argv[3];
	std::string solution_file = argv[4];

	double CR1 = 0.55;   // CR1 leads the particles towards the corresponding partial best solution found before. CR1 \in [0,1]
	double CR2 = 0.65;   // CR2 allows the particles to be influenced by the best global solution found before. CR2 \in [0,1]
	double CR3 = 0.53;  //CR3 allows the randomness degree

    load_state(solution_file, LIN);

	srand(time(0));
	Network n;

    vector<stop> set_of_stations(Number_stations);
	generation_stops(LIN, set_of_stations);//update set of stations with station ids and lines allocated to 

	read_data_stations(0, set_of_stations, Consts::stations_file); //update set of stations with id and loc
		

    vector <customer> customers_from_file = read_customers_from_file(Liste_Of_Customers_file); 

    Customer_LastMile_FirstMile(customers_from_file, set_of_stations); //customer's origin/destination distance to every station

    DiscretizationHandler discretizationHandler(discretization_info_file);
	
	read_data_stations(0, set_of_stations, Consts::stations_file); //update set of stations with id and loc

	//initialize the network with each customer's loc + 63 stops's loc + 1 depot loc and calculated their travel time matrix
	//time to depot is always 0
	//todo for xiaoyi: understanding the functionality of this function
    initialize_network(n, locations_file, Consts::stations_file, discretizationHandler);

	// exit(0);
	//calculate the travel time between one station and its following station on arc, other distance is infinity
	//without add the accelerating and decelerating time
	const vector< vector<float> > travel_time_on_arc = initialisation_of_travel_time_on_arc(set_of_stations);

	vector<Global_solution> Pop;
	for(int i=0; i <= Consts::population_size; i++){
		Global_solution glob_sol(Number_stations, customers_from_file, LIN);
		Pop.push_back(glob_sol);
		// cout<<" \n i \t"<<i<<endl;
		// for(int j: Pop[i].L[1].get_active_station_ids()){
		// 	cout<<j<<endl;
		// exit(0);
		// }
	}


	if (Pop.size() != Consts::population_size+1){
			std::stringstream err_msg;
			err_msg<<"\nLine "<<__LINE__<<": Error: Consts::population_size+1="<<
				Consts::population_size+1<<". Pop.size()="<< Pop.size() <<endl;
			cout<<err_msg.str();
			throw err_msg.str();
			exit(1);
	}

	
	//initialize each line's bus velocity to 0
	Pop.at(0).Initialise_velocity();

	//generate the number of buses on each line
	for(int line_index=0; line_index<number_of_lines; line_index++){
		Pop.at(0).L[line_index].NV = Pop.at(0).Generate_NV(line_index, Consts::initial_headway, set_of_stations);
	}

	cout<<"\n Starting -1 iteration of PSO."<<endl;

	Pop.at(0).dijkstra(set_of_stations, Pop.at(0).L, travel_time_on_arc);//calculate the shortest path from station to stations based on particle 0 solution, and store in list_of_Station, list_of_lines



	for(int stop_index_1 :Pop.at(0).list_of_stations[1][63]){
		cout<<"Pop.at(0).list_of_stations[1][12]   "<<stop_index_1<<endl;
	}
	for(int stop_index_1 :Pop.at(0).list_of_lines[1][63]){
		cout<<"Pop.at(0).list_of_lines[1][12]   "<<stop_index_1<<endl;
	}
	exit(0);


	Pop.at(0).Update_nbr_usrs();//initialize the number of buses on each line to 0
	
	cout<<"Pop.at(0).cus.size()\t"<<Pop.at(0).cus.size()<<endl;
	// exit(0);
	for(int c=0; c<Pop.at(0).cus.size(); c++){
		Pop.at(0).Customer_category(c);    //based customers' origin/destination to active stations' location to calculate the potential customer type 
		Pop.at(0).Define_customer_type(c, set_of_stations); //based on the DIStance vector solved by dijkstra, calculate actual customer's type
		if(Pop.at(0).cus[c].RT==0 || Pop.at(0).cus[c].RT==1 || Pop.at(0).cus[c].RT==4) {
			Pop.at(0).RideSharing_Users++;
		}
	}

	
	Pop.at(0).Compute_Customer_PT_trajectory(); //compute cust bus trajectory based on list_of_stations given by dijkstra algorithm 

	Pop.at(0).write_ride_sharing_results(set_of_stations);

	Pop.at(0).ATT=0;

	Pop.at(0).Eval(n, set_of_stations,0);

	// set initial G_best to index 0 
	Global_solution G_best(Number_stations, Pop.at(0).cus, Pop.at(0).L);

	// set initial P_best to index 0 
	vector<Global_solution> P_best; 
	P_best.push_back(Pop.at(0)); 

	BPSO(n, G_best, Pop, P_best, set_of_stations, LIN,
		CR1,  CR2,  CR3, travel_time_on_arc);

    cout<<"\n Finish!!!!!!!!!!!!"<<endl;

}




bool line::is_it_equal_to(const line& other) const{
	if (NV != other.NV) return false;
	if (nbr_bus != other.nbr_bus) return false;
	if (freq != other.freq) return false;
	if (waitingtime != other.waitingtime) return false;
	if (index != other.index) return false;
	if (Max_NV != other.Max_NV) return false;
	if (Min_NV != other.Min_NV) return false;
	if (Nbr_users != other.Nbr_users) return false;
	return true;
}

line::line(int NV_, vector<int> s_, vector<int> active_){
	this->NV = NV_;
	this->s = s_;
	this->active = active_;
	for (int i=0; i<s_.size(); i++){
		V_1.push_back(0);
		V_0.push_back(0);
	}
	this->nbr_bus=-1;
	this->freq=-1;
	this->waitingtime=-1;
	this->index=-1;
	this->Max_NV=-1;
	this->Min_NV=-1;
	this->Nbr_users = -1;
}