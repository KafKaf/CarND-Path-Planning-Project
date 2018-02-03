#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// Constants
const double TRAJECTORY_S_SPACING = 30.0; // spacing for trajectory polynomial
const double CAR_STEP_IN_SECONDS = 0.02; // car will visit the next point every x seconds
const double TARGET_VELOCITY_MPH = 49.5;
const double TARGET_VELOCITY_MPS = 49.5 * 0.44704;
const double LANE_WIDTH_IN_METERS = 4;
const double SAFE_LANE_CHANGE_SPEED_IN_MPH = 30;
const double SAFE_LANE_CHANGE_DISTANCE_IN_METERS = 30;
const int MIN_RUNS_BEFORE_LANE_CHANGE = 10;
const int NUM_OF_REFERENCE_POINTS = 50;
const double NORMAL_VELOCITY_CHANGE_IN_MPH = 0.224;
const double MEDIUM_VELOCITY_CHANGE_IN_MPH = 0.44;
const double HIGH_VELOCITY_CHANGE_IN_MPH = 1.0;
const double SAFE_DISTANCE_METERS = 40.0;
const double SAFE_DISTANCE_SECONDS = 3.0;
const double MEDIUM_SAFE_DISTANCE_SECONDS = 2.25;
const double NOT_SAFE_DISTANCE_SECONDS = 0.75;
const double SAFE_MPS_TO_PRESERVE_SPEED_FOR_SMALL_GAP = -3.0;
const int NUM_OF_LANES = 3;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// For converting between mph(miles per hour) and mps(meters per second)
double mph2mps(double mph) { return mph * 0.44704; }

// time to distance in seconds(assuming constant velocity) - x = v * t ==> t = x / v
double getTimeToDistance(double distance_meters, double speed_meters_per_second) {
    return distance_meters / speed_meters_per_second;
}

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Transform points to vehicle coordinates
void TransformPointsToVehicleCoordinates(vector<double> &ptsx, vector<double> &ptsy, double px, double py, double psi) {
	// pre compute
	double sin_reverse_psi = sin(-psi);
	double cos_reverse_psi = cos(-psi);

	for (int p = 0;p < ptsx.size(); p++) {
		// translation
		double x = ptsx[p] - px;
		double y = ptsy[p] - py;

		// rotation
		ptsx[p] = cos_reverse_psi * x - sin_reverse_psi * y;
		ptsy[p] = cos_reverse_psi * y + sin_reverse_psi * x;
	}
}

// Transform points to global coordinates
void TransformPointsToGlobalCoordinates(vector<double> &ptsx, vector<double> &ptsy, double px, double py, double psi) {
	// pre compute
	double sin_psi = sin(psi);
	double cos_psi = cos(psi);

	for (int p = 0;p < ptsx.size(); p++) {
		// translation
		double x = ptsx[p] + px;
		double y = ptsy[p] + py;

		// rotation
		ptsx[p] = cos_psi * x - sin_psi * y;
		ptsy[p] = cos_psi * y + sin_psi * x;
	}
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

double distance(double x, double y)
{
	return sqrt(x * x + y * y);
}

int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

// Inefficiency cost - penalizing lane velocity difference from car's target velocity, higher the difference ==> higher the cost
double inefficiencyCost(double target_velocity, double lane_velocity) {

    double diff_velocity = target_velocity - lane_velocity;

    if  (diff_velocity <= 0)
        return 0;

    return diff_velocity / target_velocity;
}

// Lane change cost - lane change is not a desired behavior unless needed, we want to penalize it, bigger the change ==> higher the cost
double laneChangeCost(int current_lane, int final_lane) {
    return pow(current_lane - final_lane, 2);
}

// Lane change times cost - we want to penalize frequent lane changes for safety and comfort, more frequent ==> higher the cost
double laneChangeTimesCost(int steps_since_last_lane_change, int current_lane, int final_lane) {
    return current_lane != final_lane ? 100 / steps_since_last_lane_change : 0.0;
}

// s distance cost - we prefer to drive in a free lane, penalize driving with car in front of us
double sDistanceCost(double s_diff) {
    if (s_diff < 100)
        return 1;

    return 0;
}

// slowing down cost - we don't like our car slowing down, penalize the lane that make us slow down
double slowingDownCost(bool is_slowing_down, int current_lane, int final_lane) {
    return (current_lane == final_lane && is_slowing_down) ? 1.0 : 0.0;
}

// Behavior planning cost function
double calculateBehaviorPlanningCost(int steps_since_last_lane_change, int current_lane, double target_velocity_mps,
                                     int final_lane, double final_lane_velocity_mps, double s_diff,
                                     bool is_slowing_down) {
    double cost1 = 5 * inefficiencyCost(target_velocity_mps, final_lane_velocity_mps);
    double cost2 = 1 * laneChangeCost(current_lane, final_lane);
    double cost3 = 1 * laneChangeTimesCost(steps_since_last_lane_change, current_lane, final_lane);
    double cost4 = 2 * sDistanceCost(s_diff);
    double cost5 = 2 * slowingDownCost(is_slowing_down, current_lane, final_lane);

    return cost1 + cost2 + cost3 + cost4 + cost5;
}

bool isCarInLane(int lane, double pred_car_d) {
    return pred_car_d > LANE_WIDTH_IN_METERS*lane && pred_car_d < LANE_WIDTH_IN_METERS*(lane+1);
}

int getCarLane(double d, double lane_width) {
    return (int)(d/lane_width);
}

double predictCarSInEndOfGeneratedTrajectory(double car_s, int prev_size, double speed_mps) {
    return car_s + (double)prev_size * CAR_STEP_IN_SECONDS * speed_mps;
}

// Acts as behavior planning layer - get lane summary statitics, returns best lane index, distance to car ahead in meters, speed of car ahead in mps.
tuple<int,double,double> getLaneSummaryStatistics(int steps_since_last_lane_change, double car_s, int prev_size,
                                                  int current_lane,
                                                  double ref_vel_mps, bool is_slowing_down,
                                                  const vector<vector<double>> &sensor_fusion) {

    double distance_to_car_ahead_meters = 1000;
    double speed_of_car_ahead_mps = 100;

    // calculate lanes speed average
    vector<int> lane_cars(NUM_OF_LANES);
    vector<double> lane_speed(NUM_OF_LANES);
    vector<double> lane_s_diff(NUM_OF_LANES);
    std::fill (lane_cars.begin(), lane_cars.end(), 0);
    std::fill (lane_speed.begin(), lane_speed.end(), 0.0);
    std::fill (lane_s_diff.begin(), lane_s_diff.end(), 1000.0);

    for (auto &sensors : sensor_fusion) {
        // car is in my lane
        double d = sensors[6];
        int check_car_lane = getCarLane(d, 4);
        if (check_car_lane < 0 || check_car_lane > 2) {
            continue;
        }
        double vx = sensors[3];
        double vy = sensors[4];
        double check_speed_mps = distance(vx,vy);
        double check_car_s = sensors[5];
        // predict the other car s into the end of trajectory
        check_car_s = predictCarSInEndOfGeneratedTrajectory(check_car_s, prev_size, check_speed_mps);
        // only care about cars that are in front of us
        double s_diff = check_car_s - car_s;
        if (s_diff > 0) {
            lane_cars[check_car_lane] += 1;
            lane_speed[check_car_lane] += check_speed_mps;
            if (s_diff < lane_s_diff[check_car_lane])
                lane_s_diff[check_car_lane] = s_diff;

            if (current_lane == check_car_lane) {
                // min distance and speed to car in front of us
                if (s_diff < distance_to_car_ahead_meters) {
                    distance_to_car_ahead_meters = s_diff;
                    speed_of_car_ahead_mps = check_speed_mps;
                }
            }
        }
    }

    // choose lane with lowest cost
    double lane_cost = 10000;
    int lane_index = current_lane;
    for (int i = 0; i < lane_cars.size(); i++) {
        int cars_in_lane = lane_cars[i];
        double avg_speed_mps = cars_in_lane == 0 ? ref_vel_mps : lane_speed[i] / cars_in_lane;
        double cost = calculateBehaviorPlanningCost(steps_since_last_lane_change, current_lane, ref_vel_mps, i,
                                                    avg_speed_mps, lane_s_diff[i], is_slowing_down);
        if (cost < lane_cost) {
            lane_cost = cost;
            lane_index = i;
        }
    }

    return make_tuple(lane_index, distance_to_car_ahead_meters, speed_of_car_ahead_mps);
}

bool isLaneChangeSafe(double car_speed, int runs_since_last_lane_change, int intended_lane, double intended_s
		, int prev_size, const vector<vector<double>> &sensor_fusion) {

	// don't turn if your speed is low or if you changed lane just a moment ago
	if (car_speed < SAFE_LANE_CHANGE_SPEED_IN_MPH && runs_since_last_lane_change < MIN_RUNS_BEFORE_LANE_CHANGE)
		return false;

	for (auto &sensors : sensor_fusion) {
		// car is in my lane
		double d = sensors[6];
		if (isCarInLane(intended_lane, d)) {
			double vx = sensors[3];
			double vy = sensors[4];
			double check_speed = distance(vx,vy);
			double check_car_s = sensors[5];

			// predict the other car s into the future of the lane change
			check_car_s = predictCarSInEndOfGeneratedTrajectory(check_car_s, prev_size, check_speed) + TRAJECTORY_S_SPACING;
			// check if other cars are far away from us
            if (fabs(intended_s - check_car_s) < SAFE_LANE_CHANGE_DISTANCE_IN_METERS)
			    return false;
		}
	}

	return true;
}

// calculate the new car velocity(mph) based on data from the car ahead
double getNewReferenceVelocityFrom(double ref_vel, double car_speed, double speed_of_car_ahead_mps, double distance_diff_meters) {

    double car_speed_mps = mph2mps(car_speed);
    double speed_diff_mps = car_speed_mps - speed_of_car_ahead_mps;

    // assuming constant RELATIVE velocity and distance
    double time_to_collision_seconds = getTimeToDistance(distance_diff_meters, speed_diff_mps);

    // need to slow down, we are too close and/or getting close fast
    if (speed_diff_mps > 0 && (distance_diff_meters < SAFE_DISTANCE_METERS || time_to_collision_seconds < SAFE_DISTANCE_SECONDS)) {
        // slow down stronger as we get closer to the car ahead
        if (time_to_collision_seconds > MEDIUM_SAFE_DISTANCE_SECONDS)
            ref_vel -= NORMAL_VELOCITY_CHANGE_IN_MPH;
        else if (time_to_collision_seconds > NOT_SAFE_DISTANCE_SECONDS)
            ref_vel -= MEDIUM_VELOCITY_CHANGE_IN_MPH;
        else
            ref_vel -= HIGH_VELOCITY_CHANGE_IN_MPH;
    // need to keep same pace, we are close, but the gap is widening
    } else if (distance_diff_meters < SAFE_DISTANCE_METERS && speed_diff_mps > SAFE_MPS_TO_PRESERVE_SPEED_FOR_SMALL_GAP) {
        // do nothing - keep same speed
    // safe to accelerate again if we are not at max speed
    } else if (ref_vel < TARGET_VELOCITY_MPH) {
        ref_vel += NORMAL_VELOCITY_CHANGE_IN_MPH;
    }

    // can't go backwards on highway
    if (ref_vel < 0)
        ref_vel = 0;

    return ref_vel;
}

// execute lane change if needed and if safe to do it
void executeLaneChangeWithPreconditions(int &lane, int &runs_since_last_lane_change, int new_lane, double car_speed,
										double car_s, int prev_size, const vector<vector<double>> &sensor_fusion) {
    // need to execute a lane change only if the best lane is different from our current lane
	if (new_lane != lane) {

        // should only make a single lane change, remember not doing a full change
		bool full_change = true;
		int lane_diff = new_lane - lane;

        // change new lane to the closest lane in the direction of change
		if (abs(lane_diff) > 1) {
			if (lane_diff > 0) {
				new_lane = lane + 1;
			} else {
				new_lane = lane - 1;
			}
			full_change = false;
		}

        // s of our vehicle in the end of the lane change
		double intended_s = car_s + TRAJECTORY_S_SPACING;
		std::cout << "Intended lane:" << new_lane << std::endl;

        // check if lane change is safe
		if (isLaneChangeSafe(car_speed, runs_since_last_lane_change, new_lane, intended_s, prev_size, sensor_fusion)) {

            lane = new_lane;
			std::cout << "Changed lane to:" << lane << std::endl;

			// initialize only when full change is complete, make it worthwhile to do it in next iteration
			if (full_change)
				runs_since_last_lane_change = 0;
		}
	}
}

double getMiddleOfLaneD(int lane) {
	return LANE_WIDTH_IN_METERS * (lane + 0.5);
}

int main() {
    uWS::Hub h;

    // Load up map values for waypoint's x,y,s and d normalized normal vectors
    vector<double> map_waypoints_x;
    vector<double> map_waypoints_y;
    vector<double> map_waypoints_s;
    vector<double> map_waypoints_dx;
    vector<double> map_waypoints_dy;

    // Waypoint map to read from
    string map_file_ = "../data/highway_map.csv";
    // The max s value before wrapping around the track back to 0
    double max_s = 6945.554;

    ifstream in_map_(map_file_.c_str(), ifstream::in);

    string line;
    while (getline(in_map_, line)) {
  	    istringstream iss(line);
  	    double x;
  	    double y;
  	    float s;
  	    float d_x;
  	    float d_y;
  	    iss >> x;
  	    iss >> y;
  	    iss >> s;
  	    iss >> d_x;
  	    iss >> d_y;
  	    map_waypoints_x.push_back(x);
  	    map_waypoints_y.push_back(y);
  	    map_waypoints_s.push_back(s);
  	    map_waypoints_dx.push_back(d_x);
  	    map_waypoints_dy.push_back(d_y);
    }

	// start in lane 1;
	int lane = 1;

	// refrence velocity
	double ref_vel = 0.0; // mph

    // runs since last lane change
    int runs_since_last_lane_change = 0;

    // if the car is in the process of slowing down
    bool slowing_down = false;

  h.onMessage([&slowing_down,&runs_since_last_lane_change,&lane,&ref_vel,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx
                      ,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

            // non executed steps size from last run
            int prev_size = previous_path_x.size();

            // set s from last generated position(future) for continuity
            if (prev_size > 0 ) {
                car_s = end_path_s;
            }

            // increment runs since last lane change
            runs_since_last_lane_change++;

            // get lane summary statistics - acts as the behavior planning layer
            tuple<int,double,double> lane_summary_stats = getLaneSummaryStatistics(runs_since_last_lane_change, car_s,
                                                                            prev_size, lane, TARGET_VELOCITY_MPS,
                                                                            slowing_down, sensor_fusion);

            int new_lane = get<0>(lane_summary_stats);
            double distance_diff_meters = get<1>(lane_summary_stats);
            double speed_of_car_ahead_mps = get<2>(lane_summary_stats);

            // safety first - get new velocity
            double new_ref_vel = getNewReferenceVelocityFrom(ref_vel, car_speed, speed_of_car_ahead_mps, distance_diff_meters);
            // are we slowing down?
            if (new_ref_vel < ref_vel)
                slowing_down = true;
            // set new reference velocity
            ref_vel = new_ref_vel;

			// execute lane change if needed and if possible
			executeLaneChangeWithPreconditions(lane, runs_since_last_lane_change, new_lane, car_speed, car_s, prev_size, sensor_fusion);

			// create a list of widely spaced (x,y) waypoints, evenly spaced at 30m
			// later we will interpolate these waypoints with a spline and fill it in with more points that control speed.
			vector<double> ptsx;
			vector<double> ptsy;

			// reference x,y,yaw states
			// either we will reference the starting point as where the car is or the previous paths end point
			double ref_x = car_x;
			double ref_y = car_y;
			double ref_yaw = deg2rad(car_yaw);

			// if previous size is almost empty, use the car as starting point
			if (prev_size < 2) {

				// Use two points that make the path tangent to the car
				double prev_car_x = car_x - cos(car_yaw);
				double prev_car_y = car_y - sin(car_yaw);

				ptsx.push_back(prev_car_x);
				ptsx.push_back(car_x);

				ptsy.push_back(prev_car_y);
				ptsy.push_back(car_y);
			}
				// use the previous path end point as starting reference
			else {
				// redefine reference state as previous end point
				ref_x = previous_path_x[prev_size-1];
				ref_y = previous_path_y[prev_size-1];

				double ref_x_prev = previous_path_x[prev_size-2];
				double ref_y_prev = previous_path_y[prev_size-2];
				ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

				// use two points(unless it's the same point if velocity is 0) that make the path tangent to the previous path end point
				ptsx.push_back(ref_x_prev);
				ptsy.push_back(ref_y_prev);

				if (ref_x_prev != ref_x && ref_y_prev != ref_y) {
					ptsx.push_back(ref_x);
					ptsy.push_back(ref_y);
				}
			}

			// get d for middle of lane
			double middle_of_lane_d = getMiddleOfLaneD(lane);

			// In Frenet add evenly 30m spaced points ahead of the starting reference
			vector<double> next_wp0 = getXY(car_s + TRAJECTORY_S_SPACING, middle_of_lane_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
			vector<double> next_wp1 = getXY(car_s + 2*TRAJECTORY_S_SPACING, middle_of_lane_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
			vector<double> next_wp2 = getXY(car_s + 3*TRAJECTORY_S_SPACING, middle_of_lane_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

			ptsx.push_back(next_wp0[0]);
			ptsx.push_back(next_wp1[0]);
			ptsx.push_back(next_wp2[0]);

			ptsy.push_back(next_wp0[1]);
			ptsy.push_back(next_wp1[1]);
			ptsy.push_back(next_wp2[1]);

			// convert to car perspective
			TransformPointsToVehicleCoordinates(ptsx, ptsy, ref_x, ref_y, ref_yaw);

			// create a spline
			tk::spline spline;

			// set x,y points to the spline
			spline.set_points(ptsx, ptsy);

			// define the actual x,y points we will use for the planner
          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

            // append previous path to current one for continuity
            next_x_vals.insert(std::begin(next_x_vals), std::begin(previous_path_x), std::end(previous_path_x));
            next_y_vals.insert(std::begin(next_y_vals), std::begin(previous_path_y), std::end(previous_path_y));

			// calculate how to break up spline points so that we travel at our desired reference velocity
			double target_x = TRAJECTORY_S_SPACING;
			double target_y = spline(target_x);
			double target_dist = distance(target_x, target_y);
			double x_add_on = 0;

			// calculate number of steps in current velocity to pass target distance
			double N = (target_dist / (CAR_STEP_IN_SECONDS * mph2mps(ref_vel)) );

			// x spacing to generate new points from
			double x_step_spacing = target_x / N;

			// add the next points from polynomial/spline trajectory
			for (int i = 1; i <= NUM_OF_REFERENCE_POINTS - previous_path_x.size(); i++) {
				double x_point = x_add_on + x_step_spacing;
				double y_point = spline(x_point);

				x_add_on = x_point;

				// transform points to global coordinates - rotate and shift
				double x_ref = x_point;
				double y_ref = y_point;

				x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
				y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

				x_point += ref_x;
				y_point += ref_y;

				next_x_vals.push_back(x_point);
				next_y_vals.push_back(y_point);
			}

            json msgJson;
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
