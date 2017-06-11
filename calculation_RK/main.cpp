#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>        
#include <sys/types.h>
#include <time.h>

using namespace std;
//===========================================================================================
//=========================        parameter     ============================================
//===========================================================================================
//--------------------------------------- spatial adjustment
#define pi 3.1415926
#define two_pi 6.2831853
double adjust_angle = 0.25*pi*1.0; //45 degree network
double adjust_cos_angle = cos(adjust_angle); //*.92
double a = 1. / adjust_cos_angle;
double b = 1. / sin(adjust_angle);
//---------------------------------------
double h = 0.005;   //time step 
//---------------------------------------
#define n_print 50   // output data once every 50 time steps
//---------------------------------------
#define tend 50.   // total simulation time
//---------------------------------------
#define Factin_min 500 //F-actin parameter for force clutch
//---------------------------------------
#define a_beta 0.65 // a/kBT with the unit of 1/pN
//---------------------------------------
double k_spring = 0.5 / a; //khmc in the paper
double inverse_gamma = 1.0 *a;//1/gamma in the paper
double Y_branch = 20. *a; // y_br in the paper
double Y_nuc = 10.*a; // y_nuc in the paper
#define Y_force 0.5*a // a small buffer above y_L to count additional subunits for the spring force
//---------------------------------------
#define L0 20 // initial Las17 number 
#define N0 1 // initial filament number
#define L2 101 // N_full in the paper
//-------------------------------- 
int phenotype = 0; // choose 0-WT, 1-LPDA (las17 pan1 Delta acidic)
double LPA_factor = 0.6; //the ratio of kbr(LPDA)/kbr(WT)
//=====================================================================
double l0 = 60.; //l_bar_max in the paper
double knuc = 0.00015;  // k_nuc_max in the paper
double k0, kbr, ksev, kdet; // the 4 fitting parameters in the paper (table 1)
//--------------------------------------------
#define konG 32.            // polymerization rate (unused parameter)
#define kcap 0.667          // capping rate (unused parameter)
//--------------------------------------------
//===========================================================================================
//=========================        variables     ============================================
//===========================================================================================
double L, F, F_br, F_mid;  // L-Las17 #, F-F-actin #, F_br-F_br in the paper, F_mid-F-actin # in the middle exerting pulling force
double lbar; // l_bar(y,t) in the paper
double Y_m, Y_Sla2 = 0.; //Y_m-y_L in the paper, Y_Sla2-y_S in the paper
double Force; // f_out in the paper
double N_att; // n_att in the paper
double BR; // B(t) in the paper (Brownian ratchet factor)
double Lmax, Fmax, tmaxL, tmaxF; // maximum of Las17, F-actin, and the time when Las17 and F-actin occur
int Imax = 0; // maximum of invagination
int current_L_x_10; // current invagination. this also gives which shape function is used currently
double norm = 1.; //normalization factor
//-------------------------------
const int N_y = 2000; // pre-load ram space for y grid
const int N_x = 300; // pre-load ram space for x grid
const int N_L_shape = 100; // total allowed shape function #
double f[N_x][N_y]; // pdf of F-actin (rho in the paper)
double f_y[N_y]; // pdf of F-actin in y direction
double f_x[N_x]; // pdf of F-actin in x direction
double f_Nedelec[N_L_shape]; // f_in in the paper (reference: Serge Dmitrieff, François Nédélec, PLoS Com. Bio., 2015)
bool Shape_Func_exist[N_L_shape]; // if the shape function for given invagination exists or not, 1: yes, 0: no. 
double convert_L = 7.85, convert_f = 565.0; // length and force scales converted into the coordinate of this simulation
int longest_invagination = 30; // if the invagination exceeds this value, an additional term will destruct the F-actin network
//-------------------------------
int ring_center_right = int(20. * b + 0.5); //r_R^in in the paper
int ring_center_left = -ring_center_right; //r_L^in in the paper
int ring_half_width = int(10. * b+0.5);   //width of the Las17 ring
double center = (ring_center_left + ring_center_right)*0.5; //center of the Las17 ring
int ix_start = ring_center_right, ix_end = ring_center_right; // initial position in x (staring and ending points) of the F-actin network
int iy_start, iy_end; // initial position in y (staring and ending points) of the F-actin network
//-------------------------------
//------------------------------------------------------------------------------
const unsigned int N_K = 8; //store the parameter values in a 8 array
double *K = (double*)malloc(sizeof(double)*N_K);
//------------------------------------------------------------------------------ 
const unsigned int N_phi = 8; // discretize rotational dimension into 8 grids
double *phi = (double*)malloc(sizeof(double)*N_phi); // corresponding 8 discrete angles
double *cos_phi = (double*)malloc(sizeof(double)*N_phi); // corresponding 8 discrete cosines
double *R = (double*)malloc(sizeof(double)*N_phi); // corresponding 8 discrete radial lengths (see eq.2 in the paper)
//------------------------------------------------------------------------------ 
const unsigned int num_figure = 10; // output 10 time frames during the simulation
bool plot_series = true; // wether output num_figure frames (true) or just one (false)
bool plot_yet[num_figure]; // wether each of the num_figure frames are printed (true) or not yet (false)
double plot_time[num_figure]; // which frames to print
//=========================================================================================
int n_char; // character variables for loading parameters and out put
char* str1 = new char[40];
char* str2 = new char[40];
char* str3 = new char[40];
//=========================================================================================
const int n_nuc_br_in = 100; // total # of pre-calculated branching, nucleation and forbiden region profiles
const int n_tot_shape = 91; // total # of pre-calculated shape function
int nuc_start[n_tot_shape][n_nuc_br_in]; // starting point of each nucleation region profile
int nuc_end[n_tot_shape][n_nuc_br_in]; // ending point of each nucleation region profile
int br_start[n_tot_shape][n_nuc_br_in]; // starting point of each branching region profile
int br_end[n_tot_shape][n_nuc_br_in]; // ending point of each branching region profile
int in_start[n_tot_shape][n_nuc_br_in]; // starting point of each forbiden region profile
int in_end[n_tot_shape][n_nuc_br_in]; // ending point of each forbiden region profile
int n_br_point[n_tot_shape]; // total spot # of each nucleation region profile
int n_nuc_point[n_tot_shape]; // total spot # of each branching region profile
int n_in_point[n_tot_shape]; // total spot # of each forbiden region profile
//=========================================================================================

//==============================================================================================================
//----------------------------------------------------------------------------------------------- initialization
//==============================================================================================================
void init() {
	int move = 2;
	F = N0*l0, L = L0, F_br = 0.;
	lbar = l0;
	Y_m = N_y - l0*move;
	Y_Sla2 = Y_m;
	for (int iy = 0; iy < N_y; iy++) {
		for (int ix = 0; ix < N_x; ix++) {
			f[ix][iy] = 0;
		}
	}
	//------------------------------------------------------
	for (int iy = N_y - l0*move; iy < N_y - l0*(move - 1); iy++) {
		for (int ix = ix_start; ix <= ix_end; ix++) {
			norm = two_pi * double(ix);
			f[ix][iy] = double(N0) / norm; //***changed
		}
	}
	//------------------------------------------------------*/
	//------------------------------------------------------*/
	iy_start = N_y - l0*move, iy_end = N_y;
	for (int iy = 0; iy < N_y; iy++) {
		f_y[iy] = 0;
		for (int ix = 0; ix < N_x; ix++) {
			norm = two_pi * double(ix);
			f_y[iy] += f[ix][iy] * norm; //***changed
		}
	}
	//------------------------------------------------------*/
	for (int ix = 0; ix < N_x; ix++) {
		f_x[ix] = 0;
		for (int iy = 0; iy < N_y; iy++) {
			norm = two_pi * double(ix);
			f_x[ix] += f[ix][iy] * norm; //***changed
		}
	}
	//------------------------------------------------------*/
	F = 0;
	for (int iy = 0; iy < N_y; iy++) {
		for (int ix = 0; ix < N_x; ix++) {
			norm = two_pi * double(ix);
			F += f[ix][iy] * norm; //***changed
		}
	}
	//------------------------------------------------------*/
	for (int ix = 0; ix < N_L_shape; ix++) {
		Shape_Func_exist[ix] = false;
	}
	Shape_Func_exist[0] = true;
	f_Nedelec[0] = 0.;
	//------------------------------------------------------*/
	int iiiii = int(Y_m+0.5);
	N_att = f_y[iiiii];
	//------------------------------------------------------*/
	for (int i = 0; i < N_phi; i++) {
		phi[i] = double(i) * 0.25 * pi + 0. * pi;
		cos_phi[i] = cos(phi[i]);
	}
	//-------------------------------------------------------
	for (int i = 0; i < num_figure; i++) {
		plot_yet[i] = false;
		plot_time[i] = 8 + i; //14
	}
	//-------------------------------------------------------
	for (int i = 0; i < n_tot_shape; i++) {
		for (int j = 0; j < n_nuc_br_in; j++) {
			nuc_start[i][j] = 0;
			nuc_end[i][j] = 0;
			br_start[i][j] = 0;
			br_end[i][j] = 0;
			in_start[i][j] = 0;
			in_end[i][j] = 0;
		}
		 n_br_point[i] = 0;
		 n_nuc_point[i] = 0;
		 n_in_point[i] = 0;
	}
	//-------------------------------------------------------
	return;
}
//==============================================================================================================
//------------------------------------------------------------------------------------------- mechanical balance
//==============================================================================================================
void force() {
	int n_try = 0;
	double Force_net = 100.;
	double Force_net_pre = Force_net;
	double Y_m_pre = Y_m;
	double Force_diff = 1.;
	ofstream myfile_test;

	bool stop = false;
	double inverse_gamma_curr = inverse_gamma;
	bool oscillation = false, still = false;
	double force_min_threshold = 1.;
	//===============================================================
	while ((fabs(Force_net) > force_min_threshold) && (n_try < 10000)/* && (oscillation == false)*/) {
		n_try++;
		Y_m_pre = Y_m;
		Force_diff = fabs(Force_net - Force_net_pre);
		Force_net_pre = Force_net;
		//------------------------------------------------------------------------
		double current_L = ((Y_Sla2 - Y_m)); //***changed

		current_L_x_10 = 0;
		double current_L_diff_min = 1000;

		for (int L_x_10 = 0; L_x_10 < N_L_shape; L_x_10++) {
			double diff = fabs(current_L - (double(L_x_10)*(convert_L) / 10.));
			if ((Shape_Func_exist[L_x_10] == true) && (diff < current_L_diff_min)) {
				current_L_x_10 = L_x_10;
				current_L_diff_min = diff;
			}
		}
		//------------------------------------------------------------------------
		double force_li;
		double Y_i = Y_Sla2 - Y_m;
		double force_small, force_large, Y_small, Y_large;
		bool minus = false;
		if (Y_i > double(current_L_x_10)*(convert_L) / 10.) {
			force_small = f_Nedelec[current_L_x_10];
			force_large = f_Nedelec[current_L_x_10 + 1];
			Y_small = double(current_L_x_10)*(convert_L) / 10.;
			Y_large = (double(current_L_x_10) + 1.)*(convert_L) / 10.;
		}
		else{
			if (current_L_x_10 - 1 < 0) {
				minus = true;
				force_small = f_Nedelec[0];
			}
			else
				force_small = f_Nedelec[current_L_x_10 - 1];
			force_large = f_Nedelec[current_L_x_10];
			Y_small = (double(current_L_x_10) - 1.)*(convert_L) / 10.;
			Y_large = double(current_L_x_10)*(convert_L) / 10.;
		}
		if (minus == true)
			force_li = force_small;
		else
			force_li = force_small + (force_large - force_small)*(Y_i - Y_small) / (Y_large - Y_small);
		//------------------------------------------------------------------------*/
		Force = 0.;
		for (int iy = 0; iy < N_y; iy++) {
			if ((double(iy) < Y_m + Y_force) && f_y[iy] > 0.)
				Force += k_spring*(double(iy) - Y_m)*double((f_y[iy]));
		}
		//------------------------------------------------------------------------*/
		//Force_net = Force + f_Nedelec[current_L_x_10]; 
		Force_net = Force + force_li;

		if ((Force_net_pre*Force_net < 0.) && (fabs(Force_net - Force_net_pre) - Force_diff< 0.01))
			oscillation = true;

		if (oscillation == true)
			inverse_gamma_curr = 0.01*inverse_gamma;

		Y_m += (inverse_gamma_curr*Force_net)*h;

		//------------------------------------
		if (F < Factin_min)
			Y_Sla2 = Y_m;
		if (Y_Sla2 < Y_m)
			Y_Sla2 = Y_m;

		//------------------------------------
		//Y_Sla2 = Y_m; //***test
	}
	if ((fabs(Force_net) > force_min_threshold)) {
		still = true;
	}
	//=======================================================================
	if (still == true) {
		myfile_test.open("force_test.txt");
		//myfile_test << setprecision(15) << Force_net << "	"  << Y_m << endl;

		double Y_high, Y_low;
		bool up_or_down = false;
		if (Y_m_pre > Y_m) {
			Y_high = Y_m_pre; Y_low = Y_m;
			up_or_down = true;
		}
		else {
			Y_high = Y_m; Y_low = Y_m_pre;
			up_or_down = false;
		}

		//----------------------------------
		cout << "===================================" << endl;
		cout << Y_m_pre << "	" << Y_m << endl;
		cout << Force_net_pre << "	" << Force_net << endl;
		cout << "===================================" << endl;
		//----------------------------------
		double dddY = 0.00001;
		int nY = int((Y_high - Y_low) / dddY);
		double Force_net_min = Force_net, Y_m_best = Y_m, Force_best = Force;
		for (int i_Y = 0; i_Y <= nY; i_Y++) {
			if (fabs(Force_net_min) < 1.)
				break;
			//Y_m = Y_low + double(i_Y)*dddY;
			if (up_or_down == true)
				Y_m += dddY;
			else
				Y_m -= dddY;
			//------------------------------------
			if (F < Factin_min)
				Y_Sla2 = Y_m;
			if (Y_Sla2 < Y_m)
				Y_Sla2 = Y_m;

			//------------------------------------
			//------------------------------------------------------------------------
			double current_L = ((Y_Sla2 - Y_m)); //***changed

			current_L_x_10 = 0;
			double current_L_diff_min = 1000;

			for (int L_x_10 = 0; L_x_10 < N_L_shape; L_x_10++) {
				double diff = fabs(current_L - (double(L_x_10)*(convert_L) / 10.));
				if ((Shape_Func_exist[L_x_10] == true) && (diff < current_L_diff_min)) {
					current_L_x_10 = L_x_10;
					current_L_diff_min = diff;
				}
			}
			//------------------------------------------------------------------------
			double force_li;
			double Y_i = Y_Sla2 - Y_m;
			double force_small, force_large, Y_small, Y_large;
			bool minus = false;
			if (Y_i > double(current_L_x_10)*(convert_L) / 10.) {
				force_small = f_Nedelec[current_L_x_10];
				force_large = f_Nedelec[current_L_x_10 + 1];
				Y_small = double(current_L_x_10)*(convert_L) / 10.;
				Y_large = (double(current_L_x_10) + 1.)*(convert_L) / 10.;
			}
			else{
				if (current_L_x_10 - 1 < 0) {
					minus = true;
					force_small = f_Nedelec[0];
				}
				else
					force_small = f_Nedelec[current_L_x_10 - 1];
				force_large = f_Nedelec[current_L_x_10];
				Y_small = (double(current_L_x_10) - 1.)*(convert_L) / 10.;
				Y_large = double(current_L_x_10)*(convert_L) / 10.;
			}
			if (minus == true)
				force_li = force_small;
			else
				force_li = force_small + (force_large - force_small)*(Y_i - Y_small) / (Y_large - Y_small);
			//------------------------------------------------------------------------*/
			Force = 0.;

			for (int iy = 0; iy < N_y; iy++) {
				if ((double(iy) < Y_m + Y_force) && f_y[iy] > 0.)
					Force += k_spring*(double(iy) - Y_m)*double((f_y[iy]));
			}
			//----------------------------------------------------------------------------
			Force_net = Force + force_li;
			if (fabs(Force_net_min) > fabs(Force_net)) {
				Force_net_min = Force_net;
				Y_m_best = Y_m;
				Force_best = Force;
			}
			myfile_test << Y_m << "	" << Force_net << "	" << force_small << "	" << force_large << "	" << force_li << "	" << Y_i << "	" << Y_small << "	" << Y_large << endl;
		}
		//cout << "stop " << Y_m_best << "	" << Force_net_min << endl;
		myfile_test.close();
		//cin.get();
		Y_m = Y_m_best;
		Force = Force_best;
		if (fabs(Force_net_min) > force_min_threshold) {
			cout << "still too large " << Force_net_min << "	" << F << endl;
			cout << Y_high << "	" << Y_low << endl;
			cin.get();
		}
	}
	//======================================================================= */
	myfile_test.close();

	//===============================================================*/
}
//==============================================================================================================
//------------------------------------------------- calculating filament length for branching (Brownian ratchet)
//==============================================================================================================
void calculate_lbar(int i) {
	if (i - l0 >= Y_m + 0.5)
		lbar = l0;
	else
		lbar = (i - Y_m) + (l0 - (i - Y_m))*(1. / BR);
	return;
}
//==============================================================================================================
//---------------------------------------------------------------------------------------------- output function 
//==============================================================================================================
void output(double t, int i_figure){
	//------------------------------------------------------
	//-------------------------------------------
	n_char = sprintf(str3, "%.6g", float(i_figure));
	strcpy(str1, "data_pdf_");
	strcat(str1, str3);
	strcat(str1, ".txt");
	//------------------------------------------------
	//------------ output rho of the i_figure'th frame 
	//------------------------------------------------
	ofstream myfilePdf;
	myfilePdf.open(str1);
	for (int ix = 0; ix < N_x; ix++) {
		for (int iy = int(Y_m) - 50; iy < int(Y_m) + 200; iy++) {
			myfilePdf << (f[ix][iy]) << "	";
		}
		myfilePdf << endl;
	}
	myfilePdf.close();
	//-----------------------------------------------------
	//-------------------------------------------
	n_char = sprintf(str3, "%.6g", float(i_figure));
	strcpy(str1, "shape_test_");
	strcat(str1, str3);
	strcat(str1, ".txt");
	//------------------------------------------------
	//------------ output rho of the i_figure'th branching, nucleation, forbidden regions, and Las17 ring, and F, L values
	//------------------------------------------------
	ofstream myfiletest;
	myfiletest.open(str1);
	myfiletest << "descriptor x_nuc y_nuc " << endl;
	for (int iy_nuc = 0; iy_nuc < n_nuc_point[current_L_x_10]; iy_nuc++) {
		for (int ix_nuc = nuc_start[current_L_x_10][iy_nuc]; ix_nuc <= nuc_end[current_L_x_10][iy_nuc]; ix_nuc++){
			myfiletest << ix_nuc << "	" << iy_nuc << endl;
		}
	}

	myfiletest << "descriptor x_br y_br " << endl;
	for (int iy_nuc = 0; iy_nuc < n_br_point[current_L_x_10]; iy_nuc++) {
		for (int ix_nuc = br_start[current_L_x_10][iy_nuc]; ix_nuc <= br_end[current_L_x_10][iy_nuc]; ix_nuc++){
			myfiletest << ix_nuc << "	" << iy_nuc << endl;
		}
	}
	myfiletest << "descriptor x_in y_in " << endl;
	for (int iy_nuc = 0; iy_nuc < n_in_point[current_L_x_10]; iy_nuc++) {
		for (int ix_nuc = in_start[current_L_x_10][iy_nuc]; ix_nuc <= in_end[current_L_x_10][iy_nuc]; ix_nuc++){
			myfiletest << ix_nuc << "	" << iy_nuc << endl;
		}
	}
	myfiletest << "descriptor convert" << endl;
	myfiletest << convert_L << endl;
	myfiletest << "descriptor center" << endl;
	myfiletest << center << endl;
	myfiletest << "descriptor x_Las17 y_Las17" << endl;
	for (int i = ring_center_right - ring_half_width; i < ring_center_right + ring_half_width; i++) {
		myfiletest << i << "	" << 0 << endl;
	}
	myfiletest << "descriptor t_rep F_rep L_rep" << endl;
	myfiletest << int(t+0.5) << "	" << F << "	" << L << endl;
	myfiletest.close();
	//--------------------------------------------------------------------
	int digit_1 = current_L_x_10 / 10;
	int digit_decimal = current_L_x_10 - digit_1 * 10;
	//digit_1 = 6, digit_decimal = 5;
	n_char = sprintf(str3, "%.6g", float(digit_1));
	strcpy(str1, "membrane\\shape\\");
	strcat(str1, str3);
	strcat(str1, ".");
	n_char = sprintf(str3, "%.6g", float(digit_decimal));
	strcat(str1, str3);
	strcat(str1, ".txt");
	//--------------------------------------------------------------------
	ifstream infile;
	std::string s;
	infile.open(str1);
	getline(infile, s, '\n');
	getline(infile, s, '	');
	getline(infile, s, '	');
	double R_start = atof(s.c_str());
	infile.close();
	//--------------------------------------------------------------------
	//-------------------------------------------
	n_char = sprintf(str3, "%.6g", float(i_figure));
	strcpy(str2, "shape_");
	strcat(str2, str3);
	strcat(str2, ".txt");
	//------------------------------------------------
	//------------ output rho of the i_figure'th shape
	//------------------------------------------------
	ofstream myfileAdd;
	myfileAdd.open(str2);
	myfileAdd << "descriptor s_test R_test Z_test" << endl;
	myfileAdd << 0 << "	" << double(ring_center_right + ring_half_width)/convert_L << "	" << 0 << endl;
	infile.open(str1);
	getline(infile, s, '\n');
	while (getline(infile, s, '	')) {
		double s_Nedelec = atof(s.c_str());
		getline(infile, s, '	');
		double R_Nedelec = atof(s.c_str());
		getline(infile, s, '	');
		double Z_Nedelec = atof(s.c_str());
		getline(infile, s, '\n');
		myfileAdd << s_Nedelec << "	" << R_Nedelec << "	" << Z_Nedelec << endl;
	}
	infile.close();
	myfileAdd.close();
	//--------------------------------------------------------------------
	return;
}
//==============================================================================================================
//------------------------------------------------------------------------------------------------------- clutch 
//==============================================================================================================
double alpha(double x, double x_0, double A){
	double xxx = exp(A * (x - x_0));
	return (((xxx - 1. / xxx) / (xxx + 1. / xxx)) + 1)*0.5;
	//return 1.;
}
//==============================================================================================================
//---------------------------------------------------------------------------------------------------- main loop
//==============================================================================================================

int main() {
	clock_t t_cpu;
	t_cpu = clock();
	//------------------------------------------------------------
	ifstream infile;
	std::string s;
	//------------------------------------------------------------
	//--------------------------- reading fitting parameter values
	//------------------------------------------------------------
	infile.open("..//parameter.txt");
	getline(infile, s, '	'); getline(infile, s, '\n');
	k0 = atof(s.c_str());
	getline(infile, s, '	'); getline(infile, s, '\n');
	kbr = atof(s.c_str());
	getline(infile, s, '	'); getline(infile, s, '\n');
	ksev = atof(s.c_str());
	getline(infile, s, '	'); getline(infile, s, '\n');
	kdet = atof(s.c_str());
	getline(infile, s, '	'); getline(infile, s, '\n');
	infile.close();
	//------------------------------------------------------------
	// ------------------------------- output file for time course
	//------------------------------------------------------------
	ofstream myfile;
	myfile.open("data_timeCourse.txt");
	myfile << "descriptor time F L n_att ivagination(nm) y_L f_out" << endl;
	double t = 0., t0 = 0.;
	int nt = (tend / h);

	//================================================================================================== 
	if (phenotype == 0) {
		K[0] = k0, K[1] = kbr, K[2] = ksev, K[3] = kdet, K[4] = konG, K[5] = kcap, K[6] = knuc, K[7] = 0;
	}
	else if (phenotype == 1) {
		K[0] = k0, K[1] = kbr*LPA_factor, K[2] = ksev, K[3] = kdet, K[4] = konG, K[5] = kcap, K[6] = knuc, K[7] = 0;
	}

	//==================================================================================================
	Lmax = 0, Fmax = 0, tmaxL = 0, tmaxF = 0;
	init(); t = 0.;
	//==================================================================================================
	//=========================================================== reading shape functions and parameters
	//==================================================================================================
	for (int digit_1 = 0; digit_1 < 10; digit_1++){
		int decimal_tot = 10;
		for (int digit_decimal = 0; digit_decimal < decimal_tot; digit_decimal++){
			int L_x_10 = digit_1 * 10 + digit_decimal;
			if (L_x_10 >= N_L_shape) {
				cout << "Too many shapes!" << endl;
				cin.get();
			}
			n_char = sprintf(str3, "%.6g", float(digit_1));
			strcpy(str1, "membrane\\parameter\\");
			strcat(str1, str3);
			strcat(str1, ".");
			n_char = sprintf(str3, "%.6g", float(digit_decimal));
			strcat(str1, str3);
			strcat(str1, ".txt");
			infile.open(str1);
			if (getline(infile, s, '	'))
				Shape_Func_exist[L_x_10] = true;
			infile.close();
			//-------------------------------------------------
			if (Shape_Func_exist[L_x_10] == true) {
				//-------------------------------------------
				infile.open(str1);
				getline(infile, s, '	');
				getline(infile, s, '\n');
				double Ri_Nedelec = atof(s.c_str());
				getline(infile, s, '	');
				getline(infile, s, '\n');
				double L_Nedelec = atof(s.c_str());
				getline(infile, s, '\n');
				getline(infile, s, '\n');
				getline(infile, s, '\n');
				getline(infile, s, '	');
				getline(infile, s, '\n');
				f_Nedelec[L_x_10] = atof(s.c_str()) * convert_f;
				infile.close();
		    } // for if (Shape_Func_exist[L_x_10] == true)
			//-------------------------------------------------
	    } //for digit_decimal
    } //for digit_1
	//==================================================================================================
	//================================================= reading branching, nucleation, forbidden regions
	//==================================================================================================
	for (int digit_1 = 0; digit_1 < 10; digit_1++){
		int decimal_tot = 10;
		for (int digit_decimal = 0; digit_decimal < decimal_tot; digit_decimal++){
			int n_shape = 10 * digit_1 + digit_decimal;
			//-----------------------------------------------------
			n_char = sprintf(str3, "%.6g", float(digit_1));
			strcpy(str1, "membrane\\nuc_br_in\\br_");
			if (digit_1 > 0)
			strcat(str1, str3);
			n_char = sprintf(str3, "%.6g", float(digit_decimal));
			strcat(str1, str3);
			strcat(str1, ".txt");
			infile.open(str1);
			while (getline(infile, s, '	')) {
			getline(infile, s, '	');
			br_start[n_shape][n_br_point[n_shape]] = int(atof(s.c_str()));
			getline(infile, s, '\n');
			br_end[n_shape][n_br_point[n_shape]] = int(atof(s.c_str()));
			n_br_point[n_shape] ++;
			}
			infile.close();
			//-----------------------------------------------------
			n_char = sprintf(str3, "%.6g", float(digit_1));
			strcpy(str1, "membrane\\nuc_br_in\\nuc_");
			if (digit_1 > 0)
				strcat(str1, str3);
			n_char = sprintf(str3, "%.6g", float(digit_decimal));
			strcat(str1, str3);
			strcat(str1, ".txt");
			infile.open(str1);
			while (getline(infile, s, '	')) {
				getline(infile, s, '	');
				nuc_start[n_shape][n_nuc_point[n_shape]] = int(atof(s.c_str()));
				getline(infile, s, '\n');
				nuc_end[n_shape][n_nuc_point[n_shape]] = int(atof(s.c_str()));
				n_nuc_point[n_shape] ++;
			}
			infile.close();
			//-----------------------------------------------------
			n_char = sprintf(str3, "%.6g", float(digit_1));
			strcpy(str1, "membrane\\nuc_br_in\\in_");
			if (digit_1 > 0)
				strcat(str1, str3);
			n_char = sprintf(str3, "%.6g", float(digit_decimal));
			strcat(str1, str3);
			strcat(str1, ".txt");
			infile.open(str1);
			while (getline(infile, s, '	')) {
				getline(infile, s, '	');
				in_start[n_shape][n_in_point[n_shape]] = int(atof(s.c_str()));
				getline(infile, s, '\n');
				in_end[n_shape][n_in_point[n_shape]] = int(atof(s.c_str()));
				n_in_point[n_shape] ++;
			}
			infile.close();
			//-----------------------------------------------------
		} //for digit_decimal
	} //for digit_1
    //==================================================================================================
				//======================================================================
				int i_print = 0;
				int number_of_lines = 0;
				bool pdf_output = false;
				//==================================================================================================

				//==================================================================================================
				for (int it = 0; it < nt; it++) {  //time loop
			    //==================================================================================================
					//-----------------------------------------------------
					i_print++;
					//-----------------------------------------------------
					if (current_L_x_10 == 0) {
						double clutch_force_factor = alpha(F, Factin_min, 0.005);
						f_Nedelec[0] = clutch_force_factor*f_Nedelec[1];//0.1
					}
					//-----------------------------------------------------*/
					//-----------------------------------------------------
					if (current_L_x_10 > longest_invagination) {
						K[7] = 1.0;
						//cout << "invagination is " << current_L_x_10 << ", too long!" << endl; cin.get();
					}
					if (current_L_x_10 > Imax)
						Imax = current_L_x_10;
					//-----------------------------------------------------
					double Force_pre = Force;
					force();
					int Y_m_int = int(Y_m + 0.5);

					//N_att = f_y[Y_m_int];
					N_att = (f_y[Y_m_int] + f_y[Y_m_int + 1] + f_y[Y_m_int - 1]) / 3.;
					BR = exp(-Force / N_att*a_beta);
					//------------------------------------------------
					if (F / BR < 1000.)
						h = 0.005;
					else
						h = 0.0002;
					//t = t0 + double(it) * h;
					t += h;
					//------------------------------------------------
					//--------------------------------- Las17 dynamics
					//------------------------------------------------
					L += (K[0] * L *L*(L2 - L) - K[3] * L*F_br - K[7] * L)*h;
					//==================================================================
					//=================== nucleation ===================================
					//==================================================================
					for (int iy_nuc = 0; iy_nuc < n_nuc_point[current_L_x_10]; iy_nuc++) {
						for (int ix_nuc = nuc_start[current_L_x_10][iy_nuc]; ix_nuc <= nuc_end[current_L_x_10][iy_nuc]; ix_nuc++){
							int ix = ix_nuc, iy = Y_m_int + iy_nuc;
							//----------------------------------------				
							double xxx = knuc * 80 * (two_pi* double(ix));
							//-------------------------------------------------------
							double increment;
							double sigma = (ring_center_right - ring_half_width)*1.0;
							double gaussian_factor = exp(-double(ix)*double(ix) / 2. / sigma / sigma);
							increment = gaussian_factor * xxx * h;
							//-----------------------------------------------------
							norm = two_pi * double(ix);
							f[ix][iy] += increment / norm;
							f_y[iy] += increment; //***changed;
							f_x[ix] += increment; //***changed;
							//----------------------------------------
						}
					}
					//==================================================================
					//=================== branching ====================================
					//==================================================================					
					F_br = 0.;
					for (int iy_br = 0; iy_br < n_br_point[current_L_x_10]; iy_br++) {
						for (int ix_br = br_start[current_L_x_10][iy_br]; ix_br <= br_end[current_L_x_10][iy_br]; ix_br++){
							int ix = ix_br, iy = Y_m_int + iy_br;
							//-------------------------------------------------------------
							double sigma = Y_branch*0.5;
							double xxx = K[1] * exp(-double(iy_br)*double(iy_br) / 2. / sigma / sigma) * double((f[ix][iy] * two_pi* double(ix)));
							//----------------------------*/
							//-------------------------------------------------------
							F_br += xxx;
							xxx *= (h*L);
							//--------------------------------------
							calculate_lbar(iy);
							//-------------------------------------------------------
							double increment = xxx / double(N_phi); //***changed;
							for (int kkk = 0; kkk <= 1 * (int(lbar) - 1); kkk++){
								int j = iy - kkk;
								double Dy = double(kkk), r = double(ix);
								//-----------------------------------------------------
								if ((j < N_y) && (j >= 0)) {
									//-----------------------------------------------------
									for (int i_phi = 0; i_phi < N_phi; i_phi++) {
										R[i_phi] = sqrt(r*r + Dy*Dy - 2.*r*Dy*cos_phi[i_phi]);
										int r_in, r_out;
										int round_R = int(R[i_phi] + 0.5);
										if (round_R < R[i_phi]){
											r_in = round_R, r_out = round_R + 1;
										}
										else{
											r_in = round_R - 1, r_out = round_R;
										}
										//------------------------------------------------
										int j_br = j - Y_m_int;
										if (j_br < 0)
											j_br = 0;
										if (r_in > in_end[current_L_x_10][j_br]) {
										if ((r_in > 0) && (r_in < N_x) ) {
											double ratio = (r_out - R[i_phi]);
											norm = two_pi * double(r_in);
											f[r_in][j] += ratio*increment / norm;
											f_y[j] += ratio*increment; //***changed;
											f_x[r_in] += ratio*increment; //***changed;
										}
										}
										if (r_out > in_end[current_L_x_10][j_br]) {
											if ((r_out > 0) && (r_out < N_x)) {
												double ratio = (R[i_phi] - r_in);
												norm = two_pi * double(r_out);
												f[r_out][j] += ratio*increment / norm;
												f_y[j] += ratio*increment; //***changed;
												f_x[r_out] += ratio*increment; //***changed;
											}
										}
										//------------------------------------------------
									} //i_phi
									//-----------------------------------------------------
								}
							}
							//----------------------
							//--------------------------------------------------------
						}//iy_nuc
					}//ix_nuc		
					//==================================================================
					//=================== severing =====================================
					//==================================================================
					for (int iy = 0; iy < N_y; iy++) {
						for (int ix = 0; ix < N_x; ix++) {
							double xxx = f[ix][iy] * K[2] * h;
							f[ix][iy] -= xxx;
							norm = two_pi * double(ix);
							f_y[iy] -= xxx* norm; //***changed
							f_x[ix] -= xxx* norm; //***changed
						}
					}
					//==============================================================================
					//===================================================================
					
					//=====================================================================
					F = 0.;
					for (int iy = 0; iy < N_y; iy++) {
						for (int ix = 0; ix < N_x; ix++) {
							norm = two_pi * double(ix);
							F += f[ix][iy] * norm; //***changed
						}
					}
					//=====================================================================
					if (L > Lmax) {
						tmaxL = t, Lmax = L;
					}
					if (F > Fmax) {
						tmaxF = t, Fmax = F;
					}
					//----------------------------------------------------
					//----------------------------------------------------
					//--------------------------------------------------------------
					if ((i_print > n_print)){
						cout << t << "	" << Y_Sla2 << "	" << Y_m << "	" << current_L_x_10 << endl;
						//------------------------------------------------
						//-------------------------- output maximum values
						//------------------------------------------------
						ofstream myfile3;
						myfile3.open("max.txt");
						myfile3 << "descriptor tFmax Fmax tLmax Lmax Imax" << endl;
						myfile3 << tmaxF << "	" << Fmax << "	" << tmaxL << "	" << Lmax << "	" << double(Imax)*(convert_L) / 10. / a*2.7 << endl;
						myfile3.close();
						//---------------------------------------------						
						myfile << t << "	" << F << "	" << L << "	" << N_att << "	" << double(current_L_x_10)*(convert_L) / 10. / a*2.7 << "	" << (Y_m - N_y) / a << "	" << Force << endl;
						i_print -= n_print;

					}
					//--------------------------------------------------------------
					if (plot_series == true){
						for (int i_figure = 0; i_figure < num_figure; i_figure++) {
							if ((plot_yet[i_figure] == false) && (fabs(plot_time[i_figure] - t) < 0.004)){
								plot_yet[i_figure] = true;
								output(t, i_figure);
							}
						}
				    }
					//--------------------------------------------------------------
				} //it

	return 0;
}