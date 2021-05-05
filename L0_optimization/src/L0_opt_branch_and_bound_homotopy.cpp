//============================================================================
// Name        : BranchandBound_l2l0M.cpp
// Author      : ramzi
// Version     :
// Copyright   : Your copyright notice
// Description : BranchandBound algorithm in C++, Ansi-style
//============================================================================
#include <iostream>
#include <vector>
#include <math.h>

#include <time.h>

#include <sstream>
#include <string>
#include <ilcplex/ilocplex.h>

#include "armadillo"
//#include "model.h"
ILOSTLBEGIN

using namespace arma;

#include "global_var.h"
#include "Util.h"

#include "Homotopy_opt_l2l0.h"
#include "Homotopy_opt_l0l2.h"
#include "Homotopy_opt_l2pl0.h"

#include "Active_set.h"

#include "optim_l2l0.h"
#include "optim_l0l2.h"
#include "optim_l2pl0.h"

#include "branch_and_bound_l2l0.h"
#include "branch_and_bound_l0l2.h"
#include "branch_and_bound_l2pl0.h"


#include <sys/ioctl.h>
#include <fcntl.h>
#include <linux/kd.h>

/************************    Principal Program    ****************************************/

int main(int argc, char** argv) {

	clock_t t_, t2_;
	T_homo = 0;
	T_test = 0;
	T_model = 0;

	vec epsilon;
	vec lambda;

	stringstream pathdata,pathdata_H, pathdata_y, pathdata_eps, pathdata_lambda, pathdata_out, pathdata_out_sol, pathdata_out_time;



// default version = 8 ;
string	version ="";
string	regle ="";
bool opt_sos=false;
bool paral = 0;

bool L2L0, L2pL0, L0L2;
	L2pL0 = 0;
	L2L0 = 0;
	L0L2 = 0;

int SNR = 3;
int k= 5;
float rho = 0.8;
int N = 100;
int Q=1000;



cout <<endl<< "Emplacement données : "<<endl;
string path = argv[1];
cout << path << endl;

cout <<endl<< "Type de probleme (1: L2L0 ; 2: L2pL0 ; 3: L0L2 ): ";
int type_p = atoi(argv[2]);
cout << type_p << endl;

cout << endl << "TimeBBMax = ";
TimeBBMax = atoi(argv[3]);
cout << TimeBBMax << endl;


if(type_p ==1)
		L2L0 = 1;
else if(type_p ==2)
		L2pL0 = 1;
else
		L0L2 = 1;


opt_with_cplex = false  ;

opt_relaxation_with_cplex = 0;
opt_relaxation_with_hom = true;
opt_relaxation_with_actset = 0;


//cout <<"version cplex : (0: no paral, 1: paral)"<<endl;
paral = 0;
		if(paral)
			version +="_paral";
//		cout <<"opt sos : (0: no , 1: yes)"<<endl;
opt_sos = 0;
			if(opt_sos)
				version +="_sos_bigM";



if (L2L0) {
	if (opt_with_cplex)
		cout<< "L2L0_MIP_Cplex"<<version;
	else if(opt_relaxation_with_cplex)
		cout<< "L2L0_BB_Rcplex";
	else if(opt_relaxation_with_hom)
		cout<< "L2L0_BB_Rhom";
	else
		cout<< "L2L0_BB_Ract";
}

if (L2pL0) {
	if (opt_with_cplex)
		cout<< "L2pL0_MIP_Cplex"<<version;
	else if(opt_relaxation_with_cplex)
		cout<< "L2pL0_BB_Rcplex";
	else if(opt_relaxation_with_hom)
		cout<< "L2pL0_BB_Rhom";
	else
			cout<< "L2pL0_BB_Ract";

}

if (L0L2) {
	if (opt_with_cplex)
		cout<< "L0L2_MIP_Cplex"<<version;
	else if(opt_relaxation_with_cplex)
		cout<< "L0L2_BB_Rcplex";
	else if(opt_relaxation_with_hom)
		cout<< "L0L2_BB_Rhom";
}


branching_rule =1;

if(branching_rule){
regle = "";
cout << regle << endl;
}

		pathdata.str("");
		pathdata_H.str("");
		pathdata_y.str("");
		pathdata_eps.str("");
		pathdata_lambda.str("");

		//pathdata << "/home/rbenmhenni/Bureau/these/Matlab/test_random_dic_gloutonvsMIP/genere_donnees/Cpp_multimedia_material_bertsimas_rand/N"<<N<<"Q"<<Q<<"SA_SNR"<< SNR << "_rho" << rho << "_K" << k << "_instance" << instance << "/";
		pathdata << path;

		cout << pathdata.str()<<endl;



		pathdata_H << pathdata.str() << "H.dat";
		pathdata_y << pathdata.str() << "y.dat";



		pathdata_eps << pathdata.str() << "alpha_bruit.dat";
		pathdata_lambda << pathdata.str() << "lambda.dat";


		A.load(pathdata_H.str(), raw_ascii);
		eyeN.eye(A.n_rows, A.n_rows);
		AA = A.t() * A;


		y.load(pathdata_y.str(), raw_ascii);
		norm2carre_y = pow(norm((y) , 2), 2);


		born_sup = 10000;
		T_model = 0;
		T_homo = 0;


		BigM = 1.1 * max(abs(A.t() * y)) / pow(norm(A.col(1)), 2);
		vec x_sol, x0;
		x0 =zeros<vec> (A.n_cols);


		pathdata_out.str("");
		pathdata_out_sol.str("");
		pathdata_out_time.str("");



		try{

		if (L2L0) {
			k= atoi(argv[4]);
			cout << endl <<" k=" << k << endl;
			t_ = clock();
			if (opt_with_cplex){
				x_sol = Cplexsolver_l2l0(&A, &y, k, BigM, opt_sos);
				pathdata_out << pathdata.str() << "L2L0_MIP_Cplex_T" << TimeBBMax <<""<< regle  << version<<"_res.txt";
				pathdata_out_sol << pathdata.str() << "Sol_L2L0_MIP_Cplex_T" << TimeBBMax <<""<< regle << version << ".dat";
				pathdata_out_time << pathdata.str() << "Time_L2L0_MIP_Cplex_T" << TimeBBMax <<""<< regle << version << ".dat";

			}else{
				x_sol = BranchandBound_l2l0(&A, &y, k, BigM);

				if(opt_relaxation_with_cplex){
					pathdata_out << pathdata.str() << "L2L0_BB_Rcplex_T" << TimeBBMax <<""<<  regle <<"_res.txt";
					pathdata_out_sol << pathdata.str() << "Sol_L2L0_BB_Rcplex_T" << TimeBBMax <<  regle <<".dat";
					pathdata_out_time << pathdata.str() << "Time_L2L0_BB_Rcplex_T" << TimeBBMax <<""<< regle  << ".dat";
				}else{
					pathdata_out << pathdata.str() << "L2L0_BB_Rhom_T" << TimeBBMax <<""<< regle <<"_res.txt";
					pathdata_out_sol << pathdata.str() << "Sol_L2L0_BB_Rhom_T" << TimeBBMax <<""<< regle <<".dat";
					pathdata_out_time << pathdata.str() << "Time_L2L0_BB_Rhom_T" << TimeBBMax <<""<< regle  << ".dat";
				}

				}
			t2_ = clock();

		} else if (L0L2) {
			epsilon.load(pathdata_eps.str(), raw_ascii);

			cout << "epsilon = " << epsilon(0) << endl;
			cut = 1;



			t_ = clock();
			if (opt_with_cplex){

				x_sol = Cplexsolver_l0l2(&A, &y, epsilon(0), BigM);
				pathdata_out << pathdata.str() << "L0L2_MIP_Cplex_T" << TimeBBMax<<""<< regle << version<<"_res.txt";
				pathdata_out_sol << pathdata.str() << "Sol_L0L2_MIP_Cplex_T" << TimeBBMax<<""<< regle << version<<".dat";
				pathdata_out_time << pathdata.str() << "Time_L0L2_MIP_Cplex_T" << TimeBBMax<<""<< regle << version<<".dat";

			}else{

				x_sol = BranchandBound_l0l2(&A, &y, epsilon(0), BigM);

				if(opt_relaxation_with_cplex){
					 pathdata_out << pathdata.str() << "L0L2_BB_Rcplex_T" << TimeBBMax<<"_res.txt";
					 pathdata_out_sol << pathdata.str() << "Sol_L0L2_BB_Rcplex_T" << TimeBBMax<<""<< regle <<".dat";
					 pathdata_out_time << pathdata.str() << "Time_L0L2_BB_Rcplex_T" << TimeBBMax<<""<< regle <<".dat";
				}else {
				     pathdata_out << pathdata.str() << "L0L2_BB_Rhom_T" << TimeBBMax<<""<< regle <<"_res.txt";
				     pathdata_out_sol << pathdata.str() << "Sol_L0L2_BB_Rhom_T" << TimeBBMax<<""<< regle <<".dat";
				     pathdata_out_time << pathdata.str() << "Time_L0L2_BB_Rhom_T" << TimeBBMax<<""<< regle <<".dat";
				}
			}

			t2_ = clock();
		}else if (L2pL0) {
			lambda.load(pathdata_lambda.str(), raw_ascii);

			t_ = clock();
			if (opt_with_cplex){

				x_sol = Cplexsolver_l2pl0(&A, &y, lambda(0), BigM, opt_sos);

				pathdata_out << pathdata.str() << "L2pL0_MIP_Cplex_T" <<  TimeBBMax<< regle << version<<"_res.txt";
				pathdata_out_sol << pathdata.str() << "Sol_L2pL0_MIP_Cplex_T" << TimeBBMax <<""<< regle<< version <<".dat";
				pathdata_out_time << pathdata.str() << "Time_L2pL0_MIP_Cplex_T" << TimeBBMax <<""<< regle<< version<<".dat";

			}else{
				x_sol = BranchandBound_l2pl0(&A, &y, lambda(0), BigM, x0);

				if(opt_relaxation_with_cplex){
					pathdata_out << pathdata.str() << "L2pL0_BB_Rcplex_T" << TimeBBMax<<""<< regle <<"_res.txt";
					pathdata_out_sol << pathdata.str() << "Sol_L2pL0_BB_Rcplex_T" << TimeBBMax<<""<< regle <<".dat";
					pathdata_out_time << pathdata.str() << "Time_L2pL0_BB_Rcplex_T" << TimeBBMax<<""<< regle <<".dat";
				}else if(opt_relaxation_with_hom){
					pathdata_out << pathdata.str() << "L2pL0_BB_Rhom_T" << TimeBBMax<<""<< regle <<"_res.txt";
					pathdata_out_sol << pathdata.str() << "Sol_L2pL0_BB_Rhom" << TimeBBMax <<""<< regle <<".dat";
					pathdata_out_time << pathdata.str() << "Time_L2pL0_BB_Rhom" << TimeBBMax <<""<< regle <<".dat";


				}else if(opt_relaxation_with_actset){
					pathdata_out << pathdata.str() << "L2pL0_BB_Ract_T" << TimeBBMax<<""<< regle <<"_nowarmstart_res.txt";
					pathdata_out_sol << pathdata.str() << "Sol_L2pL0_BB_Ract_T"  << TimeBBMax <<""<< regle <<".dat";
					pathdata_out_time << pathdata.str() << "Time_L2pL0_BB_Ract_T" << TimeBBMax <<""<< regle <<".dat";

				}

				}
			t2_ = clock();
		}


		}catch(exception e ){
			cout << "execption"<<endl;
			//continue;
		}

		//         RESULT


		ofstream fichier(pathdata_out.str().c_str(), ios::out);
		ofstream fichier_time(pathdata_out_time.str().c_str(), ios::out);

		if (fichier) {
			cout << endl << "BigM:" << BigM << endl;
			fichier  << "BigM: " << BigM << endl;

			if (opt_with_cplex) {

				cout << "Node Number Cplex: " << CplexNodenum << endl;
				cout << "val opt: " << born_sup << endl;
				cout << "temps d'execution: " << (float) (t2_ - t_)
						/ CLOCKS_PER_SEC << endl;

				fichier << "Node_Number_Cplex: " << CplexNodenum << endl;
				fichier << "val_opt: " << born_sup << endl;
				fichier << "temps_d'execution: " << (float) (t2_ - t_)
						/ CLOCKS_PER_SEC << endl;
				fichier << "CplexStatus: " << CplexStatus << endl;

				//	cout << "temps test: : " <<T_test<<endl;

			} else {

				cout << "Node_Number_BB: " << BBNodenum << endl;
				cout << "val_opt: " << born_sup << endl;
				cout << "temps_d'execution: " << (float) (t2_ - t_)
						/ CLOCKS_PER_SEC - T_model << endl;
				cout << "temps de modelisation : " << T_model << endl;

				cout << "temps homotopy: " << T_homo << endl;
				//	cout << "temps test: : " <<T_test<<endl;
				//	cout << "itmax= " <<itmax<<endl;
				cout << "Best_Node_num: "<< BestNodenum << endl;
				cout << "Best_temps_num: "<<T_best_Bsupp << endl;

				fichier << "Node_Number_BB: " << BBNodenum << endl;
				fichier << "val_opt: " << born_sup << endl;
				fichier << "temps_d'execution: " << (float) (t2_ - t_)
						/ CLOCKS_PER_SEC - T_model << endl;
				fichier << "temps_de_modelisation: " << T_model << endl;
				fichier << "temps_homotopy: " << T_homo << endl;
				fichier << "Best_Node_num: "<< BestNodenum << endl;
				fichier << "Best_temps_num: "<<T_best_Bsupp << endl;

			}
			fichier.close();


		} else
			cerr << "Impossible d'ouvrir le fichier !" << endl;

		if (fichier_time) {
			fichier_time << (float) (t2_ - t_)
									/ CLOCKS_PER_SEC - T_model << endl;

			fichier_time.close();

		} else
					cerr << "Impossible d'ouvrir le fichier time !" << endl;


cout << pathdata_out_sol.str()<<endl;
x_sol.t().print("x_sol");
x_sol.save(pathdata_out_sol.str(),raw_ascii);



}
