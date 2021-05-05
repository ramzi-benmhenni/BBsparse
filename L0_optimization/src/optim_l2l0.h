
int CplexNodenum = 0;
int CplexStatus = 0;

double T_model;

vec optim_QP_bigM(mat H, vec y, double bigM){
	clock_t t0, t1;
	t0=clock();

	mat  HH;
	vec  Hy;

	HH = H.t() * H;
	Hy = H.t() * y;

    int Q = H.n_cols;
    //cout << "Q= "<<Q<<endl;

	vec w = zeros(Q);
	IloEnv env;

			try {


					IloModel model(env);
					IloCplex cplex(model);

		//-----------variables------------
					IloNumVarArray x(env, Q, -bigM, bigM); //declaration de Q variables réelles de borne inf -bigM et borne supp bigM

					model.add(x);
					//-----------Objective------------

					IloExpr objExpr(env);

					// -2 Y' * H * X
					for (int i = 0; i < Q; i++)
						objExpr += -2 * Hy(i) * x[i];

					// X'*HH*X
					for (int i = 0; i < Q; i++)
						for (int j = 0; j < Q; j++)
							objExpr += x[i] * HH(i, j) * x[j];

					IloObjective obj(env, objExpr, IloObjective::Minimize, "OBJ");
					model.add(obj);

					cplex.setParam(IloCplex::Threads, 1); // nombre de threads
					cplex.setOut(env.getNullStream());





t1=clock();
												 T_model +=(float)(t1-t0)/CLOCKS_PER_SEC;

					cplex.solve();

					IloNumArray vals(env);
					cplex.getValues(vals,x);



					for (int i = 0; i < Q; i++)
							w[i]=vals[i];

					env.end();

					//w.print("w");

				} catch (IloException & e) {
					cerr << " ERREUR : exception = " << e << endl;
				}

	return w;

}
vec optim_withCplex_model_bigM_l2l0( mat A ,vec y, int k, double bigM, uvec Q1, uvec Q0){

	vec  Ay;
		Ay = A.t() * y;
	//	bigM = 1*max(abs(Ay)) / pow(norm(A.col(1)), 2);
    int Q = A.n_cols;
    vec x_relache = zeros(Q);
	clock_t t0, t1;


		IloEnv env;
			try {

				 t0=clock();

					IloModel model(env);
					IloCplex cplex(model);

		//-----------variables------------

					IloNumVarArray x(env, Q, -bigM, bigM); //declaration de Q variables réelles de borne inf -bigM et borne supp bigM
					IloNumVarArray b(env, Q,0,1); //** declaration de Q variables Booléennes

					//** Facultatif:: nommer les variables **//
					/*for (int i = 0; i < Q; i++) {
						std::stringstream name;
						name << "b" << i + 1;
						b[i] = IloNumVar(env,0,1, name.str().c_str());
					}
*/
					model.add(x);
					model.add(b);

					//-----------Objective------------

					IloExpr objExpr(env);

					// -2 Y' * A * X
					for (int i = 0; i < Q; i++)
						objExpr += -2 * Ay(i) * x[i];

					// X'*AA*X
					for (int i = 0; i < Q; i++)
						for (int j = 0; j < Q; j++)
							objExpr += x[i] * AA(i, j) * x[j];

					IloObjective obj(env, objExpr, IloObjective::Minimize, "OBJ");
					model.add(obj);

					//-----------Contraintes------------

					IloRange ctr_sum(env, 0, IloSum(b), k, "ctr_sum");
					model.add(ctr_sum);

					for (int i = 0; i < Q; i++) {
						std::stringstream name;
						name << "ctr_bigM_" << i + 1 << "_1";
						model.add(
								IloRange(env, -IloInfinity, -bigM * b[i] - x[i], 0,
										name.str().c_str()));

						name << "ctr_bigM_" << i + 1 << "_2";
						model.add(
								IloRange(env, 0, bigM * b[i] - x[i], IloInfinity,
										name.str().c_str()));
						// or	model.add(IloRange  (env, 0, bigM * b[i]- IloAbs(x[i]), IloInfinity,"ctr_bigM"));
					}

					for (int i = 0; i < Q0.n_rows; i++) {
					model.add(IloRange(env, 0, b[Q0(i)], 0));
					}

					for (int i = 0; i < Q1.n_rows; i++) {
										model.add(IloRange(env, 1, b[Q1(i)], 1));
										}

					cplex.setOut(env.getNullStream());
					cplex.setParam(IloCplex::Threads, 1); // nombre de threads



					 t1=clock();
							 T_model +=(float)(t1-t0)/CLOCKS_PER_SEC;



										cplex.solve();

					IloNumArray vals(env);
					cplex.getValues(vals,x);


					for (int i = 0; i < Q; i++)
						x_relache[i]=vals[i];

					env.end();

				} catch (IloException & e) {
					cerr << " ERREUR : exception = " << e << endl;
				}


return x_relache;
}


 void generate_model_with_sos_l2l0(mat H, vec y, int k, double bigM){
 	mat  HH;
 	vec  Hy;

 	HH = H.t() * H;
 		Hy = H.t() * y;
 	//	bigM = 1*max(abs(Hy)) / pow(norm(H.col(1)), 2);
     int Q = H.n_cols;
 		//-----------Affichage------------

 		cout << "K = " << k << endl;
 		cout << "SOS1 " << endl;

 		IloEnv env;
 			try {

 					IloModel model(env);
 					IloCplex cplex(model);


 		//-----------variables------------

 					IloNumVarArray x(env, Q, -bigM, bigM); //declaration de Q variables réelles de borne inf -bigM et borne supp bigM
 					IloBoolVarArray b(env, Q); //** declaration de Q variables Booléennes
 					IloNumVarArray s(env, Q, 0, 1); //

 					//** Facultatif:: nommer les variables **//
 					for (int i = 0; i < Q; i++) {
 						std::stringstream name;
 						name << "b" << i + 1;
 						b[i] = IloBoolVar(env, name.str().c_str());
 					}

 					model.add(x);
 					model.add(b);

 					//-----------Objective------------

 					IloExpr objExpr(env);

 					// -2 Y' * H * X
 					for (int i = 0; i < Q; i++)
 						objExpr += -2 * Hy(i) * x[i];

 					// X'*HH*X
 					for (int i = 0; i < Q; i++)
 						for (int j = 0; j < Q; j++)
 							objExpr += x[i] * HH(i, j) * x[j];

 					IloObjective obj(env, objExpr, IloObjective::Minimize, "OBJ");
 					model.add(obj);

 					//-----------Contraintes------------

 					IloRange ctr_sum(env, 0, IloSum(b), k, "ctr_sum");
 					model.add(ctr_sum);


 					for (int i = 0; i < Q; i++) {
 						//std::stringstream name;
 						//name << "ctr_bigM_" << i + 1 << "_1";
 						//model.add(
 						//		IloRange(env, -IloInfinity, -bigM * b[i] - x[i], 0,
 						//				name.str().c_str()));

 						//name << "ctr_bigM_" << i + 1 << "_2";
 						//model.add(
 						//		IloRange(env, 0, bigM * b[i] - x[i], IloInfinity,
 						//				name.str().c_str()));
 						std::stringstream name_eg;
 						name_eg << "ctr_egality_" << i + 1 ;
 						model.add( IloRange(env, 0, s[i]- 1 + b[i], 0, name_eg.str().c_str()) );
 						IloNumVarArray ar(env, 2, x[i], s[i]);
 						std::stringstream name_sos;
 						name_sos << "sos_" << i + 1;
 						IloSOS1 sos(env, ar, name_sos.str().c_str());
 						model.add(sos);

 						// or	model.add(IloRange  (env, 0, bigM * b[i]- IloAbs(x[i]), IloInfinity,"ctr_bigM"));
 					}

 					//------------ save Model------------
 					cplex.exportModel("model_L2L0_BigM.lp");
 					std::cout<<"modèle avec sos créé"<<endl;

 					env.end();

 				} catch (IloException & e) {
 					cerr << " ERREUR : exception = " << e << endl;
 				}

 }

void generate_model_with_bigM_l2l0(mat H, vec y, int k, double bigM){
	mat  HH;
	vec  Hy;

	HH = H.t() * H;
		Hy = H.t() * y;
	//	bigM = 1*max(abs(Hy)) / pow(norm(H.col(1)), 2);
    int Q = H.n_cols;
		//-----------Affichage------------
		cout << "bigM = " << bigM << endl;
		cout << "K = " << k << endl;

		IloEnv env;
			try {

					IloModel model(env);
					IloCplex cplex(model);


		//-----------variables------------

					IloNumVarArray x(env, Q, -bigM, bigM); //declaration de Q variables réelles de borne inf -bigM et borne supp bigM
					IloBoolVarArray b(env, Q); //** declaration de Q variables Booléennes

					//** Facultatif:: nommer les variables **//
					for (int i = 0; i < Q; i++) {
						std::stringstream name;
						name << "b" << i + 1;
						b[i] = IloBoolVar(env, name.str().c_str());
					}

					model.add(x);
					model.add(b);

					//-----------Objective------------

					IloExpr objExpr(env);

					// -2 Y' * H * X
					for (int i = 0; i < Q; i++)
						objExpr += -2 * Hy(i) * x[i];

					// X'*HH*X
					for (int i = 0; i < Q; i++)
						for (int j = 0; j < Q; j++)
							objExpr += x[i] * HH(i, j) * x[j];

					IloObjective obj(env, objExpr, IloObjective::Minimize, "OBJ");
					model.add(obj);

					//-----------Contraintes------------

					IloRange ctr_sum(env, 0, IloSum(b), k, "ctr_sum");
					model.add(ctr_sum);

					for (int i = 0; i < Q; i++) {
						std::stringstream name;
						name << "ctr_bigM_" << i + 1 << "_1";
						model.add(
								IloRange(env, -IloInfinity, -bigM * b[i] - x[i], 0,
										name.str().c_str()));

						name << "ctr_bigM_" << i + 1 << "_2";
						model.add(
								IloRange(env, 0, bigM * b[i] - x[i], IloInfinity,
										name.str().c_str()));
						// or	model.add(IloRange  (env, 0, bigM * b[i]- IloAbs(x[i]), IloInfinity,"ctr_bigM"));
					}

					//------------ save Model------------
					cplex.exportModel("model_L2L0_BigM.lp");
					std::cout<<"modèle avec BigM créé"<<endl;

					env.end();

				} catch (IloException & e) {
					cerr << " ERREUR : exception = " << e << endl;
				}

}

void generate_model_with_bigM_l2l0_NC(mat H, vec y, int k, double bigM){
	mat  HH;
	vec  Hy;

	HH = H.t() * H;
		Hy = H.t() * y;
	//	bigM = 1*max(abs(Hy)) / pow(norm(H.col(1)), 2);
    int Q = H.n_cols;
    int N = H.n_rows;
		//-----------Affichage------------
		cout << "bigM = " << bigM << endl;
		cout << "K = " << k << endl;

		IloEnv env;
			try {

					IloModel model(env);
					IloCplex cplex(model);


		//-----------variables------------

					IloNumVarArray x(env, Q, -bigM, bigM); //declaration de Q variables réelles de borne inf -bigM et borne supp bigM
					IloBoolVarArray b(env, Q); //** declaration de Q variables Booléennes
					IloNumVarArray z(env, N, -IloInfinity, IloInfinity );

					cout << "N"<< N <<endl;

					model.add(x);
					model.add(b);
					model.add(z);


					cout << "2"<<endl;
					//-----------Objective------------

					IloExpr objExpr(env);

					// -2 Y' * Z
					for (int i = 0; i < N; i++)
						objExpr += -2  * z[i] * y(i);

					cout << "3"<<endl;

					// Z'*Z
					for (int i = 0; i < N; i++)
							objExpr += z[i] * z[i];

					cout << "4"<<endl;

					IloObjective obj(env, objExpr, IloObjective::Minimize, "OBJ");
					model.add(obj);

					cout << "5"<<endl;

					//-----------Contraintes------------

					for (int i = 0; i < N; i++){
						IloExpr exp_Hx_i(env);
						for (int j = 0; j < Q; j++) {
							exp_Hx_i += H(i,j) * x[j];
						//cout << i << " " << j;
									}
						//cout <<endl;
					model.add(IloRange(env, 0, exp_Hx_i - z[i], 0 ));
					}


					IloRange ctr_sum(env, 0, IloSum(b), k, "ctr_sum");
					model.add(ctr_sum);

					for (int i = 0; i < Q; i++) {
						std::stringstream name;
						name << "ctr_bigM_" << i + 1 << "_1";
						model.add(
								IloRange(env, -IloInfinity, -bigM * b[i] - x[i], 0,
										name.str().c_str()));

						name << "ctr_bigM_" << i + 1 << "_2";
						model.add(
								IloRange(env, 0, bigM * b[i] - x[i], IloInfinity,
										name.str().c_str()));
						// or	model.add(IloRange  (env, 0, bigM * b[i]- IloAbs(x[i]), IloInfinity,"ctr_bigM"));
					}

					//------------ save Model------------
					cplex.exportModel("model_L2L0_BigM.lp");
					std::cout<<"modèle avec BigM créé"<<endl;

					env.end();

				} catch (IloException & e) {
					cerr << " ERREUR : exception = " << e << endl;
				}

}

void generate_model_with_sos_l2l0_NC(mat H, vec y, int k, double bigM){
	mat  HH;
	vec  Hy;

	HH = H.t() * H;
		Hy = H.t() * y;
	//	bigM = 1*max(abs(Hy)) / pow(norm(H.col(1)), 2);
    int Q = H.n_cols;
    int N = H.n_rows;
		//-----------Affichage------------
		cout << "K = " << k << endl;

		IloEnv env;
			try {

					IloModel model(env);
					IloCplex cplex(model);


		//-----------variables------------

					IloNumVarArray x(env, Q, -bigM, bigM); //declaration de Q variables réelles de borne inf -bigM et borne supp bigM
					IloBoolVarArray b(env, Q); //** declaration de Q variables Booléennes
					IloNumVarArray z(env, N, -IloInfinity, IloInfinity );
					IloNumVarArray s(env, Q, 0, 1); //

					cout << "N"<< N <<endl;

					model.add(x);
					model.add(b);
					model.add(z);



					//-----------Objective------------

					IloExpr objExpr(env);

					// -2 Y' * Z
					for (int i = 0; i < N; i++)
						objExpr += -2  * z[i] * y(i);



					// Z'*Z
					for (int i = 0; i < N; i++)
							objExpr += z[i] * z[i];



					IloObjective obj(env, objExpr, IloObjective::Minimize, "OBJ");
					model.add(obj);



					//-----------Contraintes------------

					for (int i = 0; i < N; i++){
						IloExpr exp_Hx_i(env);
						for (int j = 0; j < Q; j++) {
							exp_Hx_i += H(i,j) * x[j];
						//cout << i << " " << j;
									}
						//cout <<endl;
					model.add(IloRange(env, 0, exp_Hx_i - z[i], 0 ));
					}


					IloRange ctr_sum(env, 0, IloSum(b), k, "ctr_sum");
					model.add(ctr_sum);

					for (int i = 0; i < Q; i++) {
					/*	std::stringstream name;
						name << "ctr_bigM_" << i + 1 << "_1";
						model.add(
								IloRange(env, -IloInfinity, -bigM * b[i] - x[i], 0,
										name.str().c_str()));

						name << "ctr_bigM_" << i + 1 << "_2";
						model.add(
								IloRange(env, 0, bigM * b[i] - x[i], IloInfinity,
										name.str().c_str()));
						// or	model.add(IloRange  (env, 0, bigM * b[i]- IloAbs(x[i]), IloInfinity,"ctr_bigM"));
						*/
 						std::stringstream name_eg;
 						name_eg << "ctr_egality_" << i + 1 ;
 						model.add( IloRange(env, 0, s[i]- 1 + b[i], 0, name_eg.str().c_str()) );
 						IloNumVarArray ar(env, 2, x[i], s[i]);
 						std::stringstream name_sos;
 						name_sos << "sos_" << i + 1;
 						IloSOS1 sos(env, ar, name_sos.str().c_str());
 						model.add(sos);

					}

					//------------ save Model------------
					cplex.exportModel("model_L2L0_BigM.lp");
					std::cout<<"modèle avec sos créé"<<endl;

					env.end();

				} catch (IloException & e) {
					cerr << " ERREUR : exception = " << e << endl;
				}

}

vec Cplexsolver_l2l0(mat *A,vec *y,int k ,double BigM, bool sos){

    int Q = A->n_cols;
    int N = A->n_rows;
    vec x_sol = zeros(Q);


IloEnv env;
	try {
			IloModel model(env);
			IloCplex cplex(model);

			clock_t t0, t1;
				t0=clock();

		//if ( 0) {
		//	generate_modell2plusl0with_bigM(H, y, lambda, bigM);
		//	generate_modell2plusl0with_bigM_bcarre(H, y, lambda, bigM);
			if(Q<N){
				cout << "Q<N" << endl;
				if(sos)
				generate_model_with_sos_l2l0(*A, *y, k, BigM);
				else
				generate_model_with_bigM_l2l0(*A, *y, k, BigM);
			}else{
				cout << "Q>=N" << endl;
				if(sos)
				generate_model_with_sos_l2l0_NC(*A, *y, k, BigM);
				else
				generate_model_with_bigM_l2l0_NC(*A, *y, k, BigM);
			}


		t1=clock();
			T_model +=(float)(t1-t0)/CLOCKS_PER_SEC;
		//}else{

		IloObjective obj;
		IloNumVarArray var(env);
		IloRangeArray rng(env);
		cplex.importModel(model, "model_L2L0_BigM.lp", obj, var, rng);
		//------------Resolution------------

		cplex.setParam(IloCplex::Param::TimeLimit,TimeBBMax);
		//cplex.setParam(IloCplex::Param::Threads,1);
		cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap,  gap);
		if(branching_rule==4){
		cplex.setParam(IloCplex::VarSel, CPX_VARSEL_STRONG);
		cplex.setParam(IloCplex::BrDir, CPX_BRDIR_UP);
		}

		cplex.solve();

		//------------Resultat--------------
		std::cout << std::endl << "------------Resultat--------------"
				<< std::endl;
		CplexNodenum =cplex.getNnodes() ;
		CplexStatus =cplex.getStatus();
		born_sup = cplex.getObjValue();

		IloNumArray vals(env);
		cplex.getValues(vals, var);

	//	cout << "size(var) " << var.getSize() << endl;


	if(Q<N){

		for (int i = 0; i < Q; i++)
				x_sol[i]=vals[i];
	}else{
		for (int i = 0; i < Q; i++)
				x_sol[i]=vals[i+N];
	}



	//	}
		env.end();

		} catch (IloException & e) {
		cerr << " ERREUR : exception = " << e << endl;
	}

	return x_sol;
}

