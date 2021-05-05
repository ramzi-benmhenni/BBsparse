using namespace arma;


void generate_modelwith_bigM_logic(mat H, vec y, int k, double bigM){
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

						model.add( IloIfThen(env, b[i]==0, IloAbs(x[i]) == 0  ) );
						// or	model.add(IloRange  (env, 0, bigM * b[i]- IloAbs(x[i]), IloInfinity,"ctr_bigM"));
					}

					//------------ save Model------------
					cplex.exportModel("model_BigM.lp");
					std::cout<<"modèle avec BigM créé"<<endl;

					env.end();

				} catch (IloException & e) {
					cerr << " ERREUR : exception = " << e << endl;
				}

}

void generate_modell2_l0with_bigM(mat H, vec y, int k, double bigM){
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
					cplex.exportModel("model_BigM.lp");
					std::cout<<"modèle avec BigM créé"<<endl;

					env.end();

				} catch (IloException & e) {
					cerr << " ERREUR : exception = " << e << endl;
				}

}

void generate_modelwithout_bigM(mat H, vec y, int k){
	mat  HH;
	vec  Hy;

	HH = H.t() * H;
		Hy = H.t() * y;
	//	bigM = 1*max(abs(Hy)) / pow(norm(H.col(1)), 2);
    int Q = H.n_cols;
		//-----------Affichage------------

		cout << "K = " << k << endl;


try {
				IloEnv env;
				IloModel model(env);
				IloCplex cplex(model);


				std::cout<<"-----------variables------------"<<endl;

					IloNumVarArray x(env, Q,-IloInfinity,IloInfinity); //declaration de Q variables réelles
					IloBoolVarArray b(env, Q); //** declaration de Q variables Booléennes

					//** Facultatif:: nommer les variables **//
					for (int i = 0; i < Q; i++) {
						std::stringstream name;
						name << "b" << i + 1;
						b[i] = IloBoolVar(env, name.str().c_str());
					}

					model.add(x);
					model.add(b);

					std::cout<<"-----------Objective------------"<<endl;


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

					std::cout<<"-----------Contraintes------------"<<endl;

					IloRange ctr_sum(env, 0, IloSum(b), k, "ctr_sum");
					model.add(ctr_sum);

					for (int i = 0; i < Q; i++) {
						std::stringstream name;
						name << "ctr_bigM_" << i + 1 << "_1";

						model.add( IloIfThen(env, b[i]==0, IloAbs(x[i]) == 0  ) );


						// or	model.add(IloRange  (env, 0, bigM * b[i]- IloAbs(x[i]), IloInfinity,"ctr_bigM"));
					}

					//------------ save Model------------
					cplex.exportModel("model_sans_BigM.lp");
					std::cout<<"modèle sans BigM créé"<<endl;

					env.end();

				} catch (IloException & e) {
					std::cout<< " ERREUR : exception = " << e << endl;
				}

}


void generate_modell2plusl0with_bigM(mat H, vec y, double lambda, double bigM){
	mat  HH;
	vec  Hy;

	HH = H.t() * H;
		Hy = H.t() * y;
	//	bigM = 1*max(abs(Hy)) / pow(norm(H.col(1)), 2);
    int Q = H.n_cols;

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

					// -labda b
					for (int i = 0; i < Q; i++)
						objExpr += lambda * b[i];

					// X'*HH*X
					for (int i = 0; i < Q; i++)
						for (int j = 0; j < Q; j++)
							objExpr += x[i] * HH(i, j) * x[j];

					IloObjective obj(env, objExpr, IloObjective::Minimize, "OBJ");
					model.add(obj);

					//-----------Contraintes------------


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
					cplex.exportModel("modell2plusl0_BigM.lp");
					std::cout<<"modèle avec BigM créé"<<endl;

					env.end();

				} catch (IloException & e) {
					cerr << " ERREUR : exception = " << e << endl;
				}

}




void generate_modell2plusl0with_bigM_bcarre(mat H, vec y, double lambda, double bigM){
	mat  HH;
	vec  Hy;

	HH = H.t() * H;
		Hy = H.t() * y;
	//	bigM = 1*max(abs(Hy)) / pow(norm(H.col(1)), 2);
    int Q = H.n_cols;

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

					// -labda b
					for (int i = 0; i < Q; i++)
						objExpr += lambda * b[i]* b[i];

					// X'*HH*X
					for (int i = 0; i < Q; i++)
						for (int j = 0; j < Q; j++)
							objExpr += x[i] * HH(i, j) * x[j];

					IloObjective obj(env, objExpr, IloObjective::Minimize, "OBJ");
					model.add(obj);

					//-----------Contraintes------------


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
					cplex.exportModel("modell2plusl0_BigM_bcarre.lp");
					std::cout<<"modèle avec BigM créé"<<endl;

					env.end();

				} catch (IloException & e) {
					cerr << " ERREUR : exception = " << e << endl;
				}

}
