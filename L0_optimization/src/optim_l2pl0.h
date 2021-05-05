


vec optim_withCplex_model_bigM_l2pl0(Node *node, mat A ,vec y, double lambda, double bigM, uvec Q1, uvec Q0){

	vec  Ay;
		Ay = A.t() * y;
	//	bigM = 1*max(abs(Ay)) / pow(norm(A.col(1)), 2);
    int Q = A.n_cols;
    vec x_relache = zeros(Q);
    vec b_relache = zeros(Q);
	clock_t t0, t1;


		IloEnv env;
			try {

				 t0=clock();

					IloModel model(env);
					IloCplex cplex(model);

		//-----------variables------------

					IloNumVarArray x(env, Q, -bigM, bigM); //declaration de Q variables réelles de borne inf -bigM et borne supp bigM
					IloNumVarArray b(env, Q,0,1); //** declaration de Q variables Booléennes

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

					// lambda sum b
					for (int i = 0; i < Q; i++)
						objExpr += lambda * b[i];

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
					}

					for (int i = 0; i < Q0.n_rows; i++) {
					model.add(IloRange(env, 0, b[Q0(i)], 0));
					}

					for (int i = 0; i < Q1.n_rows; i++) {
										model.add(IloRange(env, 1, b[Q1(i)], 1));
										}

					//-------------Initialization------------
					/*			IloNumArray startx(env);
					for (int i = 0; i < Q; i++) {
							startx.add(node->x_relache(i));
					}

					cplex.add(x, startx,IloCplex::MIPStartAuto,"secondMIPStart");

					IloNumArray startb(env);
					for (int i = 0; i < Q; i++) {
							startb.add(node->b_relache(i));
					}

					cplex.addMIPStart(b, startb,IloCplex::MIPStartAuto,"secondMIPStart");
					startb.end();

					startx.end();
*/

					cplex.setOut(env.getNullStream());
					//cplex.setParam(IloCplex::Threads, 1); // nombre de threads
					cplex.setParam(IloCplex::Param::RootAlgorithm, 1);


					t1=clock();
					T_model +=(float)(t1-t0)/CLOCKS_PER_SEC;

					cplex.solve();

					IloNumArray vals_x(env);
					IloNumArray vals_b(env);

					cplex.getValues(vals_x,x);
					cplex.getValues(vals_b,b);


					for (int i = 0; i < Q; i++)
						x_relache[i]=vals_x[i];

					for (int i = 0; i < Q; i++)
						b_relache[i]=vals_b[i];

		//			Q1.print("Q1");
		//			Q0.print("Q0");

		//			x_relache.t().print("[optim] x_relache");
		//			b_relache.t().print("[optim] b_relache");

		//			cout << max(abs(x_relache)) << endl;
		//			cout << cplex.getObjValue() << endl;

		//			cout << pow(norm(y - A * x_relache, 2), 2) +  lambda * sum(b_relache) - pow(norm(y, 2), 2) << endl;

		//			pause_cin();

					env.end();

				} catch (IloException & e) {
					cerr << " ERREUR : exception = " << e << endl;
				}

node->b_relache = b_relache;

return x_relache;
}

void generate_model_with_bigM_l2pl0(mat H, vec y,  double lambda, double bigM){
	mat  HH;
	vec  Hy;

	HH = H.t() * H;
		Hy = H.t() * y;
	//	bigM = 1*max(abs(Hy)) / pow(norm(H.col(1)), 2);
    int Q = H.n_cols;
		//-----------Affichage------------
		cout << "bigM = " << bigM << endl;
		cout << "lambda = " << lambda << endl;

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

					// lambda sum b
					for (int i = 0; i < Q; i++)
					objExpr += lambda * b[i];

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


					IloNumArray startx(env);
					for (int i = 0; i < Q; i++) {
							startx.add(0);
					}

					cplex.addMIPStart(x, startx,IloCplex::MIPStartAuto,"secondMIPStart");


					//------------ save Model------------
					cplex.exportModel("model_L2pL0_BigM.lp");
					std::cout<<"modèle avec BigM créé"<<endl;

					env.end();

				} catch (IloException & e) {
					cerr << " ERREUR : exception = " << e << endl;
				}

}

void generate_model_with_sos_l2pl0(mat H, vec y,  double lambda, double bigM){
	mat  HH;
	vec  Hy;

	HH = H.t() * H;
		Hy = H.t() * y;
	//	bigM = 1*max(abs(Hy)) / pow(norm(H.col(1)), 2);
    int Q = H.n_cols;
		//-----------Affichage------------
		cout << "bigM = " << bigM << endl;
		cout << "lambda = " << lambda << endl;

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

					// lambda sum b
					for (int i = 0; i < Q; i++)
					objExpr += lambda * b[i];

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


					IloNumArray startx(env);
					for (int i = 0; i < Q; i++) {
							startx.add(0);
					}

					cplex.addMIPStart(x, startx,IloCplex::MIPStartAuto,"secondMIPStart");


					//------------ save Model------------
					cplex.exportModel("model_L2pL0_BigM.lp");
					std::cout<<"modèle avec BigM créé"<<endl;

					env.end();

				} catch (IloException & e) {
					cerr << " ERREUR : exception = " << e << endl;
				}

}

vec Cplexsolver_l2pl0(mat *A,vec *y, double lambda,double BigM, bool sos){

    int Q = A->n_cols;
    vec x = zeros(Q);


IloEnv env;
	try {
			IloModel model(env);
			IloCplex cplex(model);


		clock_t t0, t1;
		t0=clock();

if(sos)
			generate_model_with_sos_l2pl0(*A, *y, lambda, BigM);
else
			generate_model_with_bigM_l2pl0(*A, *y, lambda, BigM);

		t1=clock();
		T_model +=(float)(t1-t0)/CLOCKS_PER_SEC;


		IloObjective obj;
		IloNumVarArray var(env);
		IloRangeArray rng(env);
		cplex.importModel(model, "model_L2pL0_BigM.lp", obj, var, rng);
		//------------Resolution------------


		//cplex.setParam(IloCplex::Threads, 4); // nombre de threads
		cplex.setParam(IloCplex::Param::TimeLimit,TimeBBMax);
		cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap,  gap);
		if(branching_rule==4){
		cplex.setParam(IloCplex::VarSel, CPX_VARSEL_STRONG);
		cplex.setParam(IloCplex::BrDir, CPX_BRDIR_UP);
		}

		cplex.solve();




		std::cout << std::endl << "------------Resultat--------------"<< std::endl;
		CplexNodenum =cplex.getNnodes() ;
		CplexStatus =cplex.getStatus();
		born_sup = cplex.getObjValue();
		//std::cout << "Time = " << cplex.getTime() - start << " sec"<< std::endl;

		IloNumArray vals(env);
		cplex.getValues(vals, var);


		for (int i = 0; i < Q; i++){
				x(i)=vals[i];

		}
		env.end();

	} catch (IloException & e) {
		cerr << " ERREUR : exception = " << e << endl;
	}

	return x;
}
