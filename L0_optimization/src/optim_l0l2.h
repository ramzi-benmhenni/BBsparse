

vec optim_withCplex_model_bigM_l0l2( mat A, vec y, double epsilon, double bigM, uvec Q1, uvec Q0){

	vec  Ay;
		Ay = A.t() * y;
		vec yy= y.t()*y;

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

					IloExpr cnstExpr(env);

					// -2 Y' * A * X
					for (int i = 0; i < Q; i++)
						cnstExpr += -2 * Ay(i) * x[i];

					// X'*AA*X
					for (int i = 0; i < Q; i++)
						for (int j = 0; j < Q; j++)
							cnstExpr += x[i] * AA(i, j) * x[j];


					//-----------Contraintes------------

				//	IloRange ctr_sum(env, 0,cnstExpr, epsilon  , "ctr_Qp");
					IloRangeArray ctr(env);

					ctr.add(cnstExpr +yy(0) <= pow(epsilon,2) );
					model.add(ctr);


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

					for (unsigned i = 0; i < Q0.n_rows; i++) {
					model.add(IloRange(env, 0, b[Q0(i)], 0));
					}

					for (unsigned i = 0; i < Q1.n_rows; i++) {
										model.add(IloRange(env, 1, b[Q1(i)], 1));
										}



					IloObjective obj(env, IloSum(b) , IloObjective::Minimize, "OBJ");
					model.add(obj);


					cplex.setOut(env.getNullStream());
					//cplex.setParam(IloCplex::Threads, 1); // nombre de threads



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

void generate_model_l0l2_with_bigM(mat H, vec y, double epsilon, double bigM){
	mat  HH;
	vec  Hy;
	vec yy;

	HH = H.t() * H;
		Hy = H.t() * y;

	yy= y.t()*y;

	//	bigM = 1*max(abs(Hy)) / pow(norm(H.col(1)), 2);
    int Q = H.n_cols;
		//-----------Affichage------------
		cout << "bigM = " << bigM << endl;
		cout << "epsilon = " << epsilon << endl;

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

					IloObjective obj(env, IloSum(b) , IloObjective::Minimize, "OBJ");
					model.add(obj);

					std::cout<<"obj ajouté"<<endl;



					//-----------Contraintes------------

					IloExpr cnstExpr(env);

					// X'*AA*X
					for (int i = 0; i < Q; i++)
						for (int j = 0; j < Q; j++)
							cnstExpr += HH(i, j) * x[i] * x[j];

				//	cnstExpr += HH(0, 0) * x[0] * x[0];

					// -2 Y' * A * X
					for (int i = 0; i < Q; i++)
						cnstExpr += -2 * Hy(i) * x[i];


				//	IloRange ctr_Qp(env, 0, cnstExpr, epsilon , "ctr_Qp");


					IloRangeArray ctr(env);
					ctr.add(cnstExpr +yy(0) <= pow(epsilon,2) );
					//ctr.add(HH(0, 0) * x[0] * x[0] <= 0);



					std::cout<<"Contraintes créés"<<endl;
					model.add(ctr);

					std::cout<<"Contraintes ajoutés"<<endl;



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



				//	cplex.setOut(env.getNullStream());

					//------------ save Model------------
					cplex.exportModel("model_L0L2_BigM.lp");
					std::cout<<"modèle avec BigM créé"<<endl;

					cnstExpr.end();
					env.end();

				} catch (IloException & e) {
					cerr << " ERREUR : exception = " << e << endl;
					exit(0);
				}

}



vec Cplexsolver_l0l2(mat *A,vec *y,double epsilon ,double BigM){

    int Q = A->n_cols;
    vec x = zeros(Q);


IloEnv env;
	try {
			IloModel model(env);
			IloCplex cplex(model);


			clock_t t0, t1;
				t0=clock();

		generate_model_l0l2_with_bigM(*A, *y, epsilon, BigM);

		t1=clock();
				T_model +=(float)(t1-t0)/CLOCKS_PER_SEC;

		IloObjective obj;
		IloNumVarArray var(env);
		IloRangeArray rng(env);
		cplex.importModel(model, "model_L0L2_BigM.lp", obj, var, rng);
		//------------Resolution------------



		cplex.setParam(IloCplex::Param::TimeLimit,TimeBBMax);
		cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap,  gap);
		if(branching_rule==4){
		cplex.setParam(IloCplex::VarSel, CPX_VARSEL_STRONG);
		cplex.setParam(IloCplex::BrDir, CPX_BRDIR_UP);
		}
		cplex.solve();




		std::cout << std::endl << "------------Resultat--------------"<< std::endl;
		CplexNodenum =cplex.getNnodes() ;
		CplexStatus =cplex.getStatus() ;



		born_sup = cplex.getObjValue();
		//std::cout << "Time = " << cplex.getTime() - start << " sec"<< std::endl;

		IloNumArray vals(env);
		cplex.getValues(vals, var);


		for (int i = 0; i < Q; i++){
				x(i)=vals[i+Q];

		}
		env.end();

	} catch (IloException & e) {
		cerr << " ERREUR : exception = " << e << endl;
	}

	return x;
}
