
bool is_binary(vec b) {

	for (int i = 0; i < b.n_elem; i++) {
		if (b(i) > eps_binary && 1 - b(i) > eps_binary)
			return false;
	}
	return true;

}

/************************** regle de branchement ***********************************/
int branching_rules(Node *node, vec x_relache) {

	int qj = -1;
	//max_xj
	uvec Qv = conv_to<uvec>::from(node->vec_Qv);
	uvec Q1 = conv_to<uvec>::from(node->vec_Q1);
	uvec Q0 = conv_to<uvec>::from(node->vec_Q0);

	try {

		//if(branching_rule==1)
		qj = Qv(index_max(abs(x_relache(Qv))));

		if(branching_rule==2)
		qj = Qv(index_max(abs(node->score(Qv))));

		else if(branching_rule==3){
			vec score;
			score =abs( 0.5*ones(size(node->b_relache)) - node->b_relache);

		//	node->b_relache.t().print("b");
		//	score.t().print("score");

			//int a;
		//	cin>>a;
			qj = Qv(index_min(score(Qv)));

		}


		if (!binary_search(node->vec_Qv.begin(), node->vec_Qv.end(), qj))
			cout << "attention qj not in Qv";

	} catch (std::logic_error e) {
		x_relache(Qv).print();
		Qv.print("Qv");
		Q1.print("Q1");
		Q0.print("Q0");


		node->b_relache.t().print("b");

		cout << "message from branchin_rules"<< endl;
		pause_cin();

	}

	return qj;

}

/****************************    relaxation     ****************************************/
vec relaxation_l2l0(Node *node, mat *A, vec *y, int k, double BigM) {
	clock_t th, thf;
	uvec Q1, Q0, Qv;
	vec x_Qv, x, w, b,score;
	mat A_Q1, A_Qv;

	x = zeros(A->n_cols);
	score = x;
	b = zeros(A->n_cols);

	double Tau = (k - node->vec_Q1.size()) * BigM;

	Q1 = conv_to<uvec>::from(node->vec_Q1);
	Q0 = conv_to<uvec>::from(node->vec_Q0);
	Qv = conv_to<uvec>::from(node->vec_Qv);

	b(Q1)= ones(Q1.n_elem);



	if (opt_relaxation_with_cplex) {
		x = optim_withCplex_model_bigM_l2l0(*A, *y, k, BigM, Q1, Q0);

		if(Qv.n_elem)
		b(Qv)= abs(x(Qv))/BigM ;

		node->x_relache = x ;
		node->b_relache = b ;

		node->score = abs(x) ;
		return x;
	}

	A_Q1 = A->cols(Q1);
	A_Qv = A->cols(Qv);

	if (node->branch == 0) {
		node->invAAQ1 = node->invAAQ1_parent;
	} else {
		node->invAAQ1 = inv(AA.submat(Q1, Q1));//version recursive!!
	}

	w = node->invAAQ1 * (A_Q1.t() * *y);

	if (!w.is_empty())
		if (BigM < max(abs(w))) {
			//cout<<"BigM < max(abs(w)) apple a optim_QP"<<endl;
			//w.print("w_old");
			w = optim_QP_bigM(A_Q1, *y, BigM);


		}

	if (node->vec_Q1.size() < k) {

		vec score_tmp;
		th = clock();
		x_Qv = mixed_BigM_homotopytemps_l2l0(A_Q1, &w, &A_Qv, Q1, Qv, *y, BigM,
				Tau,&score_tmp);

		thf = clock();
		T_homo += (float) (thf - th) / CLOCKS_PER_SEC;
		//x_Qv=homotopy(&BB,&Bz,Tau);
		x(Qv) = x_Qv;
		score(Qv) = score_tmp;
	}
	x(Q1) = w;



	if(Qv.n_elem)
	b(Qv)= abs(x(Qv))/BigM ;

	node->x_relache = x ;
	node->b_relache = b ;
	node->score = score ;

	return x;
}

vec relaxation_bigM_l2l0(Node *node, mat *A, vec *y, int k, double BigM) {
	clock_t th, thf;
	uvec Q1, Q0, Qv;
	vec x_Qv, x, b, Bz;
	mat P, BB, s;
	mat A_Q1, A_Qv;

	x = zeros<vec> (A->n_cols);
	b = zeros(A->n_cols);

	double Tau = (k - node->vec_Q1.size()) * BigM;

	Q1 = conv_to<uvec>::from(node->vec_Q1);
	Q0 = conv_to<uvec>::from(node->vec_Q0);
	Qv = conv_to<uvec>::from(node->vec_Qv);


	b(Q1)= ones(Q1.n_elem);

	A_Q1 = A->cols(Q1);
	A_Qv = A->cols(Qv);

	/*	if(node->branch==0){
	 node->invAAQ1=node->invAAQ1_parent;
	 BB=node->BB;
	 Bz=node->Bz;

	 s=node->invAAQ1*A_Q1.t();

	 }else{*/

	node->invAAQ1 = inv(AA.submat(Q1, Q1));//version recursive!!

	s = node->invAAQ1 * A_Q1.t();
	P = -A_Q1 * s + eyeN;
	//AA=A_Qv.t()*P;

	BB = A_Qv.t() * P * A_Qv;

	Bz = A_Qv.t() * (P * (*y));

	th=clock();
	//	 node->BB=BB;
	//	 node->Bz=Bz;
	//}


	// x_Qv=homotopy_onecol(&BB,&Bz,Tau);
	  x_Qv=homotopy_onecol_l2l0(&BB,&Bz,Tau);

	 thf=clock();

	x(Qv) = x_Qv;
	x(Q1) = s * ((*y) - A_Qv * x_Qv);

	T_homo +=(float)(thf-th)/CLOCKS_PER_SEC;


	b(Qv)= abs(x(Qv))/BigM ;

	node->x_relache = x ;
	node->b_relache = b ;

	return x;
}

/****************************    Node evaluate   ***************************************/
vec evaluation_l2l0(Node *node, mat *A, vec *y, int k, double BigM) {

	vec x_relache = zeros(A->n_cols);
	/*if(node->vec_Q1.size()==k){
	 node->is_solution=true;
	 uvec Q1= conv_to <uvec>::from(node->vec_Q1);

	 // RQ: l'invertion à mettre à jour récurcivement
	 //x_relache(Q1) = solve( AA.submat(Q1,Q1),  A->cols(Q1).t() * (*y) );
	 x_relache(Q1) = inv( AA.submat(Q1,Q1)) * ( A->cols(Q1).t() * (*y) );


	 node->born_inf= pow(norm((*y)- (*A) * x_relache,2),2);
	 return x_relache;
	 }else{*/

	x_relache = relaxation_l2l0(node, A, y, k, BigM);

	if (!x_relache.empty()) {
		node->is_solution = false;
		node->is_feasible = true;
		node->born_inf = pow(norm((*y) - (*A) * x_relache, 2), 2);
	}

	vec b = abs(x_relache) / BigM;


	if (is_binary(node->b_relache)) {
		node->is_solution = true;
	}

	return x_relache;
}

/************************    BranchandBound algo   ***************************************/
vec BranchandBound_l2l0(mat *A, vec *y, int k, double BigM) {

	clock_t t0;
	t0 = clock();
	double current_time = 0;

	vector<Node> nodeList;
	Node currentnode, Upnode, downnode;
	vec x_opt = zeros<vec> (A->n_cols);
	vec x_relache;
	int qj;
	int iter = 0;

	/******** Initialisation du noeud racine *********/
	Node rootNode;
	rootNode.born_inf = -exp10(10);
	rootNode.vec_Q0.clear();
	rootNode.vec_Q1.clear();
	for (unsigned int i = 0; i < A->n_cols; i++)
		rootNode.vec_Qv.push_back(i);
	rootNode.qj_parent = -1;
	//rootNode.BB=AA;
	//rootNode.Bz=A->t()* *y;
	/**************************************************/

	nodeList.push_back(rootNode);

	while (!nodeList.empty()) {
		iter++;
		currentnode = nodeList.back();
		nodeList.pop_back();

		current_time = (float)(clock()-t0)/CLOCKS_PER_SEC ;
				if (current_time> TimeBBMax){
						BBNodenum = iter;
						return x_opt;
					}


				if(iter %300==0){
					cout<<"N node : " << iter << " times: "<< current_time <<endl;
				}

		x_relache = evaluation_l2l0(&currentnode, A, y, k, BigM);
		//x_relache.print("x_relache");
		//cout<<"currentnode.born_inf :" << currentnode.born_inf<<endl;

		if (currentnode.is_solution && currentnode.born_inf < born_sup) {
			x_opt = x_relache;
			born_sup = currentnode.born_inf;
			cout << " is solution_ bornsupp: " << born_sup << endl;
			BestNodenum = iter;
		} else {
			if (currentnode.is_feasible && currentnode.born_inf < born_sup) {
				//choisir la variable de branchement
				qj = branching_rules(&currentnode, x_relache);

				Node filsNodeUp, filsNodeDown;
				// filsNodeUp
				filsNodeUp.branch = 1;
				filsNodeUp.num = 2 * currentnode.num;
				filsNodeUp.vec_Q1 = currentnode.vec_Q1;
				filsNodeUp.vec_Q1.push_back(qj);
				filsNodeUp.vec_Q0 = currentnode.vec_Q0;
				filsNodeUp.vec_Qv = currentnode.vec_Qv;
				filsNodeUp.vec_Qv.erase(
						remove(filsNodeUp.vec_Qv.begin(),
								filsNodeUp.vec_Qv.end(), qj),
						filsNodeUp.vec_Qv.end());
				//	filsNodeUp.qj_parent=qj;
				//	filsNodeUp.invAAQ1_parent=currentnode.invAAQ1;


				// filsNodeDown
				filsNodeDown.num = 2 * currentnode.num + 1;
				filsNodeDown.vec_Q1 = currentnode.vec_Q1;
				filsNodeDown.vec_Q0 = currentnode.vec_Q0;
				filsNodeDown.vec_Q0.push_back(qj);
				filsNodeDown.vec_Qv = filsNodeUp.vec_Qv;
				//	filsNodeDown.qj_parent=qj;
				//	filsNodeDown.branch=0;
				filsNodeDown.invAAQ1_parent = currentnode.invAAQ1;

				//	uvec Q1= conv_to <uvec>::from(currentnode.vec_Q1);
				//	uvec Q0= conv_to <uvec>::from(currentnode.vec_Q0);
				//	int qj_indice = qj- ((uvec)find(Q1<qj)).n_rows-((uvec)find(Q0<qj)).n_rows;

				//Q1.print("Q1");
				//Q0.print("Q0");
				//cout<<"qj: "<< qj <<" qj_indice: "<< qj_indice<<endl;

				//int x;
				//cin >> x;
				//filsNodeDown.BB = currentnode.BB;
				//filsNodeDown.BB.shed_col(qj_indice);
				//filsNodeDown.BB.shed_row(qj_indice);IloEnv env;
				//filsNodeDown.Bz = currentnode.Bz;
				//filsNodeDown.Bz.shed_row(qj_indice);


				// add node shild to List

				nodeList.push_back(filsNodeDown);
				nodeList.push_back(filsNodeUp);

			}
		}

		//delete &currentnode;
	}

	BBNodenum = iter;
	return x_opt;
}

