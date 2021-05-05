
/****************************    relaxation     ****************************************/
vec relaxation_l2pl0(Node *node, mat *A, vec *y, double lambda, double BigM) {
	clock_t th, thf;
	uvec Q1, Q0, Qv;
	vec x_Qv, x, b, w,score;
	mat A_Q1, A_Qv;

	x = zeros(A->n_cols);
	score = x;
	b = zeros(A->n_cols);



	Q1 = conv_to<uvec>::from(node->vec_Q1);
	Q0 = conv_to<uvec>::from(node->vec_Q0);
	Qv = conv_to<uvec>::from(node->vec_Qv);

	b(Q1)=ones(Q1.n_elem);




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

			w = optim_QP_bigM(A_Q1, *y, BigM);
		}




	double local_born_supp;

	    if(Q1.n_elem > 0){
			local_born_supp = pow(norm( (*y) - A_Q1 * w, 2), 2) + lambda * Q1.n_elem  - norm2carre_y ;

			if( local_born_supp < born_sup ){


				born_sup = local_born_supp;

				x_opt.zeros();
				x_opt(Q1)=w;


				nbr_update_Bsupp ++;
				T_best_Bsupp= current_time;
				cout << " Update born_sup " << born_sup << endl;

			}
	    }



	    if (opt_relaxation_with_cplex) {
	    		x=optim_withCplex_model_bigM_l2pl0(node ,*A, *y, lambda, BigM, Q1, Q0);


	    		if(Qv.n_elem)
	    		b(Qv) = abs(x(Qv)) / BigM;

	    		node->x_relache = x;
	    		node->b_relache = b;

	    		node->score = abs(x) ;
	    		return x;
	    	}



		th = clock();

		vec score_tmp;

		if(opt_relaxation_with_hom){
		x_Qv = mixed_BigM_homotopytemps_l2pl0(A_Q1, &w, &A_Qv, Q1, Qv, *y,
				BigM, lambda/(2*BigM) ,&score_tmp);

		score(Qv) = score_tmp;

		}else{
		//x_Qv.t().print("x_hom");

		vec x_Qv0 ;
		x_Qv0.reset();

		if(warm_restart){
			x_Qv0 = node->x_relache(Qv);

			if(node->branch==1){

			}else{

			}
		}

			if(opt_relaxation_with_actset)
				x_Qv =	mixed_BigM_featuresign_l2pl0(A_Q1, &w, &A_Qv, x_Qv0 , Q1, Qv, *y, BigM, lambda/(2*BigM), node);
		//x_Qv.t().print("x_act");
				score(Qv) =x_Qv;
		}

		thf = clock();
		T_homo += (float) (thf - th) / CLOCKS_PER_SEC;

		//x_Qv=homotopy(&BB,&Bz,Tau);

		x(Qv) = x_Qv;


	//}
	x(Q1) = w;

	b(Qv) = abs(x(Qv)) / BigM;


	node->x_relache = x;
    node->b_relache = b;
    node->score = score ;





	return x;
}

/****************************    Node evaluate   *********************************/
vec evaluation_l2pl0(Node *node, mat *A, vec *y, double lambda, double BigM) {

	vec x_relache = zeros(A->n_cols);

	x_relache = relaxation_l2pl0(node, A, y, lambda, BigM);

	if (!x_relache.empty()) {
		node->is_solution = false;
		node->is_feasible = true;
		node->born_inf = pow(norm( (*y) - (*A) * node->x_relache, 2), 2) + lambda * sum(node->b_relache)  - norm2carre_y;


	}

	if (is_binary(node->b_relache)) {
	//	cout<< "node is solution" << endl;
		node->is_solution = true;
	}

	return x_relache;
}

/************************    BranchAndBound algo   ********************************/
vec BranchandBound_l2pl0(mat *A, vec *y, double lambda, double BigM, vec x0) {

	clock_t t0;
	t0 = clock();

	current_time =0;
	nbr_update_Bsupp=0;
	T_best_Bsupp=0;

	vector<Node> nodeList;
	Node currentnode, Upnode, downnode;
	x_opt = x0;
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
	rootNode.x_relache = zeros(A->n_cols);
	rootNode.b_relache = zeros(A->n_cols);

	//	rootNode.BB=AA;
	//	rootNode.Bz=A->t()* *y;
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


				if(iter %1000==0){
					cout<<"N node : " << iter << " times: "<< current_time <<endl;
				}




		x_relache = evaluation_l2pl0(&currentnode, A, y, lambda, BigM);


		if (currentnode.is_solution && currentnode.born_inf <= born_sup) {
			x_opt = x_relache;
			born_sup = currentnode.born_inf;
			cout <<"Node: "<<iter << " is solution_ bornsupp: " << born_sup << endl;
			nbr_update_Bsupp ++;
			T_best_Bsupp= current_time;
			BestNodenum = iter;

		} else {
			if (currentnode.is_feasible && currentnode.born_inf < born_sup ) {

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

				filsNodeUp.x_relache= currentnode.x_relache;
				filsNodeUp.b_relache= currentnode.b_relache;
				//	filsNodeUp.qj_parent=qj;
				//	filsNodeUp.invAAQ1_parent=currentnode.invAAQ1;


				// filsNodeDown
				filsNodeDown.branch=0;
				filsNodeDown.num = 2 * currentnode.num + 1;
				filsNodeDown.vec_Q1 = currentnode.vec_Q1;
				filsNodeDown.vec_Q0 = currentnode.vec_Q0;
				filsNodeDown.vec_Q0.push_back(qj);
				filsNodeDown.vec_Qv = filsNodeUp.vec_Qv;
				//	filsNodeDown.qj_parent=qj;
				//	filsNodeDown.branch=0;
				filsNodeDown.invAAQ1_parent = currentnode.invAAQ1;
				//filsNodeDown.invAAQv_in_parent = currentnode.invAAQv_in;
				filsNodeDown.x_relache= currentnode.x_relache;
				filsNodeDown.b_relache= currentnode.b_relache;

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

