#include <time.h>
#include <iostream>
#include <armadillo>

using namespace arma;




/****************************    relaxation     ****************************************/
vec relaxation_l0l2(Node *node,mat *A,vec *y, double epsilon, double BigM){
	clock_t th,thf;
	uvec Q1,Q0,Qv;
	vec x_Qv,x,b,w,score;
	mat A_Q1,A_Qv;


	x=zeros(A->n_cols);
	score = x;
	b = zeros(A->n_cols);

	//double Tau=(k- node->vec_Q1.size())*BigM;


	Q1= conv_to <uvec>::from( node->vec_Q1);
	Q0= conv_to <uvec>::from( node->vec_Q0);
	Qv= conv_to <uvec>::from( node->vec_Qv);

	b(Q1)= ones(Q1.n_elem);

	A_Q1=A->cols(Q1);
	A_Qv=A->cols(Qv);

if(opt_relaxation_with_cplex){

	x = optim_withCplex_model_bigM_l0l2 (*A, *y,  epsilon, BigM, Q1,Q0);

	if(norm(*y- A_Q1* x(Q1) ,2)<=epsilon){
			x(Qv).zeros();
	        node->born_inf = Q1.size();
	        node->is_solution = 1;
	        return x;
	 }


	if(Qv.n_elem){
	node->born_inf= norm(x(Qv),1)/BigM + Q1.size() ;
	node->is_feasible=1;
	}else{
	node->is_feasible=0;
	}

	return x;
}


//x=zeros(A->n_cols);


	if(node->branch==0){
		node->invAAQ1=node->invAAQ1_parent;
	}else{
		node->invAAQ1=inv(AA.submat(Q1,Q1));//version recursive!!
	}

	//node->invAAQ1=inv(AA.submat(Q1,Q1));

	w= node->invAAQ1*(A_Q1.t()* *y);
	//w.t().print("w");



	if(!w.is_empty())
	if(BigM< max(abs(w))  ){
		w= optim_QP_bigM( A_Q1 , *y,  BigM);
	}


	 if(norm(*y- A_Q1* w ,2)<=epsilon){
	        x(Q1)=w;
	        node->born_inf = Q1.size();
	        node->is_solution = 1;
	        return x;
	 }


	 vec score_tmp;

	 th=clock();
	 x_Qv=  mixed_BigM_homotopytemps_l0l2(A_Q1, &w, &A_Qv,Q1, Qv, *y, BigM, epsilon, &score_tmp);
	// x_Qv.print("x_Qv");

	 score(Qv) = score_tmp;

	 thf=clock();
	 T_homo +=(float)(thf-th)/CLOCKS_PER_SEC;

	 //x_Qv=homotopy(&BB,&Bz,Tau);

	  x(Qv)=x_Qv;
	  x(Q1)=w;


	  if(Qv.n_elem)
	  b(Qv)= abs(x(Qv))/BigM ;

	  node->x_relache = x ;
	  node->b_relache = b ;
	  node->score = score ;

	  node->born_inf=  sum(b) ;
	  node->is_feasible=true;

	  return x;
}



/****************************    Node evaluate   ***************************************/
vec evaluation_l0l2(Node *node,mat *A,vec *y, double epsilon, double BigM){

vec x_relache=zeros(A->n_cols);

		x_relache=relaxation_l0l2(node,A,y,epsilon,BigM);


return x_relache;
}

/************************    BranchandBound algo   ***************************************/
vec BranchandBound_l0l2 (mat *A,vec *y, double epsilon , double BigM){
	//cout<<"BranchandBound"<<endl;
	clock_t t0;
	t0=clock();
	double current_time=0;

	vector<Node> nodeList;
	Node currentnode, Upnode , downnode;
	vec x_opt=zeros<vec>(A->n_cols);
	vec x_relache;
	int qj;
	int iter =0;



	/******** Initialisation du noeud racine *********/
	Node rootNode;
	rootNode.born_inf=exp10(10);
	rootNode.vec_Q0.clear();
	rootNode.vec_Q1.clear();
	for(unsigned int i=0; i< A->n_cols;i++)
	rootNode.vec_Qv.push_back(i);
	rootNode.qj_parent=-1;
	//rootNode.BB=AA;
	//rootNode.Bz=A->t()* *y;
	/**************************************************/

	nodeList.push_back(rootNode);



	while(! nodeList.empty()){
		iter++;
		currentnode= nodeList.back();
		nodeList.pop_back();

		current_time = (float)(clock()-t0)/CLOCKS_PER_SEC ;
		if (current_time> TimeBBMax){
				BBNodenum = iter;
				return x_opt;
			}


		if(iter %1000==0){
			cout<<"N node : " << iter << " times: "<< current_time <<endl;
		}

		//cout<<currentnode.branch<<endl;

		x_relache=evaluation_l0l2(&currentnode,A,y,epsilon,BigM);
	//	x_relache.print("x_relache");

		if(currentnode.is_solution && currentnode.born_inf<born_sup){
			x_opt=x_relache;
			born_sup=currentnode.born_inf;

			cout<<" is solution_ bornsupp: "<<born_sup<<endl;
			BestNodenum = iter;
		}else{
			if(currentnode.is_feasible && currentnode.born_inf <= born_sup - cut){
				//choisir la variable de branchement
				qj=branching_rules(&currentnode,x_relache);

				Node filsNodeUp,filsNodeDown;
				// filsNodeUp
				filsNodeUp.branch=1;
				filsNodeUp.num=2*currentnode.num;
				filsNodeUp.vec_Q1=currentnode.vec_Q1;
				filsNodeUp.vec_Q1.push_back(qj);
				filsNodeUp.vec_Q0=currentnode.vec_Q0;
				filsNodeUp.vec_Qv=currentnode.vec_Qv;
				filsNodeUp.vec_Qv.erase(remove(filsNodeUp.vec_Qv.begin(), filsNodeUp.vec_Qv.end(), qj), filsNodeUp.vec_Qv.end());


				// filsNodeDown
				filsNodeDown.num=2*currentnode.num+1;
				filsNodeDown.vec_Q1=currentnode.vec_Q1;
				filsNodeDown.vec_Q0=currentnode.vec_Q0;
				filsNodeDown.vec_Q0.push_back(qj);
				filsNodeDown.vec_Qv=filsNodeUp.vec_Qv;
				filsNodeDown.invAAQ1_parent=currentnode.invAAQ1;

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

