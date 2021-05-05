

vec calculr(vec y,mat B,mat H,vec w,vec x,uvec ind_w_M,uvec ind_x_M,double M,int N){

	vec a1,a2;
	if (! ind_x_M.is_empty() )
		a1=H.cols(ind_x_M)*(M*sign(x(ind_x_M)));
	else
		a1=zeros(N);


	if (! ind_w_M.is_empty())
		a2=B.cols(ind_w_M)*(M*sign(w(ind_w_M)));
	else
		a2=zeros(N,1);




	return y-a2-a1;

}


double fobj_featuresign(vec x,vec w,mat B,mat H,vec y,double lambda){

	return 0.5 * pow(norm(y -H*x -B*w),2) + lambda * norm(x,1);

}

void updates_activeset_x(vec x, uvec &ind_x_0, uvec &ind_x_in, uvec &ind_x_M,double M, double eps ){

	vec Etat_var_x = ones(x.n_elem); //1 -> Xint ; 0-> (x == 0) ; 2 -> X_M

			ind_x_M = find(abs(abs(x) - M) <= eps);
			x.rows(ind_x_M) = (M * sign(x.rows(ind_x_M)));
			Etat_var_x.elem(ind_x_M) = 2 * ones(ind_x_M.n_rows);

			ind_x_0 = find(abs(x) <= eps);
			x(ind_x_0)=0* ones(ind_x_0.n_elem);
			Etat_var_x.elem(ind_x_0) = zeros(ind_x_0.n_rows);

			ind_x_in = find(Etat_var_x == 1);

}

vec mixed_BigM_featuresign_l2pl0(mat B, vec *w, mat *H,vec x, uvec uQ1, uvec uQv, vec y, double M, double lambda, Node *node) {

		double eps = exp10(-3);
		double eps_ = exp10(-10);

		double lambda0;


		vec Hy = abs(H->t() * (y - (B * *w)));


		uword I0 = Hy.index_max();
		lambda0 = Hy.at(I0);

		int N  = H->n_rows;
		int Q1 = H->n_cols;
		int Q2 = B.n_cols;

		if(x.is_empty())
			x = zeros(Q1);


		if(lambda >=lambda0)
		    return x;


		uvec ind_x_M, ind_w_M;
		uvec ind_x_in, ind_w_in;
		uvec ind_x_inp, ind_w_inp;

		uvec ind_x_0;




		vec Etat_var_x = ones(Q1); //1 -> Xint ; 0-> (x == 0) ; 2 -> X_M

		ind_x_M = find(abs(abs(x) - M) <= eps);
		x.rows(ind_x_M) = (M * sign(x.rows(ind_x_M)));
		Etat_var_x.elem(ind_x_M) = 2 * ones(ind_x_M.n_rows);

		ind_x_0 = find(abs(x) <= eps);
		x(ind_x_0)= zeros(ind_x_0.n_elem);
		Etat_var_x.elem(ind_x_0) = zeros(ind_x_0.n_rows);

		ind_x_in = find(Etat_var_x == 1);





		vec Etat_var_w = ones(Q2); //1 -> Wint ; 2 -> w_M

		ind_w_M = find(abs(abs(*w) - M) <= eps);
		w->rows(ind_w_M) = (M * sign(w->rows(ind_w_M)));

		Etat_var_w.elem(ind_w_M) = 2 * ones(ind_w_M.n_rows);
		ind_w_in = find(Etat_var_w == 1);


		mat BB = AA.submat(uQ1, uQ1);
		mat HH = AA.submat(uQv, uQv);
		mat HB = AA.submat(uQv, uQ1);
		mat BH = AA.submat(uQ1, uQv);




		mat inv_BB, R, F, S, inv_S,inv_Sp;

		clock_t t0, t1;


		t0 =clock();
		 if(! ind_w_in.is_empty()){

					inv_BB = inv(BB.submat(ind_w_in, ind_w_in));
					R =  B.cols(ind_w_in) * inv_BB;
					F = eyeN + R * -B.cols(ind_w_in).t() ;
					S = H->cols(ind_x_in).t() * F * H->cols(ind_x_in);
					inv_S = inv(S);
		 }else{
					S = HH.submat(ind_x_in,ind_x_in) ;//H->cols(ind_x_in).t() * H->cols(ind_x_in);
					inv_S=inv( S );
		 }

		 t1=clock();
					 T_test +=(float)(t1-t0)/CLOCKS_PER_SEC;

		vec r, c1, c2,d1,d2, e, grad,thetax;

		uword i1, i2, i3 ;
		int ie1, ie2, ie3 ;



		double sigma,sigma1,sigma2,sigma3;

		sigma = -1;
		sigma1 = -1;
		sigma2 = -1;
		sigma3 = -1;

		int maxit,it2,it=0;

		vec progress, progress_xMp, progress_xMn,	progress_wMp, progress_wMn;
		double min_progress_xMp, min_progress_xMn,	min_progress_wMp, min_progress_wMn, t;


		int ii1, ii2, ii3, ii4;
		uword id_min_progress_xMp, id_min_progress_xMn,	id_min_progress_wMp, id_min_progress_wMn, progress_id;
		uvec id_sort_progress;


		vec x_new, w_new, sx_temp, sw_temp;
		double progress_born=1;
		uword id_progress_born, id_lsearch, pos;

		double fobj_lsearch, fobj_temp, lsearch ;

		bool born_x, born_change;



		 thetax = sign(x);

		 		int niter = 1000;


		 		it=0;

		 r= calculr(y,B, *H, *w,x, ind_w_M, ind_x_M, M, N);


		 fobj_lsearch = fobj_featuresign(x, *w, B, *H, y, lambda);



		while(it<niter){



		    if(it>=1){
		        // Activation de variables
		        if(! ind_w_in.is_empty() ){
		            c1 = H->t()*r - HB.cols(ind_w_in)* w->rows(ind_w_in)- HH.cols(ind_x_in)*x(ind_x_in);
		            c2 = B.t()*r - BH.cols(ind_x_in)*x(ind_x_in)- BB.cols(ind_w_in)* w->rows(ind_w_in);
		        }else if(! ind_x_in.is_empty()){
		            c1 = H->t()*r - HH.cols(ind_x_in)*x(ind_x_in);
		            c2 = B.t()*r - BH.cols(ind_x_in)*x(ind_x_in);
		        }else{
		            c1 = H->t()*r;
		            c2 = B.t()*r;
		        }




		        e= abs(c1(ind_x_0))- lambda ;
		        i1 = arg_max_p(ind_x_0, e, &sigma1, &ie1, eps);


		        e= (-c1(ind_x_M) % sign(x(ind_x_M)) )+ lambda ;
		        i2 = arg_max_p(ind_x_M, e, &sigma2, &ie2, eps);


		        e= -c2(ind_w_M) % sign(w->rows(ind_w_M));
		        i3 = arg_max_p(ind_w_M, e, &sigma3, &ie3, eps);




		        if ( (sigma1<0) && (sigma2<0) && (sigma3<0) ){

		        	grad = -c1(ind_x_in) % sign(x(ind_x_in))+lambda;


		        	if(! grad.is_empty()){
		        		if(max(abs( grad )) < eps){
		        			nbr_iter += it;
		        			return x;

		        			node->ind_Qv_in = ind_x_in;
		        			node->ind_Qv_0 = ind_x_0;
		        			node->ind_Qv_M = ind_x_M;

		        			node->invAAQv_in = inv_S;


		        		}
		        		 sigma=-1;

		        	}else{

		        		nbr_iter += it;
		        		return x;

		        	}

		        }

		        //if(sigma1>0)
		        //sigma = sigma1;
		        //else
	        	//sigma = max(sigma2,sigma3);

		        sigma = max(sigma1, max(sigma2,sigma3));

		        if( (sigma == sigma1) && (ie1!= -1) ){

		        	thetax.row(i1) = sign(c1.row(i1));

		        	ind_x_inp= ind_x_in;
		        	inv_Sp=inv_S;
		        	pos =  deplacer(i1, &ind_x_0, &ind_x_in);
		        	if( pos != -1){

		        		if(! ind_w_in.is_empty()){


									inv_S = inversion_rec_add_P(inv_S,	H->cols(ind_x_inp).t() * F * H->col(i1),	H->col(i1).t() * F * H->col(i1));

						 }else{

							 	 	inv_S = inversion_rec_add_P(inv_S, 	H->cols(ind_x_inp).t()  * H->col(i1),	HH.submat(i1,i1,i1,i1) );

						 }
		        	}else{
		        		cout<<"deplacer1 fail"<<endl;


		        	}
		        		//pause_cin();


		        }else if( (sigma == sigma2) && (ie2!= -1) ){

		        	ind_x_inp= ind_x_in;
		        	inv_Sp=inv_S;
		        	if( deplacer(i2, &ind_x_M, &ind_x_in) != -1){
		        		if(! ind_w_in.is_empty()){

									inv_S = inversion_rec_add_P(inv_S,	H->cols(ind_x_inp).t() * F * H->col(i2),	H->col(i2).t() * F * H->col(i2));

						 }else{

									inv_S = inversion_rec_add_P(inv_S,	H->cols(ind_x_inp).t()  * H->col(i2),	H->col(i2).t()  * H->col(i2));

						 }

		        		r+= H->col(i2) * x.row(i2);

		        	}else	cout<<"deplacer2 fail"<<endl;



		        }else if( (sigma == sigma3) && (ie3!= -1) ){

		        	inv_Sp=inv_S;
		        	if( deplacer(i3, &ind_w_M, &ind_w_in) != -1){
		        		if(! ind_w_in.is_empty()){

									inv_BB = inv(BB.submat(ind_w_in, ind_w_in));
									R =  B.cols(ind_w_in) * inv_BB;
									F = eyeN + R * -B.cols(ind_w_in).t() ;
									S = H->cols(ind_x_in).t() * F * H->cols(ind_x_in);
									inv_S = inv(S);
						 }else{
									S =  H->cols(ind_x_in).t() * H->cols(ind_x_in);
									inv_S=inv( S );
						 }

		        	}else	cout<<"deplacer3 fail"<<endl;

		        	r+= B.col(i3) * w->row(i3);


		        }

		    }
		         maxit=1000;


		         it2=0;
		         while(it<niter){

		        	 it2++;
		        	 it++;


		        	 x_new = x;
		        	 w_new = *w;




		        	 if(! ind_w_in.is_empty()){
		        	            x_new(ind_x_in) = inv_S* (H->cols(ind_x_in).t() * F * r - lambda * thetax(ind_x_in));
		        	            w_new(ind_w_in) = R.t() *(r- H->cols(ind_x_in) * x_new(ind_x_in));
		        	 }else
		        	            x_new(ind_x_in) = inv_S* (H->cols(ind_x_in).t()*r - lambda * thetax(ind_x_in));




		        	        if( all(sign(x(ind_x_in)) == sign(x_new(ind_x_in))) || ind_x_in.is_empty() )
		        	            if (all(abs(x_new) <=M) || x_new.is_empty() )
		        	        //    	if(all(sign(w->rows(ind_w_in)) == sign(w_new(ind_w_in))) || ind_w_in.is_empty())
		        	                if (all(abs(w_new) <=M) || w_new.is_empty()){

		        	                    x=x_new;
		        	                    *w=w_new;
		        	                    break;
		        	                }



		        	        d1=(x_new(ind_x_in) - x(ind_x_in));
							d2=(w_new(ind_w_in) - w->rows(ind_w_in));



							progress_xMp = (+M  - x(ind_x_in)) /d1;
							progress_xMn = (-M  - x(ind_x_in)) /d1;
							progress_wMp = (+M  - w->rows(ind_w_in)) /d2;
							progress_wMn = (-M  - w->rows(ind_w_in)) /d2;

							progress_born = 1;


							id_min_progress_xMp = arg_min_p(ind_x_in, progress_xMp, &min_progress_xMp, &ii1, eps_);

							if(progress_born > min_progress_xMp){
								progress_born 	 = min_progress_xMp;
								id_progress_born = id_min_progress_xMp;
								born_x= 1;
							}

							id_min_progress_xMn = arg_min_p(ind_x_in, progress_xMn, &min_progress_xMn, &ii2, eps_);

							if(progress_born > min_progress_xMn){
								progress_born 	 = min_progress_xMn;
								id_progress_born = id_min_progress_xMn;
								born_x= 1;
							}

							id_min_progress_wMp = arg_min_p(ind_w_in, progress_wMp, &min_progress_wMp, &ii3, eps_);

							if(progress_born > min_progress_wMp){
								progress_born 	 = min_progress_wMp;
								id_progress_born = id_min_progress_wMp;
								born_x= 0;
							}

							id_min_progress_wMn = arg_min_p(ind_w_in, progress_wMn, &min_progress_wMn, &ii4, eps_);

							if(progress_born > min_progress_wMn){
								progress_born 	 = min_progress_wMn;
								id_progress_born = id_min_progress_wMn;
								born_x= 0;
							}





							progress = -x(ind_x_in) /d1;
							id_sort_progress = sort_index(progress);


							lsearch = 0;
							bool pasplus1=0;

							for(int i=0; i< progress.n_elem ; i++){

								t = progress(id_sort_progress(i));



								if (t <= 0)
									continue;

								if (t >= progress_born)
									break;


								//if(i1 == ind_x_in(id_sort_progress(i)) )
								//	continue;

								sx_temp= x;
								sw_temp= *w;

								sx_temp(ind_x_in)= x(ind_x_in)+ d1*t;
								sw_temp(ind_w_in) = w->rows(ind_w_in)+ d2*t;


								fobj_temp = fobj_featuresign(sx_temp, sw_temp, B, *H, y, lambda);


								if ( fobj_temp < fobj_lsearch ) {
									fobj_lsearch = fobj_temp;
									lsearch = t;
									id_lsearch = ind_x_in(id_sort_progress(i));
								}else{
									//if(pasplus1)
									if (ind_x_in(id_sort_progress(i))== i1)
																		continue;
									break;
									//pasplus1 = 1;
								}



							}

							// test obj de progress_born
								t = progress_born;

								sx_temp= x;
								sw_temp= *w;

								sx_temp(ind_x_in)= x(ind_x_in)+ d1*t;
								sw_temp(ind_w_in) = w->rows(ind_w_in)+ d2*t;


								fobj_temp = fobj_featuresign(sx_temp, sw_temp, B, *H, y, lambda);


								if ( fobj_temp <= fobj_lsearch ) {
									fobj_lsearch = fobj_temp;
									lsearch = t;
									id_lsearch = id_progress_born;
								}




							if ( (lsearch >0) && (lsearch < progress_born)){

								x(ind_x_in) = x(ind_x_in)+ ( (x_new(ind_x_in)- x(ind_x_in)) % (lsearch * ones(ind_x_in.n_elem)) ) ;
								w->rows(ind_w_in) = w->rows(ind_w_in)+ ( (w_new(ind_w_in)- w->rows(ind_w_in)) % (lsearch* ones(ind_w_in.n_elem) ) );
								thetax(ind_x_in) = sign(x(ind_x_in));


								born_change=0;
								pos = deplacer(id_lsearch, &ind_x_in, &ind_x_0);
								if( pos  != -1){


									thetax(id_lsearch) = 0;
									inv_S= inversion_rec_del_P(inv_S, pos);



									}else	cout<<"deplacer4 fail"<<endl;




							}else if((lsearch >0) && lsearch == progress_born ){

								x(ind_x_in) = x(ind_x_in)+ ( (x_new(ind_x_in)- x(ind_x_in)) % (lsearch * ones(ind_x_in.n_elem)) ) ;
								w->rows(ind_w_in) = w->rows(ind_w_in)+ ( (w_new(ind_w_in)- w->rows(ind_w_in)) % (lsearch* ones(ind_w_in.n_elem) ) );
								thetax(ind_x_in) = sign(x(ind_x_in));

								born_change=1;


								if (lsearch !=1){
									if(born_x==1){
										pos =  deplacer(id_progress_born, &ind_x_in, &ind_x_M);
										if( pos != -1){

											inv_S= inversion_rec_del_P(inv_S, pos);
											r-= H->col(id_progress_born) * x.row(id_progress_born);


										}else
											cout<<"deplacer 5 fail"<<endl;


									}else if(born_x==0)
										if( deplacer(id_progress_born, &ind_w_in, &ind_w_M) != -1){
											if(! ind_w_in.is_empty()){

														inv_BB = inv(BB.submat(ind_w_in, ind_w_in));
														R =  B.cols(ind_w_in) * inv_BB;
														F = eyeN + R * -B.cols(ind_w_in).t() ;
														S = H->cols(ind_x_in).t() * F * H->cols(ind_x_in);
														inv_S = inv(S);
											 }else{
														S =  H->cols(ind_x_in).t() * H->cols(ind_x_in);
														inv_S=inv( S );
									 			 }

												r-= B.col(id_progress_born) * w->row(id_progress_born);

										}else
											cout<<"deplacer 6 fail id:"<< id_progress_born << endl;



									w->rows(ind_w_M)= M* sign(w->rows(ind_w_M));
									x(ind_x_M)= M* sign(x(ind_x_M));

								}//else{
								//	if(it2==1)
								//	break;
								//}



							}else {

								if(it2==1){
										if( (sigma == sigma1) && (ie1!= -1) ){

														thetax.row(i1) = 0;

														pos =  deplacer(i1, &ind_x_in, &ind_x_0);
														if( pos != -1){
															inv_S= inversion_rec_del_P(inv_S, pos);
															//inv_S = inv_Sp;
														}else{
															cout<<"deplacer_1 fail"<<endl;
														}


										}else if( (sigma == sigma2) && (ie2!= -1)  ){
											if( deplacer(i2, &ind_x_in, &ind_x_M) != -1)
												inv_S=inv_Sp;
											else
												cout<<"deplacer_2 fail"<<endl;

										}else if( (sigma == sigma3)  && (ie3!= -1) ){

											if( deplacer(i3,&ind_w_in, &ind_w_M) != -1)
												inv_S=inv_Sp;
											else
												cout<<"deplacer_3 fail"<<endl;

										}

								}else{
									cout<<"lsearch ==0 it2=="<<it2<<endl;
									w->rows(ind_w_M)= M* sign(w->rows(ind_w_M));
									x(ind_x_M)= M* sign(x(ind_x_M));

								}


								}

						}


		         }

		cout<<"it : "<< it << endl;








}














/*

vec mixed_BigM_featuresign_l2pl0_old(mat B, vec *w, mat *H,vec x, uvec uQ1, uvec uQv, vec y, double M, double lambda) {

		double eps = exp10(-7);
		double lambda0;


		vec Hy = abs(H->t() * (y - (B * *w)));


		uword I0 = Hy.index_max();
		lambda0 = Hy.at(I0);

		int N  = H->n_rows;
		int Q1 = H->n_cols;
		int Q2 = B.n_cols;

		if(x.is_empty())
			x = zeros(Q1);


		if(lambda >=lambda0)
		    return x;


		uvec ind_x_M, ind_w_M;
		uvec ind_x_in, ind_w_in;
		uvec ind_x_0;


		vec Etat_var_x = ones(Q1); //1 -> Xint ; 0-> (x == 0) ; 2 -> X_M

		ind_x_M = find(abs(abs(x) - M) <= eps);
		x.rows(ind_x_M) = (M * sign(x.rows(ind_x_M)));
		Etat_var_x.elem(ind_x_M) = 2 * ones(ind_x_M.n_rows);

		ind_x_0 = find(x == 0);
		Etat_var_x.elem(ind_x_0) = zeros(ind_x_0.n_rows);

		ind_x_in = find(Etat_var_x == 1);


		vec Etat_var_w = ones(Q2); //1 -> Wint ; 2 -> w_M

		ind_w_M = find(abs(abs(*w) - M) <= eps);
		w->rows(ind_w_M) = (M * sign(w->rows(ind_w_M)));

		Etat_var_w.elem(ind_w_M) = 2 * ones(ind_w_M.n_rows);
		ind_w_in = find(Etat_var_w == 1);


		mat BB = AA.submat(uQ1, uQ1);
		mat HH = AA.submat(uQv, uQv);
		mat HB = AA.submat(uQv, uQ1);
		mat BH = AA.submat(uQ1, uQv);




		vec r, c1, c2,d1,d2, e, grad,thetax;

		uword i1, i2, i3 ;
		int ie1, ie2, ie3 ;

		r= calculr(y,B, *H, *w,x, ind_w_M, ind_x_M, M, N);


		double sigma,sigma1,sigma2,sigma3=0;
		int maxit,it,it2=0;

		vec progress, progress_xMp, progress_xMn,	progress_wMp, progress_wMn;
		double min_progress_xMp, min_progress_xMn,	min_progress_wMp, min_progress_wMn, t;


		int ii1, ii2, ii3, ii4;
		uword id_min_progress_xMp, id_min_progress_xMn,	id_min_progress_wMp, id_min_progress_wMn, progress_id;
		uvec id_sort_progress;


		vec x_new, w_new, sx_temp, sw_temp;
		double progress_born=1;
		uword id_progress_born, id_lsearch;

		double fobj_lsearch, fobj_temp, lsearch ;

		bool born_x;

		mat inv_BB, R, F, S, inv_S;

		thetax = sign(x);

		int niter = 1000;


		it=0;
		while(it<niter){

			it=it+1;

			if(it>100)
			cout << "it :" <<it<< endl;


			r= calculr(y,B, *H, *w,x, ind_w_M, ind_x_M, M, N);

		    if(it>1){
		        // Activation de variables
		        if(! ind_w_in.is_empty() ){
		            c1 = H->t()*r - HB.cols(ind_w_in)* w->rows(ind_w_in)- HH.cols(ind_x_in)*x(ind_x_in);
		            c2 = B.t()*r - BH.cols(ind_x_in)*x(ind_x_in)- BB.cols(ind_w_in)* w->rows(ind_w_in);
		        }else if(! ind_x_in.is_empty()){
		            c1 = H->t()*r - HH.cols(ind_x_in)*x(ind_x_in);
		            c2 = B.t()*r - BH.cols(ind_x_in)*x(ind_x_in);
		        }else{
		            c1 = H->t()*r;
		            c2 = B.t()*r;
		        }





		        e= abs(c1(ind_x_0))- lambda ;
		        i1 = arg_max_p(ind_x_0, e, &sigma1, &ie1, eps);


		     //   abs(c1(ind_x_0)).t().print("abs(c1(ind_x_0))");
		     //   cout <<"lambda : " << lambda<< endl;
		     //   e.t().print("e");

		        e= (-c1(ind_x_M) % sign(x(ind_x_M)) )+ lambda ;
		        i2 = arg_max_p(ind_x_M, e, &sigma2, &ie2, eps);

		      //  e.t().print("e");

		        e= -c2(ind_w_M) % sign(w->rows(ind_w_M));
		        i3 = arg_max_p(ind_w_M, e, &sigma3, &ie3, eps);

		      //  e.t().print("e");


		        if ( (sigma1<0) && (sigma2<0) && (sigma3<0) ){
		        	return x;
		        	grad = -c1(ind_x_in) % sign(x(ind_x_in))+lambda;



		        	if(! grad.is_empty()){
		        		if(max(abs( grad )) < eps)
		        			return x;
		        		else{
		        			cout<< "grad not null"<<endl;
		        			grad.t().print("grad");
		        		}

		        		 sigma=-1;

		        	}else
		        		return x;
		        }

		        sigma = max(sigma1, max(sigma2,sigma3));

		        if( (sigma == sigma1) && (ie1!= -1) ){

		        	thetax.row(i1) = sign(c1.row(i1));
		        	if(! deplacer(i1, &ind_x_0, &ind_x_in)){
		        		cout<<"deplacer1 fail"<<endl;
		        		//pause_cin();
		        	}

		        }else if( (sigma == sigma2) && (ie2!= -1) ){

		        	if(! deplacer(i2, &ind_x_M, &ind_x_in)){
		        		cout<<"deplacer2 fail"<<endl;
		        	//	pause_cin();
		        	}

		        }else if( (sigma == sigma3) && (ie3!= -1) ){

		        	if(! deplacer(i3, &ind_w_M, &ind_w_in)){
		        		cout<<"deplacer3 fail"<<endl;
		        		//pause_cin();
		        	}
		        }

		    }
		         maxit=1000;

		         it2=0;

		         while(it2<maxit){

		        	 if(it2> 10)
		        	 cout << "it2 :"<< it2 <<endl;

		        	 it2++;

		        	 r= calculr(y,B, *H, *w,x, ind_w_M, ind_x_M, M, N);

		        	 x_new = x;
		        	 w_new = *w;




		        	 if(! ind_w_in.is_empty()){

								inv_BB = inv(BB.submat(ind_w_in, ind_w_in));
								R =  B.cols(ind_w_in) * inv_BB;
								F = eyeN + R * -B.cols(ind_w_in).t() ;
								S = H->cols(ind_x_in).t() * F * H->cols(ind_x_in);
								inv_S = inv(S);




		        	            x_new(ind_x_in) = inv_S* (H->cols(ind_x_in).t() * F * r - lambda * thetax(ind_x_in));
		        	            w_new(ind_w_in) = R.t() *(r- H->cols(ind_x_in) * x_new(ind_x_in));




		        	 }else{
		        		 	 	S = H->cols(ind_x_in).t() * H->cols(ind_x_in);
		        	            inv_S=inv(S );
		        	            x_new(ind_x_in) = inv_S* (H->cols(ind_x_in).t()*r - lambda * thetax(ind_x_in));
		        	 }





		        	 fobj_lsearch = fobj_featuresign(x, *w, B, *H, y, lambda);

		        	// cout<<"fobj_featuresign : "<<fobj_lsearch<< endl;

		        	        if( all(sign(x(ind_x_in)) == sign(x_new(ind_x_in))) || ind_x_in.is_empty() )
		        	            if (all(abs(x_new) <=M) || x_new.is_empty() )
		        	            	if(all(sign(w->rows(ind_w_in)) == sign(w_new(ind_w_in))) || ind_w_in.is_empty())
		        	                if (all(abs(w_new) <=M) || w_new.is_empty()){

		        	                    x=x_new;
		        	                    *w=w_new;
		        	                    break;
		        	                }


		        	        d1=(x_new(ind_x_in) - x(ind_x_in));
							d2=(w_new(ind_w_in) - w->rows(ind_w_in));



							progress_xMp = (+M  - x(ind_x_in)) /d1;
							progress_xMn = (-M  - x(ind_x_in)) /d1;
							progress_wMp = (+M  - w->rows(ind_w_in)) /d2;
							progress_wMn = (-M  - w->rows(ind_w_in)) /d2;

							progress_born = 1;


							id_min_progress_xMp = arg_min_p(ind_x_in, progress_xMp, &min_progress_xMp, &ii1, eps);

							if(progress_born > min_progress_xMp){
								progress_born 	 = min_progress_xMp;
								id_progress_born = id_min_progress_xMp;
								born_x= 1;
							}

							id_min_progress_xMn = arg_min_p(ind_x_in, progress_xMn, &min_progress_xMn, &ii2, eps);

							if(progress_born > min_progress_xMn){
								progress_born 	 = min_progress_xMn;
								id_progress_born = id_min_progress_xMn;
								born_x= 1;
							}

							id_min_progress_wMp = arg_min_p(ind_w_in, progress_wMp, &min_progress_wMp, &ii3, eps);

							if(progress_born > min_progress_wMp){
								progress_born 	 = min_progress_wMp;
								id_progress_born = id_min_progress_wMp;
								born_x= 0;
							}

							id_min_progress_wMn = arg_min_p(ind_w_in, progress_wMn, &min_progress_wMn, &ii4, eps);

							if(progress_born > min_progress_wMn){
								progress_born 	 = min_progress_wMn;
								id_progress_born = id_min_progress_wMn;
								born_x= 0;
							}


							progress = -x(ind_x_in) /d1;
							id_sort_progress = sort_index(progress);

							lsearch = 0;

							for(int i=0; i< progress.n_elem ; i++){

								t = progress(id_sort_progress(i));


								if( (t <= 0) || (t >= progress_born) )
									continue;

								sx_temp= x;
								sw_temp= *w;

								sx_temp(ind_x_in)= x(ind_x_in)+ d1*t;
								sw_temp(ind_w_in) = w->rows(ind_w_in)+ d2*t;


								fobj_temp = fobj_featuresign(sx_temp, sw_temp, B, *H, y, lambda);


								if ( fobj_temp < fobj_lsearch ) {
									fobj_lsearch = fobj_temp;
									lsearch = t;
									id_lsearch = ind_x_in(id_sort_progress(i));
								}else
									break;

							}

							// test obj de progress_born

							t = progress_born;

							sx_temp= x;
							sw_temp= *w;

							sx_temp(ind_x_in)= x(ind_x_in)+ d1*t;
							sw_temp(ind_w_in) = w->rows(ind_w_in)+ d2*t;


							fobj_temp = fobj_featuresign(sx_temp, sw_temp, B, *H, y, lambda);


							if ( fobj_temp < fobj_lsearch ) {
								fobj_lsearch = fobj_temp;
								lsearch = t;
								id_lsearch = id_progress_born;
							}



							if ( (lsearch >0) && (lsearch < progress_born)){

								x(ind_x_in) = x(ind_x_in)+ ( (x_new(ind_x_in)- x(ind_x_in)) % (lsearch * ones(ind_x_in.n_elem)) ) ;
								w->rows(ind_w_in) = w->rows(ind_w_in)+ ( (w_new(ind_w_in)- w->rows(ind_w_in)) % (lsearch* ones(ind_w_in.n_elem) ) );
								thetax(ind_x_in) = sign(x(ind_x_in));


								//ind_x_in.t().print("ind_x_in");
								//cout << "id_lsearch : "<< id_lsearch <<endl;

								if(! deplacer(id_lsearch, &ind_x_in, &ind_x_0)){
									cout<<"deplacer4 fail"<<endl;
									//pause_cin();

								}

							}else if((lsearch >0) && lsearch == progress_born ){

								x(ind_x_in) = x(ind_x_in)+ ( (x_new(ind_x_in)- x(ind_x_in)) % (lsearch * ones(ind_x_in.n_elem)) ) ;
								w->rows(ind_w_in) = w->rows(ind_w_in)+ ( (w_new(ind_w_in)- w->rows(ind_w_in)) % (lsearch* ones(ind_w_in.n_elem) ) );
								thetax(ind_x_in) = sign(x(ind_x_in));


								if (lsearch !=1){
									if(born_x==1)
										if(! deplacer(id_progress_born, &ind_x_in, &ind_x_M))
											cout<<"deplacer 5 fail"<<endl;


									else if(born_x==0)
										if(! deplacer(id_progress_born, &ind_w_in, &ind_w_M))
											cout<<"deplacer 6 fail id:"<< id_progress_born << endl;



									w->rows(ind_w_M)= M* sign(w->rows(ind_w_M));
									x(ind_x_M)= M* sign(x(ind_x_M));

								}



							}else if(lsearch == 0)
								break;

							//	if(it2>100){
							//		cout <<"progress_born : "<< progress_born <<endl;
							//		cout <<"lsearch : "<< lsearch <<endl;
							//	}

							if ( ind_x_in.is_empty() && ind_w_in.is_empty() )
							   break;


		         }





		}



}
*/
