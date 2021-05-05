
int itmax = 0;

vec mixed_BigM_homotopytemps_l2l0(mat B, vec *w, mat *H, uvec uQ1, uvec uQv, vec y,
		double M, double Tau,vec *score) {
	clock_t t0, t1;

	double eps = exp10(-7);
	double lambda, gamma, nx;
	double gamma_xint_plus_P, gamma_xint_plus_N, gamma_x0_plus, gamma_xM_plus,
			gamma_xM_moin_N, gamma_xM_moin_P, gamma_wM_plus, gamma_wM_moin;

	uword i1, i2, i3, i4, i5, i6, i7, i8;
	int i1_, i2_, i3_, i4_, i5_, i6_, i7_, i8_;

	vec Hy = abs(H->t() * (y - (B * *w)));
	uword I0 = Hy.index_max();
	lambda = Hy.at(I0);

	mat BB = AA.submat(uQ1, uQ1);
	mat HH = AA.submat(uQv, uQv);
	mat HB = AA.submat(uQv, uQ1);
	mat BH = AA.submat(uQ1, uQv);

	int N = H->n_rows;
	int Q1 = H->n_cols;
	int Q2 = B.n_cols;

	sp_mat eyeN1;
	eyeN1 = eyeN(span(0, N - 1), span(0, N - 1));
	vec Hr, Br, res;
	mat inv_BB,R,F,P, S , inv_S;

	int niter = 1000;
	int it_non_elag= 1000;
	if(non_elag_hom)
		it_non_elag=50;

	vec x = zeros(Q1);
	*score =x;

	vec r, a1, a2, d1, d2, v1, v2, vv1, vv2, e, c, cc;
	uvec IM1, IM2, Iint2;
	uvec Iint1p, Iint1 = regspace<uvec> (0, 0);
	Iint1(0) = I0;

	vec Etat_var_x = zeros(Q1); //1 -> Xint ; 0-> (x == 0) ; 2 -> X_M
	vec Etat_var_w = ones(Q2); //1 -> Wint ; 2 -> w_M

	IM2 = find(abs(abs(*w) - M) <= exp10(-7));
	w->rows(IM2) = (M * sign(w->rows(IM2)));

	Etat_var_w.elem(IM2) = 2 * ones(IM2.n_rows);
	Iint2 = find(Etat_var_w == 1);

	a1 = zeros(N);

	d1 = zeros(Q1);
	d2 = zeros(Q2);



	if (!IM2.is_empty()) {
		a2 = B.cols(IM2) * (w->rows(IM2));
		r = y - a2;

	} else {
		a2 = zeros(N);
		r = y;
	}

	 Hr = H->t() * r;
	 Br = B.t() * r;

	inv_BB = inv(BB.submat(Iint2, Iint2));
	R = B.cols(Iint2) * inv_BB;
	P = R * -B.cols(Iint2).t();

	F = eyeN1 + P;

	if(!Iint2.is_empty())
	 S = H->cols(Iint1).t() * F * H->cols(Iint1);
	else
	 S = HH.submat(Iint1,Iint1);

	 inv_S = inv(S);


	c = Hr - HB.cols(Iint2) * w->rows(Iint2);// - HH.cols(Iint1) * x.rows(Iint1);
	d1.elem(Iint1) = inv_S * sign(c.rows(Iint1));

	uvec J, es; // arg(x == 0)
	Etat_var_x.at(I0) = 1;

	Iint1 = find(Etat_var_x == 1);
	IM1 = find(Etat_var_x == 2);
	IM2 = find(Etat_var_w == 2);
	Iint2 = find(Etat_var_w == 1);
	J = find(Etat_var_x == 0); // à optimiser

	for (int it = 0; it < niter; it++) {

		if (it > itmax)
			itmax = it;
		//J = find(Etat_var_x == 0); // à optimiser

		v1 = HH.submat(J, Iint1) * d1.rows(Iint1);

		if(!Iint2.is_empty())
		v2 = HB.submat(J,Iint2) * inv_BB * (BH.submat(Iint2, Iint1) * d1.rows(Iint1));
		//v2 = H->cols(J).t() * R * (BH.submat(Iint2, Iint1) * d1.rows(Iint1));
		else
		v2 = zeros(J.n_elem);


		e = (lambda - c.rows(J)) / (1 + v2 - v1);
		i1 = arg_min_p(J, e, &gamma_xint_plus_P, &i1_, eps);

		e = (lambda + c.rows(J)) / (1 - v2 + v1);
		i2 = arg_min_p(J, e, &gamma_xint_plus_N, &i2_, eps);

		e = -x.rows(Iint1) / d1.rows(Iint1);
		i3 = arg_min_p(Iint1, e, &gamma_x0_plus, &i3_, eps);

		e = (M * sign(x.rows(Iint1)) - x.rows(Iint1)) / d1.rows(Iint1);
		i4 = arg_min_p(Iint1, e, &gamma_xM_plus, &i4_, eps);


		gamma = min(gamma_xint_plus_P,min(gamma_xint_plus_N, min(gamma_x0_plus, gamma_xM_plus)));

		if (!IM1.is_empty()) {

			vv1 = HH.submat(IM1, Iint1) * d1.rows(Iint1);
			vv2 = HB.submat(IM1,Iint2) * inv_BB * (BH.submat(Iint2, Iint1) * d1.rows(Iint1));

			e = (lambda - c.rows(IM1)) / (1 - vv2 + vv1);
			i5 = arg_min_p(IM1, e, &gamma_xM_moin_N, &i5_, eps);

			e = (lambda + c.rows(IM1)) / (1 + vv2 - vv1);
			i6 = arg_min_p(IM1, e, &gamma_xM_moin_P, &i6_, eps);

		}

		if (!Iint2.is_empty()) {
			d2(Iint2) = - inv_BB.t() * (BH.submat(Iint2,Iint1) * d1.rows(Iint1));
		//	d2(Iint2) = -( R.t() * (H->cols(Iint1) * d1.rows(Iint1)));

			e = (M * sign(w->rows(Iint2)) - w->rows(Iint2)) / d2.rows(Iint2);
			i7 = arg_min_p(Iint2, e, &gamma_wM_plus, &i7_, eps);
			gamma = min(gamma, gamma_wM_plus);
		}

		if (!IM2.is_empty()) {

			if (!Iint2.is_empty()) {
				cc = Br.rows(IM2) - BB.submat(IM2, Iint2) * w->rows(Iint2) - BH.submat(IM2, Iint1) * x.rows(Iint1);
			} else {
				cc = Br.rows(IM2) - BH.submat(IM2, Iint1) * x.rows(Iint1);
			}

			vv1 = BH.submat(IM2, Iint1) * d1.rows(Iint1);
			vv2 = B.cols(IM2).t() * R * (BH.submat(Iint2, Iint1) * d1.rows(Iint1));

			e = (cc) / (-vv2 + vv1);
			i8 = arg_min_p(IM2, e, &gamma_wM_moin, &i8_, eps);

			gamma = min(gamma, gamma_wM_moin);

		}

		//  any condition is in force
		if (gamma == exp10(10)) {
		//	cout << "aucune condition";
			return x;
		}

		//  new solution
		x.elem(Iint1) = x.rows(Iint1) + gamma * d1.rows(Iint1);
		*w += gamma * d2;

		score->elem(Iint1) += abs(x.elem(Iint1));

		lambda = lambda - gamma;

		nx = norm(x.rows(Iint1), 1);
		if (nx >= Tau) {
			if (nx > Tau) {
				//<<"fin";
				vec sn = sign(x.rows(Iint1) - gamma * d1.rows(Iint1) / 2);
				vec gm = (nx - Tau) / (d1.rows(Iint1).t() * sn);
				x.rows(Iint1) -= d1.rows(Iint1) * gm;
				*w -= d2 * gm;
			}
			return x;
		}


		if( it>it_non_elag ){
		res = y- H->cols(Iint1) * x.rows(Iint1)- B * *w;
		if (pow(norm(res,2),2) < born_sup ){
			//cout << " itg :"<< it<<endl;
			return x;
		   }
		it_non_elag+=10;
		}

		//  lambda positive
		if (lambda < eps) {
		//	cout << "aucune condition";
			return x;
		}

		t0 =clock();

		//  mettre à jour selon le cas (gamma)
		if (gamma == gamma_x0_plus) {

			Etat_var_x.at(i3) = 0;
			deplacer(i3, &Iint1, &J);
			x.at(i3) = 0;

			if (!Iint2.is_empty()) {
				c = Hr - (HB.cols(Iint2) * w->rows(Iint2)) - (HH.cols(Iint1)
						* x.rows(Iint1));
			} else {
				c = Hr - HH.cols(Iint1) * x.rows(Iint1);
			}

			//  % update direction Iint
			d1 = zeros(Q1);
			d2 = zeros(Q2);

			inv_S= inversion_rec_del_P(inv_S, i3_);

			d1.elem(Iint1) = inv_S * sign(c.rows(Iint1));


		} else if (gamma == gamma_xint_plus_P || gamma == gamma_xint_plus_N) {

			if (gamma == gamma_xint_plus_N)
				i1 = i2;

			Etat_var_x.at(i1) = 1;
			Iint1p = Iint1;


			deplacer(i1, &J, &Iint1);

			if (!Iint2.is_empty()) {
				c = Hr - (HB.cols(Iint2) * w->rows(Iint2)) - (HH.cols(Iint1)
						* x.rows(Iint1));
			} else {
				c = Hr - HH.cols(Iint1) * x.rows(Iint1);
			}

			//  % update direction Iint
			d1 = zeros(Q1);
			d2 = zeros(Q2);

			inv_S = inversion_rec_add_P(inv_S,	H->cols(Iint1p).t() * F * H->col(i1),	H->col(i1).t() * F * H->col(i1));

			d1.elem(Iint1) = inv_S * sign(c.rows(Iint1));

		} else if (gamma == gamma_xM_plus) {
			//cout << "gamma_xM_plus"<<endl;

			Etat_var_x.at(i4) = 2;
			deplacer(i4, &Iint1, &IM1);

			x.row(i4) = (M * sign(x.row(i4)));
			a1 += H->col(i4) * x.row(i4);
			r = y - a1 - a2;
			Hr = H->t() * r;

			if (!Iint2.is_empty()) {
				c = Hr - (HB.cols(Iint2) * w->rows(Iint2)) - (HH.cols(Iint1)
						* x.rows(Iint1));
			} else {
				c = Hr - HH.cols(Iint1) * x.rows(Iint1);
			}

			//  % update direction Iint
			d1 = zeros(Q1);
			d2 = zeros(Q2);


			inv_S= inversion_rec_del_P(inv_S, i4_ );

			//S = H->cols(Iint1).t() * F * H->cols(Iint1);
			//inv_S = inv(S);

			d1.elem(Iint1) = inv_S * sign(c.rows(Iint1));

		} else if (gamma == gamma_xM_moin_N || gamma == gamma_xM_moin_P) {

			if (gamma == gamma_xM_moin_P)
				i5 = i6;
			//cout << "gamma_xM_moin_N"<<endl;

			Etat_var_x.at(i5) = 1;
			deplacer(i5, &IM1, &Iint1);

			x.row(i5) = (M * sign(x.row(i5)));
			a1 -= H->col(i5) * x.row(i5);
			r = y - a1 - a2;
			Hr = H->t() * r;

			c = Hr - (HB.cols(Iint2) * w->rows(Iint2)) - (HH.cols(Iint1)
					* x.rows(Iint1));

			d1 = zeros(Q1);
			d2 = zeros(Q2);

			S = H->cols(Iint1).t() * F * H->cols(Iint1);
			inv_S = inv(S);
			d1.elem(Iint1) = inv_S * sign(c.rows(Iint1));

		} else if (gamma == gamma_wM_plus) {
			//cout << "gamma_wM_plus it: "<<it<<endl;



			Etat_var_w.at(i7) = 2;
			deplacer(i7, &Iint2, &IM2);

			w->row(i7) = (M * sign(w->row(i7)));
			a2 += B.col(i7) * w->row(i7);
			r = y - a1 - a2;
			Hr = H->t() * r;
			Br = B.t() * r;

			inv_BB = inv(BB.submat(Iint2, Iint2));

			R =  B.cols(Iint2) * inv_BB;
			F = eyeN + R * -B.cols(Iint2).t() ;


			S = H->cols(Iint1).t() * F * H->cols(Iint1);


			if (!Iint2.is_empty()) {
				c = Hr - (HB.cols(Iint2) * w->rows(Iint2)) - (HH.cols(Iint1)
						* x.rows(Iint1));
			} else {
				c = Hr - HH.cols(Iint1) * x.rows(Iint1);
			}

			// % update direction Iint
			d1 = zeros(Q1);
			d2 = zeros(Q2);

			inv_S = inv(S);
			d1.elem(Iint1) = inv_S * sign(c.rows(Iint1));



		} else if (gamma == gamma_wM_moin) {
			//	cout << "gamma_wM_moin"<<endl;

			Etat_var_w.at(i8) = 1;
			deplacer(i8, &IM2, &Iint2);

			w->row(i8) = (M * sign(w->row(i8)));
			a2 -= B.col(i8) * w->row(i8);
			r = y - a1 - a2;
			Hr = H->t() * r;
			Br = B.t() * r;

			inv_BB = inv(BB.submat(Iint2, Iint2));
			R = B.cols(Iint2) * inv_BB;
			F = R * -B.cols(Iint2).t() + eyeN;

			if (!Iint2.is_empty()) {
				c = Hr - (HB.cols(Iint2) * w->rows(Iint2)) - (HH.cols(Iint1)
						* x.rows(Iint1));
			} else {
				c = Hr - HH.cols(Iint1) * x.rows(Iint1);
			}

			// % update direction Iint
			d1 = zeros(Q1);
			d2 = zeros(Q2);

			S = H->cols(Iint1).t() * F * H->cols(Iint1);
			inv_S = inv(S);
			d1.elem(Iint1) = inv_S * sign(c.rows(Iint1));

		}

		 t1=clock();
			 T_test +=(float)(t1-t0)/CLOCKS_PER_SEC;

	}//for

	return x;

}


vec homotopy_l2l0(mat *BB, vec *By, double Tau) {
	//clock_t t, t2;
	double temps = 0;
	mat B;
	double lambda, gamma, gamma1, gamma2, gamma3, nx;
	uword i1, i2, i3;
	int N = By->n_rows;
	//cout << "N= " << N <<endl;
	int niter = 1000;

	double eps = exp10(-10);

	vec x = zeros(N);
	vec d = zeros(N);
	vec c, v, w;
	uvec ws;

	//intialization
	uword I0 = (abs(*By)).index_max();
	lambda = abs(By->at(I0));
	//cout << "lambda " << lambda << endl;

	vec JJ = ones(N);
	uvec J, I = regspace<uvec> (0, 0);
	I(0) = I0;
	JJ.at(I0) = 0;

	for (int it = 0; it < niter; it++) {

		if (it > itmax)
			itmax = it;

		J = find(JJ == 1);
		c = *By - (BB->rows(I)).t() * x.rows(I);

		d.reset();

		B = inv(BB->submat(I, I));
		d = B * sign(c.rows(I));

		v = BB->submat(J, I) * d;

		w = (lambda - c.rows(J)) / (1 - v);

		ws = find(w >= eps);
		if (!ws.empty()) {
			gamma1 = min(w.elem(ws));
			i1 = ((uvec) J.elem(find(w == gamma1))).at(0);
		} else {
			gamma1 = exp10(10);
		}

		w = (lambda + c.rows(J)) / (1 + v);
		ws = find(w >= eps);
		if (!ws.empty()) {
			gamma2 = min(w.elem(ws));
			i2 = ((uvec) J.elem(find(w == gamma2))).at(0);
		} else {
			gamma2 = exp10(10);
		}

		w = -x.rows(I) / d;
		ws = find(w >= eps);
		if (!ws.empty()) {
			gamma3 = min(w.elem(ws));
			i3 = ((uvec) I.elem(find(w == gamma3))).at(0);
		} else {
			gamma3 = exp10(10);
		}

		gamma = min(gamma1, min(gamma2, gamma3));

		if (gamma == 0 || gamma == exp10(10))
			return x;

		x.rows(I) += gamma * d;
		lambda -= gamma;
		nx = norm(x.rows(I), 1);
		if (nx >= Tau) {
			if (nx > Tau) {
				vec sn = sign(x.rows(I) - gamma * d / 2);
				vec gm = (nx - Tau) / (d.t() * sn);
				x.rows(I) -= d * gm;
			}
			if (it > 50)
				cout << "it :" << it << endl;
			return x;
		}

		if (gamma == gamma3) {
			JJ.at(i3) = 1;
			I = find(JJ == 0);

		} else {
			if (gamma == gamma2) {
				JJ.at(i2) = 0;
				I = find(JJ == 0);
				//I.insert_rows(0,1);
				//I(0)=i2;

				//I.t().print("I_i2");
			} else {
				if (gamma == gamma1) {
					JJ.at(i1) = 0;
					I = find(JJ == 0);
					//I.insert_rows(0,1);
					//I(0)=i1;

				}
			}
		}
	}

	return x;
}


vec homotopy_onecol_l2l0 (mat *BB, vec *By, double Tau){

	//clock_t t, t2;
	//double temps=0;

	double lambda,gamma,gamma1,gamma2,gamma3,nx;
	uword i,i1,i2,i3;

	int N= By->n_rows;
	//cout << "N= " << N <<endl;
	int niter=1000;
	int add=1;

	vec x=zeros(N);
	vec d=zeros(N);
	vec c,v,w;
	mat F11inv,F22inv;
	vec v1,u1,u2,u3;
	mat B= mat(1,1);
	uvec ws,pos;
	uvec i0=uvec(1);


	//intialization
	uword I0 = (abs(*By)).index_max();
	lambda = abs(By->at(I0));
	vec JJ=ones(N);
	uvec J,Ip,I=regspace<uvec>(0,0);
	I(0)=I0;
	JJ.at(I0)=0;

	for(int it=0; it<niter ;it++){
		J=find(JJ==1);

		c= *By - (BB->rows(I)).t() *x.rows(I) ;
		d.reset();

		/******************************************/
		//d=solve( ((mat) BB->cols(I)).rows(I), sign(c.rows(I)));
		 if(it>0){

		        if (add == 1){
		        	i0(0)=i;
		        	u1 = (BB->submat(Ip,i0));
		        	u2 = B*u1;
		        	F22inv = 1/(BB->submat(i0,i0) - u1.t()*u2);
		        	u3 = F22inv(0,0)*u2;
		            F11inv = B + (u2*F22inv(0,0))*u2.t();
		            B=F11inv;
		            B.insert_cols(B.n_cols,1,1);
		            B.col(B.n_cols-1)=-u3;
		            B.insert_rows(B.n_rows,1);
		            u3.insert_rows(u3.size(),-F22inv);
		            B.row(B.n_rows-1)=-u3.t();

		           // B=inv(((mat) BB->cols(I)).rows(I));
		        }else{

		        	// find the position in new matrix
		           	pos = find(Ip==i);
		            B.insert_rows(B.n_rows,B.row(pos(0)));
		            B.insert_cols(B.n_cols,B.col(pos(0)));
		            B.shed_col(pos(0));
		            B.shed_row(pos(0));
		            F11inv =B.submat( 0, 0, B.n_rows-2, B.n_cols -2 );
		            F22inv = B(B.n_rows-1, B.n_cols-1);
		            u3 = (vec) -B.submat( 0, B.n_cols-1, B.n_rows-2, B.n_cols-1 );
		            u2 = u3/F22inv(0,0);
		            B = F11inv - u2*(u2.t()*F22inv(0,0));

		        }

		 }else{
			 B= 1/(BB->submat(I,I));
		 }

		 //B.print("B");
		 d=B* sign(c.rows(I));
		/*****************************************/
		v =BB->submat(J,I) *d;
		w = (lambda - c.rows(J)) / (1-v);

		ws= find(w>=exp10(-10)) ;
		if(!ws.empty()){
			gamma1 = min(w.elem( ws ));
			i1= ((uvec)J.elem(find(w==gamma1))).at(0);
		}else{
			gamma1 = exp10(10);
		}


		w = (lambda + c.rows(J)) / (1+v);
		ws= find(w>=exp10(-10)) ;
		if(!ws.empty()){
			gamma2 = min(w.elem( ws ));
			i2=  ((uvec)J.elem(find(w==gamma2))).at(0);
		}else{
			gamma2 = exp10(10);
		}

		w = -x.rows(I) / d ;
		ws= find(w>=exp10(-10)) ;
		if(!ws.empty()){
			gamma3 = min(w.elem( ws ));
			i3=  ((uvec)I.elem(find(w==gamma3))).at(0);
		}else{
			gamma3 = exp10(10);
		}


		gamma=min(gamma1,min(gamma2,gamma3));

		if(gamma==0 || gamma == exp10(10))
			return x;


		x.rows(I)+= gamma*d;
		lambda-= gamma;
		nx=norm(x.rows(I),1);
		if (nx>=Tau){
			if(nx>Tau){
				vec sn= sign(x.rows(I)- gamma*d/2);
				vec gm=(nx-Tau)/ (d.t()*sn);
				x.rows(I)-= d*gm;
			}


			return x;
		}


		Ip=I;

		if(gamma==gamma3){
			JJ.at(i3)=1;
			I.shed_row( ((uvec)find(I==i3))(0) );

			i=i3;
			add=0;

		}else{
			if(gamma==gamma2){
				JJ.at(i2)=0;
				I.insert_rows(I.n_rows,1);
				I(I.n_rows-1)=i2;

				add=1;
				i=i2;

			}else{
				if(gamma==gamma1){
					JJ.at(i1)=0;
					I.insert_rows(I.n_rows,1);
					I(I.n_rows-1)=i1;
					add=1;
					i=i1;

				}
			}
		}
	}

	cout<<"sort fin"<<endl;
	return x;

}


