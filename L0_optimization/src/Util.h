

uword arg_min_p(uvec J, vec e, double *val, int *i1, double eps) {


	*i1 = -1;

	*val = exp10(10);

	if(e.n_elem < 1){
		//cout << "fonction arg_min_p vecteur vide"<<endl;
		return -1;
		}



	for (uword i = 0; i < e.n_elem; i++) {
		if (e(i) >= eps && e(i) < *val) {
			*i1 = i;
			*val = e(i);
		}
	}

	if (*i1 > -1)
		return J.at(*i1);

	return -1;
}



uword arg_max_p(uvec J, vec e, double *val, int *i1, double eps) {
	*i1 = -1;

	*val = -1;

	for (uword i = 0; i < e.n_elem; i++) {
		if (e(i) >= eps && e(i) > *val) {
			*i1 = i;
			*val = e(i);
		}
	}

	if (*i1 > -1)
		return J.at(*i1);

	return 0;
}


void pause_cin() {
	int a;
	cin >> a;
}

uword deplacer(uword u1, uvec *v1, uvec *v2) {
	uword i(0);
	bool trouve = false;

	while (i < v1->n_rows) {
		if (v1->at(i) == u1) {
			trouve = true;
			break;
		}
			i++;
	}

	if (!trouve) {
		cout << "Function deplacer : erreur u1 not in v1" << endl;
		//pause_cin();
		return -1;
	}

	v1->shed_row(i);

	v2->insert_rows(v2->n_rows, 1);
	v2->at(v2->n_rows - 1) = u1;

	return i;
}

uvec vec_add_val(uvec I, uword i1) {
	I.insert_rows(I.n_rows, 1);
	I(I.n_rows - 1) = i1;
	return I;
}

mat inversion_rec_add_P(mat B, vec HPv, mat vPv) {
	mat F11inv, F22inv;
	vec v1, u2, u3;

	//i0(0) = pos;
	//Hv = (BB->submat(Ip, i0));
	u2 = B * HPv;
	F22inv = 1 / (vPv - HPv.t() * u2);
	u3 = F22inv(0, 0) * u2;
	F11inv = B + (u2 * F22inv(0, 0)) * u2.t();
	B = F11inv;
	B.insert_cols(B.n_cols, 1, 1);
	B.col(B.n_cols - 1) = -u3;
	B.insert_rows(B.n_rows, 1);
	u3.insert_rows(u3.size(), -F22inv);
	B.row(B.n_rows - 1) = -u3.t();

	return B;

}

mat inversion_rec_del_P(mat B, uword pos ) {
	mat F11inv, F22inv;
	vec v1, u2, u3;

	if(B.n_elem==1)
	{		B= mat(0,0);

	}else{

				B.insert_rows(B.n_rows, B.row(pos));
				B.insert_cols(B.n_cols, B.col(pos));
				B.shed_col(pos);
				B.shed_row(pos);
				F11inv = B.submat(0, 0, B.n_rows - 2, B.n_cols - 2);
				F22inv = B(B.n_rows - 1, B.n_cols - 1);
				u3 = (vec) -B.submat(0, B.n_cols - 1, B.n_rows - 2,	B.n_cols - 1);
				u2 = u3 / F22inv(0, 0);
				B = F11inv - u2 * (u2.t() * F22inv(0, 0));

	}

return B;

}



int norm_zero(vec x){
	return sum(abs(x)-err_x>0);
}
