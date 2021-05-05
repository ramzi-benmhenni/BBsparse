



bool opt_relaxation_with_cplex =0;
bool opt_relaxation_with_actset = 0;
bool opt_relaxation_with_hom = 0;
bool opt_with_cplex =0;
bool cut;

bool warm_restart=1;
bool non_elag_hom=0;
int branching_rule =1; // 1 -> bmax; 2 ->  score


/**** global_variable ***/
mat A;
vec y;
vec x_opt;



double norm2carre_y;
//int k;

sp_mat eyeN;
mat AA;

double BigM=0;
double born_sup= exp10(10);

int BBNodenum =0;
int BestNodenum=0;

double TimeBBMax= 1000;


long int nbr_iter=0;

int nbr_update_Bsupp=0;
double T_best_Bsupp=0;



double eps_binary = exp10(-8);
double gap = exp10(-8);

double err_x = exp10(-4);
/************************/
float T_homo=0;
float T_test=0;
double current_time = 0;

struct Node
{
	int num;
    double born_inf;
    vector<int> vec_Qv,vec_Q1,vec_Q0;

    vec x_relache;
    vec b_relache;

    vec score;


    bool is_solution=false;
    bool is_feasible=false;

    int branch;
    int qj_parent;


    mat invAAQ1;
    mat invAAQ1_parent;

    uvec ind_Qv_in_parent;
    uvec ind_Qv_0_parent;
    uvec ind_Qv_M_parent;

    mat invAAQv_in_parent;

    uvec ind_Qv_in;
    uvec ind_Qv_0;
    uvec ind_Qv_M;

    mat invAAQv_in;


   // mat BB;
   // mat Bz;
};




