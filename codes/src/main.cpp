#include "sdde_solver.h"
#include "lib.h"

unsigned seed;

int main(int argc, char const *argv[])
{

	if (argc < 3)
	{
		std::cerr << "input error \n";
		exit(2);
	}

	/*-----simulation parameters---------------------------------------------*/
	const int N = atoi(argv[1]);
	const double t_trans = atof(argv[2]);
	const double t_sim = atof(argv[3]);
	const double g = atof(argv[4]);
	const double mu = atof(argv[5]);
	const string cname = argv[6];
	const string dname = argv[7];
	const int ncluster = atoi(argv[8]);
	const int num_sim = atoi(argv[9]);
	const bool cor_to_file = atoi(argv[10]);

	const double tinitial = 0.0;
	const double dt = 0.05;
	const double dtmax = 0.1;
	/*-----------------------------------------------------------------------*/

	seed = 1234;
	double wtime = get_wall_time();

	if (num_sim > 1)
	{
		INITIALIZE_RANDOM_CLOCK(seed);
	}
	else
	{
		INITIALIZE_RANDOM_F(seed);
	}

	int index_transition = (int)(t_trans / dt);
	dim1 hi(N);
	dim1 r_global;
	dim1 R_local(ncluster);
	dim1 R_local_2(2);
	// const double sigma = 0.0;

	/*-----------------------------------------------------------------------*/
	const string subname = to_string(g) + "-" + to_string(mu);
	const string ofname = "../data/text/R-" + subname + ".txt";

	FILE *ofile;
	ofile = fopen(ofname.c_str(), "a");
	if (!fileExists(ofname))
	{
		cout << "output file for r did not open correctly \n!";
		exit(EXIT_FAILURE);
	}

	vector<int> clst = {N};
	vector<vector<int>> local_nodes_l1 = read_nodes_of_each_cluster(
		"networks/communities_"+to_string(N)+"_l1.txt", ncluster);

	vector<vector<int>> local_nodes_l2 = read_nodes_of_each_cluster(
		"networks/communities_"+to_string(N)+"_l2.txt", 2);

	vector<vector<int>> global_nodes = nodes_of_each_cluster(clst);

	Eigen::MatrixXd Cij = read_matrix(N, "networks/" + cname + ".txt");
	Eigen::MatrixXd Dij = read_matrix(N, "networks/" + dname + ".txt");

	vector<vector<int>> adjlistC = adjmat_to_adjlist(Cij);
	const double maxdelay = Dij.maxCoeff();

	for (int ens = 0; ens < num_sim; ens++)
	{

		printf("g = %10.3f, omega = %10.3f, sim = %5d \n", g, mu, ens);

		for (int ii = 0; ii < N; ii++)
			hi[ii] = RANDOM * 2 * M_PI - M_PI;

		SDDE dde(N);
		dde.set_params(tinitial, t_sim, g, maxdelay, dtmax);
		dde.set_matrices(Cij, Dij, adjlistC);
		dde.set_history(hi);
		dim1 w(N, mu);
		dde.set_initial_frequencies(w);
		dde.integrate();

		if (cor_to_file)
		{
			Eigen::MatrixXd cor = dde.get_correlation();
			const string cfname = "../data/text/c-" + subname +
								  "-" + to_string(ens);
			write_matrix_to_file(cfname, cor);
		}

		r_global = dde.interpolate_order_parameter(
			global_nodes[0], 1.0, t_sim - 5, dt, "linear");
		double R_global = mean(r_global, index_transition);

		for (int nl = 0; nl < ncluster; nl++)
		{
			dim1 r = dde.interpolate_order_parameter(
				local_nodes_l1[nl], 1.0, t_sim - 5, dt, "linear");
			R_local[nl] = mean(r, index_transition);
		}
		for (int nl = 0; nl < 2; nl++)
		{
			dim1 r = dde.interpolate_order_parameter(
				local_nodes_l2[nl], 1.0, t_sim - 5, dt, "linear");
			R_local_2[nl] = mean(r, index_transition);
		}

		fprintf(ofile, "%15.9f %15.9f %15.9f", g, mu, R_global);
		for (int nl = 0; nl < ncluster; nl++)
			fprintf(ofile, "%15.9f", R_local[nl]);
		for (int nl = 0; nl < 2; nl++)
			fprintf(ofile, "%15.9f", R_local_2[nl]);
		fprintf(ofile, "\n");
	}
	fclose(ofile);
	/*-----------------------------------------------------------------------*/
	FREE_RANDOM;
	wtime = get_wall_time() - wtime;
	display_timing(wtime, 0);

	return 0;
}
