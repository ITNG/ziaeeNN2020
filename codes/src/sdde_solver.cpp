
#include "sdde_solver.h"

/*------------------------------------------------------------*/
void SDDE::set_history(const std::vector<double> &hi)
{
    for (int i = 0; i < (nstart + 1); i++)
        pt_t_ar.push_back(-(nstart - i) / (double)nstart * maxdelay);

    // set pt_x_Var to zero for constant history
    for (int i = 0; i < N; i++)
    {
        pt_x_Var.push_back(dim1());
        for (int j = 0; j < (nstart + 1); j++)
        {
            pt_x_Var[i].push_back(0);
        }
    }

    // p_x_ar: N x nstart
    for (int i = 0; i < N; i++)
    {
        pt_x_ar.push_back(dim1());
        for (int j = 0; j < (nstart + 1); j++)
        {
            if (j == 0)
            {
                pt_x_ar[i].push_back(hi[i]);
            }
            else
                pt_x_ar[i].push_back(pt_x_ar[i][0]);
        }
    }
}
/*------------------------------------------------------------*/
void SDDE::set_params(
    double itinitial,
    double itfinal,
    double iPARk,
    double imaxdelay,
    double idtmax)
{
    tinitial = itinitial;
    tfinal = itfinal;
    PARk = iPARk;
    K_over_N = iPARk / (N + 0.0);
    maxdelay = imaxdelay;

    dtmax = idtmax;
    PARD = 0.0;
    dtmin = 1.e-4;
    dt0 = 1.e-2;
    RelTol = 1.e-5;
    AbsTol = 1.e-8;
    nstart = N+1;
    MaxIter = 10000000;
}
/*------------------------------------------------------------*/
void SDDE::set_matrices(
    const Eigen::MatrixXd &iAdj,
    const Eigen::MatrixXd &iDij,
    const std::vector<std::vector<int>> &adjlist)
{
    Cij = iAdj;
    Dij = iDij;
    adjlistC = adjlist;
}
/*------------------------------------------------------------*/
void SDDE::set_initial_frequencies(const std::vector<double> &w)
{
    pt_w = w;
}
/*------------------------------------------------------------*/
dim1 SDDE::KuramotoK1(const double t, dim2I &plag, int SIM_n)
{
    double sum = 0;
    dim1 k1(N);

    for (int i = 0; i < N; i++)
    {
        sum = 0.0;
        for (size_t j = 0; j < adjlistC[i].size(); j++)
        {
            int k = adjlistC[i][j];
            double d = Dij(i, k);
            double c = Cij(i, k);
            if (d > dtmax)
            {
                sum += c * sin(interp_x(t - d, plag[i][k], k) - pt_x_ar[i][SIM_n]);
            }
            else
            {
                sum += c * sin(pt_x_ar[k][SIM_n] - pt_x_ar[i][SIM_n]);
                // std::cout << "false" << std::endl;
            }
        }
        k1[i] = pt_w[i] + K_over_N * sum;
    }


    return k1;
}
/*------------------------------------------------------------*/
dim1 SDDE::KuramotoK2(const double t, dim2I &plag, int SIM_n, const double dt, const dim1 &k1)
{
    double sum = 0;
    dim1 k2(N);

    for (int i = 0; i < N; i++)
    {
        sum = 0.0;
        for (size_t j = 0; j < adjlistC[i].size(); j++)
        {
            int k = adjlistC[i][j];
            double d = Dij(i, k);
            double c = Cij(i, k);
            if (d > dtmax)
                sum += c * sin(interp_x(t - d, plag[i][k], k) - (pt_x_ar[i][SIM_n] + 0.5 * dt * k1[i]));
            else
                sum += c * sin((pt_x_ar[k][SIM_n] + 0.5 * dt * k1[k]) - (pt_x_ar[i][SIM_n] + 0.5 * dt * k1[i]));
        }
        k2[i] = pt_w[i] + K_over_N * sum;
    }

    return k2;
}
/*------------------------------------------------------------*/
dim1 SDDE::KuramotoK3(const double t, dim2I &plag, int SIM_n, const double dt, const dim1 &k2)
{
    double sum = 0;
    dim1 k3(N);

    for (int i = 0; i < N; i++)
    {
        sum = 0.0;
        for (size_t j = 0; j < adjlistC[i].size(); j++)
        {
            int k = adjlistC[i][j];
            double d = Dij(i, k);
            double c = Cij(i, k);
            if (d > dtmax)
                sum += c * sin(interp_x(t - d, plag[i][k], k) - (pt_x_ar[i][SIM_n] + 0.75 * dt * k2[i]));
            else
                sum += c * sin((pt_x_ar[k][SIM_n] + 0.75 * dt * k2[k]) - (pt_x_ar[i][SIM_n] + 0.75 * dt * k2[i]));
        }
        k3[i] = pt_w[i] + K_over_N * sum;
    }

    return k3;
}
/*------------------------------------------------------------*/
dim1 SDDE::KuramotoK4(const double t, dim2I &plag, int SIM_n, const double dt, const dim1 &x)
{
    double sum = 0;
    dim1 k4(N);

    for (int i = 0; i < N; i++)
    {
        sum = 0.0;
        for (size_t j = 0; j < adjlistC[i].size(); j++)
        {
            int k = adjlistC[i][j];
            double d = Dij(i, k);
            double c = Cij(i, k);
            if (d > dtmax)
                sum += c * sin(interp_x(t - d, plag[i][k], k) - x[i]);
            else
                sum += c * sin(x[k] - x[i]);
        }
        k4[i] = pt_w[i] + K_over_N * sum;
    }

    return k4;
}
/*------------------------------------------------------------*/
void SDDE::integrate()
{
    using namespace std;

    // default_random_engine generator(seed);
    // normal_distribution<double> distribution (0.0,1.0);
    long unsigned int NumberOfMinSteps = 0;
    double thresh = AbsTol / RelTol;
    double RelErr;
    double dt = dt0;
    int NumOfDiscont = 4;
    double discont[4] = {maxdelay, 2 * maxdelay, 3 * maxdelay, tfinal};
    bool TakingMinStep{false};
    long unsigned int SIM_n;
    // long unsigned int plag;
    double t = 0.0;
    int nextdsc = 0;
    bool hitdsc{false};
    double dist;
    SIM_n = nstart;
    dim2I plag(N, std::vector<long unsigned int>(SIM_n));
    // Mdim2I plag = Mdim2I::Zero(N, SIM_n);

    dim1 k1(N), k2(N), k3(N), k4(N);
    dim1 TEMPx(N), ERRORx(N);

    k1 = KuramotoK1(t, plag, SIM_n);

    while ((t <= tfinal || hitdsc) && SIM_n - nstart <= MaxIter)
    {
        t += dt * 0.5;
        k2 = KuramotoK2(t, plag, SIM_n, dt, k1);

        t += dt * 0.25;
        k3 = KuramotoK3(t, plag, SIM_n, dt, k2);

        for (int i = 0; i < N; i++)
            TEMPx[i] = pt_x_ar[i][SIM_n] + dt * 1.0 / 9.0 * (2.0 * k1[i] + 3.0 * k2[i] + 4.0 * k3[i]);

        t += dt * 0.25;
        k4 = KuramotoK4(t, plag, SIM_n, dt, TEMPx);

        for (int i = 0; i < N; i++)
            ERRORx[i] = dt / 72.0 * (-5.0 * k1[i] + 6.0 * k2[i] + 8.0 * k3[i] - 9.0 * k4[i]);

        for (int i = 0; i < N; i++)
            ERRORx[i] = ERRORx[i] / MAX(MAX(fabs(pt_x_ar[i][SIM_n]), fabs(TEMPx[i])), thresh);

        RelErr = 0.0;
        for (int i = 0; i < N; i++)
            RelErr = MAX(RelErr, fabs(ERRORx[i]));

        //decide wether  the step is accepted or not
        if (RelErr <= RelTol || TakingMinStep)
        { //if it is accepted
            for (int i = 0; i < N; i++)
            {
                // pt_x_ar[i].push_back(TEMPx[i]+sqrt(dt)*(PARD*distribution(generator)));
                pt_x_ar[i].push_back(TEMPx[i]); //!!!! noise term removed
                pt_x_Var[i].push_back(k4[i]);
            }
            k1 = k4;
            pt_t_ar.push_back(t);
            SIM_n++;

            if (SIM_n - nstart > MaxIter)
            {
                cout << "Warning: MaxIter reached! EndTime of simulation " << t << endl;
                exit(1);
            }
            dt = dt * MAX(0.5, 0.8 * pow(RelTol / RelErr, 0.3333333333333));
            //hit discontinuities
            hitdsc = false;
            if (nextdsc < NumOfDiscont)
            {
                dist = discont[nextdsc] - t;
                if (dist <= MIN(1.1 * dt, dtmax))
                {
                    dt = dist;
                    nextdsc += 1;
                    hitdsc = true;
                }
                else if (dist <= 2 * dt)
                {
                    dt = 0.5 * dist;
                }
            }
        }
        else
        { //if the step not accepted, half the time step
            t -= dt;
            dt = 0.5 * dt;
        }

        if (dt < dtmin && !hitdsc)
        {
            NumberOfMinSteps++;
            TakingMinStep = true;
            dt = dtmin;
        }
        else
        {
            TakingMinStep = false;
        }
        // stop time step enlarging too much
        if (dt > dtmax)
            dt = dtmax;
    } // end of while loop
}
/*------------------------------------------------------------*/
double SDDE::interp_x(double t, long unsigned int &n0, int i)
{
    assert(n0 >= 0);
    assert(n0 < pt_t_ar.size());
    // std::cout << t << "\t" << pt_t_ar[0]<< std::endl;
    assert(t >= pt_t_ar[0]);
    while (n0 < pt_t_ar.size() - 1 && pt_t_ar[n0 + 1] < t)
        n0++;
    while (pt_t_ar[n0] > t)
        n0--;
    return hermite_x(t, pt_t_ar[n0], pt_x_ar[i][n0], pt_x_Var[i][n0],
                     pt_t_ar[n0 + 1], pt_x_ar[i][n0 + 1], pt_x_Var[i][n0 + 1]);
}
/*------------------------------------------------------------*/
double SDDE::hermite_x(double t, double tn, double Xn, double Vn,
                       double tnp1, double Xnp1, double Vnp1)
{
    double h = tnp1 - tn;
    double s = (t - tn) / h;
    double sq1 = SQUARE(s - 1.0);
    double sq2 = SQUARE(s);
    return (1.0 + 2.0 * s) * sq1 * Xn + (3.0 - 2.0 * s) * sq2 * Xnp1 + h * s * sq1 * Vn + h * (s - 1) * sq2 * Vnp1;
}
/*------------------------------------------------------------*/
double SDDE::interpolate(const std::vector<double> &X,
                         const std::vector<double> &Y,
                         const double xnew, int &n0,
                         const std::string kind)
{
    int len = X.size();
    if (xnew > X[len - 1] || xnew < X[0])
    {
        std::cerr << "warning : data out of interval [X_min, X_max]"
                  << "\n";
        std::cout << X[0] << " < " << xnew << " < " << X[len - 1] << "\n";
        std::cout << "lensght t:" << X.size() << " length y: " << Y.size() << "\n";
        std::cout << "pt_x_ar shape is " << pt_x_ar.size() << " " << pt_x_ar[0].size() << "\n";
    }
    /* find left end of interval for interpolation */
    while (X[n0 + 1] < xnew)
        n0++;
    while (X[n0] > xnew)
        n0--;

    /* linear interpolation*/
    if (kind == "linear")
    {
        std::vector<double> xp(2);
        std::vector<double> yp(2);

        for (int i = 0; i < 2; i++)
        {
            xp[i] = X[n0 + i];
            yp[i] = Y[n0 + i];
        }
        double dydx = (yp[1] - yp[0]) / (xp[1] - xp[0]);
        return yp[0] + dydx * (xnew - xp[0]);
    }

    /* quadratic interpolation */
    if (kind == "quadratic")
    {
        std::vector<double> xp(3);
        std::vector<double> yp(3);
        std::vector<double> L(3);

        if (n0 < N - 2)
        {
            for (int i = 0; i < 3; i++)
            {
                xp[i] = X[n0 + i];
                yp[i] = Y[n0 + i];
            }
        }
        else
        {
            for (int i = 0; i < 3; i++)
            {
                xp[i] = X[n0 - 1 + i];
                yp[i] = Y[n0 - 1 + i];
            }
        }
        for (int k = 0; k < 3; k++)
        {
            double term_N = 1.0;
            double term_D = 1.0;
            for (int j = 0; j < 3; j++)
            {
                if (j != k)
                {
                    term_N *= (xnew - xp[j]);
                    term_D *= (xp[k] - xp[j]);
                }
            }
            L[k] = term_N / term_D;
        }
        double Sum = 0.0;
        for (int i = 0; i < 3; i++)
            Sum += L[i] * yp[i];

        return Sum;
    }
}
/*------------------------------------------------------------*/
double SDDE::order_parameter(const dim1 &x)
{
    int n = x.size();
    double real_R = 0.;
    double imag_R = 0.;

    for (int i = 0; i < n; i++)
    {
        real_R += cos(x[i]);
        imag_R += sin(x[i]);
    }
    real_R /= (double)n;
    imag_R /= (double)n;

    return sqrt(real_R * real_R + imag_R * imag_R);
}
/*------------------------------------------------------------*/

dim1 SDDE::order_parameter_array(const std::vector<int> &nodes)
{
    // xVec[N x nsteps] for each time nstep

    int rows = nodes.size();
    int colms = pt_x_ar[0].size();

    assert(rows > 1);

    dim1 x(rows);
    dim1 order(colms);

    for (int j = 0; j < colms; j++)
    {
        for (int i = 0; i < rows; i++)
        {
            int k = nodes[i];
            x[i] = pt_x_ar[k][j];
        }
        order[j] = order_parameter(x);
    }
    return order;
}
/*------------------------------------------------------------*/
dim1 SDDE::interpolate_order_parameter(
    const std::vector<int> &nodes,
    const double ti,
    const double tf,
    const double dt_intp,
    const std::string kind)
{
    dim1 order = order_parameter_array(nodes);

    int n_new = (int)((tf - ti) / (dt_intp + 0.0));
    // int n = pt_t_ar.size();
    pt_t_intp = linspace(ti, tf, n_new);
    dim1 order_intp(pt_t_intp.size());

    int n0 = 0;
    for (int k = 0; k < n_new; k++)
        order_intp[k] = interpolate(pt_t_ar, order, pt_t_intp[k], n0, kind);

    return order_intp;
}
/*------------------------------------------------------------*/
dim1 SDDE::get_times(const std::string kind)
{
    if (kind != "interpolated")
        return pt_t_ar;
    else
        return pt_t_intp;
}
/*------------------------------------------------------------*/
Eigen::MatrixXd SDDE::correlation(const dim1 &x)
{
    /* Calculate Kuramoto correlation*/
int n = x.size();
    Eigen::MatrixXd cor(n, n);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            cor(i, j) = cos(x[j] - x[i]);

    return cor;
}
/*------------------------------------------------------------*/
Eigen::MatrixXd SDDE::get_correlation()
{
    dim1 x(N);
    int col = pt_x_ar[0].size();

    for (int i = 0; i < N; i++)
        x[i] = pt_x_ar[i][col - 1];
    Eigen::MatrixXd cor = correlation(x);
    return cor;
}

/*------------------------------------------------------------*/
dim1 SDDE::linspace(double a, double b, int n)
{
    std::vector<double> arr;
    double step = (b - a) / (n);

    for (int i = 0; i < n; i++)
    {
        arr.push_back(a + i * step);
    }

    return arr;
}
/*------------------------------------------------------------*/
void SDDE::mean_std(const dim1 &vec, const int id, double &R, double &sR)
{
    std::vector<double> newVec(vec.begin() + id, vec.end());
    int n = newVec.size();
    double m = mean(newVec, 0);
    double s = 0;

    for (int i = 0; i < n; i++)
    {
        double tmp = newVec[i] - m;
        s += tmp * tmp;
    }
    s /= (n + 0.0);

    sR = sqrt(s);
    R = m;
}
/*------------------------------------------------------------*/