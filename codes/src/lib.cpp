#include "lib.h"

Eigen::MatrixXd read_matrix(int Node, std::string filename)
{
    using namespace Eigen;

    std::ifstream ifile(filename);

    if (fileExists(filename))
    {
        dim2 Cij(Node, dim1(Node));
        for (int i = 0; i < Node; i++)
        {
            for (int j = 0; j < Node; j++)
            {
                ifile >> Cij[i][j];
            }
        }
        ifile.close();

        Eigen::MatrixXd mat(Node, Node);
        for (int i = 0; i < Node; i++)
            mat.row(i) = Eigen::VectorXd::Map(&Cij[i][0], Cij[i].size());

        return mat;
    }
    else
    {
        std::cerr << "\n file : " << filename << " not found \n";
        exit(2);
    }
}
/*------------------------------------------------------------*/
bool fileExists(const std::string &filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}
/*------------------------------------------------------------*/
dim2 kuramoto_correlation(const dim1 &x)
{
    /* Calculate Kuramoto correlation*/
    int n = x.size();
    dim2 cor(n, dim1(n));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            cor[i][j] = cos(x[j] - x[i]);

    return cor;
}
/*------------------------------------------------------------*/
std::vector<double> arange(const double start, const double end,
                           const double step)
{
    int nstep = round((end - start) / step);
    std::vector<double> arr(nstep);

    for (int i = 0; i < nstep; i++)
        arr[i] = start + i * step;
    return arr;
}
//-----------------------------------------------------------//
double get_wall_time()
{
    /*measure real passed time */
    struct timeval time;
    if (gettimeofday(&time, NULL))
    {
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
//------------------------------------------------------------//
double get_cpu_time()
{
    /*measure cpu passed time*/
    return (double)clock() / CLOCKS_PER_SEC;
}
//------------------------------------------------------------//
void display_timing(double wtime, double cptime)
{
    int wh;      //, ch;
    int wmin;    //, cpmin;
    double wsec; //, csec;
    wh = (int)wtime / 3600;
    // ch = (int)cptime / 3600;
    wmin = ((int)wtime % 3600) / 60;
    // cpmin = ((int)cptime % 3600) / 60;
    wsec = wtime - (3600. * wh + 60. * wmin);
    // csec = cptime - (3600. * ch + 60. * cpmin);
    printf("Wall Time : %d hours and %d minutes and %.4f seconds.\n", wh, wmin, wsec);
    // printf ("CPU  Time : %d hours and %d minutes and %.4f seconds.\n",ch,cpmin,csec);
}
//------------------------------------------------------------//
void print_matrix(const dim2 &A, std::string filename)
{
    using namespace std;
    int row = A.size();
    int col = A[0].size();

    ofstream ofile;
    ofile.open(filename);
    if (ofile.is_open())
    {
        ofile.precision(9);
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                ofile << setw(20) << A[i][j];
            }
            ofile << "\n";
        }
        ofile.close();
    }
    else
    {
        std::cout << "Error opening file to print R \n";
    }
}
//------------------------------------------------------------//
void find_dominant_frequency(const dim1 &Freq, const dim1 &Pxx,
                             double &f, double &p)
{
    int n = Pxx.size();
    double max_v = -100000.0;
    int index = 0;
    for (int i = 1; i < n; i++)
    {
        if (max_v < Pxx[i])
        {
            max_v = Pxx[i];
            index = i;
        }
    }
    p = max_v;
    f = Freq[index];
    cout << f << "\t" << p << endl;
}
//------------------------------------------------------------//

vector<vector<int>> nodes_of_each_cluster(vector<int> &clusters)
{
    std::vector<std::vector<int>> local_nodes;
    int nc = clusters.size();
    local_nodes.resize(nc);

    vector<int> a = {0};
    for (int i = 0; i < nc; i++)
    {
        a.push_back(a[i] + clusters[i]);
        vector<int> c(clusters[i]);
        std::iota(c.begin(), c.end(), a[i]);

        for (int j = 0; j < clusters[i]; j++)
            local_nodes[i].push_back(c[j]);
    }

    return local_nodes;
}
//------------------------------------------------------------//

vector<vector<int>> read_nodes_of_each_cluster(
    const std::string filename,
    const int n)
{
    vector<vector<int>> adjlist;
    adjlist.resize(n);
    std::ifstream file(filename);
    if (file.is_open())
    {
        std::string line;
        int counter = 0;

        while (getline(file, line))
        {
            std::istringstream iss(line);
            int value;
            while (iss >> value)
                adjlist[counter].push_back(value);
            counter++;
        }
    }
    else
    {
        printf("%s file not found! \n", filename.c_str());
        exit(EXIT_FAILURE);
    }
    file.close();

    return adjlist;
}
//------------------------------------------------------------//

std::vector<int> read_from_file(std::string filename)
{
    std::vector<int> numbers;
    std::ifstream inputFile(filename);
    if (inputFile.good())
    {
        int current_number = 0;
        while (inputFile >> current_number)
            numbers.push_back(current_number);

        inputFile.close();
    }

    else
    {
        printf("%s file not found! \n", filename.c_str());
        exit(EXIT_FAILURE);
    }
    return numbers;
}
//------------------------------------------------------------//
std::vector<std::vector<int>> adjmat_to_adjlist(const Eigen::MatrixXd &A)
{
    int row = A.rows();
    int col = A.cols();
    std::vector<std::vector<int>> adjlist;
    adjlist.resize(row);

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            if (std::abs(A(i, j)) > 1.e-6)
                adjlist[i].push_back(j);
        }
    }

    return adjlist;
}

/*------------------------------------------------------------*/
void write_matrix_to_file(
    const string fname, Eigen::MatrixXd &m)
{
    int n = m.rows();
    
    std::ofstream oX((fname + ".bin").c_str(),
                     std::ios::out | std::ios::binary);
    if (!oX)
    {
        cout << "Could not open file for binary output!";
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            oX.write((char *)&m(i, j), sizeof(double));

    oX.close();
}

/*------------------------------------------------------------*/