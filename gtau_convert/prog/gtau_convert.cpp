//#include <alps/params.hpp>
#include <fstream>
#include <iomanip>
#include <boost/tuple/tuple.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "green_function.h"
#include "fouriertransform.h"

static std::complex<double> I(0.0, 1.0);

typedef boost::tuple<double, std::vector<double>, std::vector<std::complex<double> > > input_t;

double tau_f(double beta, int ntau, int i) { return beta / ntau * (i + 0.5); }

//forward declarations
// save gtau
void print_green_itime(std::ostream &os, const itime_green_function_t &v, const double beta);

itime_green_function_t simple_gtau_convert(matsubara_green_function_t const& gw, double beta, int ntau, std::vector<double> tail)
{
    itime_green_function_t gtau (ntau, gw.nsite(), gw.nflavor(), diagonal);
    if (tail.size() != 3) throw std::logic_error("tail should have 3 elements");
    std::cout << "tail = " << tail[0] << "/iw + " << tail[1] << "/(iw)^2 + " << (tail[2]) << "/(iw)^3" << std::endl;
    //SimpleG0FourierTransformer tr(beta, 0, 0, gw.nflavor(), eps, epssq);
    FourierTransformerFixedTail tr(beta, 1, tail);
    tr.frequency_to_time_ft(gtau, gw);
    return gtau;
}

#ifdef BUILD_PYTHON_MODULE
#include <boost/python.hpp>

namespace py = boost::python;
typedef boost::python::numeric::array pyarray;

pyarray pyconvert (pyarray &python_data, int ntau, pyarray pytail) // wrapper from numpy array to STL vector
{
    typedef std::vector<std::complex<double> > cvector_type;
    const py::tuple &shape = py::extract<py::tuple>(python_data.attr("shape"));
    int ncols = py::extract<int>(shape[0]);
    int nw = py::extract<int>(shape[1]); // Python takes care of exceptions

    if (ncols != 2) throw std::logic_error("need a 2*nfreqs array");

    // extract gw data
    cvector_type data(nw);
    cvector_type grid(nw);

    for (int w=0; w<nw; ++w) { 
        grid[w] = py::extract<std::complex<double> >(python_data[py::make_tuple(0,w)]);
        data[w] = py::extract<std::complex<double> >(python_data[py::make_tuple(1,w)]);
        }

    matsubara_green_function_t gw (nw, 1, 1, diagonal);
    std::copy(data.begin(), data.end(), gw(0,0,0));

    // get beta
    double beta = M_PI*2.0 / std::abs(grid[1] - grid[0]);
    std::cout << "beta = " << beta << std::endl;

    // extract tail
    py::tuple tail_shape = py::extract<py::tuple>(pytail.attr("shape"));
    if (py::extract<int>(tail_shape[0]) != 3) throw std::logic_error("tail should be an array of 3 numbers");
    std::vector<double> tail(3,0.0);
    tail[0] = py::extract<double>(pytail[0]);
    tail[1] = py::extract<double>(pytail[1]);
    tail[2] = py::extract<double>(pytail[2]);


    itime_green_function_t gtau(simple_gtau_convert(gw, beta, ntau, tail));
    std::ofstream out("gtau.dat");
    out.setf( std::ios_base::scientific);
    print_green_itime(out, gtau, beta);
    out.close();

    py::tuple shape_out = py::make_tuple(2,ntau);
    typedef boost::python::list pylist; 
    pylist l1;

    for(unsigned int i=0;i<gtau.ntime();++i){
        double tau = tau_f(beta, gtau.ntime(), i);
        double val = gtau(i,0,0,0);
        //l1.append(py::make_tuple(py::make_tuple(tau),py::make_tuple(val),py::make_tuple(0)));
        l1.append(py::make_tuple(tau, val, 0.0));
        }
    pyarray py_gtau(l1);
    py_gtau.transpose();
    return py_gtau;
}

BOOST_PYTHON_MODULE(pygtau_convert)
{
    using namespace boost::python;
    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
    def("convert", &pyconvert);
}


#else  
// read gw
input_t readin(const char *filename);

int main(int argc, char **argv)
{
    boost::program_options::options_description opts;
    opts.add_options()
        ("diagonal" , po::value<bool>()->default_value(true), "leading asymptotics is 1/w?") 
        ("ntau"     , po::value<int>()->default_value(256), "Number of tau points")
        ("gw"       , po::value<std::string>()->default_value("gw.dat"), "file with gw. 1st column : Matsubara, 2nd : real, 3rd : imaginary")
        ("tail"     , po::value<std::vector<double> >()->multitoken(), "3 tail coefficients")
        ("help,h",          "help");

    po::options_description cmdline_opts;
    cmdline_opts.add(opts);
    po::variables_map p;
    po::store(po::parse_command_line(argc, argv, cmdline_opts), p);
    po::notify(p);
    if (p.count("help")) { std::cout << cmdline_opts << std::endl; exit(0); }

    int ntau = p["ntau"].as<int>();
    std::vector<double> tail; 
    if (p["tail"].empty()) { tail = std::vector<double> (3,0.0); tail[0] = 1.0; }
    else tail = p["tail"].as<std::vector<double> >();

    input_t input = readin(p["gw"].as<std::string>().c_str()); 
    std::vector<double>& grid = boost::get<1>(input);
    std::vector<std::complex<double> > vals = boost::get<2>(input);

    int nw = grid.size();
    double beta = boost::get<0>(input);

    int nflavor = 1;
    int nsite = 1;

    matsubara_green_function_t gw (nw, nsite, nflavor, diagonal);
    std::copy(vals.begin(), vals.end(), gw(0,0,0));

    itime_green_function_t gtau(simple_gtau_convert(gw, beta, ntau, tail));

    std::ofstream out("gtau.dat");
    out.setf( std::ios_base::scientific);
    print_green_itime(out, gtau, beta);
    out.close();

/* // this works w alpscore
    alps::params p(argc, (const char**) argv);

    p.define<bool>("diagonal", false, "leading asymptotics is 1/w?");
    p.define<int>("ntau", 256, "Number of tau points");
    p.define<std::string>("gw", "gw.dat", "file with gw. 1st column : Matsubara, 2nd : real, 3rd : imaginary");
    p.define<std::vector<double> >("tail",  "3 tail coefficients");

    if (p.help_requested(std::cerr)){exit(0);}

    if (!p.exists("tail")) { std::vector<double> temp_tail(3,0.0); temp_tail[0] = 1.0; p["tail"] = temp_tail; }


    //std::cout << "beta = " << beta << std::endl;

    */
}

input_t readin(const char *filename)
{
    double temp=0, temp2=0; 
    unsigned int wn=0;
    std::vector<double> grid;
    std::vector<std::complex<double> > vals;
    std::ifstream input_gw(filename);
    if (!input_gw.is_open()) {std::cerr << "No such file" << std::endl; exit(1);}
    while (!input_gw.eof()){
        input_gw >> temp; // matsubara
        if (input_gw.eof()) break;
        grid.push_back(temp);
              
        input_gw >> temp;
        input_gw >> temp2;
              
        if (input_gw.eof()) break;
        vals.push_back(temp + I*temp2);
        wn++;
        }
    int wn_max=wn;
    std::cout << wn_max << " Matsubaras in " << filename << std::endl;
    double beta = M_PI*2.0 / (grid[1] - grid[0]);
    std::cout << "beta = " << beta << std::endl;

    return boost::make_tuple(beta, grid, vals); 
}

#endif

void print_green_itime(std::ostream &os, const itime_green_function_t &v, const double beta){
  os<<std::setprecision(14);
  if(!v.get_header().empty()){
	  os<<v.get_header()<<std::endl;
  }
  for(unsigned int i=0;i<v.ntime();++i){
    double tau = tau_f(beta, v.ntime(), i);
    os<<tau<<" ";
    if (v.shape()==diagonal){
      for (unsigned int k=0; k<v.nsite(); ++k)
        for(unsigned int f=0;f<v.nflavor();++f)
          os<<v(i,k,k,f)<<" ";
    }else if(v.shape()==blockdiagonal){
      for (unsigned int k=0; k<v.nsite(); k+=2)
        for (unsigned int p1=0; p1<=1; p1++)
          for (unsigned int p2=0; p2<=1; p2++)
            for(unsigned int f=0;f<v.nflavor();++f)
              os<<v(i,k+p1,k+p2,f)<<" ";
    }else if(v.shape()==nondiagonal){
      for (unsigned int s1=0; s1<v.nsite(); ++s1)
        for (unsigned int s2=0; s2<v.nsite(); ++s2)
          for(unsigned int f=0;f<v.nflavor();++f)
            os<<v(i,s1,s2,f)<<" ";
    }
    os<<std::endl;
  }
}

