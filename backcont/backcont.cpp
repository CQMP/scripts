#include<iostream>
#include<complex>
#include<fstream>
#include<boost/program_options.hpp>


enum kernel_type{
  standard,
  anomalous,
  bosonic,
  me_bosonic,
  me_anomalous
};

inline std::complex<double> fermionic_standard_kernel(const double &omega_n, const double omega){
  return 1./(std::complex<double>(0., omega_n)-omega);
}
inline std::complex<double> anomalous_kernel(const double &omega_n, const double omega){
  return -1.0 / (std::complex<double>(0., omega_n) - omega);
}
inline std::complex<double> bosonic_kernel(const double &omega_n, const double omega){
  return 1.0 / (std::complex<double>(0., omega_n) + omega);
}
inline std::complex<double> me_anomalous_kernel(const double &omega_n, const double omega){
  return -omega / (std::complex<double>(0., omega_n) - omega);
}
inline std::complex<double> me_bosonic_kernel(const double &omega_n, const double omega){
  return omega / (std::complex<double>(0., omega_n) + omega);
}

inline std::complex<double> kernel(const double &omega_n, const double omega, const kernel_type k_type){
  if(k_type==standard) return fermionic_standard_kernel(omega_n, omega);
  else if(k_type==anomalous) return anomalous_kernel(omega_n, omega);
  else if(k_type==bosonic) return bosonic_kernel(omega_n, omega);
  else if(k_type==me_bosonic) return me_bosonic_kernel(omega_n, omega);
  else if(k_type==me_anomalous) return me_anomalous_kernel(omega_n, omega);
  else throw std::logic_error("kernel not implemented.");
}

int main(int argc, char**argv){
  namespace po = boost::program_options;
  std::string spectra_filename;
  std::string orig_filename;
  std::string output_omega_filename;
  std::string output_diff_filename;
  std::string output_tau_filename;
  std::string kernel_name;
  int n_matsubara, n_tau;
  kernel_type k_type=standard;
  double beta;
  bool multiply_m1divpi=false;
  
  po::options_description desc("Allowed options");
  desc.add_options()
  ("help", "show this help")
  ("beta", po::value<double>(&beta), "inverse temperature")
  ("n_matsubara", po::value<int>(&n_matsubara)->default_value(-1), "number of matsubara frequencies")
  ("n_tau", po::value<int>(&n_tau)->default_value(20000), "number of imaginary time points")
  ("imag_freq_file", po::value<std::string>(&orig_filename)->default_value("G_omega_av.dat"), "input G(i omega_n) to maxent")
  ("real_freq_file", po::value<std::string>(&spectra_filename)->default_value("spectra.dat"), "output A=-1/pi*ImG(omega) from maxent")
  ("output_freq_file", po::value<std::string>(&output_omega_filename)->default_value("G_omega_back.dat"), "backcontinued output G(omega) with errors")
  ("diff_freq_file", po::value<std::string>(&output_diff_filename)->default_value("G_omega_diff.dat"), "difference to input file")
  ("output_tau_file", po::value<std::string>(&output_tau_filename)->default_value("G_tau_back.dat"), "backcontinued output G(tau) with errors")
  ("kernel", po::value<std::string>(&kernel_name)->default_value("standard"), "kernel type: standard, anomalous, ...")
  ("multiply_m1divpi", "if not specified: scales results by -pi, as required if converting ImG to A. If specified: standard Kramers Kronig (required for Sigma/Anomalous/etc backcont)")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout<<desc;
    return 1;
  }

  //toggle between continuation of G to A and continuation of any other quantity with standard Kramers Kronig relation  
  if (vm.count("multiply_m1divpi")) {
    multiply_m1divpi=true;
  }

  std::ifstream orig_file(orig_filename.c_str());
  if(!orig_file.good()) throw std::invalid_argument("imag freq file: "+orig_filename+" could not be opened. specify with --imag_freq_file");
  std::ifstream spectra_file(spectra_filename.c_str());
  if(!spectra_file.good()) throw std::invalid_argument("real freq file: "+spectra_filename+" could not be opened. specify with --real_freq_file");
  if(!vm.count("beta")) throw std::runtime_error("you need to specify the inverse temperature with --beta.");
  if(vm.count("kernel")){
    if(kernel_name==std::string("standard")){
      k_type=standard;
      std::cout<<"using standard kernel."<<std::endl;
    }else if(kernel_name==std::string("anomalous")){
      k_type=anomalous;
      std::cout<<"using anomalous kernel."<<std::endl;
    }else if(kernel_name==std::string("bosonic")){
      k_type=bosonic;
      std::cout<<"using bosonic kernel."<<std::endl;
    }else if(kernel_name==std::string("me_bosonic")){
      k_type=me_bosonic;
      std::cout<<"using maxent's bosonic kernel."<<std::endl;
    }else if(kernel_name==std::string("me_anomalous")){
      k_type=me_anomalous;
      std::cout<<"using maxent's anomalous kernel."<<std::endl;
    }else{
      throw std::runtime_error("kernel type not recognized.");
    }
  }
  
  std::vector<std::complex<double> > imag_freq_data;
  std::vector<std::complex<double> > imag_freq_error;
  std::vector<double > real_freq_data;
  std::vector<double > real_freq_freq;
  do{
    double dummy, imag_freq_data_real, imag_freq_data_imag, imag_freq_sigma_real, imag_freq_sigma_imag;
    orig_file>>dummy>>imag_freq_data_real>>imag_freq_data_imag>>imag_freq_sigma_real>>imag_freq_sigma_imag>>std::ws;
    imag_freq_data.push_back(std::complex<double>(imag_freq_data_real, imag_freq_data_imag));
    imag_freq_error.push_back(std::complex<double>(imag_freq_sigma_real, imag_freq_sigma_imag));
  }while(orig_file.good());
  do{
    double frequency, value, defaultm;
    spectra_file>>frequency>>value>>defaultm>>std::ws;
    real_freq_data.push_back(value);
    real_freq_freq.push_back(frequency);
  }while(spectra_file.good());
  
  std::cout<<"read in files: "<<imag_freq_data.size()<<" matsubara freqs and "<<real_freq_data.size()<<" real frequency points."<<std::endl;
  if(n_matsubara ==-1) n_matsubara=imag_freq_data.size();
  if(real_freq_data[0]+real_freq_data.back() > 1.e-4) std::cerr<<"problem with spectra: does not go to zero at boundary?";
  std::cout<<real_freq_data[0]<<" "<<real_freq_data.back()<<std::endl;
  
  //back-continue to the imaginary axis
  if(k_type==standard){
    std::vector<double > imag_time_back(n_tau,0.);
    std::ofstream gtau_file(output_tau_filename.c_str());
    gtau_file.precision(14);
    for(int i=0;i<n_tau;++i){
      double tau=i/(double)n_tau*beta;
      imag_time_back[i]=0.;
      for(int w=1;w<real_freq_freq.size()-1;++w){
        double freq =real_freq_freq[w];
        double value=real_freq_data[w];
        double delta=(real_freq_freq[w+1]-real_freq_freq[w-1])/2.;
        double kernel=-std::exp(-freq*tau)/(std::exp(-freq*beta)+1);
        if(!std::isnan(kernel))
          imag_time_back[i]+=kernel*value*delta;
        //std::cout<<freq<<" "<<value<<" "<<delta<<" "<<std::exp(-freq*tau)/(std::exp(-freq*beta)+1)*value*delta<<" "<<imag_time_back[i]<<std::endl;
      }
      double kernel1=-std::exp(-real_freq_freq[0]*tau    )/(std::exp(-real_freq_freq[0]    *beta)+1);
      double kernel2=-std::exp(-real_freq_freq.back()*tau)/(std::exp(-real_freq_freq.back()*beta)+1);
      if(!std::isnan(kernel1))
        imag_time_back[i]+=kernel1*real_freq_data[0]*(real_freq_freq[1]-real_freq_freq[0])/2.;
      if(!std::isnan(kernel2))
        imag_time_back[i]+=kernel2*real_freq_data.back()*(real_freq_freq.back()-real_freq_freq[real_freq_freq.size()-2])/2.;
      if(multiply_m1divpi) imag_time_back[i]*=-1./M_PI;
      gtau_file<<tau<<" "<<imag_time_back[i]<<std::endl;
    }
  }
  
  std::vector<std::complex<double> > imag_freq_data_back(n_matsubara);
  std::ofstream gomega_file(output_omega_filename.c_str());
  gomega_file.precision(14);
  for(int n=0;n<n_matsubara;++n){
    double omega_n;
    if(k_type==standard) omega_n=(2.*n+1)*M_PI/beta;
    else omega_n=(2.*n)*M_PI/beta;
    imag_freq_data_back[n]=0.;
    for(int w=1;w<real_freq_freq.size()-1;++w){
      double freq =real_freq_freq[w];
      double value=real_freq_data[w];
      double delta=(real_freq_freq[w+1]-real_freq_freq[w-1])/2.;
      std::complex<double> kernel_val=kernel(omega_n, freq,k_type);
      imag_freq_data_back[n]+=kernel_val*value*delta;
    }
    std::complex<double> kernel1=kernel(omega_n, real_freq_freq[0],k_type); 
    std::complex<double> kernel2=kernel(omega_n, real_freq_freq.back(),k_type);
    imag_freq_data_back[n]+=kernel1*real_freq_data[0]*(real_freq_freq[1]-real_freq_freq[0])/2.;
    imag_freq_data_back[n]+=kernel2*real_freq_data.back()*(real_freq_freq.back()-real_freq_freq[real_freq_freq.size()-2])/2.;
    if(multiply_m1divpi) imag_freq_data_back[n]*=-1./M_PI;
    gomega_file<<omega_n<<" "<<imag_freq_data_back[n].real()<<" "<<imag_freq_data_back[n].imag()<<std::endl;
  }
  
  std::ofstream gomega_diff_file(output_diff_filename.c_str());
  for(int n=0;n<n_matsubara;++n){
    double diff_real=imag_freq_data[n].real()-imag_freq_data_back[n].real();
    double diff_imag=imag_freq_data[n].imag()-imag_freq_data_back[n].imag();
    gomega_diff_file<<(2.*n+1)*M_PI/beta<<" "<<diff_real<<" "<<diff_imag<<" "<<imag_freq_error[n].real()<<" "<<imag_freq_error[n].imag()<<std::endl;
  }
}
