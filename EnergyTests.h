#ifndef _ENERGY_TESTS_H_
#define _ENERGY_TESTS_H_

#include <Eigen/Dense>
#include <Eigen/Core>
#include <cstdlib>

namespace JIXIE{
template <class T>
class EnergyTest{
public:
  typedef Eigen::Matrix<T,Eigen::Dynamic,1> TVect;
  LagrangianForces<T>& lf;
  std::string output_directory;
  int Nrf;
  EnergyTest(const std::string& output_directory_input,LagrangianForces<T>& lf_input,const int Nrf_input):
  output_directory(output_directory_input),lf(lf_input),Nrf(Nrf_input){}

  void RandomPerturbation(TVect& dx){
    for(int i=0;i<dx.size();i++)
      dx(i)=-(T)1+(T)2*(T)rand()/(T)RAND_MAX;
  }

  virtual void RefinementTest(const TVect& x){
    std::string error_force_filename(output_directory+std::string("/")+std::string("refinement_force_error.dat"));
    std::string error_derivative_filename(output_directory+std::string("/")+std::string("refinement_derivative_error.dat"));
    std::string epsilon_filename(output_directory+std::string("/")+std::string("refinement_epsilon.dat"));

    TVect f=TVect::Zero(x.size());
    TVect delta_x(x.size());
    TVect df(x.size());
    RandomPerturbation(delta_x);
    lf.AddForce(f,x);

    TVect force_error(Nrf),epsilon(Nrf),derivative_error(Nrf);
    for(int rf=0;rf<Nrf;rf++){
      epsilon(rf)=(T).1/(T)(rf+1);
      //energy/force consistency check
      T PE_plus,PE_minus;
      PE_plus=lf.PotentialEnergy(x+(T).5*epsilon(rf)*delta_x);
      PE_minus=lf.PotentialEnergy(x-(T).5*epsilon(rf)*delta_x);
      force_error(rf)=fabs((PE_plus-PE_minus)/epsilon(rf)-f.dot(delta_x));
      //force/derivative consistency check
      df=TVect::Zero(x.size());
      lf.AddForce(df,x+(T).5*epsilon(rf)*delta_x,(T)1/epsilon(rf));
      lf.AddForce(df,x-(T).5*epsilon(rf)*delta_x,-(T)1/epsilon(rf));
      lf.AddForceDifferential(df,x,delta_x,(T)-1);
      derivative_error(rf)=df.norm();
    }
    FILE_IO::Write_DAT_File(error_force_filename,force_error);
    FILE_IO::Write_DAT_File(error_derivative_filename,derivative_error);
    FILE_IO::Write_DAT_File(epsilon_filename,epsilon);
  }

};
}
#endif
