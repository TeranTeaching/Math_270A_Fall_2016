#ifndef _LAGRANGIAN_FORCES_H_
#define _LAGRANGIAN_FORCES_H_

#include <Eigen/Dense>
#include <Eigen/Core>
#include "SymmetricTridiagonal.h"

namespace JIXIE{
template <class T>
class LagrangianForces{
  typedef Eigen::Matrix<T,Eigen::Dynamic,1> TVect;
public:
  LagrangianForces(){}
  virtual T PotentialEnergy(const TVect& x){return (T)0;}
  virtual void AddForce(TVect& force,const TVect& x,T scale=(T)1){}
  virtual void AddForceDerivative(SymmetricTridiagonal<T>& A,const TVect& x,T scale=(T)1){}
  virtual void AddForceDifferential(TVect& df,const TVect& x,const TVect& dx,T scale=(T)1){

  }
};

template <class T>
class ConstitutiveModel{
public:
  ConstitutiveModel(){}
  virtual T EnergyDensity(const T& F){return (T)0;}
  virtual void P(T& P, const T& F){P=(T)0;}
  virtual void dPdF(T& dPdF,const T F){dPdF=(T)0;}
};

template <class T>
class LinearElasticity:public ConstitutiveModel<T>{
public:
  T k;
  LinearElasticity(const T k_input):k(k_input){}
  T EnergyDensity(const T& F){return (T).5*k*(F-(T)1)*(F-(T)1);}
  void P(T& P, const T& F){P=k*(F-(T)1);}
  void dPdF(T& dPdF,const T F){dPdF=k;}
};

template <class T>
class NeoHookean:public ConstitutiveModel<T>{
public:
  T k;
  NeoHookean(const T k_input):k(k_input){}
  T EnergyDensity(const T& F){
    if(F<0){
      std::cout<< "Negative F: F = " << F << std::endl;
      exit(1);
    }
    T logF=log(F);
    return (T).5*k*logF*logF;
  }
  void P(T& P, const T& F){
    if(F<0){
      std::cout<< "Negative F: F = " << F << std::endl;
      exit(1);
    }
    T logF=log(F);
    P=k*logF/F;
  }
  void dPdF(T& dPdF,const T F){
    if(F<0){
      std::cout<< "Negative F: F = " << F << std::endl;
      exit(1);
    }
    dPdF=k/(F*F)*(1-log(F));
  }
};

template <class T>
class FEMHyperelasticity:public LagrangianForces<T>{
  typedef Eigen::Matrix<T,Eigen::Dynamic,1> TVect;
  typedef Eigen::Matrix<T,2,2> TMat2;
  int N;
  T a,dX;
  ConstitutiveModel<T>& cons_model;
  std::vector<int> constrained_nodes;
public:
  FEMHyperelasticity(const T a_input, const T dX_input,const int N_input,ConstitutiveModel<T>& cons_input):N(N_input),
  a(a_input),dX(dX_input),cons_model(cons_input){constrained_nodes.resize(0);}

  T F(const TVect& x,const int e)const{
    return (x(e+1)-x(e))/dX;
  }

  T PotentialEnergy(const TVect& x){
    T PE=(T)0;
    for(int e=0;e<N-1;e++)
      PE+=cons_model.EnergyDensity(F(x,e))*dX;
    return PE;
  }

  void AddForce(TVect& force,const TVect& x,T scale=(T)1){
    for(int e=0;e<N-1;e++){
      T P;cons_model.P(P,F(x,e));
      force(e)-=scale*P;
      force(e+1)+=scale*P;
    }
  }

  void AddForceDerivative(SymmetricTridiagonal<T>& A,const TVect& x,T scale){
    for(int e=0;e<N-1;e++){
      T dPdF;cons_model.dPdF(dPdF,F(x,e));
      TMat2 element_stiffness;
      element_stiffness << -scale*dPdF/dX,scale*dPdF/dX,scale*dPdF/dX,-scale*dPdF/dX;
      for(int i=0;i<2;i++){
        for(int j=i;j<2;j++){
          A(e+i,e+j)+=element_stiffness(i,j);}}
    }
  }

  void AddForceDifferential(TVect& result,const TVect& x,const TVect& dx,T scale){
    SymmetricTridiagonal<T> dfdx(N);
    dfdx.SetToZero();
    for(int e=0;e<N-1;e++){
      T dPdF;cons_model.dPdF(dPdF,F(x,e));
      TMat2 element_stiffness;
      element_stiffness << -scale*dPdF/dX,scale*dPdF/dX,scale*dPdF/dX,-scale*dPdF/dX;
      for(int i=0;i<2;i++){
        for(int j=i;j<2;j++){
          dfdx(e+i,e+j)+=element_stiffness(i,j);}}
    }
    TVect df(N);df=TVect::Zero(N);
    dfdx.Multiply(dx,df);
    result+=scale*df;
  }
};

}
#endif
