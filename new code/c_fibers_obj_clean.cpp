#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <cmath>
#include <chrono>
#include <random>
#include<vector>
#include<iostream>
#include<algorithm>

// Double Typedefs
typedef Eigen::VectorXd Vector;
typedef Eigen::Ref<Vector> RefVector;
typedef Eigen::Ref<const Vector> CstRefVector;
typedef Eigen::MatrixXd Matrix;
typedef Eigen::Ref<Matrix, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>> RefMatrix;

typedef Eigen::Triplet<float> Trip;
typedef Eigen::SparseMatrix<float> SparseM;

typedef Eigen::Triplet<double> Trip_d;
typedef Eigen::SparseMatrix<double> SparseMd;
// Float typedefs
/*
typedef Eigen::VectorXf Vector;
typedef Eigen::Ref<Vector> RefVector;
typedef Eigen::Ref<const Vector> CstRefVector;
typedef Eigen::MatrixXf Matrix;
typedef Eigen::Ref<Matrix, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>> RefMatrix;
*/
typedef Eigen::VectorXi IVector;
typedef Eigen::Ref<IVector> IRefVector;
typedef Eigen::Ref<const IVector> CstIRefVector;


struct SpecialParameter{
  float psi = -1.0;
  float Lx = -1, Ly=-1, Lz=-1;
  float shear = 0.0;
  int NPartsPerBatch;
  int NBatches;
};

enum class geometry{rpy, triply_periodic, single_wall};

class CManyFibers{
    float a, ds, dt, k_bend, M0, kBT, eta, Lp;
    int num_parts;
    //PSEParameters par, par_Mhalf;
    geometry geom;
    SpecialParameter parStar;
    bool clamp; // fibers clamped at one end or not
    Vector T_fix; // Ghost vector for clamped fibers (all use this same 3x1 for now)
    // Solver parameters
    float impl;
    float impl_c;
    float alpha;// = 2.0*impl*M0;
    float scale_00; // = M0/(1.0+alpha);
    float scale_10; // = M0;
    float scale; // scale for impl_hydro block
    bool PC_wall = false; // use wall corrections in PC or not
    bool perm_set = false;
    IVector perm;
    IVector i_perm;
    
    // Fiber configurations
    Vector X_0_h, X_0_p, X_0_m;
    Matrix T_h, U_h, V_h;
    Matrix T_p, U_p, V_p;
    Matrix T_m, U_m, V_m;
    bool parametersSet = false;
    //tbb::global_control tbbControl;
  public:
    // CFibers():
    //   tbbControl(tbb::global_control::max_allowed_parallelism, 3){     
    // }
  
    void setParametersPSE(float psi){
      parStar.psi = psi;
    }
      
    void setParameters(int DomainInt, int Nfib, int NblobPerFib, float a, float ds, float dt, float k_bend, float M0, float impl_c, float kBT, float eta, float Lp, bool clamp, Vector& T_fix){
      // TODO: Put the list of parameters into a structure
      this->a = a;
      this->ds = ds;
      this->dt = dt;
      this->k_bend = k_bend;
      this->M0 = M0;
      this->kBT = kBT;
      this->eta = eta;
      this->Lp = Lp;
      this->clamp = clamp;
      this->T_fix = T_fix;
      this->parametersSet = true;
      // solver parameters
      this->impl_c = impl_c;
      this->impl = (impl_c*(dt*k_bend/(ds*ds*ds)));
      this->alpha = 2.0*impl*M0;
      this->scale_00 = M0/(1.0+alpha);
      this->scale_10 = M0;
      float K_scale = 2.0*ds;
      this->scale = (M0/(1+alpha))/(K_scale);
      int numParts = NblobPerFib*Nfib;
      this->num_parts = numParts;
      
      if(DomainInt == 0){
          geom = geometry::rpy;
          parStar.NBatches = Nfib;
          parStar.NPartsPerBatch = NblobPerFib;
          PC_wall = false;
          std::cout << "using batch hydro (NO INTERFIBER HYDRO)\n";
      }
      else if(DomainInt == 1){
          geom = geometry::rpy;
          parStar.NBatches = 1;
          parStar.NPartsPerBatch = numParts;
          PC_wall = false;
          std::cout << "using batch hydro (NO INTERFIBER HYDRO)\n";
      }
      else if(DomainInt == 2){
          geom = geometry::single_wall;
          parStar.NBatches = Nfib;
          parStar.NPartsPerBatch = NblobPerFib;
          PC_wall = true;
          std::cout << "using batch hydro (NO INTERFIBER HYDRO)\n";
      }
      else if(DomainInt == 3){
          geom = geometry::single_wall;
          parStar.NBatches = 1;
          parStar.NPartsPerBatch = numParts;
          PC_wall = true;
          std::cout << "using batch hydro (NO INTERFIBER HYDRO)\n";
      }
      else if(DomainInt == 4){
          //geom = geometry::doubly_periodic;
          std::cout << "doubly periodic not implemeneted yet\n";
      }
      else if(DomainInt == 5){
          geom = geometry::triply_periodic;
      }
      else{
          std::cout << "using built in rpy implementation\n";
      }
      
      if(DomainInt <= 5){
        // TODO: set up solver parameters

        // libmobility::Parameters par;
        // par.hydrodynamicRadius = {a};
        // par.viscosity = eta;
        // par.temperature = 0.5;
        // par.numberParticles = numParts;
        // parStar.Lx = Lp; parStar.Ly = Lp; parStar.Lz = Lp;
        // if(DomainInt == 3 and parStar.psi == -1.0){
        //     std::cout << "you need to set psi before initializing\n";
      }
      // TODO: set up solver
      //   solver = createSolver(geom, parStar);
      //   solver->initialize(par);
    }
    
    void update_T_fix(Vector& T_fix){
      this->T_fix = T_fix;
    }
    
    std::vector<float> single_fiber_Pos(Matrix& T, Vector& X_0){
      const int N_lk = T.cols();
      Eigen::Vector3d T_j, Sum;
      Sum[0] = X_0[0];
      Sum[1] = X_0[1];
      Sum[2] = X_0[2];
      
      int size = 3*N_lk + 3;
      std::vector<float> pos;
      pos.reserve(size);
      
      pos.push_back(Sum[0]);
      pos.push_back(Sum[1]);
      pos.push_back(Sum[2]);
      
      // Rotate the frame
      for(int j = 0; j < N_lk; ++j){
          // Set components of the frame
          T_j = T.col(j);
          Sum += ds*T_j;
          pos.push_back(Sum[0]);
          pos.push_back(Sum[1]);
          pos.push_back(Sum[2]);
      }
      return pos;
  
    }
    
    
    template<class AMatrix, class BVector>
    std::vector<float> multi_fiber_Pos(AMatrix& T, BVector& X_0){
      const int N_lk = T.cols();
      const int N_fib = T.rows()/3;
      
      int size = N_fib*(3*N_lk + 3);
      std::vector<float> pos;
      std::vector<float> pos_j;
      pos.reserve(size);
      Matrix T_j(3,N_lk);
      Vector X_0_j(3);
      for(int j = 0; j < N_fib; ++j){
          for(int k = 0; k < N_lk; ++k){
              T_j(0,k) = T(3*j+0,k);
              T_j(1,k) = T(3*j+1,k);
              T_j(2,k) = T(3*j+2,k);
              //
              X_0_j(0) = X_0(3*j+0);
              X_0_j(1) = X_0(3*j+1);
              X_0_j(2) = X_0(3*j+2);
          }
          pos_j = single_fiber_Pos(T_j, X_0_j);
          pos.insert( pos.end(), pos_j.begin(), pos_j.end() );
      }
      return pos;
  
    }
    
    template<class AMatrix>
    std::vector<float> end_to_end_distance(AMatrix& T){
      const int N_lk = T.cols();
      const int N_fib = T.rows()/3;
      
      std::vector<float> e_2_e;
      e_2_e.reserve(N_fib);
      for(int j = 0; j < N_fib; ++j){
          Vector T_sum = Vector::Zero(3);
          for(int k = 0; k < N_lk; ++k){
              T_sum(0) += ds*T(3*j+0,k);
              T_sum(1) += ds*T(3*j+1,k);
              T_sum(2) += ds*T(3*j+2,k);
          }
          e_2_e.push_back(T_sum.norm());
      }
      return e_2_e;
    }

private:

};

using namespace pybind11::literals;
namespace py = pybind11;

//TODO: Fill python documentation here
PYBIND11_MODULE(c_fibers_obj, m) {
    m.doc() = "Fibers code";
    py::class_<CManyFibers>(m, "CManyFibers").
      def(py::init()).
      def("setParameters", &CManyFibers::setParameters,
	  "Set parameters for the module").
      def("update_T_fix", &CManyFibers::update_T_fix,
      "update ghost tangent").
    //   def("RHS_and_Midpoint", &CManyFibers::RHS_and_Midpoint,
	//   "Generate the RHS for the solve and the midpoint positions").
      def("multi_fiber_Pos", &CManyFibers::multi_fiber_Pos<RefMatrix&, RefVector& >,
	  "Get the blob positions").
      def("end_to_end_distance", &CManyFibers::end_to_end_distance<RefMatrix&>,
	  "end_to_end_distance");
}
    //   def("frame_rot", &CManyFibers::frame_rot<RefMatrix&, RefVector&, RefMatrix& >,
	//   "Rotate the fibers and frame").
    //   def("apply_A_x_Banded_PC", &CManyFibers::apply_A_x_Banded_PC,
	//   "apply A*PC(X)").
    //   def("apply_Banded_PC", &CManyFibers::apply_Banded_PC,
    //   "apply the PC").
    //   def("Outer_Product_Transpose", &CManyFibers::Outer_Product_Transpose<RefMatrix&, RefVector& >,
	//   "Calculate the left outer product matrix transpose").
    //   def("Outer_Product_Tens", &CManyFibers::Outer_Product_Tens<RefMatrix&, RefVector& >,
	//   "Calculate the left outer product matrix").
    //   def("Stresslet_RFD", &CManyFibers::Stresslet_RFD,
	//   "Calculate the RFD of the stress").
    //   def("Stresslet_Strat", &CManyFibers::Stresslet_Strat,
	//   "Calculate the Strat integral of the stress").
    //   def("Stresslet_Correct", &CManyFibers::Stresslet_Correct,
	//   "Correct the stress").
    //   def("Kill_Variance", &CManyFibers::Kill_Variance,
	//   "Kill_Variance").
    //   def("Stresslet_KsF", &CManyFibers::Stresslet_KsF,
	//   "Calculate the K_s*F");