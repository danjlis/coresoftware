//CosmicSpray
#include "CosmicSpray.h"

#include <g4main/PHG4InEvent.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4Utils.h>
#include <g4main/PHG4ParticleGeneratorBase.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>      // for PHDataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <iostream>
#include <TSystem.h>
#include "TROOT.h"
#include "TF2.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>  // for operator<<, endl, basic_ostream
#include <memory>    // for allocator_traits<>::value_type

#include <vector>  // for vector, vector<>::const_iterator

double CosmicSpray::_o_z_max;
double CosmicSpray::_o_z_min;
double CosmicSpray::_o_x_max;
double CosmicSpray::_o_x_min;
double CosmicSpray::_mean_x;
double CosmicSpray::_x_rad;
double CosmicSpray::_z_rad;
double CosmicSpray::_y_fix;
double CosmicSpray::_y_top_fix;
double CosmicSpray::_y_bot_fix;
double CosmicSpray::_x_max;
double CosmicSpray::_x_min;
double CosmicSpray::_z_max;
double CosmicSpray::_z_min;
double CosmicSpray::_offset;
double CosmicSpray::_gun_e;


class PHCompositeNode;
class PHG4Particle;
class PHG4ParticleGeneratorBase;

CosmicSpray::CosmicSpray(const std::string &name = "COSMICS", const double &gun_e = 4., const double &y_top_fix = 253., const double &y_bot_fix = 43., const double &z_max = 304.91, const double &x_max = 264.71, const double &x_min = 183.3, const int &debug = 1)
  : PHG4ParticleGeneratorBase(name)
{
  _gun_e = gun_e;
  _z_max = 304.91;
  _z_min = -1*z_max;
  _x_max = x_max;
  _x_min = x_min;
  _mean_x = (_x_max + _x_min)/2;
  _x_rad = (_x_max - _x_min)/2;
  _z_rad = 0;
  _debug = debug;
  _R_earth = 637100000;
  _d_earth = 1000000;
  _R_det = _x_max;
  _f_vertex = new TF1("_f_vertex","0.0007*cos(63*x/48) + 0.0007", -1*TMath::Pi()/2, TMath::Pi()/2);
  _f_momentum = new TF2("_f_momentum","1 - TMath::Power(R/(R+d),2)*(1 - TMath::Power(cos(var[0])*cos(var[1]),2))",-1*TMath::Pi()/2, TMath::Pi()/2, -1*TMath::Pi()/2, TMath::Pi()/2)
  return;
}

void CosmicSpray::set_energy(const double &e){
  _gun_e = e;
}
void CosmicSpray::add_particle(const std::string &name, const unsigned int num)
{
  _particle_names.clear();
  _particle_names.push_back(std::make_pair(name, num));
  return;
}

int CosmicSpray::process_event(PHCompositeNode *topNode)
{
  // set_vertex
  if(_debug) std::cout<<"Processing Event"<<std::endl;
  std::string pdgname;
  for (unsigned int i = 0; i < _particle_names.size(); ++i){
    pdgname = _particle_names[i].first;
    if(_debug) std::cout<<"Particle added: "<<pdgname << std::endl;
  }
  gRandom->SetSeed(0);
  int pdgcode = get_pdgcode(pdgname);
  int trackid = 0;
  double gun_t = 0.0;
  double gun_x, gun_y, gun_z;
  double gun_px, gun_py, gun_pz;
  double y_diff;
  TRandom *randomGen = new TRandom();
  
  randomGen->SetSeed();
  gun_z = randomGen->Uniform(_z_min, _z_max);
  double theta_ran;
  _f_vertex->GetRandom(theta_ran);
  gun_x = _x_max * cos(theta_ran);
  gun_y = _x_max * sin(theta_ran);

  double tr, pr;
  _f_momentum->SetRange(-1*(TMath::Pi()/2 - theta_ran), TMath::Pi());
  _f_momentum->GetRandom2(tr, pr);
  
  if(_debug) std::cout<<"Vertex z and t: "<<gun_z<<" / " <<gun_x<<" / "<<theta_ran<<std::endl;
  
  if (fabs(gun_z) <= _z_max ) {
    gun_y = _y_top_fix * sin(theta_ran);
    
  }    
  else {
    gun_y = _y_top_fix * sin(theta_ran) - (_offset/( _y_top_fix * sin(theta_ran) ))*(fabs(gun_z) - _z_max) ;
    if (gun_y < 0) continue;
  }
  y_diff = gun_y - _y_bot_fix;
  
  if(_debug) std::cout<<"Vertex: "<<gun_x<<" / " <<gun_y<<" / "<<gun_z<<std::endl;
  
  double tr, pr;
  tr = randomGen->Uniform()*TMath::Pi()/2;
  if(_debug) std::cout<<"E/T: "<<_gun_e<<" / " <<tr<<std::endl;
  
  pr = (randomGen->Uniform() - 0.5) * TMath::Pi() * 2;
  double py = -1 * _gun_e * cos(tr);
  double px = _gun_e*sin(tr)*cos(pr);
  double pz = _gun_e*sin(tr)*sin(pr);
  double ddx =  y_diff *sin(tr)*cos(pr);
  double ddz = y_diff *sin(tr)*sin(pr);
  
  if(_debug) std::cout<<"dx/dy/dz: "<<ddx<<" / " <<y_diff<<" / "<<ddz<<std::endl;
  
  if (((gun_x + ddx) > _x_max) || ((gun_x + ddx) < _x_min) || ((gun_z + ddz) > _z_max) || ((gun_z - ddz) < _z_min)) continue;
  
  gun_px = px;
  gun_py = py;
  gun_pz = pz;
  
  TVector3 gun_p(gun_px, gun_py, gun_pz);
  TVector3 dir = gun_p.Unit();
  //if(_debug) std::cout<<"Dir: "<<dx<<" / " <<dy<< " / " <<dz<<std::endl;
  if(_debug) std::cout<<"Norm Dir: "<<dir.X()<<" / " <<dir.Y()<<" / "<<dir.Z()<<std::endl;

  if(_debug) std::cout<<"Momentum: "<<gun_px<<" / "<<gun_py<<" / "<<gun_pz<<std::endl;
  if(_debug)std::cout<<"total mom: "<<_gun_e<<std::endl;
  if(_debug)std::cout<<"Before adding vertex"<<std::endl;
  _InEvent = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");

  int vtxindex = _InEvent->AddVtx(gun_x, gun_y ,gun_z,gun_t);
  if(_debug)std::cout<<"After adding vertex"<<std::endl;

  PHG4Particle *particle = new PHG4Particlev2();
  particle->set_track_id(trackid);
  if(_debug)std::cout<<"track_id: "<<trackid<<std::endl;
  particle->set_vtx_id(vtxindex);
  if(_debug)std::cout<<"vtxindex: "<<vtxindex<<std::endl;
  particle->set_parent_id(0);
  if(_debug)std::cout<<"parent_id "<<std::endl;
  particle->set_name(pdgname);
  if(_debug)std::cout<<"pdgname: "<<pdgname<<std::endl;
  particle->set_pid(pdgcode);
  if(_debug)std::cout<<"pdgcode: "<<pdgcode<<std::endl;
  particle->set_px(gun_px);
  if(_debug)std::cout<<"px"<<std::endl;
  particle->set_py(gun_py);
  if(_debug)std::cout<<"py"<<std::endl;
  particle->set_pz(gun_pz);
  if(_debug)std::cout<<"pz"<<std::endl;
  particle->set_e(_gun_e);
  if(_debug)std::cout<<"ene"<<std::endl;

  _InEvent->AddParticle(vtxindex, particle);
  if(_debug)std::cout<<"left COSMICS"<<std::endl;

  return 0;
}
/*
double CosmicSpray::fun_energy_theta(double *val, double *par){

  double energy = val[0];
  double theta = val[1];

  double E_0 = 4.03;
  double n = 3.08;
  double R = 174;

  double D = sqrt(R*R*cos(theta)*cos(theta) + 2*R + 1) - R*cos(theta);

  double result = TMath::Power(E_0 + energy,-1*n)*TMath::Power(D,-1*n +1);
  return result;
}
*/
