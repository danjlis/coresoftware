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
#include "TF3.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>  // for operator<<, endl, basic_ostream
#include <memory>    // for allocator_traits<>::value_type

#include <vector>  // for vector, vector<>::const_iterator

double CosmicSpray::_y_top_fix;
double CosmicSpray::_y_bot_fix;
double CosmicSpray::_x_max;
double CosmicSpray::_x_min;
double CosmicSpray::_z_max;
double CosmicSpray::_R_earth;
double CosmicSpray::_d_earth;
double CosmicSpray::_gun_e;


class PHCompositeNode;
class PHG4Particle;
class PHG4ParticleGeneratorBase;


double CosmicSpray::GetGamma0(double &energy){
  double mass = 105; // MeV
  double gamma0 = sqrt(1 + TMath::Power(energy/(mass),2));
  return gamma0;
}

double CosmicSpray::GetBeta(double &gamma){
  double b = sqrt(1-TMath::Power((1/gamma),2));
  return b;
}

double CosmicSpray::EnergyAngularDistribution(double *val, double *par){
  double en = TMath::Power(10, val[1]);
  double th = val[0];
  double p;
  double R = 600.;
  double c = 29979245800; //cm/s
  double tau = 0.0000022; // seconds
  double h_0 = 1000000; //cm
  double g = GetGamma0(en);
  double b = GetBeta(g);
  double ld = c*tau*g*b;
  double norm = 1/((sqrt(R*R*cos(th)*cos(th) + 2*R + 1) - R*cos(th))*ld);
  double f = exp(-h_0*(sqrt(R*R*cos(th)*cos(th) + 2*R + 1) - R*cos(th))/(ld));
//std::cout<<"g: "<<g<<"b: "<<b<<"ld: "<<ld<<"norm: "<<norm<<"f: "<<f<<std::endl;
  p = fabs(norm*f);
  return p;
}

CosmicSpray::CosmicSpray(const std::string &name = "COSMICS", const std::string &detector = "FULL", const int &debug = 1)
  : PHG4ParticleGeneratorBase(name)
{
  _x_max = 264.71;
  _x_min = 183.3;
  _R_earth = 637100000;
  _d_earth = 1000000;
  _z_max = 304.91;


  if (detector == "HCALSECTOR"){
    _y_top_fix = 200;
    _y_bot_fix = 0;
  }
  else if(detector == "EMCAL"){
    _x_max =116.;
    _x_min = 92.;
    _y_top_fix = _x_max;
    _y_bot_fix= -1*_x_max;
  } 
  else {
    _y_top_fix = _x_max;
    _y_bot_fix = -1*_x_max;
  }
  _f_vertex = new TF1("_f_vertex","(1/TMath::Power(TMath::Pi(),2))*((TMath::Pi()/2)-fabs(x)+TMath::Pi()/2-sin(fabs(x)+TMath::Pi()/2)*cos(fabs(x)+TMath::Pi()/2))",-1*TMath::Pi(), TMath::Pi());

  _f_momentum = new TF2("_f_momentum",CosmicSpray::EnergyAngularDistribution, 0, TMath::Pi()/2,2, 5);

  _debug = debug;
  return;
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
  TRandom *randomGen = new TRandom();
  randomGen->SetSeed();

  gun_z = randomGen->Uniform(-1*_z_max, _z_max);
  double theta_ran;
  theta_ran = _f_vertex->GetRandom();
  gun_x = _x_max * sin(theta_ran);
  gun_y = _x_max * cos(theta_ran);
  if(_debug) std::cout<<"Vertex: "<<gun_x<<" / " <<gun_y<<" / "<<gun_z<<std::endl;

  double er, tr,t_lim, pr, p_max;
  t_lim = 0;
  if (theta_ran > TMath::Pi()/2){
    t_lim = fabs(theta_ran) - TMath::Pi()/2;
    _f_momentum->SetRange(t_lim, TMath::Pi()/2, 2,5);
  }
  else if (theta_ran < -1*TMath::Pi()/2){
    t_lim = fabs(theta_ran) - TMath::Pi()/2;
    _f_momentum->SetRange(t_lim, TMath::Pi()/2, 2,5);
  }
  else {
    t_lim = TMath::Pi()/2 - fabs(theta_ran);
    _f_momentum->SetRange(0, TMath::Pi()/2, 2,5);
  }
  _f_momentum->GetRandom2(tr, er);
  _gun_e = TMath::Power(10, er - 3);

  if (theta_ran > TMath::Pi()/2){
    p_max = TMath::ACos(sin(t_lim)/sin(tr));
  }
  else if (theta_ran < -1*TMath::Pi()/2){
    p_max = TMath::ACos(sin(t_lim)/sin(tr));
  }
  else {
    if (tr < t_lim) p_max = TMath::Pi();
    else p_max = TMath::Pi()/2 + TMath::ASin(sin(t_lim)/sin(tr));
  }

  pr = (randomGen->Uniform() - 0.5 )*2*p_max;

  if (theta_ran > 0){
    if (pr > 0) pr = pr + TMath::Pi() - p_max;
    else if (pr < 0) pr = pr - TMath::Pi() + p_max;
  }

  gun_px = _gun_e*sin(tr)*sin(pr);
  gun_py = -1*fabs(_gun_e*cos(tr));
  gun_pz = _gun_e*sin(tr)*cos(pr);


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
