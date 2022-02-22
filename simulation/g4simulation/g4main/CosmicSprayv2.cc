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

double CosmicSpray::_y_fix;
double CosmicSpray::_y_top_fix;
double CosmicSpray::_y_bot_fix;
double CosmicSpray::_x_max;
double CosmicSpray::_x_min;
double CosmicSpray::_z_max;
double CosmicSpray::_z_min;
double CosmicSpray::_gun_e;
double CosmicSpray::_offset;

class PHCompositeNode;
class PHG4Particle;
class PHG4ParticleGeneratorBase;

CosmicSpray::CosmicSpray(const std::string &name = "COSMICS", const double &y_top_fix = 253., const double &y_bot_fix = 43., const double &z_max = 304.91, const double &x_max = 264.71, const double &x_min = 183.3, const int &debug = 1)
  : PHG4ParticleGeneratorBase(name)
{
  _gun_e = 4;
  _offset = 100;
  _y_top_fix = y_top_fix;
  _y_bot_fix = y_bot_fix;
  _y_fix = y_top_fix - y_bot_fix;
  _z_max = z_max + _offset;
  _z_min = -1*z_max - _offset;
  _x_max = x_max + _offset;
  _x_min = x_min - _offset;
  _debug = debug;
  _ftop = new TF2("ftop",CosmicSpray::TopDistributionFunction,x_min,x_max,-1*z_max,z_max);
  _fbottom = new TF2("fbottom",CosmicSpray::BottomDistributionFunction,-1*_gun_e,_gun_e,-1*_gun_e,_gun_e);
  return;
}

void CosmicSpray::set_gun_energy(const double &e){
  _gun_e = e;
  _fbottom = new TF2("fbottom",CosmicSpray::BottomDistributionFunction,-1*_gun_e,_gun_e,-1*_gun_e,_gun_e);

}

void CosmicSpray::set_offset_distance(const double &o){
  _offset = o;
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
  double gun_x = 0, gun_y = 0, gun_z = 0;
  double gun_px = 0, gun_py = 0, gun_pz = 0;
  double botx = _gun_e;
  double botz = _gun_e;
  double dy = 0;
  
  TRandom *randomGen = new TRandom();
  while (sqrt(botx*botx+botz*botz) > _gun_e){
    randomGen->SetSeed();
    gun_x = randomGen->Uniform(_x_min, _x_max);
    gun_z = randomGen->Uniform(_z_min, _z_max);
    
    if(_debug) std::cout<<"Vertex: "<<gun_x<<" / " <<gun_z<<std::endl;
    _fbottom->SetParameters(gun_x, gun_z);

    double dx, dz;
    _fbottom->GetRandom2(dx, dz);
    dy = -1 *sqrt(_gun_e*_gun_e - dx*dx - dz*dz);
    TVector3 xy(dx, dy, 0);
    TVector3 zy(0, dy, dz);
    TVector3 y(0,-1,0);
    double ax = y.Angle(xy);
    double az = y.Angle(zy);
    
    double ddx = _y_fix * tan(ax);
    double ddz = _y_fix * tan(az);
    if (((gun_x + ddx) > (_x_max-_offset)) || ((gun_x + ddx) < (_x_min + _offset)) || ((gun_z + ddz) > (_z_max - _offset)) || ((gun_z - ddz) < (-1*_z_max + _offset))) continue;
    
    botx = dx;
    botz = dz;    
  }
  TVector3 dir(botx, dy, botz);
  
  //if(_debug) std::cout<<"Dir: "<<dx<<" / " <<dy<< " / " <<dz<<std::endl;
  if(_debug) std::cout<<"Norm Dir: "<<dir.X()<<" / " <<dir.Y()<<" / "<<dir.Z()<<std::endl;
  
  gun_y = _y_top_fix;

  gun_px = dir.X();
  gun_py = dir.Y();
  gun_pz = dir.Z();
  double tot_p = sqrt(gun_px*gun_px+gun_py*gun_py+gun_pz*gun_pz);
  
  if(_debug) std::cout<<"Momentum: "<<gun_px<<" / "<<gun_py<<" / "<<gun_pz<<std::endl;
  if(_debug)std::cout<<"total mom: "<<tot_p<<std::endl;
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

double CosmicSpray::TopDistributionFunction(double *val, double *par) {
  double a1 = atan2(abs(_x_min - val[0]), _y_fix);
  double a2 = atan2(abs(_z_max - val[1]), _y_fix);
  double b1 = atan2(abs(_x_max - val[0]), _y_fix);
  double b2 = atan2(abs(_z_max + val[1]), _y_fix);
  double r1 = a1+b1+(1/2)*(sin(2*a1)+sin(2*b1));
  double r2 = a2+b2+(1/2)*(sin(2*a2)+sin(2*b2));
  return r1*r2/4;
}

double CosmicSpray::BottomDistributionFunction(double *val, double *par) {
  //double a1 = atan2(_x_min - (val[0]-par[0]), _y_fix);
  //double a2 = atan2(-1*_z_max - (val[1]-par[1]), _y_fix);
  //double b1 = atan2(_x_max - (val[0]-par[0]), _y_fix);
  //double b2 = atan2(_z_max - (val[1]-par[1]), _y_fix);
  
  double a1 = TMath::ASin((sqrt(val[0]*val[0] + val[1]*val[1]))/ _gun_e);
  double p = cos(a1)*cos(a1);
  return p;
}
