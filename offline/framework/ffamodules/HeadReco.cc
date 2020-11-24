#include "HeadReco.h"

#include <ffaobjects/EventHeader.h>
#include <ffaobjects/EventHeaderv1.h>
#include <ffaobjects/RunHeader.h>
#include <ffaobjects/RunHeaderv1.h>
#include <ffaobjects/FlagSave.h>
#include <ffaobjects/FlagSavev1.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>          // for PHIODataNode
#include <phool/PHNode.h>                // for PHNode
#include <phool/PHNodeIterator.h>        // for PHNodeIterator
#include <phool/PHObject.h>              // for PHObject
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TSystem.h>

#include <HepMC/GenEvent.h>
#include <HepMC/HeavyIon.h>              // for HeavyIon

#include <cstdlib>
#include <iostream>

using namespace std;

HeadReco::HeadReco(const std::string &name):
  SubsysReco(name)
{}	   

// the nodes need to be created here since at least one input manager uses
// the event header. Creating them in InitRun() will be too late
int HeadReco::Init(PHCompositeNode *topNode)
{
  
PHNodeIterator iter(topNode);
PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *> (iter.findFirst("PHCompositeNode","RUN"));
RunHeader *runheader = new RunHeaderv1();
PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(runheader,"RunHeader","PHObject");
runNode->addNode(newNode);


FlagSave *flgsv = new FlagSavev1();
newNode = new PHIODataNode<PHObject>(flgsv, "Flags", "PHObject");
runNode->addNode(newNode);

PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *> (iter.findFirst("PHCompositeNode","DST"));
EventHeader *eventheader = new EventHeaderv1();
newNode = new PHIODataNode<PHObject>(eventheader,"EventHeader","PHObject");
dstNode->addNode(newNode);

return Fun4AllReturnCodes::EVENT_OK;
}

int HeadReco::InitRun(PHCompositeNode *topNode)
{
recoConsts *rc = recoConsts::instance();
RunHeader *runheader = findNode::getClass<RunHeader>(topNode,"RunHeader");
runheader->set_RunNumber(rc->get_IntFlag("RUNNUMBER"));
return Fun4AllReturnCodes::EVENT_OK;
}

int
HeadReco::process_event(PHCompositeNode *topNode)
{

  Fun4AllServer *se = Fun4AllServer::instance();
  EventHeader *evtheader = findNode::getClass<EventHeader>(topNode,"EventHeader");
  PHHepMCGenEventMap *genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

  if (genevtmap)
  {
  for (PHHepMCGenEventMap::ReverseIter iter = genevtmap->rbegin(); iter != genevtmap->rend(); ++iter)
  {
  PHHepMCGenEvent *genevt = iter->second;
  int embed_flag = genevt->get_embedding_id();
  if (embed_flag == 0) // should be foreground event
  {

  HepMC::GenEvent *hepmcevt = genevt->getEvent();

  if (hepmcevt)
  {
  HepMC::HeavyIon *hi = hepmcevt->heavy_ion();
  if (hi)
  {
  evtheader->set_ImpactParameter(hi->impact_parameter());
  evtheader->set_EventPlaneAngle(hi->event_plane_angle());
  evtheader->set_eccentricity(hi->eccentricity());
  evtheader->set_ncoll(hi->Ncoll());
  evtheader->set_npart(hi->Npart_targ()+hi->Npart_proj());
  }
  }
  }
  }
  }
  evtheader->set_RunNumber(se->RunNumber());
  evtheader->set_EvtSequence(se->EventNumber());
  if (Verbosity() > 0)
  {
  evtheader->identify();
  }

return Fun4AllReturnCodes::EVENT_OK;
}

int HeadReco::EndRun(const int runno)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  FlagSave* flagsave = findNode::getClass<FlagSave>(se->topNode(), "Flags");
  if (flagsave)
    {
      recoConsts *rc = recoConsts::instance();
      cout << "Saving recoConst Flags: " << endl;
      flagsave->FillFromPHFlag(rc);
      flagsave->identify();
    }
  return Fun4AllReturnCodes::EVENT_OK;
}
