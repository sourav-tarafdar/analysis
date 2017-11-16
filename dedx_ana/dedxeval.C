# include "dedxeval.h"

#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <g4hough/SvtxVertexMap.h>
#include <g4hough/SvtxVertex.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrack.h>
#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxCluster.h>
#include <g4hough/SvtxHitMap.h>
#include <g4hough/SvtxHit.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4TruthInfoContainer.h>


#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>
#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxTruthEval.h>
#include <g4eval/SvtxVertexEval.h>
#include <g4eval/SvtxHitEval.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TMath.h>
#include <math.h>

#include <iostream>
#include <set>
#include <cmath>
#include <cassert>
#include <vector>

using namespace std;
using std::vector;
using std::pair;

dEDxEval::dEDxEval(const string &name , const string &filename, const string &trackmapname,
		unsigned int nlayers_maps,
		unsigned int nlayers_intt,
		unsigned int nlayers_tpc) :
  SubsysReco("dEDxEval"), _ievent(0), _tr_g4dedx(NULL), _filename(filename), _tfile (NULL)
{
  verbosity = 0;
}

int dEDxEval::Init(PHCompositeNode *topNode) {
  
  _ievent = 0;

  _tfile = new TFile(_filename.c_str(), "RECREATE");
  
  _tr_g4dedx = new TTree("tr_g4dedx","dedx from G4");
  _tr_g4dedx->Branch("event", &event, "event/I");


  _tr_g4dedx->Branch("g4_layer", &g4_layer, "g4_layer/I");
  _tr_g4dedx->Branch("g4_trkid", &g4_trkid, "g4_trkid/I");
  _tr_g4dedx->Branch("g4hitid", &g4hitid, "g4hitid/F");
  _tr_g4dedx->Branch("g4x1", &g4x1, "g4x1/F");
  _tr_g4dedx->Branch("g4x0", &g4x0, "g4x0/F");
  _tr_g4dedx->Branch("g4y0", &g4y0, "g4y0/F");
  _tr_g4dedx->Branch("g4y1", &g4y1, "g4y1/F");
  _tr_g4dedx->Branch("g4z0", &g4z0, "g4z0/F");
  _tr_g4dedx->Branch("g4z1", &g4z1, "g4z1/F");
  _tr_g4dedx->Branch("g4dx", &g4dx, "g4dx/F");
  _tr_g4dedx->Branch("g4dx_tot", &g4dx_tot, "g4dx_tot/F");
  _tr_g4dedx->Branch("g4x", &g4x, "g4x/F");
  _tr_g4dedx->Branch("g4y", &g4y, "g4y/F");
  _tr_g4dedx->Branch("g4z", &g4z, "g4z/F");
  _tr_g4dedx->Branch("g4e", &g4e, "g4e/F");
  _tr_g4dedx->Branch("g4e_tot", &g4e_tot, "g4e_tot/F");
  _tr_g4dedx->Branch("g4px", &g4px, "g4px/F");
  _tr_g4dedx->Branch("g4py", &g4py, "g4py/F");
  _tr_g4dedx->Branch("g4pz", &g4pz, "g4pz/F");
  
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int dEDxEval::InitRun(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int dEDxEval::process_event(PHCompositeNode *topNode) {
  
  cout << "dEDxEvaluator::process_event - Event = " << _ievent << endl;
  
  if (!_svtxevalstack) {
    cout << " evalstk 1 :" << endl;
    _svtxevalstack = new SvtxEvalStack(topNode);
  }else {
    cout << "evalstk 2 :" << endl;
    _svtxevalstack->next_event(topNode);
  }
  

  //Fill out tree ---------

  fill_tree(topNode);

  ++_ievent;
  return Fun4AllReturnCodes::EVENT_OK;

}



int dEDxEval::End(PHCompositeNode *topNode) {
  
  _tfile->cd();

  _tr_g4dedx->Write();

  // _tr_dedx->Write();

  _tfile->Close();

  delete _tfile;

  delete _svtxevalstack;

  return Fun4AllReturnCodes::EVENT_OK;

}

// Start Filling the Tree with cluster variables and  truth info

void dEDxEval::fill_tree(PHCompositeNode *topNode) {
  
  reset_variables();

  event = _ievent;
  //SvtxEvalStack  _svtxevalstack(topNode);

  //SvtxClusterEval* clustereval = _svtxevalstack->get_cluster_eval();
  SvtxTruthEval*     trutheval = _svtxevalstack->get_truth_eval();
  //SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
  //SvtxHitMap* hitmap = findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");
  //SvtxHitEval*         hiteval = _svtxevalstack->get_hit_eval();
  
  g4dx_tot = 0;
  g4e_tot = 0;
   
  std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits();
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
	 iter != g4hits.end();
	 ++iter) {
            
      PHG4Hit *g4hit = *iter;
      PHG4Particle *g4particle = trutheval->get_particle(g4hit);
      
       g4hitid =  g4hit->get_hit_id();
       g4x0 =  g4hit->get_x(0);
       g4x1 =  g4hit->get_x(1);
       g4y0 =  g4hit->get_y(0);
       g4y1 =  g4hit->get_y(1);
       g4z0 =  g4hit->get_z(0);
       g4z1 =  g4hit->get_z(1);
       g4x = g4hit->get_avg_x();
       g4y = g4hit->get_avg_y();
       g4z = g4hit->get_avg_z();
       g4_trkid = g4hit->get_trkid();
       g4_layer =  g4hit->get_layer();
       float x1_tr = g4hit->get_x(0) - g4hit->get_x(1) ;
       float y1_tr = g4hit->get_y(0) - g4hit->get_y(1) ;
       float z1_tr  =  g4hit->get_z(0) - g4hit->get_z(1) ;
       g4dx = sqrt( x1_tr*x1_tr + y1_tr*y1_tr + z1_tr*z1_tr  );
       g4e = g4hit->get_edep();
       g4trkid1 = g4hit->get_trkid();
       //cout << " hit id :" << g4hitid << endl;
       if(g4particle) {
	 g4px      = g4particle->get_px();
	 g4py      = g4particle->get_py();
	 g4pz      = g4particle->get_pz();
       }
       
       _tr_g4dedx->Fill();
  } // g4 hits iterator

  
  
  return;

}

void dEDxEval::reset_variables(){
  
  event = -9999;
  
  g4_trkid = -9999;
  g4_layer = -9999.;
  g4hitid = -9999.;
  g4x1 = -9999.;
  g4x0 = -9999.;
  g4y1 = -9999.;
  g4y0 = -9999.;
  g4z1 = -9999.;
  g4z0 = -9999.;
  g4dx = -9999.;
  g4e = -9999.;
  g4px = -9999.;
  g4py = -9999.;
  g4pz = -9999.;
  g4e_tot = -9999.;
  g4dx_tot = -9999.;
  g4x = -9999.;
  g4y = -9999.;
  g4z = -9999.;
  /* 
  layer = -9999.;
  hitid = -9999.;
  dx = -9999.;
  e = -9999.;
  x = -9999.;
  y = -9999.;
  z = -9999.;
  */
    
}
