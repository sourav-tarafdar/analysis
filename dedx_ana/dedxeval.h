#ifndef __DEDXEVAL_H__
#define __DEDXEVAL_H__

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

using std::vector;
using std::pair;
class PHCompositeNode;

class SvtxEvalStack;
class TFile;
class TNtuple;
class TTree;

class dEDxEval : public SubsysReco 
{
 
public :

  dEDxEval(const std::string &name = "DEDXEVAL",
	   const std::string &filename = "dEDxEval.root",
	   const std::string &trackmapname = "SvtxTrackMap",
	   unsigned int nlayers_maps = 3,
	   unsigned int nlayers_intt = 4,
	   unsigned int nlayers_tpc = 40);
  virtual ~dEDxEval(){}
  
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void fill_tree(PHCompositeNode *topNode);
  

  void reset_variables();

  private:

  unsigned int _ievent;

  // eval stack
  SvtxEvalStack* _svtxevalstack;
  
  //TTree *_tr_dedx;
  TTree *_tr_g4dedx;
  
  std::string _filename;
  TFile *_tfile;

  // Variables for comparing succesive elements of Hitmap container

  int g4trkid1, g4trkid2;

  // Variables that goes into tree
  int event,  g4_trkid, g4_layer;
  float g4x1, g4x0, g4y1, g4y0, g4z0, g4z1, g4x, g4y, g4z, g4hitid, g4dx, g4e, g4px, g4py, g4pz, g4e_tot, g4dx_tot;

  int  trkid, layer;
  float x1, x0, y1, y0, z0, z1, x, y, z, hitid, dx, e, px, py, pz, e_tot, dx_tot;
  float x2, y2, z2;
  /*
  vector <pair< int, float>> x1_layer;
  vector <pair< int, float>> y1_layer;
  vector <pair< int, float>> z1_layer;
  vector <pair< int, float>> e_layer;
  */
  vector <vector<float>> x_layer;
  vector <vector<float>> y_layer;
  vector <vector<float>> z_layer; 
  
  //vector<vector<float>>dx_layer; 
};

#endif //__DEDXEVAL_H__

