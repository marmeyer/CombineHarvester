#include "CombineHarvester/CombinePdfs/interface/MorphFunctions.h"
#include "CombineHarvester/CombineTools/interface/Algorithm.h"
#include "CombineHarvester/CombineTools/interface/AutoRebin.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"
#include "CombineHarvester/CombineTools/interface/CardWriter.h"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/HTTSM2017/interface/HttSystematics_SMRun2.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TF1.h"
#include "TH2.h"
#include "TString.h"
#include "boost/algorithm/string/predicate.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/regex.hpp"
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <math.h>

using namespace std;
using boost::starts_with;
namespace po = boost::program_options;

int main(int argc, char **argv) {
  typedef vector<string> VString;
  typedef vector<pair<int, string>> Categories;
  using ch::syst::bin_id;
  using ch::JoinStr;

  // Define program options
  string output_folder = "sm_run2";
  string base_path = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/HTTSM2017/shapes";
  string input_folder_em = "Vienna/";
  string input_folder_et = "Vienna/";
  string input_folder_mt = "Vienna/";
  string input_folder_tt = "Vienna/";
  string chan = "all";
  string postfix = "-ML";
  bool regional_jec = true;
  bool ggh_wg1 = true;
  bool auto_rebin = false;
  bool manual_rebin = false;
  bool real_data = false;
  bool jetfakes = true;
  bool embedding = false;
  bool verbose = false;
  string stxs_signals = "stxs_stage0"; // "stxs_stage0" or "stxs_stage1"
  string categories = "stxs_stage0"; // "stxs_stage0", "stxs_stage1" or "gof"
  string gof_category_name = "gof";
  int era = 2016; // 2016 or 2017
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()
      ("base_path", po::value<string>(&base_path)->default_value(base_path))
      ("input_folder_em", po::value<string>(&input_folder_em)->default_value(input_folder_em))
      ("input_folder_et", po::value<string>(&input_folder_et)->default_value(input_folder_et))
      ("input_folder_mt", po::value<string>(&input_folder_mt)->default_value(input_folder_mt))
      ("input_folder_tt", po::value<string>(&input_folder_tt)->default_value(input_folder_tt))
      ("postfix", po::value<string>(&postfix)->default_value(postfix))
      ("channel", po::value<string>(&chan)->default_value(chan))
      ("auto_rebin", po::value<bool>(&auto_rebin)->default_value(auto_rebin))
      ("regional_jec", po::value<bool>(&regional_jec)->default_value(regional_jec))
      ("ggh_wg1", po::value<bool>(&ggh_wg1)->default_value(ggh_wg1))
      ("real_data", po::value<bool>(&real_data)->default_value(real_data))
      ("manual_rebin", po::value<bool>(&manual_rebin)->default_value(manual_rebin))
      ("verbose", po::value<bool>(&verbose)->default_value(verbose))
      ("output_folder", po::value<string>(&output_folder)->default_value(output_folder))
      ("stxs_signals", po::value<string>(&stxs_signals)->default_value(stxs_signals))
      ("categories", po::value<string>(&categories)->default_value(categories))
      ("gof_category_name", po::value<string>(&gof_category_name)->default_value(gof_category_name))
      ("jetfakes", po::value<bool>(&jetfakes)->default_value(jetfakes))
      ("embedding", po::value<bool>(&embedding)->default_value(embedding))
      ("era", po::value<int>(&era)->default_value(era));
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  // Define the location of the "auxiliaries" directory where we can
  // source the input files containing the datacard shapes
  std::map<string, string> input_dir;
  input_dir["mt"] = base_path + "/" + input_folder_mt + "/";
  input_dir["et"] = base_path + "/" + input_folder_et + "/";
  input_dir["tt"] = base_path + "/" + input_folder_tt + "/";
  input_dir["em"] = base_path + "/" + input_folder_em + "/";

  // Define channels
  VString chns;
  if (chan.find("mt") != std::string::npos)
    chns.push_back("mt");
  if (chan.find("et") != std::string::npos)
    chns.push_back("et");
  if (chan.find("tt") != std::string::npos)
    chns.push_back("tt");
  if (chan.find("em") != std::string::npos)
    chns.push_back("em");
  if (chan == "all")
    chns = {"mt", "et", "tt", "em"};

  // Define background processes
  map<string, VString> bkg_procs;
  VString bkgs, bkgs_em;
  bkgs = {"W", "ZTT", "QCD", "ZL", "ZJ", "TTT", "TTL", "TTJ", "VVJ", "VVT", "VVL", "WH125", "ZH125"};
  bkgs_em = {"W", "ZTT", "QCD", "ZL", "TT", "VV", "ST", "WH125", "ZH125"};
  if(embedding){
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "ZTT"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "TTT"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "VVT"), bkgs.end());
    bkgs = JoinStr({bkgs,{"EMB"}});
    bkgs_em.erase(std::remove(bkgs_em.begin(), bkgs_em.end(), "ZTT"), bkgs_em.end());
    bkgs_em = JoinStr({bkgs_em,{"EMB"}});
  }
  if(jetfakes){
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "QCD"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "W"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "VVJ"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "TTJ"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "ZJ"), bkgs.end());
    bkgs = JoinStr({bkgs,{"jetFakes"}});
  }

  std::cout << "[INFO] Considerung the following processes:\n";
  if (chan.find("em") != std::string::npos) {
    std::cout << "For em channel : \n";
    for (unsigned int i=0; i < bkgs_em.size(); i++) std::cout << bkgs_em[i] << std::endl;
  }
  if (chan.find("mt") != std::string::npos || chan.find("et") != std::string::npos || chan.find("tt") != std::string::npos) {
    std::cout << "For et,mt,tt channels : \n";
    for (unsigned int i=0; i < bkgs.size(); i++) std::cout << bkgs[i] << std::endl;
  }
  bkg_procs["et"] = bkgs;
  bkg_procs["mt"] = bkgs;
  bkg_procs["tt"] = bkgs;
  bkg_procs["em"] = bkgs_em;

  // Define categories
  map<string, Categories> cats;
  // STXS stage 0 categories (optimized on ggH and VBF)
  if(categories == "stxs_stage0"){
    cats["et"] = {
        { 1, "et_ggh"},
        { 2, "et_qqh"},
        {11, "et_w"},
        {12, "et_ztt"},
        {13, "et_tt"},
        {14, "et_ss"},
        {15, "et_zll"},
        {16, "et_misc"},
    };
     cats["mt"] = {
        { 1, "mt_ggh"},
        { 2, "mt_qqh"},
        {11, "mt_w"},
        {12, "mt_ztt"},
        {13, "mt_tt"},
        {14, "mt_ss"},
        {15, "mt_zll"},
        {16, "mt_misc"},
    };
     cats["tt"] = {
        { 1, "tt_ggh"},
        { 2, "tt_qqh"},
        {12, "tt_ztt"},
        {16, "tt_misc"},
        {17, "tt_noniso"},
    };
     cats["em"] = {
        { 1, "em_ggh"},
        { 2, "em_qqh"},
        {12, "em_ztt"},
        {13, "em_tt"},
        {14, "em_ss"},
        {16, "em_misc"},
        {18, "em_st"},
        {19, "em_vv"},
    };
  }
  // STXS stage 1 categories (optimized on STXS stage 1 splits of ggH and VBF)
  else if(categories == "stxs_stage1"){
    cats["et"] = {
        { 1, "et_ggh_unrolled"},
        { 2, "et_qqh_unrolled"},
        {11, "et_w"},
        {12, "et_ztt"},
        {13, "et_tt"},
        {14, "et_ss"},
        {15, "et_zll"},
        {16, "et_misc"},
    };
     cats["mt"] = {
        { 1, "mt_ggh_unrolled"},
        { 2, "mt_qqh_unrolled"},
        {11, "mt_w"},
        {12, "mt_ztt"},
        {13, "mt_tt"},
        {14, "mt_ss"},
        {15, "mt_zll"},
        {16, "mt_misc"},
    };
     cats["tt"] = {
        { 1, "tt_ggh_unrolled"},
        { 2, "tt_qqh_unrolled"},
        {12, "tt_ztt"},
        {16, "tt_misc"},
        {17, "tt_noniso"},
    };
     cats["em"] = {
        { 1, "em_ggh_unrolled"},
        { 2, "em_qqh_unrolled"},
        {12, "em_ztt"},
        {13, "em_tt"},
        {14, "em_ss"},
        {16, "em_misc"},
        {18, "em_st"},
        {19, "em_vv"},
    };
  }
  else if(categories == "gof"){
    cats["et"] = {
        { 100, gof_category_name.c_str() },
    };
    cats["mt"] = {
        { 100, gof_category_name.c_str() },
    };
    cats["tt"] = {
        { 100, gof_category_name.c_str() },
    };
    cats["em"] = {
        { 100, gof_category_name.c_str() },
    };
  }
  else throw std::runtime_error("Given categorization is not known.");

  // Specify signal processes and masses
  vector<string> sig_procs;
  // STXS stage 0: ggH and VBF processes
  if(stxs_signals == "stxs_stage0") sig_procs = {"ggH", "qqH"};
  // STXS stage 1: Splits of ggH and VBF processes
  // References:
  // - https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGFiducialAndSTXS
  // - https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG2
  else if(stxs_signals == "stxs_stage1") sig_procs = {
      // ggH
      "ggH_0J", "ggH_1J_PTH_0_60", "ggH_1J_PTH_60_120", "ggH_1J_PTH_120_200",
      "ggH_1J_PTH_GT200", "ggH_GE2J_PTH_0_60", "ggH_GE2J_PTH_60_120",
      "ggH_GE2J_PTH_120_200", "ggH_GE2J_PTH_GT200", "ggH_VBFTOPO_JET3VETO",
      "ggH_VBFTOPO_JET3",
      // VBF
      "qqH_VBFTOPO_JET3VETO", "qqH_VBFTOPO_JET3", "qqH_REST",
      "qqH_PTJET1_GT200", "qqH_VH2JET"};
  else throw std::runtime_error("Given STXS signals are not known.");
  vector<string> masses = {"125"};

  // Create combine harverster object
  ch::CombineHarvester cb;

  // Add observations and processes
  std::string era_tag;
  if (era == 2016) era_tag = "Run2016";
  else if (era == 2017) era_tag = "Run2017";
  else std::runtime_error("Given era is not implemented.");

  for (auto chn : chns) {
    cb.AddObservations({"*"}, {"htt"}, {era_tag}, {chn}, cats[chn]);
    cb.AddProcesses({"*"}, {"htt"}, {era_tag}, {chn}, bkg_procs[chn], cats[chn],
                    false);
    cb.AddProcesses(masses, {"htt"}, {era_tag}, {chn}, sig_procs, cats[chn],
                    true);
  }

  // Add systematics
  ch::AddSMRun2Systematics(cb, jetfakes, embedding, regional_jec, ggh_wg1, era);

  // Extract shapes from input ROOT files
  for (string chn : chns) {
    cb.cp().channel({chn}).backgrounds().ExtractShapes(
        input_dir[chn] + "htt_" + chn + ".inputs-sm-" + era_tag + postfix + ".root",
        "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC");
    cb.cp().channel({chn}).process(sig_procs).ExtractShapes(
        input_dir[chn] + "htt_" + chn + ".inputs-sm-" + era_tag + postfix + ".root",
        "$BIN/$PROCESS$MASS", "$BIN/$PROCESS$MASS_$SYSTEMATIC");
  }

  // Delete processes with 0 yield
  cb.FilterProcs([&](ch::Process *p) {
    bool null_yield = !(p->rate() > 0.0);
    if (null_yield) {
      std::cout << "[INFO] Removing process with null yield: \n ";
      std::cout << ch::Process::PrintHeader << *p << "\n";
      cb.FilterSysts([&](ch::Systematic *s) {
        bool remove_syst = (MatchingProcess(*p, *s));
        return remove_syst;
      });
    }
    return null_yield;
  });

  // Replacing observation with the sum of the backgrounds (Asimov data)
  // useful to be able to check this, so don't do the replacement
  // for these
  if (!real_data) {
    for (auto b : cb.cp().bin_set()) {
      std::cout << "[INFO] Replacing data with asimov in bin " << b << "\n";
      cb.cp().bin({b}).ForEachObs([&](ch::Observation *obs) {
        obs->set_shape(cb.cp().bin({b}).backgrounds().GetShape() +
                           cb.cp().bin({b}).signals().GetShape(),
                       true);
      });
    }
  }

  // Rebin manually
  if (manual_rebin) {
    vector<double> binning = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    for (auto b : cb.cp().bin_set()) {
      std::cout << "[INFO] Rebinning by hand for bin: " << b << std::endl;
      cb.cp().bin({b}).VariableRebin(binning);
    }
  }

  // At this point we can fix the negative bins
  std::cout << "[INFO] Fixing negative bins.\n";
  cb.ForEachProc([](ch::Process *p) {
    if (ch::HasNegativeBins(p->shape())) {
      auto newhist = p->ClonedShape();
      ch::ZeroNegativeBins(newhist.get());
      p->set_shape(std::move(newhist), false);
    }
  });

  cb.ForEachSyst([](ch::Systematic *s) {
    if (s->type().find("shape") == std::string::npos)
      return;
    if (ch::HasNegativeBins(s->shape_u()) ||
        ch::HasNegativeBins(s->shape_d())) {
      auto newhist_u = s->ClonedShapeU();
      auto newhist_d = s->ClonedShapeD();
      ch::ZeroNegativeBins(newhist_u.get());
      ch::ZeroNegativeBins(newhist_d.get());
      s->set_shapes(std::move(newhist_u), std::move(newhist_d), nullptr);
    }
  });

  // Perform auto-rebinning
  if (auto_rebin) {
    for (auto b : cb.cp().bin_set()) {
      TString bstr = b;
      auto shape = cb.cp().bin({b}).backgrounds().GetShape();
      auto min = shape.GetBinLowEdge(1);
      auto max = shape.GetBinLowEdge(shape.GetNbinsX()+1);
      std::cout << "[INFO] Rebin " << b << " with boundaries (" << min << ", " << max << ")\n";
      if ( max == 1.0 ) {
          cb.cp().bin({b}).VariableRebin({min, max});
      }
      else if ( min == 0.125 && max == 8.0 ) {
          cb.cp().bin({b}).VariableRebin({
      0.125,  0.15 ,  0.2  ,  0.25 ,  0.3  ,  0.35 ,  0.4  ,  0.45 ,  0.5  ,  0.55 ,  0.6, 0.65 ,  0.7  ,  0.75 ,  0.8  ,  0.85 ,  0.9  ,  0.95 ,  1.   ,
      1.025,  1.075, 1.125,  1.175,  1.225,  1.275,  1.325,  1.375,  1.425,  1.475,  1.525,  1.575, 1.625,  1.675,  1.725,  1.775,  1.825,  1.875,
      1.9  ,  1.95 ,  2.   ,  2.05 , 2.1, 2.15 ,  2.2  ,  2.25 ,  2.3  ,  2.35 ,  2.4  ,  2.45 ,  2.5  ,  2.55 ,  2.6  , 2.65, 2.7  ,  2.75 ,
      2.775,  2.825,  2.875,  2.925,  2.975,  3.025,  3.075,  3.125, 3.175,  3.225,  3.275,  3.325,  3.375,  3.425,  3.475,  3.525,  3.575,  3.625,
      3.65 ,  3.7  ,  3.75 ,  3.8  ,  3.85 ,  3.9  ,  3.95 ,  4.   ,  4.05 ,  4.1  , 4.15, 4.2  ,  4.25 ,  4.3  ,  4.35 ,  4.4  ,  4.45 ,  4.5 ,
      4.525,  4.575,  4.625, 4.675,  4.725,  4.775,  4.825,  4.875,  4.925,  4.975,  5.025,  5.075,  5.125, 5.175,  5.225,  5.275,  5.325,  5.375,
      5.4  ,  5.45 ,  5.5  ,  5.55 ,  5.6  , 5.65, 5.7  ,  5.75 ,  5.8  ,  5.85 ,  5.9  ,  5.95 ,  6.   ,  6.05 ,  6.1  ,  6.15 , 6.2, 6.25 ,
      6.275,  6.325,  6.375,  6.425,  6.475,  6.525,  6.575,  6.625,  6.675, 6.725,  6.775,  6.825,  6.875,  6.925,  6.975,  7.025,  7.075,  7.125,
      7.15 , 7.2, 7.25 ,  7.3  ,  7.35 ,  7.4  ,  7.45 ,  7.5  ,  7.55 ,  7.6  ,  7.65 ,  7.7  , 7.75, 7.8  ,  7.85 ,  7.9  ,  7.95 ,  8.
                  });
      }
      else if ( min == 0.125 && max == 4.5 ) {
          cb.cp().bin({b}).VariableRebin({
      0.125,  0.15 ,  0.2  ,  0.25 ,  0.3  ,  0.35 ,  0.4  ,  0.45 ,  0.5  ,  0.55 ,  0.6, 0.65 ,  0.7  ,  0.75 ,  0.8  ,  0.85 ,  0.9  ,  0.95 ,  1.   ,
      1.025,  1.075, 1.125,  1.175,  1.225,  1.275,  1.325,  1.375,  1.425,  1.475,  1.525,  1.575, 1.625,  1.675,  1.725,  1.775,  1.825,  1.875,
      1.9  ,  1.95 ,  2.   ,  2.05 , 2.1, 2.15 ,  2.2  ,  2.25 ,  2.3  ,  2.35 ,  2.4  ,  2.45 ,  2.5  ,  2.55 ,  2.6  , 2.65, 2.7  ,  2.75 ,
      2.775,  2.825,  2.875,  2.925,  2.975,  3.025,  3.075,  3.125, 3.175,  3.225,  3.275,  3.325,  3.375,  3.425,  3.475,  3.525,  3.575,  3.625,
      3.65 ,  3.7  ,  3.75 ,  3.8  ,  3.85 ,  3.9  ,  3.95 ,  4.   ,  4.05 ,  4.1  , 4.15, 4.2  ,  4.25 ,  4.3  ,  4.35 ,  4.4  ,  4.45 ,  4.5
                  });
      }
      else if ( min == 0.2 && max == 7.4 ) {
          cb.cp().bin({b}).VariableRebin({
     0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,  0.45,  0.5 ,  0.55,  0.6 ,  0.65,  0.7 ,  0.75, 0.8  , 0.85,  0.9 ,  0.95,  1.  ,
     1.05,  1.1 ,  1.15,  1.2 ,  1.25,  1.3 ,  1.35, 1.4 ,  1.45,  1.5 ,  1.55,  1.6 , 1.65,  1.7 ,  1.75,  1.8 ,
     1.85,  1.9 ,  1.95, 2.  ,  2.05,  2.1 ,  2.15,  2.2 ,  2.25,  2.3 ,  2.35,  2.4 , 2.45,  2.5 ,  2.55, 2.6 ,
     2.65,  2.7 ,  2.75,  2.8 ,  2.85,  2.9 ,  2.95,  3.  ,  3.05,  3.1 ,  3.15, 3.2 ,  3.25,  3.3 ,  3.35,  3.4 ,
     3.45,  3.5 ,  3.55,  3.6 ,  3.65,  3.7 ,  3.75, 3.8 ,  3.85,  3.9 ,  3.95,  4.  ,  4.05,  4.1 ,  4.15,  4.2 ,
     4.25,  4.3 ,  4.35, 4.4 ,  4.45,  4.5 ,  4.55,  4.6 ,  4.65,  4.7 ,  4.75,  4.8 ,  4.85,  4.9 ,  4.95, 5.  ,
     5.05,  5.1 ,  5.15,  5.2 ,  5.25,  5.3 ,  5.35,  5.4 ,  5.45,  5.5 ,  5.55, 5.6 ,  5.65,  5.7 ,  5.75,  5.8 ,
     5.85,  5.9 ,  5.95,  6.  ,  6.05,  6.1 ,  6.15, 6.2 ,  6.25,  6.3 ,  6.35,  6.4 ,  6.45,  6.5 ,  6.55,  6.6 ,
     6.65,  6.7 ,  6.75, 6.8 ,  6.85,  6.9 ,  6.95,  7.  ,  7.05,  7.1 ,  7.15,  7.2 ,  7.25,  7.3 ,  7.35, 7.4
                  });
      }
      else if ( min == 0.2 && max == 4.2 ) {
          cb.cp().bin({b}).VariableRebin({
      0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,  0.45,  0.5 ,  0.55,  0.6 ,  0.65,  0.7 ,  0.75, 0.8 ,  0.85,  0.9 ,  0.95,  1.  ,
      1.05,  1.1 ,  1.15,  1.2 ,  1.25,  1.3 ,  1.35, 1.4 ,  1.45,  1.5 ,  1.55,  1.6 ,  1.65,  1.7 ,  1.75,  1.8 ,
      1.85,  1.9 ,  1.95, 2.  ,  2.05,  2.1 ,  2.15,  2.2 ,  2.25,  2.3 ,  2.35,  2.4 ,  2.45,  2.5 ,  2.55, 2.6 ,
      2.65,  2.7 ,  2.75,  2.8 ,  2.85,  2.9 ,  2.95,  3.  ,  3.05,  3.1 ,  3.15, 3.2 ,  3.25,  3.3 ,  3.35,  3.4 ,
      3.45,  3.5 ,  3.55,  3.6 ,  3.65,  3.7 ,  3.75, 3.8 ,  3.85,  3.9 ,  3.95,  4.  ,  4.05,  4.1 ,  4.15,  4.2
                  });
      }
      else {
          throw std::runtime_error("[ERROR] Given configuration is not valid for rebinning.");
      }
    }
  }

  // Merge bins and set bin-by-bin uncertainties
  /*
  auto bbb = ch::BinByBinFactory()
                 .SetAddThreshold(0.00)
                 .SetMergeThreshold(0.5)
                 .SetFixNorm(false);
  bbb.MergeBinErrors(cb.cp().backgrounds());
  bbb.AddBinByBin(cb.cp().backgrounds(), cb);
  */

  // This function modifies every entry to have a standardised bin name of
  // the form: {analysis}_{channel}_{bin_id}_{era}
  ch::SetStandardBinNames(cb, "$ANALYSIS_$CHANNEL_$BINID_$ERA");

  // Write out datacards. Naming convention important for rest of workflow. We
  // make one directory per chn-cat, one per chn and cmb. In this code we only
  // store the individual datacards for each directory to be combined later.
  string output_prefix = "output/";
  ch::CardWriter writer(output_prefix + output_folder + "/$TAG/$MASS/$BIN.txt",
                        output_prefix + output_folder +
                            "/$TAG/common/htt_input_" + era_tag + ".root");

  // We're not using mass as an identifier - which we need to tell the
  // CardWriter
  // otherwise it will see "*" as the mass value for every object and skip it
  //    writer.SetWildcardMasses({});
  if (verbose)
    writer.SetVerbosity(1);

  // Write datacards combined and per channel
  writer.WriteCards("cmb", cb);

  for (auto chn : chns) {
    writer.WriteCards(chn, cb.cp().channel({chn}));
  }

  if (verbose)
    cb.PrintAll();

  cout << "[INFO] Done producing datacards.\n";
}
