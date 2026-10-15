/*
* ProduceRoIMap.cpp
*
*  Created on: 12 Oct 2019
*      Author: Valerio Formato
*/

// Example headers
#include "ProduceRoIMap.h"
#include "TriggerPattern.h"
#include "UlegGeomPsdRoi.h"

// HerdSoftware headers
#include "herd/analysis/dataobjects/Line.h"
#include "herd/analysis/dataobjects/MCTruth.h"
#include "herd/analysis/dataobjects/FitHitsColl.h"
#include "herd/analysis/dataobjects/FitParticleHitsColl.h"
#include "herd/analysis/dataobjects/PsdHitsColl.h"

#include "herd/analysis/dataobjects/CaloHits.h"
#include "herd/analysis/dataobjects/CaloGeoParams.h"


// EA headers
#include "utils/SmartLog.h"

// Root headers
// #include "TTree.h"
#include "TNtuple.h"



// C/C++ standard headers
#include <iterator>
#include <numeric>

RegisterAlgorithm(ProduceRoIMap);

ProduceRoIMap::ProduceRoIMap(const std::string &name) : Algorithm{name}, _threshold{0}, _maxDelta{1}, _patternSize{3}, _nHits{3}, _layer_start{0}, _force_full_spanning{false}, _disable_filter{false}, _ignoreIDs{{}}, _require_only_one_side{false}, _matHitsCollName{"fitMatHitsCollMC"}, _store_patterns_name{""}, _suppress_patterns_name{""}
, _dumpPatterns{false}, _requireCalo{false}, _storePatternToTileList{false}
{
  DefineParameter("threshold", _threshold);
  DefineParameter("maxDelta", _maxDelta);
  DefineParameter("nHits", _nHits);
  DefineParameter("patternSize", _patternSize);
  DefineParameter("force_full_spanning", _force_full_spanning);
  DefineParameter("disable_filter", _disable_filter);
  DefineParameter("startLayer", _layer_start);
  // DefineParameter("ignore_plusminus", _ignore_plusminus);
  DefineParameter("ignoreIDs", _ignoreIDs);
  DefineParameter("require_only_one_side", _require_only_one_side);
  DefineParameter("matHitsCollName", _matHitsCollName);
  DefineParameter("store_patterns_name", _store_patterns_name);
  DefineParameter("suppress_patterns_name", _suppress_patterns_name);
  DefineParameter("dumpPatterns", _dumpPatterns);
  DefineParameter("requireCalo", _requireCalo);
  DefineParameter("storePatternToTileList", _storePatternToTileList);


}

bool ProduceRoIMap::Initialize() {
  const std::string routineName("ProduceRoIMap::Initialize");


  // _treeTriggers = std::make_shared<TTree>(std::string("TriggerMats"+GetName()).c_str(), "Trigger signals");
  // _tupleTriggersEven = std::make_shared<TNtuple>(std::string("TriggersEven_"+GetName()).c_str(), "Trigger patterns", "eTrue:cosTheta:iLayer:iMat:patternID:oldID");
  // _tupleTriggers = std::make_shared<TNtuple>(std::string("Triggers_"+GetName()).c_str(), "Trigger patterns", "eTrue:cosTheta:iSide:iLayer:iMat:patternID:oldID:xPSD0:yPSD0:maxECalo:nHitsPsd:totMult");
  // _tupleFirstTriggersEven = std::make_shared<TNtuple>(std::string("FirstTriggersEven_"+GetName()).c_str(), "Trigger patterns", "eTrue:cosTheta:iLayer:iMat:patternID:oldID");
  // _tupleFirstTriggers = std::make_shared<TNtuple>(std::string("FirstTriggers_"+GetName()).c_str(), "Trigger patterns", "eTrue:cosTheta:iSide:iLayer:iMat:patternID:oldID:xPSD0:yPSD0:eDepPsd:maxECalo:nHitsPsd:width");
  // _tupleFirstLayer = std::make_shared<TNtuple>(std::string("FirstLayer_"+GetName()).c_str(), "Trigger patterns", "eTrue:cosTheta:iLayer:iMat:patternID:oldID:xPSD0:yPSD0:zPSD0:x:y:z:px:py:pz");
  // _tupleTriggersOdd = std::make_shared<TNtuple>(std::string("TriggersOdd_"+GetName()).c_str(), "Trigger patterns", "eTrue:cosTheta:iLayer:iMat:patternID:oldID");
  // _tupleFirstTriggersOdd = std::make_shared<TNtuple>(std::string("FirstTriggersOdd_"+GetName()).c_str(), "Trigger patterns", "eTrue:cosTheta:iLayer:iMat:patternID:oldID");
  // _tupleMats = std::make_shared<TNtuple>(std::string("Mats_"+GetName()).c_str(), "Triggered mats", "eTrue:cosTheta:iLayer:iMat:patternID:oldID:xPSD0:yPSD0:zPSD0:x:y:z:px:py:pz");
  // _tupleFirstMats = std::make_shared<TNtuple>(std::string("FirstMats_"+GetName()).c_str(), "Triggered mats", "eTrue:cosTheta:iLayer:iMat:patternID:oldID");

  // Setup the filter
  // SetFilterStatus(FilterStatus::ENABLED); // don't set up for now
  // _nAccept = _nReject = 0
  _patternsAndTiles=std::make_shared<TNtuple>(("patternsAndTiles"+GetName()).c_str(),"temp", "fitSide:fitLayer:mat0:mat1:mat2:iPsdSide:iPsdLayer:tileID:tileAngleFromMat");


  return true;
}

bool ProduceRoIMap::Process() {
  const std::string routineName("ProduceRoIMap::Process");
  return true;
}

bool ProduceRoIMap::Finalize() {
  const std::string routineName("ProduceRoIMap::Finalize");


  try {
    _globStore = GetDataStoreManager()->GetGlobalDataStore("globStore");
    if (_globStore->IsPresent("fitDigitizedGeoParamsColl")) {
      _fitGeoParamsColl = _globStore->GetObject<Herd::FitGeoParamsColl>("fitDigitizedGeoParamsColl");
    } else {
      _fitGeoParamsColl = _globStore->GetObject<Herd::FitGeoParamsColl>("fitGeoParamsColl");
    }
    _psdGeoParamsColl = _globStore->GetObject<Herd::PsdGeoParamsColl>("psdGeoParamsColl");
    _caloGeoParams = _globStore->GetObject<Herd::CaloGeoParams>("caloGeoParams");
    // _ScdGeoParamsColl = _globStore->GetObject<Herd::ScdGeoParamsColl>("scdGeoParamsColl");
    // _evStore = GetDataStoreManager()->GetEventDataStore("evStore");
  }
  catch (Retrieval::Exception &exc) {
    COUT(ERROR) << "Event data store not found." << ENDL;
    return false;
  }

    Herd::SidesArray<std::vector<TriggerPattern>> triggersColl;

    auto _patternPsdMap = std::make_shared<PatternVetoListMap<PsdVetoList,Herd::PsdGeoParamsColl>>(PatternVetoListMap<PsdVetoList,Herd::PsdGeoParamsColl>());
    auto _patternCaloMap = std::make_shared<PatternVetoListMap<CaloVetoList,Herd::CaloGeoParams>>(PatternVetoListMap<CaloVetoList,Herd::CaloGeoParams>());
    // int totalMultiplicity=0;
    // int triggeredSides=0;
    // // std::cout << "algo: " << GetName() << '\n';

    for (auto side : Herd::RefFrame::Sides) {
      auto &fitGeoParamsVector = (*_fitGeoParamsColl)[side];
      // auto &fitGeoParamsVector = (*_FitGeoParamsColl)[*side_selected];
      // auto &fitMatHitsSide = (*fitMatHitsColl)[side];
      // auto &fitMatHitsSide = (*fitMatHitsColl)[*side_selected];
      auto nLayers = fitGeoParamsVector.size();
      auto geoParamsEven = fitGeoParamsVector[0];
      auto geoParamsOdd = fitGeoParamsVector[1];
      auto primSegAxis = Herd::RefFrame::ToAxis(geoParamsEven.SegmentationDirections().first); //Z/Y
      // std::cout << "even seg: " << primSegAxis << '\n';
      // std::cout << geoParamsEven.MatPosition(0) << '\n';
      // std::cout << geoParamsEven.MatPosition(1) << '\n';
      auto secSegAxis = Herd::RefFrame::ToAxis(geoParamsOdd.SegmentationDirections().first); //H/X
      auto nMatsEven = geoParamsEven.NMats()[primSegAxis];
      // std::cout << nMatsEven << '\n';
      auto nMatsOdd = geoParamsOdd.NMats()[secSegAxis];

      std::vector<int> nMatsPair{nMatsEven,nMatsOdd};
      FitHitsArray hitsArrayEven(nLayers/2, std::vector<int>(nMatsEven,1)); // init the hit arrays to 1
      FitHitsArray hitsArrayOdd(nLayers/2, std::vector<int>(nMatsOdd,1)); // init the hit arrays to 1
      auto hitsArraysPair = std::vector<FitHitsArray>{hitsArrayEven,hitsArrayOdd};

      std::vector<TriggerPattern> vTriggers;
      // auto vTriggersPair = std::vector<std::vector<TriggerPattern>>{vTriggersEven,vTriggersOdd};

      for (size_t iLayer = _layer_start; iLayer < nLayers-2*(_patternSize-1); iLayer++) {
        int parity = (iLayer%2);
        // auto &triggersLayer=triggerArraysPair[parity];
        // auto &vTriggers=vTriggersPair[parity];
        auto &hitsArray=hitsArraysPair[parity];
        auto nMats=nMatsPair[parity];
        // std::vector<TriggerPattern> vTriggers;
        for (size_t iMat = 0; iMat < nMats; iMat++) {
          // std::cout << "mat: " << iMat << '\n';
          // std::cout << hitsArray[iLayer][iMat] << '\n';
          int hasHit = hitsArray[iLayer/2][iMat];
          // if (hasHit || _patternSize > _nHits) {
          TriggerPattern iPattern = TriggerPattern(side,iLayer,iMat,!hasHit);
          iPattern.FindPatterns(_nHits-hasHit,_patternSize,_maxDelta,hitsArray,vTriggers,_force_full_spanning);
          // }
        } // loop over mats
      } // loop over layers
      // std::cout << side << '\n';
      // std::cout << "found patterns: " << vTriggers.size() << '\n';
      // for (auto patt : vTriggers) std::cout << patt << '\n';
      std::sort(vTriggers.begin(),vTriggers.end()); // sort triggers, just in case
      if (!_ignoreIDs.empty()) { // pattern suppression by ID, loaded in the datacard
        // prune patters if their IDs match a member of _ignore_pattern_IDs
        auto idChecker = [ignoreIDs=_ignoreIDs](auto &pattern) {
          return std::any_of(ignoreIDs.begin(),ignoreIDs.end(),[&pattern](auto &id){return id==pattern.PatternID();});
        };
        auto it = std::remove_if(vTriggers.begin(),vTriggers.end(), idChecker);
        vTriggers.erase(it,vTriggers.end());
      }

      if (_requireCalo) {
        auto cubeListIsEmpty = [&](const TriggerPattern& triggerPattern){
          auto &cubeList = _patternCaloMap->ConstructRetrieve(triggerPattern, *_fitGeoParamsColl, *_caloGeoParams);
          return (cubeList.NTotalCubes()==0);
        };
        auto it = std::remove_if(vTriggers.begin(),vTriggers.end(),cubeListIsEmpty);
        vTriggers.erase(it,vTriggers.end());
      }


      for (auto triggerPattern : vTriggers) {

        // auto roi = UlegGeomRoi(triggerPattern,*_FitGeoParamsColl);
        auto &cubeList = _patternCaloMap->ConstructRetrieve(triggerPattern, *_fitGeoParamsColl, *_caloGeoParams);
        auto &tileList = _patternPsdMap->ConstructRetrieve(triggerPattern, *_fitGeoParamsColl, *_psdGeoParamsColl);
    }


      triggersColl[side]=vTriggers;
      // totalMultiplicity+=vTriggers.size();
    } //loop over sides

    if (_dumpPatterns) {

    COUT(INFO) << "ULEG Pattern -> tiles in RoI dump: " << ENDL;
    for (auto fitSide : Herd::RefFrame::Sides) {
      std::cout << "#################################################################################" << '\n';
      std::cout << "FIT side: " << fitSide << "; " << triggersColl[fitSide].size() << " patterns" << '\n';
      for (auto triggerPattern : triggersColl[fitSide]) {
        std::cout << "Trigger pattern: " << triggerPattern.FullID() << '\n';
        auto &cubeList = _patternCaloMap->ConstructRetrieve(triggerPattern, *_fitGeoParamsColl, *_caloGeoParams);
        auto &tileList = _patternPsdMap->ConstructRetrieve(triggerPattern, *_fitGeoParamsColl, *_psdGeoParamsColl);
        std::cout << tileList;
        // std::cout << cubeList << '\n';
        std::cout << "Cubes in RoI: " << cubeList.NTotalCubes() << '\n';
        std::cout << '\n';
        // if (cubeList.NTotalCubes()==0) {
        //   std::cout << "EMPTY EMPTY EMPTY EMPTY" << '\n';
        // }
      }
    }
  }

  if (_storePatternToTileList) {
    for (auto fitSide : Herd::RefFrame::Sides) {
      auto fitGeoParamsVec = (*_fitGeoParamsColl)[fitSide];
      for (auto triggerPattern : triggersColl[fitSide]) {
        auto fitLayer = triggerPattern.Layer();
        auto mats = triggerPattern.Mats();
        auto fitGeoParams = fitGeoParamsVec[fitLayer];
        auto posFitFirstMat = fitGeoParams.MatPosition(mats[0]);
        // std::cout << "ping" << '\n';
        auto normalUnitVec = Herd::Vec3D(fitGeoParams.NormalDirection());
        // std::cout << "pong" << '\n';
        auto transUnitVec = Herd::Vec3D(fitGeoParams.SegmentationDirections().first);
        auto &tileList = _patternPsdMap->ConstructRetrieve(triggerPattern, *_fitGeoParamsColl, *_psdGeoParamsColl);
        auto &cubeList = _patternCaloMap->ConstructRetrieve(triggerPattern, *_fitGeoParamsColl, *_caloGeoParams);
        // std::cout << "pang" << '\n';
        auto nCubes = cubeList.NTotalCubes();
        for (auto iPsdSide : Herd::RefFrame::Sides) {
          if (tileList.TilesPerSide()[iPsdSide]>0){
          auto psdGeoParamsVector = (*_psdGeoParamsColl)[iPsdSide];
          auto tilesSide = tileList.TileIDs()[iPsdSide];
          for (size_t iPsdLayer = 0; iPsdLayer<tilesSide.size(); iPsdLayer++) {
            auto psdGeoParams = psdGeoParamsVector[iPsdLayer];
            // std::cout << "peng" << '\n';
            // std::cout << tileList.TilesPerSide()[iPsdSide] << '\n';
            auto tileIDs = tilesSide[iPsdLayer];

            if (tileIDs.size()>0) {
              for (auto tileID : tileIDs) {

              // for (auto psdTile : tileList.TileIDs()[iPsdSide][iPsdLayer]) {
                // std::cout << tileID << '\n';
                auto posPsdTile = psdGeoParams.Position(tileID);
                auto psdTileRelativePos = posPsdTile - posFitFirstMat;
                // std::cout << "pung" << '\n';
                auto tileAngleFromMat = std::atan2(psdTileRelativePos*transUnitVec,psdTileRelativePos*normalUnitVec);
                // std::cout << "pung" << '\n';
                _patternsAndTiles->Fill(static_cast<int>(fitSide),fitLayer,mats[0],mats[1],mats[2],static_cast<int>(iPsdSide),iPsdLayer,tileID,tileAngleFromMat,nCubes);
                // std::cout << "pung" << '\n';


              }
            }
            }
          }
          /* code */

        }

      }
    }
  }
    _globStore->AddObject("patternPsdMap",_patternPsdMap);
    _globStore->AddObject("patternCaloMap",_patternCaloMap);
    _globStore->AddObject(_patternsAndTiles->GetName(),_patternsAndTiles);

  return true;
}
