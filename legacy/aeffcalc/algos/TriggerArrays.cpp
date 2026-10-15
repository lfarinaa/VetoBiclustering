/*
* TriggerArrays.cpp
*
*  Created on: 12 Oct 2019
*      Author: Valerio Formato
*/

// Example headers
#include "TriggerArrays.h"
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
#include <chrono>

RegisterAlgorithm(TriggerArrays);

TriggerArrays::TriggerArrays(const std::string &name) : Algorithm{name}, _threshold{0}, _maxDelta{1}, _patternSize{3}, _nHits{3}, _layer_start{0}, _force_full_spanning{false}, _ignore_pattern_IDs{{}}, _require_only_one_side{false}, _matHitsCollName{"fitMatHitsCollMC"}, _store_patterns_name{""}, _suppress_patterns_name{""}, _require_calo_in_roi{false},
_sidesToStudy{}, _removeUpwardsPatterns{false}, _removeAlmostUpwardsPatterns{false} {
  DefineParameter("threshold", _threshold);
  DefineParameter("maxDelta", _maxDelta);
  DefineParameter("nHits", _nHits);
  DefineParameter("patternSize", _patternSize);
  DefineParameter("force_full_spanning", _force_full_spanning);
  // DefineParameter("disable_filter", _disable_filter);
  DefineParameter("layer_start", _layer_start);
  // DefineParameter("ignore_plusminus", _ignore_plusminus);
  DefineParameter("ignore_pattern_IDs", _ignore_pattern_IDs);
  DefineParameter("sidesToStudy", _sidesToStudy);
  DefineParameter("require_only_one_side", _require_only_one_side);
  DefineParameter("matHitsCollName", _matHitsCollName);
  DefineParameter("store_patterns_name", _store_patterns_name);
  DefineParameter("suppress_patterns_name", _suppress_patterns_name);
  DefineParameter("require_calo_in_roi", _require_calo_in_roi);
  DefineParameter("removeUpwardsPatterns", _removeUpwardsPatterns);
  DefineParameter("removeAlmostUpwardsPatterns", _removeAlmostUpwardsPatterns);


}

bool TriggerArrays::Initialize() {
  const std::string routineName("TriggerArrays::Initialize");

  try {
    _globStore = GetDataStoreManager()->GetGlobalDataStore("globStore");
    _FitGeoParamsColl = _globStore->GetObject<Herd::FitGeoParamsColl>("fitGeoParamsColl");
    _PsdGeoParamsColl = _globStore->GetObject<Herd::PsdGeoParamsColl>("psdGeoParamsColl");
    // _ScdGeoParamsColl = _globStore->GetObject<Herd::ScdGeoParamsColl>("scdGeoParamsColl");
    caloGeoParams = _globStore->GetObject<Herd::CaloGeoParams>("caloGeoParams");

    _evStore = GetDataStoreManager()->GetEventDataStore("evStore");
    if (!_globStore->IsPresent("patternCaloMap")) {
      auto _patternCaloMapToBeAdded = std::make_shared<PatternVetoListMap<CaloVetoList,Herd::CaloGeoParams>>(PatternVetoListMap<CaloVetoList,Herd::CaloGeoParams>());
      _globStore->AddObject("patternCaloMap",_patternCaloMapToBeAdded);
      // } else {
    } else {
      // std::cout << "Found pattern map in algo " + GetName() << '\n';
    }


    _patternCaloMap = _globStore->GetObject<PatternVetoListMap<CaloVetoList,Herd::CaloGeoParams>>("patternCaloMap");
  }
  catch (Retrieval::Exception &exc) {
    COUT(ERROR) << "Event data store not found." << ENDL;
    return false;
  }
  // _treeTriggers = std::make_shared<TTree>(std::string("TriggerMats"+GetName()).c_str(), "Trigger signals");
  // _tupleTriggersEven = std::make_shared<TNtuple>(std::string("TriggersEven_"+GetName()).c_str(), "Trigger patterns", "eTrue:cosTheta:iLayer:iMat:patternID:oldID");
  _tupleFirstTriggers = std::make_shared<TNtuple>(std::string("TriggerArrays"+GetName()+"FirstTriggers").c_str(), "Trigger patterns", "eTrue:cosPolar:azimuth:iSideSelected:cosIncidenceAngle:iSide:iLayer:iMat:patternID:oldID:eDepPsd:maxECalo:totalECalo:nHitsPsd:width");
  // _tupleFirstTriggersEven = std::make_shared<TNtuple>(std::string("FirstTriggersEven_"+GetName()).c_str(), "Trigger patterns", "eTrue:cosTheta:iLayer:iMat:patternID:oldID");
  _tupleTriggers = std::make_shared<TNtuple>(std::string("TriggerArrays"+GetName()+"AllTriggers").c_str(), "Trigger patterns", "eTrue:cosTheta:iSide:iLayer:iMat:patternID:oldID:maxECalo:totalECalo:nHitsPsd:totMult");
  _tupleTriggerEvents = std::make_shared<TNtuple>(std::string("TriggerArrays"+GetName()+"Events").c_str(), "Trigger patterns", "eTrue:cosPolar:azimuth:iSideSelected:triggeredSides:maxECalo:totalECalo:nHitsPsd:totalMultiplicity");
  _tupleSides = std::make_shared<TNtuple>(std::string("TriggerArrays"+GetName()+"Sides").c_str(), "Trigger patterns", "eTrue:cosPolar:azimuth:iSideSelected:cosIncidenceAngle:iSide:triggeredSides:nPatternsInSide:EPsd:maxECalo:totalECalo:nHitsPsd");
  // _tupleFirstLayer = std::make_shared<TNtuple>(std::string("FirstLayer_"+GetName()).c_str(), "Trigger patterns", "eTrue:cosTheta:iLayer:iMat:patternID:oldID:xPSD0:yPSD0:zPSD0:x:y:z:px:py:pz");
  // _tupleTriggersOdd = std::make_shared<TNtuple>(std::string("TriggersOdd_"+GetName()).c_str(), "Trigger patterns", "eTrue:cosTheta:iLayer:iMat:patternID:oldID");
  // _tupleFirstTriggersOdd = std::make_shared<TNtuple>(std::string("FirstTriggersOdd_"+GetName()).c_str(), "Trigger patterns", "eTrue:cosTheta:iLayer:iMat:patternID:oldID");
  // _tupleMats = std::make_shared<TNtuple>(std::string("Mats_"+GetName()).c_str(), "Triggered mats", "eTrue:cosTheta:iLayer:iMat:patternID:oldID:xPSD0:yPSD0:zPSD0:x:y:z:px:py:pz");
  // _tupleFirstMats = std::make_shared<TNtuple>(std::string("FirstMats_"+GetName()).c_str(), "Triggered mats", "eTrue:cosTheta:iLayer:iMat:patternID:oldID");

  // Setup the filter
  // SetFilterStatus(FilterStatus::ENABLED); // don't set up for now
  _nAccept = _nReject = 0;


  // if (_restrictSidesToConversionSide){
  //   auto conversionSide = *(_evStore->GetObject<Herd::Side>("conversionSide"));
  //   _sides = std::vector<Herd::RefFrame::Side>(conversionSide);
  // }

  return true;
}

bool TriggerArrays::Process() {
  const std::string routineName("TriggerArrays::Process");
  try {
    // auto psdHitsColl = _evStore->GetObject<Herd::PsdHitsColl>("psdHitsCollMC");
    auto fitMatHitsColl = _evStore->GetObject<Herd::FitHitsColl>(_matHitsCollName);
    if (!fitMatHitsColl) {
      COUT(ERROR) << "Fit mat hits not present; need to generate at digitization time!" << GetEventLoopProxy()->GetCurrentEvent() << ENDL;
      return false;
    }
    SetFilterResult(FilterResult::REJECT);

    // auto side_selected = _evStore->GetObject<Herd::RefFrame::Side>("SelectedFace");
    auto _sides = std::vector<Herd::RefFrame::Side>();

    if (_sidesToStudy.empty()) {
      if (_evStore->IsPresent("sidesToStudy")) {
        _sides = *(_evStore->GetObject<std::vector<Herd::RefFrame::Side>>("sidesToStudy"));
      } else {
        _sides = std::vector<Herd::RefFrame::Side>(Herd::RefFrame::Sides.begin(),Herd::RefFrame::Sides.end());
      }
    } else {
      for (auto sideNumber : _sidesToStudy) {
        _sides.push_back(Herd::RefFrame::Side(sideNumber));
      }
    }
    // std::cout << _sides.size() << '\n';

    // retrieve MC truth info for studies:
    auto mcTruth = _evStore->GetObject<Herd::MCTruth>("mcTruth");
    const auto &primaryPart = mcTruth->primaries.at(0);
    auto initPos = primaryPart.InitialPosition();
    auto momentum = primaryPart.InitialMomentum();
    auto trajectory = primaryPart.InitialTrajectory();
    auto eTrue = momentum.Mag();
    auto MCInteraction = primaryPart.PrimaryDisappearance();
    auto azimuth = momentum.AzimuthAngle();
    Herd::SidesArray<float> cosIncidenceAngles;
    cosIncidenceAngles.fill(-2);
    for (auto side : Herd::RefFrame::Sides) {
      auto normalAxis = Herd::RefFrame::ToAxis(Herd::RefFrame::NormalDirection(side));
      auto cosIncidenceAngle=momentum[Herd::RefFrame::ToCoo(normalAxis)]/eTrue;
      cosIncidenceAngles[side]=cosIncidenceAngle;
    }
    auto &cosPolar = cosIncidenceAngles[Herd::RefFrame::Side::TOP];

    auto selectedSide = _evStore->GetObject<Herd::RefFrame::Side>("selectedSide");
    auto iSideSelected = static_cast<int>(*selectedSide);


    // auto x=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000;

    // Get PSD info for studies:
    Herd::Point crossPsd = initPos;
    float distToOrigin=0;
    // float minDist=1e3;
    bool foundPsdCrossing=false;
    Herd::RefFrame::Side crossingSide=Herd::RefFrame::Side::NONE;
    for (auto side : Herd::RefFrame::Sides) {
      const auto &psdGeoParamsVector = (*_PsdGeoParamsColl)[side];
      auto normalAxisToSide = Herd::RefFrame::ToAxis(Herd::RefFrame::NormalDirection(side));
      Herd::Vec3D normalVector = psdGeoParamsVector.front().NormalDirection();
      if (momentum*normalVector<0) {
        for (size_t iLayer=0;iLayer<psdGeoParamsVector.size();iLayer++) {
          const auto &psdGeoParams = psdGeoParamsVector[iLayer];
          auto nTiles = psdGeoParams.NTotElements();
          for (size_t iTile =0; iTile<nTiles;iTile++) {
            tilePrism=psdGeoParams.ElementShape(iTile);
            maybeIntersections=tilePrism.Intersection();
            foundPsdCrossing=foundPsdCrossing||maybeIntersection.has_values();
          }
          // auto normalAxis = Herd::RefFrame::ToAxis(psdGeoParams.NormalDirection());
          // float normalPos = psdGeoParams.Position(0)[Herd::RefFrame::ToCoo(normalAxis)];
          // try {
          //   if (auto tempCrossPsd = trajectory.Intersection(normalAxis,normalPos)){

          //     auto distThisCrossing = (tempCrossPsd.value()-initPos).Mag();
          //     if (distThisCrossing>distToOrigin && tempCrossPsd.value().Mag()<crossPsd.Mag() && side!=crossingSide) {
          //       // very bad. Remind Nicola to add dimensions to geoParams objects
          //       distToOrigin=distThisCrossing;
          //       crossingSide=side;
          //       crossPsd=tempCrossPsd.value();
          //       foundPsdCrossing=true;
          //     }
          //   }
          // } catch (Herd::Line::Degenerate &exc) {
          //   continue; // why do I even need to be doing this...
          // }
          // for (size_t iTile = 0; iTile < psdGeoParams.NTotElements(); iTile++) {
          //   auto distThisSide = (tempCrossPsd-initPos).Mag(); // Distance to primary origin
          //   // if ((tempCrossPsd-psdGeoParams.Position(iTile)).Mag()<minDist) {
          //
          //   // }
          //   minDist=std::min((tempCrossPsd-psdGeoParams.Position(iTile)).Mag(),minDist);
          //   if ((tempCrossPsd-psdGeoParams.Position(iTile)).Mag()<10) {
          //     if (distThisSide<distToOrigin) {
          //       distToOrigin=distThisSide;
          //       crossPsd=tempCrossPsd;
          //       foundPsdCrossing=true;
          //     }
          //   }
          // }
        }
      }
    } // loop over sides


    ////////////////////////////////////////////////////////////////////////////////
    // Calo information for studies

    auto caloHits = _evStore->GetObject<Herd::CaloHits>("caloHitsMC");
    // auto caloGeoParams = _globStore->GetObject<Herd::CaloGeoParams>("caloGeoParams");

    std::array<float,6> eCalo = {0.};
    for (auto &caloHit : *caloHits) { //loop over hits in CALO, assign EDep to the corresponding caloRegion
      auto volID=caloHit.VolumeID();
      auto hitIndexes=caloGeoParams->CooIndexes(volID); // note that CooIndexes returns a uint
      // also the indexes increase towards the bottom of the CALO
      int xReduced=hitIndexes[Herd::RefFrame::Coo::X]-10; // remap the 21 cubes 0..20 to -10..10
      int yReduced=hitIndexes[Herd::RefFrame::Coo::Y]-10;
      int zReduced=hitIndexes[Herd::RefFrame::Coo::Z]-10;

      int caloRegion;
      if (std::abs(zReduced)>=8) {
        if (zReduced<0){
          caloRegion = 2; //zTop
        } else {
          caloRegion = 1; //zBot
        }
      } else if (std::abs(xReduced)<=7 && std::abs(yReduced)<=7 && (std::abs(xReduced)+std::abs(yReduced))<=11){
        caloRegion = 0; //core
      } else if (xReduced>=4) {
        caloRegion = 3; //xPos
      } else if (yReduced>=0) {
        caloRegion = 4; //yPos
      } else {
        caloRegion = 5; //yNeg
      }
      auto hitEDep = caloHit.EDep();
      eCalo[caloRegion]+=hitEDep;
    }
    auto maxECalo=*std::max_element(eCalo.begin()+2,eCalo.end()); //Skips the first 2 elements of ECalo, corresponding to the core+bot regions
    auto totalECalo=std::accumulate(eCalo.begin(),eCalo.end(),0.);

    ////////////////////////////////////////////////////////////////////////////////
    // Psd information for studies

    auto psdHitsColl = _evStore->GetObject<Herd::PsdHitsColl>("psdHitsCollMC");
    float EPsd=0.;
    int nHitsPsd=0;
    // auto side=*side_selected;
    for (auto side : Herd::RefFrame::Sides) {
      for (auto iLayer = 0; iLayer < (*_PsdGeoParamsColl)[side].size(); ++iLayer) {
        for (const auto &hit : (*psdHitsColl)[side][iLayer]) {
          EPsd += hit.EDep();
          nHitsPsd++;
        }
      }
    }

    ////////////////////////////////////
    // Suppress patterns found in a TriggerPattern object in the evStore
    Herd::SidesArray<std::vector<TriggerPattern>> *suppressPatterns;
    if (!_suppress_patterns_name.empty()) {
      suppressPatterns = _evStore->GetObject<Herd::SidesArray<std::vector<TriggerPattern>>>(_suppress_patterns_name);
    }

    Herd::SidesArray<std::vector<TriggerPattern>> triggersColl;
    int totalMultiplicity=0;
    Herd::SidesArray<int> triggerSides;
    triggerSides.fill(0);

    // std::cout << "algo: " << GetName() << '\n';
    for (auto side : _sides) {
      // std::cout << "side: " << side << '\n';

      auto &fitGeoParamsVector = (*_FitGeoParamsColl)[side];
      // auto &fitGeoParamsVector = (*_FitGeoParamsColl)[*side_selected];
      auto &fitMatHitsSide = (*fitMatHitsColl)[side];
      // auto &fitMatHitsSide = (*fitMatHitsColl)[*side_selected];
      auto nLayers = fitGeoParamsVector.size();
      auto geoParamsEven = fitGeoParamsVector[0];
      auto geoParamsOdd = fitGeoParamsVector[1];
      auto primSegAxis = Herd::RefFrame::ToAxis(geoParamsEven.SegmentationDirections().first);
      auto secSegAxis = Herd::RefFrame::ToAxis(geoParamsOdd.SegmentationDirections().first);
      auto nMatsEven = geoParamsEven.NMats()[primSegAxis];
      auto nMatsOdd = geoParamsOdd.NMats()[secSegAxis];

      std::vector<int> nMatsPair{nMatsEven,nMatsOdd};
      FitHitsArray hitsArrayEven(nLayers/2, std::vector<int>(nMatsEven,0)); // init the hit arrays to 0
      FitHitsArray hitsArrayOdd(nLayers/2, std::vector<int>(nMatsOdd,0)); // init the hit arrays to 0
      auto hitsArraysPair = std::vector<FitHitsArray>{hitsArrayEven,hitsArrayOdd};
      for (size_t iLayer = 0; iLayer < nLayers; iLayer++) {
        // std::cout << iLayer << "\n";
        auto &fitMatHitsLayer = fitMatHitsSide[iLayer];
        bool parity = (iLayer%2);
        for (auto &hit : fitMatHitsLayer) { // fill both arrays by looking at the mats data
          if (hit.EDep()>_threshold) { //threshold is applied at mat level. To do: apply it at SiPM level, or wherever it should be applied
            auto iMat = hit.VolumeID();
            hitsArraysPair[parity][iLayer/2][iMat]=1;
          }
          // the line above uses dirty int-division trick. needs explicit -1?
        } // loop over mats
      } // loop over layers
      // todo: find better way of handling this codint pattern
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
          if (hasHit || _patternSize > _nHits) {
            auto iPattern = TriggerPattern(side,iLayer,iMat,!hasHit);
            // TriggerPattern iPattern = TriggerPattern(side,iLayer,iMat,!hasHit};
            iPattern.FindPatterns(_nHits-hasHit,_patternSize,_maxDelta,hitsArray,vTriggers,_force_full_spanning);
          }
        } // loop over mats
      } // loop over layers
      // std::cout << side << '\n';
      // std::cout << "found patterns: " << vTriggers.size() << '\n';
      // for (auto patt : vTriggers) std::cout << patt << '\n';
      std::sort(vTriggers.begin(),vTriggers.end()); // sort triggers, just in case
      // remove duplicated patterns, if there are any:
      vTriggers.erase(std::unique(vTriggers.begin(),vTriggers.end(),[=](auto &p1, auto &p2){return p1.IsEquivalent(p2,_patternSize);}),vTriggers.end()); // unique-erase: erase duplicates; needed if allowing holes
      if (!_suppress_patterns_name.empty()) { // pattern suppression by collection
        auto it = std::remove_if(vTriggers.begin(),vTriggers.end(),[suppressSide=(*suppressPatterns)[side]](auto &pattern){return std::any_of(suppressSide.begin(),suppressSide.end(),[&pattern](auto &pattToSuppress){return pattToSuppress==pattern;});});
        vTriggers.erase(it,vTriggers.end()); // delete matches in other container
      }
      if (!_ignore_pattern_IDs.empty()) { // pattern suppression by ID, loaded in the datacard
        // prune patters if their IDs match a member of _ignore_pattern_IDs
        auto it = std::remove_if(vTriggers.begin(),vTriggers.end(),[ignore_pattern_IDs=_ignore_pattern_IDs](auto &pattern) {return std::any_of(ignore_pattern_IDs.begin(),ignore_pattern_IDs.end(),[&pattern](auto &id){return id==pattern.PatternID();});});
        vTriggers.erase(it,vTriggers.end());
      }

      if (_require_calo_in_roi) {
        // auto x=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000;
        // std::cout << "a " << side << (unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000 << " " << GetName() << '\n';

        auto cubeListIsEmpty = [&](const TriggerPattern& triggerPattern){


          // auto x=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000;


          // auto geomRoi = UlegGeomRoi(triggerPattern,*_FitGeoParamsColl);
          // auto hasCubes=geomRoi.MinimumCubes(*caloGeoParams,1);

          auto &caloList = _patternCaloMap->ConstructRetrieve(triggerPattern, *_FitGeoParamsColl, *caloGeoParams);

          // auto y=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000-x;
          // std::cout << y << '\n';
          // return !hasCubes;

          return (caloList.NTotalCubes()==0);
        };
        auto it = std::remove_if(vTriggers.begin(),vTriggers.end(),cubeListIsEmpty);
        vTriggers.erase(it,vTriggers.end());
        // auto y=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000-x;
        // if (y>0) {
        // std::cout << side << " " << y << " " << vTriggers.size() << '\n';
        // for (auto triggerPattern : vTriggers) std::cout << triggerPattern.PatternID() << '\n';
        // }


      }
      // std::cout << "found patterns: " << vTriggers.size() << '\n';
      if (_removeUpwardsPatterns) {
        auto isUpwardsPointing = [&fitGeoParamsColl=*_FitGeoParamsColl](const TriggerPattern& triggerPattern){return triggerPattern.PointsUpwards(fitGeoParamsColl);};
        auto it = std::remove_if(vTriggers.begin(),vTriggers.end(),isUpwardsPointing);
        vTriggers.erase(it,vTriggers.end());
      }
      // std::cout << "found patterns: " << vTriggers.size() << '\n';

      if (_removeAlmostUpwardsPatterns) {
        auto isAlmostUpwardsPointing = [&fitGeoParamsColl=*_FitGeoParamsColl](const TriggerPattern& triggerPattern){return triggerPattern.PointsAlmostUpwards(fitGeoParamsColl);};
        auto it = std::remove_if(vTriggers.begin(),vTriggers.end(),isAlmostUpwardsPointing);
        vTriggers.erase(it,vTriggers.end());
      }

      // cross-match triggers and store which layers are active (have a matching trigger in the other direction)
      std::vector<int> vTriggeredLayers;
      // for (auto layer : vTriggeredLayers) std::cout << layer << '\n';
      for (size_t i0=0; i0 < vTriggers.size(); i0++) {
        if (vTriggers[i0].Layer()%2==0) {
          // if this trigger is placed on a (side,layer) already detected as active, skip
          if (std::find(vTriggeredLayers.begin(),vTriggeredLayers.end(),vTriggers[i0].Layer())!=vTriggeredLayers.end()) continue;
          for (size_t i1=i0+1; i1 < vTriggers.size(); i1++) {
            // if this trigger is placed on a (side,layer) already detected as active, skip
            if (std::find(vTriggeredLayers.begin(),vTriggeredLayers.end(),vTriggers[i1].Layer())!=vTriggeredLayers.end()) continue;
            if (vTriggers[i0].Layer()==vTriggers[i1].Layer()-1) { // if the layers belong to the same double layer

              // int nCubes=0;
              bool hasCubes = false;
              if (_require_calo_in_roi) {
                auto &caloList0 = _patternCaloMap->ConstructRetrieve(vTriggers[i0], *_FitGeoParamsColl, *caloGeoParams);
                auto &caloList1 = _patternCaloMap->ConstructRetrieve(vTriggers[i1], *_FitGeoParamsColl, *caloGeoParams);

                // auto intersection = caloList0.Intersection(caloList1);
                // nCubes = intersection.NTotalCubes();
                hasCubes = caloList0.IntersectionHasCubes(caloList1,1);
                // if (nCubes==0) std::cout << "it happened" << hasCubes << '\n';
              }
              if (!_require_calo_in_roi || hasCubes) {
                vTriggeredLayers.push_back(vTriggers[i0].Layer());
                vTriggeredLayers.push_back(vTriggers[i1].Layer());
                break;
              }
            }
          }
        }
      }
      // for (auto layer : vTriggeredLayers) std::cout << layer << '\n';
      // if (!vTriggeredLayers.empty()) {
      // prune patterns in unmatched layers using remove_if-erase:
      auto it = std::remove_if(vTriggers.begin(),vTriggers.end(),[&vTriggeredLayers](auto &pattern) {return std::none_of(vTriggeredLayers.begin(),vTriggeredLayers.end(),[&pattern](auto &layer){return layer==pattern.Layer();});});
      vTriggers.erase(it,vTriggers.end());
      // } else {
      // vTriggers.clear();
      // }
      if (!vTriggeredLayers.empty()) {
        triggerSides[side]+=1;
        
      } else {
        // std::cout << "No triggered layers in side " << side << std::endl;
      }
      triggersColl[side]=vTriggers;
      totalMultiplicity+=vTriggers.size();
    } //loop over sides
    // for (auto side : _sides) {
    // std::cout << side << '\n';
    //   for (auto patt : triggersColl[side]) std::cout << patt << '\n';
    //
    // }
    auto triggeredSides=std::accumulate(triggerSides.begin(),triggerSides.end(),0);
    // store all patterns, to be read by another algorithm
    if (!_store_patterns_name.empty()) {
      auto triggerPatterns = std::make_unique<Herd::SidesArray<std::vector<TriggerPattern>>>();
      *triggerPatterns = triggersColl;
      // the patterns are stored post-pruning. I will fix this someday, maybe. It's only relevant if you allow for holes in the trigger patterns
      _evStore->AddObject(_store_patterns_name, std::move(triggerPatterns));
    }

    if (_require_only_one_side && triggeredSides>1) {
      _nReject++;
      SetFilterResult(FilterResult::REJECT);
      goto end;
    }
    if (triggeredSides>0) {
    // if (triggersColl[*side_selected].size()>0) {
      _nAccept++;
      SetFilterResult(FilterResult::ACCEPT);

      auto firstTriggers = std::make_unique<Herd::SidesArray<std::vector<TriggerPattern>>>();
      for (auto side : _sides) {
        bool firstEven = true;
        bool firstOdd = true;
        if (triggersColl[side].size()>0) {
          _tupleSides->Fill(eTrue,cosPolar,azimuth,iSideSelected,cosIncidenceAngles[side],static_cast<int>(side),triggeredSides,triggersColl[side].size(),EPsd,maxECalo,totalECalo,nHitsPsd);

          // if (iSideSelected!=static_cast<int>(side)) std::cout << "just happened" << '\n';
          // std::cout << _sides.size() << '\n';
          // std::cout << iSideSelected << '\n';
          // std::cout << (int)side << '\n';
        for (auto pattern : triggersColl[side]) {
          // std::cout << pattern << '\n';
          _tupleTriggers->Fill(eTrue,cosIncidenceAngles[side],static_cast<int>(pattern.Side()),pattern.Layer(),pattern.Mats().front(),std::stoi(pattern.PatternID()),pattern.OldID(),maxECalo,totalECalo,nHitsPsd,totalMultiplicity);
          // if (side==*side_selected) {
            if (firstEven && !(pattern.Layer()%2)) {
              _tupleFirstTriggers->Fill(eTrue,cosPolar,azimuth,iSideSelected,cosIncidenceAngles[side],static_cast<int>(pattern.Side()),pattern.Layer(),pattern.Mats().front(),std::stoi(pattern.PatternID()),pattern.OldID(),EPsd,maxECalo,totalECalo,nHitsPsd,pattern.Width());
              firstEven = false;
              // if (GetFilterStatus()==FilterStatus::ENABLED) {
              // auto triggerPatternEven = std::make_unique<TriggerPattern>(pattern);
              (*firstTriggers)[side].push_back(pattern);

              // _evStore->AddObject("evenTrigger"+GetName(), std::move(triggerPatternEven));
              // }
            } else if (firstOdd && pattern.Layer()%2) {
              _tupleFirstTriggers->Fill(eTrue,cosPolar,azimuth,iSideSelected,cosIncidenceAngles[side],static_cast<int>(pattern.Side()),pattern.Layer(),pattern.Mats().front(),std::stoi(pattern.PatternID()),pattern.OldID(),EPsd,maxECalo,totalECalo,nHitsPsd,pattern.Width());
              firstOdd = false;
              (*firstTriggers)[side].push_back(pattern);
              // if (GetFilterStatus()==FilterStatus::ENABLED) {
              // auto triggerPatternOdd = std::make_unique<TriggerPattern>(pattern);
              // }
            }
          // }
        }
      }
      }
      _evStore->AddObject("firstTriggers"+GetName(), std::move(firstTriggers));

      // if (GetFilterStatus()==FilterStatus::ENABLED) {
      auto triggersFirstLayers = std::make_unique<Herd::SidesArray<std::vector<std::vector<TriggerPattern>>>>();
      for (auto side : _sides) {
        if (triggersColl[side].size()>0) {
          auto firstTriggeredLayer=triggersColl[side].front().Layer(); // get first triggered layer in this side, if any
          (*triggersFirstLayers)[side].push_back({}); // init the things even if there are no elements
          (*triggersFirstLayers)[side].push_back({});
          for (auto pattern : triggersColl[side]) {
            auto diff=pattern.Layer()-firstTriggeredLayer;
            if (diff<2) { // for the first and second layers, put the triggers in this object:
              (*triggersFirstLayers)[side][diff].push_back(pattern);
            }
          }
        }
      }
      // publish the triggers in the first valid layer
      _evStore->AddObject("triggersFirstLayers"+GetName(), std::move(triggersFirstLayers));
      // }
      _tupleTriggerEvents->Fill(eTrue,cosPolar,azimuth,iSideSelected,triggeredSides,maxECalo,totalECalo,nHitsPsd,totalMultiplicity);


    } else {
      _nReject++;
      SetFilterResult(FilterResult::REJECT);
    }


    // to do:
    // - clarify where layers are 0-13 and where 0-6. maybe move to objects 0-13 with variable length
    // - maintain layer/mat structure in the triggers vector prior to writing to tuples
    // (for readability reasons)
    // - improve printed diagnostic on algorithm end (?)
    // auto y=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000-x;

    // std::cout << GetName()+": " << y << '\n';

  }
  catch (Retrieval::Exception &exc) {
    COUT(ERROR) << exc.what() << ENDL;
    return false;
  }
  end:
  // if (GetFilterResult()==FilterResult::REJECT) _nReject++;
  // if (_disable_filter) SetFilterResult(FilterResult::ACCEPT); // only for testing purposes
  return true;
}

bool TriggerArrays::Finalize() {
  const std::string routineName("TriggerArrays::Finalize");

  _globStore->AddObject(_tupleTriggers->GetName(), _tupleTriggers);
  _globStore->AddObject(_tupleTriggerEvents->GetName(), _tupleTriggerEvents);
  _globStore->AddObject(_tupleSides->GetName(), _tupleSides);
  _globStore->AddObject(_tupleFirstTriggers->GetName(), _tupleFirstTriggers);

  COUT(INFO) << "Selected events: " << _nAccept << "/" << _nReject + _nAccept << " ("
  << (float)_nAccept / (float)(_nReject + _nAccept) * 100. << "\%)" << ENDL;
  // if (_force_full_spanning) COUT(INFO) << "Patterns fully spanning" << ENDL;

  return true;
}
