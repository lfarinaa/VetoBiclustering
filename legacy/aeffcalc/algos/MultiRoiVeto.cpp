/*
* MultiRoiVeto.cpp
*
*  Created on: 11 Feb 2019
*      Author: Nicola Mori
*/

// Reject the event that depose less than 50% of photon energy in the calorimeter

// Example headers
#include "MultiRoiVeto.h"
#include "TriggerPattern.h"
#include "UlegGeomPsdRoi.h"


// HerdSoftware headers
#include "herd/analysis/dataobjects/CaloHits.h"
#include "herd/analysis/dataobjects/MCTruth.h"
#include "herd/analysis/dataobjects/CaloGeoParams.h"
#include "herd/analysis/dataobjects/FitGeoParamsColl.h"
#include "herd/analysis/dataobjects/PsdGeoParamsColl.h"
#include "herd/analysis/dataobjects/FitParticleHitsColl.h"
#include "herd/analysis/dataobjects/PsdHitsColl.h"
#include "herd/analysis/dataobjects/Parallelogram.h"

#include "herd/analysis/dataobjects/TrackInfoForCalo.h"


// EA headers
#include "utils/SmartLog.h"
// C/C++ standard headers
#include <numeric>
#include <cmath>

//Root headers
#include "TNtuple.h"

RegisterAlgorithm(MultiRoiVeto);

MultiRoiVeto::MultiRoiVeto(const std::string &name) : _threshold{0.333e-3}, _fullPsdCheckTh{1e10}, _triggerAlgoName{GetName()}, Algorithm{name},
_sidesToStudy{} {

  DefineParameter("threshold", _threshold);
  DefineParameter("fullPsdCheckTh", _fullPsdCheckTh);
  // DefineParameter("sidesToStudy", _sidesToStudy);
  DefineParameter("triggerAlgoName", _triggerAlgoName);
}

bool MultiRoiVeto::Initialize() {
  const std::string routineName("MultiRoiVeto::Initialize");

  try {
    _evStore = GetDataStoreManager()->GetEventDataStore("evStore");
    _globStore = GetDataStoreManager()->GetGlobalDataStore("globStore");
    if (!_globStore->IsPresent("patternPsdMap")) {
      auto _patternPsdMapToBeAdded = std::make_shared<PatternVetoListMap<PsdVetoList,Herd::PsdGeoParamsColl>>(PatternVetoListMap<PsdVetoList,Herd::PsdGeoParamsColl>());
      _globStore->AddObject("patternPsdMap",_patternPsdMapToBeAdded);
      // } else {
    }
    _patternPsdMap = _globStore->GetObject<PatternVetoListMap<PsdVetoList,Herd::PsdGeoParamsColl>>("patternPsdMap");
  }
  catch (Retrieval::Exception &exc) {
    COUT(ERROR) << "Event data store not found." << ENDL;
    return false;
  }

  // // Setup the filter
  // SetFilterStatus(FilterStatus::ENABLED);
  _nAccept = _nReject = 0;

  // _tupleE=std::make_shared<TNtuple>(GetName().c_str(),"temp", "eTrue:cosTheta:maxERoI:iLayer:nVetoes");
  _tupleVetoEvents=std::make_shared<TNtuple>(("MultiRoiVeto"+GetName()+"Events").c_str(),"temp", "eTrue:cosPolar:azimuth:iSideSelected:maxERoI:sidesPassingRoi:totalECalo:disHeight:allSidesClean:passFullPsdVeto");
  _tupleVetoSides=std::make_shared<TNtuple>(("MultiRoiVeto"+GetName()+"Sides").c_str(),"temp", "eTrue:cosPolar:azimuth:cosIncidenceAngle:iSideSelected:maxERoI:maxEPsd:totalCombinations:iSide:nTiles:totalECalo");
  // _tupleAllHits=std::make_shared<TNtuple>((GetName()+"AllHits").c_str(),"temp", "eTrue:cosTheta:eDep:iLayer");
  return true;
}

bool MultiRoiVeto::Process() {
  const std::string routineName("MultiRoiVeto::Process");
  try {
    // auto caloHits = _evStore->GetObject<Herd::CaloHits>("caloHitsMC");
    // auto mcTruth = _evStore->GetObject<Herd::MCTruth>("mcTruth");
    // retrieve MC truth info for studies:
    auto mcTruth = _evStore->GetObject<Herd::MCTruth>("mcTruth");
    const auto &primaryPart = mcTruth->primaries.at(0);
    auto initPos = primaryPart.InitialPosition();
    auto momentum = primaryPart.InitialMomentum();
    auto trajectory = primaryPart.InitialTrajectory();
    auto eTrue = momentum.Mag();
    auto MCInteraction = primaryPart.PrimaryDisappearance();
    auto disHeight = MCInteraction.intPosition.Mag();

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

    // std::cout << MCInteraction.processName << '\n';
    // std::cout << MCInteraction.Products().size() << '\n';
    // if (MCInteraction.Products().size()) {
    //   for (auto particle : MCInteraction.Products()) {
    //     std::cout << particle.PDGCode() << '\n';
    //   }
    // }
    // std::cout << MCInteraction.Products() << '\n';





    auto psdGeoParamsColl = _globStore->GetObject<Herd::PsdGeoParamsColl>("psdGeoParamsColl");
    auto caloGeoParams = _globStore->GetObject<Herd::CaloGeoParams>("caloGeoParams");
    auto psdHitsColl = _evStore->GetObject<Herd::PsdHitsColl>("psdHitsCollMC");
    auto fitGeoParamsColl = _globStore->GetObject<Herd::FitGeoParamsColl>("fitGeoParamsColl");

    // auto triggerPatternEven = _evStore->GetObject<TriggerPattern>("evenTrigger");
    // auto triggerPatternOdd = _evStore->GetObject<TriggerPattern>("oddTrigger");
    if (!_evStore->IsPresent("triggersFirstLayers"+_triggerAlgoName)) goto end;
    auto x=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000;


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



    auto triggersFirstLayers = _evStore->GetObject<Herd::SidesArray<std::vector<std::vector<TriggerPattern>>>>("triggersFirstLayers"+_triggerAlgoName);


    // int vetoedPatterns=0;
    // auto testPos = Herd::Vec3D(1,1,0);
    // auto geop = (*psdGeoParamsColl)[Herd::RefFrame::Side::TOP][0];
    //
    // for (size_t i = 0; i < geop.NTotElements(); i++) {
    //   /* code */
    //   auto pos = geop.Position(i);
    //   if (std::fabs(pos[Herd::RefFrame::Coo::X])<5 && std::fabs(pos[Herd::RefFrame::Coo::Y])<5) std::cout << i << " " << pos << '\n';
    //   // if (((*psdGeoParamsColl)[Herd::RefFrame::Side::TOP][0].Position(i)*testPos*testPos).Mag()<10) std::cout << "0" << '\n';
    //   // std::cout <<  << '\n';
    // }
    // geop = (*psdGeoParamsColl)[Herd::RefFrame::Side::TOP][1];
    // for (size_t i = 0; i < geop.NTotElements(); i++) {
    //   /* code */
    //   auto pos = geop.Position(i);
    //   if (std::fabs(pos[Herd::RefFrame::Coo::X])<5 && std::fabs(pos[Herd::RefFrame::Coo::Y])<5) std::cout << i << " " << pos << '\n';
    //   // std::cout << (*psdGeoParamsColl)[Herd::RefFrame::Side::TOP][1].Position(i) << '\n';
    // }
    // // std::cout << "ID: " << (*psdGeoParamsColl)[Herd::RefFrame::Side::TOP][0].VolumeID(testPos) << '\n';
    // // std::cout << "ID: " << (*psdGeoParamsColl)[Herd::RefFrame::Side::TOP][1].VolumeID(testPos) << '\n';
    auto passRoiVeto = Herd::SidesArray<int>();
    passRoiVeto.fill(0);
    // auto passFullPsdVeto = Herd::SidesArray<int>();
    // passFullPsdVeto.fill(0);
    std::vector<float> maxERoIVector;
    int sidesTriggered=0;
    float maxEPsd = 0.;
    // bool passFullPsdVeto;
    for (auto side : Herd::RefFrame::Sides) {
      // auto normalAxis = Herd::RefFrame::ToAxis(Herd::RefFrame::NormalDirection()Direction(side));
      // auto cosIncidenceAngle = momentum[Herd::RefFrame::ToCoo(normalAxis)]/eTrue;
      // auto side_selected = _evStore->GetObject<Herd::RefFrame::Side>("SelectedFace");
      // auto side=*side_selected;

      int totalCombinations = 0;
      float maxERoI = 0.;
      int totalTiles=0;
      if ((*triggersFirstLayers)[side].size()>0){
        sidesTriggered++;
        totalCombinations = (*triggersFirstLayers)[side][0].size()*(*triggersFirstLayers)[side][1].size();
        PsdVetoList psdTileListEvenUnion(*psdGeoParamsColl);
        PsdVetoList psdTileListOddUnion(*psdGeoParamsColl);
        for (auto triggerPatternEven : (*triggersFirstLayers)[side][0]) {
          // auto geomRoiEven = UlegGeomRoi(triggerPatternEven,*fitGeoParamsColl);
          // auto psdTileListEven = PsdVetoList(geomRoiEven,*psdGeoParamsColl);
          auto &psdTileListEven = _patternPsdMap->ConstructRetrieve(triggerPatternEven,*fitGeoParamsColl,*psdGeoParamsColl);
          psdTileListEvenUnion = psdTileListEvenUnion.Union(psdTileListEven);
        } // loop over even triggers
        for (auto triggerPatternOdd : (*triggersFirstLayers)[side][1]) {
          // auto geomRoiOdd = UlegGeomRoi(triggerPatternOdd,*fitGeoParamsColl);
          // auto psdTileListOdd = PsdVetoList(geomRoiOdd,*psdGeoParamsColl);
          auto &psdTileListOdd = _patternPsdMap->ConstructRetrieve(triggerPatternOdd,*fitGeoParamsColl,*psdGeoParamsColl);
          psdTileListOddUnion = psdTileListOddUnion.Union(psdTileListOdd);
        } // loop over odd triggers
        auto psdTileListJoint = psdTileListEvenUnion.Intersection(psdTileListOddUnion);
        // auto veto = psdTileListJoint.Veto(*psdHitsColl,_threshold);
        // passRoiVeto[side]=!veto;


        totalTiles = psdTileListJoint.NTotalTiles();
        // float maxEPsd =0.;
        // float maxESide =0.;
        // float totalERoI = 0.;
        // float totalEPsd =0.;
        // float totalESide =0.;
        ////
        // auto caloCubeListEven = CaloVetoList(geomRoiEven,*caloGeoParams);
        // auto caloCubeListOdd = CaloVetoList(geomRoiOdd,*caloGeoParams);
        // auto caloCubeListJoint = caloCubeListEven.Intersection(caloCubeListOdd);
        ////
        // auto triggerLayer = (*triggerPatternEven).Layer();
        for (auto psdSide : Herd::RefFrame::Sides) {
          auto psdHitsSide = (*psdHitsColl)[psdSide];
          for (size_t iPsdLayer=0; iPsdLayer<psdHitsSide.size();iPsdLayer++) {
            auto psdHits = psdHitsSide[iPsdLayer];
            for (auto psdHit : psdHits) {
              auto isInRoi = true;
              maxEPsd = std::max(maxEPsd,psdHit.EDep());
              if (totalCombinations<=_fullPsdCheckTh) {
                isInRoi = psdTileListJoint.IsInRoi(psdHit.VolumeID(),psdSide,iPsdLayer);
              }
              if (isInRoi) {
                maxERoI = std::max(maxERoI,psdHit.EDep());
              }
            } // loop over psd hits in the layer
          } // loop over psd layers
        } // loop over psd sides
        passRoiVeto[side]=maxERoI<_threshold;
        _tupleVetoSides->Fill(eTrue,cosPolar,azimuth,cosIncidenceAngles[side],iSideSelected,maxERoI,maxEPsd,totalCombinations,static_cast<int>(side),totalTiles,totalECalo);
        maxERoIVector.push_back(maxERoI);
      }
    } // loop over triggers sides
    // } // loop over both patterns
    bool passFullPsdVeto = maxEPsd<_threshold;
    auto sidesPassingRoi = std::accumulate(passRoiVeto.begin(),passRoiVeto.end(),0);
    int allSidesClean = sidesTriggered == sidesPassingRoi;
    auto minMaxERoI = *std::min_element(maxERoIVector.begin(),maxERoIVector.end());
    if (sidesPassingRoi) {
      _tupleVetoEvents->Fill(eTrue,cosPolar,azimuth,iSideSelected,minMaxERoI,sidesPassingRoi, totalECalo, disHeight, allSidesClean, passFullPsdVeto);
      _nAccept++;
      SetFilterResult(FilterResult::ACCEPT);
    } else {
      _nReject++;
      SetFilterResult(FilterResult::REJECT);
    }

    auto y=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000-x;
    // std::cout << GetName()+": " << y << '\n';
    // _tupleE->Fill(eTrue,cosIncidenceAngle,maxERoI,triggerPatternEven->Layer(), vetoedPatterns);
  }
  catch (Retrieval::Exception &exc) {
    COUT(ERROR) << exc.what() << ENDL;
    return false;
  }

  end:
  return true;
}

bool MultiRoiVeto::Finalize() {
  const std::string routineName("MultiRoiVeto::Finalize");

  COUT(INFO) << "Selected events: " << _nAccept << "/" << _nReject + _nAccept << " ("
  << (float)_nAccept / (float)(_nReject + _nAccept) * 100. << "\%)" << ENDL;

  try {
    // _globStore->AddObject(_tupleE->GetName(), _tupleE);
    _globStore->AddObject(_tupleVetoSides->GetName(), _tupleVetoSides);
    _globStore->AddObject(_tupleVetoEvents->GetName(), _tupleVetoEvents);
    // _globStore->AddObject(_tupleAllHits->GetName(), _tupleAllHits);
  } catch (Retrieval::Exception &) {
    COUT(ERROR) << "Global data store not found." << ENDL;
    return false;
  } catch (Insertion::Exception &) {
    COUT(ERROR) << "Can't add object " << GetName() << " to the data store." << ENDL;
  }
  return true;
}
