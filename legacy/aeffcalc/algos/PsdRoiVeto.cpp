/*
* PsdRoiVeto.cpp
*
*  Created on: 11 Feb 2019
*      Author: Nicola Mori
*/

// Reject the event that depose less than 50% of photon energy in the calorimeter

// Example headers
#include "PsdRoiVeto.h"
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

RegisterAlgorithm(PsdRoiVeto);

PsdRoiVeto::PsdRoiVeto(const std::string &name) : _threshold{1e-3}, _triggerAlgoName{GetName()}, Algorithm{name},
_sidesToStudy{} {

  DefineParameter("threshold", _threshold);
  // DefineParameter("sidesToStudy", _sidesToStudy);
  DefineParameter("triggerAlgoName", _triggerAlgoName);

}

bool PsdRoiVeto::Initialize() {
  const std::string routineName("PsdRoiVeto::Initialize");

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

  // if (_sidesToStudy.empty()) {
  //   _sides = std::vector<Herd::RefFrame::Side>(Herd::RefFrame::Sides.begin(),Herd::RefFrame::Sides.end());
  // } else {
  //   for (auto sideNumber : _sidesToStudy) {
  //     _sides.push_back(Herd::RefFrame::Side(sideNumber));
  //   }
  // }

  // _tupleE=std::make_shared<TNtuple>(GetName().c_str(),"temp", "eTrue:cosTheta:maxERoI:iLayer:nVetoes");
  _tupleVetoEvents=std::make_shared<TNtuple>(("PsdRoiVeto"+GetName()+"Events").c_str(),"temp", "eTrue:cosPolar:azimuth:iSideSelected:maxERoI:maxEPsd:maxESide:totalERoI:totalEPsd:totalESide");
  _tupleVetoSides=std::make_shared<TNtuple>(("PsdRoiVeto"+GetName()+"Sides").c_str(),"temp", "eTrue:cosPolar:azimuth:iSideSelected:cosIncidenceAngle:triggerSide:maxERoI:maxEPsd:maxESide:totalERoI:totalEPsd:totalESide:iLayer");
  // _tupleVeto->Fill(eTrue,cosPolar,azimuth,iSideSelected,cosIncidenceAngles[side],static_cast<int>(side),maxERoI,maxEPsd,maxESide,totalERoI,totalEPsd,totalESide,triggerLayer);

  // _tupleAllHits=std::make_shared<TNtuple>((GetName()+"AllHits").c_str(),"temp", "eTrue:cosTheta:eDep:iSide:iLayer:iTile");
  // _tupleAllHitsInTrigger=std::make_shared<TNtuple>((GetName()+"AllHitsInTrigger").c_str(),"temp", "eTrue:cosIncidenceAngle:triggerSide:eDep:iSide:iLayer:iTile");
  return true;
}

bool PsdRoiVeto::Process() {
  const std::string routineName("PsdRoiVeto::Process");
  try {
    // auto caloHits = _evStore->GetObject<Herd::CaloHits>("caloHitsMC");
    // auto mcTruth = _evStore->GetObject<Herd::MCTruth>("mcTruth");
    // auto side_selected = _evStore->GetObject<Herd::RefFrame::Side>("SelectedFace");
    // retrieve MC truth info for studies:
    auto mcTruth = _evStore->GetObject<Herd::MCTruth>("mcTruth");
    const auto &primaryPart = mcTruth->primaries.at(0);
    auto initPos = primaryPart.InitialPosition();
    auto momentum = primaryPart.InitialMomentum();
    auto trajectory = primaryPart.InitialTrajectory();
    auto eTrue = momentum.Mag();
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
    // std::cout << iSideSelected << '\n';

    // auto MCInteraction = primaryPart.PrimaryDisappearance();
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

    // for (auto psdSide : Herd::RefFrame::Sides) {
    //   auto psdHitsSide = (*psdHitsColl)[psdSide];
    //   for (size_t iPsdLayer=0; iPsdLayer<psdHitsSide.size();iPsdLayer++) {
    //     auto psdHits = psdHitsSide[iPsdLayer];
    //     for (auto psdHit : psdHits) {
    //       _tupleAllHits->Fill(eTrue,cosIncidenceAngles[psdSide],psdHit.EDep(),static_cast<int>(psdSide),iPsdLayer,psdHit.VolumeID());
    //     }
    //   }
    // }


    if (!_evStore->IsPresent("firstTriggers"+_triggerAlgoName)) goto end;
    auto firstTriggers = _evStore->GetObject<Herd::SidesArray<std::vector<TriggerPattern>>>("firstTriggers"+_triggerAlgoName);
    int vetoedPatterns=0;
    auto passVeto = Herd::SidesArray<int>();
    passVeto.fill(0);

    std::vector<float> maxERoIVector;
    std::vector<float> maxEPsdVector;
    std::vector<float> maxESideVector;
    std::vector<float> totalERoIVector;
    std::vector<float> totalEPsdVector;
    std::vector<float> totalESideVector;


    for (auto side : Herd::RefFrame::Sides) {
      if ((*firstTriggers)[side].size()>0){
        auto normalAxis = Herd::RefFrame::ToAxis(Herd::RefFrame::NormalDirection(side));
        auto cosIncidenceAngle = momentum[Herd::RefFrame::ToCoo(normalAxis)]/eTrue;
        // auto cosPolar = momentum[Herd::RefFrame::Coo::Z]/eTrue;
        // auto azimuth = momentum.AzimuthAngle();
        auto triggerPatternEven = (*firstTriggers)[side][0];
        auto triggerPatternOdd = (*firstTriggers)[side][1];



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

        // for (auto triggerPattern : {*triggerPatternEven,*triggerPatternOdd}) {

        float maxERoI = 0.;
        float maxEPsd =0.;
        float maxESide =0.;
        float totalERoI = 0.;
        float totalEPsd =0.;
        float totalESide =0.;
        ////
        // auto geomRoiEven = UlegGeomRoi(triggerPatternEven,*fitGeoParamsColl);
        // auto geomRoiOdd = UlegGeomRoi(triggerPatternOdd,*fitGeoParamsColl);
        // auto psdTileListEven = PsdVetoList(geomRoiEven,*psdGeoParamsColl);
        // auto psdTileListOdd = PsdVetoList(geomRoiOdd,*psdGeoParamsColl);
        auto &psdTileListEven = _patternPsdMap->ConstructRetrieve(triggerPatternEven,*fitGeoParamsColl,*psdGeoParamsColl);

        auto &psdTileListOdd = _patternPsdMap->ConstructRetrieve(triggerPatternOdd,*fitGeoParamsColl,*psdGeoParamsColl);

        auto psdTileListJoint = psdTileListEven.Intersection(psdTileListOdd);
        // auto caloCubeListEven = CaloVetoList(geomRoiEven,*caloGeoParams);
        // auto caloCubeListOdd = CaloVetoList(geomRoiOdd,*caloGeoParams);
        // auto caloCubeListJoint = caloCubeListEven.Intersection(caloCubeListOdd);
        ////
        auto veto = psdTileListJoint.Veto(*psdHitsColl,_threshold);
        passVeto[side]=!veto;
        // if (!veto) {
        //   std::cout << side << '\n';
        //   std::cout << psdTileListJoint << '\n';
        // }
        auto triggerLayer = triggerPatternEven.Layer();
        auto triggerSide = triggerPatternEven.Side();

        // if (veto) vetoedPatterns++;
        for (auto psdSide : Herd::RefFrame::Sides) {
          auto psdHitsSide = (*psdHitsColl)[psdSide];
          auto isInSide = triggerPatternEven.Side()==psdSide; // both are in same side
          //   // std::cout << "psd side: " << side << '\n';
          for (size_t iPsdLayer=0; iPsdLayer<psdHitsSide.size();iPsdLayer++) {


            auto psdHits = psdHitsSide[iPsdLayer];
            for (auto psdHit : psdHits) {
              // _tupleAllHitsInTrigger->Fill(eTrue,cosIncidenceAngles[side],static_cast<int>(side),psdHit.EDep(),static_cast<int>(psdSide),iPsdLayer,psdHit.VolumeID());

              if (isInSide) {
                maxESide = std::max(maxERoI,psdHit.EDep());
                totalESide +=psdHit.EDep();
              }

              auto isInRoi = psdTileListJoint.IsInRoi(psdHit.VolumeID(),psdSide,iPsdLayer);
              if (isInRoi) {
                maxERoI = std::max(maxERoI,psdHit.EDep());
                totalERoI +=psdHit.EDep();
              }
              maxEPsd = std::max(maxEPsd,psdHit.EDep());
              totalEPsd +=psdHit.EDep();
            } // loop over psd hits in the layer
          } // loop over psd layers
        } // loop over sides in the PSD
        _tupleVetoSides->Fill(eTrue,cosPolar,azimuth,iSideSelected,cosIncidenceAngles[side],static_cast<int>(side),maxERoI,maxEPsd,maxESide,totalERoI,totalEPsd,totalESide,triggerLayer);
        // } // loop over both patterns
        maxERoIVector.push_back(maxERoI);
        maxEPsdVector.push_back(maxEPsd);
        maxESideVector.push_back(maxESide);
        totalERoIVector.push_back(totalERoI);
        totalEPsdVector.push_back(totalEPsd);
        totalESideVector.push_back(totalESide);
      }
    }
    auto sidesPassing = std::accumulate(passVeto.begin(),passVeto.end(),0);

    auto minMaxERoI = *std::min_element(maxERoIVector.begin(),maxERoIVector.end());
    auto minMaxEPsd = *std::min_element(maxEPsdVector.begin(),maxEPsdVector.end());
    auto minMaxESide = *std::min_element(maxESideVector.begin(),maxESideVector.end());
    auto minTotalERoI = *std::min_element(totalERoIVector.begin(),totalERoIVector.end());
    auto minTotalEPsd = *std::min_element(totalEPsdVector.begin(),totalEPsdVector.end());
    auto minTotalESide = *std::min_element(totalESideVector.begin(),totalESideVector.end());
    if (sidesPassing) {
      _tupleVetoEvents->Fill(eTrue,cosPolar,azimuth,iSideSelected,minMaxERoI,minMaxEPsd,minMaxESide,minTotalERoI,minTotalEPsd,minTotalESide);
      _nAccept++;
      SetFilterResult(FilterResult::ACCEPT);
    } else {
      _nReject++;
      SetFilterResult(FilterResult::REJECT);
    }
    // _tupleE->Fill(eTrue,cosIncidenceAngle,maxERoI,triggerPatternEven->Layer(), vetoedPatterns);
  }
  catch (Retrieval::Exception &exc) {
    COUT(ERROR) << exc.what() << ENDL;
    return false;
  }
  end:
  return true;
}

bool PsdRoiVeto::Finalize() {
  const std::string routineName("PsdRoiVeto::Finalize");

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
