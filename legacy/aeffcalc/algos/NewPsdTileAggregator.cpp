/*
 * PsdGeometricDigitizerAlgo.cpp
 *
 *  Created on: 10 Jan 2019
 *      Author: Nicola Mori
 */

// HerdSoftware headers
#include "algorithms/digitization/PsdGeometricDigitizerAlgo.h"
#include "herd/analysis/dataobjects/Hit.h"

// C/C++ standard headers
#include <algorithm>
#include <cassert>

namespace Herd {

RegisterAlgorithm(PsdGeometricDigitizerAlgo);

PsdGeometricDigitizerAlgo::PsdGeometricDigitizerAlgo(const std::string &name)
    : Algorithm(name), _aggrPsdGeoParamsColl(std::make_shared<PsdGeoParamsColl>()) {
  _topAggrFactor[RefFrame::Axis::X] = _topAggrFactor[RefFrame::Axis::Y] = _topAggrFactor[RefFrame::Axis::Z] = 1;
  _latAggrFactor[RefFrame::Axis::X] = _latAggrFactor[RefFrame::Axis::Y] = _latAggrFactor[RefFrame::Axis::Z] = 1;
  DefineParameter("topAggrFactorX", _topAggrFactor[RefFrame::Axis::X]);
  DefineParameter("topAggrFactorY", _topAggrFactor[RefFrame::Axis::Y]);
  DefineParameter("latAggrFactorXY", _latAggrFactor[RefFrame::Axis::X]);
  DefineParameter("latAggrFactorZ", _latAggrFactor[RefFrame::Axis::Z]);

  DeclareConsumedObject("psdHitsCollMC", ObjectCategory::EVENT, "evStore");
  DeclareConsumedObject("psdParticleHitsCollMC", ObjectCategory::EVENT, "evStore");
  DeclareConsumedObject("psdGeoParamsColl", ObjectCategory::GLOBAL, "globStore");

  DeclareProducedObject("psdHitsCollAggregated", ObjectCategory::EVENT, "psdHitsCollMC", "evStore");
  DeclareProducedObject("psdParticleHitsCollAggregated", ObjectCategory::EVENT, "psdParticleHitsCollMC", "evStore");
  DeclareProducedObject("psdGeoParamsCollAggregated", ObjectCategory::GLOBAL, "psdGeoParamsColl", "globStore");

  // Set geometry helper containers
  _tileMainAggrAxis[RefFrame::Side::TOP] = RefFrame::Axis::Y;
  _tileSecAggrAxis[RefFrame::Side::TOP] = RefFrame::Axis::X;
  _tileMainAggrAxis[RefFrame::Side::Xneg] = _tileMainAggrAxis[RefFrame::Side::Xpos] =
      _tileMainAggrAxis[RefFrame::Side::Yneg] = _tileMainAggrAxis[RefFrame::Side::Ypos] = RefFrame::Axis::Z;
  _tileSecAggrAxis[RefFrame::Side::Xneg] = _tileSecAggrAxis[RefFrame::Side::Xpos] = RefFrame::Axis::Y;
  _tileSecAggrAxis[RefFrame::Side::Yneg] = _tileSecAggrAxis[RefFrame::Side::Ypos] = RefFrame::Axis::X;
}

bool PsdGeometricDigitizerAlgo::Initialize() {
  const std::string routineName("PsdGeometricDigitizerAlgo::Initialize");

  // Set Y aggregation factor for Y (for X it has been set as the sideAggrFactorXY parameter)
  _latAggrFactor[RefFrame::Axis::Y] = _latAggrFactor[RefFrame::Axis::X];

  try {
    // Retrieve data stores once to avoid having to retrieve them for each event.
    _evStore = GetDataStoreManager()->GetEventDataStore("evStore");
    _globStore = GetDataStoreManager()->GetGlobalDataStore("globStore");

    // Retrieve PSD geometric parameters
    _origPsdGeoParamsColl = _globStore->GetObject<PsdGeoParamsColl>("psdGeoParamsColl");

    AxesArray<unsigned int> nAggrElements;
    AxesArray<float> aggrSizes;
    // Initialize the aggregated parameters equal to the original ones
    (*_aggrPsdGeoParamsColl) = (*_origPsdGeoParamsColl);
    // Compute aggregated geometry parameters
    for (auto side : RefFrame::Sides) { // Loop on all sides
      if ((*_origPsdGeoParamsColl)[side].size() != 0) {
        auto nLayers = (*_origPsdGeoParamsColl)[side].size();
        for (decltype(nLayers) iLayer = 0; iLayer < nLayers; ++iLayer) { // Loop on all layers
          PsdGeoParams &origLayerParams = (*_origPsdGeoParamsColl)[side][iLayer];
          PsdGeoParams &aggrLayerParams = (*_aggrPsdGeoParamsColl)[side][iLayer];
          auto segmDirs = origLayerParams.SegmentationDirections();
          std::pair<RefFrame::Axis, RefFrame::Axis> segmAxes;
          segmAxes.first = RefFrame::ToAxis(segmDirs.first);
          segmAxes.second = RefFrame::ToAxis(segmDirs.second);
          int nOrigPerAggr = 1; // Number of original elements in each aggregated element
          auto aggrFactor = (side == RefFrame::Side::TOP ? _topAggrFactor : _latAggrFactor);
          for (auto axis : RefFrame::Axes) { // Loop on all axes
            // Set the aggregation factor to 1 if the number of elements is 1 (e.g. for bars along the longitudinal
            // axis)
            if (origLayerParams.NElements()[axis] == 1) {
              aggrFactor[axis] = 1;
            }
            // Check that aggregation factors are exact divisors of the number of elements (if number of elements is not
            // already 1)
            if (origLayerParams.NElements()[axis] % aggrFactor[axis] != 0) {
              COUT(ERROR) << "Layer " << iLayer << ": aggregation factor for side " << side << " along " << axis << " ("
                          << aggrFactor[axis] << ") is not a divisor of the number of " << side << " elements along "
                          << axis << " (" << origLayerParams.NElements()[axis] << ")" << ENDL;
              return false;
            }
            // Check that no gap is present along the aggregation direction
            if (aggrFactor[axis] != 1) {
              int nearbyElem = (axis == segmAxes.first ? 1 : origLayerParams.NElements()[segmAxes.first]);
              if (origLayerParams.Sizes()[axis] <
                  std::abs(origLayerParams.Position(0)[RefFrame::ToCoo(axis)] -
                           origLayerParams.Position(nearbyElem)[RefFrame::ToCoo(axis)])) {
                COUT(ERROR) << "Gap along the " << axis << " axis for " << side
                            << " PSD; cannot aggregate PSD elements along " << axis << ENDL;
                return false;
              }
            }
            nAggrElements[axis] = origLayerParams.NElements()[axis] / aggrFactor[axis];
            aggrSizes[axis] = origLayerParams.Sizes()[axis] * aggrFactor[axis];
            nOrigPerAggr *= aggrFactor[axis];
          }
          // Set the aggregated geometry parameters
          aggrLayerParams.SetNElements(nAggrElements);
          aggrLayerParams.SetSizes(aggrSizes);
          // Set the position of each aggregated element as the mean position of the component original elements
          // Compute the sum of coordinates for each aggregated element
          std::vector<Point> sumCoo(aggrLayerParams.NTotElements(), Point(0., 0., 0.));
          for (size_t origVolID = 0; origVolID < origLayerParams.NTotElements(); ++origVolID) {
            // Retrieve the indexes of the current PSD element
            CooArray<unsigned char> aggrIndexes = _AggregatedIndexes(origLayerParams.CooIndexes(origVolID), aggrFactor);
            // Compute the aggregated volume ID
            //    volID   =   index on secondary axis   *  elements along primary axis  +   index on primary axis
            unsigned int aggrVolID = aggrIndexes[RefFrame::ToCoo(segmAxes.second)] * nAggrElements[segmAxes.first] +
                                     aggrIndexes[RefFrame::ToCoo(segmAxes.first)];
            assert(aggrVolID < aggrLayerParams.NTotElements());
            // Accumulate the coordinates of the original element on the aggregated element
            auto &origPoint = origLayerParams.Position(origVolID);
            for (auto coo : RefFrame::Coos) {
              sumCoo[aggrVolID][coo] += origPoint[coo];
            }
          }
          // Compute the position as the mean of the coordinates
          float nAggrElemInverse = 1. / (nOrigPerAggr);
          for (size_t aggrVolID = 0; aggrVolID < aggrLayerParams.NTotElements(); ++aggrVolID) {
            sumCoo[aggrVolID] *= nAggrElemInverse;
            aggrLayerParams.SetPosition(aggrVolID, sumCoo[aggrVolID]);
          }
        }
      }
    }

    // Publish new aggregated objects
    _globStore->AddObject("psdGeoParamsCollAggregated", _aggrPsdGeoParamsColl);
    _globStore->SetAlias("psdGeoParamsCollAggregated", "psdGeoParamsColl");
  } catch (Retrieval::Exception &exc) {
    COUT(ERROR) << exc.what() << ENDL;
    return false;
  } catch (Insertion::Exception &) {
    COUT(ERROR) << "Can't add the PSD aggregated geometry parameters object to the globaldata store." << ENDL;
    return false;
  }
  return true;
}

bool PsdGeometricDigitizerAlgo::Process() {
  static const std::string routineName("PsdGeometricDigitizerAlgo::Process");

  // Retrieve hits and particle hits
  observer_ptr<PsdHitsColl> psdHitsColl;
  try {
    psdHitsColl = _evStore->GetObject<PsdHitsColl>("psdHitsCollMC");
  } catch (Retrieval::Exception &) {
    COUT(ERROR) << "PSD hits collection \"psdHitsCollMC\" not found." << ENDL;
    return false;
  }
  observer_ptr<PsdParticleHitsColl> psdParticleHitsColl;
  try {
    psdParticleHitsColl = _evStore->GetObject<PsdParticleHitsColl>("psdParticleHitsCollMC");
  } catch (Retrieval::NotAvailable &) {
  } catch (Retrieval::NotFound &) {
  }

  if (!_origPsdGeoParamsColl) {
    COUT(ERROR) << "Hits are available but geometry parameters object is missing" << ENDL;
    return false;
  }

  // Digitize
  for (auto side : RefFrame::Sides) {
    auto &aggrFactor = (side == RefFrame::Side::TOP ? _topAggrFactor : _latAggrFactor);
    _Aggregate((*_origPsdGeoParamsColl)[side], (*psdHitsColl)[side],
               (psdParticleHitsColl ? &((*psdParticleHitsColl)[side]) : nullptr), aggrFactor,
               (*_aggrPsdGeoParamsColl)[side], _aggrPsdHitsColl[side], _aggrPsdParticleHitsColl[side]);
  }

  try {
    _evStore->AddObject("psdHitsCollAggregated", _aggrPsdHitsColl);
    _evStore->SetAlias("psdHitsCollAggregated", "psdHitsCollMC");
    if (psdParticleHitsColl) {
      _evStore->AddObject("psdParticleHitsCollAggregated", _aggrPsdParticleHitsColl);
      _evStore->SetAlias("psdParticleHitsCollAggregated", "psdParticleHitsCollMC");
    }
    _globStore->SetAlias("psdGeoParamsCollAggregated", "psdGeoParamsColl");
  } catch (Insertion::Exception &exc) {
    COUT(ERROR) << "While adding data objects to store: " << exc.what() << ENDL;
    return false;
  }
  return true;
}

std::vector<std::string> PsdGeometricDigitizerAlgo::FreeObjects(const std::vector<std::string> &objs,
                                                                Memory::Status memStatus) {
  std::vector<std::string> freed;
  if (memStatus == Memory::Status::EXHAUSTED) {
    for (auto &obj : objs) {
      if (obj == "psdHitsCollAggregated") {
        PsdHitsColl().swap(_aggrPsdHitsColl);
        freed.push_back(obj);
      } else if (obj == "psdParticleHitsCollAggregated") {
        PsdParticleHitsColl().swap(_aggrPsdParticleHitsColl);
        freed.push_back(obj);
      }
    }
  }
  return freed;
}

void PsdGeometricDigitizerAlgo::_Aggregate(const std::vector<PsdGeoParams> &originalGeoParams,
                                           const PsdHits &originalHits,
                                           observer_ptr<const PsdParticleHits> originalParticleHits,
                                           const AxesArray<int> &aggrFactors,
                                           const std::vector<PsdGeoParams> &aggregatedGeoParams,
                                           PsdHits &aggregatedHits, PsdParticleHits &aggregatedParticleHits) {
  if (originalHits.size() == 0) {
    return;
  }

  aggregatedHits.resize(originalHits.size());
  bool particleHits = false;
  if (originalParticleHits) {
    aggregatedParticleHits.resize(originalParticleHits->size());
    particleHits = true;
  }
  for (size_t iLayer = 0; iLayer < aggregatedHits.size(); ++iLayer) {

    // Retrieve objects for this layer
    auto &layerOrigHits = originalHits[iLayer];
    auto &layerAggrHits = aggregatedHits[iLayer];
    observer_ptr<const std::vector<std::vector<ParticleHit>>> layerOrigParticleHits{nullptr};
    observer_ptr<std::vector<std::vector<ParticleHit>>> layerAggrParticleHits{nullptr};
    if (particleHits) {
      layerOrigParticleHits = &((*originalParticleHits)[iLayer]);
      layerAggrParticleHits = &(aggregatedParticleHits[iLayer]);
    }
    auto &layerOrigGeoParams = originalGeoParams[iLayer];

    // Clean aggregated hits
    layerAggrHits.resize(0);
    if (particleHits) {
      layerAggrParticleHits->resize(0);
    }

    // Retrieve segmentation information
    auto segmDirs = aggregatedGeoParams[iLayer].SegmentationDirections();
    std::pair<RefFrame::Axis, RefFrame::Axis> segmAxes;
    segmAxes.first = RefFrame::ToAxis(segmDirs.first);
    segmAxes.second = RefFrame::ToAxis(segmDirs.second);

    auto actualAggrFactors = aggrFactors;
    // Fix aggregation factors for bars (i.e. set aggregation factor to 1 for the secondary segmentation direction)
    if (layerOrigGeoParams.NElements()[segmAxes.second] == 1) {
      actualAggrFactors[segmAxes.second] = 1;
    }

    for (size_t iHit = 0; iHit < layerOrigHits.size(); ++iHit) {
      auto &origHit = layerOrigHits[iHit];
      auto aggrIndexes = _AggregatedIndexes(layerOrigGeoParams.CooIndexes(origHit.VolumeID()), actualAggrFactors);
      unsigned int aggrVolID =
          aggrIndexes[RefFrame::ToCoo(segmAxes.second)] * aggregatedGeoParams[iLayer].NElements()[segmAxes.first] +
          aggrIndexes[RefFrame::ToCoo(segmAxes.first)];
      assert(aggrVolID < aggregatedGeoParams[iLayer].NTotElements());
      auto hitIter = std::find_if(layerAggrHits.begin(), layerAggrHits.end(),
                                  [aggrVolID](const Hit &hit) { return (hit.VolumeID() == aggrVolID); });
      if (hitIter == layerAggrHits.end()) {
        layerAggrHits.push_back(Hit(0., aggrVolID));
        if (particleHits) {
          layerAggrParticleHits->push_back(std::vector<ParticleHit>());
        }
        hitIter = layerAggrHits.end() - 1;
      }
      hitIter->SetEDep(hitIter->EDep() + origHit.EDep());
      if (particleHits) {
        auto &aggrHit = (*layerAggrParticleHits)[(size_t)(hitIter - layerAggrHits.begin())];
        aggrHit.insert(aggrHit.end(), (*layerOrigParticleHits)[iHit].begin(), (*layerOrigParticleHits)[iHit].end());
      }
    }
  }
}

CooArray<unsigned char> PsdGeometricDigitizerAlgo::_AggregatedIndexes(const CooArray<unsigned char> &origIndexes,
                                                                      const AxesArray<int> &aggrFactors) {
  CooArray<unsigned char> aggrIndexes;
  for (auto &axis : RefFrame::Axes) {
    aggrIndexes[RefFrame::ToCoo(axis)] = origIndexes[RefFrame::ToCoo(axis)] / aggrFactors[axis];
  }
  return aggrIndexes;
}

} // namespace Herd
