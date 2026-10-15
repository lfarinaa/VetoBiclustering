/*
 * TriggerArrays.h
 *
 *  Created on: 12 Oct 2019
 *      Author: Valerio Formato
 */

#ifndef TriggerArrays_H_
#define TriggerArrays_H_

#include "algorithm/Algorithm.h"
#include "herd/analysis/dataobjects/FitGeoParamsColl.h"
#include "herd/analysis/dataobjects/PsdGeoParamsColl.h"
#include "herd/analysis/dataobjects/ScdGeoParamsColl.h"
#include "herd/analysis/dataobjects/CaloGeoParams.h"
#include "UlegGeomPsdRoi.h"



using namespace EA;

class TNtuple;
class TTree;


// clang-format off
/*! @brief An algorithm that finds the Fit hits closest to the MC truth incoming direction
 * * <B>Needed event objects:</B>
 *   name        |     type    |  store  | optional | description
 * --------------|-------------|---------|----------|---------------------------------------------------------------------------------------------------
 * mcTruth       | MCTruth     | evStore |    no    | Hits in the Fit
 * FitHitsCollMC | FitHitsColl | evStore |    no    | Hits in the Fit
 *
 * * <B>Needed global objects:</B>
 *   name        |     type    |  store  | optional | description
 * --------------|-------------|---------|----------|---------------------------------------------------------------------------------------------------
 * FitGeoParamsColl | FitGeoParamsColl | evStore |    no    | Hits in the Fit
 *
 * * <B>Published event objects:</B>
 *
 *   name               |     type    |  store  | optional | description
 * ---------------------|-------------|---------|----------|---------------------------------------------------------------------------------------------------
 * FitPrimaryHitsCollMC | FitHitsColl | evStore |    no    | Primary hits in the Fit
 *
 * <B>Parameters:</B>
 *  name       |    type     | default value | description
 * ------------|-------------|---------------|-----------------------------------------------------------------------------------------
 *  distance   |   float     |  0.1 (cm)     | The maximum distance from the MCTruth impact point
 */
// clang-format on

class TriggerArrays : public Algorithm {
public:
  /*! @brief Constructor.
   *
   * @param name The name of the algorithm object.
   */
  TriggerArrays(const std::string &name);

  /*! @brief Initializes the cut.
   *
   * @return true if initialization is done without errors, false otherwise.
   */
  bool Initialize();

  /*! @brief Discard events below a given energy deposit threshold.
   *
   * @return true if no error occurs during processing, false otherwise.
   */
  bool Process();

  /*! @brief Prints the selection efficiency.
   *
   * @return true if no error occurs during finalization, false otherwise.
   */
  bool Finalize();

private:
  float _threshold;
  int _maxDelta;
  int _nHits;
  int _patternSize;
  int _layer_start;
  bool _force_full_spanning;
  bool _disable_filter;
  std::vector<std::string> _ignore_pattern_IDs;
  std::vector<int> _sidesToStudy;
  // std::vector<Herd::RefFrame::Side> _sides;
  // bool _ignore_plusminus;
  bool _require_only_one_side;
  bool _require_calo_in_roi;
  bool _removeUpwardsPatterns;
  bool _removeAlmostUpwardsPatterns;
  std::string _matHitsCollName;
  std::string _store_patterns_name;
  std::string _suppress_patterns_name;

  // Utility variables
  observer_ptr<GlobalDataStore> _globStore; // Pointer to the event data store
  observer_ptr<EventDataStore> _evStore; // Pointer to the event data store
  observer_ptr<Herd::FitGeoParamsColl> _FitGeoParamsColl; // Pointer to the FitGeoParams collection
  observer_ptr<Herd::PsdGeoParamsColl> _PsdGeoParamsColl; // Pointer to the PsdGeoParams collection
  observer_ptr<Herd::CaloGeoParams> caloGeoParams; // Pointer to the PsdGeoParams collection
  // observer_ptr<Herd::ScdGeoParamsColl> _ScdGeoParamsColl; // Pointer to the PsdGeoParams collection
  // std::shared_ptr<TTree> _treeTriggers; // Objects to be pushed on global store must be held by a shared_ptr
  // std::shared_ptr<TNtuple> _tupleTriggersEven; // Objects to be pushed on global store must be held by a shared_ptr
  // std::shared_ptr<TNtuple> _tupleFirstTriggersEven; // Objects to be pushed on global store must be held by a shared_ptr
  // std::shared_ptr<TNtuple> _tupleTriggersOdd; // Objects to be pushed on global store must be held by a shared_ptr
  // std::shared_ptr<TNtuple> _tupleTriggersOdd; // Objects to be pushed on global store must be held by a shared_ptr
  std::shared_ptr<TNtuple> _tupleTriggers; // Objects to be pushed on global store must be held by a shared_ptr
  std::shared_ptr<TNtuple> _tupleSides; // Objects to be pushed on global store must be held by a shared_ptr
  std::shared_ptr<TNtuple> _tupleTriggerEvents; // Objects to be pushed on global store must be held by a shared_ptr
  // std::shared_ptr<TNtuple> _tupleFirstTriggersOdd; // Objects to be pushed on global store must be held by a shared_ptr
  // std::shared_ptr<TNtuple> _tupleFirstTriggersOdd; // Objects to be pushed on global store must be held by a shared_ptr
  std::shared_ptr<TNtuple> _tupleFirstTriggers; // Objects to be pushed on global store must be held by a shared_ptr
  std::shared_ptr<TNtuple> _tupleFirstLayer; // Objects to be pushed on global store must be held by a shared_ptr
  std::shared_ptr<TNtuple> _tupleMats; // Objects to be pushed on global store must be held by a shared_ptr
  // observer_ptr<PatternToCaloVetoListMap> _patternCaloMap;
  observer_ptr<PatternVetoListMap<CaloVetoList,Herd::CaloGeoParams>> _patternCaloMap;
  int _nAccept, _nReject;

};

#endif /* TriggerArrays_H_ */
