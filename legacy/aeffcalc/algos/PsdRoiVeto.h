/*
 * PsdRoiVeto.h
 *
 *  Created on: 11 Feb 2019
 *      Author: Nicola Mori
 */

#ifndef PsdRoiVeto_H_
#define PsdRoiVeto_H_

#include "algorithm/Algorithm.h"

#include "herd/analysis/dataobjects/FitGeoParamsColl.h"
#include "herd/analysis/dataobjects/PsdGeoParamsColl.h"
#include "herd/analysis/dataobjects/ScdGeoParamsColl.h"
#include "UlegGeomPsdRoi.h"

using namespace EA;

class TNtuple;

/*! @brief A cut/filter class that rejects events below a given energy deposit in the calorimeter.
 *
 * * <B>Needed event objects:</B>
 *
 *   name       |     type    |  store  | optional | description
 * -------------|-------------|---------|----------|---------------------------------------------------------------------------------------------------
 * caloHits     |   CaloHits   | evStore |    no    | Hits in the CALO.
 *
 *
 * <B>Parameters:</B>
 *
 *  name       |    type     | default value | description
 * ------------|-------------|---------------|-----------------------------------------------------------------------------------------
 *  minEDep    |   float     |      5.       | The minimum value for energy deposit [GeV].
 */
class PsdRoiVeto : public Algorithm {
public:
  /*! @brief Constructor.
   *
   * @param name The name of the algorithm object.
   */
  PsdRoiVeto(const std::string &name);

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

  // Algorithm parameters
  // float _minEDep;


  // Utility variables
  observer_ptr<EventDataStore> _evStore; // Pointer to the event data store
  observer_ptr<GlobalDataStore> _globStore; // Pointer to the global data store
  // int _nAccept, _nReject;
  float _threshold;
  std::string _triggerAlgoName;
  std::vector<int> _sidesToStudy;
  std::vector<Herd::RefFrame::Side> _sides;
  observer_ptr<PatternVetoListMap<PsdVetoList,Herd::PsdGeoParamsColl>> _patternPsdMap;
  // std::shared_ptr<TNtuple> _tupleE;
  // std::shared_ptr<TNtuple> _tupleAllHits;
  // std::shared_ptr<TNtuple> _tupleAllHitsInTrigger;
  std::shared_ptr<TNtuple> _tupleVetoEvents;
  std::shared_ptr<TNtuple> _tupleVetoSides;
  // std::shared_ptr<TNtuple> _tupleInRoI;
  // std::shared_ptr<TNtuple> _tupleZ;
  int _nAccept, _nReject;

};

#endif /* PsdRoiVeto_H_ */
