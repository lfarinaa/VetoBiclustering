/*
 * TriggerPattern.h
 *
 *  Created on: 12 Oct 2019
 *      Author: Valerio Formato
 */

#ifndef TriggerPattern_H_
#define TriggerPattern_H_


// C/C++ standard headers
#include <iterator>
#include <numeric>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>


// #include "algorithm/Algorithm.h"
#include "herd/analysis/dataobjects/FitGeoParamsColl.h"
#include "herd/analysis/dataobjects/PsdGeoParamsColl.h"
#include "herd/analysis/dataobjects/PsdHitsColl.h"
// #include "herd/analysis/dataobjects/ScdGeoParamsColl.h"
#include "herd/analysis/dataobjects/Line.h"
#include "herd/analysis/dataobjects/FitHitsColl.h"

//
//
// using namespace EA;
//
// class TNtuple;
// class TH1F;


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


// void calculate_patterns_recursion(std::vector<std::vector<int>> &vPatts, std::vector<int> vDeltas, std::vector<int> pattern, int nDeltas) {
//   // helper function to handle the recursion
//   for (auto delta : vDeltas){
//     pattern.push_back(delta);
//     nDeltas?calculate_patterns_recursion(vPatts,vDeltas,pattern,nDeltas-1):vPatts.push_back(pattern);
//     pattern.pop_back();
//   }
// }
//
// std::vector<std::vector<int>> calculate_patterns(int _pattSize, int _maxDelta) { // calculate all patterns of a given length and max lateral distance
//   const auto nDeltas=2*_maxDelta+1;
//   // there are (nDeltas)^_pattSize total patterns:
//   std::vector<int> vDeltas(nDeltas);
//   std::iota(vDeltas.begin(),vDeltas.end(),-_maxDelta);
//   std::vector<std::vector<int>> vPatts;
//   std::vector<int> pattern;
//   calculate_patterns_recursion(vPatts,vDeltas,pattern,_pattSize-1);
//   // std::cout << nPatts << " " << vPatts.size() << '\n';
//   return vPatts;
// }

// std::vector<std::vector<int>> find_patterns(std::vector<std::vector<int>> hitsarray, int pattSize, int maxDelta, int nWildcards) { // calculate all patterns of a given length and max lateral distance
//   const auto nDeltas=2*maxDelta+1;
//   // there are (nDeltas)^_pattSize total patterns:
//   std::vector<int> vDeltas(nDeltas);
//   std::string
//   std::iota(vDeltas.begin(),vDeltas.end(),-maxDelta);
//   std::vector<std::vector<int>> vPatts;
//   std::vector<int> pattern;
//   find_patterns_recursion(vPatts,vDeltas,pattern,pattSize-1);
//   // std::cout << nPatts << " " << vPatts.size() << '\n';
//   return vPatts;
// }

// std::vector<std::pair<int,int>> find_pattern_bruteforce(std::vector<std::vector<int>> hitsArray, std::vector<int> pattern) {
//   int nLayers = hitsArray.size();
//   int nMats = hitsArray[0].size();
//   std::vector<std::pair<int,int>> depthLocationVec;
//   int totalDelta=0;
//   for (size_t iPattern = 0; iPattern < pattern.size(); iPattern++) { // loop over patterns
//     totalDelta+=pattern[iPattern];
//     for (size_t iRow = 0; iRow < nLayers-pattern.size(); iRow++) { // loop over layers of mat hits array
//       if (std::max(0,-totalDelta) < std::min(nMats,nMats-totalDelta)) { // loop over hits
//         for (size_t iCol = std::max(0,-totalDelta); iCol < std::min(nMats,nMats-totalDelta); iCol++) {
//           hitsArray[iRow][iCol]+=hitsArray[iRow+iPattern+1][iCol+totalDelta];
//         }
//       } // loop over hits
//     } // loop over layers of mat hits array
//   } // loop over patterns
//   for (size_t iRow = 0; iRow < nLayers-pattern.size(); iRow++) { // loop over layers of mat hits array
//     for (size_t iCol = 0; iCol < nMats; iCol++) { // loop over hits
//       if (hitsArray[iRow][iCol]==pattern.size()+1) { // if enough matching hits, store matching pattern
//         std::pair<int,int> depthLocation={iRow,iCol};
//         depthLocationVec.push_back(depthLocation);
//       }
//     } // loop over hits
//   } // loop over layers of mat hits array
//   return depthLocationVec;
// }

//
// find_patterns (hitsArray, position, maxdelta, nhits, wildcards, patternsStorage)
// int nLayers = hitsArray.size();
// int nMats = hitsArray[0].size();
// for loop over layers
// for loop over mats
// std::string partialpattern = "";
//
// //
// find_pattern_recursion(hitsarray, position,nhits,nwildcards,partialpattern, patternsstorage) {
//   for
// }
//
//
//
//


// int patternID(std::vector<int> pattern) { // transform pattern vector to int code
//   int id=0;
//   for (size_t i = 0; i < pattern.size(); i++) {
//     id+=(5+pattern[i])*std::pow(10,i); // fine for low max delta and 3 hits
//   }
//   return id;
// }

// int subtriggerID (int depth, int location, std::vector<int> pattern) {
//   int patternID=patternID(pattern);
//   std::string idString = std::to_string(depth*100+location)+std::to_string(patternID);
//   int
// }

// void findPatterns(Int_t nHits,Int_t nRows,Int_t dM,Int_t nLayers,Int_t nMats,Bool_t* FIThits,Int_t iLayer,Int_t mat0,Int_t iPattern,vector<Int_t>& patterns);
// const Int_t kHit   = 5; // code for hit found with DeltaM=0


//
// Find and label all patterns of type <nHits>-hits-in-<nRows>-rows present in the <FIThits> array of L0 triggers in a FIT module (one tracking direction only)
// The FIT module consists of <nLayers> layers and <nMats> mats per layer
// We look for <nHits> hits in <nRows> consecutive layers, allowing a displacement of <dM> mats between consecutive layers
// The found patterns are labelled as <LMMabcde...> with
// L:  layer index of the first hit (or hole) of the pattern (one digit)
// MM: mat   index of the first hit (or hole) of the pattern (two digits)
// a, b, c, ...: code indicating mat difference between one layer and the previous one, with:
//   5 : no displacement (used by default for most upstream hit of a given pattern
//   4, 6: one-mat    displacement to left and right, respectively
//   3, 7: two-mats   displacement to left and right, respectively
//   2, 8: three-mats displacement to left and right, respectively
//   1, 9: four-mats  displacement to left and right, respectively
//
// void FITTrigger()
// {
//   // Configure the patterns searched for
//   const Int_t nHits = 3;         // how many hits in a row am I looking for
//   const Int_t nRows = nHits;     // in how many rows (>= nHits)
//   const Int_t  dM   = 1;         // maximum displacement in mat between consecutive layers to be considered part of the same pattern
//
//
//   // Initialize, define and configure environment
//   const Int_t nLayers = 7;      // number of layers of the FIT module
//   const Int_t nMats   = 10;     // number of mats in the layers (<=4)
//   Bool_t FIThits[nLayers*nMats] =    // level-0 hits pattern in the FIT module
//   {0,0,0,0,0,0,0,0,0,0,
//     1,0,0,0,0,0,0,0,0,0,
//     0,0,1,0,1,0,0,0,1,0,
//     0,1,0,1,0,1,0,0,0,1,
//     0,0,0,0,0,1,0,0,0,0,
//     1,0,1,0,0,0,0,0,0,0,
//     0,1,0,0,0,0,0,0,0,1};
//
//     // check for inconsistent demands
//     if(nRows<1 || nHits<1 || dM<0 || dM>4)
//     {
//       cout << "Error: nRows and nHits must be larger than 0 and dM cannot be negative or larger than 4";
//       return;
//     }
//     if(nRows<nHits || nRows>nLayers)
//     {
//       cout << "Error: nRows (" << nRows<< ") must be between nHits (" << nHits << ") and nLayers (" << nLayers << ")" << endl;
//       return;
//     }
//
//     // look for all the patterns in FIThits starting at all mats of all layers (except for the last <nRows> ones)
//     vector<TriggerPattern> allPatterns;
//     for(Int_t iLayer=0;iLayer<nLayers-nRows+1;iLayer++)
//     for(Int_t iMat=0;iMat<nMats;iMat++)
//     {
//       TriggerPattern iPattern=TriggerPattern(iLayer,iMat,0); // initial value of pattern
//       vector<Int_t> fPattern;   // vector to keep the all patterns that initiate at mat <iMat> of layer <iLayer>
//
//       if(FIThits[iLayer*nMats+iMat])
//       {
//         iPattern+= Int_t(kHit*TMath::Power(10,nRows-1));
//         findPatterns(nHits-1,nRows-1,dM,nLayers,nMats,FIThits,iLayer+1,iMat,iPattern,fPattern);
//       }
//       else
//       {
//         if(nHits==nRows) continue; // this happens when no more holes are allowed
//         findPatterns(nHits,nRows-1,dM,nLayers,nMats,FIThits,iLayer+1,iMat,iPattern,fPattern);
//       }
//
//       // keep all patterns in a general vector
//       for(Int_t ipat=0;ipat<fPattern.size();ipat++)
//       allPatterns.push_back(fPattern[ipat]);
//     }
//
//
//     // Report on found patterns
//     Int_t npatterns = allPatterns.size();
//     if(npatterns>0)
//     {
//       cout << npatterns << " patterns ("<< nHits << " hits in " << nRows << " rows) found: ";
//       for(Int_t ipat=0;ipat<allPatterns.size()-1;ipat++)
//       cout << allPatterns[ipat] << ", ";
//       if(allPatterns.size())
//       cout << allPatterns[allPatterns.size()-1] << endl;
//     }
//     else
//     cout << "No patterns ("<< nHits << " hits in " << nRows << " rows) found" << endl;
//   }

// std::vector<std::pair<int,int>> find_pattern_bruteforce(std::vector<std::vector<int>> hitsArray, std::vector<int> pattern, int nLayers, int nMats) {
//   std::vector<std::pair<int,int>> depthLocationVec;
//   int totalDelta=0;
//   for (size_t iPattern = 0; iPattern < pattern.size(); iPattern++) { // loop over patterns
//     totalDelta+=pattern[iPattern];
//     for (size_t iRow = 0; iRow < nLayers-pattern.size(); iRow++) { // loop over layers of mat hits array
//       if (std::max(0,-totalDelta) < std::min(nMats,nMats-totalDelta)) { // loop over hits
//         for (size_t iCol = std::max(0,-totalDelta); iCol < std::min(nMats,nMats-totalDelta); iCol++) {
//           hitsArray[iRow][iCol]+=hitsArray[iRow+iPattern+1][iCol+totalDelta];
//         }
//       } // loop over hits
//     } // loop over layers of mat hits array
//   } // loop over patterns
//   for (size_t iRow = 0; iRow < nLayers-pattern.size(); iRow++) { // loop over layers of mat hits array
//     for (size_t iCol = 0; iCol < nMats; iCol++) { // loop over hits
//       if (hitsArray[iRow][iCol]==pattern.size()+1) { // if enough matching hits, store matching pattern
//         std::pair<int,int> depthLocation={iRow,iCol};
//         depthLocationVec.push_back(depthLocation);
//       }
//     } // loop over hits
//   } // loop over layers of mat hits array
//   return depthLocationVec;
// }

using FitHitsArray = std::vector<std::vector<int>>;

class TriggerPattern{
public:
  // single hit init:
  TriggerPattern() :  _side{}, _layer{}, _mats{}, _holes{} {}
  TriggerPattern(Herd::RefFrame::Side side, int layer, int mat, int hole) :  _side{side}, _layer{layer}, _mats{mat}, _holes{hole} {}
  // full init:
  TriggerPattern(Herd::RefFrame::Side side, int layer, std::vector<int> mats, std::vector<int> holes) :  _side{side}, _layer{layer}, _mats{mats}, _holes{holes} {}
  // copy init:
  // TriggerPattern(const TriggerPattern &p1) : _side{p1.Side()}, _layer{p1.Layer()}, _mats{p1.Mats()}, _holes{p1.Holes()} {}
  // destroy:
  // virtual ~TriggerPattern() = default;

  // getters for priv vars:
  Herd::RefFrame::Side Side() const {return _side;}
  int Layer() const {return _layer;}
  std::vector<int> Mats() const {return _mats;};
  std::vector<int> Holes() const {return _holes;};

  // query tools:
  int Length() const {return _mats.size();} // total mats in pattern
  int NHits() const { // total hits (= mats-holes)
    int sum=0;
    for (size_t i = 0; i < _mats.size(); i++) {
      if (!_holes[i]) sum+=1;
    }
    return sum;
  }
  int NHoles() const {return (_mats.size()-NHits());} // total holes
  int Width() const {
    return *std::max_element(_mats.begin(),_mats.end())-*std::min_element(_mats.begin(),_mats.end())+1;
  }

  // add a mat to the pattern
  void append(int mat, int isHole) {
    _mats.push_back(mat);
    _holes.push_back(isHole);
  }

  // cut to size
  void cut(int size) {
    for (size_t i=_mats.size();i>size;i--) {
      _mats.pop_back();
      _holes.pop_back();
    }
  }

  bool IsEquivalent(TriggerPattern otherP, int firstHits) const {
    auto p1=TriggerPattern(*this); // copy so they can be cut
    auto p2=TriggerPattern(otherP); // copy so they can be cut
    p1.cut(firstHits);
    p2.cut(firstHits);
    if ((p1.Side() == p2.Side()) && (p1.Layer() == p2.Layer()) && (p1.Length() == p2.Length()) && (p1.Holes() == p2.Holes())) {
      for (size_t i = 0; i < p1.Length(); i++) {
        if (!p1.Holes()[i]) {
          if (p1.Mats()[i] != p2.Mats()[i]) {
            return false;
          }
        }
      }
      return true;
    } else {
      return false;
    }
  }

  bool PointsUpwards(Herd::FitGeoParamsColl fitGeoParamsColl) const {
    auto fitGeoParamsLayer = fitGeoParamsColl[_side][_layer];
    auto firstMatHeight = fitGeoParamsLayer.MatPosition(_mats.front())[Herd::RefFrame::Coo::Z];
    auto lastMatHeight = fitGeoParamsLayer.MatPosition(_mats.back())[Herd::RefFrame::Coo::Z];
    auto isUp = lastMatHeight > firstMatHeight;
    // std::cout << FullID() << " " << isUp << '\n';
    return isUp;
    // std::cout << _side << " " << _layer << " " << _mats.front() << " " << _mats.back() << '\n';
    // return _side!=Herd::RefFrame::Side::TOP&&!(_layer%2)&&_mats.front()>_mats.back();
  }
  bool PointsAlmostUpwards(Herd::FitGeoParamsColl fitGeoParamsColl) const {
    auto fitGeoParamsLayer = fitGeoParamsColl[_side][_layer];
    auto firstMatHeight = fitGeoParamsLayer.MatPosition(_mats.front())[Herd::RefFrame::Coo::Z];
    auto lastMatHeight = fitGeoParamsLayer.MatPosition(_mats.back())[Herd::RefFrame::Coo::Z];
    auto isQuasiUp = lastMatHeight = firstMatHeight;
    return isQuasiUp;
    // return _side!=Herd::RefFrame::Side::TOP&&!(_layer%2)&&_mats.front()==_mats.back();
  }



  void FindPatterns(int nHitsRemaining, int patternSize, int dM, FitHitsArray FITHits, std::vector<TriggerPattern>& patterns, bool force_full_spanning) {
    // find patterns of nHits hits in patternSize layers with mat distances up to dM;
    // starting in layer iLayer and mat mat0
    // then store them in vector patterns
    int mat0 = Mats().back();
    int iLayer = Layer()/2+Length();
    // std::cout << "Looking for patterns in layer " << iLayer << '\n';
    // std::cout << "partial: " << (*this) << '\n';
    // std::cout << "start on mat " << mat0 << '\n';
    int minMat=std::max(0,mat0-dM);
    int maxMat=std::min((int)(FITHits[iLayer].size()-1),mat0+dM);
    // std::cout << "look at mats in (" << minMat << ".." << maxMat << ")" << '\n';
    // if (minMat>maxMat) std::cout << "WARNING: bad limits in loop" << '\n';
    for (int iMat = minMat; iMat <= maxMat; iMat++) {
      TriggerPattern newPattern = (*this);
      int hasHit = FITHits[iLayer][iMat];
      newPattern.append(iMat,!hasHit);
      int lengthRemaining=patternSize-newPattern.Length();
      // std::cout << "mat " << iMat << " hits: " << hasHit << '\n';
      if(hasHit || lengthRemaining>nHitsRemaining) {
        bool check_span = force_full_spanning?(iLayer==(int)(FITHits.size()-1)):true;
        // bool check_span = force_full_spanning?(iLayer==(int)(FITHits.size()-1)||iMat==0||iMat==(int)(FITHits[iLayer].size()-1)):true;
        if(lengthRemaining<=0&&check_span) {
          // newPattern.cut(patternSize); //uncomment when done with >3 hit reporting
          patterns.push_back(newPattern);
          // std::cout << "found full pattern:" << '\n';
          // std::cout << newPattern << '\n';
        } else {
          // std::cout << "Trying to extend partial pattern:" << '\n';
          // std::cout << newPattern << '\n';
          // std::cout << lengthRemaining << check_span << '\n';
          newPattern.FindPatterns(nHitsRemaining-hasHit,patternSize,dM,FITHits,patterns,force_full_spanning);
        }
      }
    }
  }

  // various IDs of the pattern
  // The found patterns are labelled as <LMMabcde...> with
  // L:  layer index of the first hit (or hole) of the pattern (one digit)
  // MM: mat   index of the first hit (or hole) of the pattern (two digits)
  // a, b, c, ...: code indicating mat difference between one layer and the previous one, with:
  //   5 : no displacement (used by default for most upstream hit of a given pattern
  //   4, 6: one-mat    displacement to left and right, respectively
  //   3, 7: two-mats   displacement to left and right, respectively
  //   2, 8: three-mats displacement to left and right, respectively
  //   1, 9: four-mats  displacement to left and right, respectively
  //   0: no hits in this mat, but allowed due to nHits<nMats
  std::string PatternID() const { // TODO: find unique mapping (maps,holes) -> ID; have to pay some price for it
    // todo: add holes to pattern
    std::string pattID = _holes.front()?"0":"5"; // to do: should we allow patterns to start on a hole?
    for (size_t i = 1; i < _mats.size(); i++) {
      if (_holes[i]) {
        pattID += "0";
      } else {
        pattID += std::to_string(5+_mats[i]-_mats[i-1]);
      }
    }
    return pattID;
  }
  std::string SideLayerMatID() const { // LMM, where L is the first layer and MM the first mat
    auto sideID = Herd::RefFrame::SideName[static_cast<int>(_side)];
    auto layerID = std::to_string(_layer);
    auto matID = std::to_string(_mats.front());
    if (matID.size()<2) matID = "0"+matID; // zero-pad if single-digit
    return sideID+layerID+matID;
  }
  std::string FullID() const {
    return SideLayerMatID()+PatternID();
  }
  int OldID() const {
    // no support for holes
    int pattID = 0;
    for (size_t i = 1; i < _mats.size(); i++) {
      pattID *= 10;
      pattID += (_mats[i]-_mats[i-1]);
    }
    pattID=std::fabs(pattID);
    return pattID;
  }

private:
  Herd::RefFrame::Side _side;
  int _layer;
  std::vector<int> _mats;
  std::vector<int> _holes;
};

// overloaded operators:
inline bool operator==(const TriggerPattern &p1, const TriggerPattern &p2) {
  // patterns are "equal" if the contain the same active mats, regardless of inactive
  // this way we remove duplicates with sort-erase
  if ((p1.Side() == p2.Side()) && (p1.Layer() == p2.Layer()) && (p1.Length() == p2.Length()) && (p1.Holes() == p2.Holes())) {
    for (size_t i = 0; i < p1.Length(); i++) {
      if (!p1.Holes()[i]) {
        if (p1.Mats()[i] != p2.Mats()[i]) {
          return false;
        }
      }
    }
    return true;
  } else {
    return false;
  }
}

inline bool operator<(const TriggerPattern &p1, const TriggerPattern &p2) {
  // no support for p1 and p2 being different length
  if (p1.Side() != p2.Side()) {
    return p1.Side() < p2.Side();
  } else if (p1.Layer()!=p2.Layer()) { // compare origin layers
    return p1.Layer() < p2.Layer();
  // } else if (p1.Length() != p2.Length())
  } else if (p1.Holes() != p2.Holes()) { // if same layer, compare positions of holes, later holes are sorted earlier
    return p1.Holes() < p2.Holes();
    // if same hole positions:
  } else {
    for (size_t i = 0; i < p1.Length(); i++) { // compare mat IDs in non-hole layers
      if (!p1.Holes()[i]) {
        if (p1.Mats()[i]!=p2.Mats()[i]) {
          return p1.Mats()[i] < p2.Mats()[i];
        }
      }
    }
    for (size_t i = 0; i < p1.Length(); i++) { // compare mat IDs in yes-hole layers
      if (p1.Holes()[i]) {
        if (p1.Mats()[i]!=p2.Mats()[i]) {
          return p1.Mats()[i] < p2.Mats()[i];
        }
      }
    }
    // to do: select the representative pattern from the holes family here
    return false; // they are equal
  }
}

inline std::ostream &operator<<(std::ostream &out, const TriggerPattern &p1) {
  out << "{First layer=" << p1.Layer() << ", Mats=(";
  for (size_t iMat=0; iMat<p1.Length()-1; iMat++) {
    out << p1.Mats()[iMat];
    if (p1.Holes()[iMat]) out << "*";
    out << ",";
  }
  out << p1.Mats()[p1.Length()-1] << ")} " << '\n';
  out << "ID: " << p1.FullID();
  return out;
}

template<>
struct std::hash<TriggerPattern> {
    std::size_t operator()(TriggerPattern const& pattern) const noexcept {
      // std::cout << "hashing " + pattern.FullID() << '\n';
      return std::hash<std::string>{}(pattern.FullID());
      // std::cout << "hash generated: " + std::to_string(hash) << '\n';
    }
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Search recursively for hits satisfying the condition of belonging to a pattern defined by
// <nHits>-hits-in-<nLayers>-layers, allowing a maximum displacement of <dM> mats betweeen consecutive layers
// <FIThits>  = <nLayers> x <nMats> matrix of level-0 hits in the FIT module
// <iPattern> = partial pattern until the current level
// patterns   = vector to store all found complete patterns

// using MatsArray <ContainedClass> = std::vector<std::vector<ContainedClass>>>;


//


// int main() {
//   int _nHits = 3;
//   int _patternSize = 3;
//   int _dM = 1;
//
//   FitHitsArray hitsArray = // 7 layers x 10 mats
//    {{0,0,0,0,0,0,0,0,0,0},
//     {1,0,0,0,0,1,0,0,0,0},
//     {0,0,1,0,1,0,0,0,1,0},
//     {0,1,0,1,0,1,0,0,0,0},
//     {0,0,0,0,0,1,0,0,0,1},
//     {1,0,1,0,0,0,1,0,0,1},
//     {0,1,0,0,0,0,0,1,0,1}};
//   Herd::RefFrame::Side side=Herd::RefFrame::Side::TOP;
//   int nLayers = hitsArray.size();
//   int nMats = hitsArray[0].size();
//   std::vector<TriggerPattern> vTriggers;
//   bool force_full_spanning=true;
//
//   for (size_t iLayer = 0; iLayer < nLayers-(_patternSize-1); iLayer++) {
//     for (size_t iMat = 0; iMat < nMats; iMat++) {
//       // std::cout << "mat: " << iMat << '\n';
//       // std::cout << hitsArray[iLayer][iMat] << '\n';
//       int hasHit = hitsArray[iLayer][iMat];
//       if (hasHit || _patternSize>_nHits) {
//         TriggerPattern iPattern = TriggerPattern(side,2*iLayer,iMat,!hasHit);
//         findPatterns(_nHits-hasHit,_patternSize,_dM,hitsArray,iPattern,vTriggers,force_full_spanning);
//       }
//     }
//   }
//   std::cout << "Raw patterns: " << vTriggers.size() << '\n';
//   for (auto patt : vTriggers) {
//     std::cout << patt << '\n';
//     std::cout << std::stoi(patt.PatternID()) << '\n';
//     std::cout << patt.OldID() << '\n';
//   }
//   if (_patternSize>_nHits && vTriggers.size()>0) {
//     std::cout << "Sorted and cleaned patterns: " << '\n';
//     std::sort(vTriggers.begin(),vTriggers.end());
//     vTriggers.erase(std::unique(vTriggers.begin(),vTriggers.end()),vTriggers.end());
//     for (auto patt : vTriggers) {
//       std::cout << patt << '\n';
//       std::cout << std::stoi(patt.PatternID()) << '\n';
//       std::cout << patt.OldID() << '\n';
//     }
//   }
//   return 0;
// }



#endif /* TriggerPattern_H_ */
