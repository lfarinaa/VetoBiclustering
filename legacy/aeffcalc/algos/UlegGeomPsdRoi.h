/*
* UlegGeomPsdRoi.h
*
*  Created on: 12 Oct 2019
*      Author: Luis Fariña
*/

#ifndef UlegGeomPsdRoi_H_
#define UlegGeomPsdRoi_H_


// C/C++ standard headers
#include <iterator>
#include <numeric>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <unordered_map>


// #include "algorithm/Algorithm.h"
#include "herd/analysis/dataobjects/FitGeoParamsColl.h"
#include "herd/analysis/dataobjects/PsdGeoParamsColl.h"
#include "herd/analysis/dataobjects/CaloGeoParams.h"
#include "herd/analysis/dataobjects/PsdHitsColl.h"
// #include "herd/analysis/dataobjects/ScdGeoParamsColl.h"
#include "herd/analysis/dataobjects/Line.h"
#include "herd/analysis/dataobjects/FitHitsColl.h"
#include "herd/analysis/dataobjects/CaloHits.h"
#include "TriggerPattern.h"




class UlegGeomRoi {

private:
  Herd::RefFrame::Side _side;
  Herd::Vec3D _fitNormalVec;
  Herd::Vec3D _fitTransVec;
  Herd::Vec3D _fitLongVec;
  Herd::Point _crossPoint;
  float _angleLeft;
  float _angleRight;
  float _firstMatHeight;
  float _lastMatHeight;
  // IDs of the tiles in the patern's RoI.

  void BuildFrom(const TriggerPattern &triggerPattern, const Herd::FitGeoParamsColl &fitGeoParamsColl) {
    // auto x=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000;

    // std::cout << "Building" << triggerPattern << '\n';

    // void BuildFrom(TriggerPattern triggerPattern, Herd::FitGeoParamsColl &fitGeoParamsColl, Herd::PsdGeoParamsColl &psdGeoParamsColl) {

    _angleLeft = 6;
    _angleRight = -6;
    _side = triggerPattern.Side();
    auto firstLayer = triggerPattern.Layer();
    auto lastLayer = triggerPattern.Layer() + 2*(triggerPattern.Length() - 1);
    auto firstMat = triggerPattern.Mats().front();
    auto lastMat = triggerPattern.Mats().back();

    auto firstFitGeoParams = fitGeoParamsColl[_side][firstLayer];
    auto lastFitGeoParams = fitGeoParamsColl[_side][lastLayer];
    auto fitNormalDir = firstFitGeoParams.NormalDirection();
    auto fitNormalCoo = Herd::RefFrame::ToCoo(fitNormalDir);
    _fitNormalVec = Herd::Vec3D(fitNormalDir);
    auto fitTransDir = firstFitGeoParams.SegmentationDirections().first;
    auto fitLongDir = firstFitGeoParams.SegmentationDirections().second;
    // auto mainSegCoo = Herd::RefFrame::ToCoo(fitTransDir);
    // auto secSegCoo = Herd::RefFrame::ToCoo(fitLongDir);
    _fitTransVec = Herd::Vec3D(fitTransDir);
    _fitLongVec = Herd::Vec3D(fitLongDir);
    auto matDimensions = Herd::Vec3D(firstFitGeoParams.MatSize()[Herd::RefFrame::Axis::X],firstFitGeoParams.MatSize()[Herd::RefFrame::Axis::Y],firstFitGeoParams.MatSize()[Herd::RefFrame::Axis::Z]);
    // matDimensions[secSegCoo]=0;
    auto deltaMat = lastMat - firstMat;
    auto firstMatPos = firstFitGeoParams.MatPosition(firstMat);
    auto lastMatPos = lastFitGeoParams.MatPosition(lastMat);
    _firstMatHeight = firstMatPos*_fitNormalVec-0.5*firstFitGeoParams.MatSize()[Herd::RefFrame::ToAxis(fitNormalCoo)];
    _lastMatHeight = lastMatPos*_fitNormalVec+0.5*lastFitGeoParams.MatSize()[Herd::RefFrame::ToAxis(fitNormalCoo)];
    // std::cout << "First mat at " << firstMatPos << '\n';

    // this part requires the swap that follows. There are probably more elegant ways of writing it. Maybe a minmax_element() if more legible?
    auto leftFaceFirstMat=firstMatPos-matDimensions/2.*_fitTransVec*_fitTransVec;
    // std::cout << "left1: " << leftEdgeFirstMat << '\n';
    auto rightFaceFirstMat=firstMatPos+matDimensions/2.*_fitTransVec*_fitTransVec;
    // std::cout << "right1: " << rightEdgeFirstMat << '\n';
    auto leftFaceLastMat=lastMatPos-matDimensions/2.*_fitTransVec*_fitTransVec;
    // std::cout << "left2: " << leftEdgeLastMat << '\n';
    auto rightFaceLastMat=lastMatPos+matDimensions/2.*_fitTransVec*_fitTransVec;
    if (!Herd::RefFrame::IsPositive(fitTransDir)) {
      std::swap(leftFaceFirstMat,rightFaceFirstMat);
      std::swap(leftFaceLastMat,rightFaceLastMat);
    }
    // std::cout << "right2: " << rightEdgeLastMat << '\n';
    auto halfMatThickness=matDimensions/2.*_fitNormalVec*_fitNormalVec;
    std::vector<Herd::Vec3D> vUpDown{-1.*halfMatThickness,halfMatThickness};


    Herd::Line leftLine, rightLine;
    for (auto leftUpDown : vUpDown) {
      for (auto rightUpDown : vUpDown) {
        // std::cout << leftUpDown << '\n';
        // for (auto leftEdgeLastMat : {leftFaceLastMat-halfMatThickness,leftFaceLastMat+halfMatThickness}) {
        // for (auto rightEdgeLastMat : {rightFaceLastMat-halfMatThickness,rightFaceLastMat+halfMatThickness}) {
        auto leftEdgeFirstMat=leftFaceFirstMat+leftUpDown;
        // std::cout << leftEdgeFirstMat << '\n';
        auto rightEdgeFirstMat=rightFaceFirstMat+rightUpDown;
        auto leftEdgeLastMat=leftFaceLastMat+leftUpDown;
        auto rightEdgeLastMat=rightFaceLastMat+rightUpDown;

        auto leftLineTemp = Herd::Line(leftEdgeFirstMat,leftEdgeFirstMat-rightEdgeLastMat);
        // std::cout << "left: " << leftLineTemp.Origin() << " + t" << leftLineTemp.Direction() << '\n';
        auto rightLineTemp = Herd::Line(rightEdgeFirstMat,rightEdgeFirstMat-leftEdgeLastMat);
        // std::cout << "right: " << rightLineTemp.Origin() << " + t" << rightLineTemp.Direction() << '\n';
        auto angleLeftTemp = std::atan2(leftLineTemp.Direction()*_fitTransVec,leftLineTemp.Direction()*_fitNormalVec);
        // std::cout << "cosL: " << angleLeftTemp << '\n';
        auto angleRightTemp = std::atan2(rightLineTemp.Direction()*_fitTransVec,rightLineTemp.Direction()*_fitNormalVec);
        // std::cout << "cosR: " << angleRightTemp << '\n';
        // std::cout << (angleLeft<angleRight) << '\n';
        if (angleLeftTemp<_angleLeft) {
          leftLine=leftLineTemp;
          _angleLeft=angleLeftTemp;
        }
        if (angleRightTemp>_angleRight) {
          rightLine=rightLineTemp;
          _angleRight=angleRightTemp;
        }
        // vLeftLines.push_back(leftLine);
        // vRightLines.push_back(rightLine);
        // }
        // }
      }
    }
    // auto leftLine = *std::min_element(vLeftLines.begin(),vLeftLines.end(),[_fitTransVec,_fitNormalVec](Herd::Line l1, Herd::Line l2) {return l1.Direction()*_fitTransVec/(l1.Direction()*_fitTransVec*_fitTransVec+l1.Direction()*_fitNormalVec*_fitNormalVec).Mag()<l2.Direction()*_fitTransVec/(l2.Direction()*_fitTransVec*_fitTransVec+l2.Direction()*_fitNormalVec*_fitNormalVec).Mag();});
    // auto rightLine = *std::max_element(vRightLines.begin(),vRightLines.end(),[_fitTransVec,_fitNormalVec](Herd::Line l1, Herd::Line l2) {return l1.Direction()*_fitTransVec/(l1.Direction()*_fitTransVec*_fitTransVec+l1.Direction()*_fitNormalVec*_fitNormalVec).Mag()<l2.Direction()*_fitTransVec/(l2.Direction()*_fitTransVec*_fitTransVec+l2.Direction()*_fitNormalVec*_fitNormalVec).Mag();});
    // std::cout << "left: " << leftLine.Origin() << " + t" << leftLine.Direction() << '\n';
    // std::cout << "right: " << rightLine.Origin() << " + t" << rightLine.Direction() << '\n';
    _crossPoint = leftLine.Intersection(rightLine).value();
    // auto y=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000-x;

    // std::cout << "geomRoi: " << y << '\n';
    // _angleLeft = leftLine.Direction()*_fitTransVec/(leftLine.Direction()*_fitTransVec*_fitTransVec+leftLine.Direction()*_fitNormalVec*_fitNormalVec).Mag();
    // std::cout << "cosL: " << _angleLeft << '\n';
    // _angleRight = rightLine.Direction()*_fitTransVec/(rightLine.Direction()*_fitTransVec*_fitTransVec+rightLine.Direction()*_fitNormalVec*_fitNormalVec).Mag();


    // std::cout << "cosR: " << _angleRight << '\n';
    // std::cout << (angleLeft<angleRight) << '\n';
    // _dCosPhiLR = angleRight-angleLeft;
  }
public:
  UlegGeomRoi(const TriggerPattern& triggerPattern, const Herd::FitGeoParamsColl& fitGeoParamsColl) : _side{}, _crossPoint{}, _fitLongVec{}, _fitTransVec{}, _angleLeft{}, _angleRight{} {BuildFrom(triggerPattern, fitGeoParamsColl);}
  // TriggerPattern(Herd::RefFrame::Side side, int layer, int mat, int hole) :  _side{side}, _layer{layer}, _mats{mat}, _holes{hole} {}
  // copy init:
  // UlegGeomRoi(const UlegGeomRoi &roi) : _tilesInRoi{roi.TileIDs()} {}
  // destroy:
  // virtual ~UlegGeomRoi() = default;

  // Herd::Side _side;
  Herd::RefFrame::Side Side() const {return _side;}
  // Herd::Vec3D _fitTransVec;
  Herd::Vec3D TransVec() const {return _fitTransVec;}
  Herd::Vec3D NormalVec() const {return _fitNormalVec;}
  Herd::Vec3D LongVec() const {return _fitLongVec;}
  // Herd::Vec3D _fitLongVec;
  // Herd::Point _crossPoint;
  Herd::Point CrossPoint() const {return _crossPoint;}
  std::pair<float,float> Angles() const {return std::pair<float,float>(_angleLeft,_angleRight);}

  float AngleFromCrossPoint(Herd::Point testPos) const {
    auto vecDirCross = testPos-_crossPoint;
    // vecDirCross /= vecDirCross.Mag();
    auto angle = std::atan2(vecDirCross*_fitTransVec,vecDirCross*_fitNormalVec);
    // auto anglePos = vecDirCross*_fitTransVec/(vecDirCross*_fitTransVec*_fitTransVec+vecDirCross*_fitNormalVec*_fitNormalVec).Mag();
    return angle;
  }
  // float _angleLeft;
  // float _angleRight;
  bool IsHigher(Herd::Point testPos) const {
    auto posHeight = testPos*_fitNormalVec;
    return (posHeight>_firstMatHeight);
  }

  bool IsLower(Herd::Point testPos) const {
    auto posHeight = testPos*_fitNormalVec;
    return (posHeight<_lastMatHeight);
  }

  bool IsInDirectRoi(Herd::Point testPos) const {
    auto anglePos = AngleFromCrossPoint(testPos);
    auto isInDirectRoi = (anglePos>_angleLeft && anglePos<_angleRight);
    // auto reverseAnglePos = anglePos - std::copysign(3.14159265,anglePos); // angle in (-pi,pi). The opposite is obtained by adding/subtracting
    // auto isInReverseRoi = (reverseAnglePos>_angleLeft && reverseAnglePos<_angleRight);
    return (isInDirectRoi && IsHigher(testPos));
  }

  bool IsInReverseRoi(Herd::Point testPos) const {
    auto anglePos = AngleFromCrossPoint(testPos);
    // auto isInDirectRoi = (anglePos>_angleLeft && anglePos<_angleRight);
    auto reverseAnglePos = anglePos - std::copysign(3.14159265,anglePos); // angle in (-pi,pi). The opposite is obtained by adding/subtracting
    auto isInReverseRoi = (reverseAnglePos>_angleLeft && reverseAnglePos<_angleRight);
    return (isInReverseRoi && IsLower(testPos));
  }



  bool IsInRoi(Herd::Point testPos) const {
    auto anglePos = AngleFromCrossPoint(testPos);
    auto isInDirectRoi = (anglePos>_angleLeft && anglePos<_angleRight);
    auto reverseAnglePos = anglePos - std::copysign(3.14159265,anglePos); // angle in (-pi,pi). The opposite is obtained by adding/subtracting
    auto isInReverseRoi = (reverseAnglePos>_angleLeft && reverseAnglePos<_angleRight);
    return ((isInDirectRoi && IsHigher(testPos)) || (isInReverseRoi && IsLower(testPos)));
  }


  void LowerCrossPoint(float dist) {
    _crossPoint-=dist*_fitNormalVec;
  }
  void OpenAnglesByFactor(float factor) {
    _angleLeft*=factor;
    _angleRight*=factor;
  }
  void OpenAnglesByConst(float factor) {
    _angleLeft+=std::copysign(factor,_angleLeft);
    _angleRight+=std::copysign(factor,_angleRight);
  }

  bool MinimumCubes(const Herd::CaloGeoParams& caloGeoParams, const int minCubes) const {
    // auto cubeSize = caloGeoParams.CubeSize();
    // auto vPlusMinus = {-cubeSize/2.,cubeSize/2.};
    int nCubes=0;
    for (size_t cubeID = 0; cubeID < caloGeoParams.NCubes(); cubeID++) {
      auto caloHitPos=caloGeoParams.Position(cubeID);
      // bool isInRoi = false;
      // for (auto dX : vPlusMinus) {
      //   for (auto dY : vPlusMinus) {
      //     for (auto dZ : vPlusMinus) {
      //       auto dPos=Herd::Point(dX,dY,dZ);
      //       auto testPos = caloHitPos + dPos;
      //       isInRoi=isInRoi || geomRoi.IsInRoi(testPos);
      //     }
      //   }
      // }
      // uncomment to check cube centers, but code is slow then. Better add early stopping or smarter stuff
      auto isInRoi=IsInRoi(caloHitPos);
      if (isInRoi) {
        nCubes++;
        if (nCubes>minCubes) return true;
      }
    } // loop over cubes
    return false;
  }

};

class PsdVetoList {
private:
  Herd::SidesArray<std::vector<std::vector<int>>> _tilesInRoi;

  void CalculateTileList(const UlegGeomRoi& geomRoi, const Herd::PsdGeoParamsColl &psdGeoParamsColl) {
    for (auto psdSide : Herd::RefFrame::Sides) {
      // std::cout << "psd psdSide: " << psdSide << '\n';
      auto psdGeoParamsVector = psdGeoParamsColl[psdSide];
      for (size_t iPsdLayer=0; iPsdLayer<psdGeoParamsVector.size();iPsdLayer++) {
        std::vector<int> tilesInLayer;
        auto psdGeoParams = psdGeoParamsVector[iPsdLayer];
        auto tileDimensions = Herd::Vec3D(psdGeoParams.Sizes()[Herd::RefFrame::Axis::X],psdGeoParams.Sizes()[Herd::RefFrame::Axis::Y],psdGeoParams.Sizes()[Herd::RefFrame::Axis::Z]);
        auto psdNormalVec = Herd::Vec3D(psdGeoParams.NormalDirection());
        auto psdMainVec = Herd::Vec3D(psdGeoParams.SegmentationDirections().first);
        auto psdSecVec = Herd::Vec3D(psdGeoParams.SegmentationDirections().second);
        for (size_t psdTileID = 0; psdTileID < psdGeoParams.NTotElements(); psdTileID++) { // rewrite with an std magic (foreach?) with early return
          auto psdHitPos=psdGeoParams.Position(psdTileID);
          bool isInRoi = false;
          for (auto dN : {-1*tileDimensions/2.*psdNormalVec*psdNormalVec,tileDimensions/2.*psdNormalVec*psdNormalVec}) {
            for (auto dM : {-1*tileDimensions/2.*psdMainVec*psdMainVec,tileDimensions/2.*psdMainVec*psdMainVec}) {
              for (auto dS : {-1*tileDimensions/2.*psdSecVec*psdSecVec,tileDimensions/2.*psdSecVec*psdSecVec}) {
                auto dPos = dN + dM + dS;
                auto testPos = psdHitPos + dPos;
                isInRoi = isInRoi || geomRoi.IsInDirectRoi(testPos);
                // isInRoi = isInRoi || (geomRoi.IsInRoi(testPos) && (twoSided || geomRoi.IsHigher(testPos)));
              }
            }
          }
          if (isInRoi) {
            tilesInLayer.push_back(psdTileID);
          }
        } // loop over tiles
        _tilesInRoi[psdSide].push_back(tilesInLayer);
      }
    }
  }



public:
  PsdVetoList(const UlegGeomRoi& geomRoi, const Herd::PsdGeoParamsColl& psdGeoParamsColl) : _tilesInRoi{} {CalculateTileList(geomRoi, psdGeoParamsColl);}
  // TriggerPattern(Herd::RefFrame::Side side, int layer, int mat, int hole) :  _side{side}, _layer{layer}, _mats{mat}, _holes{hole} {}
  // copy init:
  // PsdVetoList(const PsdVetoList &tileList) : _tilesInRoi{tileList.TileIDs()} {}
  PsdVetoList() : _tilesInRoi{} {}

  PsdVetoList(const Herd::PsdGeoParamsColl& psdGeoParamsColl) : _tilesInRoi{} {
    for (auto side : Herd::RefFrame::Sides) {
      for (auto psdGeoParams : psdGeoParamsColl[side]) {
        _tilesInRoi[side].push_back({}); // init to shape determined by the geoparams
      }
    }
  }
  PsdVetoList(const Herd::SidesArray<std::vector<std::vector<int>>>&& tileList) : _tilesInRoi{std::move(tileList)} {}
  // destroy:
  // virtual ~PsdVetoList() = default;

  Herd::SidesArray<std::vector<std::vector<int>>> TileIDs() const {return _tilesInRoi;}

  bool IsInRoi(int psdTileID, Herd::RefFrame::Side psdSide, int iPsdLayer) const {
    return std::find(_tilesInRoi[psdSide][iPsdLayer].begin(),_tilesInRoi[psdSide][iPsdLayer].end(),psdTileID)!=_tilesInRoi[psdSide][iPsdLayer].end();
  }

  bool Veto(Herd::Hit psdHit, Herd::RefFrame::Side psdSide, int iPsdLayer, float energyThreshold) const {
    return (IsInRoi(psdHit.VolumeID(), psdSide, iPsdLayer) && psdHit.EDep()>energyThreshold);
  }

  bool Veto(Herd::PsdHitsColl psdHitsColl, float energyThreshold) const {
    for (auto psdSide : Herd::RefFrame::Sides) {
      for (size_t iPsdLayer=0; iPsdLayer<psdHitsColl[psdSide].size(); iPsdLayer++) {
        for (auto psdHit : psdHitsColl[psdSide][iPsdLayer]) {
          if (IsInRoi(psdHit.VolumeID(), psdSide, iPsdLayer) && psdHit.EDep()>energyThreshold) return true;
        }
      }
    }
    return false;
  }

  int NTotalTiles() const {
    int acc=0;
    for (auto tilesSide : _tilesInRoi) {
      for (auto tilesLayer : tilesSide) {
        acc += tilesLayer.size();
      }
    }
    return acc;
  }

  Herd::SidesArray<int> TilesPerSide() const {
    Herd::SidesArray<int> tiles;
    for (auto psdSide : Herd::RefFrame::Sides) {
      int acc=0;
      for (auto tilesLayer : _tilesInRoi[psdSide]) {
        acc += tilesLayer.size();
      }
      tiles[psdSide]=acc;
    }
    return tiles;
  }

  PsdVetoList Intersection(const PsdVetoList& other) {
    Herd::SidesArray<std::vector<std::vector<int>>> intersection;
    auto &otherTilesList = other._tilesInRoi;
    // Herd::SidesArray<std::vector<std::vector<int>>> otherTilesList = other.TileIDs();
    for (auto psdSide : Herd::RefFrame::Sides) {
      for (size_t iPsdLayer = 0; iPsdLayer < _tilesInRoi[psdSide].size(); iPsdLayer++) {
        std::vector<int> layerIntersection;
        auto &layerThis=_tilesInRoi[psdSide][iPsdLayer];
        auto &layerOther=otherTilesList[psdSide][iPsdLayer];
        std::set_intersection(layerThis.begin(), layerThis.end(), layerOther.begin(), layerOther.end(), std::back_inserter(layerIntersection));
        intersection[psdSide].push_back(std::move(layerIntersection));
      }
    }
    return PsdVetoList(std::move(intersection));
  }

  PsdVetoList Union(const PsdVetoList& other) {
    Herd::SidesArray<std::vector<std::vector<int>>> vUnion;
    auto &otherTilesList = other._tilesInRoi;
    for (auto psdSide : Herd::RefFrame::Sides) {
      for (size_t iPsdLayer = 0; iPsdLayer < _tilesInRoi[psdSide].size(); iPsdLayer++) {
        std::vector<int> layerUnion;
        auto &layerThis=_tilesInRoi[psdSide][iPsdLayer];
        auto &layerOther=otherTilesList[psdSide][iPsdLayer];
        std::set_union(layerThis.begin(), layerThis.end(), layerOther.begin(), layerOther.end(), std::back_inserter(layerUnion));
        vUnion[psdSide].push_back(std::move(layerUnion));
      }
    }
    return PsdVetoList(std::move(vUnion));
  }

  float Overlap(const PsdVetoList& other) {
    float sizeOfIntersection = Intersection(other).NTotalTiles();
    float minSize = std::min(NTotalTiles(),other.NTotalTiles());
    return sizeOfIntersection/minSize;
  }

  float Jaccard(const PsdVetoList& other) {
    float sizeOfIntersection = Intersection(other).NTotalTiles();
    float sizeOfUnion = Union(other).NTotalTiles();
    return sizeOfIntersection/sizeOfUnion;
  }
};

class CaloVetoList {
private:
  std::vector<int> _cubesInRoi;
  void CalculateCubeList(const UlegGeomRoi& geomRoi, const Herd::CaloGeoParams& caloGeoParams) noexcept {
    // auto x=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000;


    auto cubeSize = caloGeoParams.CubeSize();
    auto vPlusMinus = {-cubeSize/2.,cubeSize/2.};
    for (size_t cubeID = 0; cubeID < caloGeoParams.NCubes(); cubeID++) {
      auto caloHitPos=caloGeoParams.Position(cubeID);
      // bool isInRoi = false;
      // for (auto dX : vPlusMinus) {
      //   for (auto dY : vPlusMinus) {
      //     for (auto dZ : vPlusMinus) {
      //       auto dPos=Herd::Point(dX,dY,dZ);
      //       auto testPos = caloHitPos + dPos;
      //       isInRoi=isInRoi || geomRoi.IsInRoi(testPos);
      //     }
      //   }
      // }
      // uncomment to check cube centers, but code is slow then. Better add early stopping or smarter stuff
      auto isInRoi=geomRoi.IsInRoi(caloHitPos);
      // if (geomRoi.Side()==Herd::RefFrame::Side::TOP) {
        // std::cout << "caloHitPos: " << caloHitPos << isInRoi << '\n';
      // }
      if (isInRoi) {
        _cubesInRoi.push_back(cubeID);
      }
    } // loop over cubes
    // auto y=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000-x;
    // std::cout << "caloList: " << y << '\n';
  }
public:
  CaloVetoList(const UlegGeomRoi& geomRoi, const Herd::CaloGeoParams& caloGeoParams) : _cubesInRoi{} {CalculateCubeList(geomRoi, caloGeoParams);}
  CaloVetoList() : _cubesInRoi{} {}
    // copy init:
    // CaloVetoList(const CaloVetoList &cubeList) : _cubesInRoi{cubeList.CubeIDs()} {}
    // CaloVetoList(const std::vector<int>& cubeList) : _cubesInRoi{cubeList} {}
    CaloVetoList(const std::vector<int>&& cubeList) : _cubesInRoi{std::move(cubeList)} {}
    // destroy:
    // virtual ~CaloVetoList() = default;

    std::vector<int> CubeIDs() const {return _cubesInRoi;}
    int NTotalCubes() const {return _cubesInRoi.size();}
    bool MinimumCubes(const UlegGeomRoi& geomRoi, const Herd::CaloGeoParams& caloGeoParams, const int minCubes) const {
      // auto cubeSize = caloGeoParams.CubeSize();
      // auto vPlusMinus = {-cubeSize/2.,cubeSize/2.};
      int nCubes=0;
      for (size_t cubeID = 0; cubeID < caloGeoParams.NCubes(); cubeID++) {
        auto caloHitPos=caloGeoParams.Position(cubeID);
        // bool isInRoi = false;
        // for (auto dX : vPlusMinus) {
        //   for (auto dY : vPlusMinus) {
        //     for (auto dZ : vPlusMinus) {
        //       auto dPos=Herd::Point(dX,dY,dZ);
        //       auto testPos = caloHitPos + dPos;
        //       isInRoi=isInRoi || geomRoi.IsInRoi(testPos);
        //     }
        //   }
        // }
        // uncomment to check cube centers, but code is slow then. Better add early stopping or smarter stuff
        auto isInRoi=geomRoi.IsInRoi(caloHitPos);
        if (isInRoi) {
          nCubes++;
          if (nCubes>minCubes) return true;
        }
      } // loop over cubes
      return false;
    }

    bool IsInRoi(int caloCubeID) const {
      return std::find(_cubesInRoi.begin(),_cubesInRoi.end(),caloCubeID)!=_cubesInRoi.end();
    }

    bool Veto(Herd::Hit caloHit, float energyThreshold) const {
      return (IsInRoi(caloHit.VolumeID()) && caloHit.EDep()>energyThreshold);
    }

    bool Veto(Herd::CaloHits caloHits, float energyThreshold) const {
      for (auto caloHit : caloHits) {
        if (IsInRoi(caloHit.VolumeID()) && caloHit.EDep()>energyThreshold) return true;
      }
      return false;
    }

    CaloVetoList Intersection(const CaloVetoList& other) noexcept {
      // auto x=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000;


      std::vector<int> intersection;
      auto &otherCubesList = other._cubesInRoi;
      std::set_intersection(_cubesInRoi.begin(), _cubesInRoi.end(), otherCubesList.begin(), otherCubesList.end(), std::back_inserter(intersection));
      // auto y=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000-x;
      // std::cout << "caloIntersection: " << y << '\n';
      return CaloVetoList(std::move(intersection));
      // return CaloVetoList(intersection);
    }

    bool IntersectionHasCubes(const CaloVetoList& other, const int minCubes) const {
      // auto cubeSize = caloGeoParams.CubeSize();
      // auto vPlusMinus = {-cubeSize/2.,cubeSize/2.};
      // auto x=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000;
      int nCubesInIntersection=0;
      auto &otherCubesList = other._cubesInRoi;
      size_t minCubeID=0;
      for (auto Cube : _cubesInRoi) {
        for (size_t jCube=minCubeID; jCube < otherCubesList.size(); jCube++) {
          if (Cube == otherCubesList[jCube]) {
            nCubesInIntersection++;
            minCubeID=jCube+1;
          } else if (Cube < otherCubesList[jCube]) break;
          if (nCubesInIntersection>=minCubes) {
            // auto y=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000-x;

            // std::cout << "caloIntersectionHasCubes: " << y << '\n';

            return true;
          }
        }
      }
      // std::set_intersection(_cubesInRoi.begin(), _cubesInRoi.end(), otherCubesList.begin(), otherCubesList.end(), [](){nCubes>=minCubes?nCubes++:return});
      // auto y=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000-x;
      // std::cout << "caloIntersectionHasCubes: " << y << '\n';
      // return CaloVetoList(std::move(intersection));
      // return CaloVetoList(intersection);
      return false;
    }

    CaloVetoList Union(const CaloVetoList& other) noexcept {
      // auto x=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000;

      std::vector<int> vUnion;
      auto &otherCubesList = other._cubesInRoi;
      std::set_union(_cubesInRoi.begin(), _cubesInRoi.end(), otherCubesList.begin(), otherCubesList.end(), std::back_inserter(vUnion));

      // auto y=(unsigned long long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()%1000000-x;
      // std::cout << "caloUnion: " << y << '\n';
      // return CaloVetoList(vUnion);
      return CaloVetoList(std::move(vUnion));
    }
};

inline std::ostream &operator<<(std::ostream &out, const PsdVetoList &psdVetoList) {
  auto tilesPerSide = psdVetoList.TilesPerSide();
  out << "Tiles in RoI: " << psdVetoList.NTotalTiles() << "\n";
  if (psdVetoList.NTotalTiles()>0) {
    for (auto psdSide : Herd::RefFrame::Sides) {
      out << "PSD side " << psdSide << ": " << tilesPerSide[psdSide] << " tiles" << ": \n";
      if (tilesPerSide[psdSide]>0) {
        auto tilesSide = psdVetoList.TileIDs()[psdSide];
        for (size_t iPsdLayer=0; iPsdLayer<tilesSide.size(); iPsdLayer++) {
          auto tileIDs = tilesSide[iPsdLayer];
          if (tileIDs.size()>0) {
            out << "Layer " << iPsdLayer << ": ";
            std::string separator = "(";
            for (auto tileID : tileIDs) {
              out << separator << tileID;
              separator = ",";
            }
            out << ")\n";
          }
        }
      }
    }
  }
  return out;
}

inline std::ostream &operator<<(std::ostream &out, const CaloVetoList &caloVetoList) {
  out << "Cubes in RoI: " << caloVetoList.NTotalCubes() << "\n";
    if (caloVetoList.NTotalCubes()>0) {
        auto cubeIDs = caloVetoList.CubeIDs();
          std::string separator = "(";
          for (auto cubeID : cubeIDs) {
            out << separator << cubeID;
            separator = ",";
          }
          out << ")\n";
        }
  return out;
}


template<typename SubdetectorVetoList, typename GeoContainer>
class PatternVetoListMap {
private:
  std::unordered_map<TriggerPattern,SubdetectorVetoList> _map;
public:
  PatternVetoListMap() : _map{} {}

  SubdetectorVetoList& ConstructRetrieve(TriggerPattern pattern, Herd::FitGeoParamsColl& fitGeoParamsColl, GeoContainer& geo) {
    if (!_map.count(pattern)) { // .contains() in C++20
      auto geomRoi = UlegGeomRoi(pattern,fitGeoParamsColl);
      _map[pattern] = SubdetectorVetoList(geomRoi, geo);
    }
    return _map.at(pattern);
  }
};

// template<typename SubdetectorVetoList>
// class PatternVetoListMap {
// private:
//   std::unordered_map<TriggerPattern,SubdetectorVetoList> _map;
//   std::shared_ptr<void> _geo;
// public:
//   template<typename GeoContainer>
//   PatternVetoListMap(GeoContainer& geo) : _map{}, _geo(std::make_shared<GeoContainer>(geo)) {}
//
//   SubdetectorVetoList& ConstructRetrieve(TriggerPattern pattern, Herd::FitGeoParamsColl& fitGeoParamsColl) {
//     if (!_map.count(pattern)) {
//       auto geomRoi = UlegGeomRoi(pattern,fitGeoParamsColl);
//       _map[pattern] = SubdetectorVetoList(geomRoi, *static_cast<GeoContainer*>(_geo.get()));
//     }
//     return _map.at(pattern);
//   }
// };


#endif /* UlegGeomPsdRoi_H_ */
