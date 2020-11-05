//////////////////////////////////////////////////////////
// TRACK Class
// J.C. Zamora, cardona@if.usp.br
// University of Sao Paulo, 2020
//////////////////////////////////////////////////////////
#ifndef Track_H
#define Track_H

#include <TObject.h>
#include "TVector3.h"
#include "TMath.h"

#include <numeric>
#include <algorithm>
#include <iostream>


class Track : public TObject {
public:
  Track();
  ~Track();


  void HelloWorld();
  void SetTrackVertex(TVector3 vertex);
  TVector3 GetTrackVertex() { return fTrackVertex;};
  void SetTrackRange(Double_t range);
  double GetTrackRange() {return fTrackRange;};
  void SetTrackAngle(double angle);
  double GetTrackAngle() {return fTrackAngle;};
  void SetTrackCharge(double charge);
  double GetTrackCharge() {return fTrackCharge;};

protected:
  TVector3 fTrackVertex;
  double fTrackRange;
  double fTrackAngle;
  double fTrackCharge;

private:

public:
  ClassDef(Track,0)
};

#endif
