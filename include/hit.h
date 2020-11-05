//////////////////////////////////////////////////////////
// HIT Class
// J.C. Zamora, cardona@if.usp.br
// University of Sao Paulo, 2020
//////////////////////////////////////////////////////////
#ifndef Hit_H
#define Hit_H

#include <TObject.h>
#include "TVector3.h"
#include "TMath.h"

#include <numeric>
#include <algorithm>
#include <iostream>


class Hit : public TObject {
public:
  Hit();
  ~Hit();


  void HelloWorld();
  void SetHit(int hitID, TVector3 vpos, double charge);
  TVector3 GetHitPosition(){return fHitPosition;};
  double GetHitCharge(){return fHitCharge;};

protected:
  TVector3 fHitPosition;
  double fHitCharge;

private:

public:
  ClassDef(Hit,0)
};

#endif
