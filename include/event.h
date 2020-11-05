//////////////////////////////////////////////////////////
// EVENT Class
// J.C. Zamora, cardona@if.usp.br
// University of Sao Paulo, 2020
//////////////////////////////////////////////////////////
#ifndef Event_H
#define Event_H

#include <TObject.h>
#include "TVector3.h"
#include "TMath.h"

#include <numeric>
#include <algorithm>
#include <iostream>


class Event : public TObject {
public:
  Event();
  ~Event();


  void HelloWorld();
  void SetEventCharge(double charge);
  double GetEventCharge(){return fEventCharge;};
  void SetEventNumTracks(int ntracks);
  int GetEventNumTracks(){return fEventNumTracks;};
  void SetEventRvertex(TVector3 rvertex);
  TVector3 GetEventRvertex(){return fEventRvertex;};

protected:

  double fEventCharge;
  int fEventNumTracks;
  TVector3 fEventRvertex;

private:

public:
  ClassDef(Event,0)
};

#endif
