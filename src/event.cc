//////////////////////////////////////////////////////////
// EVENT Class
// J.C. Zamora, cardona@if.usp.br
// University of Sao Paulo, 2020
//////////////////////////////////////////////////////////
#define Event_cc

#include "event.h"
#include<iostream>

using namespace std;

ClassImp(Event)

Event::Event() {}

Event::~Event() {}


/*
 * HelloWorld function called to print out the message
 */
void Event::HelloWorld()
{
  cout << "Hello, Event!" << endl;
}
