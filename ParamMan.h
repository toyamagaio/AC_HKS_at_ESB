/*
  ParamMan.h
*/

#ifndef ParamMan_hh
#define ParamMan_hh 1

#include <iomanip>
#include <sstream>
#include <string>
#include <map>
#include "define.h"
//#include "ParamMan.hh"

class ParamMan
{
private:
//   std::string ParamFileName;
  const char* ParamFileName;
  char* OutputFileName;

public:
//   explicit ParamMan( const std::string & filename );
  ParamMan();
  explicit ParamMan( const char* filename );
  ~ParamMan(){}

private:
  struct ACparam{
    int cid;
    int ud;
    double adcOffset[nAC*2];
    double adcGain[nAC*2];
    double tdcOffset[nAC*2];
    double tdcGain[nAC*2];
  };

  struct WCparam{
    int cid;
    int ud;
    double adcOffset[nWC*2];
    double adcGain[nWC*2];
    double tdcOffset[nWC*2];
    double tdcGain[nWC*2];
  };

  struct ToFparam{
    int cid;
    int ud;
    double adcOffset[nToF*2];
    double adcGain[nToF*2];
    double tdcOffset[nToF*2];
    double tdcGain[nToF*2];
  };

  ACparam AC;
  WCparam WC;
  ToFparam ToF;

public:
  //void SetFileName( const char* ParamFileName & filename )
  void SetFileName( const char* filename )
  { ParamFileName=filename; } 
  void SetAdcOffset( int cid, int seg, int tb, double adcOffset );
  void SetAdcGain(   int cid, int seg, int tb, double adcGain   );
  void SetAdcPeak(   int cid, int seg, int tb, double adcPeak   );
  void SetNpeTune(   int cid, int seg, int tb, double npe       );
  void SetTdcOffset( int cid, int seg, int tb, double tdcOffset );
  void SetTdcGain(   int cid, int seg, int tb, double tdcGain   );
  void SetTimeTune(  int cid, int seg, int tb, double tdcGain   );
  
public:
  bool Initialize( void );
  double GetAdcOffset( int cid, int seg, int tb );
  double GetAdcGain(   int cid, int seg, int tb );
  double GetTdcOffset( int cid, int seg, int tb );
  double GetTdcGain(   int cid, int seg, int tb );
  double npe(          int cid, int seg, int tb, double adc);
  double time(         int cid, int seg, int tb, double tdc);
  void WriteToFile( const char* OutputFileName );
};


#endif
