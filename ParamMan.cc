/*
  ParamMan.cc
*/


#include "ParamMan.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <map>
using namespace std;

ParamMan::ParamMan( const char* filename )
  : ParamFileName(filename)
{
  ToF.cid=1;
  for(int i=0; i<nToF*2; i++){
    ToF.adcOffset[i]=8000.00;   ToF.adcGain[i]=0.0350;
    ToF.tdcOffset[i]=-9999.00;   ToF.tdcGain[i]=0.0250;
  }

  AC.cid=2;
  for(int i=0; i<nAC*2; i++){
    AC.adcOffset[i]=8000.00;   AC.adcGain[i]=0.0350;
    AC.tdcOffset[i]=-9999.00;   AC.tdcGain[i]=0.0250;
  }

  WC.cid=3;
  for(int i=0; i<nWC*2; i++){
    WC.adcOffset[i]=8000.00;   WC.adcGain[i]=0.0350;
    WC.tdcOffset[i]=-9999.00;   WC.tdcGain[i]=0.0250;
  }
}
////////////////////////////////
//ParamMan::~ParamMan()
//{
//}
///////////////////////////////////
bool ParamMan::Initialize( void )
{
  static const std::string funcname = "ParamMan::Initialize";
  FILE *fp;
  char str[MaxChar];
  
  if((fp=fopen(ParamFileName,"r"))==0){
    std::cerr << "[" << funcname << "]: file open fail" << std::endl;
    return false;
    //exit(-1);
  }
  
  int cid,seg,at,tb;
  double p0, p1;
  
  while(fgets(str,MaxChar,fp)!=0){
    if(str[0]=='#') continue;
    else if(sscanf(str,"%d %d %d %d %lf %lf",&cid,&seg,&at,&tb,
		   &p0,&p1)==6){
      if(at==0){ /* ADC */
	    if(cid==ToF.cid)
	      { ToF.adcOffset[seg-1+nToF*tb]=p0;
	        ToF.adcGain[seg-1+nToF*tb]=p1;}
	    else if(cid==AC.cid)
	      { AC.adcOffset[seg-1+nAC*tb]=p0;
	        AC.adcGain[seg-1+nAC*tb]=p1;}
	    else if(cid==WC.cid)
	      { WC.adcOffset[seg-1+nWC*tb]=p0;
	        WC.adcGain[seg-1+nWC*tb]=p1;}
      }

      else if(at==1){ /* TDC */
	    if(cid==ToF.cid)
	      { ToF.tdcOffset[seg-1+nToF*tb]=p0;
	        ToF.tdcGain[seg-1+nToF*tb]=p1;}
	    else if(cid==AC.cid)
	      { AC.tdcOffset[seg-1+nAC*tb]=p0;
	        AC.tdcGain[seg-1+nAC*tb]=p1;}
	    else if(cid==WC.cid)
	      { WC.tdcOffset[seg-1+nWC*tb]=p0;
	        WC.tdcGain[seg-1+nWC*tb]=p1;}
      }

 	  else{
 	    std::cerr << "[" << funcname << "]: new fail (A) "
 		      << " Cid=" << std::setw(2) << cid
 		      << " at="  << std::setw(2) << at
 		      << " tb=" << std::setw(2) << tb << std::endl;
 	  }
	
    }   /* if(sscanf...) */
  }       /* while(fgets...) */
  fclose(fp);
  
//   std::cout << "[" << funcname << "]: Initialization finished" << std::endl;
  return true;
}

///////////////////////////////////
void ParamMan::SetAdcOffset( int cid, int seg, int tb,
				 double adcOffset )
{
  static const std::string funcname = "ParamMan::SetAdcOffset";
  
  if(cid==AC.cid)    AC.adcOffset[seg-1+nAC*tb]=adcOffset;
  else if(cid==ToF.cid)    ToF.adcOffset[seg-1+nToF*tb]=adcOffset;
  else if(cid==WC.cid)     WC.adcOffset[seg-1+nWC*tb]=adcOffset;
  else   cerr << "[" << funcname << "]: unknown id" << endl;
}
///////////////////////////////////
void ParamMan::SetAdcGain( int cid, int seg, int tb,
			       double adcGain )
{
  static const std::string funcname = "ParamMan::SetAdcGain";
  
  if(cid==AC.cid)          AC.adcGain[seg-1+nAC*tb]=adcGain;
  else if(cid==ToF.cid)    ToF.adcGain[seg-1+nToF*tb]=adcGain;
  else if(cid==WC.cid)     WC.adcGain[seg-1+nWC*tb]=adcGain;
  else   cerr << "[" << funcname << "]: unknown id" << endl;
}

///////////////////////////////////
void ParamMan::SetAdcPeak( int cid, int seg, int tb,
			       double adcPeak )
{
  static const std::string funcname = "ParamMan::SetAdcPeak";
  
  if(adcPeak==0) return;
  else if(cid==AC.cid){
        AC.adcGain[seg-1+nAC*tb]=1.0/(adcPeak-AC.adcOffset[seg-1+nAC*tb]);
         }
  else if(cid==ToF.cid){
        ToF.adcGain[seg-1+nToF*tb]=1.0/(adcPeak-ToF.adcOffset[seg-1+nToF*tb]);
         }
  else if(cid==WC.cid){
        WC.adcGain[seg-1+nWC*tb]=1.0/(adcPeak-WC.adcOffset[seg-1+nWC*tb]);
         }
  else   cerr << "[" << funcname << "]: unknown id" << endl;
}

///////////////////////////////////
void ParamMan::SetNpeTune( int cid, int seg, int tb,
				double npe )
{
  static const std::string funcname = "ParamMan::SetNpeTune";
  
  if(cid==AC.cid )
    AC.adcGain[seg-1+nAC*tb]/=npe;
  else if(cid==ToF.cid )
    ToF.adcGain[seg-1+nToF*tb]/=npe;
  else if(cid==WC.cid )
    WC.adcGain[seg-1+nWC*tb]/=npe;
  else   cerr << "[" << funcname << "]: unknown id" << endl;
}

///////////////////////////////////

double ParamMan::GetAdcOffset( int cid, int seg, int tb )
{
  static const std::string funcname = "ParamMan::GetAdcOffset";

  if(cid==AC.cid )          return AC.adcOffset[seg-1+nAC*tb];
  else if(cid==ToF.cid )    return ToF.adcOffset[seg-1+nToF*tb];
  else if(cid==WC.cid )     return WC.adcOffset[seg-1+nWC*tb];
  else   cerr << "[" << funcname << "]: unknown id" << endl;

  return -1.;
}

///////////////////////////////////
double ParamMan::GetAdcGain( int cid, int seg, int tb )
{
  static const std::string funcname = "ParamMan::GetAdcGain";

  if(cid==AC.cid )          return AC.adcGain[seg-1+nAC*tb];
  else if(cid==ToF.cid )    return ToF.adcGain[seg-1+nToF*tb];
  else if(cid==WC.cid )     return WC.adcGain[seg-1+nWC*tb];
  else   cerr << "[" << funcname << "]: unknown id" << endl;

  return -1.;
}

///////////////////////////////////
double ParamMan::npe( int cid, int seg, int tb, double adc )
{
  static const std::string funcname = "ParamMan::npe";

  if(cid==AC.cid )
    return AC.adcGain[seg-1+nAC*tb]*(adc-AC.adcOffset[seg-1+nAC*tb]);
  else if(cid==ToF.cid )
    return ToF.adcGain[seg-1+nToF*tb]*(adc-ToF.adcOffset[seg-1+nToF*tb]);
  else if(cid==WC.cid )
    return WC.adcGain[seg-1+nWC*tb]*(adc-WC.adcOffset[seg-1+nWC*tb]);
  else   cerr << "[" << funcname << "]: unknown id" << endl;

  return -1.;
}

///////////////////////////////////
double ParamMan::time( int cid, int seg, int tb, double tdc )
{
  static const std::string funcname = "ParamMan::npe";

  if(cid==AC.cid )
    return AC.tdcGain[seg-1+nAC*tb]*(tdc-AC.tdcOffset[seg-1+nAC*tb]);
  else if(cid==ToF.cid )
    return ToF.tdcGain[seg-1+nToF*tb]*(tdc-ToF.tdcOffset[seg-1+nToF*tb]);
  else if(cid==WC.cid )
    return WC.tdcGain[seg-1+nWC*tb]*(tdc-WC.tdcOffset[seg-1+nWC*tb]);
  else   cerr << "[" << funcname << "]: unknown id" << endl;

  return -1.;
}

///////////////////////////////////
void ParamMan::SetTdcOffset( int cid, int seg, int tb,
				 double tdcOffset )
{
  static const std::string funcname = "ParamMan::SettdcOffset";
  
  if(cid==AC.cid)    AC.tdcOffset[seg-1+nAC*tb]=tdcOffset;
  else if(cid==ToF.cid)    ToF.tdcOffset[seg-1+nToF*tb]=tdcOffset;
  else if(cid==WC.cid)     WC.tdcOffset[seg-1+nWC*tb]=tdcOffset;
  else   cerr << "[" << funcname << "]: unknown id" << endl;
}
///////////////////////////////////
void ParamMan::SetTdcGain( int cid, int seg, int tb,
			       double tdcGain )
{
  static const std::string funcname = "ParamMan::SettdcGain";
  
  if(cid==AC.cid)          AC.tdcGain[seg-1+nAC*tb]=tdcGain;
  else if(cid==ToF.cid)    ToF.tdcGain[seg-1+nToF*tb]=tdcGain;
  else if(cid==WC.cid)     WC.tdcGain[seg-1+nWC*tb]=tdcGain;
  else   cerr << "[" << funcname << "]: unknown id" << endl;
}
///////////////////////////////////
void ParamMan::SetTimeTune( int cid, int seg, int tb,
				double time )
{
  static const std::string funcname = "ParamMan::SetNpeTune";
  
  if(cid==AC.cid )
    AC.tdcOffset[seg-1+nAC*tb]+=time/AC.tdcGain[seg-1+nAC*tb];
  else if(cid==ToF.cid )
    ToF.tdcOffset[seg-1+nToF*tb]+=time/ToF.tdcGain[seg-1+nToF*tb];
  else if(cid==WC.cid )
    WC.tdcOffset[seg-1+nWC*tb]+=time/WC.tdcGain[seg-1+nWC*tb];
  else   cerr << "[" << funcname << "]: unknown id" << endl;
}

///////////////////////////////////
double ParamMan::GetTdcOffset( int cid, int seg, int tb )
{
  static const std::string funcname = "ParamMan::GetTdcOffset";

  if(cid==AC.cid )          return AC.tdcOffset[seg-1+nAC*tb];
  else if(cid==ToF.cid )    return ToF.tdcOffset[seg-1+nToF*tb];
  else if(cid==WC.cid )     return WC.tdcOffset[seg-1+nWC*tb];
  else   cerr << "[" << funcname << "]: unknown id" << endl;

  return -1.;
}

///////////////////////////////////
double ParamMan::GetTdcGain( int cid, int seg, int tb )
{
  static const std::string funcname = "ParamMan::GetTdcGain";

  if(cid==AC.cid )          return AC.tdcGain[seg-1+nAC*tb];
  else if(cid==ToF.cid )    return ToF.tdcGain[seg-1+nToF*tb];
  else if(cid==WC.cid )     return WC.tdcGain[seg-1+nWC*tb];
  else   cerr << "[" << funcname << "]: unknown id" << endl;

  return -1.;
}


///////////////////////////////////
void ParamMan::WriteToFile(const char* OutputFileName)   //wrinting param file
{
  ofstream fout;
  if( fout.is_open() ) fout.close();
  fout.open(OutputFileName, ios::out|ios::trunc);
  fout.setf(ios_base::fixed);
//   fout.open(name.str().c_str(), std::ios::out|std::ios::trunc);
//   fout.setf(std::ios_base::fixed);
  fout << "#" << endl
       << "#  "  << OutputFileName << endl
       << "#" << endl;
  int at=0; // ADC
  fout << "#ADC#" << endl;
  fout << "# CID SEG AT  TB      Offs        Gain" << endl;
  fout << "# ToF"<< endl;
  for(int k=0; k<2; k++){//k:u,d
    for(int i=0; i<nToF; i++)
      fout << std::setw(4) << CID_ToF
	   << std::setw(4) << i+1
	   << std::setw(4) << at
	   << std::setw(4) << k
	   << std::setw(13) << std::setprecision(6)
	   << ToF.adcOffset[i+nToF*k]
	   << std::setw(11) << std::setprecision(6)
	   << ToF.adcGain[i+nToF*k] << endl;
  }
  fout << "# AC"<< endl;
  for(int k=0; k<2; k++){//k:u,d
    for(int i=0; i<nAC; i++)
      fout << std::setw(4) << CID_AC
	   << std::setw(4) << i+1
	   << std::setw(4) << at
	   << std::setw(4) << k
	   << std::setw(13) << std::setprecision(6)
	   << AC.adcOffset[i+nAC*k]
	   << std::setw(11) << std::setprecision(6)
	   << AC.adcGain[i+nAC*k] << endl;
  }
  fout << "# WC"<< endl;
  for(int k=0; k<2; k++){//k:u,d
    for(int i=0; i<nWC; i++)
      fout << std::setw(4) << CID_WC
	   << std::setw(4) << i+1
	   << std::setw(4) << at
	   << std::setw(4) << k
	   << std::setw(13) << std::setprecision(6)
	   << WC.adcOffset[i+nWC*k]
	   << std::setw(11) << std::setprecision(6)
	   << WC.adcGain[i+nWC*k] << endl;
  }
  at=1; // TDC
  fout << "#TDC#" << endl;
  fout << "# CID SEG AT  TB      Offs        Gain" << endl;
  fout << "# ToF"<< endl;
  for(int k=0; k<2; k++){//k:u,d
    for(int i=0; i<nToF; i++)
      fout << std::setw(4) << CID_ToF
	   << std::setw(4) << i+1
	   << std::setw(4) << at
	   << std::setw(4) << k
	   << std::setw(13) << std::setprecision(6)
	   << ToF.tdcOffset[i+nToF*k]
	   << std::setw(11) << std::setprecision(6)
	   << ToF.tdcGain[i+nToF*k] << endl;
  }
  fout << "# AC"<< endl;
  for(int k=0; k<2; k++){//k:u,d
    for(int i=0; i<nAC; i++)
      fout << std::setw(4) << CID_AC
	   << std::setw(4) << i+1
	   << std::setw(4) << at
	   << std::setw(4) << k
	   << std::setw(13) << std::setprecision(6)
	   << AC.tdcOffset[i+nAC*k]
	   << std::setw(11) << std::setprecision(6)
	   << AC.tdcGain[i+nAC*k] << endl;
  }
  fout << "# WC"<< endl;
  for(int k=0; k<2; k++){//k:u,d
    for(int i=0; i<nWC; i++)
      fout << std::setw(4) << CID_WC
	   << std::setw(4) << i+1
	   << std::setw(4) << at
	   << std::setw(4) << k
	   << std::setw(13) << std::setprecision(6)
	   << WC.tdcOffset[i+nWC*k]
	   << std::setw(11) << std::setprecision(6)
	   << WC.tdcGain[i+nWC*k] << endl;
  }
  if(fout.is_open()) fout.close();
  cout << OutputFileName << " was written"<<endl;
}
