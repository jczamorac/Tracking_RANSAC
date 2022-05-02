//////////////////////////////////////////////////////////
// LMedS Class
// J.C. Zamora, cardona@if.usp.br
// University of Sao Paulo, 2020
//////////////////////////////////////////////////////////
#define LMedS_cc

#include "lmeds.h"
#include<iostream>

using namespace std;

ClassImp(LMedS)

LMedS::LMedS() {

  fLMedSMaxIteration = 500;
	fLMedSMinPoints = 30;
	fLMedSThreshold = 15;
  fRandSamplMode = 0;
}

LMedS::~LMedS() {

  delete Rand;
}

void LMedS::Init(vector<double> v1, vector<double> v2, vector<double> v3, vector<double> v4)
{

  Reset();
  vX = v1;
	vY = v2;
	vZ = v3;
	vQ = v4;

	fOriginalCloudSize = vX.size();
	double TotalCharge=0;
	for(unsigned int i=0; i< vQ.size(); i++){
		TotalCharge += vQ[i];
	}
	fTotalCharge = TotalCharge;

  Rand =new TRandom();

  //std::cout << "/* LMedS inicializa  */" << '\n';

}

void LMedS::Reset()
{

	vX.clear();
	vY.clear();
	vZ.clear();
	vQ.clear();
  Vs.Clear();
  Ps.Clear();
  errorsVec.clear();
  //delete Rand;

}


void LMedS::Solve(double dist_thres, double Nminpoints, int Nintera)
{

    fLMedSThreshold = dist_thres;
    fLMedSMinPoints = Nminpoints;
    fLMedSMaxIteration = Nintera;

      //std::cout << "numero de puntos  "<<vX.size()<< '\n';
    std::vector<int> remainIndex;
    for (size_t i = 0; i < vX.size(); i++)
    remainIndex.push_back(i);

  	TVector3 V1, V2;
  	std::vector< int> inliners;
    inliners.clear();
    std::vector< std::pair <double,int> >  IdxMod1;
    std::vector< std::pair <double,int> > IdxMod2;




  	  for(int i=0;i<fLMedSMaxIteration;i++){


        if(remainIndex.size()<fLMedSMinPoints) break;

        std::vector< int> Rsamples = RandSam(remainIndex,fRandSamplMode);  //random sampling
        EstimModel(Rsamples); //estimate the linear model


        //std::vector<int> inlIdxR;
        int nbInliers = 0;

        for (auto j = remainIndex.begin(); j != remainIndex.end(); ++j){

  	      double error = EstimError(*j); //error of each point relative to the model
          error = error*error;

  	      if(error<(fLMedSThreshold*fLMedSThreshold)){
  	        	//inlIdxR.push_back(*j);
              nbInliers++;
              errorsVec.push_back(error);
  	        	}
  		     }

         double med = GetMedian(errorsVec);
         errorsVec.clear();


  	    if(nbInliers>fLMedSMinPoints){
          //getting the best models
          double scale = med/nbInliers;
          IdxMod1.push_back( std::make_pair(scale, Rsamples[0]) );
          IdxMod2.push_back( std::make_pair(scale, Rsamples[1]) );
  			} //if a cluster was found

  		}//for Lmeds interactions


      //sort clusters
      sort(IdxMod1.begin(), IdxMod1.end());
      sort(IdxMod2.begin(), IdxMod2.end());

      remainIndex.clear(); // track remaining points
      for (size_t i = 0; i < vX.size(); i++)
      remainIndex.push_back(i);

      //extract inliers using the models
      for (int i = 0; i < IdxMod1.size(); ++i)
      {
        std::vector<int> ModInx = {IdxMod1[i].second, IdxMod2[i].second};
        EstimModel(ModInx);
        std::vector<int> inlIdxR;
        ModInx.clear();

        if(remainIndex.size()<fLMedSMinPoints) break;

        int counter = 0;

        for (auto j = remainIndex.begin(); j != remainIndex.end(); ++j)
        {
            double error = EstimError(*j);

            if((error*error)<(fLMedSThreshold*fLMedSThreshold))
            {
              inlIdxR.push_back(*j);
              counter++;
            }

        }

        if(counter>fLMedSMinPoints){
          TVector3 v1, v2;
  	      double chi2=Fit3D(inlIdxR,v1,v2);
          SetCluster(inlIdxR, IdxMod1[i].first, chi2,v1,v2);
          v1.Clear();
          v2.Clear();
        }
        std::vector<int> tempRemain;
        std::set_difference(remainIndex.begin(), remainIndex.end(), inlIdxR.begin(), inlIdxR.end(),
        std::inserter(tempRemain, tempRemain.begin()));
        remainIndex = tempRemain;
        inlIdxR.clear();
        tempRemain.clear();
      }

      IdxMod1.clear();
      IdxMod2.clear();
      remainIndex.clear();


}

vector<double> LMedS::GetPDF(const std::vector<int>  samplesIdx){

  size_t pclouds = samplesIdx.size();
  double Tcharge = 0;
  for(int i=0;i<pclouds;i++) Tcharge += vQ[samplesIdx[i]];

  SetAvCharge(Tcharge/pclouds);
  std::vector<double> w;
  if(Tcharge>0)
  for(int i=0;i<pclouds;i++) w.push_back(vQ[samplesIdx[i]]/Tcharge);

  return w;
}

vector<int> LMedS::RandSam(vector<int> indX, Int_t mode)
{
  size_t pclouds = indX.size();
  std::vector<double> Proba = GetPDF(indX);
  int p1,p2;
  double w1,w2;
  vector<int> ranpair;
  ranpair.resize(2);

  if(mode==0){
    //-------Uniform sampling
    p1=(int)(gRandom->Uniform(0,pclouds));

     do{
       p2=(int)(gRandom->Uniform(0,pclouds));
     } while(p2==p1);

     ranpair[0] = indX[p1];
     ranpair[1] = indX[p2];
  }

  if(mode==1){
  //--------Gaussian sampling
    double dist = 0;
    double sigma = 30.0;
    double y = 0;
    double gauss = 0;
    int counter = 0;
    p1=(int)(gRandom->Uniform(0,pclouds));
    TVector3 P1 ={vX[indX[p1]],vY[indX[p1]],vZ[indX[p1]]};
    do{
      p2=(int)(gRandom->Uniform(0,pclouds));
      TVector3 P2 ={vX[indX[p2]],vY[indX[p2]],vZ[indX[p2]]};
      TVector3 dif = P2-P1;
      dist = dif.Mag();
      gauss = 1.0*exp(-1.0*pow(dist/sigma,2.0));
      y = (gRandom->Uniform(0,1));
      counter++;
      if(counter>20 && p2!=p1) break;
      } while(p2==p1 || y>gauss);

      ranpair[0] = indX[p1];
      ranpair[1] = indX[p2];
  }

  if(mode==2){
    //-------Weighted sampling
    bool cond = false;
    int counter = 0;
    p1=(int)(gRandom->Uniform(0,pclouds));
    do{
      counter++;
      if(counter>30 && p2!=p1) break;
      p2=(int)(gRandom->Uniform(0,pclouds));
      cond = false;
      double TwiceAvCharge = 2*GetAvCharge();
      if(Proba.size()==pclouds){
        w2 = gRandom->Uniform(0,TwiceAvCharge);
        if(Proba[p2]>=w2) cond = true;
      }else{
        w2 = 1;
        cond = true;
      }
    } while(p2==p1 || cond==false);

    ranpair[0] = indX[p1];
    ranpair[1] = indX[p2];
  }

  if(mode==3){
    //-------Weighted sampling + Gauss dist.
    bool cond = false;
    double dist = 0;
    double sigma = 30.0;
    double y = 0;
    double gauss = 0;
    int counter = 0;
    p1=(int)(gRandom->Uniform(0,pclouds));
    TVector3 P1 ={vX[indX[p1]],vY[indX[p1]],vZ[indX[p1]]};
    do{
      p2=(int)(gRandom->Uniform(0,pclouds));
      TVector3 P2 ={vX[indX[p2]],vY[indX[p2]],vZ[indX[p2]]};
      TVector3 dif = P2-P1;
      dist = dif.Mag();
      gauss = 1.0*exp(-1.0*pow(dist/sigma,2));
      y = (gRandom->Uniform(0,1));
      counter++;
      if(counter>30 && p2!=p1) break;

      cond = false;
      double TwiceAvCharge = 2*GetAvCharge();
      if(Proba.size()==pclouds){
        w2 = gRandom->Uniform(0,TwiceAvCharge);
        if(Proba[p2]>=w2) cond = true;
        }else{
        w2 = 1;
        cond = true;
      }

    } while(p2==p1 || cond==false || y>gauss);

    ranpair[0] = indX[p1];
    ranpair[1] = indX[p2];
  }

  return ranpair;

}


void LMedS::EstimModel(const std::vector<int>  samplesIdx)
{

  //line from two points
  TVector3 Po1 = {vX[samplesIdx[0]], vY[samplesIdx[0]], vZ[samplesIdx[0]]};
  TVector3 Po2 = {vX[samplesIdx[1]], vY[samplesIdx[1]], vZ[samplesIdx[1]]};

  Vs = Po2 - Po1;
  Ps = Po1;

}

double LMedS::EstimError(int i)
{
    //distance point to line
    TVector3 newPoint = {vX[i], vY[i], vZ[i]};
    TVector3 vec = Ps - newPoint;
    TVector3 nD = Vs.Cross(vec);

	  double dist = nD.Mag()/Vs.Mag();


    return  dist;
}


void LMedS::SetCluster(const std::vector<int> samplesIdx, const double cost, const double Chi2, TVector3 CP1, TVector3 CP2)
{

    Cluster cstr;
    cstr.ClusterIndex = samplesIdx;
    cstr.ClusterSize = samplesIdx.size();
    cstr.ClusterStrength = cost;
    cstr.ClusterChi2 = Chi2;
    cstr.ClusterFitP1 = CP1;
    cstr.ClusterFitP2 = CP2;
    cluster_vector.push_back(cstr);
}

double LMedS::GetMedian(std::vector<double> errvec)
{
  size_t vsize = errvec.size();

  if (vsize == 0)
  {
    return 0;
  }
  else
  {
    sort(errvec.begin(), errvec.end());
    if (vsize % 2 == 0)
    {
      return (errvec[vsize / 2 - 1] + errvec[vsize / 2]) / 2;
    }
    else
    {
      return errvec[vsize / 2];
    }
  }

}

double LMedS::Fit3D( vector<int> inliners, TVector3& V1, TVector3& V2)
{
    //------3D Line Regression
    //----- adapted from: http://fr.scribd.com/doc/31477970/Regressions-et-trajectoires-3D
    int R, C;
    double Q;
    double Xm,Ym,Zm;
    double Xh,Yh,Zh;
    double a,b;
    double Sxx,Sxy,Syy,Sxz,Szz,Syz;
    double theta;
    double K11,K22,K12,K10,K01,K00;
    double c0,c1,c2;
    double p,q,r,dm2;
    double rho,phi;

    Q=Xm=Ym=Zm=0.;
		double total_charge=0;
    Sxx=Syy=Szz=Sxy=Sxz=Syz=0.;

    for (auto i : inliners)
    {
        Q+=vQ[i]/10.;
        Xm+=vX[i]*vQ[i]/10.;
        Ym+=vY[i]*vQ[i]/10.;
        Zm+=vZ[i]*vQ[i]/10.;
        Sxx+=vX[i]*vX[i]*vQ[i]/10.;
        Syy+=vY[i]*vY[i]*vQ[i]/10.;
        Szz+=vZ[i]*vZ[i]*vQ[i]/10.;
        Sxy+=vX[i]*vY[i]*vQ[i]/10.;
        Sxz+=vX[i]*vZ[i]*vQ[i]/10.;
        Syz+=vY[i]*vZ[i]*vQ[i]/10.;
    }
    //vTrackCharge.push_back(total_charge);

    Xm/=Q;
    Ym/=Q;
    Zm/=Q;
    Sxx/=Q;
    Syy/=Q;
    Szz/=Q;
    Sxy/=Q;
    Sxz/=Q;
    Syz/=Q;
    Sxx-=(Xm*Xm);
    Syy-=(Ym*Ym);
    Szz-=(Zm*Zm);
    Sxy-=(Xm*Ym);
    Sxz-=(Xm*Zm);
    Syz-=(Ym*Zm);

    theta=0.5*atan((2.*Sxy)/(Sxx-Syy));

    K11=(Syy+Szz)*pow(cos(theta),2)+(Sxx+Szz)*pow(sin(theta),2)-2.*Sxy*cos(theta)*sin(theta);
    K22=(Syy+Szz)*pow(sin(theta),2)+(Sxx+Szz)*pow(cos(theta),2)+2.*Sxy*cos(theta)*sin(theta);
    K12=-Sxy*(pow(cos(theta),2)-pow(sin(theta),2))+(Sxx-Syy)*cos(theta)*sin(theta);
    K10=Sxz*cos(theta)+Syz*sin(theta);
    K01=-Sxz*sin(theta)+Syz*cos(theta);
    K00=Sxx+Syy;

    c2=-K00-K11-K22;
    c1=K00*K11+K00*K22+K11*K22-K01*K01-K10*K10;
    c0=K01*K01*K11+K10*K10*K22-K00*K11*K22;


    p=c1-pow(c2,2)/3.;
    q=2.*pow(c2,3)/27.-c1*c2/3.+c0;
    r=pow(q/2.,2)+pow(p,3)/27.;


    if(r>0) dm2=-c2/3.+pow(-q/2.+sqrt(r),1./3.)+pow(-q/2.-sqrt(r),1./3.);
    if(r<0)
    {
        rho=sqrt(-pow(p,3)/27.);
        phi=acos(-q/(2.*rho));
        dm2=min(-c2/3.+2.*pow(rho,1./3.)*cos(phi/3.),min(-c2/3.+2.*pow(rho,1./3.)*cos((phi+2.*TMath::Pi())/3.),-c2/3.+2.*pow(rho,1./3.)*cos((phi+4.*TMath::Pi())/3.)));
    }

    a=-K10*cos(theta)/(K11-dm2)+K01*sin(theta)/(K22-dm2);
    b=-K10*sin(theta)/(K11-dm2)-K01*cos(theta)/(K22-dm2);

    Xh=((1.+b*b)*Xm-a*b*Ym+a*Zm)/(1.+a*a+b*b);
    Yh=((1.+a*a)*Ym-a*b*Xm+b*Zm)/(1.+a*a+b*b);
    Zh=((a*a+b*b)*Zm+a*Xm+b*Ym)/(1.+a*a+b*b);

    V1.SetXYZ(Xm,Ym,Zm);
    V2.SetXYZ(Xh,Yh,Zh);

    return(fabs(dm2/Q));
}
