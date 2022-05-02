//////////////////////////////////////////////////////////
// RANSAC 2D circle Class
// J.C. Zamora, cardona@if.usp.br
// University of Sao Paulo, 2020
//////////////////////////////////////////////////////////
#define RANSAC2Dcir_cc

#include "ransac_2Dcircle.h"
#include<iostream>

using namespace std;

ClassImp(Ransac_circle)

Ransac_circle::Ransac_circle() {

  fRANSACMaxIteration = 500;
	fRANSACMinPoints = 30;
	fRANSACThreshold = 15;
  fRandSamplMode = 1;
}

Ransac_circle::~Ransac_circle() {

  delete Rand;
}

void Ransac_circle::Init(vector<double> v1, vector<double> v2, vector<double> v3, vector<double> v4)
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

  //std::cout << "/* Ransac_circle inicializa  */" << '\n';

}

void Ransac_circle::Reset()
{

	vX.clear();
	vY.clear();
	vZ.clear();
	vQ.clear();
  Mcir.Clear();
  //delete Rand;

}


void Ransac_circle::Solve(double dist_thres, double Nminpoints, int Nintera)
{

    fRANSACThreshold = dist_thres;
    fRANSACMinPoints = Nminpoints;
    fRANSACMaxIteration = Nintera;

  //std::cout << "numero de puntos  "<<vX.size()<< '\n';
  std::vector<int> remainIndex;
  for (size_t i = 0; i < vX.size(); i++)
    remainIndex.push_back(i);

  TVector3 V1, V2;
  std::vector< int> inliners;
  inliners.clear();
  std::vector< std::pair <double,int> >  IdxMod1;
  std::vector< std::pair <double,int> > IdxMod2;
  std::vector< std::pair <double,int> > IdxMod3;




  for(int i=0;i<fRANSACMaxIteration;i++){


    if(remainIndex.size()<fRANSACMinPoints) break;

    std::vector< int> Rsamples = RandSam(remainIndex,fRandSamplMode);  //random sampling
    EstimModel(Rsamples); //estimate the linear model


    //std::vector<int> inlIdxR;
    int nbInliers = 0;
    double weight = 0;

    for (auto j = remainIndex.begin(); j != remainIndex.end(); ++j){

      double error = EstimError(*j); //error of each point relative to the model
      error = error*error;

      if(error<(fRANSACThreshold*fRANSACThreshold)){
          //inlIdxR.push_back(*j);
          nbInliers++;
          weight +=error;
          }
       }


    if(nbInliers>fRANSACMinPoints){
      //getting the best models
      double scale = weight/nbInliers;
      IdxMod1.push_back( std::make_pair(scale, Rsamples[0]) );
      IdxMod2.push_back( std::make_pair(scale, Rsamples[1]) );
      IdxMod3.push_back( std::make_pair(scale, Rsamples[2]) );
    } //if a cluster was found

  }//for RANSAC interactions


  //sort clusters
  sort(IdxMod1.begin(), IdxMod1.end());
  sort(IdxMod2.begin(), IdxMod2.end());
  sort(IdxMod3.begin(), IdxMod3.end());

  remainIndex.clear(); // track remaining points
  for (size_t i = 0; i < vX.size(); i++)
    remainIndex.push_back(i);

  //extract inliers using the models
  for (int i = 0; i < IdxMod1.size(); ++i)
  {
    std::vector<int> ModInx = {IdxMod1[i].second, IdxMod2[i].second, IdxMod3[i].second};
    EstimModel(ModInx);
    std::vector<int> inlIdxR;
    //inlIdxR.clear();
    ModInx.clear();

    if(remainIndex.size()<fRANSACMinPoints) break;

    int counter = 0;

    for (auto j = remainIndex.begin(); j != remainIndex.end(); ++j)
    {
        double error = EstimError(*j);


        if((error*error)<(fRANSACThreshold*fRANSACThreshold))
        {
          inlIdxR.push_back(*j);
          counter++;
        }

    }

    if(counter>fRANSACMinPoints){
      TVector3 v1, v2;
      double chi2=Fit2DCircle(inlIdxR,v2);
      //double chi2=1;
      SetCluster(inlIdxR, IdxMod1[i].first, chi2,Mcir,v2);
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
  IdxMod3.clear();
  remainIndex.clear();


}

vector<double> Ransac_circle::GetPDF(const std::vector<int>  samplesIdx){

  size_t pclouds = samplesIdx.size();
  double Tcharge = 0;
  for(int i=0;i<pclouds;i++) Tcharge += vQ[samplesIdx[i]];

  SetAvCharge(Tcharge/pclouds);
  std::vector<double> w;
  if(Tcharge>0)
  for(int i=0;i<pclouds;i++) w.push_back(vQ[samplesIdx[i]]/Tcharge);

  return w;
}

vector<int> Ransac_circle::RandSam(vector<int> indX, Int_t mode)
{

  size_t pclouds = indX.size();
  std::vector<double> Proba = GetPDF(indX);
  int p1,p2, p3;
  double w1,w2;
  vector<int> ranpair;
  ranpair.resize(3);

  if(mode==0){
    //-------Uniform sampling
    p1=(int)(gRandom->Uniform(0,pclouds));

	   do{
       p2=(int)(gRandom->Uniform(0,pclouds));
     } while(p2==p1);

     do{
       p3=(int)(gRandom->Uniform(0,pclouds));
     } while(p3==p1 || p3==p2);

     ranpair[0] = indX[p1];
     ranpair[1] = indX[p2];
     ranpair[2] = indX[p3];
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

      counter = 0;
    do{
      p3=(int)(gRandom->Uniform(0,pclouds));
      TVector3 P3 ={vX[indX[p3]],vY[indX[p3]],vZ[indX[p3]]};
      TVector3 dif = P3-P1;
      dist = dif.Mag();
      gauss = 1.0*exp(-1.0*pow(dist/sigma,2.0));
      y = (gRandom->Uniform(0,1));
      counter++;
      if(counter>20 && p3!=p1) break;
    } while((p2==p1 || p3==p2 ) || y>gauss);

      ranpair[0] = indX[p1];
      ranpair[1] = indX[p2];
      ranpair[2] = indX[p3];
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

    cond = false;
    counter = 0;
    do{
      counter++;
      if(counter>30 && p3!=p1) break;
      p3=(int)(gRandom->Uniform(0,pclouds));
      cond = false;
      double TwiceAvCharge = 2*GetAvCharge();
      if(Proba.size()==pclouds){
        w2 = gRandom->Uniform(0,TwiceAvCharge);
        if(Proba[p3]>=w2) cond = true;
      }else{
        w2 = 1;
        cond = true;
      }
    } while(p3==p1 || p2==p3 || cond==false);

    ranpair[0] = indX[p1];
    ranpair[1] = indX[p2];
    ranpair[2] = indX[p3];
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

    cond = false;
    counter = 0;
    do{
      p3=(int)(gRandom->Uniform(0,pclouds));
      TVector3 P3 ={vX[indX[p3]],vY[indX[p3]],vZ[indX[p3]]};
      TVector3 dif = P3-P1;
      dist = dif.Mag();
      gauss = 1.0*exp(-1.0*pow(dist/sigma,2));
      y = (gRandom->Uniform(0,1));
      counter++;
      if(counter>30 && p3!=p1 &&p3!=p2) break;

      cond = false;
      double TwiceAvCharge = 2*GetAvCharge();
      if(Proba.size()==pclouds){
        w2 = gRandom->Uniform(0,TwiceAvCharge);
        if(Proba[p3]>=w2) cond = true;
        }else{
        w2 = 1;
        cond = true;
      }

    } while(p3==p1 || p3==p2 || cond==false || y>gauss);


    ranpair[0] = indX[p1];
    ranpair[1] = indX[p2];
    ranpair[2] = indX[p3];
  }

  return ranpair;

}


void Ransac_circle::EstimModel(const std::vector<int>  samplesIdx)
{

  //circle from three points
  //for the moment forget about the Z coordinate
  TVector3 Po1 = {vX[samplesIdx[0]], vY[samplesIdx[0]], 0};
  TVector3 Po2 = {vX[samplesIdx[1]], vY[samplesIdx[1]], 0};
  TVector3 Po3 = {vX[samplesIdx[2]], vY[samplesIdx[2]], 0};

  double a = 2 * (Po2.X() - Po1.X());
  double b = 2 * (Po2.Y() - Po1.Y());
  double c = Po2.X()*Po2.X() + Po2.Y()*Po2.Y() - Po1.X()*Po1.X() - Po1.Y()*Po1.Y();
	double d = 2 * (Po3.X() - Po2.X());
  double e = 2 * (Po3.Y() - Po2.Y());
  double f = Po3.X()*Po3.X() + Po3.Y()*Po3.Y() - Po2.X()*Po2.X() - Po2.Y()*Po2.Y();

	double xc = (b*f - e*c)/(b*d - e*a);
	double yc = (d*c - a*f)/(b*d - e*a);
	double r = sqrt((xc - Po1.X())*(xc - Po1.X()) + (yc - Po1.Y())*(yc - Po1.Y()));

  Mcir.SetXYZ(xc,yc,r);


}

double Ransac_circle::EstimError(int i)
{
    //distance between point to center
    TVector3 newPoint = {vX[i], vY[i], vZ[i]};
    TVector3 vec = Mcir - newPoint;


	  double dist = fabs(sqrt(vec.X()*vec.X() + vec.Y()*vec.Y()) - Mcir.Z());


    return  dist;
}


void Ransac_circle::SetCluster(const std::vector<int> samplesIdx, const double cost, const double Chi2, TVector3 CP1, TVector3 CP2)
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


double Ransac_circle::Fit2DCircle( vector<int> inliners, TVector3& V1)
{

    //  2D circle fit, due to Taubin, based on the journal article
    // G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
    //             Space Curves Defined By Implicit Equations, With
    //             Applications To Edge And Range Image Segmentation",
    //             IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
    //----- adapted from: https://people.cas.uab.edu/~mosya/cl/CPPcircle.html

  int i,iter,IterMAX=99;
  int Niliers = 0;
  double Xi,Yi,Zi;
  double Xm,Ym,Zm;
  double Mz,Mxy,Mxx,Myy,Mxz,Myz,Mzz,Cov_xy,Var_z;
  double A0,A1,A2,A22,A3,A33;
  double Dy,xnew,x,ynew,y;
  double DET,Xcenter,Ycenter;
  double Q;
  double a,b, r;

  for (auto i : inliners)
  {
      Q+=vQ[i]/10.;
      Xm+=vX[i]*vQ[i]/10.;
      Ym+=vY[i]*vQ[i]/10.;
      Zm+=vZ[i]*vQ[i]/10.;
      Niliers++;

  }

  Xm/=Q;
  Ym/=Q;
  Zm/=Q;


  Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;

  for (auto i : inliners)
  {
      Xi = vX[i] - Xm;   //  centered x-coordinates
      Yi = vY[i] - Ym;   //  centered y-coordinates
      Zi = Xi*Xi + Yi*Yi;

      Mxy += Xi*Yi;
      Mxx += Xi*Xi;
      Myy += Yi*Yi;
      Mxz += Xi*Zi;
      Myz += Yi*Zi;
      Mzz += Zi*Zi;
  }
  Mxx /= Niliers;
  Myy /= Niliers;
  Mxy /= Niliers;
  Mxz /= Niliers;
  Myz /= Niliers;
  Mzz /= Niliers;


//      computing coefficients of the characteristic polynomial

  Mz = Mxx + Myy;
  Cov_xy = Mxx*Myy - Mxy*Mxy;
  Var_z = Mzz - Mz*Mz;
  A3 = 4.0*Mz;
  A2 = -3.0*Mz*Mz - Mzz;
  A1 = Var_z*Mz + 4.0*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
  A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
  A22 = A2 + A2;
  A33 = A3 + A3 + A3;

//    finding the root of the characteristic polynomial
//    using Newton's method starting at x=0
//     (it is guaranteed to converge to the right root)

  for (x=0.,y=A0,iter=0; iter<IterMAX; iter++)  // usually, 4-6 iterations are enough
  {
        Dy = A1 + x*(A22 + A33*x);
      xnew = x - y/Dy;
      if ((xnew == x)||(!isfinite(xnew))) break;
      ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
      if (abs(ynew)>=abs(y))  break;
      x = xnew;  y = ynew;
  }

//       computing paramters of the fitting circle

  DET = x*x - x*Mz + Cov_xy;
  Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2.0;
  Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2.0;




//       assembling the output


  a = Xcenter + Xm;
  b = Ycenter + Ym;
  r = sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz);
  V1.SetXYZ(a,b,r);

  return fabs(Mz - r*r);

}
