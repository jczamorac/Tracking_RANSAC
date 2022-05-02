#include <sstream>
#include <string>
#include <fstream>
#include <vector>

#include "ransac.h"
#include "mlesac.h"
#include "lmeds.h"
#include "ransac_2Dcircle.h"
using namespace std;





std::vector< std::vector<float> *> pts;

inline void Read_Line(std::istream& in)
{
	while(in) {
		std::string str;
		std::getline(in,str,'\n');
		if(str.length()==0) continue;

		double x0,y0,z0;
		std::stringstream ss;
		ss<<str;

		ss>>x0; ss>>y0; ss>>z0;

		std::vector<float>* p = new std::vector<float>(3);
		pts.push_back(p);
		(*p)[0]=(float)x0;
		(*p)[1]=(float)y0;
    (*p)[2]=(float)z0;

	}
	//std::cout<<"Read Line Done!"<<std::endl;
}


int main(int argc, const char* argv[])
    {


			if(argc<3) {
    		std::cout<<argv[0]<<" infilename outfilename"<<std::endl;
    		return 1;
    	}

    	std::ifstream ifile(argv[1]);
    	Read_Line(ifile);
    	ifile.close();
    	std::ofstream ofile(argv[2]);

      // ***** Part 1: Set input data.
    	//cout<<pts.size()<<endl;
			vector<double> vx;
			vector<double> vy;
			vector<double> vz;
			vector<double> vq;
      const size_t len = pts.size();
    	for(int i=0;i<len;i++){
    		vx.push_back( (*pts[i])[0] );
    		vy.push_back( (*pts[i])[1] );
        vz.push_back( (*pts[i])[2] );
        vq.push_back( 1.0 );
				//std::cout <<vx[i]<<"  "<<vy[i]<< '\n';
    	}// ofile<<(*pts[i])[0]<<"   "<<(*pts[i])[1]<<std::endl;



			// ***** Part 2: Solve.
			/*Ransac *myransac = new Ransac();
			myransac->Init(vx,vy,vz,vq);
			myransac->Solve(10.0,10,500);
			Ransac::AllClusters myClusters = myransac->GetClusters();
			*/
			/*
			Mlesac *mymlesac = new Mlesac();
			mymlesac->Init(vx,vy,vz,vq);
			mymlesac->Solve(10.0,10,500);
			Mlesac::AllClusters myClusters = mymlesac->GetClusters();
			*/
			/*LMedS *mylmeds = new LMedS();
			mylmeds->Init(vx,vy,vz,vq);
			mylmeds->Solve(10.0,10,500);
			LMedS::AllClusters myClusters = mylmeds->GetClusters();
			*/
			Ransac_circle *myransac = new Ransac_circle();
			myransac->Init(vx,vy,vz,vq);
			myransac->Solve(7.0,50,2500);
			Ransac_circle::AllClusters myClusters = myransac->GetClusters();


			// ***** Part 3: Get the clusters
    	int Nclusters = myClusters.size();
			std::cout << "number de clusters "<<Nclusters << '\n';

			for (int i = 0; i < Nclusters; i++) {
      	size_t clustersize = myClusters[i].ClusterSize;
      	std::vector<int> indicesCluster = myClusters[i].ClusterIndex;
      	double costo =  myClusters[i].ClusterStrength;
      	double Chi2 = myClusters[i].ClusterChi2;
      	TVector3 punto1 = myClusters[i].ClusterFitP1;
      	TVector3 punto2 = myClusters[i].ClusterFitP2;
      	std::cout << i<<"  "<<punto1.X()<<" "<<punto1.Y()<<"  "<<punto1.Z()<<"  "<<clustersize<< '\n';
      	std::cout << i<<"  "<<punto2.X()<<" "<<punto2.Y()<<"  "<<punto2.Z()<< '\n';
      	for(int j =0; j<clustersize; j++){
        	ofile<<vx[indicesCluster[j]]<<"  "<<vy[indicesCluster[j]]<<"  "<<vz[indicesCluster[j]]<<"  "<<i<<"  "<<costo<<"  "<<Chi2<<std::endl;
      	}

    }


		// ***** Part 4: Perform the Clustering


    return 0;
}
