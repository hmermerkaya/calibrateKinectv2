/*
 * sample.cpp
 *
 *  Created on: Jul 10, 2018
 *      Author: hamit
 */

//#include <TRint.h>
//#include <TEveManager.h>
#include <TROOT.h>
#include <TMath.h>
#include <Math/GenVector/Transform3D.h>
#include <Math/GenVector/Translation3D.h>
#include <Math/GenVector/EulerAngles.h>
#include <TVector3.h>
#include <TRotation.h>
#include <TMinuit.h>
#include <TVirtualFitter.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <highgui.h>
#include <algorithm>
#include <TMatrixD.h>

#include "k2g.h"
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/registration/transforms.h>

#include <chrono>

// extra headers for writing out ply file
#include <pcl/console/print.h>
#include <pcl/console/parse.h>
#include <pcl/console/time.h>
#include <pcl/io/ply_io.h>
#include <pcl/io/png_io.h>

using namespace std;

typedef ROOT::Math::Translation3D Translation;
typedef ROOT::Math::Transform3D Tranformation;
typedef ROOT::Math::Rotation3D Rotation;
typedef ROOT::Math::Transform3D::Vector Vector3;


cv::Mat euler2rot( float thetax, float thetay, float thetaz)
{
  cv::Mat rotationMatrix(3,3,CV_32F);

  float x = thetax;
  float y = thetay;
  float z = thetaz;

  // Assuming the angles are in radians.
  float ch = cos(z);
  float sh = sin(z);
  float ca = cos(y);
  float sa = sin(y);
  float cb = cos(x);
  float sb = sin(x);

  float m00, m01, m02, m10, m11, m12, m20, m21, m22;

  m00 = ch * ca;
  m01 = sh*sb - ch*sa*cb;
  m02 = ch*sa*sb + sh*cb;
  m10 = sa;
  m11 = ca*cb;
  m12 = -ca*sb;
  m20 = -sh*ca;
  m21 = sh*sa*cb + ch*sb;
  m22 = -sh*sa*sb + ch*cb;

  rotationMatrix.at<float>(0,0) = m00;
  rotationMatrix.at<float>(0,1) = m01;
  rotationMatrix.at<float>(0,2) = m02;
  rotationMatrix.at<float>(1,0) = m10;
  rotationMatrix.at<float>(1,1) = m11;
  rotationMatrix.at<float>(1,2) = m12;
  rotationMatrix.at<float>(2,0) = m20;
  rotationMatrix.at<float>(2,1) = m21;
  rotationMatrix.at<float>(2,2) = m22;

  return rotationMatrix;
}
TMatrixD eulerAnglesToRotationMatrix(double thetax, double thetay, double thetaz) {

    double elemenstR_x[9]= { 1,         0,                  0                   ,
                            0,         cos(thetax),    -sin(thetax),
                            0,         sin(thetax),     cos(thetax)
    					   };

    TMatrixD R_x(3,3);
    R_x.SetMatrixArray(elemenstR_x);

    double elemenstR_y[9]={ cos(thetay),    0,      sin(thetay)  ,
    						0,                     1,                        0 ,
							-sin(thetay),   0,      cos(thetay)
    					};
    TMatrixD R_y(3,3);
    R_y.SetMatrixArray(elemenstR_y);

    double elemenstR_z[9]={ cos(thetaz),    -sin(thetaz),    0,
                    		sin(thetaz),    cos(thetaz),     0,
							0,                     0,                      1
     };
    TMatrixD R_z(3,3);
    R_z.SetMatrixArray(elemenstR_z);

    return R_z*R_y*R_x;
}
TMatrixD euler2TMatrixD( double thetax, double thetay, double thetaz)
{
  TMatrixD rotationMatrix(3,3);

  double x = thetax;
  double y = thetay;
  double z = thetaz;

  // Assuming the angles are in radians.
  double ch = cos(z);
  double sh = sin(z);
  double ca = cos(y);
  double sa = sin(y);
  double cb = cos(x);
  double sb = sin(x);

 // float m00, m01, m02, m10, m11, m12, m20, m21, m22;

 double m[9]={ch * ca, sh*sb - ch*sa*cb, ch*sa*sb + sh*cb, sa, ca*cb, -ca*sb, -sh*ca, sh*sa*cb + ch*sb, -sh*sa*sb + ch*cb};

//  m00 = ch * ca;
//  m01 = sh*sb - ch*sa*cb;
//  m02 = ch*sa*sb + sh*cb;
//  m10 = sa;
//  m11 = ca*cb;
//  m12 = -ca*sb;
//  m20 = -sh*ca;
//  m21 = sh*sa*cb + ch*sb;
//  m22 = -sh*sa*sb + ch*cb;
//
//  rotationMatrix(0,0) = m00;
//  rotationMatrix(0,1) = m01;
//  rotationMatrix(0,2) = m02;
//  rotationMatrix(1,0) = m10;
//  rotationMatrix(1,1) = m11;
//  rotationMatrix(1,2) = m12;
//  rotationMatrix(2,0) = m20;
//  rotationMatrix(2,1) = m21;
//  rotationMatrix(2,2) = m22;

  rotationMatrix.SetMatrixArray(m);

  return rotationMatrix;
}

std::vector<std::pair<cv::Point3d,cv::Point>> pointsUVvectorIRAndHD;
std::pair<cv::Mat, cv::Mat> pointsCoordsMatIRAndHD;
std::pair<TMatrixD, TMatrixD>  pointsCoordsTMatrixIRAndHD;
//cv::Mat UVs;

std::vector<cv::Point3d>  pointvectorIR;
void loadPointsFromFiles(const std::string &, std::vector<std::pair<cv::Point3d,cv::Point>> &);




void printFinalResults(double par[]){
//	  for (int i = 0; i <9; ++i) {
//
//
//		     std::cout<<par[i]<<" ";
//
//	    }

       // par[0]=0; par[1]=0; par[2]=0; par[3]=0; par[4]=0; par[5]=0; par[6]=1081.37; par[7]=959.5; par[8]=539.5;
	    TMatrixD rot = euler2TMatrixD(par[0], par[1], par[2]);
	    TMatrixD transl= TMatrixD(1,3);
	    transl[0][0]=par[3];  transl[0][1]=par[4];  transl[0][2]=par[5];
		//pointsCoordsTMatrixIRAndHD.first.Print();
	//rot.Print();
		TMatrixD rotated=  rot*pointsCoordsTMatrixIRAndHD.first;
		rotated.T();

		//rotated.Print();
		//trans.T().Print();

	    std::vector<double> uvsVec;
		double fx_hd=1081.37,fy_hd=1081.37,cx_hd=959.5,cy_hd=539.5;

		for(int i = 0; i < rotated.GetNrows(); i++) {

		//		double x= fx_hd * (rotated(i,0) + transl(0,0)) / (rotated(i,2) + transl(0,2))  + cx_hd;
		//		double y= fy_hd * (rotated(i,1) + transl(0,1)) / (rotated(i,2) + transl(0,2))  + cy_hd;
				double x= par[6] * (rotated(i,0) + transl(0,0)) / (rotated(i,2) + transl(0,2))  + par[7];
				double y= par[6] * (rotated(i,1) + transl(0,1)) / (rotated(i,2) + transl(0,2))  + par[8];
				uvsVec.push_back(x);uvsVec.push_back(y);
	         //  std::cout<<"x and y: "<<x <<" "<<y<<std::endl;

		}


		TMatrixD uvsFound= TMatrixD(rotated.GetNrows(),2);
		uvsFound.SetMatrixArray(uvsVec.data());
		uvsFound.Print();
		TMatrixD uvs= TMatrixD(pointsCoordsTMatrixIRAndHD.second);
		//uvs.T();
		uvs.Print();

		TMatrixD deltaUVs=uvsFound-uvs;
		deltaUVs*=1/sqrt(1);
		double chi2 = deltaUVs.Sqr().Sum();
		std::cout<<"chi2= "<<chi2<<std::endl;


}


void fcnchi2_2(int &, double *, double & sum, double * par, int ) {

  pointsCoordsTMatrixIRAndHD.first.ResizeTo(0,0);
  pointsCoordsTMatrixIRAndHD.first.ResizeTo(0,0);
	std::vector<double> vecOfCoordinatesIR;
	std::vector<double> vecOfCoordinatesHD;
     int i=0;
	for (auto &it:pointsUVvectorIRAndHD) {
		double c = it.first.x;
		double r = it.first.y;
		double depth_val=it.first.z/1000.;

		if (!isnan(depth_val) && depth_val >= 0.001)
		  {
			double x_ir = (c + 0.5 - par[10]) * (1./par[9]) * depth_val;
			double y_ir = (r + 0.5 - par[11]) * (1./par[9]) * depth_val;
			double z_ir = depth_val;
			//std::cout<<"Coords "<<c<<" "<<r<<" "<<x_ir<<" "<<y_ir<<" "<<z_ir<<std::endl;
//			double x_hd = (c + 0.5 - cx_ir) * fx_ir * depth_val;
//			double y_ir = (r + 0.5 - cy_ir) * fy_ir * depth_val;
//			double z_ir = depth_val;
			//cv::Mat tmpmatir(1,3,CV_32F);

			vecOfCoordinatesIR.push_back(x_ir);vecOfCoordinatesIR.push_back(y_ir);vecOfCoordinatesIR.push_back(z_ir);

			vecOfCoordinatesHD.push_back(it.second.x);vecOfCoordinatesHD.push_back(it.second.y);
			i++;

		  }


	}
	pointsCoordsTMatrixIRAndHD.first.ResizeTo(i,3);
	pointsCoordsTMatrixIRAndHD.first.SetMatrixArray(vecOfCoordinatesIR.data());
	pointsCoordsTMatrixIRAndHD.first.T();
	pointsCoordsTMatrixIRAndHD.second.ResizeTo(i, 2 );
	pointsCoordsTMatrixIRAndHD.second.SetMatrixArray(vecOfCoordinatesHD.data());


	TMatrixD rot = euler2TMatrixD(par[0], par[1], par[2]);
    TMatrixD transl= TMatrixD(1,3);
    transl[0][0]=par[3];  transl[0][1]=par[4];  transl[0][2]=par[5];
//rot.Print();
	TMatrixD rotated =  rot*pointsCoordsTMatrixIRAndHD.first;
	rotated.T();

	//rotated.Print();
	//trans.T().Print();

    std::vector<double> uvsVec;
	double fx_hd=1081.37,fy_hd=1081.37,cx_hd=959.5,cy_hd=539.5;

	for(int i = 0; i < rotated.GetNrows(); i++) {

	//		double x= fx_hd * (rotated(i,0) + transl(0,0)) / (rotated(i,2) + transl(0,2))  + cx_hd;
	//		double y= fy_hd * (rotated(i,1) + transl(0,1)) / (rotated(i,2) + transl(0,2))  + cy_hd;
			double x= par[6] * (rotated(i,0) + transl(0,0)) / (rotated(i,2) + transl(0,2))  + par[7];
			double y= par[6] * (rotated(i,1) + transl(0,1)) / (rotated(i,2) + transl(0,2))  + par[8];
			uvsVec.push_back(x);uvsVec.push_back(y);
         //  std::cout<<"x and y: "<<x <<" "<<y<<std::endl;

	}


	TMatrixD uvsFound= TMatrixD(rotated.GetNrows(),2);
	uvsFound.SetMatrixArray(uvsVec.data());
//	uvsFound.Print();
	TMatrixD uvs= TMatrixD(pointsCoordsTMatrixIRAndHD.second);
	//uvs.T();
//	uvs.Print();

	TMatrixF deltaUVs=uvsFound-uvs;
	deltaUVs*=1/sqrt(1);
	double chi2 = deltaUVs.Sqr().Sum();

//    std::vector<double> arruvsFound, arruvs;
//    uvsFound.GetMatrix2Array(arruvsFound.data());
//    uvs.GetMatrix2Array(arruvs.data());
//    double tot = 0;
//    for (int i =0; i<arruvsFound.size();i+=2){
//
//    	tot = tot + fabs(pow(arruvsFound[i],2) + pow(arruvsFound[i+1],2)  -  pow(arruvs[i],2) - pow(arruvs[i+1],2) );
//
//    }










   sum=chi2;





}
void fcnchi2(int &, double *, double & sum, double * par, int ) {
//	ROOT::Math::EulerAngles eulerangles(par[0], par[1], par[2]);
//
//	Tranformation trans(Rotation(eulerangles), Translation(par[3], par[4], par[5]));
//	Vector3 v3(1,3,4);


	//TMatrixF rot = eulerAnglesToRotationMatrix(par[0], par[1], par[2]);
	TMatrixD rot = euler2TMatrixD(par[0], par[1], par[2]);

    TMatrixD transl= TMatrixD(1,3);
    transl[0][0]=par[3];  transl[0][1]=par[4];  transl[0][2]=par[5];
	//pointsCoordsTMatrixIRAndHD.first.Print();
//rot.Print();
	TMatrixD rotated=  rot*pointsCoordsTMatrixIRAndHD.first;
	rotated.T();

	//rotated.Print();
	//trans.T().Print();

    std::vector<double> uvsVec;
	double fx_hd=1081.37,fy_hd=1081.37,cx_hd=959.5,cy_hd=539.5;

	for(int i = 0; i < rotated.GetNrows(); i++) {

	//		double x= fx_hd * (rotated(i,0) + transl(0,0)) / (rotated(i,2) + transl(0,2))  + cx_hd;
	//		double y= fy_hd * (rotated(i,1) + transl(0,1)) / (rotated(i,2) + transl(0,2))  + cy_hd;
			double x= par[6] * (rotated(i,0) + transl(0,0)) / (rotated(i,2) + transl(0,2))  + par[7];
			double y= par[6] * (rotated(i,1) + transl(0,1)) / (rotated(i,2) + transl(0,2))  + par[8];
			uvsVec.push_back(x);uvsVec.push_back(y);
         //  std::cout<<"x and y: "<<x <<" "<<y<<std::endl;

	}


	TMatrixD uvsFound= TMatrixD(rotated.GetNrows(),2);
	uvsFound.SetMatrixArray(uvsVec.data());
//	uvsFound.Print();
	TMatrixD uvs= TMatrixD(pointsCoordsTMatrixIRAndHD.second);
	//uvs.T();
//	uvs.Print();

	TMatrixD deltaUVs=uvsFound-uvs;
	deltaUVs*=1/sqrt(1);
	double chi2 = deltaUVs.Sqr().Sum();






//	double U_depth= fx_hd*pt.x/pt.z + cx_hd;
//	double V_depth= fy_hd*pt.y/pt.z + cy_hd;





   sum=chi2;





}
void minimize() {

		//loadPointsFromFiles("./", pointvectorIRAndHD);


	   TVirtualFitter *min = TVirtualFitter::Fitter(0,9);
	   // min->SetObjectFit(gr);
	    min->SetFCN(fcnchi2);

	    Double_t arglist[10];
	    arglist[0] = 1;
	    int ierflg = 1982;
	    min->ExecuteCommand("SET PRINT",arglist,ierflg);

	    double pStart[9] = { 0, 0.0, 0,
	    					 0,  0.0, -0.0,
							 1081.37, 960., 540.
						//	 364.963, 512./2, 424./2.
							};

	    double pStep[9]= {  0.0001, 0.0001, 0.0001,
                			0.001, 0.001, 0.001,
							0.05000, 0.05000, 0.5000
					//		0.5000, 0.5000, 0.5000
						 };


	    min->SetParameter(0, "par0", pStart[0], pStep[0], -TMath::Pi(), TMath::Pi() );
	    min->SetParameter(1, "par1", pStart[1], pStep[1], -TMath::Pi(), TMath::Pi() );
	    min->SetParameter(2, "par2", pStart[2], pStep[2], -TMath::Pi(), TMath::Pi() );
	    min->SetParameter(3, "par3", pStart[3], pStep[3], -2,2);
	    min->SetParameter(4, "par4", pStart[4], pStep[4], -2,2);
	    min->SetParameter(5, "par5", pStart[5], pStep[5], -2,2);
	    min->SetParameter(6, "par6", pStart[6], pStep[6],  1000, 1200);
	 	min->SetParameter(7, "par7", pStart[7], pStep[7],  850, 1000);
	 	min->SetParameter(8, "par8", pStart[8], pStep[8],  450, 600);
//	 	min->SetParameter(9, "par9", pStart[9], pStep[9],  320, 380);
//		min->SetParameter(10, "par10", pStart[10], pStep[10],  220, 290);
//		min->SetParameter(11, "par11", pStart[11], pStep[11],  180, 250);



	 	 min->FixParameter(6);
	 	 min->FixParameter(7);
	 	 min->FixParameter(8);
//	 	 min->FixParameter(9);
//	 	 min->FixParameter(10);
//	 	 min->FixParameter(11);

	   // min->SetParameter(2,"b2",pStart[2],0.01,0,0);

	    arglist[0] = 5000000000000; // number of function calls
	    arglist[1] = 0.001; // tolerance
//	    min->ExecuteCommand("MIGRAD",arglist,2);
//	    min->ExecuteCommand("HESSE", arglist,2);
	    min->ExecuteCommand("MINOS",arglist,2);
	    int nvpar,nparx;
	    double amin,edm, errdef;
	    min->GetStats(amin,edm,errdef,nvpar,nparx);
	   //min->PrintResults(1,amin);
	    double errFit[9], parFit[9];
	    std::cout<<"Parameters: \n";
	    for (int i = 0; i <9; ++i) {
	      errFit[i] = min->GetParError(i);
	      parFit[i]=min->GetParameter(i);
	      std::cout<<parFit[i]<<" ";

	    }
	    std::cout<<std::endl;
        TMatrixD res= euler2TMatrixD(parFit[0], parFit[1], parFit[2]);
        res.Print();

      //  printFinalResults(parFit);

        TMatrixD matr(4,4);
        matr.UnitMatrix();
        matr(0,3)=1.5;
        matr(1,3)=1.5;
        matr(2,3)=-0.3;
        res.ResizeTo(4,4);
        res(0,3)=parFit[3];
	    res(1,3)=parFit[4];
	    res(2,3)=parFit[5];
	    res(3,3)=1;
	    res=matr*res.Invert();

        std::ofstream fs("par.txt");
        fs<<"TVector\n";
        fs<<res(0,3)<<	"\n";
        fs<<res(1,3)<<	"\n";
        fs<<res(2,3)<<	"\n";
        fs<<"\n";

        fs<<"RMatrix\n";
        fs<<res(0,0)<<" "<<res(0,1)<<" "<<res(0,2)<<"\n"
        	<<res(1,0)<<" "<<res(1,1)<<" "<<res(1,2)<<"\n"
			<<res(2,0)<<" "<<res(2,1)<<" "<<res(2,2)<<"\n";
        fs<<"\n";
        fs<<"Camera Intrinsics: focal height width\n";
        fs<<parFit[6]<<" "<<parFit[8]*2<<" "<<parFit[7]*2;



}

const int grayscale_slider_max = 255;
int grayscale=12;
int hdgrayscale=170;
cv::Mat src, src_gray, dst;
cv::RNG rng(12345);

void on_trackbar( int, void* )
{

	cv::Mat canny_output;
	std::vector<std::vector<cv::Point> > contours;
	std::vector<cv::Vec4i> hierarchy;
	//cv::threshold(src_gray, dst, grayscale, grayscale_slider_max, cv::THRESH_BINARY);
	//cv::Canny( dst, canny_output, grayscale, grayscale*2, 3 );
	//cv::Canny( src_gray, canny_output, grayscale, grayscale*2, 3 );
	cv::threshold(src_gray, dst, grayscale, grayscale_slider_max, cv::THRESH_BINARY);
	cv::findContours( dst, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, cv::Point(0, 0) );
	std::vector<std::vector<cv::Point> >  filtered_contours(contours.size());

	  // copy only positive numbers:
	 // auto it = std::copy_if(foo.begin(), foo.end(), bar.begin(), [](int i){return !(i<0);} );
	//  bar.res//ize(std::distance(bar.begin(),it));  // shrink container to new size

	auto it = std::copy_if (contours.begin(), contours.end(), filtered_contours.begin(),
			[]( std::vector<cv::Point>  i){ return cv::contourArea(i) > 200 &&  cv::contourArea(i) < 2000 ;} );
	filtered_contours.resize(std::distance(filtered_contours.begin(),it));
	std::vector<std::vector<cv::Point> >  filtered_hulls;
	 for( size_t i = 0; i < filtered_contours.size(); i++ ) {
	                // approximate contour with accuracy proportional
	                // to the contour perimeter
		std::vector<cv::Point> approx;

		cv::approxPolyDP (cv::Mat(filtered_contours[i]), approx, 0.1 * sqrt(cv::contourArea(filtered_contours[i])) , true);

		std::vector<cv::Point> hull;

		cv::convexHull( cv::Mat(approx), hull, false );
		if (approx.size()  == 8 &&  hull.size() == 6) filtered_hulls.push_back(hull);


	  }
	// Draw contours
	 cv::Mat drawing = cv::Mat::zeros( dst.size(), CV_8UC3 );
	 for( int i = 0; i< filtered_hulls.size(); i++ )
	    {
	      cv::Scalar color = cv::Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
	      drawContours( drawing, filtered_hulls, i, color, 1, 8, hierarchy, 0, cv::Point() );
	     for (auto it:filtered_hulls[i] ) cv::circle(drawing, it, 1, cv::Scalar(0, 0, 255) , 1, 8);
	    }

		if (filtered_hulls.size()==1){

			std::vector<int> contDist(filtered_hulls[0].size());
			std::transform(filtered_hulls[0].begin(), filtered_hulls[0].end(), contDist.begin(), [](cv::Point po){
				return pow(po.x,2)+pow(po.y,2);


			});
		//	int index=contorDist.index(min(contoursDist))
		   int min= *std::min_element(contDist.begin(),contDist.end());
		    auto it=std::find(contDist.begin(),contDist.end(), min);
		    auto pos = it - contDist.begin();
			std::vector<cv::Point> rotatedConvexhull(filtered_hulls[0].size());
		//	std::rotate(filtered_hulls[0].begin(), filtered_hulls[0].begin()+int(pos),rotatedConvexhull.begin());
			std::rotate_copy(filtered_hulls[0].begin(), filtered_hulls[0].begin()+int(pos), filtered_hulls[0].end(), rotatedConvexhull.begin() );
		//	std::rotate_copy()

			cv::imshow("GrayScale", drawing );
			std::cout<<"start printing out\n";
			for (auto &it:filtered_hulls[0]) std::cout<< it<<std::endl;
			std::cout<<"////////\n";
			for (auto &it:rotatedConvexhull) std::cout<< it<<std::endl;

			for (auto &it:contDist) std::cout<< it<<std::endl;
			std::cout<<"end printing out\n";
		}


	 //Show in a window
	 //cv::namedWindow( "GrayScale", CV_WINDOW_AUTOSIZE );
	 //cv::resize(drawing, drawing, cv::Size(drawing.cols/2, drawing.rows/2)); // to half size or even smaller

	// cv::imshow("GrayScale", drawing );
	 cv::waitKey(0);


// cv::threshold(src,dst, grayscale, grayscale_slider_max, cv::THRESH_BINARY);

 //cv::imshow( "GrayScale", canny_output );
}

void on_trackbar2( int, void* ) {

 cv::threshold(src,dst, grayscale, grayscale_slider_max, cv::THRESH_BINARY);
 cv::namedWindow( "GrayScale", CV_WINDOW_AUTOSIZE );
 cv::imshow( "GrayScale", dst );
 int c=cv::waitKey(1);
 //cv::waitKey(0);

}







///Recording 3d model and image

int maxIr=0xFF, minIr=0;



void findMinMax(const cv::Mat &ir)
  {
    for(size_t r = 0; r < (size_t)ir.rows; ++r)
    {
      const uint16_t *it = ir.ptr<uint16_t>(r);

      for(size_t c = 0; c < (size_t)ir.cols; ++c, ++it)
      {
        minIr = std::min(minIr, (int) * it);
        maxIr = std::max(maxIr, (int) * it);
      }
    }
}
void convertIr(const cv::Mat &ir, cv::Mat &grey)
  {
    maxIr=0xFFFF;
    minIr=0;
    cv::Ptr<cv::CLAHE> clahe;
    clahe = cv::createCLAHE(1.5, cv::Size(32, 32));
    const float factor = 255.0f / (maxIr - minIr);
    cv::Mat ir_scaled;
    cv::resize(ir, ir_scaled, cv::Size(), 1, 1, cv::INTER_CUBIC);

    grey.create(ir_scaled.rows, ir_scaled.cols, CV_8U);

    		//  #pragma omp parallel for
    for(size_t r = 0; r < (size_t)ir_scaled.rows; ++r)
    {
      const uint16_t *itI = ir_scaled.ptr<uint16_t>(r);
      uint8_t *itO = grey.ptr<uint8_t>(r);

      for(size_t c = 0; c < (size_t)ir_scaled.cols; ++c, ++itI, ++itO)
      {
        *itO = std::min(std::max(*itI - minIr, 0) * factor, 255.0f);
      }
    }
    clahe->apply(grey, grey);
  }


//struct K2G_generator{
//public:
//	K2G_generator(Processor freenectprocessor, bool mirroring, char ** argv): freenectprocessor_(freenectprocessor), mirroring_(mirroring), argv_(argv),n_(0){}
//    K2G * operator ()(){return new K2G(freenectprocessor_, mirroring_, argv_[n_++ + 2]);}
//private:
//	unsigned int n_;
//	Processor freenectprocessor_;
//	bool mirroring_;
//	char ** argv_;
//};

struct K2G_generator{
public:
	K2G_generator(Processor freenectprocessor, bool mirroring, const std::vector<std::string> &argv): freenectprocessor_(freenectprocessor), mirroring_(mirroring), argv_(argv),n_(0){}
    K2G * operator ()(){return new K2G(freenectprocessor_, mirroring_, argv_.at(n_++));}
private:
	unsigned int n_;
	Processor freenectprocessor_;
	bool mirroring_;
	std::vector<std::string> argv_;
};

//struct PlySaver{
//
//  PlySaver(std::vector<boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>>> & clouds,  bool binary, bool use_camera, std::vector<K2G *> & kinects):
//           binary_(binary), use_camera_(use_camera), clouds_(clouds), kinects_(kinects){}
//
//  std::vector<boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>>> & clouds_;
//  std::vector<K2G *> & kinects_;
//  bool binary_;
//  bool use_camera_;
//};

struct PlySaver{

  PlySaver(std::vector<boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>>> & clouds, const std::vector<cv::Point> &rotatedConvexhull,  bool binary, bool use_camera, std::vector<K2G *> & kinects):
           binary_(binary), use_camera_(use_camera), clouds_(clouds), kinects_(kinects), rotatedConvexhull_(rotatedConvexhull){}

  std::vector<boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>>> & clouds_;
  const std::vector<cv::Point> & rotatedConvexhull_;
  std::vector<K2G *> & kinects_;
  bool binary_;
  bool use_camera_;
};


void
KeyboardEventOccurred(const pcl::visualization::KeyboardEvent &event, void * data)
{
  std::string pressed;
  pressed = event.getKeySym();
  PlySaver * s = (PlySaver*)data;
  if(event.keyDown())
  {
    if(pressed == "s")
    {

      pcl::PLYWriter writer;
      std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();
      std::string now = std::to_string((long)std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch()).count());
      cv::Mat color,depth, bigmat, bigmat_scaled, bigmat_grey;
      for(size_t i = 0; i < s->kinects_.size(); ++i){
      	writer.write ("cloud_"+ std::to_string(i) + "_" + now + ".ply", *(s->clouds_[i]), s->binary_, s->use_camera_);
        s->kinects_[i]->get(color, depth, bigmat, false);
      	//cv::imwrite("color_" + std::to_string(i) + "_" + now + ".jpg", color);
        // cv::imwrite("bigmat_" + std::to_string(i) + "_" + now + ".exr", bigmat);
         // cv::imwrite("depth_" + std::to_string(i) + "_" + now + ".exr", depth);

         //cv::resize(bigmat, bigmat_scaled, cv::Size(), 2.0, 2.0, cv::INTER_CUBIC);
         //findMinMax(bigmat);
         //convertIr(bigmat, bigmat_grey);
         cv::imwrite(s->kinects_[i]->getSerial()+"_"+std::to_string(i)+".png", color);
         bigmat.convertTo(bigmat, CV_16UC1);
         findMinMax(bigmat);
         convertIr(bigmat, bigmat_grey);


         cv::imwrite("ir_" + std::to_string(i) + ".png", bigmat_grey);
       ///cv::imwrite("depth_" + std::to_string(i) + ".exr", bigmat);
       //  cv::imwrite("depth_"+std::to_string(i) + ".exr", depth);

          depth.convertTo(depth, CV_16UC1);

          cv::imwrite(s->kinects_[i]->getSerial()+"_"+std::to_string(i)+ "_depth.png", depth);
      //   cv::imwrite(s->kinects_[i]->getSerial()+"_"+std::to_string(i)+ ".exr", depth);
       // cv::FileStorage file("bigmat_" + std::to_string(i) + "_" + now + ".xml", cv::FileStorage::WRITE);


         //file << bigmat;
        //std::vector<int> compression_params;
        //compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
       // compression_params.push_back(9);
     //   cv::imwrite("color_" + std::to_string(i) + "_" + now + ".png", color,compression_params);
      //   pcl::io::saveRgbPNGFile ("color_" + std::to_string(i) + "_" + now + ".png", color,512,424);

          for (auto &it:s->rotatedConvexhull_) std::cout<< it<<std::endl;

      }
      std::cout << "saved " << "cloud and color " + now << std::endl;
    }
    if(pressed == "x")
    {
        for(auto & k : s->kinects_){
        	k->storeParameters();
        }
        std::cout << "stored calibration parameters" << std::endl;
    }
  }
}

void saveToFiles(const string &filename, const std::vector<cv::Point3d>  &rotatedConvexhullIR,const std::vector<cv::Point> & rotatedConvexhullHD  ) {
//	std::ofstream fs(filename);
	cv::FileStorage store(filename, cv::FileStorage::WRITE);
	cv::write(store,"keypointsIR",rotatedConvexhullIR);
	cv::write(store,"keypointsHD",rotatedConvexhullHD);

	store.release();



}

void loadPointsFromFiles(const string &dirname, std::vector<std::pair<cv::Point3d, cv::Point>> &pointsvectorIRAndHD ) {
	using namespace boost::filesystem;

	try
	  {
	    if (exists(dirname))    // does p actually exist?
	    {
	      if (is_regular_file(dirname))        // is p a regular file?
	        cout << dirname << " size is " << file_size(dirname) << '\n';

	      else if (is_directory(dirname))      // is p a directory?
	      {
	        cout << dirname << " is a directory containing:\n";

	        typedef std::vector<path> vec;             // store paths,
	        vec v;                                // so we can sort them later

	        copy(directory_iterator(dirname), directory_iterator(), back_inserter(v));

	        sort(v.begin(), v.end());             // sort, since directory iteration
	        for (const auto& it:v)   {
	        	cout <<"file "<< it.string() <<endl;             // is not ordered on some file systems
	        	std::vector<cv::Point3d> keypointvectorIR;
	        	std::vector<cv::Point> keypointvectorHD;
	        	cv::FileStorage store(it.string(), cv::FileStorage::READ);
	        	cv::FileNode n1 = store["keypointsIR"];
	        	cv::FileNode n2 = store["keypointsHD"];
	        	cv::read(n1,keypointvectorIR);
	        	cv::read(n2,keypointvectorHD);
	        	assert(keypointvectorIR.size()==keypointvectorHD.size());
	        	for ( short i=0; i< keypointvectorIR.size();i++)


	        	pointsvectorIRAndHD.push_back(std::make_pair(keypointvectorIR[i], keypointvectorHD[i]));

	        	//for (const auto &it: keypointvectorIR) std::cout<<it.x<<" "<<it.y<<" "<<it.z<<std::endl;
	        	store.release();

	        }


	      }
	      else
	        cout << dirname << " exists, but is neither a regular file nor a directory\n";
	    }
	    else
	      cout << dirname << " does not exist\n";
	  }

	  catch (const filesystem_error& ex)
	  {
	    cout << ex.what() << '\n';
	  }

	  for (const auto &it: pointsvectorIRAndHD) {
		  std::cout<<it.first.x<<" "<<it.first.y<<" "<<it.first.z<<std::endl;
		  std::cout<<it.second.x<<" "<<it.second.y<<std::endl;
	  }
}



void convertUVs2Coordinates (const std::vector<std::pair<cv::Point3d, cv::Point>> & pointsUVsvectorIRAndHD, std::pair<cv::Mat, cv::Mat> & pointsCoordsCVMatrixIRAndHD) {

	pointsCoordsCVMatrixIRAndHD.first=cv::Mat();
	pointsCoordsCVMatrixIRAndHD.second=cv::Mat();

	double fx_ir=364.963,fy_ir=364.963,cx_ir=257.487,cy_ir=199.466;
	double fx_hd=1081.37,fy_hd=1081.37,cx_hd=959.5,cy_hd=539.5;

	for (auto &it: pointsUVsvectorIRAndHD) {
		double c = it.first.x;
		double r = it.first.y;
		double depth_val=it.first.z/1000.;

		if (!isnan(depth_val) && depth_val >= 0.001)
		  {
			double x_ir = (c + 0.5 - cx_ir) * (1./fx_ir) * depth_val;
			double y_ir = (r + 0.5 - cy_ir) * (1./fy_ir) * depth_val;
			double z_ir = depth_val;
			//std::cout<<"Coords "<<c<<" "<<r<<" "<<x_ir<<" "<<y_ir<<" "<<z_ir<<std::endl;
//			double x_hd = (c + 0.5 - cx_ir) * fx_ir * depth_val;
//			double y_ir = (r + 0.5 - cy_ir) * fy_ir * depth_val;
//			double z_ir = depth_val;
			cv::Mat tmpmatir(1,3,CV_32F);
			tmpmatir.at<double>(0,0)=x_ir;
			tmpmatir.at<double>(0,1)=y_ir;
			tmpmatir.at<double>(0,2)=z_ir;



			pointsCoordsCVMatrixIRAndHD.first.push_back(tmpmatir);
			cv::Mat tmpmathd(1,2,CV_32S);

			tmpmathd.at<int>(0,0)=it.second.x;
			tmpmathd.at<int>(0,1)=it.second.y;


			pointsCoordsCVMatrixIRAndHD.second.push_back(tmpmathd);

		  }


	}

	pointsCoordsCVMatrixIRAndHD.first=pointsCoordsCVMatrixIRAndHD.first.t();
	cv::Mat tmpsones= cv::Mat::ones(1,pointsCoordsCVMatrixIRAndHD.first.cols , CV_32F);
	pointsCoordsCVMatrixIRAndHD.first.push_back(tmpsones);
	pointsCoordsCVMatrixIRAndHD.second=pointsCoordsCVMatrixIRAndHD.second.t();
}

void convertUVs2Coordinates (const std::vector<std::pair<cv::Point3d, cv::Point>> & pointsUVsvectorIRAndHD, std::pair<TMatrixD, TMatrixD> & pointsCoordsTMatrixIRAndHD  ) {


	double fx_ir=364.963,fy_ir=364.963,cx_ir=512./2,cy_ir=424./2;
	//double fx_ir=364.963,fy_ir=364.963,cx_ir=257.487,cy_ir=199.466;
	double fx_hd=1081.37,fy_hd=1081.37,cx_hd=959.5,cy_hd=539.5;

	std::vector<double> vecOfCoordinatesIR;
	std::vector<double> vecOfCoordinatesHD;
     int i=0, j=0;
	for (auto &it: pointsUVsvectorIRAndHD) {
		double c = it.first.x;
		double r = it.first.y;
		double depth_val=it.first.z/1000.;

		if (!isnan(depth_val) && depth_val >= 0.001 /*&& j%4==0*/)
		  {
			double x_ir = (c + 0.5 - cx_ir) * (1./fx_ir) * depth_val;
			double y_ir = (r + 0.5 - cy_ir) * (1./fy_ir) * depth_val;
			double z_ir = depth_val;
			//std::cout<<"Coords "<<c<<" "<<r<<" "<<x_ir<<" "<<y_ir<<" "<<z_ir<<std::endl;
//			double x_hd = (c + 0.5 - cx_ir) * fx_ir * depth_val;
//			double y_ir = (r + 0.5 - cy_ir) * fy_ir * depth_val;
//			double z_ir = depth_val;
			//cv::Mat tmpmatir(1,3,CV_32F);

			vecOfCoordinatesIR.push_back(x_ir);vecOfCoordinatesIR.push_back(y_ir);vecOfCoordinatesIR.push_back(z_ir);

			vecOfCoordinatesHD.push_back(it.second.x);vecOfCoordinatesHD.push_back(it.second.y);
			i++;

		  }

		j++;
	}

	pointsCoordsTMatrixIRAndHD.first.ResizeTo(i,3);
	pointsCoordsTMatrixIRAndHD.first.SetMatrixArray(vecOfCoordinatesIR.data());
	pointsCoordsTMatrixIRAndHD.first.T();
	pointsCoordsTMatrixIRAndHD.second.ResizeTo(i, 2 );
	pointsCoordsTMatrixIRAndHD.second.SetMatrixArray(vecOfCoordinatesHD.data());
//	pointsCoordsTMatrixIRAndHD.second.T();


}


int main(int argc, char **argv)
{
//	 src = cv::imread("/home/hamit/pcl-sw/libfreenect2pclgrabber/build/ir_0.png");
//   // src = cv::imread("/home/hamit/pcl-sw/libfreenect2pclgrabber/build/291296634347_0.png");
//	 if( !src.data ) { printf("Error loading src1 \n"); return -1; }
//	 cv::cvtColor( src, src_gray, CV_BGR2GRAY );
//	 cv::blur( src_gray, src_gray, cv::Size(3,3) );
//
////     cv::imshow("GrayScale", src);
////	 cv::waitKey(0);
//	 cv::namedWindow("GrayScale", 1);
//	 char TrackbarName[50];
//	 sprintf( TrackbarName, "Max: %d", grayscale_slider_max );
//	 grayscale=0;
//	 cv::createTrackbar( TrackbarName, "GrayScale", &grayscale, grayscale_slider_max, on_trackbar );
//	 on_trackbar( grayscale, 0 );
//	 cv::waitKey(0);
//
//return 0;






	  std::cout << "Syntax is: " << argv[0] << " [-processor 0|1|2] -processor options 0,1,2,3 correspond to CPU, OPENCL, OPENGL, CUDA respectively" << std::endl;
	   std::cout << "followed by the kinect2 serials and a last 1 for mirroring" << std::endl;
	   std::cout << "Press \'s\' to store both clouds and color images." << std::endl;
	   std::cout << "Press \'x\' to store both calibrations." << std::endl;

	   Processor freenectprocessor= OPENGL;
	   int processor;
	   pcl::console::parse_argument (argc, argv, "-processor", processor);
	   freenectprocessor = static_cast<Processor>(processor);
	   bool mirroring=false;
	   pcl::console::parse_argument (argc, argv, "-mirror", mirroring);
	   if(mirroring)
	  	   	std::cout << "mirroring enabled" << std::endl;
	   std::vector<std::string> serials;
	   if ( pcl::console::parse_multiple_arguments(argc, argv, "-serials", serials)==0) {

		   PCL_ERROR("No Kinect v2 serial given\n");
		   exit(-1);

      }
	  // bool loadFilesANdCalibrate=false;
	   int isloadFilesAndCalibrate = pcl::console::find_argument(argc, argv, "-isloadAndSave");

	 //  serials={"291296634347"};
	   if(argc < 3){
	     std::cout << "Wrong syntax! specify at least processor and one serial" << std::endl;
	     return -1;
	   }

	   std::cout<<"isloadFilesAndCalibrate "<<isloadFilesAndCalibrate<<std::endl;
	   if (isloadFilesAndCalibrate > 0){
		   loadPointsFromFiles("/home/hamit/calibrateKinectv2/bin/test3",pointsUVvectorIRAndHD);
		   convertUVs2Coordinates(pointsUVvectorIRAndHD, pointsCoordsTMatrixIRAndHD);
		  // UVs = cv::Mat(pointsCoordsMatIRAndHD.second.rows, pointsCoordsMatIRAndHD.second.cols, CV_32S);
		   pointsCoordsTMatrixIRAndHD.first.Print();
		   pointsCoordsTMatrixIRAndHD.second.Print();

		   minimize();
		   return 0;

	   }



	   // Set te processor to input value or to default OPENGL
//	   Processor freenectprocessor = OPENGL;
//	   freenectprocessor = static_cast<Processor>(atoi(argv[1]));

	   // check if mirroring is enabled
// 	   bool mirroring = atoi(argv[argc - 1]) == 1 ? true : false;
//	   if(mirroring)
//	   	std::cout << "mirroring enabled" << std::endl;

	   // Count number of input serials which represent the number of active kinects
	  // int kinect2_count = mirroring ? argc - 3 : argc - 2;
	   int kinect2_count = static_cast<int>(serials.size());
	   std::cout << "loading " << kinect2_count << " devices" << std::endl;
	  // kinect2_count=1;
	   // Initialize container structures
	   std::vector<K2G *> kinects(kinect2_count);
	   std::vector<boost::shared_ptr<pcl::PointCloud<pcl::PointXYZRGB>>> clouds(kinect2_count);
	   std::vector<cv::Point> rotatedConvexhull;
	  // boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer ("3D viewer"));
	  // viewer->setBackgroundColor (0, 0, 0);

	   // Generate all the kinect grabbers
	   //K2G_generator kinect_generator(freenectprocessor, mirroring, argv);
	   K2G_generator kinect_generator(freenectprocessor, mirroring, serials);

	   std::generate(kinects.begin(), kinects.end(), kinect_generator);
	   // Initialize clouds and viewer viewpoints
	   for(size_t i = 0; i < kinect2_count; ++i)
	   {
	   	clouds[i] = kinects[i]->getCloud();
	   	kinects[i]->printParameters();

	   	// clouds[i]->sensor_orientation_.w() = 0.0;
	   	// clouds[i]->sensor_orientation_.x() = 1.0;
	   	// clouds[i]->sensor_orientation_.y() = 0.0;
	   	// clouds[i]->sensor_orientation_.z() = 0.0;
	  //   std::cout<<" cloud point: "<<clouds[i]->points[10006].z;
	   //	viewer->addPointCloud<pcl::PointXYZRGB>(clouds[i], pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB>(clouds[i]), "sample cloud_" + i);
	  // 	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "sample cloud_" + i);
	   }

	   // Add keyboards callbacks
	   PlySaver ps(clouds, rotatedConvexhull, true, false, kinects);
	  // viewer->registerKeyboardCallback(KeyboardEventOccurred, (void*)&ps);
	   std::cout << "starting cycle" << std::endl;

	   std::chrono::high_resolution_clock::time_point tnow, tpost, tpost_prev;
	   tpost_prev=std::chrono::high_resolution_clock::now();
	   std::vector<pcl::visualization::Camera> cam;


	   cv::Mat color[kinect2_count], color_grey[kinect2_count],  color_hd[kinect2_count], color_hd_grey[kinect2_count], depth[kinect2_count], irmat[kinect2_count], irmat_scaled[kinect2_count], irmat_grey[kinect2_count];



	bool isShowGrayScale = true;
	bool isShowHDColor = false;
	bool exit = false;
	bool createdIRWindow=false, createdHdWindow=false;
	std::vector<cv::Point3d> rotatedConvexhullIR;
	std::vector<cv::Point> rotatedConvexhullHD;
	unsigned int fileIndex=0;


	//cv::resizeWindow("HDGrayScale", 900,900);


	while( !exit /*!viewer->wasStopped()*/ ){

	    //   viewer->spinOnce ();

	    // cout<< viewer->getViewerPose().linear();
	    // viewer->getCameras(cam);



	     tnow = std::chrono::high_resolution_clock::now();

	     for(size_t i = 0; i < kinect2_count; ++i) {
	     //	clouds[i] = kinects[i]->updateCloud(clouds[i]);
	        kinects[i]->get(color_hd[i], color[i], depth[i], irmat[i]);


	     }


	     tpost = std::chrono::high_resolution_clock::now();
	    // std::cout << "<<< delta: " << std::chrono::duration_cast<std::chrono::milliseconds>(tpost.time_since_epoch()).count()-std::chrono::duration_cast<std::chrono::milliseconds>(tpost_prev.time_since_epoch()).count()<<std::endl;

	     tpost_prev=tpost;

	     for(size_t i = 0; i < kinect2_count; ++i) {
	    //  pcl::transformPointCloud(*clouds[i], *clouds[i], transform_1 );




	    	 cv::namedWindow("GrayScale", CV_WINDOW_AUTOSIZE);
	    	 cv::namedWindow("HDGrayScale", CV_WINDOW_FREERATIO);

	         if (isShowGrayScale) {

				 char trackbarNameIR[50];
				 sprintf( trackbarNameIR, "TrackbarGrayScale");

				 if (!createdIRWindow) {
					 cv::moveWindow("GrayScale", 600,400);
					 createdIRWindow=true;
				 }

				 //cv2.resizeWindow(windowName, 900,900)
				 grayscale=65;
				 cv::createTrackbar( trackbarNameIR, "GrayScale", &grayscale, 255, NULL);

//				cv::Mat canny_output;
				std::vector<std::vector<cv::Point> > contours;
				std::vector<cv::Vec4i> hierarchy;
//				irmat[i].convertTo(irmat[i], CV_16UC1);
//				findMinMax(irmat[i]);
//				convertIr(irmat[i], irmat_grey[i]);
				// cv::cvtColor(irmat_grey[i] ,  irmat_grey[i], CV_BGR2GRAY );
				// cv::blur( irmat_grey[i] ,  irmat_grey[i], cv::Size(3,3) );
				 cv::cvtColor(color[i] , color_grey[i], CV_BGR2GRAY );


			    int grayScale = cv::getTrackbarPos( trackbarNameIR, "GrayScale");
			//	cv::threshold(irmat_grey[i], irmat_grey[i], grayScale, grayscale_slider_max, cv::THRESH_BINARY);
				cv::threshold(color_grey[i], color_grey[i], grayScale, grayscale_slider_max, cv::THRESH_BINARY);


    ///Find chessboard board
//
//				cv::namedWindow( "GrayScale", 1 );
//				cv::Size boardSize(5,7);
//				std::vector<std::vector<cv::Point2f> > imagePoints(1);
//				  bool found = cv::findChessboardCorners(irmat_grey[i], boardSize, imagePoints[0]);
//				  if(!found)
//				    {
//				    std::cerr << "Could not find chess board!" << std::endl;
//
//				    }
//
//				  cv::cvtColor(irmat_grey[i],irmat_grey[i] , CV_GRAY2BGR);
//				cv::drawChessboardCorners(irmat_grey[i], boardSize, cv::Mat(imagePoints[0]), found );
//
//				cv::imshow("GrayScale", irmat_grey[i] );
//
//				if (cv::waitKey(30)  > 0) {
//								   cv::destroyAllWindows();
//								   isShowGrayScale=false;
//				}

				// Find the figure provided
				cv::findContours(color_grey[i] , contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, cv::Point(0, 0) );
				std::vector<std::vector<cv::Point> >  filtered_contours(contours.size());



				auto it = std::copy_if (contours.begin(), contours.end(), filtered_contours.begin(),
						[]( std::vector<cv::Point>  i){ return cv::contourArea(i) > 100 &&  cv::contourArea(i) < 3000 ;} );
				filtered_contours.resize(std::distance(filtered_contours.begin(),it));
				std::vector<std::vector<cv::Point> >  filtered_hulls;
				for( size_t i = 0; i < filtered_contours.size(); i++ ) {
				                // approximate contour with accuracy proportional
				                // to the contour perimeter
					std::vector<cv::Point> approx;

					cv::approxPolyDP (cv::Mat(filtered_contours[i]), approx, 0.1 * sqrt(cv::contourArea(filtered_contours[i])) , true);

					std::vector<cv::Point> hull;

					cv::convexHull( cv::Mat(approx), hull, false );
					if (approx.size()  == 8 &&  hull.size() == 6) filtered_hulls.push_back(hull);


				  }
				// if (filtered_hulls.size()!=1) break;
				// Draw contours

//				for( int i = 0; i< filtered_hulls.size(); i++ )
//				    {
//				      cv::Scalar color = cv::Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
//				      cv::drawContours( drawing, filtered_hulls, i, color, 1, 8, hierarchy, 0, cv::Point() );
//				      for (auto it:filtered_hulls[i] ) cv::circle(drawing, it, 1, cv::Scalar(0, 0, 255) , 1, 8);
//				    }

				// cv::imshow("GrayScale",irmat_grey[i]);
				// cv::namedWindow( "GrayScale", CV_WINDOW_AUTOSIZE );

				if (filtered_hulls.size()==1){
//					cv::imshow("GrayScale", drawing );
//					std::cout<<"start printing out\n";
//					for (auto &it:filtered_hulls[0]) std::cout<< it<<std::endl;
//					std::cout<<"end printing out\n";
					cv::Mat drawing = cv::Mat::zeros( color_grey[i].size(), CV_8UC3 );
					for( int i = 0; i< filtered_hulls.size(); i++ )
					{
					  cv::Scalar color = cv::Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
					  cv::drawContours( drawing, filtered_hulls, i, color, 1, 8, hierarchy, 0, cv::Point2d() );
					  for (auto it:filtered_hulls[i] ) cv::circle(drawing, it, 1, cv::Scalar(0, 0, 255) , 1, 8);
					}

//					std::vector<int> contDist(filtered_hulls[0].size());
//					std::transform(filtered_hulls[0].begin(), filtered_hulls[0].end(), contDist.begin(), [](cv::Point po){
//						return pow(po.x,2)+pow(po.y,2);
//
//
//					});
					std::vector<int> contDist;
					std::for_each(filtered_hulls[0].begin(), filtered_hulls[0].end(), [&contDist](cv::Point po){

									contDist.push_back(pow(po.x,2)+pow(po.y,2));


								});

//					for (auto it:contDist )  std::cout<< it<<" ";
//					std::cout<<std::endl;
//					for (auto it:contDist2 ) std::cout<< it<<" ";
//					std::cout<<std::endl;

					int min= *std::min_element(contDist.begin(),contDist.end());
					auto it=std::find(contDist.begin(),contDist.end(), min);
					auto pos = it - contDist.begin();
					std::vector<cv::Point> rotatedConvexhull(filtered_hulls[0].size());
					rotatedConvexhullIR.clear();rotatedConvexhullIR.resize(filtered_hulls[0].size());
					std::rotate_copy(filtered_hulls[0].begin(), filtered_hulls[0].begin()+int(pos), filtered_hulls[0].end(), rotatedConvexhull.begin() );

					std::transform(rotatedConvexhull.begin(),rotatedConvexhull.end(), rotatedConvexhullIR.begin(), [&depth, i](cv::Point po){
						return cv::Point3d(po.x, po.y, depth[i].at<double>(po.x, po.y));


					});

					cv::imshow("GrayScale", drawing );
					cv::imshow("HDGrayScale",color_hd[i]);
//					std::cout<<"start printing out\n";
//					for (auto &it:filtered_hulls[0]) std::cout<< it<<std::endl;
//					std::cout<<"////////\n";
		//		for (auto &it:rotatedConvexhullIR) std::cout<< it<<std::endl;
//
//					for (auto &it:contDist) std::cout<< it<<std::endl;
//					std::cout<<"end printing out\n";



				}



                char c=(char)cv::waitKey(30);
				if ( c =='q') {
					//cv::destroyWindow("GrayScale");
					//cv::destroyWindow("HDGrayScale");
					cv::destroyAllWindows();
					isShowGrayScale=false;
					isShowHDColor=false;
					exit=true;
				} else if  ( c =='s') {
					//cv::destroyWindow("GrayScale");
					cv::destroyAllWindows();
					isShowGrayScale=false;
					isShowHDColor=true;
					createdIRWindow=false;

				}




	         }


			if (!isShowHDColor) continue;
			hdgrayscale=65;
			cv::namedWindow("HDGrayScale", CV_WINDOW_FREERATIO);


	 		for (;isShowHDColor;) {

				//trackbar for HD image
				char trackbarNameHD[50];

				sprintf( trackbarNameHD, "TrackbarHDColor");

				if (!createdHdWindow) {
					cv::moveWindow("HDGrayScale", 600,400);

					createdHdWindow=true;
				}
				cv::createTrackbar( trackbarNameHD, "HDGrayScale", &hdgrayscale, 255, NULL);

				cv::Mat canny_output;
				std::vector<std::vector<cv::Point> > contours;
				std::vector<cv::Vec4i> hierarchy;

				cv::cvtColor(color_hd[i] , color_hd_grey[i], CV_BGR2GRAY );
				//cv::blur( color_grey[i] ,  color_grey[i], cv::Size(3,3) );
				int grayScale = cv::getTrackbarPos( trackbarNameHD, "HDGrayScale");
				cv::threshold(color_hd_grey[i], color_hd_grey[i], grayScale, grayscale_slider_max, cv::THRESH_BINARY);


	///Find chessboard board
//
//				cv::namedWindow( "GrayScale", 1 );
//				cv::Size boardSize(5,7);
//				std::vector<std::vector<cv::Point2f> > imagePoints(1);
//				  bool found = cv::findChessboardCorners(irmat_grey[i], boardSize, imagePoints[0]);
//				  if(!found)
//				    {
//				    std::cerr << "Could not find chess board!" << std::endl;
//
//				    }
//
//				  cv::cvtColor(irmat_grey[i],irmat_grey[i] , CV_GRAY2BGR);
//				cv::drawChessboardCorners(irmat_grey[i], boardSize, cv::Mat(imagePoints[0]), found );
//
//				cv::imshow("GrayScale", irmat_grey[i] );
//
//				if (cv::waitKey(30)  > 0) {
//								   cv::destroyAllWindows();
//								   isShowGrayScale=false;
//				}

				// Find the figure provided
				cv::findContours(color_hd_grey[i] , contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, cv::Point(0, 0) );
				std::vector<std::vector<cv::Point> >  filtered_contours(contours.size());



				auto it = std::copy_if (contours.begin(), contours.end(), filtered_contours.begin(),
						[]( std::vector<cv::Point>  i){ return cv::contourArea(i) > 1000 &&  cv::contourArea(i) < 50000 ;} );
				filtered_contours.resize(std::distance(filtered_contours.begin(),it));
				std::vector<std::vector<cv::Point> >  filtered_hulls;
				for( size_t i = 0; i < filtered_contours.size(); i++ ) {
								// approximate contour with accuracy proportional
								// to the contour perimeter
					std::vector<cv::Point> approx;

					cv::approxPolyDP (cv::Mat(filtered_contours[i]), approx, 0.1 * sqrt(cv::contourArea(filtered_contours[i])) , true);

					std::vector<cv::Point> hull;

					cv::convexHull( cv::Mat(approx), hull, false );
					if (approx.size()  == 8 &&  hull.size() == 6) filtered_hulls.push_back(hull);


				  }
				// if (filtered_hulls.size()!=1) break;
				// Draw contours

//				for( int i = 0; i< filtered_hulls.size(); i++ )
//				    {
//				      cv::Scalar color = cv::Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
//				      cv::drawContours( drawing, filtered_hulls, i, color, 1, 8, hierarchy, 0, cv::Point2f() );
//				      for (auto it:filtered_hulls[i] ) cv::circle(drawing, it, 1, cv::Scalar(0, 0, 255) , 1, 8);
//				    }

				// cv::imshow("GrayScale",irmat_grey[i]);
				// cv::namedWindow( "GrayScale", CV_WINDOW_AUTOSIZE );

				if (filtered_hulls.size()==1){
//					cv::imshow("GrayScale", drawing );
//					std::cout<<"start printing out\n";
//					for (auto &it:filtered_hulls[0]) std::cout<< it<<std::endl;
//					std::cout<<"end printing out\n";
					cv::Mat drawing = cv::Mat::zeros( color_hd_grey[i].size(), CV_8UC3 );
					//std::cout<<"color_grey[i].size() "<<color_grey[i].size()<<std::endl;
					for( int i = 0; i< filtered_hulls.size(); i++ )
					{
					  cv::Scalar color = cv::Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
					  cv::drawContours( drawing, filtered_hulls, i, color, 1, 8, hierarchy, 0, cv::Point() );
					  for (auto it:filtered_hulls[i] ) cv::circle(drawing, it, 1, cv::Scalar(0, 0, 255) , 1, 8);
					}

					std::vector<int> contDist(filtered_hulls[0].size());
					std::transform(filtered_hulls[0].begin(), filtered_hulls[0].end(), contDist.begin(), [](cv::Point po){
						return pow(po.x,2)+pow(po.y,2);


					});
					int min= *std::min_element(contDist.begin(),contDist.end());
					auto it=std::find(contDist.begin(),contDist.end(), min);
					auto pos = it - contDist.begin();
					//std::vector<cv::Point2f> rotatedConvexhullHD(filtered_hulls[0].size());
					rotatedConvexhullHD.clear();rotatedConvexhullHD.resize(filtered_hulls[0].size());
					std::rotate_copy(filtered_hulls[0].begin(), filtered_hulls[0].begin()+int(pos), filtered_hulls[0].end(), rotatedConvexhullHD.begin() );

					cv::imshow("HDGrayScale", drawing );
					//cv::resizeWindow("HDGrayScale", 1000,600);
					//cv::resizeWindow("HDGrayScale",960,540);

//					std::cout<<"start printing out\n";
//					for (auto &it:filtered_hulls[0]) std::cout<< it<<std::endl;
//					std::cout<<"////////\n";
//                	for (auto &it:rotatedConvexhullHD) std::cout<< it<<std::endl;
//
//					for (auto &it:contDist) std::cout<< it<<std::endl;
//					std::cout<<"end printing out\n";



				}


				char c=(char)cv::waitKey(30);
				if ( c == 's' ) {
					fileIndex++;
					std::string filename="IR_Color_Points_"+std::to_string(fileIndex)+".xml";
					saveToFiles(filename, rotatedConvexhullIR, rotatedConvexhullHD);
					cv::destroyWindow("HDGrayScale");
					isShowGrayScale=true;
					isShowHDColor=false;
					createdHdWindow=false;

				}
				else if ( c == 'q' ){
//					fileIndex++;
//					std::string filename="IR_Color_Points_"+std::to_string(fileIndex)+".xml";
//					saveToFiles(filename, rotatedConvexhullIR, rotatedConvexhullHD);
					cv::destroyWindow("HDGrayScale");
					isShowHDColor=false;
					isShowGrayScale=true;
				//	exit=true;

				//	loadPointsFromFiles("./",pointvectorIRAndHD);

				}








	 		}

	     }
	   }

	   // Close all kinect grabbers
	   for(auto & k : kinects)
	   	k->shutDown();



	 return 0;
}


