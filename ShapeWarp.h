//	Version 2.0
//  ShapeWarp.hpp
//  LinearWarp
//
//  Created by t_wangju on 24/10/2016.
//  The 2.0 Version Completed on 21/11/2016.
//  Copyright © 2016 t_wangju. All rights reserved.
//

#ifndef ShapeWarp_H
#define ShapeWarp_H

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include <map>

#include <GL/glew.h>
#include <GL/glut.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/SparseExtra>
#include <opencv2/opencv.hpp>

#include "PredictedPath.h"

using namespace PredictedPath;
using namespace Eigen;

typedef std::vector<Sample> Curve;

class ShapeWarp{
	enum ChooseKind:int{CENTRE,EXIST,INSERT,NONE};
	enum PoissonKind:int{WHOLE,BLOCK};
	enum SamplerKind:int{DISTANCE,SHAPE};
	enum AlgorithmKind:int{POISSONIMAGEEDITING,TRANSFINITEINTERPOLATION};
	//Apwws up 1:not needed to flush all edge each time
	enum FlushEdge:int{NO,LEFT,BOTTOM,RIGHT,TOP,ALL};
	const double HandleRadius=5.0;
	
public:
	ShapeWarp(double xmin, double ymin, double xmax,double ymax, int cols, int rows, int dims);
	~ShapeWarp();
	int Choose(double x, double y);
	int Set(double x, double y);
	int Delete(double x, double y);
	int Flush();
	inline void SetDelta(double x){
		this->delta=x;
	}
	inline const double GetDelta(){
		return delta;
	}
	inline void ChangeAlgorithm(){
		switch (AK) {
			case POISSONIMAGEEDITING:
				AK=TRANSFINITEINTERPOLATION;
				break;
			case TRANSFINITEINTERPOLATION:
				AK=POISSONIMAGEEDITING;
				break;
			default:
				return;
		}
		return;
	}
	Curve ReadControlPoints();
	//users may forget to release Matrix or Load a Matrix which has been released;
	//I can solve this problem,but it is a little complex...
	inline const double*** const ReadAMatrix()const{return const_cast<const double***>(AMatrix);}
	inline int Cols(){return cols;}
	inline int Rows(){return rows;}
	inline Sample MIApex(){return MiApex;}
	
private:

	FlushEdge FE = ALL;
	MatrixXd AInverse;
	
	AlgorithmKind AK=POISSONIMAGEEDITING;
	double*** AMatrix;
	const int cols,rows,dims;// key matrix;
	Sample MiApex;// Point's centre coordinate;
	
	ShapeWarp(const ShapeWarp &sw);
	ShapeWarp& operator=(const ShapeWarp &sw) const ;
	double delta;
	Sample CMApex;
	double xmove;// space the centre point coordinate by the movement it is.
	double ymove;
	
	Curve LBound,BBound,RBound,TBound;// always anticlockwise order;
	std::vector<double> LSect,BSect,RSect,TSect;
	mutable Curve BCurve,RCurve,TCurve,LCurve;// keep points' number unchange;
	
	double dis(Sample s1, Sample s2) const;
	double dis(Sample s1, Sample s2, Sample s) const;
	
	double GaussCE(double distance,double delta);
	int Liquidfy(Sample from ,Sample to, int ox, int oy, double delta);
	
	ChooseKind CK=NONE;
	Sample* choose=NULL;
	
	Curve CurveCut(Curve &c, const Sample s, double &totaldis) const;
	Curve CurveCut(Curve &c, double distance) const ;
	Curve Sampler(Curve c, const Curve &b, const std::vector<double> &sects) const;
	
	// all the functions below use a FULL CURVE except CHOOSE;
	bool ChooseCurve(Sample s, Curve &c);
	bool InsertCurve(Sample s, Curve &bound, Sample end, std::vector<double> &sect);
	// use CHOOSE behind INSERT to make sure the pointer is valued;
	Curve RefineCurve(const Curve &bound, const std::vector<double> &sect) const;
	
	int Bind();
	int Refine();
	int Smooth();
	double Cross(Sample v1, Sample v2, Sample u1, Sample u2) const ;
	
	// these may be needed to make sure every point can be move;
	std::vector<double> XMove;
	std::vector<double> YMove;
	// question:
	// how to compute a point in grid G';
	// when model a curve edge make the F:G->G';
	// make a f:p->p' when F:G->G' is easy;
	// but compute out the f':p'->p may be a very difficult challenge;
	std::vector<Sample> Basis;
	
private:
	int PoissonRefine();
	int UnFlip();
};

#endif /*  ShapeWarp_H */
