//  Version 2.0
//  ShapeWarp.cpp
//  LinearWarp
//
//  Created by t_wangju on 24/10/2016.
//  Copyright © 2016 t_wangju. All rights reserved.
//

#include "ShapeWarp.h"
#include "PredictedPath.h"
#include <sys/timeb.h>
#include <stdio.h>

Curve KCurve(Curve c){
	return PredictedPath::SolvePath(c,10,40);
}

ShapeWarp::ShapeWarp(double xmin, double ymin, double xmax,double ymax, int icols, int irows, int idims):cols(icols+1),rows(irows+1),dims(idims){
	printf("ShapeWarp Construct\n"); 
	Sample LBApex(xmin,ymin);
	Sample LTApex(xmin,ymax);
	Sample RBApex(xmax,ymin);
	Sample RTApex(xmax,ymax);
	
	LBound.clear();
	LBound.push_back(LTApex);
	BBound.clear();
	BBound.push_back(LBApex);
	RBound.clear();
	RBound.push_back(RBApex);
	TBound.clear();
	TBound.push_back(RTApex);
	
	LSect.clear();
	BSect.clear();
	RSect.clear();
	TSect.clear();
	LSect.push_back(irows);
	BSect.push_back(icols);
	RSect.push_back(irows);
	TSect.push_back(icols);
	
	AMatrix=new double** [cols];
	for(int i=0;i<cols;i++){
		AMatrix[i]=new double* [rows];
		for(int j=0;j<rows;j++){
			AMatrix[i][j]=new double [dims];
		}
	}
	delta=sqrt(cols*rows)/4;
	Bind();

	printf("Init\n");
	int n=cols*rows;
	int m=(cols-2)*(rows-2);
	int l=n-m;
	MatrixXd Relationship(n,n);
	Relationship.setZero();
	MatrixXd Topology(n,n);
	Topology.setZero();

	int q=m;
	for(int i=1;i<rows-1;i++){
		for(int j=1;j<cols-1;j++){
			int I=(i-1)*(cols-2)+j-1;
			Topology(I,I)=4;
			if(i==1){
				Topology(I,m+j)=-1;
				Topology(I,I+cols-2)=-1;
			}else if(i==rows-2){
				Topology(I,m+cols-1+rows-1+cols-1-j)=-1;
				Topology(I,I-cols+2)=-1;
			}else{
				Topology(I,I+cols-2)=-1;
				Topology(I,I-cols+2)=-1;
			}
			if(j==cols-2){
				Topology(I,m+cols-1+i)=-1;
				Topology(I,I-1)=-1;
			}else if(j==1){
				Topology(I,m+cols-1+rows-1+cols-1+rows-1-i)=-1;
				Topology(I,I+1)=-1;
			}else{
				Topology(I,I+1)=-1;
				Topology(I,I-1)=-1;	
			}
			std::cout<<i<<" "<<j<<std::endl; 
		}
	}
	
	for(int i=m;i<n;i++){
		Topology(i,i)=1;
	}
	
	MatrixXd A00=Topology.block(0,0,m,m);
	MatrixXd A01=Topology.block(0,m,m,l);
	MatrixXd A10=Topology.block(m,0,l,m);
	MatrixXd A11=Topology.block(m,m,l,l);

	std::cout<<"Block over"<<std::endl;
	AInverse.resize(m,l);
	AInverse.setZero();
	std::cout<<"AInverse Resize over"<<std::endl;
	AInverse=-1*(A00.inverse())*(A01)*(A11.inverse());

	PoissonRefine();
	MiApex.x=AMatrix[cols/2][rows/2][0];
	MiApex.y=AMatrix[cols/2][rows/2][1];
	CMApex.x=AMatrix[cols/2][rows/2][0];
	CMApex.y=AMatrix[cols/2][rows/2][1];
	xmove=0;
	ymove=0;	
}

ShapeWarp::ShapeWarp(const ShapeWarp & sw):cols(sw.cols),rows(sw.rows),dims(sw.dims){
	
	LBound=sw.LBound;
	BBound=sw.BBound;
	RBound=sw.RBound;
	TBound=sw.TBound;
	LSect=sw.LSect;
	BSect=sw.BSect;
	RSect=sw.RSect;
	TSect=sw.TSect;
	
	AMatrix=new double** [cols];
	for(int i=0;i<cols;i++){
		AMatrix[i]=new double* [rows];
		for(int j=0;j<rows;j++){
			AMatrix[i][j]=new double [dims];
			for(int k=0;k<dims;k++){
				AMatrix[i][j][k]=sw.AMatrix[i][j][k];
			}
		}
	}
	
	Bind();
	AInverse=sw.AInverse;
	PoissonRefine();
	MiApex.x=AMatrix[cols/2][rows/2][0];
	MiApex.y=AMatrix[cols/2][rows/2][1];
	CMApex.x=sw.CMApex.x;
	CMApex.y=sw.CMApex.y;
	xmove=sw.xmove;
	ymove=sw.ymove;
}

ShapeWarp::~ShapeWarp(){
	for(int i=0;i<rows;i++){
		for(int j=0;j<cols;j++){
			delete [] AMatrix[i][j];
		}
		delete [] AMatrix[i];
	}
	delete [] AMatrix;
}

int ShapeWarp::Bind(){
	printf("Bind\n");
	FE=ALL;
	std::vector<Curve> LPCurves,RPCurves,BPCurves,TPCurves;
	// re constructing four edges;
	// should judge whether the sect number and bound number are same;

	//For Speed Up This should be rewrite;
	//And may be not needed predicted strock to complete smooth edge curve;
	//And a vector Curve may have huge improvement in Curve Intergal;
	if(FE==ALL||FE==BOTTOM){
		BCurve=BBound;
		BCurve.push_back(*RBound.begin());
		BCurve=RefineCurve(BCurve, BSect);
		for(int i=0,j=0,k=0;i<cols-1;i++,k++){
			AMatrix[i][j][0]=BCurve[k].x;
			AMatrix[i][j][1]=BCurve[k].y;
		}
	}
	if(FE==ALL||FE==RIGHT){
		RCurve=RBound;
		RCurve.push_back(*TBound.begin());
		RCurve=RefineCurve(RCurve, RSect);
		for(int i=cols-1,j=0,k=0;j<rows-1;j++,k++){
			AMatrix[i][j][0]=RCurve[k].x;
			AMatrix[i][j][1]=RCurve[k].y;
		}
	}
	if(FE==ALL||FE==TOP){
		TCurve=TBound;
		TCurve.push_back(*LBound.begin());
		TCurve=RefineCurve(TCurve, TSect);
		for(int i=cols-1,j=rows-1,k=0;i>0;i--,k++){
			AMatrix[i][j][0]=TCurve[k].x;
			AMatrix[i][j][1]=TCurve[k].y;
		}
	}
	if(FE==ALL||FE==LEFT){
		LCurve=LBound;
		LCurve.push_back(*BBound.begin());
		LCurve=RefineCurve(LCurve, LSect);
		for(int i=0,j=rows-1,k=0;j>0;j--,k++){
			AMatrix[i][j][0]=LCurve[k].x;
			AMatrix[i][j][1]=LCurve[k].y;
		}
	}
	return 0;
}

int ShapeWarp::PoissonRefine(){
	printf("PoissonRefine\n");
	// int n=cols*rows;
	// SparseMatrix<double> AM(n, n);
	// MatrixXd rv(n,dims),bv(n,dims),xv(n,dims);
	// AM.reserve(VectorXd::Constant(n, 5));
	
	// for(int i=0;i<cols;i++){
	// 	for(int j=0;j<rows;j++){
	// 		int I=i+j*cols;
			
	// 		if(i==0||j==0||i==cols-1||j==rows-1){
	// 			AM.insert(I,I)=1;
	// 			for(int k=0;k<dims;k++){
	// 				bv(I,k)=AMatrix[i][j][k];
	// 			}
	// 		}else{
	// 			for(int k=0;k<dims;k++){
	// 				bv(I,k)=0;
	// 			}
	// 			AM.insert(I, I-1) = -1.0;
	// 			AM.insert(I, I+1) = -1.0;
	// 			AM.insert(I, I-cols) = -1.0;
	// 			AM.insert(I, I+cols) = -1.0;
	// 			AM.insert(I, I) = 4.0;
	// 		}
	// 	}
	// }
	// BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver; // faster
	// solver.preconditioner().setDroptol(0.01);
	// solver.setTolerance(1.0e-5);
	// solver.setMaxIterations(10);
	// rv=solver.compute(AM).solve(bv);
	
	// for(int i=0;i<cols;i++){
	// 	for(int j=0;j<rows;j++){
	// 		for(int k=0;k<dims;k++){
	// 			AMatrix[i][j][k]=static_cast<double>(rv(j*cols+i, k));
	// 		}
	// 	}
	// }
	// return 0;
	int m=(cols-2)*(rows-2);
	int q=2*cols+2*rows-4;
	VectorXd answerxv(m);
	VectorXd answeryv(m);
	VectorXd xv(q);
	VectorXd yv(q);
	std::cout<<m<<" "<<q<<std::endl;

	int p=0;
	for(int i=0,j=0;j<cols-1;j++){
		xv(p)=AMatrix[i][j][0];
		yv(p)=AMatrix[i][j][1];
		p++;
	}
	for(int i=0,j=cols-1;i<rows-1;i++){
		xv(p)=AMatrix[i][j][0];
		yv(p)=AMatrix[i][j][1];
		p++;
	}
	for(int i=rows-1,j=cols-1;j>0;j--){
		xv(p)=AMatrix[i][j][0];
		yv(p)=AMatrix[i][j][1];
		p++;
	}
	for(int i=rows-1,j=0;i>0;i--){
		xv(p)=AMatrix[i][j][0];
		yv(p)=AMatrix[i][j][1];
		p++;
	}

	answerxv=AInverse*xv;
	answeryv=AInverse*yv;
	
	for(int i=1;i<rows-1;i++){
		for(int j=1;j<cols-1;j++){
			int I=(i-1)*(cols-2)+j-1;
			AMatrix[i][j][0]=static_cast<double>(answerxv(I));
			AMatrix[i][j][1]=static_cast<double>(answeryv(I));	
		}
	}
	
	return 0;

}

int ShapeWarp::Refine(){
	for(int i=1;i<cols-1;i++){
		for(int j=1;j<rows-1;j++){
			double di=(double)i/(cols-1);
			double dj=(double)j/(rows-1);
			AMatrix[i][j][0]=	(1-dj)*AMatrix[i][0][0]+
								dj*AMatrix[i][rows-1][0]+
								(1-di)*AMatrix[0][j][0]+
								di*AMatrix[cols-1][j][0]-
								(1-di)*(1-dj)*AMatrix[0][0][0]-
								di*dj*AMatrix[cols-1][rows-1][0]-
								di*(1-dj)*AMatrix[cols-1][0][0]-
								(1-di)*dj*AMatrix[0][rows-1][0];
			
			AMatrix[i][j][1]=	(1-dj)*AMatrix[i][0][1]+
								dj*AMatrix[i][rows-1][1]+
								(1-di)*AMatrix[0][j][1]+
								di*AMatrix[cols-1][j][1]-
								(1-di)*(1-dj)*AMatrix[0][0][1]-
								di*dj*AMatrix[cols-1][rows-1][1]-
								di*(1-dj)*AMatrix[cols-1][0][1]-
								(1-di)*dj*AMatrix[0][rows-1][1];;
		}
	}
	return 0;
}

double ShapeWarp::GaussCE(double distance,double delta){
	double ans=(distance)/(delta)*(distance)/(delta);
	ans=exp(-ans);
	return ans;
}

int ShapeWarp::Liquidfy(Sample from, Sample to, int ox, int oy, double strongness){
	double dx=to.x-from.x;
	double dy=to.y-from.y;
	for(int i=1;i<cols-1;i++){
		for(int j=1;j<rows-1;j++){
			Sample temp(AMatrix[i][j][0],AMatrix[i][j][1]);
			double distance=std::sqrt((i-ox)*(i-ox)+(j-oy)*(j-oy));//index distance;
//			double distance=dis(temp, MiApex)/10;//real distance;
			double ce=GaussCE(distance, strongness);
			
			AMatrix[i][j][0]=AMatrix[i][j][0]+ce*dx;
			AMatrix[i][j][1]=AMatrix[i][j][1]+ce*dy;
			
		}
	}
	
	return 0;
}

double ShapeWarp::Cross(Sample v1, Sample v2, Sample u1, Sample u2) const {
	double vx=v2.x-v1.x;
	double vy=v2.y-v1.y;
	double ux=u2.x-u1.x;
	double uy=u2.y-u1.y;
	return vx*uy-ux*vy;
}

double ShapeWarp::dis(Sample s1, Sample s2) const {
	return std::sqrt((s1.x-s2.x)*(s1.x-s2.x)+(s1.y-s2.y)*(s1.y-s2.y));
}
double ShapeWarp::dis(Sample s1, Sample s2, Sample s) const {
	if(Cross(s, s1, s1, s2)*Cross(s, s2, s1, s2)<0){
		double h=dis(s1, s2);
		double area=std::abs((s.x*s1.y+s1.x*s2.y+s2.x*s.y)-
							 (s.y*s1.x+s1.y*s2.x+s2.y*s.x));
		return area/h;
	}else{
		return dis(s,s1)<dis(s,s2)?dis(s,s1):dis(s,s2);
	}
}

// choose a point: return true;
bool ShapeWarp::ChooseCurve(Sample s, Curve &bound){
	for(auto i=bound.begin();i!=bound.end();i++){
		if(dis(*i, s)<HandleRadius){
			choose=&(*i);
			return true;
		}
	}
	return false;
}

// insert a point: return true, need full curve;
bool ShapeWarp::InsertCurve(Sample s, Curve &bound, Sample end, std::vector<double> &sect){
	
	Curve tbound=bound;
	tbound.push_back(end);
	Curve curve=KCurve(tbound);
	
	auto c=sect.begin();
	for(auto b=tbound.begin();b!=tbound.end()-1;b++,c++){
		double td;
		Curve temp=CurveCut(curve, *(b+1), td);
		double t1=0;
		for(auto i=temp.begin();i!=temp.end()-1;i++){
			if(dis(*i,*(i+1),s)<HandleRadius){
				t1+=dis(*i,s);
				double s1=t1/td*(*c)+1e-5;
				double s2=(*c)-s1+1e-5;
				tbound.insert(b+1, s);
				*c=s2;
				sect.insert(c, s1);
				tbound.erase(tbound.end()-1);
				bound=tbound;
				return true;
			}else{
				t1+=dis(*i,*(i+1));
			}
		}
	}
	return false;
}
Curve ShapeWarp::RefineCurve(const Curve &bound, const std::vector<double> &sect) const {
	if(bound.size()-sect.size()!=1){
		exit(1);
	}
	Curve answer=KCurve(bound);
	answer=Sampler(answer, bound, sect);
	return answer;
}

// Cut a Curve to two Curves by the endpoint s or distance d;
Curve ShapeWarp::CurveCut(Curve &c, const Sample s, double &totaldis) const {
	Curve answer,surplus;
	totaldis=0;
	answer.clear();
	surplus.clear();
	auto iter=c.begin();
	if(dis(*c.begin(), s)<2*HandleRadius){
		answer.push_back(*c.begin());
		answer.push_back(s);
		c.insert(c.begin(), s);
		return answer;
	}
	for(iter=c.begin();iter!=c.end()-1;iter++){
		if(dis(*iter, s)<1e-5){
			answer.push_back(*iter);
			break;
		}else if(dis(*(iter+1),s)<1e-5){
			totaldis+=dis(*iter, *(iter+1));
			answer.push_back(*iter);
			continue;
		}else if(dis(*iter, *(iter+1), s)<HandleRadius){// end point s is on this segments;
			answer.push_back(*iter);
			answer.push_back(s);
			surplus.push_back(s);
			totaldis+=dis(*iter, s);
			iter++;
			break;
		}else{
			totaldis+=dis(*iter, *(iter+1));
			answer.push_back(*iter);
		}
	}
	
	for(;iter!=c.end();iter++){
		surplus.push_back(*iter);
	}
	c.clear();
	c=surplus;
	return answer;
}
Curve ShapeWarp::CurveCut(Curve &c, double d) const {
	double totaldis=0,length=1;
	Curve answer,surplus;
	answer.clear();
	surplus.clear();
	auto iter=c.begin();
	while(iter!=c.end()-1){
		length=dis(*iter, *(iter+1));
		if(totaldis+length>d){
			answer.push_back(*iter);
			double partdis=d-totaldis;
			double rate=partdis/length;
			double tx=(1-rate)*iter->x+rate*(iter+1)->x;
			double ty=(1-rate)*iter->y+rate*(iter+1)->y;
			Sample temp(tx, ty);
			answer.push_back(temp);
			surplus.push_back(temp);
			totaldis+=dis(temp, *iter);
			iter++;
			break;
		}else if(totaldis+length==d){
			answer.push_back(*iter++);
			answer.push_back(*iter);
			totaldis+=length;
			break;
		}else{
			answer.push_back(*iter);
			iter++;
			totaldis+=length;
		}
	}
	for(;iter!=c.end();iter++){
		surplus.push_back(*iter);
	}
	c=surplus;
	return answer;
}

// Sampler a Curve with Sects points as below;
// |<----11--->|<---8--->|<-----14----->|  (length)
// |-----------|---------|--------------|
// |<---5.5--->|<---2--->|<----3.5----->|  (sects)
// the Sampler answer is:
// |*--*--*--*--*--*-|--*----*--|--*----*----*----*| (answer);
// sum of '-' equal to length & sum of '*' equal to sects+1;
// b: whole bound; c: whole curve; sects: each segments' number;
Curve ShapeWarp::Sampler(Curve c, const Curve &b, const std::vector<double> &sects) const {
	if(b.size()-sects.size()!=1){
		exit(1);
	}
	
	Curve answer;
	answer.clear();
	double residue=0;
	for(int i=0;i<b.size()-1;i++){
		double td;
		Curve tc=CurveCut(c, b[i+1], td);
		double sect=sects[i]-residue;
		int segs=(int)sect;
		double next=sect-segs;
		double stepdis=td/sects[i];
		
		Curve ttc=CurveCut(tc, residue*stepdis);
		residue=next<1e-5?0:1-next;
		
		for(int k=0;k<=segs;k++){
			answer.push_back(*tc.begin());
			CurveCut(tc, stepdis);
		}
	}
	return answer;
}

int ShapeWarp::Choose(double x, double y){
	FE=NO;
	Sample cho(x,y);
	choose=NULL;
	if(dis(cho, CMApex)<HandleRadius){
		CK=CENTRE;
		choose=&CMApex;
		return 0;
	}
	
	if(ChooseCurve(cho, BBound)){
		if(FE==NO){
			FE=BOTTOM;
		}else{
			FE=ALL;
			return  0;
		}
		CK=EXIST;
	}
	if(ChooseCurve(cho, RBound)){
		if(FE==NO){
			FE=RIGHT;
		}else{
			FE=ALL;
			return 0;
		}
		CK=EXIST;
	}
	if(ChooseCurve(cho, TBound)){
		if(FE==NO){
			FE=TOP;
		}else{
			FE=ALL;
			return 0;
		}
		CK=EXIST;
	}
	if(ChooseCurve(cho, LBound)){
		if(FE==NO){
			FE=LEFT;
		}else{
			FE=ALL;
			return 0;
		}
		CK=EXIST;
	}

	if(CK==EXIST)
		return 0;

	if(InsertCurve(cho, BBound, *RBound.begin(), BSect)){
		ChooseCurve(cho, BBound);
		FE=BOTTOM;
		FE=ALL;
		CK=EXIST;
		return 0;
	}
	if(InsertCurve(cho, RBound, *TBound.begin(), RSect)){
		ChooseCurve(cho, RBound);
		FE=RIGHT;
		FE=ALL;
		CK=EXIST;
		return 0;
	}
	if(InsertCurve(cho, TBound, *LBound.begin(), TSect)){
		ChooseCurve(cho, TBound);
		FE=TOP;
		FE=ALL;
		CK=EXIST;
		return 0;
	}
	if(InsertCurve(cho, LBound, *BBound.begin(), LSect)){
		ChooseCurve(cho, LBound);
		FE=LEFT;
		FE=ALL;
		CK=EXIST;
		return 0;
	}
	
	choose=NULL;
	CK=NONE;
	return 1;
}
int ShapeWarp::Set(double x, double y){
	Sample temp(x, y);
	if(dis(temp, MiApex)<2*HandleRadius){
		x=MiApex.x;
		y=MiApex.y;
	}
	for(auto i=BBound.begin();i!=BBound.end();i++){
		if(choose==&*i){
			continue;
		}
		if(dis(temp,*i)<2*HandleRadius){
			x=i->x;
			y=i->y;
		}
	}
	for(auto i=RBound.begin();i!=RBound.end();i++){
		if(choose==&*i){
			continue;
		}
		if(dis(temp,*i)<2*HandleRadius){
			x=i->x;
			y=i->y;
		}
	}
	for(auto i=TBound.begin();i!=TBound.end();i++){
		if(choose==&*i){
			continue;
		}
		if(dis(temp,*i)<2*HandleRadius){
			x=i->x;
			y=i->y;
		}
	}
	for(auto i=LBound.begin();i!=LBound.end();i++){
		if(choose==&*i){
			continue;
		}
		if(dis(temp,*i)<2*HandleRadius){
			x=i->x;
			y=i->y;
		}
	}
	switch (CK) {
		case CENTRE:
			xmove=x-MiApex.x;
			ymove=y-MiApex.y;
			CMApex.x=x;
			CMApex.y=y;
			Flush();
			return 0;
		case EXIST:
			choose->x=x;
			choose->y=y;
			Flush();
			return 0;
		case NONE:
			return 0;
		default:
			break;
	}
	return 1;
}
int ShapeWarp::Delete(double x, double y){
	Sample cho(x,y);
	auto bsi=BSect.begin();
	for(auto iter=BBound.begin();iter!=BBound.end();iter++,bsi++){
		if(dis(*iter, cho)<HandleRadius){
			if(iter==BBound.begin()){
				std::cout<<"can't delete a vertex"<<std::endl;
				return 0;
			}
			BBound.erase(iter);
			double newsegs=*bsi+*(bsi-1);
			*(bsi-1)=newsegs;
			BSect.erase(bsi);
			Flush();
			return 0;
		}
	}
	auto rsi=RSect.begin();
	for(auto iter=RBound.begin();iter!=RBound.end();iter++,rsi++){
		if(dis(*iter, cho)<HandleRadius){
			RBound.erase(iter);
			double newsegs=*rsi+*(rsi-1);
			*(rsi-1)=newsegs;
			RSect.erase(rsi);
			Flush();
			return 0;
		}
	}
	auto tsi=TSect.begin();
	for(auto iter=TBound.begin();iter!=TBound.end();iter++,tsi++){
		if(dis(*iter, cho)<HandleRadius){
			TBound.erase(iter);
			double newsegs=*tsi+*(tsi-1);
			*(tsi-1)=newsegs;
			TSect.erase(tsi);
			Flush();
			return 0;
		}
	}
	auto lsi=LSect.begin();
	for(auto iter=LBound.begin();iter!=LBound.end();iter++,lsi++){
		if(dis(*iter, cho)<HandleRadius){
			LBound.erase(iter);
			double newsegs=*lsi+*(lsi-1);
			*(lsi-1)=newsegs;
			LSect.erase(lsi);
			Flush();
			return 0;
		}
	}
	return 1;
}
int ShapeWarp::Flush(){
	if(Bind()) return 1;

	switch (AK) {
		case POISSONIMAGEEDITING:
			PoissonRefine();
			break;
		case TRANSFINITEINTERPOLATION:
			Refine();
			break;
		default:
			return 1;
	}
	MiApex.x=AMatrix[cols/2][rows/2][0];
	MiApex.y=AMatrix[cols/2][rows/2][1];
	CMApex.x=MiApex.x+xmove;
	CMApex.y=MiApex.y+ymove;
	if(Liquidfy(MiApex, CMApex, cols/2, rows/2 ,delta)) return 1;
	// if(AK==POISSONIMAGEEDITING)
	// 	Smooth();
	return 0;
}

Curve ShapeWarp::ReadControlPoints(){
	Curve answer;
	for(auto i=BBound.begin();i!=BBound.end();i++){
		answer.push_back(*i);
	}
	for(auto i=RBound.begin();i!=RBound.end();i++){
		answer.push_back(*i);
	}
	for(auto i=TBound.begin();i!=TBound.end();i++){
		answer.push_back(*i);
	}
	for(auto i=LBound.begin();i!=LBound.end();i++){
		answer.push_back(*i);
	}
	answer.push_back(CMApex);
	return answer;
}

int ShapeWarp::Smooth(){
	int n=cols*rows;
	int d=sqrt(n)/4;
//	d=d-delta/d;
	SparseMatrix<double> AM(n, n);
	MatrixXd rv(n,dims),bv(n,dims),xv(n,dims);
	AM.reserve(VectorXd::Constant(n, 5));
	
	for(int i=0;i<cols;i++){
		for(int j=0;j<rows;j++){
			int I=i+j*cols;
			
			if(i==0||j==0||i==cols-1||j==rows-1){
				AM.insert(I,I)=1;
				for(int k=0;k<dims;k++){
					bv(I,k)=AMatrix[i][j][k];
				}
			}else if(i>d&&j>d&&i<cols-1-d&&j<rows-1-d){
				AM.insert(I,I)=1;
				for(int k=0;k<dims;k++){
					bv(I,k)=AMatrix[i][j][k];
				}
			}else{
				for(int k=0;k<dims;k++){
					bv(I,k)=0;
				}
				AM.insert(I, I-1) = -1.0;
				AM.insert(I, I+1) = -1.0;
				AM.insert(I, I-cols) = -1.0;
				AM.insert(I, I+cols) = -1.0;
				AM.insert(I, I) = 4.0;
			}
		}
	}
	BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver; // faster
	solver.preconditioner().setDroptol(0.01);
	solver.setTolerance(1.0e-5);
	solver.setMaxIterations(10);
	rv=solver.compute(AM).solve(bv);
	
	for(int i=0;i<cols;i++){
		for(int j=0;j<rows;j++){
			for(int k=0;k<dims;k++){
				AMatrix[i][j][k]=static_cast<double>(rv(j*cols+i, k));
			}
		}
	}
	return 0;
}

int ShapeWarp::UnFlip(){
	
	return 0;
}
