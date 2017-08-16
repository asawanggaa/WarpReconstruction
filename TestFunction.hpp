//
//  TestFunction.hpp
//  LinearWarp
//
//  Created by t_wangju on 07/11/2016.
//  Copyright © 2016 t_wangju. All rights reserved.
//

#ifndef TestFunction_hpp
#define TestFunction_hpp

#include <stdio.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/gl.h>

#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "ShapeWarp.h"

#define HandleRadius 5
void DrawHandle(int x, int y, int Radius){
	double theta;
	
	glBegin(GL_LINE_LOOP);
	{
		for (int i=0; i<36; i++) {
			theta=i*2*M_PI/36;
			double dx=Radius*cos(theta)+x;
			double dy=Radius*sin(theta)+y;
			glVertex2d(dx, dy);
		}
	}
	glEnd();
}
void DrawTex(ShapeWarp &psw, GLuint tex){
	glDisable(GL_BLEND);
	glReadBuffer(GL_FRONT);
	glDrawBuffer(GL_BACK);
	
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(1, 1, 1, 1);
	glColor3f(1.0, 1.0, 1.0);
	glBindTexture(GL_TEXTURE_2D, tex);
	glEnable(GL_TEXTURE_2D);
	glBegin(GL_QUADS);
	{
		for(int i=0;i<psw.Cols()-1;i++){
			for(int j=0;j<psw.Rows()-1;j++){
				glTexCoord2f(i/(double)(psw.Cols()-1), j/(double)(psw.Rows()-1));
				glVertex2i(psw.ReadAMatrix()[i][j][0], psw.ReadAMatrix()[i][j][1]);
				glTexCoord2f((i+1)/(double)(psw.Cols()-1), j/(double)(psw.Rows()-1));
				glVertex2i(psw.ReadAMatrix()[i+1][j][0], psw.ReadAMatrix()[i+1][j][1]);
				glTexCoord2f((i+1)/(double)(psw.Cols()-1), (j+1)/(double)(psw.Rows()-1));
				glVertex2i(psw.ReadAMatrix()[i+1][j+1][0], psw.ReadAMatrix()[i+1][j+1][1]);
				glTexCoord2f(i/(double)(psw.Cols()-1), (j+1)/(double)(psw.Rows()-1));
				glVertex2i(psw.ReadAMatrix()[i][j+1][0], psw.ReadAMatrix()[i][j+1][1]);
				
			}
		}
	}
	glEnd();
	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);

}
void DrawLine(std::vector<Sample> input){
	glBegin(GL_LINE_STRIP);
	{
		for(auto i=input.begin();i!=input.end();i++){
			glVertex2d(i->x, i->y);
		}
	}
	glEnd();
}
void DrawConstrain(ShapeWarp &psw){
	glColor3f(1.0f, 0.0f, 0.0f);
	Curve answer=psw.ReadControlPoints();
	for(auto i=answer.begin();i!=answer.end();i++){
		DrawHandle(i->x, i->y, HandleRadius);
	}
	glColor3f(0.0f, 1.0f, 0.0f);
	DrawHandle(psw.MIApex().x, psw.MIApex().y, HandleRadius);
}
void DrawShapeWarp(ShapeWarp &psw){
	glColor3f(1.0f, 0.0f, 0.0f);
	for(int i=0;i<psw.Cols();i++){
		for(int j=0;j<psw.Rows();j++){
			DrawHandle(psw.ReadAMatrix()[i][j][0], psw.ReadAMatrix()[i][j][1], 5);
		}
	}
}
#endif /* TestFunction_hpp */
