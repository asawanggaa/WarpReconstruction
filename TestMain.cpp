/* Using standard C++ output libraries */
#include <stdio.h>
#include <cstdlib>
#include <iostream>

#include <GL/glew.h>
#include <GL/glut.h>
#include <eigen3/Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <cmath>
#include <vector>
#include <set>
#include <map>


#include "PredictedPath.h"
#include "ShapeWarp.h"
#include "TestFunction.hpp"

// using namespace PredictedPath;
using namespace Eigen;
using namespace std;

const int cols=64;
const int rows=64;

ShapeWarp PSW(100,100,700,700,cols,rows,2);

std::vector<Sample> example;

std::vector<Sample> Example(int k, int mod){
	std::vector<Sample> answer;
	for(int i=0;i<k;i++){
		Sample temp(random()%mod,random()*i%mod);
		answer.push_back(temp);
	}
	return answer;
}
GLuint tex=0;
GLubyte data[800][800][4];

void texinit(const char * path){
	printf("TexInit\n");
	glShadeModel(GL_SMOOTH);
	glEnable(GL_COLOR_MATERIAL);
	glMatrixMode(GL_PROJECTION);
	
	cv::Mat img = cv::imread(path);
	cv::flip(img, img, 0);

	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
	glPixelStorei(GL_UNPACK_ROW_LENGTH, img.step/img.elemSize());

	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &tex);
	glBindTexture(GL_TEXTURE_2D, tex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img.cols, img.rows, 0, GL_BGR, GL_UNSIGNED_BYTE, img.ptr());
	
	glBindTexture(GL_TEXTURE_2D, 0);
}

void init(int argc, char * argv[], int x, int y, int width, int height){
	printf("Init\n");
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA);
	glutInitWindowPosition(x, y);
	glutInitWindowSize(width, height);
	glutCreateWindow("LinearWarp");
	gluOrtho2D(x, x+width, y, y+height);
	
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glReadBuffer(GL_FRONT);
	glDrawBuffer(GL_BACK);
	int xmin=100;
	int xmax=width-100;
	int ymin=100;
	int ymax=height-100;
	
	glFlush();
	glutSwapBuffers();
}

void DisplayFunc(void){
	printf("DisplayFunc\n");
	glDrawBuffer(GL_FRONT_AND_BACK);
	glClearColor(1, 1, 1, 1);
	glClear(GL_COLOR_BUFFER_BIT);
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, tex);
	glBegin(GL_QUADS);
	{
		glTexCoord2f(0.0f, 0.0f);    glVertex2i(100, 100);
		glTexCoord2f(1.0f, 0.0f);    glVertex2i(700, 100);
		glTexCoord2f(1.0f, 1.0f);    glVertex2i(700, 700);
		glTexCoord2f(0.0f, 1.0f);    glVertex2i(100, 700);
	}
	glEnd();
	glFlush();
	glutSwapBuffers();
	glBindTexture(GL_TEXTURE_2D, 0);
}
void MouseClickFunc(int button, int state, int x, int y){
	if(button==GLUT_LEFT_BUTTON&&state==GLUT_DOWN){
		glClear(GL_COLOR_BUFFER_BIT);
		PSW.Choose(x, 800-y);
		PSW.Flush();
		DrawTex(PSW, tex);
		DrawConstrain(PSW);
		glFlush();
		glutSwapBuffers();
	}
	if(button==GLUT_LEFT_BUTTON&&state==GLUT_UP){
		DrawTex(PSW, tex);
		DrawConstrain(PSW);
		glFlush();
		glutSwapBuffers();
	}
}
void MouseMoveFunc(int x, int y){
	glClear(GL_COLOR_BUFFER_BIT);
	PSW.Set(x, 800-y);
	PSW.Flush();
	glColor3f(0.0f, 0.0f, 0.0f);
	for(int i=0;i<PSW.Cols();i++){
		for(int j=0;j<PSW.Rows();j++){
			DrawHandle(PSW.ReadAMatrix()[i][j][0], PSW.ReadAMatrix()[i][j][1], HandleRadius);
		}
	}
	DrawTex(PSW, tex);
	DrawConstrain(PSW);
	glFlush();
	glutSwapBuffers();
	
}
void KeyboardFunc(unsigned char key, int x, int y){
	y=800-y;
	switch(key){
		case 't':
		case 'T':
			PSW.ChangeAlgorithm();
			PSW.Flush();
			break;
		case 's':
		case 'S':
			break;
		case 'd':
		case 'D':
			PSW.Delete(x, y);
			PSW.Flush();
			break;
		case 'q':
		case 'Q':
			exit(0);
		case 'k':
		case 'K':
			glClear(GL_COLOR_BUFFER_BIT);
			glColor3f(0.0f, 0.0f, 0.0f);
			for(int i=0;i<PSW.Cols();i++){
				for(int j=0;j<PSW.Rows();j++){
					DrawHandle(PSW.ReadAMatrix()[i][j][0], PSW.ReadAMatrix()[i][j][1], HandleRadius);
				}
			}
			DrawConstrain(PSW);
			glFlush();
			glutSwapBuffers();
			return;
		case '+':
			PSW.SetDelta(PSW.GetDelta()*1.1);
			PSW.Flush();
			break;
		case '-':
			PSW.SetDelta(PSW.GetDelta()/1.1);
			PSW.Flush();
			break;
		default:
			return;
	}
	DrawTex(PSW, tex);
	DrawConstrain(PSW);
	glFlush();
	glutSwapBuffers();
}

int main(int argc, char* argv[]) {
	// insert code here...
	printf("Program begin\n");
	init(argc, argv, 0, 0, 800, 800);
	texinit("/home/asawang/Desktop/dancing.jpg");
	glutDisplayFunc(DisplayFunc);
	glutMouseFunc(MouseClickFunc);
	glutMotionFunc(MouseMoveFunc);
	glutKeyboardFunc(KeyboardFunc);
	glutMainLoop();
	return 0;
}