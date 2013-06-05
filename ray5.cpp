#include <iostream>
#include<gl/glut.h>
#include "ver5.h"
#define L -500
#define R 500
#define B -500
#define T 500
using namespace std; 

rtScene SCENE("scene5.txt");
double l,r,b,t;
void display();

void init()
{
	 glMatrixMode(GL_PROJECTION);
	 glLoadIdentity();
	 gluOrtho2D(L,R,B,T);
	 SCENE.SceneRender();	 
}

void display()
{
	int count=0;
	glClearColor(0,0,0,0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glColor3f(1,0,1);
	glPointSize(1);
	SCENE.GetSceneDimentions(l,r,b,t);
	glBegin(GL_POINTS);
		for(int i=b;i<t;i++)
			for(int j=l;j<r;j++)
			{
				SCENE.SceneRayValue(i+t,j+r);
				glVertex2d(i,j);
			}
	glEnd();
	glFlush();
}
int main(int agrc,char **agrv)
{
	
	glutInit(&agrc,agrv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB); //for animation programs GLUT_DOUBLE , GLUT_DEPTH for 3d apps
	glutInitWindowSize(500,500);			//screen resolution (width,height)
	glutCreateWindow("Ray Tracer");
	glutInitWindowPosition(50,25);
	glutDisplayFunc(display);		//display function
	init();							//put all initialization commands in this funtion
	glutMainLoop();

	return 0;
}