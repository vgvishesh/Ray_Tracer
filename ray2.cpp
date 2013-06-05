#include <iostream>
#include<gl/glut.h>
#include "ver2.h"
#define L -250
#define R 250
#define B -250
#define T 250
#define F 50
//#define Z 6	// z coordinate for eye...
using namespace std; 

rtCamera CAM;
double Z=10,N=1,l,b,r,t;
void display();

void Trace(rtCamera &CAM)
{
	Sphere S(0,0,0,60);
	CAM.rtCameraCastRay(S);
}

void MakeMapping()
{
	// dimentions for near plane
	b=l=-100;
	t=r=100;
}


void init(rtCamera &CAM)
{
	 glMatrixMode(GL_PROJECTION);
	 glLoadIdentity();
	 //glOrtho(0,10,-10,10,1,50);
	 gluOrtho2D(L,R,B,T);
	 MakeMapping();
	 CAM.rtCameraWindow(l,r,b,t,N,F);
	 CAM.rtCameraSet(10,10,Z,0,0,0,0,1,0);
	 CAM.rtCameraSetRays();
	 Trace(CAM);
	 /*glMatrixMode(GL_MODELVIEW);
	 glLoadIdentity();
	 gluLookAt(0,0,10,0,0,0,0,1,0);*/
	 //glViewport(0,0,50,50);
}

void specialKey(int key,int mx,int my)
{
	switch(key)
	{
		case GLUT_KEY_UP:Z++;
						init(::CAM);
						cout<<"\nz:"<<Z;
						break;
		case GLUT_KEY_DOWN:Z--;
						init(::CAM);
						cout<<"\nz:"<<Z;
						break;
		case GLUT_KEY_RIGHT:N++;
						init(::CAM);
						cout<<"\nn:"<<N;
						break;
		case GLUT_KEY_LEFT:N--;
						init(::CAM);
						cout<<"\nn:"<<N;
						break;
	}
	display();
}

void display()
{
	int count=0;
	glClearColor(0,0,0,0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glColor3f(1,0,1);
	glPointSize(1);
	glBegin(GL_POINTS);
		for(int i=b;i<t;i++)
			for(int j=l;j<r;j++)
				if(::CAM.rtvalue(i+t,j+r))
				{
					glColor3f(1,0,1);
					glVertex2d(i,j);
					count++;
				}
				else
				{
					glColor3f(0,0,0);
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
	glutSpecialFunc(specialKey);
	glutInitWindowPosition(50,25);
	glutDisplayFunc(display);		//display function
	init(::CAM);							//put all initialization commands in this funtion
	glutMainLoop();

	return 0;
}