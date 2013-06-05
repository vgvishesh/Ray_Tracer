#include <iostream>
#include<gl/glut.h>


void init()
{
	 glMatrixMode(GL_PROJECTION);
	 glLoadIdentity();
	 glOrtho(0,10,-10,10,1,50);
	 glMatrixMode(GL_MODELVIEW);
	 glLoadIdentity();
	 gluLookAt(0,0,10,0,0,0,0,1,0);
	 glViewport(0,0,50,50);
}

void display()
{
	glClearColor(0,0,0,0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glColor3f(1,0,1);
	glutSolidSphere(5,50,50);
	glFlush();
}

int main(int agrc,char **agrv)
{
	glutInit(&agrc,agrv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB); //for animation programs GLUT_DOUBLE , GLUT_DEPTH for 3d apps
	glutInitWindowSize(500,500);			//screen resolution (width,height)
	glutCreateWindow("3D bezeir curves");
	glutInitWindowPosition(50,25);
	glutDisplayFunc(display);		//display function
	init();							//put all initialization commands in this funtion
	glutMainLoop();

	return 0;
}