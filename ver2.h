#ifndef VER2_H
#define VER2_H

#include<iostream>
#include<math.h>
using namespace std;

class rtCamera;
class rtRay;
class Sphere;
class rtPoint
{
	friend class rtCamera;
	friend class rtRay;
	friend class Sphere;
private:
	double xx,yy,zz;
public:
	rtPoint():xx(0),yy(0),zz(0){};
	void rtPointSet(double x,double y,double z)
	{
		xx=x;yy=y;zz=z;
	}
	void rtPointSet(rtPoint p)
	{
		xx=p.xx;
		yy=p.yy;
		zz=p.zz;
	}
};

class rtRay
{
	friend class Sphere;
	friend class rtCamera;
private:
	rtPoint source,destination;//,dir;
	bool hit;
public:
	rtRay():source(),destination(){};
	void rtRaySetSource(rtPoint s)
	{
		source.rtPointSet(s);
	}
	/*void rtRaySetDir()
	{
		dir.rtPointSet(destination.xx-source.xx,destination.yy-source.yy,destination.zz-source.zz);
	}*/
	void rtRaySetDestination(rtPoint d)
	{
		destination.rtPointSet(d);
	}
};

class Sphere
{
private:
	rtPoint center;
	double radius;
public:
	Sphere():center(),radius(1){};
	Sphere(double x,double y,double z,double r)
	{
		center.rtPointSet(x,y,z);
		radius=r;
	}
	bool Intersect(rtRay R)
	{
		double dx,dy,dz,X,Y,Z,A,B,C,t1,t2,D;
		dx=R.destination.xx-R.source.xx;
		dy=R.destination.yy-R.source.yy;
		dz=R.destination.zz-R.source.zz;

		X=R.source.xx-center.xx;
		Y=R.source.yy-center.yy;
		Z=R.source.zz-center.zz;

		A= pow(dx,2)+pow(dy,2)+pow(dz,2);
		B= dx*X + dy*Y + dz*Z;
		C= pow(X,2) + pow(Y,2) + pow(Z,2) - pow(radius,2);
		/*cout<<"\n\n"<<R.destination.xx<<" "<<R.destination.yy<<" "<<R.destination.zz;
		cout<<"\n"<<B<<"\t"<<A<<"\t"<<C<<"\t";*/

		D=pow(B,2) - A*C;
		//cout<<"\n"<<D<<"\t";
		if(D>=0)
		{
			t1= (-B + sqrt(D))/A;
			t2= (-B - sqrt(D))/A;
			//cout<<t1<<"\t"<<t2;
			if(t1>=1 || t2>=1)
			 return 1;
		}
		return 0;
	}
};

class rtVec
{
	friend class rtCamera;
private:
	double x,y,z;
public:
	void rtVecSet(double dx,double dy,double dz)
	{
		x=dx;
		y=dy;
		z=dz;
	}
	void rtVecSet(rtVec in)
	{
		x=in.x;
		y=in.y;
		z=in.z;
	}
	rtVec cross(rtVec v)
	{
		rtVec c;
		c.x = y*v.z + (-1)*z*v.y;
		c.y = z*v.x + (-1)*x*v.z;
		c.z = x*v.y + (-1)*y*v.x;
		return (c);
	}
	double dot(rtVec v)
	{
		rtVec d;
		double sum;
		d.x=x*v.x;
		d.y=y*v.y;
		d.z=z*v.z;
		sum = d.x+d.y+d.z;
		return (sum);
	}
	rtVec set_diff(rtVec v)
	{
		rtVec d;
		d.x=x-v.x;
		d.y=y-v.y;
		d.z=z-d.z;
		return (d);
	}
	void flip()
	{
		x=-x;
		y=-y;
		z=-z;
	}
	void normalize()
	{
		double dem;
		double sum;
		sum=x*x+y*y+z*z;
		dem=sqrt(sum);
		x=x/dem;
		y=y/dem;
		z=z/dem;
	}
	rtVec(double xx,double yy,double zz)
	{
		x=xx;
		y=yy;
		z=zz;
	}
	rtVec(rtVec &v)
	{
		x=v.x;
		y=v.y;
		z=v.z;
	}
	//point rtVec_to_point();
	void scalar_mul(double c)
	{
		x=x*c;
		y=y*c;
		z=z*c;
	}
	rtVec()
	{
		x=y=z=0;
	}
};

class rtCamera
{
private:
	rtPoint eye, lookat, up;
	rtVec n,u,v;
	double near,far,left,right,top,bottom;
	double width,height;
	rtRay **rays;
public:
	void rtCameraSet(double a=0,double b=0,double c=10,double d=0,double e=0,double f=0,double g=0,double h=1,double i=0)
	{
		eye.rtPointSet(a,b,c);
		lookat.rtPointSet(d,e,f);
		up.rtPointSet(g,h,i);
		rtVec UP(up.xx,up.yy,up.zz);
		n.rtVecSet(eye.xx-lookat.xx,eye.yy-lookat.yy,eye.zz-lookat.zz);
		n.normalize();
		u.rtVecSet(UP.cross(n));
		//u.normalize();
		v.rtVecSet(n.cross(u));		
	}
	void rtCameraWindow(double a,double b,double c,double d,double n,double f)
	{
		left=a;
		right=b;
		bottom=c;
		top=d;
		near=n;				//on z axis
		far=f;				//on z axis
		/*glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D(left,right,bottom,top);*/
	}
	void rtCameraSetRays()
	{
		width=right-left;
		height=top-bottom;
		rays=new rtRay*[(int)height];
		rtPoint p;
		for(int i=0;i<height;i++)
			rays[i]=new rtRay[(int)width];
		//cout<<bottom<<" "<<height;
		for(int i=0,b=bottom;i<height;i++,b++)
			for(int j=0,l=left;j<width;j++,l++)
			{
				rays[i][j].rtRaySetSource(eye);
				p.rtPointSet(eye.xx -near*n.x + u.x*l +v.x*b, 
							 eye.yy -near*n.y + u.y*l +v.y*b, 
							 eye.zz -near*n.z + u.z*l +v.z*b);					//z axis point of the ray...
				rays[i][j].rtRaySetDestination(p);
			}
	}		
	void rtCameraCastRay(Sphere S)
	{
		for(int i=0;i<height;i++)
			for(int j=0;j<width;j++)
				rays[i][j].hit=S.Intersect(rays[i][j]);
	}
	bool rtvalue(int i,int j)
	{
		return rays[i][j].hit;
	}
};
#endif