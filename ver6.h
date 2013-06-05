#ifndef VER4_H
#define VER4_H

#include<iostream>
#include<math.h>
#include<fstream>
#include<string>
#include<gl/glut.h>
using namespace std;

class rtCamera;
class rtRay;
class rtLight;
class Sphere;
class rtVec;

class rtPoint
{
	friend class rtCamera;
	friend class rtRay;
	friend class Sphere;
	friend class rtScene;
private:
	double xx,yy,zz;
public:
	rtPoint():xx(0),yy(0),zz(0){/*cout<<"vi";*/};
	void rtPointSet(double x,double y,double z)
	{
		xx=x;yy=y;zz=z;
		//cout<<xx<<" "<<yy<<" "<<zz<<"\n";
	}
	void rtPointSet(rtPoint p)
	{
		xx=p.xx;
		yy=p.yy;
		zz=p.zz;
	}
	void rtPointGetPoint(double &x,double &y,double &z)
	{
		x=this->xx;
		y=this->yy;
		z=this->zz;
	}
};

class rtVec
{
	friend class rtCamera;
	friend class Sphere;
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
	void rtVecSet(rtPoint P)
	{
		P.rtPointGetPoint(this->x,this->y,this->z);
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
		d.z=z-v.z;
		return (d);
	}
	rtVec add(rtVec v)
	{
		rtVec d;
		d.x=x+v.x;
		d.y=y+v.y;
		d.z=z+v.z;
		return (d);
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
	void rtVecGetVec(double &x,double &y,double &z)
	{
		x=this->x;
		y=this->y;
		z=this->z;
	}
};

class color3
{
	friend class rtCamera;
	friend class Sphere;
private:
	double r,g,b;
public:
	color3()
	{
		r=0.0;
		g=0.0;
		b=0.0;
	}
	void set(double rr,double gg,double bb)
	{
		r=rr;
		g=gg;
		b=bb;
	}
	void set(color3 c)
	{
		r=c.r;
		g=c.g;
		b=c.b;
	}
	//void set_color_runtime();
	void getColor(double &r,double &g,double &b)
	{
		r=this->r;
		g=this->g;
		b=this->b;
	}
	void coloradd(int con,color3 c)
	{
		r=con*c.r;
		g=con*c.g;
		b=con*c.b;
	}
};


class rtRay
{
	friend class Sphere;
	friend class rtCamera;
private:
	rtPoint source,destination,hitPoint;//,dir;
	rtVec dir;
	color3 rayColor;
	double *lambert,*phong;
	int nLP;	// just an interger count of no. of values in arrays of lambert and phong , which will be equal to no. of light sources in the scene
	bool hit;
public:
	rtRay():source(),destination(),rayColor(){};
	void rtRaySetSource(rtPoint s)
	{
		source.rtPointSet(s);
	}
	void rtRayGetSource(double &x,double &y,double &z)
	{
		source.rtPointGetPoint(x,y,z);
	}
	void rtRaySetDir()
	{
		double dx,dy,dz;
		dx=destination.xx-source.xx;
		dy=destination.yy-source.yy;
		dz=destination.zz-source.zz;
		dir.rtVecSet(dx,dy,dz);
	}
	void rtRayGetDir(double &x,double &y,double &z)
	{
		dir.rtVecGetVec(x,y,z);
	}
	void rtRaySetDestination(rtPoint d)
	{
		destination.rtPointSet(d);
	}
	void rtRayGetDestination(double &x,double &y,double &z)
	{
		destination.rtPointGetPoint(x,y,z);
	}		
	void rtRaySetRayColor(color3 C)
	{
		rayColor.set(C);
	}
	void rtRaySetHitPoint(double t)
	{
		double x,y,z;
		x=source.xx*(1-t)+destination.xx*t;
		y=source.yy*(1-t)+destination.yy*t;
		z=source.zz*(1-t)+destination.zz*t;
		hitPoint.rtPointSet(x,y,z);
	}	
	void rtRayGetHitPoint(double &x,double &y,double &z)
	{
		hitPoint.rtPointGetPoint(x,y,z);
	}
	void rtRaySetCoeffSize(int nl)
	{
		nLP=nl;
		lambert=new double[nLP];
		phong=new double[nLP];
	}
	void rtDirNormalize()
	{
		dir.normalize();
	}
	void rtRaySetLambert(rtVec s,rtVec m,int i)
	{
		double dp=s.dot(m);
		lambert[i]=dp>0?dp:0;
	}
	void rtRaySetPhong(rtVec s,rtVec m,int i)
	{
		rtVec v(dir);
		v.scalar_mul(-1);
		rtVec h(v.add(s));
		h.normalize();
		double dp=h.dot(m);
		phong[i]=dp>0?dp:0;
	}
	double rtRayGetLambert(int i)
	{
		return lambert[i];
	}
	double rtRayGetPhong(int i)
	{
		return phong[i];
	}
};

class Material
{
	friend class Sphere;
private:
	double diff_r,diff_g,diff_b;
	double spec_r,spec_g,spec_b,phong_exp;		//1<=phong_exp<=200
	double amb_r,amb_g,amb_b;
public:
	Material()
	{
		//defalut material is black plastic
		amb_r=0;amb_g=0;amb_b=0;
		diff_r=0.01;diff_g=0.01;diff_b=0.01;
		spec_r=0.50;spec_g=0.50;spec_b=0.50;
		phong_exp=32;
		//cout<<"set";
	}
	void MaterialSet(double Par,double Pag,double Pab,double Pdr,double Pdg,double Pdb,double Psr,double Psg,double Psb,double PhCo)
	{
		amb_r=Par;amb_g=Pag;amb_b=Pab;
		diff_r=Pdr;diff_g=Pdg;diff_b=Pdb;
		spec_r=Psr;spec_g=Psg;spec_b=Psb;
		phong_exp=PhCo;
	}
	void MaterialGetAmbCoeff(double &Par,double &Pag,double &Pab)
	{
		Par=amb_r;
		Pag=amb_g;
		Pab=amb_r;
	}
	void MaterialGetDiffCoeff(double &Pdr,double &Pdg,double &Pdb)
	{
		Pdr=diff_r;
		Pdg=diff_g;
		Pdb=diff_r;
	}
	void MaterialGetSpecCoeff(double &Psr,double &Psg,double &Psb,double &PhCo)
	{
		Psr=spec_r;
		Psg=spec_g;
		Psb=spec_b;
		PhCo=phong_exp;
	}
};

class rtLight
{
	friend class Sphere;
private:
	rtPoint LightPosition;
	color3 ambient_intensity;
	color3 diffuse_intensity;
	color3 specular_intensity;
public:
	void rtLightSet(double x,double y,double z,double Iar,double Iag,double Iab,double Idr,double Idg,double Idb,double Isr,double Isg,double Isb)
	{
		LightPosition.rtPointSet(x,y,z);
		ambient_intensity.set(Iar,Iag,Iab);
		diffuse_intensity.set(Idr,Idg,Idb);
		specular_intensity.set(Isr,Isg,Isb);		
	}
	void rtLightGetPosition(double &x,double &y,double &z)
	{
		LightPosition.rtPointGetPoint(x,y,z);
	}
	void rtLightGetAmbient(double &r,double &g,double &b)
	{
		ambient_intensity.getColor(r,g,b);
	}
	void rtLightGetDiffuse(double &r,double &g,double &b)
	{
		diffuse_intensity.getColor(r,g,b);
	}
	void rtLightGetSpecular(double &r,double &g,double &b)
	{
		specular_intensity.getColor(r,g,b);
	}
};

class Plane
{
private:
	rtPoint point;
	rtVec normal;
	color3 PlaneColor;
	Material surface;
public:
	Plane():point(),normal(),PlaneColor(),surface(){};
	void PlaneSet(double x,double y,double z,double a,double b,double c,double rr,double gg,double bb,
					double Par,double Pag,double Pab,double Pdr,double Pdg,double Pdb,double Psr,double Psg,double Psb,double PhCo)
	{
		normal.rtVecSet(x,y,z);
		normal.normalize();
		point.rtPointSet(a,b,c);
		PlaneColor.set(rr,gg,bb);
		surface.MaterialSet(Par, Pag, Pab, Pdr, Pdg, Pdb, Psr, Psg, Psb, PhCo);
	}
	double Intersect(rtRay R)
	{
		double sx,sy,sz,nx,ny,nz,num,den,dx,dy,dz,D,t;
		rtVec B;
		normal.rtVecGetVec(nx,ny,nz);
		R.rtRayGetSource(sx,sy,sz);
		R.rtRayGetDir(dx,dy,dz);
		
		B.rtVecSet(point);
		D=(-1)*normal.dot(B);
		
		num= nx*sx + ny*sy + nz*sz + D;
		den= nx*dx + ny*dy + nz*dz;
		
		t=(-1)*num/den;
		
		if(t>1)
			return t;
		return 0;
	}
	void ColorRay(rtRay &R,rtLight *L,int num_Lights)
	{
		R.rtDirNormalize();
		color3 lghtCol;
		double Lx,Ly,Lz,Hx,Hy,Hz,r=1,g=1,b=1;
		double Par, Pag, Pab, Pdr, Pdg, Pdb, Psr, Psg, Psb, PhCo;
		// since this is a plane we already have normal, no need to calculate it explicitly as done in case of sphere.
		R.rtRayGetHitPoint(Hx,Hy,Hz);
		for(int i=0;i<num_Lights;i++)
		{
			L[i].rtLightGetPosition(Lx,Ly,Lz);
			rtVec s(Lx-Hx,Ly-Hy,Lz-Hz);
			R.rtRaySetLambert(s,normal,i);
			R.rtRaySetPhong(s,normal,i);
		}

		surface.MaterialGetAmbCoeff(Par, Pag, Pab);
		surface.MaterialGetDiffCoeff(Pdr, Pdg, Pdb);
		surface.MaterialGetSpecCoeff(Psr, Psg, Psb, PhCo);

		for(int i=0;i<num_Lights;i++)
		{
			double Iar, Iag, Iab, Idr, Idg, Idb, Isr, Isg, Isb;
			L[i].rtLightGetAmbient(Iar,Iag,Iab);
			L[i].rtLightGetDiffuse(Idr,Idg,Idb);
			L[i].rtLightGetSpecular(Isr,Isg,Isb);

			r*=Iar*Par + Idr*Pdr*R.rtRayGetLambert(i) + Isr*Psr*pow(R.rtRayGetPhong(i),PhCo);
			g*=Iag*Pag + Idg*Pdg*R.rtRayGetLambert(i) + Isg*Psg*pow(R.rtRayGetPhong(i),PhCo);
			b*=Iab*Pab + Idb*Pdb*R.rtRayGetLambert(i) + Isb*Psb*pow(R.rtRayGetPhong(i),PhCo);
		}

		double j,k,l;
		PlaneColor.getColor(j,k,l);
		lghtCol.set(r*j,g*k,b*l);
		R.rtRaySetRayColor(lghtCol);
	}
};

class Sphere
{
	friend class rtScene;
private:
	rtPoint center;
	double radius;
	color3 color;
	Material surface;
public:
	Sphere():center(),radius(1),color(),surface(){};
	Sphere(double x,double y,double z,double r,double rr,double gg,double bb,double Par,double Pag,double Pab,
			double Pdr,double Pdg,double Pdb,double Psr,double Psg,double Psb,double PhCo)
	{
		center.rtPointSet(x,y,z);
		color.set(rr,gg,bb);
		radius=r;
		surface.MaterialSet(Par, Pag, Pab, Pdr, Pdg, Pdb, Psr, Psg, Psb, PhCo);
	}
	void SphereSet(double x,double y,double z,double r,double rr,double gg,double bb,double Par,double Pag,double Pab,double Pdr,double Pdg,double Pdb,double Psr,double Psg,double Psb,double PhCo)
	{
		//cout<<" "<< y<<" "<< z<<" "<< r<<" "<< rr<<" "<< gg<<" "<< bb<<" "<< Par<<" "<< Pag<<" "<<Pdr<<" "<< Pdg<<" "<< Pdb<<" "<< Psr<<" "<< Psg<<" "<< Psb<<" "<< PhCo;
		center.rtPointSet(x,y,z);
		color.set(rr,gg,bb);
		radius=r;
		surface.MaterialSet(Par, Pag, Pab, Pdr, Pdg, Pdb, Psr, Psg, Psb, PhCo);
	}
	double Intersect(rtRay R)
	{
		double dx,dy,dz,X,Y,Z,A,B,C,t1,t2,D;
		dx=R.dir.x;
		dy=R.dir.y;
		dz=R.dir.z;

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
				if(t1>1 && t2>1)
					return t1>=t2?t2:t1;
				else if(t1>=1)
					return t1;
				else
					return t2;
		}
		return 0;
	}

	// gets the color of the ray for only 1 light source...
	void ColorRay(rtRay &R,rtLight *L,int num_Lights)
	{
		color3 lghtCol;
		double r=1,g=1,b=1;
		R.rtDirNormalize();				// will used in vector v for calculation of phong coefficient....
		rtVec m(R.hitPoint.xx-center.xx,R.hitPoint.yy-center.yy,R.hitPoint.zz-center.zz);	//normal vector at point where ths ray hit the surfae of thi sphere
		m.normalize();
		for(int i=0;i<num_Lights;i++)
		{
			rtVec s(L[i].LightPosition.xx-R.hitPoint.xx,L[i].LightPosition.yy-R.hitPoint.yy,L[i].LightPosition.zz-R.hitPoint.zz);
			R.rtRaySetLambert(s,m,i);
			R.rtRaySetPhong(s,m,i);
		}

		for(int i=0;i<num_Lights;i++)
		{
			r*=(L[i].ambient_intensity.r*surface.amb_r +
				L[i].diffuse_intensity.r*surface.diff_r*R.lambert[i] +
				L[i].specular_intensity.r*surface.spec_r*pow(R.phong[i],surface.phong_exp));

			g*=(L[i].ambient_intensity.g*surface.amb_g +
				L[i].diffuse_intensity.g*surface.diff_g*R.lambert[i] +
				L[i].specular_intensity.g*surface.spec_g*pow(R.phong[i],surface.phong_exp));

			b*=(L[i].ambient_intensity.b*surface.amb_b +
				L[i].diffuse_intensity.b*surface.diff_b*R.lambert[i] +
				L[i].specular_intensity.b*surface.spec_b*pow(R.phong[i],surface.phong_exp));
		}
		
		lghtCol.set(r*color.r,g*color.g,b*color.b);
		R.rtRaySetRayColor(lghtCol);
	}
};

class rtCamera
{
	friend class rtScene;
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
		near=n;				
		far=f;				
	}
	void rtCameraSetRays()
	{
		width=right-left;
		height=top-bottom;
		rays=new rtRay*[(int)height];
		rtPoint p;
		for(int i=0;i<height;i++)
			rays[i]=new rtRay[(int)width];
		for(int i=0,b=bottom;i<height;i++,b++)
			for(int j=0,l=left;j<width;j++,l++)
			{
				rays[i][j].rtRaySetSource(eye);
				p.rtPointSet(eye.xx -near*n.x + u.x*b +v.x*l, 
							 eye.yy -near*n.y + u.y*b +v.y*l, 
							 eye.zz -near*n.z + u.z*b +v.z*l);						//z axis point of the ray...
				rays[i][j].rtRaySetDestination(p);
				rays[i][j].rtRaySetDir();										// to set the ray direction..
			}
	}			
	void rtCameraCastRay(Sphere *S,int num_Sphere,rtLight *L,int num_Lights,Plane *P,int num_Planes)
	{
		for(int i=0;i<height;i++)
		{
			for(int j=0;j<width;j++)
			{
				double min=far,tmp;
				int object=-1,pl=-1;
				for(int k=0;k<num_Sphere;k++)
				{
					tmp=S[k].Intersect(rays[i][j]);
					if(tmp>=1 && tmp<min)
					{
						min=tmp;
						object=k;
					}
				}
				for(int k=0;k<num_Planes;k++)
				{
					tmp=P[k].Intersect(rays[i][j]);
					if(tmp>1 && tmp<min)
					{
						min=tmp;
						pl=k;
					}
				}
				if(pl!=-1)
				{
					rays[i][j].rtRaySetHitPoint(min);
					rays[i][j].rtRaySetCoeffSize(num_Lights);
					P[pl].ColorRay(rays[i][j],L,num_Lights);
				}
				else if(object!=-1)
				{
					rays[i][j].rtRaySetHitPoint(min);
					rays[i][j].rtRaySetCoeffSize(num_Lights);
					S[object].ColorRay(rays[i][j],L,num_Lights);
				}				
			}
		}
	}
	void rtvalue(int i,int j)
	{
		glColor3f(rays[i][j].rayColor.r,rays[i][j].rayColor.g,rays[i][j].rayColor.b);
	}
};

class rtScene
{
private:
	Sphere *Sph_obj;
	Plane *Planes;
	//Sphere Sph_obj[3];
	int nS,nL,nP;
	rtCamera CAM;
	rtLight *Lights;
public:
	rtScene(char *fname)
	{
		FILE *fp;
		fp=fopen(fname,"r");
		double a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s;
		char ch;


		//.........read camera.........//
		fscanf(fp,"%lf,%lf,%lf,%lf,%lf,%lf\n",&a,&b,&c,&d,&e,&f);					// 6 values
		this->CAM.rtCameraWindow(a,b,c,d,e,f);
		fscanf(fp,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",&a,&b,&c,&d,&e,&f,&g,&h,&i);	// 9 values
		this->CAM.rtCameraSet(a,b,c,d,e,f,g,h,i);
		//.............................//

		//.......read Lights.............//
		fscanf(fp,"%c %d\n",&ch,&nL);
		cout<<nL;
		if(ch=='l')
		{
			this->Lights=new rtLight[nL];
			for(int z=0;z<nL;z++)
			{
				fscanf(fp,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",&a,&b,&c,&d,&e,&f,&g,&h,&i,&j,&k,&l);	//12 values
				cout<<"\n"<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<f<<" "<<g<<" "<<h<<" "<<i<<" "<<j<<" "<<k<<" "<<l;
				this->Lights[z].rtLightSet(a,b,c,d,e,f,g,h,i,j,k,l);
			}
		}
		//..............................//


		//.....read spheres.............//
		fscanf(fp,"%c %d\n",&ch,&nS);
		cout<<nS;
		if(ch=='s')
		{
			Sph_obj=new Sphere[nS];
			//Sph_obj=(class Sphere *) malloc(sizeof(Sphere)*nS);
			for(int z=0;z<nS;z++)
			{
				//cout<<Sph_obj[i].center.xx;
				fscanf(fp,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",&a,&b,&c,&d,&e,&f,&g,&h,&i,&j,&k,&l,&m,&n,&o,&p,&q);
				cout<<"\n"<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<f<<" "<<g<<" "<<h<<" "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<" "<<n<<" "<<" "<<o<<" "<<" "<<p<<" "<<" "<<q<<" ";
				Sph_obj[z].SphereSet(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q);
			}
		}
		//.............................//

		//.........read Planes.........//
		fscanf(fp,"%c %d\n",&ch,&nP);
		if(ch=='p')
		{
			Planes=new Plane[nP];
			for(int z=0;z<nP;z++)
			{
				fscanf(fp,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",&a,&b,&c,&d,&e,&f,&g,&h,&i,&j,&k,&l,&m,&n,&o,&p,&q,&r,&s);
				cout<<"\n"<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<f<<" "<<g<<" "<<h<<" "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<" "<<n<<" "<<" "<<o<<" "<<" "<<p<<" "<<" "<<q<<" "<<r<<" "<<s<<" ";
				Planes[z].PlaneSet(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s);
			}
		}
		//..................................//
	}
	void SceneRender()
	{
		CAM.rtCameraSetRays();
		CAM.rtCameraCastRay(Sph_obj,nS,Lights,nL,Planes,nP);
	}
	void GetSceneDimentions(double &l,double &r,double &b,double &t)
	{
		l=CAM.left;
		r=CAM.right;
		b=CAM.bottom;
		t=CAM.top;
	}		
	void SceneRayValue(int i,int j)
	{
		CAM.rtvalue(i,j);
	}
};
#endif