#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=2000;
const int Ly=1;

const int Q=5;
const double W0=1.0/3;

const double C=0.5; // C<0.707 celdas/click
const double TresC2=3*C*C;
const double AUX0=1-TresC2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

const double D=0.6;

class LatticeBoltzmann{
private:
  double w[Q];
  int V[2][Q];  //V[alpha][i]  alpha=0 es x, alpha=1 es y
  double f[Lx][Ly][Q],fnew[Lx][Ly][Q]; // f[ix][iy][i]
  double densidad[Lx];
public:
  LatticeBoltzmann(void);
  double rho(int ix,int iy,bool UseNew);
  double Jx(int ix,int iy);
  double Jy(int ix,int iy);
  double feq(int i,double rho0,double Jx0,double Jy0);
  void Inicie(double rho0,double Jx0,double Jy0);
  void ImponerCampos(int ix,int iy,double & rho0,double & Jx0,double & Jy0,int t);
  void Colisione(int t);
  void Adveccione(void);
  void Imprimase(char const * NombreArchivo, int t);
  void ImprimaUnaLinea(char const * NombreArchivo,int t);
  double getdensidad(int ix){return densidad[ix];};
};
LatticeBoltzmann::LatticeBoltzmann(void){
  w[0]=W0;
  w[1]=w[2]=w[3]=w[4]=1.0/6;

  V[0][0]=0;  
  V[1][0]=0;

  V[0][1]=1;    V[0][2]=0;    V[0][3]=-1;   V[0][4]=0;  
  V[1][1]=0;    V[1][2]=1;    V[1][3]=0;    V[1][4]=-1;  
}
double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew)
      suma+=fnew[ix][iy][i];
    else
      suma+=f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::Jx(int ix,int iy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    suma+=V[0][i]*f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::Jy(int ix,int iy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    suma+=V[1][i]*f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::feq(int i,double rho0,double Jx0,double Jy0){
  if(i==0)
    return AUX0*rho0;
  else
    return w[i]*(TresC2*rho0+3*(V[0][i]*Jx0+V[1][i]*Jy0));
}
void LatticeBoltzmann::ImponerCampos(int ix,int iy,double & rho0,double & Jx0,double & Jy0,int t){
  double A=1,lambda=1000,omega=2*M_PI*C/lambda;
  if(ix==0 && iy==0)
    rho0=A*sin(omega*t);
}

void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0){
	int ix,iy,i;
	for(iy=0;iy<Ly;iy++){
		for(ix=0;ix<Lx;ix++){
			for(i=0;i<Q;i++){
				f[ix][iy][i]=feq(i,rho0,Jx0,Jy0);
				}
			}
	  }
}
void LatticeBoltzmann::Colisione(int t){
	int ix,iy,i; double rhonuevo, rho0, Jx0, Jy0;
	for(iy=0;iy<Ly;iy++){
		for(ix=0;ix<Lx;ix++){
			rho0=rho(ix,iy,false);  Jx0=Jx(ix,iy);  Jy0=Jy(ix,iy);
			ImponerCampos(ix,iy,rho0,Jx0,Jy0,t);
			if(ix==Lx-1){
				fnew[ix][iy][1]=D*f[ix][iy][3];
				fnew[ix][iy][3]=D*f[ix][iy][1];
				fnew[ix][iy][2]=D*f[ix][iy][4];
				fnew[ix][iy][4]=D*f[ix][iy][2];
			}else{
				for(i=0;i<Q;i++){
					fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(i,rho0,Jx0,Jy0);
				}
			}
			rhonuevo=rho(ix, iy, true);  densidad[ix]=rhonuevo;
		}
	}
}
void LatticeBoltzmann::Adveccione(void){
	int ix,iy,i;
	for(iy=0;iy<Ly;iy++){
		for(ix=0;ix<Lx;ix++){
			for(i=0;i<Q;i++){
				f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
			}
		}
	}
}
void LatticeBoltzmann::Imprimase(char const * NombreArchivo, int t){
  ofstream MiArchivo(NombreArchivo); double rho0,Jx0,Jy0;
  for(int iy=0;iy<Ly;iy++){
    for(int ix=0;ix<Lx-1;ix++){
      rho0=rho(ix,iy,true);   Jx0=Jx(ix,iy);  Jy0=Jy(ix,iy); //densidad[ix][t]=rho0;
      ImponerCampos(ix,iy,rho0,Jx0,Jy0,t);
      //MiArchivo<<ix<<" "<<iy<<" "<<rho0<<endl;
      MiArchivo<<ix<<" "<<rho0<<endl;
    }
    MiArchivo<<endl;
  }
  MiArchivo.close();
}
void LatticeBoltzmann::ImprimaUnaLinea(char const* NombreArchivo,int t){
  ofstream MiArchivo(NombreArchivo); double rho0,Jx0,Jy0;
  int ix=Lx/2;
  for(int iy=0;iy<Ly;iy++){
    rho0=rho(ix,iy,true);   Jx0=Jx(ix,iy);  Jy0=Jy(ix,iy);
    ImponerCampos(ix,iy,rho0,Jx0,Jy0,t);
    MiArchivo<<iy<<" "<<rho0<<endl;
  }
  MiArchivo.close();
}
//---------------- Funciones Globales --------

int main(void){
  LatticeBoltzmann Ondas;
  int i, ix, t,tmax;
  
  double densidad_max[Lx];
  double densidad_min[Lx];
  
  double Amax_pmax=0.5;
  double Amin_pmax=0.5;
  
  double SWR;
  double Cabsorcion;
  
  for(ix=0;ix<=Lx;ix++){
	  densidad_max[ix]=0.0;
	  densidad_min[ix]=0.0;
  }
  
  double rho0=0,Jx0=0,Jy0=0;

  //Inicie
  Ondas.Inicie(rho0,Jx0,Jy0);
  tmax=80000;
  for(t=0;t<tmax;t++){
    Ondas.Colisione(t);
    Ondas.Adveccione();
  }
    
  tmax=2200;
  for(t=0;t<=tmax;t++){
    Ondas.Colisione(t);
    Ondas.Adveccione();
    for(ix=0;ix<=Lx-2;ix++){
		if(Ondas.getdensidad(ix) > densidad_max[ix]) densidad_max[ix]=Ondas.getdensidad(ix);
		if(Ondas.getdensidad(ix) < densidad_min[ix]) densidad_min[ix]=Ondas.getdensidad(ix);
	}
  }
  
  for(ix=0;ix<=Lx-2;ix++){
	  if(densidad_max[ix] > Amax_pmax)  Amax_pmax=densidad_max[ix];
	  if(densidad_max[ix] < Amin_pmax)  Amin_pmax=densidad_max[ix];
  }
  
  //cout << Amin_pmax << " " << Amax_pmax << endl;
  
  SWR = Amax_pmax/Amin_pmax;
  Cabsorcion=(4.0*SWR)/((SWR+1)*(SWR+1));
  
  cout << SWR << " " << Cabsorcion << endl;
  
  //Ondas.Imprimase("Ondas.dat",t);
  //for(ix=0;ix<=Lx-2;ix++) densidad_80000_mas_2200[ix]=Ondas.getdensidad(ix);
  //for(ix=0;ix<=Lx-2;ix++) cout << ix << " " << densidad_80000[ix] << " " << densidad_80000_mas_2200[ix] << endl;
  
  return 0;
}

//Valores de:
//D=0.6
//SRW= 9.39545 
//Cabsorcion= 0.347769
