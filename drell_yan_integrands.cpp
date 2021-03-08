#include <iostream>
#include <tgmath.h>


using namespace std;




class flows {

    private:
        double p1[4];
        double p2[4];

    public:
        flows(){};
        flows(double* np1, double* np2){for(int i=0; i<4; i++){p1[i]=np1[i]; p2[i]=np2[i];}};

        double norm(double vx, double vy, double vz){return pow(pow(vx,2.)+pow(vy,2.)+pow(vz,2.),0.5);};

        double sp(double vx, double vy, double vz, double wx, double wy, double wz){return vx*wx+vy*wy+vz*wz;}

        double flow1(double kx, double ky, double kz, double lx, double ly, double lz, bool var){
            if(var==0){
                return (p1[0]+p2[0])/(norm(kx,ky,kz)+norm(kx-p1[1]-p2[1],ky-p1[2]-p2[2],kz-p1[3]-p2[3]));
                };
            if(var==1){
                return (p1[0]+p2[0])/(norm(lx,ly,lz)+norm(lx-p1[1]-p2[1],ly-p1[2]-p2[2],lz-p1[3]-p2[3]));
                };
            return 0;
        }
        double flow2(double kx, double ky, double kz, double lx, double ly, double lz, bool var){
            if(var==0){
            return (p1[0]+p2[0])/(norm(kx,ky,kz)+norm(kx-lx,ky-ly,kz-lz)+norm(lx+p1[1]+p2[1],ly+p1[2]+p2[2],lz+p1[3]+p2[3]));
                }
            if(var==1){
            return (p1[0]+p2[0])/(norm(lx,ly,lz)+norm(kx-lx,ky-ly,kz-lz)+norm(kx+p1[1]+p2[1],ky+p1[2]+p2[2],kz+p1[3]+p2[3]));
                }
            return 0;
        }
        double flow3(double kx, double ky, double kz, double lx, double ly, double lz){
            return -(p1[0]+p2[0])/(-norm(lx,ly,lz)+norm(kx+lx,ky+ly,kz+lz)+norm(kx+p1[1]+p2[1],ky+p1[2]+p2[2],kz+p1[3]+p2[3]));
        }
        double flow4(double kx, double ky, double kz, double lx, double ly, double lz){
            return -(p1[0]+p2[0])/(norm(kx+lx,ky+ly,kz+lz)+norm(kx+lx+p1[1]+p2[1],ky+ly+p1[2]+p2[2],kz+lz+p1[3]+p2[3]));
        }
        double flow5(double kx, double ky, double kz, double lx, double ly, double lz){
            return -(p1[0]+p2[0])/(norm(kx,ky,kz)+norm(lx,ly,lz)+norm(kx+lx+p1[1]+p2[1],ky+ly+p1[2]+p2[2],kz+lz+p1[3]+p2[3]));
        }
        double flow6(double kx, double ky, double kz, double lx, double ly, double lz){
            return -(p1[0]+p2[0])/(norm(lx,ly,lz)+norm(kx-lx+p1[1]+p2[1],ky-ly+p1[2]+p2[2],kz-lz+p1[3]+p2[3]));
        }

        double flow7(double kx, double ky, double kz, double lx, double ly, double lz){
            return (p1[0]+p2[0])/(norm(lx,ly,lz)+norm(kx-lx+p1[1]+p2[1],ky-ly+p1[2]+p2[2],kz-lz+p1[3]+p2[3])-norm(kx,ky,kz));
        }

        double jacques1(double kx, double ky, double kz, double lx, double ly, double lz, bool var){
            if(var==0){
                return 2*norm(kx,ky,kz);
            }
            if(var==1){
                return 2*norm(lx,ly,lz);
            }
            return 0;
        }

        double jacques2(double kx, double ky, double kz, double lx, double ly, double lz, bool var){
            if(var==0){
                return norm(kx,ky,kz)+norm(kx-lx,ky-ly,kz-lz)+norm(lx,ly,lz);
            }
            if(var==1){
                return norm(lx,ly,lz)+norm(kx-lx,ky-ly,kz-lz)+norm(kx,ky,kz);
            }
        return 0;
        }

        double jacques3(double kx, double ky, double kz, double lx, double ly, double lz){
            return abs(-norm(lx,ly,lz)+norm(kx+lx,ky+ly,kz+lz)+norm(kx,ky,kz));
        }

        double jacques4(double kx, double ky, double kz, double lx, double ly, double lz){
            return norm(kx+lx,ky+ly,kz+lz)+norm(kx+lx,ky+ly,kz+lz);
        }       

        double jacques5(double kx, double ky, double kz, double lx, double ly, double lz){
            return norm(kx,ky,kz)+norm(lx,ly,lz)+norm(kx+lx,ky+ly,kz+lz);
        }

        double jacques6(double kx, double ky, double kz, double lx, double ly, double lz){
            return norm(lx,ly,lz)+norm(kx-lx,ky-ly,kz-lz);
        }

        double jacques7(double kx, double ky, double kz, double lx, double ly, double lz){
            return abs(norm(lx,ly,lz)+norm(kx-lx,ky-ly,kz-lz)-norm(kx,ky,kz));
        }

};



class integrands{

    private:
        double p1[4];
        double p2[4];
        double p10;
        double p1x;
        double p1y;
        double p1z;
        double p20;
        double p2x;
        double p2y;
        double p2z;
        double MUV;
        flows f;
        double sigma;
        double constant;

    public:
        integrands(double* np1, double* np2, double nMUV, flows nf, double nsigma){for(int i=0; i<4; i++){p1[i]=np1[i]; p2[i]=np2[i];}; f=nf;
        p10=p1[0];
        p1x=p1[1];
        p1y=p1[2];
        p1z=p1[3];
        p20=p2[0];
        p2x=p2[1];
        p2y=p2[2];
        p2z=p2[3];
        MUV=nMUV;
        sigma=nsigma;
        constant=0.389379*pow(10,9.)*pow(4*2*acos(0.0)/132.507,2.)*(4*2*acos(0.0)*0.118/9)/(pow(2*2*acos(0.0),5.)*4*8*pow(p10*p20-p1x*p2x-p1y*p2y-p1z*p2z,3.));
        };

        double Power(double x, double c){return pow(x,c);}
        double Sqrt(double x){return pow(x,0.5);}
        double Abs(double x){return abs(x);}

        double dt_num (double k0, double kx, double ky, double kz, double l0, double lx, double ly, double lz){
            return 512*(Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 256*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 256*Power(l0*p10 - lx*p1x - ly*p1y - lz*p1z,2)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 256*(k0*l0 - kx*lx - ky*ly - kz*lz)*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 
                256*(Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 256*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                Power(k0*p20 - kx*p2x - ky*p2y - kz*p2z,2) - 512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 256*Power(k0*p10 - kx*p1x - ky*p1y - kz*p1z,2)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 512*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 
                256*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 256*(k0*l0 - kx*lx - ky*ly - kz*lz)*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 
                256*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 
                256*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                512*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 256*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                Power(l0*p20 - lx*p2x - ly*p2y - lz*p2z,2) + 512*Power(k0*l0 - kx*lx - ky*ly - kz*lz,2)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                512*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 512*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                256*(Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                512*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2));
        }

        double dt_num_UV_left(double k0, double kx, double ky, double kz, double l0, double lx, double ly, double lz){
            return (32*(Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*p10*p20)/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) - 
                (32*l0*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*p20)/Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) - 
                (16*p10*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*p20)/Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) + 
                (16*l0*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*p20)/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) + 
                (16*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*Power(p20,2))/Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) - 
                (96*(Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*(-(kx*p1x) - ky*p1y - kz*p1z)*(-(kx*p2x) - ky*p2y - kz*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) + 
                (96*(-(kx*lx) - ky*ly - kz*lz)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(-(kx*p2x) - ky*p2y - kz*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) + 
                (48*(-(kx*p1x) - ky*p1y - kz*p1z)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(-(kx*p2x) - ky*p2y - kz*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) - 
                (48*(-(kx*lx) - ky*ly - kz*lz)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(-(kx*p2x) - ky*p2y - kz*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) - 
                (48*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*Power(-(kx*p2x) - ky*p2y - kz*p2z,2))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) - 
                (32*l0*p10*(l0*p20 - lx*p2x - ly*p2y - lz*p2z))/Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) + 
                (16*Power(p10,2)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z))/Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) + 
                (96*(-(kx*lx) - ky*ly - kz*lz)*(-(kx*p1x) - ky*p1y - kz*p1z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) - 
                (48*Power(-(kx*p1x) - ky*p1y - kz*p1z,2)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) + 
                (96*Power(MUV,2)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) + 
                (128*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) - 
                (48*Power(MUV,2)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) - 
                (64*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) - 
                (16*p10*p20*(l0*p20 - lx*p2x - ly*p2y - lz*p2z))/Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) + 
                (48*(-(kx*p1x) - ky*p1y - kz*p1z)*(-(kx*p2x) - ky*p2y - kz*p2z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) - 
                (96*Power(-(kx*lx) - ky*ly - kz*lz,2)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) + 
                (48*(-Power(kx,2) - Power(ky,2) - Power(kz,2))*(Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) + 
                (32*Power(l0,2)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) - 
                (16*(Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) - 
                (16*l0*p10*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) + 
                (48*(-(kx*lx) - ky*ly - kz*lz)*(-(kx*p1x) - ky*p1y - kz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) - 
                (48*(-Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) - 
                (48*Power(MUV,2)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) - 
                (48*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) - 
                (16*l0*p20*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) + 
                (48*(-(kx*lx) - ky*ly - kz*lz)*(-(kx*p2x) - ky*p2y - kz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) - 
                (48*(-Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) - 
                (48*Power(MUV,2)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) - 
                (48*(l0*p20 - lx*p2x - ly*p2y - lz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) + 
                (16*l0*p10*(Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5) - 
                (48*(-(kx*lx) - ky*ly - kz*lz)*(-(kx*p1x) - ky*p1y - kz*p1z)*(Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) - 
                (48*Power(MUV,2)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),2.5) - 
                (64*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)))/
                Power(Power(kx,2) + Power(ky,2) + Power(kz,2) + Power(MUV,2),1.5);
        }

        double dt_num_UV_right(double k0, double kx, double ky, double kz, double l0, double lx, double ly, double lz){
            return (32*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*p10*p20)/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) - 
                (32*k0*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*p20)/Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) - 
                (16*p10*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*p20)/Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) + 
                (16*k0*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*p20)/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) + 
                (16*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*Power(p20,2))/Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) - 
                (32*k0*p10*(k0*p20 - kx*p2x - ky*p2y - kz*p2z))/Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) + 
                (16*Power(p10,2)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z))/Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) + 
                (96*Power(MUV,2)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) + 
                (128*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) + 
                (96*(-(kx*lx) - ky*ly - kz*lz)*(-(lx*p1x) - ly*p1y - lz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) - 
                (48*Power(-(lx*p1x) - ly*p1y - lz*p1z,2)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) - 
                (48*Power(MUV,2)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) - 
                (64*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) - 
                (16*p10*p20*(k0*p20 - kx*p2x - ky*p2y - kz*p2z))/Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) + 
                (96*(-(kx*lx) - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(-(lx*p2x) - ly*p2y - lz*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) - 
                (96*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(-(lx*p1x) - ly*p1y - lz*p1z)*(-(lx*p2x) - ly*p2y - lz*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) + 
                (48*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(-(lx*p1x) - ly*p1y - lz*p1z)*(-(lx*p2x) - ly*p2y - lz*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) - 
                (48*(-(kx*lx) - ky*ly - kz*lz)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(-(lx*p2x) - ly*p2y - lz*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) + 
                (48*(-(lx*p1x) - ly*p1y - lz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*(-(lx*p2x) - ly*p2y - lz*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) - 
                (48*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*Power(-(lx*p2x) - ly*p2y - lz*p2z,2))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) - 
                (96*Power(-(kx*lx) - ky*ly - kz*lz,2)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) + 
                (48*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(-Power(lx,2) - Power(ly,2) - Power(lz,2))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) + 
                (32*Power(k0,2)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) - 
                (16*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) - 
                (16*k0*p10*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) - 
                (48*(-Power(lx,2) - Power(ly,2) - Power(lz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) - 
                (48*Power(MUV,2)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) - 
                (48*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) + 
                (48*(-(kx*lx) - ky*ly - kz*lz)*(-(lx*p1x) - ly*p1y - lz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) - 
                (16*k0*p20*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) - 
                (48*(-Power(lx,2) - Power(ly,2) - Power(lz,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) - 
                (48*Power(MUV,2)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) - 
                (48*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) + 
                (48*(-(kx*lx) - ky*ly - kz*lz)*(-(lx*p2x) - ly*p2y - lz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) + 
                (16*k0*p10*(Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) - 
                (48*Power(MUV,2)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5) - 
                (64*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),1.5) - 
                (48*(-(kx*lx) - ky*ly - kz*lz)*(-(lx*p1x) - ly*p1y - lz*p1z)*(Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)))/
                Power(Power(lx,2) + Power(ly,2) + Power(lz,2) + Power(MUV,2),2.5);
        }

        double t_channel_num(double k0, double kx, double ky, double kz, double l0, double lx, double ly, double lz){
            return -512*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 1024*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2));
        }

        double t_channel_num_derivative(double k0, double kx, double ky, double kz, double l0, double lx, double ly, double lz){
            return (-256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*p20 - 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*p20 - 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*p10*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 
                1024*k0*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 
                512*l0*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 
                256*k0*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 
                512*k0*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 
                256*l0*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 
                256*k0*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                256*k0*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*p10*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                512*k0*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*l0*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                256*k0*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*p20*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                512*k0*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*l0*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                256*k0*(l0*p20 - lx*p2x - ly*p2y - lz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*p10*(Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                512*k0*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                256*l0*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                256*k0*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)))/
                (Power(k0 + l0,2) - Power(kx + lx,2) - Power(ky + ly,2) - Power(kz + lz,2)) - 
                (2*(k0 + l0)*(-512*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 1024*(k0*l0 - kx*lx - ky*ly - kz*lz)*
                (k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2))))/
                Power(Power(k0 + l0,2) - Power(kx + lx,2) - Power(ky + ly,2) - Power(kz + lz,2),2);
        }

        double u_channel_num(double k0, double kx, double ky, double kz, double l0, double lx, double ly, double lz){
            return 128*(Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 128*(k0*l0 - kx*lx - ky*ly - kz*lz)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 64*Power(l0*p10 - lx*p1x - ly*p1y - lz*p1z,2)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 64*(Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 
                128*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*Power(k0*p20 - kx*p2x - ky*p2y - kz*p2z,2) - 
                128*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                128*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 64*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 64*(k0*l0 - kx*lx - ky*ly - kz*lz)*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                64*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                128*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 64*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                Power(l0*p20 - lx*p2x - ly*p2y - lz*p2z,2) + 
                64*Power(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2),2)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                128*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*l0 - kx*lx - ky*ly - kz*lz)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 128*Power(k0*l0 - kx*lx - ky*ly - kz*lz,2)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 64*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                128*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 128*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 128*Power(k0*p10 - kx*p1x - ky*p1y - kz*p1z,2)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 64*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                64*(k0*l0 - kx*lx - ky*ly - kz*lz)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                128*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 64*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                64*(k0*l0 - kx*lx - ky*ly - kz*lz)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 128*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                128*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                128*Power(k0*p20 - kx*p2x - ky*p2y - kz*p2z,2)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                64*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 64*(k0*l0 - kx*lx - ky*ly - kz*lz)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 128*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                64*(Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                128*Power(k0*p10 - kx*p1x - ky*p1y - kz*p1z,2)*(Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                64*(k0*l0 - kx*lx - ky*ly - kz*lz)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                128*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                128*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                128*(k0*l0 - kx*lx - ky*ly - kz*lz)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                64*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                64*(k0*l0 - kx*lx - ky*ly - kz*lz)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2));
        }

        double s_channel_num(double k0, double kx, double ky, double kz, double l0, double lx, double ly, double lz){
            return -512*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 1024*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 1024*Power(k0*p10 - kx*p1x - ky*p1y - kz*p1z,2)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 
                512*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 1024*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 
                256*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 256*Power(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2),2)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 1024*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                Power(k0*p20 - kx*p2x - ky*p2y - kz*p2z,2) - 
                512*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*Power(k0*p20 - kx*p2x - ky*p2y - kz*p2z,2) - 
                512*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*Power(k0*p20 - kx*p2x - ky*p2y - kz*p2z,2) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 512*Power(k0*p10 - kx*p1x - ky*p1y - kz*p1z,2)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 256*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                512*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                512*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                512*Power(k0*p10 - kx*p1x - ky*p1y - kz*p1z,2)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                512*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 256*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 2048*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                768*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 512*Power(k0*p20 - kx*p2x - ky*p2y - kz*p2z,2)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 512*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                512*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*Power(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z,2) - 
                512*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*Power(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z,2) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                512*Power(k0*p10 - kx*p1x - ky*p1y - kz*p1z,2)*(Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                512*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                256*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                1024*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                256*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                256*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                256*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                768*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                256*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                256*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*Power(Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2),2);
        }

        double st_channel_num(double k0, double kx, double ky, double kz, double l0, double lx, double ly, double lz){
            return 512*(Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 256*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 256*Power(l0*p10 - lx*p1x - ly*p1y - lz*p1z,2)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 256*(k0*l0 - kx*lx - ky*ly - kz*lz)*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 
                256*(Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 
                256*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*Power(k0*p20 - kx*p2x - ky*p2y - kz*p2z,2) - 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 
                256*Power(k0*p10 - kx*p1x - ky*p1y - kz*p1z,2)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                512*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 256*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 256*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 256*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 256*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                Power(l0*p20 - lx*p2x - ly*p2y - lz*p2z,2) + 
                512*Power(k0*l0 - kx*lx - ky*ly - kz*lz,2)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 512*(k0*l0 - kx*lx - ky*ly - kz*lz)*
                Power(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z,2) - 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                256*(Power(l0,2) - Power(lx,2) - Power(ly,2) - Power(lz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                256*(k0*l0 - kx*lx - ky*ly - kz*lz)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2));
        }

        double bub_num(double k0, double kx, double ky, double kz, double l0, double lx, double ly, double lz){
            return -1024*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 512*(k0*l0 - kx*lx - ky*ly - kz*lz)*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2));
        }

        double bub_num_derivative(double k0, double kx, double ky, double kz, double l0, double lx, double ly, double lz){
            return -((-Sqrt(Power(lx,2) + Power(ly,2) + Power(lz,2)) + p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (512*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (-(kx*lx) - ky*ly - kz*lz + Sqrt(Power(lx,2) + Power(ly,2) + Power(lz,2))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z) - 
                1024*(-(kx*lx) - ky*ly - kz*lz + Sqrt(Power(lx,2) + Power(ly,2) + Power(lz,2))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p1x) - ky*p1y - kz*p1z + p10*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z) + 
                256*(Sqrt(Power(lx,2) + Power(ly,2) + Power(lz,2))*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z) - 
                256*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (Sqrt(Power(lx,2) + Power(ly,2) + Power(lz,2))*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                256*(-(kx*p1x) - ky*p1y - kz*p1z + p10*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (Sqrt(Power(lx,2) + Power(ly,2) + Power(lz,2))*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                512*(-(kx*lx) - ky*ly - kz*lz + Sqrt(Power(lx,2) + Power(ly,2) + Power(lz,2))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p1x) - ky*p1y - kz*p1z + p10*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(Sqrt(Power(lx,2) + Power(ly,2) + Power(lz,2))*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                512*(-(kx*lx) - ky*ly - kz*lz + Sqrt(Power(lx,2) + Power(ly,2) + Power(lz,2))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (Sqrt(Power(lx,2) + Power(ly,2) + Power(lz,2))*p20 - lx*p2x - ly*p2y - lz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                512*(-(kx*lx) - ky*ly - kz*lz + Sqrt(Power(lx,2) + Power(ly,2) + Power(lz,2))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p1x) - ky*p1y - kz*p1z + p10*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                256*(Sqrt(Power(lx,2) + Power(ly,2) + Power(lz,2))*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2))))/
                (8.*Power(-Power(kx - lx,2) - Power(ky - ly,2) - Power(kz - lz,2) + 
                Power(-Sqrt(Power(lx,2) + Power(ly,2) + Power(lz,2)) + p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2),2)*
                Abs(Sqrt(Power(lx,2) + Power(ly,2) + Power(lz,2))*
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2)*
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))) - 
                ((Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2)) + p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (256*(-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (-(lx*p1x) - ly*p1y - lz*p1z + p10*(Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2)) + p10 + 
                p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z) + 
                512*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (-(kx*lx) - ky*ly - kz*lz + (p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2)) + p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z) - 
                1024*(-(kx*p1x) - ky*p1y - kz*p1z + p10*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*lx) - ky*ly - kz*lz + (p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2)) + p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z) - 
                256*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (-(lx*p2x) - ly*p2y + p20*(Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2)) + p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - lz*p2z) + 
                256*(-(kx*p1x) - ky*p1y - kz*p1z + p10*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (-(lx*p2x) - ly*p2y + p20*(Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2)) + p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - lz*p2z) - 
                256*(-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (-(lx*p1x) - ly*p1y - lz*p1z + p10*(Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2)) + p10 + 
                p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                512*(-(kx*p1x) - ky*p1y - kz*p1z + p10*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*lx) - ky*ly - kz*lz + (p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2)) + p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                512*(-(kx*lx) - ky*ly - kz*lz + (p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2)) + p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (-(lx*p2x) - ly*p2y + p20*(Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2)) + p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - lz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (-(lx*p1x) - ly*p1y - lz*p1z + p10*(Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2)) + p10 + 
                p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                512*(-(kx*p1x) - ky*p1y - kz*p1z + p10*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*lx) - ky*ly - kz*lz + (p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2)) + p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2))))/
                (8.*Power(-Power(lx,2) - Power(ly,2) - Power(lz,2) + 
                Power(Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2)) + p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2),2)*
                Abs(Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2))*
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2)*
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))) + 
                (-256*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*p20*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2)) + 
                256*p20*(-(kx*p1x) - ky*p1y - kz*p1z + p10*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2)) + 
                256*p10*(-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + 
                Power(kz - p1z - p2z,2))) - kz*p2z) + 
                512*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + 
                Power(kz - p1z - p2z,2))) - kz*p2z) - 
                1024*(-(kx*p1x) - ky*p1y - kz*p1z + p10*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + 
                Power(kz - p1z - p2z,2))) - kz*p2z) - 
                256*p10*(-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*p20*(-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                512*(-(kx*p1x) - ky*p1y - kz*p1z + p10*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                512*(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + 
                Power(kz - p1z - p2z,2))) - kz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*p10*(-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                512*(-(kx*p1x) - ky*p1y - kz*p1z + p10*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)))/
                (16.*(-Power(lx,2) - Power(ly,2) - Power(lz,2) + 
                Power(Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2)) + p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                Abs(Sqrt(Power(-kx + lx,2) + Power(-ky + ly,2) + Power(-kz + lz,2))*
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2)*
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))));

        }

        double t_channel_num_2(double k0, double kx, double ky, double kz, double l0, double lx, double ly, double lz){
            return -1024*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 512*(k0*l0 - kx*lx - ky*ly - kz*lz)*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (l0*p20 - lx*p2x - ly*p2y - lz*p2z)*(p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2));
        }

        double t_channel_num_derivative_2(double k0, double kx, double ky, double kz, double l0, double lx, double ly, double lz){
            return (-2*(k0 + l0)*(-1024*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (l0*p10 - lx*p1x - ly*p1y - lz*p1z)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z) - 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (k0*p20 - kx*p2x - ky*p2y - kz*p2z) + 256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (k0*p10 - kx*p1x - ky*p1y - kz*p1z)*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*
                (Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z) - 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p20 - kx*p2x - ky*p2y - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p20 - lx*p2x - ly*p2y - lz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                512*(k0*l0 - kx*lx - ky*ly - kz*lz)*(k0*p10 - kx*p1x - ky*p1y - kz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                256*(Power(k0,2) - Power(kx,2) - Power(ky,2) - Power(kz,2))*(l0*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2))))/
                Power(Power(k0 + l0,2) - Power(kx + lx,2) - Power(ky + ly,2) - Power(kz + lz,2),2);
        }

        double bub_UV_num(double k0, double kx, double ky, double kz, double l0, double lx, double ly, double lz){
            return (((3/(32.*Power(1 + Power(lx,2) + Power(ly,2) + Power(lz,2),2)) - 
                (3*(Sqrt(1 + Power(lx,2) + Power(ly,2) + Power(lz,2)) - p10 - p20 + 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))/
                (32.*Power(1 + Power(lx,2) + Power(ly,2) + Power(lz,2),2.5)))*
                (512*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (-(kx*lx) - ky*ly - kz*lz + Sqrt(1 + Power(lx,2) + Power(ly,2) + Power(lz,2))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z) - 
                1024*(-(kx*lx) - ky*ly - kz*lz + Sqrt(1 + Power(lx,2) + Power(ly,2) + Power(lz,2))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p1x) - ky*p1y - kz*p1z + p10*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z) + 
                256*(Sqrt(1 + Power(lx,2) + Power(ly,2) + Power(lz,2))*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z) - 
                256*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (Sqrt(1 + Power(lx,2) + Power(ly,2) + Power(lz,2))*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                256*(-(kx*p1x) - ky*p1y - kz*p1z + p10*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (Sqrt(1 + Power(lx,2) + Power(ly,2) + Power(lz,2))*p20 - lx*p2x - ly*p2y - lz*p2z) + 
                512*(-(kx*lx) - ky*ly - kz*lz + Sqrt(1 + Power(lx,2) + Power(ly,2) + Power(lz,2))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p1x) - ky*p1y - kz*p1z + p10*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(Sqrt(1 + Power(lx,2) + Power(ly,2) + Power(lz,2))*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                512*(-(kx*lx) - ky*ly - kz*lz + Sqrt(1 + Power(lx,2) + Power(ly,2) + Power(lz,2))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*(-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (Sqrt(1 + Power(lx,2) + Power(ly,2) + Power(lz,2))*p20 - lx*p2x - ly*p2y - lz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                512*(-(kx*lx) - ky*ly - kz*lz + Sqrt(1 + Power(lx,2) + Power(ly,2) + Power(lz,2))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-(kx*p1x) - ky*p1y - kz*p1z + p10*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) - 
                256*(Sqrt(1 + Power(lx,2) + Power(ly,2) + Power(lz,2))*p10 - lx*p1x - ly*p1y - lz*p1z)*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2))))/
                Abs(Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2)*
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) + 
                (2*(-1/(32.*Power(1 + Power(lx,2) + Power(ly,2) + Power(lz,2),1.5)) + 
                (3*(Sqrt(1 + Power(lx,2) + Power(ly,2) + Power(lz,2)) - p10 - p20 + 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))/
                (64.*Power(1 + Power(lx,2) + Power(ly,2) + Power(lz,2),2)))*
                (-256*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*p20*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))\
                + 256*p20*(-(kx*p1x) - ky*p1y - kz*p1z + 
                p10*(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))\
                + 256*p10*(-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z) + 
                512*(Power(p10,2) - Power(p1x,2) - Power(p1y,2) - Power(p1z,2))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z) - 
                1024*(-(kx*p1x) - ky*p1y - kz*p1z + p10*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z) - 
                256*p10*(-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*p20*(-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                512*(-(kx*p1x) - ky*p1y - kz*p1z + p10*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) + 
                512*(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (-(kx*p2x) - ky*p2y + p20*(p10 + p20 - 
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))) - kz*p2z)*
                (p10*p20 - p1x*p2x - p1y*p2y - p1z*p2z) - 
                256*p10*(-Power(kx,2) - Power(ky,2) - Power(kz,2) + 
                Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2))*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2)) + 
                512*(-(kx*p1x) - ky*p1y - kz*p1z + p10*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))*
                (p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)))*
                (Power(p20,2) - Power(p2x,2) - Power(p2y,2) - Power(p2z,2))))/
                Abs(Power(p10 + p20 - Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2)),2)*
                Sqrt(Power(kx - p1x - p2x,2) + Power(ky - p1y - p2y,2) + Power(kz - p1z - p2z,2))))/2.;
        }



        double dt_integrand(double kx, double ky, double kz, double lx, double ly, double lz, int tag){
            double e_k=Power(Power(kx,2)+Power(ky,2)+Power(kz,2),0.5);
            double e_kl=Power(Power(kx-lx,2)+Power(ky-ly,2)+Power(kz-lz,2),0.5);
            double e_l=Power(Power(lx,2)+Power(ly,2)+Power(lz,2),0.5);

            if(tag==1){
                double dt_virtual_left=dt_num(e_k, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*e_k*2*(e_l-p10-p20))*
                    (e_k-e_l-e_kl)*(e_k-e_l+e_kl)*(e_k-p10-p20+e_k)*(-p10-p20))+
                    dt_num(e_l+e_kl, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*(e_l-p10-p20)*2*(e_kl))*
                    (e_l+e_kl-e_k)*(e_l+e_kl+e_k)*(e_l+e_kl-p10-p20-e_k)*(e_l+e_kl-p10-p20+e_k))+
                    dt_num(p10+p20+e_k, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*(e_l-p10-p20)*2*(e_k))*
                    (p10+p20)*(p10+p20+e_k+e_k)*(p10+p20+e_k-e_l-e_kl)*(p10+p20+e_k-e_l+e_kl)
                    );

                return dt_virtual_left;
            };

            if(tag==2){
                double dt_virtual_right=dt_num(e_k, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_k*2*(e_k-p10-p20)*2*e_l)*
                    (e_k-e_l-e_kl)*(e_k-e_l+e_kl)*(-p10-p20)*(-p10-p20+2*e_l))+
                    dt_num(e_k, kx, ky, kz, e_k+e_kl, lx, ly, lz)/
                    (abs(2*e_k*2*(e_k-p10-p20)*2*(e_kl))*
                    (e_k+e_kl-p10-p20-e_l)*(e_k+e_kl-p10-p20+e_l)*(e_k+e_kl-e_l)*(e_k+e_kl+e_l))+
                    dt_num(e_k, kx, ky, kz, e_l+p10+p20, lx, ly, lz)/
                    (abs(2*e_k*2*(e_k-p10-p20)*2*(e_l))*
                    (p10+p20)*(p10+p20+2*e_l)*(e_k-e_l-p10-p20-e_kl)*(e_k-e_l-p10-p20+e_kl)
                    );

                return dt_virtual_right;
            };

            if(tag==3){
                double dt_real_right=dt_num(e_k, kx, ky, kz, -e_l+p10+p20, lx, ly, lz)/
                    (abs(2*e_k*2*(p20+p10-e_l-e_k)*2*(e_l))*
                    (-p10-p20)*(-p10-p20+2*e_k)*(p10+p20)*(p10+p20-2*e_l)
                    );

                return dt_real_right;
            };

            if(tag==4){
                double dt_real_left=dt_num(p10+p20-e_k, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*(-e_l-e_k+p10+p20)*2*(e_k))*
                    (-p10-p20)*(-p10-p20+2*e_l)*(p10+p20)*(p10+p20-2*e_k)
                    );

                return dt_real_left;
            };

            if(tag==5){
                double dt_left_uv=dt_num_UV_left(0, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*(e_l-p10-p20))
                    );

                return dt_left_uv;
            };

            if(tag==6){
                double dt_right_uv=dt_num_UV_right(e_k, kx, ky, kz, 0, lx, ly, lz)/
                    (abs(2*e_k*2*(e_k-p10-p20))
                    );

                return dt_right_uv;
            };

            return 0;
        }

        double t_channel_integrand(double kx, double ky, double kz, double lx, double ly, double lz, int tag){
            double e_k=Power(Power(kx,2)+Power(ky,2)+Power(kz,2),0.5);
            double e_kl=Power(Power(kx+lx,2)+Power(ky+ly,2)+Power(kz+lz,2),0.5);
            double e_l=Power(Power(lx,2)+Power(ly,2)+Power(lz,2),0.5);

            if(tag==1){
                double usual=t_channel_num(-e_l+e_kl, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*(-e_l+e_kl+p10+p20)*2*e_kl)*
                    pow((-e_l+e_kl+e_k)*(-e_l+e_kl-e_k),2.)
                    );

                return usual;
            };

            if(tag==2){
                double weird=-t_channel_num(e_k, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*(e_k+p10+p20)*4*pow(e_k,3.))*(e_k+e_l+e_kl)*(e_k+e_l-e_kl))+
                    t_channel_num_derivative(e_k, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*(e_k+p10+p20)*4*pow(e_k,2.))
                    );

                return weird;
            };

            return 0;
        }

        double u_channel_integrand(double kx, double ky, double kz, double lx, double ly, double lz, int tag){
            double e_k=Power(Power(kx,2)+Power(ky,2)+Power(kz,2),0.5);
            double e_kl=Power(Power(kx+lx,2)+Power(ky+ly,2)+Power(kz+lz,2),0.5);
            double e_l=Power(Power(lx,2)+Power(ly,2)+Power(lz,2),0.5);

            if(tag==1){
                double usual=u_channel_num(-e_l+e_kl, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*e_kl*2*(-e_l+e_kl+p10+p20))*
                    (-e_l+e_kl+e_k)*(-e_l+e_kl-e_k)*(p10+p20)*(2*e_kl+p10+p20)
                    );

                return usual;
            };

            if(tag==2){
                double weird1=u_channel_num(e_k, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*e_k*2*(e_k+p10+p20))*
                    (e_k+e_l+e_kl)*(e_k+e_l-e_kl)*(e_k+e_l+p10+p20+e_kl)*(e_k+e_l+p10+p20-e_kl)
                    );

                return weird1;
            };

            if(tag==3){
                double weird2=u_channel_num(-e_l+e_kl, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*e_kl*2*(e_kl+p10+p20))*
                    (-e_l+e_kl-e_k)*(-e_l+e_kl+e_k)*(-e_l+e_kl+p10+p20-e_k)*(-e_l+e_kl+p10+p20+e_k)
                    );

                return weird2;
            };

            if(tag==4){
                double weird3=u_channel_num(e_k, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_k*2*e_l*2*(e_k+e_l+p10+p20))*
                    (p10+p20)*(2*e_k+p10+p20)*(e_k+e_l-e_kl)*(e_k+e_l+e_kl)
                    );

                return weird3;
            };

            if(tag==5){
                double weird4=u_channel_num(-e_k, kx, ky, kz, e_k+e_kl, lx, ly, lz)/
                    (abs(2*e_k*2*e_kl*2*(e_kl+p10+p20))*
                    (p10+p20)*(p10+p20-2*e_k)*(e_k+e_kl-e_l)*(e_k+e_kl+e_l)
                    );

                return weird4;
            };

            if(tag==6){
                double weird5=u_channel_num(e_k, kx, ky, kz, -e_k-p10-p20+e_kl, lx, ly, lz)/
                    (abs(2*e_k*2*(e_k+p10+p20)*2*e_kl)*
                    (-p10-p20)*(-p10-p20+2*e_kl)*(-e_k-p10-p20+e_kl-e_l)*(-e_k-p10-p20+e_kl+e_l)
                    );

                return weird5;
            };
            
            return 0;
        };


        double s_channel_integrand(double kx, double ky, double kz, double lx, double ly, double lz){
            double e_k=Power(Power(kx,2)+Power(ky,2)+Power(kz,2),0.5);
            double e_kl=Power(Power(kx+lx,2)+Power(ky+ly,2)+Power(kz+lz,2),0.5);
            double e_l=Power(Power(lx,2)+Power(ly,2)+Power(lz,2),0.5);

            double usual=s_channel_num(e_k, kx, ky, kz, e_l, lx, ly, lz)/
                (abs(2*e_l*2*(e_k-e_l+p10+p20)*2*e_k)*
                pow((p10+p20)*(2*e_k+p10+p20),2.)
                );

            return usual;
        };

        double bub_integrand(double kx, double ky, double kz, double lx, double ly, double lz, int tag){
            double e_k=Power(Power(kx,2)+Power(ky,2)+Power(kz,2),0.5);
            double e_kl=Power(Power(kx+lx,2)+Power(ky+ly,2)+Power(kz+lz,2),0.5);
            double e_l=Power(Power(lx,2)+Power(ly,2)+Power(lz,2),0.5);

            if(tag==1){
                double virtual_cut=bub_num_derivative(0, kx, ky, kz, 0, lx, ly, lz);

                return virtual_cut;
            };

            if(tag==2){
                double real_cut=bub_num(p10+p20-e_k, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*(p10+p20-e_k-e_l)*2*e_k)*
                    pow((p10+p20)*(p10+p20-2*e_k),2.)
                    );

                return real_cut;
            };

            if(tag==3){
                double UV_cut=bub_UV_num(0, kx, ky, kz, 0, lx, ly, lz);

                return UV_cut;
            };

            return 0;
        };

        double st_channel_integrand(double kx, double ky, double kz, double lx, double ly, double lz, int tag){
            double e_k=Power(Power(kx,2)+Power(ky,2)+Power(kz,2),0.5);
            double e_kl=Power(Power(kx+lx,2)+Power(ky+ly,2)+Power(kz+lz,2),0.5);
            double e_l=Power(Power(lx,2)+Power(ly,2)+Power(lz,2),0.5);

            if(tag==1){
                double usual=st_channel_num(-e_k, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_k*2*(-e_k+e_l+p10+p20)*2*e_l)*
                    (p10+p20)*(p10+p20+2*e_l)*(p10+p20)*(p10+p20-2*e_k)
                );

                return usual;
            };

            if(tag==2){
                double weird=st_channel_num(-e_k, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*e_k*2*(-e_k+p10+p20))*
                    (-e_k+e_l+p10+p20-e_kl)*(-e_k+e_l+p10+p20+e_kl)*(p10+p20)*(p10+p20+2*e_l)
                    );

                return weird;
            };

            return 0;
        };

        double t_channel_integrand_2(double kx, double ky, double kz, double lx, double ly, double lz, int tag){
            double e_k=Power(Power(kx,2)+Power(ky,2)+Power(kz,2),0.5);
            double e_kl=Power(Power(kx+lx,2)+Power(ky+ly,2)+Power(kz+lz,2),0.5);
            double e_l=Power(Power(lx,2)+Power(ly,2)+Power(lz,2),0.5);

            if(tag==1){
                double usual=t_channel_num_2(-e_l+e_kl, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*(-e_l+e_kl+p10+p20)*2*e_kl)*
                    pow((-e_l+e_kl+e_k)*(-e_l+e_kl-e_k),2.)
                    );

                return usual;
            };

            if(tag==2){
                double weird=-t_channel_num_2(e_k, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*(e_k+p10+p20)*4*pow(e_k,3.))*(e_k+e_l+e_kl)*(e_k+e_l-e_kl))+
                    t_channel_num_derivative_2(e_k, kx, ky, kz, e_l, lx, ly, lz)/
                    (abs(2*e_l*2*(e_k+p10+p20)*4*pow(e_k,2.))
                    );

                return weird;
            };

            return 0;
        }

        double dt_LU(double kx, double ky, double kz, double lx, double ly, double lz){
            double t_v_left=f.flow1(kx,ky,kz,lx,ly,lz,1);
            double t_v_right=f.flow1(kx,ky,kz,lx,ly,lz,0);
            double t_r_left=f.flow2(kx,ky,kz,lx,ly,lz,1);
            double t_r_right=f.flow2(kx,ky,kz,lx,ly,lz,0);

            double exp_v_left=exp(-pow(t_v_left,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);
            double exp_v_right=exp(-pow(t_v_right,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);
            double exp_r_left=exp(-pow(t_r_left,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);
            double exp_r_right=exp(-pow(t_r_right,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);

            double j_v_left=pow(t_v_left,6)/f.jacques1(kx,ky,kz,lx,ly,lz,1);
            double j_v_right=pow(t_v_right,6)/f.jacques1(kx,ky,kz,lx,ly,lz,0);
            double j_r_left=pow(t_r_left,6)/f.jacques2(kx,ky,kz,lx,ly,lz,1);
            double j_r_right=pow(t_r_right,6)/f.jacques2(kx,ky,kz,lx,ly,lz,0);

            return constant*(j_v_left*exp_v_left*dt_integrand(t_v_left*kx,t_v_left*ky,t_v_left*kz,t_v_left*lx,t_v_left*ly,t_v_left*lz,1)+
                      j_v_left*exp_v_left*dt_integrand(t_v_left*kx,t_v_left*ky,t_v_left*kz,t_v_left*lx,t_v_left*ly,t_v_left*lz,5)+
                      j_v_right*exp_v_right*dt_integrand(t_v_right*kx,t_v_right*ky,t_v_right*kz,t_v_right*lx,t_v_right*ly,t_v_right*lz,2)+
                      j_v_right*exp_v_right*dt_integrand(t_v_right*kx,t_v_right*ky,t_v_right*kz,t_v_right*lx,t_v_right*ly,t_v_right*lz,6)+
                      j_r_left*exp_r_left*dt_integrand(t_r_left*kx,t_r_left*ky,t_r_left*kz,t_r_left*lx,t_r_left*ly,t_r_left*lz,3)+
                      j_r_right*exp_r_right*dt_integrand(t_r_right*kx,t_r_right*ky,t_r_right*kz,t_r_right*lx,t_r_right*ly,t_r_right*lz,4)
                      );

        }

        double t_channel_LU(double kx, double ky, double kz, double lx, double ly, double lz){
            double t_v=f.flow3(kx,ky,kz,lx,ly,lz);
            double t_r=f.flow1(kx,ky,kz,lx,ly,lz,0);

            double exp_v=exp(-pow(t_v,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);
            double exp_r=exp(-pow(t_r,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);

            double j_v=pow(t_v,6)/f.jacques3(kx,ky,kz,lx,ly,lz);
            double j_r=pow(t_r,6)/f.jacques1(kx,ky,kz,lx,ly,lz,0);

            return constant*(j_v*exp_v*t_channel_integrand(t_v*kx,t_v*ky,t_v*kz,t_v*lx,t_v*ly,t_v*lz,1)+
                    j_r*exp_r*t_channel_integrand(t_r*kx,t_r*ky,t_r*kz,t_r*lx,t_r*ly,t_r*lz,2)
                    );
        }

        double u_channel_LU(double kx, double ky, double kz, double lx, double ly, double lz){
            double t_u=f.flow3(kx,ky,kz,lx,ly,lz);
            double t_w1=f.flow1(kx,ky,kz,lx,ly,lz,0);
            double t_w2=f.flow4(kx,ky,kz,lx,ly,lz);
            double t_w3=f.flow5(kx,ky,kz,lx,ly,lz);
            double t_w4=f.flow4(kx,ky,kz,lx,ly,lz);
            double t_w5=f.flow1(kx,ky,kz,lx,ly,lz,0);

            double exp_u=exp(-pow(t_u,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);
            double exp_w1=exp(-pow(t_w1,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);
            double exp_w2=exp(-pow(t_w2,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);
            double exp_w3=exp(-pow(t_w3,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);
            double exp_w4=exp(-pow(t_w4,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);
            double exp_w5=exp(-pow(t_w5,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);

            double j_u=pow(t_u,6)/f.jacques3(kx,ky,kz,lx,ly,lz);
            double j_w1=pow(t_w1,6)/f.jacques1(kx,ky,kz,lx,ly,lz,0);
            double j_w2=pow(t_w2,6)/f.jacques4(kx,ky,kz,lx,ly,lz);
            double j_w3=pow(t_w3,6)/f.jacques5(kx,ky,kz,lx,ly,lz);
            double j_w4=pow(t_w4,6)/f.jacques4(kx,ky,kz,lx,ly,lz);
            double j_w5=pow(t_w5,6)/f.jacques1(kx,ky,kz,lx,ly,lz,0);

            return constant*(j_u*exp_u*u_channel_integrand(t_u*kx,t_u*ky,t_u*kz,t_u*lx,t_u*ly,t_u*lz,1)+
                    j_w1*exp_w1*u_channel_integrand(t_w1*kx,t_w1*ky,t_w1*kz,t_w1*lx,t_w1*ly,t_w1*lz,2)+
                    j_w2*exp_w2*u_channel_integrand(t_w2*kx,t_w2*ky,t_w2*kz,t_w2*lx,t_w2*ly,t_w2*lz,3)+
                    j_w3*exp_w3*u_channel_integrand(t_w3*kx,t_w3*ky,t_w3*kz,t_w3*lx,t_w3*ly,t_w3*lz,4)+
                    j_w4*exp_w4*u_channel_integrand(t_w4*kx,t_w4*ky,t_w4*kz,t_w4*lx,t_w4*ly,t_w4*lz,5)+
                    j_w5*exp_w5*u_channel_integrand(t_w5*kx,t_w5*ky,t_w5*kz,t_w5*lx,t_w5*ly,t_w5*lz,6)
                    );
        }

        double s_channel_LU(double kx, double ky, double kz, double lx, double ly, double lz){
            double t=f.flow7(kx,ky,kz,lx,ly,lz);

            double exp_s=exp(-pow(t,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);

            double j=pow(t,6)/f.jacques7(kx,ky,kz,lx,ly,lz);

            return constant*j*exp_s*s_channel_integrand(t*kx,t*ky,t*kz,t*lx,t*ly,t*lz);
        }

        double bub_LU(double kx, double ky, double kz, double lx, double ly, double lz){
            double t_v=f.flow1(kx,ky,kz,lx,ly,lz,0);
            double t_r=f.flow2(kx,ky,kz,lx,ly,lz,1);

            double exp_v=exp(-pow(t_v,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);
            double exp_r=exp(-pow(t_r,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);

            double j_v=pow(t_v,6)/f.jacques1(kx,ky,kz,lx,ly,lz,0);
            double j_r=pow(t_r,6)/f.jacques2(kx,ky,kz,lx,ly,lz,1);

            return constant*(j_v*exp_v*bub_integrand(t_v*kx,t_v*ky,t_v*kz,t_v*lx,t_v*ly,t_v*lz,1)+
                    j_v*exp_v*bub_integrand(t_v*kx,t_v*ky,t_v*kz,t_v*lx,t_v*ly,t_v*lz,3)+
                    j_r*exp_r*bub_integrand(t_r*kx,t_r*ky,t_r*kz,t_r*lx,t_r*ly,t_r*lz,2)
                    );
        }

        double st_channel_LU(double kx, double ky, double kz, double lx, double ly, double lz){
            double t_u=f.flow3(kx,ky,kz,lx,ly,lz);
            double t_w=f.flow1(kx,ky,kz,lx,ly,lz,0);

            double exp_u=exp(-pow(t_u,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);
            double exp_w=exp(-pow(t_w,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);

            double j_u=pow(t_u,6)/f.jacques3(kx,ky,kz,lx,ly,lz);
            double j_w=pow(t_w,6)/f.jacques1(kx,ky,kz,lx,ly,lz,0);

            return constant*(j_u*exp_u*st_channel_integrand(t_u*kx,t_u*ky,t_u*kz,t_u*lx,t_u*ly,t_u*lz,1)+
                    j_w*exp_w*st_channel_integrand(t_w*kx,t_w*ky,t_w*kz,t_w*lx,t_w*ly,t_w*lz,2)
                    );
        }

        double t_channel_2_LU(double kx, double ky, double kz, double lx, double ly, double lz){
            double t_u=f.flow3(kx,ky,kz,lx,ly,lz);
            double t_w=f.flow1(kx,ky,kz,lx,ly,lz,0);

            double exp_u=exp(-pow(t_u,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);
            double exp_w=exp(-pow(t_w,2.)/pow(sigma,2))/(Sqrt(2*acos(0.0))*sigma);

            double j_u=pow(t_u,6)/f.jacques3(kx,ky,kz,lx,ly,lz);
            double j_w=pow(t_w,6)/f.jacques1(kx,ky,kz,lx,ly,lz,0);
            
            return constant*(j_u*exp_u*t_channel_integrand_2(t_u*kx,t_u*ky,t_u*kz,t_u*lx,t_u*ly,t_u*lz,1)+
                    j_w*exp_w*t_channel_integrand_2(t_w*kx,t_w*ky,t_w*kz,t_w*lx,t_w*ly,t_w*lz,2)
                    );
        }



};






int main(){

double p1[4]={1,0,0,1};
double p2[4]={1,0,0,-1};

flows f=flows(p1,p2);

cout<<f.flow1(0.1,0.23,0.3,0.37,0.5,0.8,0)<<endl;
cout<<f.flow1(0.1,0.23,0.3,0.37,0.5,0.8,1)<<endl;
cout<<f.flow2(0.1,0.23,0.3,0.37,0.5,0.8,0)<<endl;
cout<<f.flow2(0.1,0.23,0.3,0.37,0.5,0.8,1)<<endl;
cout<<f.flow3(0.1,0.23,0.3,0.37,0.5,0.8)<<endl;
cout<<f.flow4(0.1,0.23,0.3,0.37,0.5,0.8)<<endl;
cout<<f.flow5(0.1,0.23,0.3,0.37,0.5,0.8)<<endl;
cout<<f.flow6(0.1,0.23,0.3,0.37,0.5,0.8)<<endl;
cout<<f.flow7(0.1,0.23,0.3,0.37,0.5,0.8)<<endl;

cout<<"----"<<endl;

cout<<f.jacques1(0.1,0.23,0.3,0.37,0.5,0.8,0)<<endl;
cout<<f.jacques1(0.1,0.23,0.3,0.37,0.5,0.8,1)<<endl;
cout<<f.jacques2(0.1,0.23,0.3,0.37,0.5,0.8,0)<<endl;
cout<<f.jacques2(0.1,0.23,0.3,0.37,0.5,0.8,1)<<endl;
cout<<f.jacques3(0.1,0.23,0.3,0.37,0.5,0.8)<<endl;
cout<<f.jacques4(0.1,0.23,0.3,0.37,0.5,0.8)<<endl;
cout<<f.jacques5(0.1,0.23,0.3,0.37,0.5,0.8)<<endl;
cout<<f.jacques6(0.1,0.23,0.3,0.37,0.5,0.8)<<endl;
cout<<f.jacques7(0.1,0.23,0.3,0.37,0.5,0.8)<<endl;

integrands integrand(p1,p2,1,f,1);

cout<<"-----"<<endl;

cout<<integrand.dt_num(0.1,0.23,0.3,0.37,0.5,0.8,1.1,0.2)<<endl;
cout<<integrand.dt_num_UV_left(0.1,0.23,0.3,0.37,0.5,0.8,1.1,0.2)<<endl;
cout<<integrand.dt_num_UV_right(0.1,0.23,0.3,0.37,0.5,0.8,1.1,0.2)<<endl;
cout<<integrand.dt_num_UV_right(0.1,0.23,0.3,0.37,0.5,0.8,1.1,0.2)<<endl;
cout<<integrand.t_channel_num(0.1,0.23,0.3,0.37,0.5,0.8,1.1,0.2)<<endl;
cout<<integrand.t_channel_num_derivative(0.1,0.23,0.3,0.37,0.5,0.8,1.1,0.2)<<endl;
cout<<integrand.u_channel_num(0.1,0.23,0.3,0.37,0.5,0.8,1.1,0.2)<<endl;
cout<<integrand.s_channel_num(0.1,0.23,0.3,0.37,0.5,0.8,1.1,0.2)<<endl;
cout<<integrand.st_channel_num(0.1,0.23,0.3,0.37,0.5,0.8,1.1,0.2)<<endl;
cout<<integrand.bub_num(0.1,0.23,0.3,0.37,0.5,0.8,1.1,0.2)<<endl;

cout<<"-----"<<endl;
cout<<integrand.dt_integrand(0.23,0.3,0.37,0.8,1.1,0.2,1)<<endl;
cout<<integrand.dt_integrand(0.23,0.3,0.37,0.8,1.1,0.2,2)<<endl;
cout<<integrand.dt_integrand(0.23,0.3,0.37,0.8,1.1,0.2,3)<<endl;
cout<<integrand.dt_integrand(0.23,0.3,0.37,0.8,1.1,0.2,4)<<endl;
cout<<integrand.dt_integrand(0.23,0.3,0.37,0.8,1.1,0.2,5)<<endl;
cout<<integrand.dt_integrand(0.23,0.3,0.37,0.8,1.1,0.2,6)<<endl;

cout<<"-----"<<endl;
cout<<integrand.t_channel_integrand(0.23,0.3,0.37,0.8,1.1,0.2,1)<<endl;
cout<<integrand.t_channel_integrand(0.23,0.3,0.37,0.8,1.1,0.2,2)<<endl;

cout<<"-----"<<endl;
cout<<integrand.u_channel_integrand(0.23,0.3,0.37,0.8,1.1,0.2,1)<<endl;
cout<<integrand.u_channel_integrand(0.23,0.3,0.37,0.8,1.1,0.2,2)<<endl;
cout<<integrand.u_channel_integrand(0.23,0.3,0.37,0.8,1.1,0.2,3)<<endl;
cout<<integrand.u_channel_integrand(0.23,0.3,0.37,0.8,1.1,0.2,4)<<endl;
cout<<integrand.u_channel_integrand(0.23,0.3,0.37,0.8,1.1,0.2,5)<<endl;
cout<<integrand.u_channel_integrand(0.23,0.3,0.37,0.8,1.1,0.2,6)<<endl;

cout<<"-----"<<endl;
cout<<integrand.s_channel_integrand(0.23,0.3,0.37,0.8,1.1,0.2)<<endl;

cout<<"-----"<<endl;
cout<<integrand.bub_integrand(0.23,0.3,0.37,0.8,1.1,0.2,1)<<endl;
cout<<integrand.bub_integrand(0.23,0.3,0.37,0.8,1.1,0.2,2)<<endl;

cout<<"-----"<<endl;
cout<<integrand.st_channel_integrand(0.23,0.3,0.37,0.8,1.1,0.2,1)<<endl;
cout<<integrand.st_channel_integrand(0.23,0.3,0.37,0.8,1.1,0.2,2)<<endl;

cout<<"-----"<<endl;
cout<<integrand.t_channel_integrand_2(0.23,0.3,0.37,0.8,1.1,0.2,1)<<endl;
cout<<integrand.t_channel_integrand_2(0.23,0.3,0.37,0.8,1.1,0.2,2)<<endl;

cout<<"-----"<<endl;
cout<<integrand.dt_LU(0.23,0.3,0.37,0.8,1.1,0.2)<<endl;
cout<<integrand.t_channel_LU(0.23,0.3,0.37,0.8,1.1,0.2)<<endl;
cout<<integrand.u_channel_LU(0.23,0.3,0.37,0.8,1.1,0.2)<<endl;
cout<<integrand.s_channel_LU(0.23,0.3,0.37,0.8,1.1,0.2)<<endl;
cout<<integrand.bub_LU(0.23,0.3,0.37,0.8,1.1,0.2)<<endl;
cout<<integrand.st_channel_LU(0.23,0.3,0.37,0.8,1.1,0.2)<<endl;
cout<<integrand.t_channel_2_LU(0.23,0.3,0.37,0.8,1.1,0.2)<<endl;


return 0;

}
