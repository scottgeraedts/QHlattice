#include "weierstrass.h"
weierstrass::weierstrass(double theta0, double alpha0, int N, vector<complex<double> > zero){
	pi=3.1519265359;
    ii=complex<double>(0,1); pre = pow(10, -10);
    zeros=zero; Nphi=N; alpha=alpha0; theta=theta0;
    omega1=0.5*sqrt(1.*pi*N/sin(theta0))*alpha0;
    omega2=0.5*sqrt(1.*pi*N/sin(theta0))/alpha0*exp(theta0*complex<double>(0,1));
    
    run_weierstrass();
}
weierstrass::weierstrass(){
	pi=3.1519265359;
}	
weierstrass::weierstrass(complex<double> w1, complex<double> w2, vector<complex<double> > zero){
	pi=3.1519265359;
    ii=complex<double>(0,1); pre = pow(10, -15);
    int make_int(double);
    zeros=zero;
    Nphi=make_int( real(2.*(conj(w1)*w2-w1*conj(w2))/pi/ii) );
    alpha=sqrt(norm(2.*w1)/norm(2.*w2));
    theta=real(-ii*log(alpha*alpha*w2/w1));
    omega1=w1; omega2=w2;
    
    run_weierstrass();
}
void weierstrass::run_weierstrass(){
    tau=omega2/omega1;
    geteta_1(); geteta_2();
    Gbar=(weta1*conj(omega2)-weta2*conj(omega1))/(-0.5*pi*ii*double(Nphi));
}
void weierstrass::geteta_1(){
    complex<double> eta_tmp, eta, err;
    eta=complex<double> (1.0/6.0, 0.0);
    err=complex<double> (1.0,1.0);
    int n=1;
    while (real(err)>pre || imag(err)>pre || real(err)<-pre || imag(err)<-pre) {
        eta_tmp=eta;
        eta+=1.0/pow(sin((double)n*pi*tau), 2);
        n++;
        err=eta-eta_tmp;
    }
    eta=eta*pi*pi/complex<double>(2,0)/omega1;
    weta1=eta;
}
void weierstrass::geteta_2(){
    complex<double> eta_tmp, eta, err;
    eta=complex<double> (1.0/6.0, 0.0);
    err=complex<double> (1.0,1.0);
    int n=1;
    while (real(err)>pre || imag(err)>pre || real(err)<-pre || imag(err)<-pre) {
        eta_tmp=eta;
        eta+=1.0/pow(sin((double)n*pi/tau), 2);
        n++;
        err=eta-eta_tmp;
    }
    eta=eta*pi*pi/complex<double>(2,0)/omega2;
    weta2=eta;
}
complex<double> weierstrass::wsigma(complex<double> z){
    complex<double> sigma_tmp, sigma, err, Oij1, Oij2;
    complex<double> omega_1, omega_2;
    
    bool optimizetau=true;
    if (optimizetau==true) {
        int sl2z[2][2];
        sl2z[0][0]=1;
        sl2z[0][1]=10;
        sl2z[1][0]=0;
        sl2z[1][1]=1;
        omega_1 = (double)sl2z[0][0]*omega1+(double)sl2z[0][1]*omega2;
        omega_2 = (double)sl2z[1][0]*omega1+(double)sl2z[1][1]*omega2;
    }
    else if (optimizetau==false){
        omega_1 = omega1;
        omega_2 = omega2;
    }
    
    sigma=z; err=sigma+ii;
    int n=1;
    while (norm(err)>pre+1 || norm(err)<1-pre) {
        sigma_tmp=sigma;
        for (int i=-n; i<=n; i++) {
            Oij1=2.0*i*omega_1+2.0*(+n)*omega_2;
            Oij2=2.0*i*omega_1+2.0*(-n)*omega_2;
            sigma*=(1.0-z/Oij1)*exp(z/Oij1+z*z/2.0/Oij1/Oij1);
            sigma*=(1.0-z/Oij2)*exp(z/Oij2+z*z/2.0/Oij2/Oij2);
        }
        for (int i=-n+1; i<=n-1; i++) {
            Oij1=2.0*(+n)*omega_1+2.0*i*omega_2;
            Oij2=2.0*(-n)*omega_1+2.0*i*omega_2;
            sigma*=(1.0-z/Oij1)*exp(z/Oij1+z*z/2.0/Oij1/Oij1);
            sigma*=(1.0-z/Oij2)*exp(z/Oij2+z*z/2.0/Oij2/Oij2);
        }
        n++;

        err=sigma/sigma_tmp;
    }
    return sigma;
}
complex<double> weierstrass::wsigma_gaussian(complex<double>(z)){
    complex<double> wf;
    wf = complex<double> (1,0);
    for (int i=0; i<Nphi; i++) {
        wf*=wsigma(z-zeros[i])*exp(-0.5*pow(z-zeros[i], 2)*Gbar/(1.));
    }
    wf*=exp(-0.5*norm(z));
    
    return wf;
}


int make_int(double x){
    int y=(int)x;
    
    if (abs(x-y)<0.5) {
        y=y;
    }
    else if (abs(x-(y+1))<0.5) {
        y++;
    }
    else if(abs(x-(y-1))<0.5){
        y--;
    }
    return y;
}

