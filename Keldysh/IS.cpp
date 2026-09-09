#include "SFS.h"
 #include <complex>
 #include <iostream>
 #define pi 3.141592653589793
 #define icom (complex<double>(0, 1.))

 using namespace std;

  double E_IS_calc(double *Is, complex<double> *G,complex<double> *Del, double*q, double I)
  {
   // double  Xi1;
    double w, h, dGmax;
    complex<double> t1,t2,wm, Fx, F1x, Gx, G1x;
    double ro, D, E;
	ro = 1; E = 0;

     complex<double> *Fi, *Fi1, *Fi_old, *Ei, *F, *F1, *G1;
     Fi=new complex<double>[N];
     Fi1=new complex<double>[N];
     Fi_old = new complex<double>[N];
	 F = new complex<double>[N];
	 F1 = new complex<double>[N];
	 G1 = new complex<double>[N];
	 Ei = new complex<double>[N];

	 D = 2 * pi * 1 * 1;
	 for (int i = 1; i < N - 1; i++)
	 {
		 Is[i] = 0;
		 Ei[i] = abs(Del[i])*abs(Del[i]) * log(T);
		//cout <<fixed<< Del[i] << " "<<Is[i]<< " "<<Ei[i]<<endl;
	 }

    for(int iw=0; iw<w_obrez; iw++)
    {
            w=pi*T*(2*iw+1);
      //  wm=w+ icom*H;
        for(int i=0; i<N; i++) G[i]=1.+icom*0.1;

        dGmax=1;
         while (dGmax > epsG)
            {
            // calc Fi(w), Fi(-w)
			 for (int i = 0; i<N; i++)
				 G1[i] = -conj(G[i]);
                Prog( Fi, G, Del, w, q, I);
                Prog( Fi1, G1, Del, -w, q, I);
            // recalc G
                Gcalc( G, &dGmax, Fi, Fi1, Del, w,q);
            }
            //for (int i=0; i<N; i++) G1[i]= -conj(G[i]);
            //cout<< w<<" "<<dGmax<<" "<<G[N-N_S2]<<" "<<G1[N-N_S2]<<" "<<Fi[N-N_S2]<<" "<<Fi1[N-N_S2]<<endl;

		 // E calculation

		 //for (int i = 1; i < N - 1; i++)
		 //{
			// wm = get_wm(i, w);
			// F[i] = Fi[i] / sqrt(wm*wm + Fi[i] * conj(Fi1[i]));
			// wm = get_wm(i, -w);
			// F1[i] = Fi1[i] / sqrt(wm*wm + Fi1[i] * conj(Fi[i]));
		 //}

		 //for (int i = 1; i < N - 2; i++)
			//// if (Layer(i)==0)
			// {
			//	 h = get_h_out(i);
			//	 if (h != 0)
			//	 {
			//		 Fx = (F[i + 1] - F[i]) / h;
			//		 F1x = (F1[i + 1] - F1[i]) / h;
			//		 Gx = (G[i + 1] - G[i]) / h;
			//		 G1x = (G1[i + 1] - G1[i]) / h;
			//		 wm = get_wm(i, w);
			//		 Ei[i] = Ei[i] + pi* T*(abs(Del[i])*abs(Del[i]) / abs(w) - Del[i] * conj(F[i]) - conj(Del[i])*F[i] + D / 2.*(Fx * conj(F1x) + Gx*Gx) + 2. * wm*(1. - G[i]));
			//		 wm = get_wm(i, -w);
			//		 Ei[i] = Ei[i] + pi* T*(abs(Del[i])*abs(Del[i]) / abs(w) - Del[i] * conj(F1[i]) - conj(Del[i])*F1[i] + D / 2.*(F1x * conj(Fx) + G1x*G1x) + 2. * wm*(1. - G1[i]));
			//		// if (iw == 0)
			//		//	 cout << i << " GGG " << Fx << " " << F1x << " " << Gx << " " << G1x << endl;
			//	 }
			// }

		 // I calculation

        for (int i=1; i<N-2; i++)
            {     //wm= get_wm(i,w);
                  h=get_h(i);
				  ro = Roi[Layer(i)];
				  wm = get_wm(i, w);
				  t1 = (wm*wm + (Fi[i])*conj(Fi1[i]));
                  Is[i]+=real(-icom*T*pi/h*(conj(Fi1[i]))*(Fi[i+1]-Fi[i])/t1/ro);
				  wm = get_wm(i, -w);
                  t1 = (wm*wm+(Fi1[i])*conj(Fi[i]));
                  Is[i]+=real(-icom*T*pi/h*(conj(Fi[i]))*(Fi1[i+1]-Fi1[i])/t1/ro);
                //cout<<i<<"  "<<get_wm(i,w)<<endl;
            }

        cout<<w<<"  "<<Is[100]<<endl<<endl<<endl;

        }

	 //double E0, E1, E2;
	 //E0 = 0; E1 = 0; E2 = 0;
	 //E = 0;
	 //for (int i = 1; i < N - 2; i++)
	 //{
		// if (Layer(i) == 0)
		//	 E0 = E0 + real(Ei[i]) * get_h_out(i);
		// if (Layer(i) == 0)
		//	 E1 = E1 + imag(Ei[i]) * get_h_out(i);
		//// if (Layer(i) == 2)
		////	 E2 = E2 + real(Ei[i]) * get_h_out(i);
		// E=E+ real(Ei[i]) * get_h_out(i);
	 //}
	// cout << E0 << " " << E1 << " " << E2 << " " << E0+E1+E2 << " " << endl;
	 return(E);

      delete [] Fi;
      delete [] Fi1;
	  delete[] F;
	  delete[] F1;
	  delete[] G1;
	  delete[] Ei;

  }

