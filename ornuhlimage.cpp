
//
// single run modified tomlinson model
//

// includes

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "ornuhlimage.h"

// Lazy stuff

using namespace std;

typedef unsigned int uint;

// objects

int main()
{
	//input values
	double barref = 5.0e-20; // 4.0e-20
	double kapparef = 0.0612245;
	
	double spring1 = 3.0 * 1e0;			// 2.0 //1e-4 for harmonic
	//double spring1 = atof(argv[1]);			// 2.0 //1e-4 for harmonic
	double spring2 = 1.0 * 1e0;	    	// 2.0 //1e-4 for harmonic
	double supvel = 1.0;				// 1.0
	double latcon = 2.5e-10;			// 2.5e-10
	double barr1 = barref;				// 
	double barr2 = 0.5 * barref;
	double kappa1 = kapparef;
	double kappa2 = 0.5 * kapparef;
	double align = 1.0;
	double nu2 = 0.382653;
	double nu4 = 5.809e17;
	
	double xmass = 5e-24;
	double xdamp = 1.0e14; //1e-9

	double lmass = 3.0e-24;
	double ldamp = 1.0e13;
    double t0 = xdamp/spring1; 
    //double t0 = 2e-11;
    //double t0 = atof(argv[2]);
	
	double temp = 300; // 1e20
	
	double modif = 0.5;	// 0.5 gives reliable time step dep. 1.0 should be ok
    double tstep = modif * 3e-14;	
	uint tsteps = 1.0/modif * 4.0e4;	// has to be even beucasue lazyness

    vector <uint> antsteps;
    ifstream reads ("./antimes.csv");
    
    uint shift = round(0.025e-9 / tstep);
    if (reads.is_open())
    {
        uint atstep;

        while(reads >> atstep)
        {
            if (atstep > shift)
            {antsteps.push_back(atstep-shift);}
        }
    }

	uint ttoa = ceil(latcon/(tstep));	// timesteps to minima

	string lfile = "lout.csv";
	string tomfile = "tomout.csv";
    string pfile = "parameters.dat";
    string tfile = "time.csv";

    uint ncores = 1;
    uint nfiles = 10; 

    uint runs = 10;

        #pragma omp parallel
        {
            #pragma omp for nowait
    for (uint m = 0; m < nfiles; m++) // yes, I'm reversing my loop index concention here... shit happens when you party nude
    {
	    string xfile = "xout" + to_string(m+0) + ".csv";
	    //string tfile = "time" + to_string(m) + ".csv";
        
        vector <double> publicpos;
        vector <double> publicf;
        publicpos.reserve(runs*tsteps);

            vector <double>privatepos;
            vector <double>privatef;
            if (runs > ncores)
                privatepos.reserve(tsteps*ceil((double)runs/ncores));
            else
                privatepos.reserve(tsteps);

            for (uint l = 0; l < runs; l++)
            {
                    tomlin afm(spring1,spring2,supvel,latcon,align,barr1,barr2,kappa1,
                           kappa2,nu2,nu4,temp,tstep,tsteps,xmass,lmass,xdamp,ldamp,
                           t0, tfile,xfile,lfile,tomfile,"",pfile);
   
                afm.setposx(0e-9); 
                for ( uint k = 0; k < tsteps; k++ )
                    {
                    	//if (k == 20.0*ttoa)
                    	//{	
                    	//	afm.tpause();
                    	//	cout << "resuming at tstep / t = " << k << " / " << k*tstep << endl;
                    //}
                    	//	
                    	//if (k == 0.25*tsteps)
                    //{
                    	//	afm.treverse();
                    	//	cout << "reversing at tstep / t = " << k << " / " << k*tstep << endl;
                    //}
                    	//if (k == 0.7*tsteps) //25 and 70 gives more or less match up
                    //{
                    	//	afm.treverse();
                    	//	cout << "reversing at tstep / t = " << k << " / " << k*tstep << endl;
                    //}
                    	//if (k == 20*ttoa)
                    	//{
                    	//	cout << "pausing at tstep / t = " << k << " / " << k*tstep << endl;
                    	//	afm.tpause();
                    	//}
                    	//else if (k == 0.9*tsteps)
                    	//{
                    	//	afm.tpause();
                    	//	cout << "resuming at tstep / t = " << k << " / " << k*tstep << endl;
                    	//}

                    	afm.rk4();
                    	//afm.calckin();
                    	//afm.calcpot();
                        afm.calcfric();

                    	//afm.pushvals();
                    	afm.pushvalslite();
                    	afm.inctime();
                    }

                auto posptr = afm.getposs();
                privatepos.insert(privatepos.end(),posptr->begin(),posptr->end());
                //auto fptr = afm.getfrics();
                //privatef.insert(privatef.end(),fptr->begin(),fptr->end());
                //afm.writedata();    // REMOVE WHEN PARALLEL
            }
	    // pragma omp critical
            publicpos.insert(publicpos.end(),privatepos.begin(),privatepos.end());
            //publicf.insert(publicf.end(),privatef.begin(),privatef.end());
        
        //afm.printins();
        //afm.printouts();
        //afm.printavgs();
    
        //// write at significant x incrrease
        //vector <double> ts;
        //vector <double> xs; 
        //ts.reserve(tsteps + runs/10);
        //xs.reserve(tsteps + runs/10);
        
        //uint d0 = 10;
        //int shift = 1;

        //for (uint k = 0; k < runs; k++) 
        //{
        //    uint t0 = k*tsteps;
        //    double xscale = *max_element(publicpos.begin()+t0, publicpos.begin() + (k+1)*tsteps - 1);
        //    double threshold = 1.5e-3 * xscale;


        //    for (uint l = 0; l < tsteps; l += d0 + shift)
        //    {
        //        uint next = l + 2*shift;

        //        if (t0 + next > publicpos.size()) // out of bounds check
        //            break;
        //
        //        double x1 = publicpos[t0+l];
        //        double x2 = publicpos[t0+next];

        //        //cout << l << " " << next << " " << next - l << endl 
        //        //     << x1 << " " << x2 << " " << x2 - x1 << endl;

        //        if (abs(x2 - x1) > threshold)
        //        {
        //            ts.push_back(t0+l);
        //            xs.push_back(publicpos[t0+l]);
        //            
        //            shift = 1;
        //        }
        //        else
        //        {
        //            shift *= 2;
        //        }
        //    }
        //}
        
        //ofstream xstream, tstream ;
        
        //tstream.open(tfile);
        //for (uint k = 0; k < ts.size() - 1; k++)
        //    tstream << ts[k] << endl;

        //xstream.open(xfile);
        //for (uint k = 0; k < xs.size() - 1; k++)
        //    xstream << xs[k] << endl
        //            << xs[k+1] << endl;
        //xstream.close();
        //tstream.close();
        
        //// write at constant distance
        
        uint writeat = 500;
        vector <double> times;
        times.reserve(round(tsteps/writeat));

        for (uint k = 0; k < tsteps; k++)
            times.push_back(k*tstep);

        ofstream xstream, tstream;

        tstream.open(tfile);
        for (uint k = 0; k < times.size() - 1; k++)	// -1 because consecutive points written, avoids out of bound
        {
            if (k % writeat == 0) 
            {
                tstream << setprecision(16) << times[k] << endl << times[k+1] << endl;
            }
        }

        xstream.open(xfile);
        for (uint k = 0; k < publicpos.size() - 1; k++)
        {
            if (k % writeat == 0) 
            {
                xstream << setprecision(16) << publicpos[k] << endl << publicpos[k+1] << endl;
            }
        }

        tstream.close();
        xstream.close();
       
        // write slipping points 

        //vector <uint> slips;
        //uint slipmode = 1;
        //uint adjacent = 10;
        //uint ttoa = ceil(latcon/(tstep*supvel)); // timesteps between lattice spacings
        ////uint stride = 1000;
        //uint stride = ttoa/4.0;
        //findslips(slipmode,adjacent,stride,&times,&publicf,&slips);
        //
        //ofstream sstream;

        //sstream.open("slips.dat");
        //for (auto &k : slips)
        //    sstream << k*tstep << endl;
       
        //sstream.close();

        //// write at predefined grid
        
        //ofstream xstream; //tstream; 
        
        //vector <double> times;
        //for (auto &el : antsteps)
        //    times.push_back(el*tstep);
	//
	//tstream.open(tfile);
        //for (auto &el : times)
        //    tstream << el << endl; 

        //xstream.open(xfile);
	//for (uint k = 0; k < runs; k++)
	//{
        //for (auto &el : antsteps)
        //    xstream << setprecision(16) << publicpos[el+k*tsteps] << endl << publicpos[el+1+k*tsteps] << endl;
	//}
        ////tstream.close();
        //xstream.close();
        }
    }
}

	
