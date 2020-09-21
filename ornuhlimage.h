
//
// header file for modified tomlinson model
//


// includes

#include <fstream>
#include <random>
#include <numeric> 
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm> 
//#include <algorithm>

// Lazy stuff

using namespace std;

typedef unsigned int uint;

// constants

static const double pi = atan(1.0)*4;
static const double kB = 1.38064852e-23;

// objects

class tomlin
{
	// sub classes
	class xobj
	{
		// dynamic variables
		double pos = 0;
		double vel = 0;
		double acc = 0;
		double kin = 0;

		// parameters
		const double mass;
		const double damp;

		// derived constants
		const double rmass = 1.0/mass;

		// containers
		vector <double> poss;
		vector <double> vels;
		vector <double> accs;

		// functions
		void calckin(){kin = 0.5*mass*pow(vel,2);}

		friend tomlin;

		xobj (double mass, double damp, uint tsteps)
			: mass(mass), damp(damp)
		{
			poss.reserve(tsteps);
			vels.reserve(tsteps);
			accs.reserve(tsteps);
		}
	};

	class lobj
	{
		// dynamic variables
		double pos = 0;
		double vel = 0;
		double acc = 0;
		double kin = 0;

		// parameters
		const double mass;
		const double damp;
        const double t0;

		// derived constants
		const double rmass = 1.0/mass;
		const double rt0 = 1.0/t0;

		// containers
		vector <double> poss;
		vector <double> vels;
		vector <double> accs;
		vector <double> slips;

		// functions
		void calckin(){kin = 0.5*mass*pow(vel,2);}

		friend tomlin;

		lobj (double mass, double damp, double t0, uint tsteps)
			: mass(mass), damp(damp), t0(t0)
		{
			poss.reserve(tsteps);
			vels.reserve(tsteps);
			accs.reserve(tsteps);
		}
	};

	// model parameters
	const double spring1;
	const double spring2;
	const double supvel;
	const double latcon;
	const double align;
	const double barr1;	
	const double barr2;	
	const double kappa1;
	const double kappa2;
	const double nu2;
	const double nu4;

	// dyanamic varaibles
	double kin = 0;
	double pot = 0;
	double fric = 0;
	double time = 0;
	double suppos = 0;
	double temp;

	// non physical parameters
	const double tstep;
	const uint tsteps;
	bool reverse = false;
	bool pause = false;

	// derived constants
	const double latcona = latcon;
	const double latconb = latcona/align;
	const double rlatcona = 1.0/latcon;
	const double rlatcona2pi = 2*pi*rlatcona;
	const double rlatconb2pi = 2*pi*latcon/align;

	// objects
	xobj x;
	lobj l;

	// containers
	vector <double> pots;
	vector <double> kins;
	vector <double> frics;
	vector <double> times;
	vector <double> poss;	
	vector <double> lessnoiset;
	vector <double> lessnoisef;
	vector <double> lessnoisesuppos;

	// functions
    void xvel(double);
	void xacc(double);
	void ldot(double);
	void testxacc();
	void testldot();
	pair<double,double> langevin();
	pair<double,double> testlangevin();

	// files
	const string tfile;
	const string xfile;
	const string lfile;
	const string tomfile;
	const string avgfile;
	const string pfile;

	// need for speed
	const double oneosix = 1.0/6.0;

	public:

	// functions
	void rk4l(double);
	void rk4x();
	void rk4();
	void eulmar();
	void calcpot();
	void calcfric();
	void inctime();

	void noisered(uint,uint);

	// in line functions
	void calckin(){x.calckin(); l.calckin(); kin = x.kin + l.kin;}
	void treverse(){reverse = !reverse;}
	void tpause(){pause = !pause;}

	// accessors
	void setposx(double c){x.pos = c;}
	void setvelx(double c){x.vel = c;}
	void setaccx(double c){x.acc = c;}
	void setposl(double c){l.pos = c;}
	void setvell(double c){l.vel = c;}
	void setaccl(double c){l.acc = c;}
	void settime(double c){time = c;}
	void setsuppos(double c){suppos = c;}

	double getposx(uint t){return x.poss[t];}
	double getvelx(){return x.vel;}
	double getaccx(){return x.acc;}
	double getposl(){return l.pos;}
	double getvell(){return l.vel;}
	double getaccl(){return l.acc;}
	uint getposxsize(){return x.poss.size();}

	double gettime(double t){return times[t];}	
	double getkin(){return kin;}
	double getpot(){return pot;}
	double getfric(double t){return frics[t];}
	double getsuppos(){return suppos;} 
	double getrnfric(double t) {return lessnoisef[t];}
	double getrntime(double t) {return lessnoiset[t];}
	double getrnsuppos(double t) {return lessnoisesuppos[t];}

	vector <double>* getposs() {return &x.poss;}
	vector <double>* getfrics() {return &frics;}
	vector <double>* gettimes() {return &times;}
	vector <double>* getrnfrics() {return &lessnoisef;}
	vector <double>* getrntimes() {return &lessnoiset;}
	vector <double>* getrnsupposs() {return &lessnoisesuppos;}

	// io stuff
	void pushvals();
	void pushvalslite();
	void writedata();
	void writedatalite();

	// debugging
	void printins();
	void printouts();
	void printavgs();
	void justkicks();

	// constructor
	tomlin(double spring1, double spring2, double supvel, double latcon, double align, double barr1, 
		   double barr2, double kappa1, double kappa2, double nu2, double nu4, 
		   double temp, double tstep, uint tsteps, double xmass, double lmass, 
		   double xdamp, double ldamp, double t0, string tfile, string xfile, string lfile, 
		   string tomfile, string avgfile, string pfile)
		 : spring1(spring1), spring2(spring2), supvel(supvel), latcon(latcon), align(align), barr1(barr1),
	       barr2(barr2), kappa1(kappa1), kappa2(kappa2), nu2(nu2), nu4(nu4), temp(temp),
		   tstep(tstep), tsteps(tsteps), x(xmass,xdamp,tsteps), 
		   l(lmass,ldamp,t0,tsteps), tfile(tfile), xfile(xfile), 
		   lfile(lfile), tomfile(tomfile), avgfile(avgfile), pfile(pfile)
		{
		   kins.reserve(tsteps);
		   pots.reserve(tsteps);
		   frics.reserve(tsteps);
		   times.reserve(tsteps);
		   poss.reserve(tsteps);
		}
};

// out of line function members

void tomlin::calcpot()
{
	pot = 0.5*spring1*pow(x.pos-suppos*time,2) + 
		  barr1*(1-cos(rlatcona2pi*(x.pos))); 
}
	
void tomlin::xacc(double kick)
{
    // force part
    double force = spring1 * x.pos;
    
    x.acc = force/x.mass;
    //x.acc = - (x.rmass*force);
	
    // acceleration part
    //x.acc = x.acc - x.vel * x.damp + kick;
}

void tomlin::xvel(double kick)
{
    // force part
    double force = 2*pi*barr1*rlatcona*sin(2*pi*x.pos*rlatcona); 
    x.vel = (-force - spring1*(x.pos - suppos))/(x.damp*x.mass) + kick/(x.damp*x.mass);
    
    //x.vel = -spring1*(x.pos*x.rmass)/x.damp + kick*x.rmass/x.damp;
   
    //x.vel = -(force);
}

void tomlin::ldot(double kick)
{
    l.vel =  -l.pos*l.rt0/l.damp + kick;

    //l.vel = 0;
}


void tomlin::calcfric()
{
	fric = spring1*(suppos-x.pos);
}

void tomlin::pushvals()
{
	x.poss.push_back(x.pos);
	x.vels.push_back(x.vel);
	x.accs.push_back(x.acc);
	//l.poss.push_back(l.pos);
	//l.vels.push_back(l.vel);
	//l.accs.push_back(l.acc);

	kins.push_back(kin);
	pots.push_back(pot);
	frics.push_back(fric);
	times.push_back(time);
	poss.push_back(suppos);
}

void tomlin::pushvalslite()
{
	x.poss.push_back(x.pos);
	x.vels.push_back(x.vel);
	times.push_back(time);
}

void tomlin::inctime()
{
	time += tstep;
	
	if (reverse)
	{
		suppos -= supvel*tstep;
	}
	else if (pause)
	{
		suppos = 1.0*suppos;
	}
	else
	{
		suppos += supvel*tstep;
	}
}

void tomlin::writedata()
{   
	ofstream xstream, lstream, tomstream, tstream, avgstream, pstream;

	xstream.open(xfile);
		xstream << "position, velocity, acceleration" << endl;
		for (uint k = 0; k < tsteps; k++)
			xstream << setprecision(16) << x.poss[k] << "," << x.vels[k] << "," << x.accs[k] << endl;
	xstream.close();
	
	//lstream.open(lfile);
	//	lstream << "position, velocity, acceleration" << endl;
	//	for (uint k = 0; k < tsteps; k++)
	//		lstream << setprecision(16) << l.poss[k] << "," << l.vels[k] << "," << l.accs[k] << endl;
	//lstream.close();
	
	tomstream.open(tomfile);
		tomstream << "kinetic, potential, friction" << endl;
		for (uint k = 0; k < tsteps; k++)
			tomstream << setprecision(16) << setprecision(16) << kins[k] << "," << pots[k] << "," << frics[k] << endl;
	tomstream.close();

	tstream.open(tfile);
		tstream << "time, displacement" << endl;
		for (uint k = 0; k < tsteps; k++)
			tstream << setprecision(16) << times[k] << "," << poss[k] << endl;
	tstream.close();

	avgstream.open(avgfile);
		avgstream << "time, avged friction " << endl;
		uint size = lessnoisef.size();
		for (uint k = 0; k < size; k++)
			avgstream << setprecision(16) << lessnoiset[k] << "," << lessnoisef[k] << "," << lessnoisesuppos[k] << endl;
	avgstream.close();
	
	pstream.open(pfile);
		pstream << setprecision(16) 
                << "in order: spring1, spring2, t0, supvel, align, latcona, latconb, barr1, barr2, kappa1, kappa2, nu2, nu4, tstep, tsteps, xmass, lmass, xdamp, ldamp" << endl
                << spring1 << endl << spring2<< endl << l.t0 << endl << suppos 
                << endl << align << endl << latcona << endl << latconb  
                << endl << barr1  << endl << barr2  << endl << kappa1 
                << endl << kappa2 << endl << nu2 << endl << nu4  
                << endl << tstep  << endl << tsteps << endl << x.mass  
                << endl << l.mass  << endl << x.damp  << endl << l.damp  << endl;
	pstream.close();
} 

void tomlin::writedatalite()
{   
	ofstream xstream, tomstream, tstream;

	xstream.open(xfile);
		xstream << "position" << endl;
		for (uint k = 0; k < tsteps; k++)
			xstream << setprecision(16) << x.poss[k] << endl;
	xstream.close();
} 

void tomlin::printins()
{
	cout << "the following parameters was used: " << endl
		 << "spring " << spring1 << endl << "supvel " << suppos << endl 
		 << "latcon " << latcon << endl << "align "  << align  << endl 
		 << "barr1 "  << barr1  << endl << "barr2 "  << barr2  << endl 
		 << "kappa1 " << kappa1 << endl << "kappa2 " << kappa2 << endl 
		 << "nu2 "	  << nu2	<< endl << "nu4 "	 << nu4	   << endl 
		 << "tstep "  << tstep  << endl << "tsteps " << tsteps << endl 
		 << "xmass "  << x.mass  << endl << "lmass "  << l.mass  << endl 
		 << "xdamp "  << x.damp  << endl << "ldamp "  << l.damp  << endl
		 << "tfile "  << tfile  << endl << "xfile "  << xfile  << endl
		 << "lfile "  << lfile  << endl << "tomfile " << tomfile << endl;
}

void tomlin::printouts()
{
	cout << "and the exit values were: " << endl
		 << "time " << time << endl << "suppos " << suppos << endl << "kin " 
		 << kin << endl		 << "pot " << pot << endl << "tot en " <<  kin + pot 
		 << endl << "xpos " << x.pos << endl << "xvel " << x.vel << endl 
		 << "xacc " << x.acc << endl << "lpos " << l.pos << endl << "lvel " 
		 << l.vel << endl << "ldot " << l.acc << endl;

}

void tomlin::printavgs()
{
	double avgkin = accumulate(kins.begin(),kins.end(),0.0) / kins.size();

	cout << "avg kin " << avgkin << endl;
	
	double avgfric = accumulate(frics.begin(),frics.end(),0.0) / frics.size();

	cout << "avg fric " << avgfric << endl;
}

void tomlin::rk4l(double kick)
{
    double k1,k2,k3,k4;
    double oldpos;

    double t0 = time;
    double t1 = time + 0.5*tstep;
    double t2 = time + tstep;
    

    oldpos = l.pos;
    ldot(kick);

    k1 = l.vel;

    l.pos = oldpos*0.5*tstep*l.vel;
    ldot(kick);

    k2 = l.vel;
    
    l.pos = oldpos*0.5*tstep*l.vel;
    ldot(kick);

    k3 = l.vel;

	l.pos = oldpos+tstep*k3;
	ldot(kick);
	
	k4 = l.vel;

	l.pos = oldpos + tstep*oneosix*(k1 + 2*k2 + 2*k3 + k4);
}

void tomlin::rk4x()	// two degree of freedom RK4 algorithm
{
	// preliminaries
	
    auto kicks = langevin();	// random kicks from langevin dynamics, see below
	double xkick = kicks.first;
	double lkick = kicks.second;

    rk4l(lkick);

	// container for old positions, velocities, and accelerations
	
	struct oldvals {double pos; double vel;};
	struct k {double vel; double acc;};
	oldvals oldx;
	double oldt;

	// calculate accelerations in current configuration
	xacc(xkick); 

	// k1

	// save old values, at this point this is the ones from the beginning of 
	// the timestep

	oldx.pos = x.pos;
   	oldx.vel = x.vel;
	//oldl.acc = l.acc;
	
	// calculate the RK4 k1 constants, as: 
	// * velocity is the start of timestep velocity since we're in k1, so no time 
	//   has elapsed
	// * acceleration is the acceleration as calculated from the start of timestep 
	//   configuration langevin thermostated since no time has elapsed

	k k1x = {x.vel, x.acc};
	
	// k2
	
	// the RK4 k2 procedure is very similar to that of k1. but for k2 time is advanced
	// by half a timestep. this means we euler forward the dynamic variables half a 
	// timestep (advancing from the start of timestep configuration and using k1 as 
	// the present configuration) before we compute the new acceleration from the 
	// obtained half timestep configuration, from with we then calculate the k2 
	// constants identically to how we did for k1

	x.pos = oldx.pos+0.5*tstep*k1x.vel;
	x.vel = oldx.vel+0.5*tstep*k1x.acc;
	
	xacc(xkick);

	k k2x = {x.vel,x.acc};
	
	// k3
	
	// the RK4 k3 prodecure is identical to that of k2, just the k2 configuration is 
	// used as the present configuration in the euler step

	x.pos = oldx.pos+0.5*tstep*k2x.vel;
	x.vel = oldx.vel+0.5*tstep*k2x.acc;
	//l.vel = oldl.vel+0.5*tstep*k2l.acc;
	
	xacc(xkick);
	
	k k3x = {x.vel,x.acc};
	
	// k4
	
	// the RK4 k4 prodecure is identical to that of k3, but the present configuration 
	// used here is one full timestep from the start of timestep configuration

	time = oldt + tstep;
	x.pos = oldx.pos+tstep*k3x.vel;
	x.vel = oldx.vel+tstep*k3x.acc;
	
	xacc(xkick);
	
	k k4x = {x.vel,x.acc - x.vel * x.damp + xkick};

	// obtained approximation -- just sum everything up according to RK4 formula

	x.pos = oldx.pos + tstep*oneosix*(k1x.vel + 2*k2x.vel + 2*k3x.vel + k4x.vel);
	x.vel = oldx.vel + tstep*oneosix*(k1x.acc + 2*k2x.acc + 2*k3x.acc + k4x.acc);
}

void tomlin::rk4() // two degree of freedom RK4 algorithm
{    
    // preliminaries
	
    double oldxpos;
    //double oldlpos;
    double oldt;
    
    auto kicks = langevin();	// random kicks from langevin dynamics, see below
	double xkick = kicks.first;
	//double lkick = kicks.second;


    xvel(xkick);
    //ldot(lkick);

	// k1

	// save old values, at this point this is the ones from the beginning of 
	// the timestep
    
    oldt = time;
    oldxpos = x.pos;
    //oldlpos = l.pos;
    	
	// calculate the RK4 k1 constants, as: 
	// * velocity is the start of timestep velocity since we're in k1, so no time 
	//   has elapsed
	// * acceleration is the acceleration as calculated from the start of timestep 
	//   configuration langevin thermostated since no time has elapsed

	double k1x = x.vel;
	//double k1l = l.vel;
	
	// k2
	
	// the RK4 k2 procedure is very similar to that of k1. but for k2 time is advanced
	// by half a timestep. this means we euler forward the dynamic variables half a 
	// timestep (advancing from the start of timestep configuration and using k1 as 
	// the present configuration) before we compute the new acceleration from the 
	// obtained half timestep configuration, from with we then calculate the k2 
	// constants identically to how we did for k1

	x.pos = oldxpos+0.5*tstep*k1x;
	//l.pos = oldlpos+0.5*tstep*k1l;
	
    xvel(xkick);
    //ldot(lkick);

	double k2x = x.vel;
	//double k2l = l.vel;
	
	// k3
	
	// the RK4 k3 prodecure is identical to that of k2, just the k2 configuration is 
	// used as the present configuration in the euler step
	
    x.pos = oldxpos+0.5*tstep*k2x;
	//l.pos = oldlpos+0.5*tstep*k2l;
	
    xvel(xkick);
    //ldot(lkick);

	double k3x = x.vel;
	//double k3l = l.vel;

	// k4
	
	// the RK4 k4 prodecure is identical to that of k3, but the present configuration 
	// used here is one full timestep from the start of timestep configuration

	time = oldt + tstep;
	x.pos = oldxpos+tstep*k3x;
	//l.pos = oldlpos+tstep*k3l;
	
    xvel(xkick);
    //ldot(lkick);
	
	double k4x = x.vel;
	//double k4l = l.vel;

	// obtained approximation -- just sum everything up according to RK4 formula

    
	//cout << oldxpos << " " << tstep*oneosix*(k1x + 2*k2x + 2*k3x + k4x) << endl;
	x.pos = oldxpos + tstep*oneosix*(k1x + 2*k2x + 2*k3x + k4x);
	//l.pos = oldlpos + tstep*oneosix*(k1l + 2*k2l + 2*k3l + k4l);
}

pair <double, double> tomlin::langevin()
{
	random_device rd;
	mt19937 gen(rd());

	double mean = 0;
	//
    //double xstd = sqrt(2*kB*temp*x.damp/tstep);
	//double lstd = 0.9*xstd;

	//normal_distribution<double> xkick(mean,xstd);
	//normal_distribution<double> lkick(mean,lstd);

	//pair <double, double> out = {xkick(gen),lkick(gen)};
	
	double xstd = 1.0;
	double lstd = 1.0;
    
    normal_distribution<double> xkick(mean,xstd);
	//normal_distribution<double> lkick(mean,lstd);
	
    double c1 = sqrt(2*kB*temp*x.mass*x.damp/tstep);
    //double c2 = sqrt(2*kb*temp*l.damp/tstep);
	pair <double, double> out = {c1*xkick(gen),0};

	return out;
}

void findslips (uint mode, uint adj, uint stride, vector <double>* tvals, vector <double>* fvals, vector <uint>* slips)
{

    vector <uint> slipsish;

    // method I: identify slipping regions by slope, then fins slipping points
    // by interval halving
    if (mode == 1)
    {
        uint curr = 0;
        uint next = stride;
        uint end = fvals->size();
        bool climbing = true;
        while (next < end)
        {
            //uint next = start + round(adj+exp(logstep*loops));

            //cout << curr << " " << next << " " << end << " " << tvals->size() << " " << fvals->size() << endl;
            double x1 = tvals->at(curr);
            double x2 = tvals->at(next);

            double y1 = fvals->at(curr);
            double y2 = fvals->at(next);

            double slope = (y2 - y1) / (x2 - x1);

            if (slope > 0.0)
                climbing = true;

            if (slope < 0.0 && climbing)
            {
                climbing = false;
                slipsish.push_back(curr);
            }
            curr = next;
            next = curr + stride;

            //cout << x1 << " " << x2 << " " << y1 << " " << y2 << " " << x2-x1 << " " << y2-y1 << " " << slope << endl;
        }

        //cout << "slipsish has size: " << slipsish.size() << endl;

        //for (auto &el : slipsish)
        //	cout << el << endl;
        //cout << endl;
        
        slips->reserve(round(end/( (double) stride )));
        
        for (uint k = 0; k < slipsish.size() - 1; k++)
        {
            curr = slipsish[k];
            next = slipsish[k]+stride;
            //cout << curr << " " << next << " " << end << endl << " " << tvals->size() << " " << fvals->size() << endl;

            if (next >= end)     // slightly ugly termination
                break;

            double zero = 1e-11;
            double x1 = 0;
            x1 = tvals->at(curr);

            bool left = true;						// slip is to the left
            uint distance = abs((int)next - (int)curr);		// this is stridem but whatever...
            uint half = 0;	
            half = round( half + distance/2.0 );
            double xhalf = tvals->at(half);

            //cout << "curr is: " << curr << " next is: " << next << " distance is: " << distance << endl;

            //cout << "distance greater than adj? " << ((distance > adj) ? 1 : 0) << " " << endl;	
            while (distance > adj)
            {
                if (left)
                {
                    if (x1 - xhalf > zero)
                    {
                        left = true;
                    }
                    else
                    {
                        left = false;
                        swap(curr,half);
                    }
                }
                else
                {
                    if (x1 - xhalf > zero)
                    {
                        left = false;
                    }
                    else
                    {
                        left = true;
                        swap(curr,half);
                    }
                }
                
                distance = abs((int)half-(int)curr);
                
                x1 += tvals->at(curr);
                half = round( half + distance/2.0 );
                xhalf += tvals->at(half);
            }

            if(left)
            {
                slips->push_back(curr);
            }
            else
            {
                slips->push_back(next);
            }
        }
    }
    // method II find slipping rebion by difference tolerance, then find 
    else if (mode == 2)
    {
        uint curr = 0;
        uint next = stride;
        uint end = fvals->size();
        bool climbing;
        while (next < end)
        {
            //uint next = start + round(adj+exp(logstep*loops));

            double x1 = tvals->at(curr);
            double x2 = tvals->at(next);

            double y1 = fvals->at(curr);
            double y2 = fvals->at(next);

            double slope = (y2 - y1) / (x2 - x1);

            if (slope > 0.0)
                climbing = true;

            if (slope < 0.0 && climbing)
            {
                climbing = false;
                slipsish.push_back(next);
            }
            curr = next;
            next = curr + stride;
        }

        //cout << "slipsish has size: " << slipsish.size() << endl;

        //for (auto &el : slipsish)
        //	cout << el << endl;
        //cout << endl;
        
        slips->reserve(round(end/( (double) stride )));
        
        for (uint k = 0; k < slipsish.size() - 1; k++)
        {

            uint curr = slipsish[k];
            uint prev = curr - stride;
            uint next = curr + stride;
            uint distance = stride;

            //bool leftofmin;

            double ycurr = fvals->at(curr);
            double yprev = fvals->at(prev);
            double ynext = fvals->at(next);

            int p;

            if (ycurr - ynext > 0)
            {
                //leftofmin = true;
                p = 1;
            }
            else 
            {
                //leftofmin = false;
                p = -1;
            }

            
            while (distance > adj)
            {
                distance = round(distance / 2.0);
                
                next = round(curr + p*distance);

                ycurr = fvals->at(curr);
                ynext = fvals->at(next);
                
                if (ycurr - ynext > 0)
                {
                    p = 1;
                }
                else
                {
                    p = -1;
                }

                curr = next;
            }

            slips->push_back(curr);
        }
    }
}
