
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
	
    class lambdaobj
 	{
		// dynamic variables
		double pos = 0;
		double vel = 0;
		double acc = 0;

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

		// functions

		friend tomlin;

		lambdaobj (double mass, double damp, double t0, uint tsteps) 
            : mass(mass), damp(damp), t0(t0) 
		{
			poss.reserve(tsteps);
			vels.reserve(tsteps);
 		}
	};

	// model parameters
	const double spring;
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
    lambdaobj l;

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
	void xacc();
	void hacc();
	void tacc();
    void ldot();
	pair <double,double> langevin();

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
	void rk4();
	void eulmar();
	//void calcpot();
	void calcfric();
	void inctime();

	void noisered(uint,uint);

	// in line functions
	void calckin(){x.calckin(); kin = x.kin;}
	void treverse(){reverse = !reverse;}
	void tpause(){pause = !pause;}

	// accessors
	void setposx(double c){x.pos = c;}
	void setvelx(double c){x.vel = c;}
	void setaccx(double c){x.acc = c;}
	void settime(double c){time = c;}
	void setsuppos(double c){suppos = c;}

	double getposx(uint t){return x.poss[t];}
	double getvelx(){return x.vel;}
	double getaccx(){return x.acc;}
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
	void writedata();

	// debugging
	void printins();
	void printouts();
	void printavgs();
	void justkicks();

	// constructor
	tomlin(double spring, double supvel, double latcon, double align, double barr1, 
		   double barr2, double kappa1, double kappa2, double nu2, double nu4, 
		   double temp, double tstep, uint tsteps, double xmass, double xdamp, 
           double lmass, double ldamp, double t0, string tfile, string xfile, 
           string lfile, string tomfile, string avgfile, string pfile)
		 : spring(spring), supvel(supvel), latcon(latcon), align(align), barr1(barr1),
	       barr2(barr2), kappa1(kappa1), kappa2(kappa2), nu2(nu2), nu4(nu4), temp(temp),
		   tstep(tstep), tsteps(tsteps), x(xmass,xdamp,tsteps), l(lmass,ldamp,t0,tsteps),
		   tfile(tfile), xfile(xfile), lfile(lfile), tomfile(tomfile), avgfile(avgfile), pfile(pfile)
		{
		   kins.reserve(tsteps);
		   pots.reserve(tsteps);
		   frics.reserve(tsteps);
		   times.reserve(tsteps);
		   poss.reserve(tsteps);
		}
};

// out of line function members

void tomlin::tacc()
{
	double force = spring*(x.pos-suppos) + rlatcona2pi*barr1*sin(rlatcona2pi*(x.pos));
	
    x.acc = - (x.rmass*force);
}

void tomlin::hacc()
{
	x.acc = - x.rmass*spring*(x.pos);
}

void tomlin::ldot()
{
    l.vel = l.rt0 * l.pos; 
}

void tomlin::calcfric()
{
	fric = spring*(suppos-x.pos);
}

void tomlin::pushvals()
{
	x.poss.push_back(x.pos);
	x.vels.push_back(x.vel);
	x.accs.push_back(x.acc);
	l.poss.push_back(x.pos);
	l.vels.push_back(x.vel);

	kins.push_back(kin);
	pots.push_back(pot);
	frics.push_back(fric);
	times.push_back(time);
	poss.push_back(suppos);
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
	
	lstream.open(lfile);
		lstream << "position, velocity, acceleration" << endl;
		for (uint k = 0; k < tsteps; k++)
			lstream << setprecision(16) << l.poss[k] << "," << l.vels[k] << endl;
    lstream.close();
	
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
                << "spring " << spring << endl << "supvel " << suppos << endl 
				<< "align " << align << endl   
				<< "latcona " << latcona << endl << "latconb "  << latconb  << endl 
				<< "barr1 "  << barr1  << endl << "barr2 "  << barr2  << endl 
				<< "kappa1 " << kappa1 << endl << "kappa2 " << kappa2 << endl 
				<< "nu2 "	  << nu2	<< endl << "nu4 "	 << nu4	   << endl 
				<< "tstep "  << tstep  << endl << "tsteps " << tsteps << endl
				<< "xmass "  << x.mass  << endl << "lmass "  << l.mass  << endl 
				<< "xdamp "  << x.damp  << endl << "ldamp "  << l.damp  << endl;
	pstream.close();
} 

void tomlin::printins()
{
	cout << "the following parameters was used: " << endl
		 << "spring " << spring << endl << "supvel " << suppos << endl 
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
		 << l.vel << endl;

}

void tomlin::printavgs()
{
	double avgkin = accumulate(kins.begin(),kins.end(),0.0) / kins.size();

	cout << "avg kin " << avgkin << endl;
	
	double avgfric = accumulate(frics.begin(),frics.end(),0.0) / frics.size();

	cout << "avg fric " << avgfric << endl;
}

void tomlin::rk4()	// two degree of freedom RK4 algorithm
{
	// preliminaries

	// container for old positions, velocities, and accelerations
	
	struct oldvals {double pos; double vel; double acc;};
	struct k {double vel; double acc;};
	oldvals oldx;
	oldvals oldl;
	double oldt;

	// calculate accelerations in current configuration
	hacc();
    ldot();

    auto kicks = langevin();
	double xkick = kicks.first;	// random kicks from langevin dynamics, see below
	double lkick = kicks.second;	// random kicks from langevin dynamics, see below

	// k1

	// save old values, at this point this is the ones from the beginning of 
	// the timestep

	oldt = time;
	oldx.pos = x.pos;
   	oldx.vel = x.vel;
	oldx.acc = x.acc;
	oldl.pos = x.pos;
   	oldl.vel = x.vel;
	
	// calculate the RK4 k1 constants, as: 
	// * velocity is the start of timestep velocity since we're in k1, so no time 
	//   has elapsed
	// * acceleration is the acceleration as calculated from the start of timestep 
	//   configuration langevin thermostated since no time has elapsed

	k k1x = {x.vel, x.acc - x.vel*x.damp + xkick};
	k k1l = {l.vel, -l.vel*l.damp + lkick};
    
    //cout << k1x.vel << " " << k1x.acc << " " << k1l.vel << " " << k1x.acc << endl;

	// k2
	
	// the RK4 k2 procedure is very similar to that of k1. but for k2 time is advanced
	// by half a timestep. this means we euler forward the dynamic variables half a 
	// timestep (advancing from the start of timestep configuration and using k1 as 
	// the present configuration) before we compute the new acceleration from the 
	// obtained half timestep configuration, from with we then calculate the k2 
	// constants identically to how we did for k1

	time = oldt + 0.5*tstep;
	x.pos = oldx.pos+0.5*tstep*k1x.vel;
	x.vel = oldx.vel+0.5*tstep*k1x.acc;
	l.pos = oldl.pos+0.5*tstep*k1l.vel;
	//l.vel = oldl.vel+0.5*tstep*k1l.acc;
	
	hacc();
    ldot();

	k k2x = {x.vel, x.acc - x.vel*x.damp + xkick};
	k k2l = {l.vel, -l.vel*l.damp + lkick};
	
	// k3
	
	// the RK4 k3 prodecure is identical to that of k2, just the k2 configuration is 
	// used as the present configuration in the euler step

	time = oldt + 0.5*tstep;
	x.pos = oldx.pos+0.5*tstep*k2x.vel;
	x.vel = oldx.vel+0.5*tstep*k2x.acc;
	l.pos = oldl.pos+0.5*tstep*k2l.vel;
	//l.vel = oldl.vel+0.5*tstep*k2l.acc;
	
	hacc();
    ldot();
	
	k k3x = {x.vel, x.acc - x.vel*x.damp + xkick};
	k k3l = {l.vel, -l.vel*l.damp + lkick};
	
	// k4
	
	// the RK4 k4 prodecure is identical to that of k3, but the present configuration 
	// used here is one full timestep from the start of timestep configuration

	time = oldt + tstep;
	x.pos = oldx.pos+tstep*k3x.vel;
	x.vel = oldx.vel+tstep*k3x.acc;
	l.pos = oldl.pos+0.5*tstep*k3l.vel;
	//l.vel = oldl.vel+0.5*tstep*k3l.acc;
	
	hacc();
    ldot();
	
	k k4x = {x.vel, x.acc - x.vel*x.damp + xkick};
	k k4l = {l.vel, -l.vel*l.damp + lkick};

	// obtained approximation -- just sum everything up according to RK4 formula

	x.pos = oldx.pos + tstep*oneosix*(k1x.vel + 2*k2x.vel + 2*k3x.vel + k4x.vel);
	x.vel = oldx.vel + tstep*oneosix*(k1x.acc + 2*k2x.acc + 2*k3x.acc + k4x.acc);
	l.pos = oldl.pos + tstep*oneosix*(k1l.vel + 2*k2l.vel + 2*k3l.vel + k4l.vel);
	//l.vel = oldl.vel + tstep*oneosix*(k1l.acc + 2*k2l.acc + 2*k3l.acc + k4l.acc);

    //x.acc = xkick;
	// reset the timestep (time is ticked at the end of the main loop)
	
	time = oldt;
}

pair <double,double> tomlin::langevin()
{
	random_device rd;
	mt19937 gen(rd());

	double mean = 0;
	double xstd = sqrt(2*x.mass*kB*temp*x.damp/tstep);
	double lstd = sqrt(2*l.mass*kB*temp*l.damp/tstep);

	normal_distribution<double> xkick(mean,xstd);
	normal_distribution<double> lkick(mean,lstd);

	pair <double,double> out = {xkick(gen)*x.rmass, lkick(gen)*l.rmass};

	return out;
}


