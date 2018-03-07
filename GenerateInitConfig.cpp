//
//  GenerateInitConfig.cpp
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 12/6/12.
//  Copyright (c) 2012 Ryohei Seto and Romain Mari. All rights reserved.
//

#include <stdlib.h> // necessary for Linux
#include "GenerateInitConfig.h"
#include "Simulation.h"
#ifndef USE_DSFMT
#define RANDOM ( rand_gen.rand() ) // RNG uniform [0,1]
#endif
#ifdef USE_DSFMT
#define RANDOM ( dsfmt_genrand_close_open(&rand_gen) ) // RNG uniform [0,1]
#endif
using namespace std;

int GenerateInitConfig::generate(int rand_seed_, bool magnetic_config_)
{
	setParameters();
	rand_seed = rand_seed_;
	magnetic_config = magnetic_config_;
	sys.set_np(np);
	sys.friction = false;
	sys.repulsiveforce = false;
	sys.p.interaction_range = 2.5;
	sys.p.lub_max_gap = 0.5;
	sys.allocateRessources();
	sys.setBoxSize(lx, ly, lz);
	sys.setSystemVolume(2*a2);
	sys.in_predictor = false;
	sys.p.integration_method = 0;
	putRandom();
	double inflate_ratio = 1.03;
	for (int i=0; i<np;i++) {
		sys.radius[i] *= inflate_ratio;
	}
	sys.setInteractions_GenerateInitConfig();
	grad = new vec3d [np];
	prev_grad = new vec3d [np];
	step_size = 20;
	gradientDescent();
	step_size /= 4.;
	gradientDescent();
	step_size /= 4.;
	gradientDescent();
	// step_size /= 4.;
	// gradientDescent();
	// step_size /= 4.;
	// gradientDescent();
	// step_size /= 4.;
	// gradientDescent();
	// step_size /= 4.;
	// gradientDescent();
	// step_size /= 4.;
	// gradientDescent();
	int count = 0;
	double energy = 0;
	ofstream fout;
	fout.open("energy_decay.dat");
	double energy_previous = 0;
	double diff_energy = 99999;
	do {
		energy = zeroTMonteCarloSweep();
		cerr << energy << endl;
		diff_energy = energy-energy_previous;
		energy_previous = energy;
		if (count % 100 == 0) {
			fout << count << ' ' << energy << ' ' <<  diff_energy << endl;
		}
		count ++;
		cerr << "diff = " << abs(diff_energy) << endl;
	} while (abs(diff_energy) > 1e-10);
	//deflate
	for (int i=0; i < np; i++) {
		//if (i < np1) {
			sys.radius[i] = a_12[i];
		//} else {
		//	sys.radius[i] = a2;
		//}
	}
	position.resize(np);
	radius.resize(np);
	for (int i=0; i<sys.get_np(); i++) {
		position[i] = sys.position[i];
		radius[i] = sys.radius[i];
	}
	outputPositionData();
	delete [] grad;
	delete [] prev_grad;
	return 0;
}

void GenerateInitConfig::outputPositionData()
{
	ofstream fout;
	ofstream fout_yap;
	ostringstream ss_posdatafilename;
	if (sys.twodimension) {
		ss_posdatafilename << "D2";
	} else {
		ss_posdatafilename << "D3";
	}
	ss_posdatafilename << "N" << np;
	ss_posdatafilename << "VF" << volume_fraction;
	if (disperse_type == 'm') {
		ss_posdatafilename << "Mono";
	} else if (disperse_type == 'b') {
		ss_posdatafilename << "Bidi" << a2 << "_" << vf_ratio;
    } else if (disperse_type == 'p') {
        ss_posdatafilename << "Poly" << a2 << "_" << vf_ratio;
    } else {
		cerr << "disperse_type is wrong." << endl;
		exit(1);
	}
	if (sys.twodimension) {
		if (lx_lz == 1) {
			ss_posdatafilename << "Square"; // square
		} else {
			ss_posdatafilename << "L" << (int)(10*lx_lz) << "_" << 10;
		}
	} else {
		if (lx_lz == 1 && ly_lz == 1) {
			ss_posdatafilename << "Cubic"; //
		} else {
			ss_posdatafilename << "L" << (int)(10*lx_lz) << "_" << (int)(10*ly_lz) << "_" << 10;
		}
	}
	vector<double> magnetic_susceptibility;
	for (int i=0; i< np; i++) {
		if (i < np1) {
			magnetic_susceptibility.push_back(1);
		} else {
			magnetic_susceptibility.push_back(-1);
		}
	}
	
	ss_posdatafilename << "_" << rand_seed << ".dat";
	cerr << ss_posdatafilename.str() << endl;
	fout.open(ss_posdatafilename.str().c_str());
	ss_posdatafilename << ".yap";
	fout_yap.open(ss_posdatafilename.str().c_str());
	fout << "# np1 np2 vf lx ly lz vf1 vf2 disp" << endl;
	fout << "# " << np1 << ' ' << np2 << ' ' << volume_fraction << ' ';
	fout << lx << ' ' << ly << ' ' << lz << ' ';
	fout << volume_fraction1 << ' ' << volume_fraction2 << ' ' << 0 << endl;
	for (int i = 0; i < np; i++) {
		fout << position[i].x << ' ';
		fout << position[i].y << ' ';
		fout << position[i].z << ' ';
		fout << radius[i] << ' ';
		if (magnetic_config) {
			fout << 0 << ' ';
			fout << 0 << ' ';
			fout << 0 << ' ';
			fout << magnetic_susceptibility[i];
		}
		fout << endl;
		fout_yap << "r ";
		fout_yap << radius[i] << endl;
		fout_yap << "c ";
		fout_yap << position[i].x << ' ';
		fout_yap << position[i].y << ' ';
		fout_yap << position[i].z << endl;
	}
	fout_yap << "@ 4 \n";
	fout_yap << "y 4 \n";
	for (int k=0; k<sys.nb_interaction; k++) {
		unsigned short i,j;
		sys.interaction[k].get_par_num(i, j);
		vec3d d_pos = position[i]-position[j];
		if (d_pos.norm() < 10){
			fout_yap << "l ";
			fout_yap << position[i].x << ' ';
			fout_yap << position[i].y << ' ';
			fout_yap << position[i].z << ' ';
			fout_yap << position[j].x << ' ';
			fout_yap << position[j].y << ' ';
			fout_yap << position[j].z << endl;
		}
	}
	fout.close();
}

double GenerateInitConfig::computeGradient()
{
	for(int i=0; i<np; i++) {
		grad[i].reset();
	}
	unsigned short i,j;
	double r, rcont;
	double amp, amp2;
	double energy = 0;
	for (int k=0; k<sys.nb_interaction; k++) {
		if (sys.interaction[k].is_contact()) {
			sys.interaction[k].get_par_num(i, j);
			r = sys.interaction[k].get_r();
			rcont = sys.interaction[k].get_ro();
			const vec3d &nr_vec = sys.interaction[k].nvec;
			amp = (1/rcont-1/r); // negative
			amp2 = 4*amp/rcont;
			grad[i] -= r*nr_vec*amp2;
			grad[j] += r*nr_vec*amp2;
			energy += 2*r*amp*amp;
		}
	}
	return energy;
}

void GenerateInitConfig::moveAlongGradient(vec3d *g, int dir)
{
	double grad_norm;
	double gradient_power = 0.5;
	vec3d step;
	grad_norm = 0;
	for (int i=0; i<np;i++) {
		grad_norm += g[i].sq_norm();
	}
	if (grad_norm != 0) {
		double rescale = pow(grad_norm, gradient_power);
		for (int i=0; i<np;i++) {
			step = -dir*g[i]*step_size/rescale;
			sys.displacement(i, step);
		}
		sys.checkNewInteraction();
		sys.updateInteractions();
	}
}

void GenerateInitConfig::storeGradient()
{
	for (int i=0; i<np; i++) {
		prev_grad[i] = grad[i];
	}
}

double GenerateInitConfig::gradientDescent()
{
	double old_running_energy;
	double running_energy;
	double relative_en;
	long long int steps = 0;
	cerr << endl << " Gradient Descent..." << endl;
	storeGradient();
	running_energy = computeGradient();
	cerr << "  Starting Energy " << running_energy/np << endl;
	do {
		old_running_energy = running_energy;
		moveAlongGradient(grad, 1);
		storeGradient();
		running_energy = computeGradient();
		relative_en = (old_running_energy-running_energy)/(old_running_energy+running_energy);
		if (steps%100 == 0) {
			cerr << "    Steps = " << steps << " :::   Energy : " << running_energy/np << endl;
		}
		steps++;
	} while(relative_en > 1e-6);
	if (relative_en < 0) {
		cerr << "    Steps = " << steps;
		cerr << " :::   Last Step Upwards. Old Energy : " << old_running_energy/np;
		cerr << " New Energy : " << running_energy/np;
		cerr << " Relative : " << relative_en << endl;
		cerr << "      Reverting last step..." << endl;
		moveAlongGradient(prev_grad, -1);
		return old_running_energy;
	}
	if (relative_en > 0
		&& relative_en < 1e-6) {
		cerr << "    Steps = " << steps;
		cerr << " :::   Stuck: too slow (Relative energy difference : " << relative_en  << endl;
		cerr << "  Ending Energy "<< running_energy << endl<< endl;
		return running_energy;
	}
	return running_energy;
}

void GenerateInitConfig::putRandom()
{
	sys.allocatePositionRadius();
#ifndef USE_DSFMT
	rand_gen.seed(rand_seed);
#endif
#ifdef USE_DSFMT
	dsfmt_init_gen_rand(&rand_gen, rand_seed) ; // hash of time and clock trick from MersenneTwister code v1.0 by Richard J. Wagner
#endif
	for (int i=0; i < np; i++) {
		sys.position[i].x = lx*RANDOM;
		sys.position[i].z = lz*RANDOM;
		if (sys.twodimension) {
			sys.position[i].y = ly_half;
		} else {
			sys.position[i].y = ly*RANDOM;
		}
		if (i < np1) {
			sys.radius[i] = a1;
		} else {
			sys.radius[i] = a2;
		}
	}
}

void GenerateInitConfig::updateInteractions(int i)
{
	vector <Interaction*> inter_list;
	for (auto && inter : sys.interaction_list[i]){
		inter_list.push_back(inter);
	}

	for (auto &il : inter_list) {
		if (il->is_active()) {
			bool desactivated = false;
			il->updateState(desactivated);
			if (desactivated) {
				sys.deactivated_interaction.push(il->get_label());
			}
		}
	}
}

int GenerateInitConfig::overlapNumber(int i)
{
	int overlaps = 0;
	for (auto&& inter : sys.interaction_list[i]){
		if (inter->is_overlap()) {
			overlaps++;
		}
	}
	return overlaps;
}

double GenerateInitConfig::particleEnergy(int i)
{
	double energy = 0;
	for (auto&& inter : sys.interaction_list[i]){
		if (inter->is_overlap()) {
			double r = inter->get_r();
			double rcont = inter->get_ro();
			double amp = inter->get_a_reduced()*(1/rcont-1/r); // negative
			energy += r*amp*amp;
		}
	}
	return energy;
}

double GenerateInitConfig::zeroTMonteCarloSweep()
{
	int steps = 0;
	int init_overlaps = 0;
	double init_energy = 0;
	for(int i=0; i<np; i++) {
	 	init_overlaps += overlapNumber(i);
	}
	for(int i=0; i<np; i++) {
	 	init_energy += particleEnergy(i);
	}
	double dx = 0.04;
	while (steps < np) {
		int moved_part = (int)(RANDOM*np);
		//		int overlap_pre_move = overlapNumber(moved_part);
		double energy_pre_move = particleEnergy(moved_part);
		vec3d trial_move;
		if (sys.twodimension) {
			trial_move = randUniformCircle(dx);
		} else {
			trial_move = randUniformSphere(dx);
		}
		trial_move *= RANDOM;
		sys.displacement(moved_part, trial_move);
		updateInteractions(moved_part);
		//		int overlap_post_move = overlapNumber(moved_part);
		double energy_post_move = particleEnergy(moved_part);
		//		if( overlap_pre_move <= overlap_post_move ){
		if (energy_pre_move < energy_post_move) {
			sys.displacement(moved_part, -trial_move);
			updateInteractions(moved_part);
		}
		//sys.checkNewInteraction();
		steps ++;
	}
	sys.boxset.update();
	sys.checkNewInteraction();
	sys.updateInteractions();
	int final_overlaps = 0;
	double final_energy = 0;
	for(int i=0; i<np; i++) {
		final_overlaps += overlapNumber(i);
	}
	for(int i=0; i<np; i++) {
	 	final_energy += particleEnergy(i);
	}
	cerr << " MC sweep : init energy " << init_energy/np << " final energy " << final_energy/np;
	cerr << " init overlaps " << init_overlaps << " final overlaps " << final_overlaps << endl;
	return final_energy/np;
}

vec3d GenerateInitConfig::randUniformSphere(double r)
{
	double z = 2*RANDOM-1;
	double phi = 2*M_PI*RANDOM;
	double sin_theta = sqrt(1-z*z);
	return vec3d(r*sin_theta*cos(phi), r*sin_theta*sin(phi), r*z);
}

vec3d GenerateInitConfig::randUniformCircle(double r)
{
	double phi = 2*M_PI*RANDOM;
	return vec3d(r*cos(phi), 0, r*sin(phi));
}

template<typename T> T readStdinDefault(T default_value, string message)
{
	string input;
	T value;
	cerr << message << "[" << default_value << "]: ";
	getline(cin, input);
	if (!input.empty()) {
		istringstream stream( input );
		stream >> value;
	} else {
		value = default_value;
	}
	cerr << value << endl;
	return value;
}

void GenerateInitConfig::setParameters()
{
	/*
	 *  Read parameters from standard input
	 *
	 */
    MTRand mar;
    int np1_temp = 5000, np2_temp = 5000,i;
    double std1 = 0.0,std2= 0.0,r1_temp[4999],r2_temp[4999],r1_eq,r2_eq;
    double *r1,*r2;
    double sum =0.0, l_mu1 = 0, l_mu2=0, l_std1=0, l_std2 = 0;
	np = readStdinDefault(500, "number of particle");
	int dimension = readStdinDefault(2, "dimension (2 or 3)");
	if (dimension == 2) {
		sys.twodimension = true;
	} else {
		sys.twodimension = false;
	}
	if (sys.twodimension) {
		volume_fraction = readStdinDefault(0.78, "volume_fraction");
	} else {
		volume_fraction = readStdinDefault(0.5, "volume_fraction");
	}
	lx_lz = readStdinDefault(1.0 , "Lx/Lz [1]: "); // default value needs to be float number.
	if (!sys.twodimension) {
		ly_lz = readStdinDefault(1.0 , "Ly/Lz [1]: "); // default value needs to be float number.
	}
	disperse_type = readStdinDefault('b' , "(m)onodisperse or (b)idisperse or (p)oly");
	a1 = 1;
	a2 = 1;
    std2 = 0.0;
	volume_fraction1 = volume_fraction; // mono
    vf_ratio = 1;   // mono
	if (disperse_type == 'b') {
		cerr << "a1 = 1.0" << endl;
		do {
			a2 = readStdinDefault(1.4 , "a2 (a2>a1)");
		} while (a2 < a1);
		do {
			vf_ratio = readStdinDefault(0.5, "volume fraction ratio of smaller particle");
		} while (vf_ratio < 0 || vf_ratio > 1);
	}
    else if (disperse_type == 'p')
    {
        poly_type = readStdinDefault('n' , "type of polydispersity: (n)ormal or (l)og-normal");
        polyn_means = readStdinDefault('t' , "number of means: (o)ne or (t)wo");
        if (polyn_means == 't')
        {
            cerr << "a1 = 1.0" << endl;
            do {
                a2 = readStdinDefault(1.4 , "a2 (a2>a1)");
            } while (a2 < a1);
            do {
                vf_ratio = readStdinDefault(0.5, "volume fraction ratio of smaller particle");
            } while (vf_ratio < 0 || vf_ratio > 1);
            
            do {
                std1 = readStdinDefault(0.1, "standard deviation of normal distribution of smaller particle  (<0.25)");
                std2 = readStdinDefault(0.1, "standard deviation of normal distribution of bigger particle (a2 - 3*std2 < 1 + 3*std1)");
            } while (std1 > 0.25 && (a2-3*std2) < (1.0+3*std1));
            
        }
        else if (polyn_means == 'o')
        {
            cerr << "a1 = 1.0" << endl;
            do {
                std1 = readStdinDefault(0.1, "standard deviation of particle (<0.25)");
            } while (std1 > 0.25);
            
        }
    }
    
    else {
		vf_ratio = 1;
	}
    
    if (poly_type == 'l')    //log-normal scales
    {
        l_mu1 = log(1.0/(sqrt(pow(std1,2) + 1.0)));
        l_std1 = sqrt(log(pow(std1,2)/1.0 + 1));
        l_mu2 = log(pow(a2,2)/(sqrt(pow(std2,2) + pow(a2,2))));
        l_std2 = sqrt(log(pow(std2,2)/pow(a2,2) + 1));
        printf("%le %le %le %le", l_mu1, l_std1, l_mu2, l_std2);
    }
    
	//rand_seed = readStdinDefault(1, "random seed");
	/*
	 *  Calculate parameters
	 */
	volume_fraction1 = volume_fraction*vf_ratio;
	if (disperse_type == 'b') {
		volume_fraction2 = volume_fraction-volume_fraction1;
	}
    else if (disperse_type == 'p')
    {
        if (polyn_means == 't')
        {
            volume_fraction2 = volume_fraction-volume_fraction1;
            for (i=0;i<np1_temp;i++)
            {
                if (poly_type == 'n')
                    r1_temp[i] = ceil(100*mar.randNorm(1.0,std1))/100;
                else if (poly_type == 'l')
                    r1_temp[i] = ceil(100*exp(l_mu1 + l_std1*mar.randNorm(0.0,1.0)))/100;
                sum = sum + pow(r1_temp[i],3);
            }
            r1_eq = pow(sum/np1_temp,1.0/3);
            sum = 0.0;
            for (i=1;i<np2_temp;i++)
            {
                if (poly_type == 'n')
                    r2_temp[i] = ceil(100*mar.randNorm(a2,std2))/100;
                else if (poly_type == 'l')
                    r2_temp[i] = ceil(100*exp(l_mu2 + l_std2*mar.randNorm(0.0,1.0)))/100;
                sum = sum + pow(r2_temp[i],3);
            }
            r2_eq = pow(sum/np2_temp,1.0/3);
            //        cerr << "r2eq = " << np1+np2 << endl;
            cerr << "r1_eq : r2_eq " << r1_eq << ":" << r2_eq << endl;
        }
        else {
            volume_fraction2 = 0;
            for (i=0;i<np1_temp;i++)
            {
                if (poly_type == 'n')
                    r1_temp[i] = ceil(100*mar.randNorm(1.0,std1))/100;
                else if (poly_type == 'l')
                    r1_temp[i] = ceil(100*exp(l_mu1 + l_std1*mar.randNorm(0.0,1.0)))/100;
                sum = sum + pow(r1_temp[i],3);
              
            }
//              printf("volume fraction 1 = %f \n",volume_fraction1);
            r1_eq = pow(sum/np1_temp,1.0/3);
            r2_eq = r1_eq;
            cerr << "r1_eq : r2_eq " << r1_eq << ":" << r2_eq << endl;
            printf("std 1 = %f \n",std1);
        }
        
    }
    else {
		volume_fraction2 = 0;
	}
	double total_volume;
	double pvolume1, pvolume2;
	if (sys.twodimension) {
		pvolume1 = M_PI*a1*a1;
		pvolume2 = M_PI*a2*a2;
	} else {
		pvolume1 = (4.0/3)*M_PI*a1*a1*a1;
		pvolume2 = (4.0/3)*M_PI*a2*a2*a2;
	}
    
    a_12 = (double *) calloc(np, sizeof(double));
    
    if (disperse_type == 'p')
    {
        // if (polyn_means == 't')
        //{
                if (sys.twodimension) {
            pvolume1 = M_PI*r1_eq*r1_eq;
            pvolume2 = M_PI*r2_eq*r2_eq;
        } else {
            pvolume1 = (4.0/3)*M_PI*r1_eq*r1_eq*r1_eq;
            pvolume2 = (4.0/3)*M_PI*r2_eq*r2_eq*r2_eq;
        }
        
        if (np > 0)
        {
            total_volume = np/(volume_fraction1/pvolume1+volume_fraction2/pvolume2);
            printf("total volume %f \n",total_volume);
            np1 = floor(volume_fraction1*total_volume/pvolume1 + 0.5);
            np2 = np - np1;
            sum = 0.0;
            r1 = (double *) calloc(np1, sizeof(double));
            for (i=0;i<np1;i++)
            {
                if (poly_type == 'n')
                    r1[i] = ceil(100*mar.randNorm(1.0,std1))/100;
                else if (poly_type == 'l')
                    r1[i] = ceil(100*exp(l_mu1 + l_std1*mar.randNorm(0.0,1.0)))/100;
                sum = sum + 4.0/3*M_PI*pow(r1[i],3);
            }
            volume_fraction1 = sum/total_volume;
            sum = 0.0;
            r2 = (double *) calloc(np2, sizeof(double));
            for (i=0;i<np2;i++)
            {
                if (poly_type == 'n')
                    r2[i] = ceil(100*mar.randNorm(a2,std2))/100;
                else if (poly_type == 'l')
                    r2[i] = ceil(100*exp(l_mu2 + l_std2*mar.randNorm(0.0,1.0)))/100;
                sum = sum + 4.0/3*M_PI*pow(r2[i],3);
            }
            volume_fraction2 = sum/total_volume;
            
            //double np1_tmp = volume_fraction1*total_volume/pvolume1;
            
            double pvolume = np1*pvolume1+np2*pvolume2;
            if (sys.twodimension) {
                lz = sqrt(pvolume/(lx_lz*volume_fraction));
                lx = lz*lx_lz;
                ly = 0;
            } else {
                lz = pow(pvolume/(lx_lz*ly_lz*volume_fraction), 1.0/3);
                lx = lz*lx_lz;
                ly = lz*ly_lz;
            }
        } else
        {
            lz = readStdinDefault(10, "lz");
            lx = lz*lx_lz;
            ly = lz*ly_lz;
            double pvolume1_ = lx*ly*lz*volume_fraction1;
            double pvolume2_ = lx*ly*lz*volume_fraction2;
            np1 = (int)(pvolume1_/pvolume1+0.5);
            np2 = (int)(pvolume2_/pvolume2+0.5);
            np = np1 + np2;
        }
        
        for(i=0;i<np;i++)
        {
            if (i<np1)
                a_12[i]=r1[i];
            else
                a_12[i]=r2[i-np1];
        }

        //}
    }
    
    else
    {
        if (np > 0) {
            total_volume = np/(volume_fraction1/pvolume1+volume_fraction2/pvolume2);
            double np1_tmp = volume_fraction1*total_volume/pvolume1;
            if (np1_tmp-(int)np1_tmp <= 0.5) {
                np1 = (int)np1_tmp;
            } else {
			np1 = (int)np1_tmp+1;
            }
            np2 = np-np1;
            double pvolume = np1*pvolume1+np2*pvolume2;
            if (sys.twodimension) {
                lz = sqrt(pvolume/(lx_lz*volume_fraction));
                lx = lz*lx_lz;
                ly = 0;
            } else {
                lz = pow(pvolume/(lx_lz*ly_lz*volume_fraction), 1.0/3);
                lx = lz*lx_lz;
                ly = lz*ly_lz;
            }
        } else {
            lz = readStdinDefault(10, "lz");
            lx = lz*lx_lz;
            ly = lz*ly_lz;
            double pvolume1_ = lx*ly*lz*volume_fraction1;
            double pvolume2_ = lx*ly*lz*volume_fraction2;
            np1 = (int)(pvolume1_/pvolume1+0.5);
            np2 = (int)(pvolume2_/pvolume2+0.5);
            np = np1 + np2;
        }
        for(i=0;i<np;i++)
        {
            if (i<np1)
                a_12[i]=a1;
            else
                a_12[i]=a2;
        }
    }
    
    
	lx_half = lx/2;
	ly_half = ly/2;
	lz_half = lz/2;
	cerr << "np = " << np1+np2 << endl;
	cerr << "np1 : np2 " << np1 << ":" << np2 << endl;
	cerr << "vf1 = " << volume_fraction << endl;
	cerr << "vf2 = " << volume_fraction2 << endl;
	cerr << "box =" << lx << ' ' << ly << ' ' << lz << endl;
}
