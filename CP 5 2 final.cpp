#include <iostream>
#include "nr3.h"
#include "interp_1d.h"
#include "quadrature.h"
#include "romberg.h"
#include <time.h>
#include <cstdlib>
#include "ran.h"
#include <string>

using namespace std;

struct Model {
	//dimensions of spin array
	int m;
	//total number of spins
	int N;
	//dimensions
	int dim;
	//temperature
	Doub T;
	//Boltzmann constant
	Doub boltz;
	//file to write data to
	string file;
	//array of spins
	MatInt spin;
	//energy array for all local energies
	MatInt loc_nrg;
	//average energy of state
	Doub e_avg;
	//average magentization of state
	Doub mag_avg;


	Model(int m_in, Doub T_in, string fileName) { m = m_in;
		N = m_in * m_in; T = T_in; file = fileName; boltz = 1.0 / T_in;}
	~Model() { cout << "destructing Model...\n";}


	//creates spin and energy array with correct dimensions
	void initialize(){
		e_avg = 0.0;
		mag_avg = 0.0;
		//all spins set to 1 initially
		spin.assign(m,m,1);
		//energies all set to 0, though will be modified immediately
		loc_nrg.assign(m,m,0);
	}

	//updates array of local energy
	void update(){
		int tot = 0;
		int spin_sum = 0;

		//goes through every row and collum
		for(int r = 0; r < m; r++){
			for(int c = 0; c < m; c++){

				int spinVal = spin[r][c];
				spin_sum += spinVal;

				//finding neighbors while accounting
				//for boundary conditions

				//index to left
				int rleft = r - 1;
				if(rleft == -1){
					rleft = m-1;
				}

				//index to right
				int rright = r + 1;
				if(rright == m){
					rright = 0;
				}

				//index above
				int cup = c - 1;
				if(cup == -1){
					cup = m-1;
				}

				//index below
				int cdown = c + 1;
				if(cdown == m){
					cdown = 0;
				}

				//calculating local energy related to total energy
				int tot_comp = -spinVal*(spin[rright][c] + spin[r][cdown]);
				tot += tot_comp;

				//need other two neighbors spin to get local energy
				int other_comp = -spinVal*(spin[rleft][c] + spin[r][cup]);

				//sums all nearest neighbor interatctions * -1 (local energy)
				loc_nrg[r][c] = tot_comp + other_comp;
			}
		}
		//once all local energies have been evaluated, gives total energy of state
		e_avg = Doub(tot)/N;
		//once all spins have been summed, can calculate mag average
		mag_avg = Doub(spin_sum)/N;
	}


	//print spin array
	void print_spin_array() {
		for(int r = 0; r < m; r++){
			for(int c = 0; c < m; c++){
				cout << " " << spin[r][c];
			}
			cout << endl;
		}
	}

	//print energy array
	void print_energy_array() {
		for(int r = 0; r < m; r++){
			for(int c = 0; c < m; c++){
				cout << " " << loc_nrg[r][c];
			}
			cout << endl;
		}
	}

	 void test() {
	 	cout << "mag average" << mag_avg << endl;
	 	cout << "e average" << e_avg << endl;
	 }

	 Doub standardDev(Doub x) {
	 	Doub sqrd = (x*x)/N;
	 	Doub meanSqrd = sqrd/N;
	 	cout << sqrd << " " << meanSqrd << endl;
	 	Doub stDv = sqrt(sqrd - meanSqrd);
	 	cout << stDv << endl;
	 	return stDv;
	 }


	 //selects random spins and potentially flips (more likely if
	 //change in energy is greater)
	void sweep(int num_sweeps, int equib) {
		//seed random number generator
		time_t timer = time(0);
		Doub ran_seed = Doub(timer);
		Ran rand(ran_seed);

		ofstream data;
		data.open(file.c_str());
		data << "#sweep    magentization   standardDevMag     energy   standardDevEnergy" << endl;

		Doub magSum = 0.0;
		Doub nrgSum = 0.0;
		Doub pastEquib = 0.0;

		for(int i = 0; i < num_sweeps; i++){


			for(int l = 0; l < N; l++){

				//picks random row and colum
				int r = int(rand.doub()*m);
				int c = int(rand.doub()*m);

				//gets local energy for that spin
				int en = loc_nrg[r][c];

				//dE is local energy val flipped - local energy
				int dE = Doub(-en - en);

				//need random number
				Doub rand_num = rand.doub();

				//value to compare to 
				Doub comparison = exp(-dE*(1.0/T));

				//if random number is smaller than e f'n, flip
				if(rand_num < comparison){
					//cout << "flipping happened\n";
					int prevVal = spin[r][c];
					int newVal = -prevVal;

					spin[r][c] = newVal;
				}
			}
			update();
			if(i > equib){
				magSum = mag_avg*N;
				nrgSum = e_avg*N;

				data << i << " " << mag_avg << " " << standardDev(magSum)
					<< " " << e_avg << " " << standardDev(nrgSum) <<endl;
			}

		}

		data << num_sweeps << " " << mag_avg << " " << e_avg << endl;
		data.close();
	}
};



int main() {

	int dim = 256;

	Doub T;
	cout << "temperature?\n";
	cin >> T;

	string fileName;
	cout << "file name?\n";
	cin >> fileName;

	int swept;
	cout << "num sweeps?\n";
	cin >> swept;

	int equib;
	cout << "sweeps to reach equilibrium?\n";
	cin >> equib;

	Model mod(dim, T, fileName);
	mod.initialize();
	mod.update();
	mod.sweep(swept, equib);



	return 0;
}