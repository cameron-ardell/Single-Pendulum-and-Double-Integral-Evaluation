#include <iostream>
#include "nr3.h"
#include "quadrature.h"
#include "interp_1d.h"
#include "romberg.h"

using namespace std;


struct Integrand {
  Doub theta_o;
	
  //Constructor
  Integrand(Doub theta_in) { theta_o = theta_in; }

  //Destructor
  ~Integrand() { cout << "Destroying Integrand..." << endl; }

  //Set theta_o (since it'll be changing each integral calculated)
  Doub set_theta_o(Doub theta_in) {return theta_o = theta_in; };

  //actual function from assignment (without coefficients)
  Doub operator()(Doub theta) { return 1.0 / ( sqrt(cos(theta) -cos(theta_o))); }

};


int main() {
  Doub pi = 3.14159265358979;

  //highest value of theta_o I'm going to use
  Doub limit = pi / 2.0;

  //lower limit of integration (a is always 0)
  Doub a = 0.0;
  //upper limit of integration changes
  Doub b;

  //create output file
  ofstream file;
  file.open("dataprob1.dat");
  file << "#Theta_0         T/T_0" << endl;


  //establish integrand object
  Integrand integ(0);

  //need variable to hold answer to integral
  //also need variable to be T/ T_o, which is just the integral
  //times sqrt(2)/pi
  Doub answer;
  Doub actual_answer;

  //iterates for different theta_o's up to limit, i being theta_o
  for(Doub i = 0.01; i < limit; i += 0.01) {

    //need to convert i to radians for function since cos
    //is in radians
    // b = i * 180.0 / pi;
    b = i;

    integ.set_theta_o(b);

    //need to do midpoint on integrand in order to get integral,
    //can't include upper limit since it is undefined there
    Midsqu <Integrand> integrator(integ, a, b);

    //upper limit of integration will be current theta_o value
    answer = qromo(integrator);

    //get T/T_o and then write it all to a file
    actual_answer = answer * sqrt(2) / pi;
    file << setw(12) << i << setw(10) << actual_answer << endl;
  }

  file.close();


  return 0;
}
