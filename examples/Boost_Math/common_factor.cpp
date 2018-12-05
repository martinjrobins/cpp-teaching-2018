#include <boost/math/common_factor.hpp>
#include <iostream>

int main()
{
   using std::cout;
   using std::endl;

   cout << "The GCD and LCM of 6 and 15 are "
   << boost::math::gcd(6, 15) << " and "
   << boost::math::lcm(6, 15) << ", respectively."
   << endl;
}
