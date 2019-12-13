#ifndef ROUND_TRIP_HPP
#define ROUND_TRIP_HPP

#include "fwd_decl.hpp"

namespace beagle
{
  namespace test
  {
    void generateEuropeanMarketQuotes( const beagle::valuation::FiniteDifferenceDetails& fdDetails,
                                       beagle::dbl_vec_t& expiries,
                                       beagle::dbl_vec_vec_t& strikesColl,
                                       beagle::dbl_vec_vec_t& pricesColl );
  }
}



#endif