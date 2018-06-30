#include "pricer.hpp"

namespace beagle
{
  namespace valuation
  {
    Pricer::Pricer( void )
    { }

    Pricer::~Pricer( void )
    { }

    namespace mixins
    {
      OptionValueCollectionProvider::~OptionValueCollectionProvider( void )
      { }
    }
  }
}