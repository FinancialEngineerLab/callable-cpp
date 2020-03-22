#include "product.hpp"

namespace beagle
{
  namespace product
  {
    Product::~Product(void)
    { }

    namespace mixins
    {
      Option::~Option(void)
      { }

      Bond::~Bond(void)
      { }

      double Bond::standardFaceValue(void) const
      {
        return 100.0;
      }
    }
  }
}