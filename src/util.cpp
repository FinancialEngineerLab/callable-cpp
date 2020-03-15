#include "util.hpp"
#include <Eigen/Dense>

namespace beagle
{
  namespace util
  {
    beagle::dbl_t pi( void )
    {
      static const beagle::dbl_t PI = 3.14159265358979323846;
      return PI;
    }

    beagle::dbl_t rootTwo( void )
    {
      static const beagle::dbl_t rootTwo = 1.41421356237;
      return rootTwo;
    }

    beagle::dbl_t standardNormal( beagle::dbl_t arg )
    {
      return std::exp( -.5*arg*arg ) / std::sqrt( 2 * pi() );
    }

    beagle::dbl_t cumulativeStandardNormal( beagle::dbl_t arg )
    {
      return .5 * ( 1. + std::erf( arg / rootTwo() ) );
    }

    beagle::dbl_t bsCall( beagle::dbl_t strike,
                   beagle::dbl_t forward,
                   beagle::dbl_t expiry,
                   beagle::dbl_t vol )
    {
      beagle::dbl_t moneyness = std::log( forward / strike );
      beagle::dbl_t totalDev = vol * std::sqrt( expiry );
      beagle::dbl_t dOne = moneyness / totalDev + .5 * totalDev;
      beagle::dbl_t dTwo = dOne - totalDev;

      return forward * cumulativeStandardNormal( dOne ) - strike * cumulativeStandardNormal( dTwo );
    }

    beagle::dbl_t bsVega( beagle::dbl_t strike,
                   beagle::dbl_t forward,
                   beagle::dbl_t expiry,
                   beagle::dbl_t vol )
    {
      beagle::dbl_t moneyness = std::log( forward / strike );
      beagle::dbl_t totalDev = vol * std::sqrt( expiry );
      beagle::dbl_t dOne = moneyness / totalDev + .5 * totalDev;

      return forward * standardNormal( dOne ) * std::sqrt( expiry );
    }

    void tridiagonalSolve( beagle::dbl_vec_t& rhs,
                           beagle::dbl_vec_t& diag,
                           beagle::dbl_vec_t& upper,
                           beagle::dbl_vec_t& lower )
    {
      int size = rhs.size();
      for (int i = 0; i < size - 1; ++i)
      {
        diag[i+1] -= lower[i+1] / diag[i] * upper[i];
        rhs[i+1] -= lower[i+1] / diag[i] * rhs[i];
      }

      rhs[size-1] /= diag[size-1];

      for (int i = 0; i < size - 1; ++i)
      {
        int j = size - 2 - i;
        rhs[j] = (rhs[j] - upper[j] * rhs[j+1]) / diag[j];
      }

      //int size = rhs.size();

      //diag[1] -= lower[1] / diag[0] * upper[0];
      //upper[1] -= lower[1] / diag[0] * lower[0];
      //rhs[1] -= lower[1] / diag[0] * rhs[0];

      //lower[size - 2] -= upper[size - 2] / diag[size - 1] * upper[size - 1];
      //diag[size - 2] -= upper[size - 2] / diag[size - 1] * lower[size - 1];
      //rhs[size - 2] -= upper[size - 2] / diag[size - 1] * rhs[size - 1];

      //for (int i = 1; i < size - 2; ++i)
      //{
      //  diag[i + 1] -= lower[i + 1] / diag[i] * upper[i];
      //  rhs[i + 1] -= lower[i + 1] / diag[i] * rhs[i];
      //}

      //rhs[size - 2] /= diag[size - 2];

      //for (int i = 1; i < size - 2; ++i)
      //{
      //  int j = size - 2 - i;
      //  rhs[j] = (rhs[j] - upper[j] * rhs[j + 1]) / diag[j];
      //}

      //rhs[0] = (rhs[0] - upper[0] * rhs[1] - lower[0] * rhs[2]) / diag[0];
      //rhs[size - 1] = (rhs[size - 1] - upper[size - 1] * rhs[size - 3] - lower[size - 1] * rhs[size - 2]) / diag[size - 1];
    }

    void inverseMatrix( beagle::dbl_vec_vec_t& inputMat )
    {
      beagle::dbl_vec_vec_t::size_type numRows = inputMat.size();
      for (const auto row : inputMat)
      {
        // The input must be a square matrix in order to perform inversion
        if (row.size() != numRows)
          throw(std::string("The input matrix is not a square matrix!"));
      }

      using dbl_mat_t = Eigen::MatrixXd;
      dbl_mat_t theMatrix(numRows, numRows);
      for (beagle::dbl_vec_vec_t::size_type row=0; row<numRows; ++row)
      {
        for (beagle::dbl_vec_vec_t::size_type col=0; col<numRows; ++col)
          theMatrix(row, col) = inputMat[row][col];
      }

      theMatrix = theMatrix.inverse();
      for (beagle::dbl_vec_vec_t::size_type row=0; row<numRows; ++row)
      {
        for (beagle::dbl_vec_vec_t::size_type col=0; col<numRows; ++col)
          inputMat[row][col] = theMatrix(row, col);
      }
    }
  }
}
