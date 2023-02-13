#include <iostream>

#include "Hansen/HansenProblem.hpp"
#include "Hansen/HansenProblemFamily.hpp"

#include "Hill/HillProblem.hpp"
#include "Hill/HillProblemFamily.hpp"

#include "Shekel/ShekelProblem.hpp"
#include "Shekel/ShekelProblemFamily.hpp"

#include "Other/OptSqConstrProblem.hpp"

#include "Grishagin/grishagin_function.hpp"
#include "Grishagin/GrishaginProblemFamily.hpp"
#include "Grishagin/GrishaginConstrainedProblem.hpp"
#include "Grishagin/GrishaginConstrainedProblemFamily.hpp"

#include "GKLS/GKLSProblem.hpp"
#include "GKLS/GKLSProblemFamily.hpp"
#include "GKLS/GKLSConstrainedProblem.hpp"
#include "GKLS/GKLSConstrainedProblemFamily.hpp"

using std::cout;

int main()
{
  try
  {
    TGKLSProblem gkls;
    cout << "\n  GKLS Problem 1" << std::endl;
    cout << "GKLSProblem ( 0.5, 0.5 ) = " << gkls.ComputeFunction({ 0.5, 0.5 }) << std::endl;
    cout << "GKLSProblem Derivatives ( 0.5, 0.5 ) = {" <<
      gkls.ComputeFunctionDerivatives({ 0.5, 0.5 })[0] << ", " <<
      gkls.ComputeFunctionDerivatives({ 0.5, 0.5 })[1] << "}" << std::endl;
    cout << std::endl;

    TGKLSProblemFamily gklsFam;
    cout << "  GKLS Problem Family" << std::endl;
    for (int i = 0; i < gklsFam.GetFamilySize(); i++)
      cout << "GKLSProblem [" << i + 1 << "] ( 0.5, 0.5 ) = " <<
      gklsFam[i]->ComputeFunction({ 0.5, 0.5 }) << std::endl;

    {
      // create a problem with 30% fraction of feasible domain with respect to the whole search domain
      TGKLSConstrainedProblem gklsConst(cptInFeasibleDomain, 0.3, 0, 1);
      cout << "  GKLS Constrained Problem 1" << std::endl;
      // get the global minimizer coordinates
      std::vector<double> y = gklsConst.GetOptimumPoint();
      std::cout << "min GKLSConst(" << y[0] << ", " << y[1] << ") = " << gklsConst.ComputeFunction(y) << std::endl;

      cout << "OptimumValue = " << gklsConst.GetOptimumValue() << std::endl;
      cout << "OptimumPoint = {" << gklsConst.GetOptimumPoint()[0] << ", " << gklsConst.GetOptimumPoint()[1] << "}" << std::endl;
      cout << "Dimension = " << gklsConst.GetDimension() << std::endl;
      cout << "ConstraintsNumber = " << gklsConst.GetConstraintsNumber() << std::endl;
      cout << "ComputeFunction [0.5, 0.5] = " << gklsConst.ComputeFunction({ 0.5, 0.5 }) << std::endl;
      int index = 0;
      cout << "cctAllConstraints ConstraintsNumber [0.5, 0.5] = {" <<
        gklsConst.ComputeConstraints({ 0.5, 0.5 }, cctAllConstraints,
          index)[0] << ", " << gklsConst.ComputeConstraints({ 0.5, 0.5 }, cctAllConstraints,
            index)[1] << "}" << std::endl;
      cout << "ComputeConstraint 0 = " << gklsConst.ComputeConstraint(index, { 0.5, 0.5 })
        << std::endl;
    }

    TGKLSConstrainedProblemFamily gklsConstFam;
    cout << "  GKLS Constrained Problem Family" << std::endl;
    for (int i = 0; i < gklsConstFam.GetFamilySize(); i++)
      cout << "GKLSProblem [" << i + 1 << "] ( 0.5, 0.5 ) = " <<
      gklsConstFam[i]->ComputeFunction({ 0.5, 0.5 }) << std::endl;

    TGrishaginProblem grish;
    cout << "  Grishagin Problem 1" << std::endl;
    cout << "GrishaginProblem ( 0.5, 0.5 ) = " << grish.ComputeFunction({ 0.5, 0.5 }) << std::endl;
    cout << "GrishaginProblem Derivatives ( 0.5, 0.5 ) = {" <<
      grish.ComputeFunctionDerivatives({ 0.5, 0.5 })[0] << ", " <<
      grish.ComputeFunctionDerivatives({ 0.5, 0.5 })[1] << "}" << std::endl;
    cout << std::endl;

    TGrishaginProblemFamily grishFam;
    cout << "  Grishagin Problem Family" << std::endl;
    for (int i = 0; i < grishFam.GetFamilySize(); i++)
      cout << "GrishaginProblem [" << i + 1 << "] ( 0.5, 0.5 ) = " <<
      grishFam[i]->ComputeFunction({ 0.5, 0.5 }) << std::endl;

    {
      // create a problem with 30% fraction of feasible domain with respect to the whole search domain
      GrishaginConstrainedProblem grishConst(cptInFeasibleDomain, 0.3, 0, 1);
      cout << "  Grishagin Constrained Problem 1" << std::endl;
	  // get the global minimizer coordinates
	  std::vector<double> y = grishConst.GetOptimumPoint();
      std::cout << "min grishConst(" << y[0] << ", " << y[1] << ") = " << grishConst.ComputeFunction(y) << std::endl;

      cout << "OptimumValue = " << grishConst.GetOptimumValue() << std::endl;
      cout << "OptimumPoint = {" << grishConst.GetOptimumPoint()[0] << ", " << grishConst.GetOptimumPoint()[1] << "}" << std::endl;
      cout << "Dimension = " << grishConst.GetDimension() << std::endl;
      cout << "ConstraintsNumber = " << grishConst.GetConstraintsNumber() << std::endl;
      cout << "ComputeFunction [0.5, 0.5] = " << grishConst.ComputeFunction({ 0.5, 0.5 }) << std::endl;
      int index = 0;
      cout << "cctAllConstraints ConstraintsNumber [0.5, 0.5] = {" <<
        grishConst.ComputeConstraints({ 0.5, 0.5 }, cctAllConstraints,
          index)[0] << ", " << grishConst.ComputeConstraints({ 0.5, 0.5 }, cctAllConstraints,
            index)[1] << "}" << std::endl;
      cout << "ComputeConstraint 0 = " << grishConst.ComputeConstraint(index, { 0.5, 0.5 })
        << std::endl;
    }

    TGrishaginConstrainedProblemFamily grishConstFam;
    cout << "  Grishagin Constrained Problem Family" << std::endl;
    for (int i = 0; i < grishConstFam.GetFamilySize(); i++)
      cout << "GrishaginProblem [" << i + 1 << "] ( 0.5, 0.5 ) = " <<
      grishConstFam[i]->ComputeFunction({ 0.5, 0.5 }) << std::endl;

    THansenProblem0 hans;
    cout << "  Hansen Problem 0" << std::endl;
    cout << "HansenProblem ( 0.0 ) = " << hans.ComputeFunction({ 0 }) << std::endl;
    cout << "HansenProblem Derivatives ( 0.0 ) = " <<
      hans.ComputeFunctionDerivatives({ 0 })[0] << std::endl;
    cout << std::endl;

    THansenProblemFamily hansFam;
    cout << "  Hansen Problem Family" << std::endl;
    for (int i = 0; i < hansFam.GetFamilySize(); i++)
      cout << "HansenProblem [" << i << "] (0) = " <<
      hansFam[i]->ComputeFunction({ 0 }) << std::endl;

    THillProblem hill;
    cout << "  Hil Problem 0" << std::endl;
    cout << "HilProblem ( 0.6 ) = " << hill.ComputeFunction({ 0.6 }) << std::endl;
    cout << "HilProblem Derivatives ( 0.6 ) = " <<
      hill.ComputeFunctionDerivatives({ 0.6 })[0] << std::endl;
    cout << std::endl;

    THillProblemFamily hillFam;
    cout << "  Hill Problem Family" << std::endl;
    for (int i = 0; i < hillFam.GetFamilySize(); i++)
      cout << "HillProblem [" << i << "] (" << 0.5 + double(i) / 2000.0 << ") = " <<
      hillFam[i]->ComputeFunction({ 0.5 + double(i) / 2000.0 }) << std::endl;

    TShekelProblem shekel;
    cout << "  Shekel Problem 0" << std::endl;
    cout << "ShekelProblem ( 0.6 ) = " << shekel.ComputeFunction({ 0.6 }) << std::endl;
    cout << "ShekelProblem Derivatives ( 0.6 ) = " <<
      shekel.ComputeFunctionDerivatives({ 0.6 })[0] << std::endl;
    cout << std::endl;

    TShekelProblemFamily shekelFam;
    cout << "  Shekel Problem Family" << std::endl;
    for (int i = 0; i < shekelFam.GetFamilySize(); i++)
      cout << "ShekelProblem [" << i << "] (" << 0.5 + double(i) / 2000.0 << ") = " <<
      shekelFam[i]->ComputeFunction({ 0.5 + double(i) / 2000.0 }) << std::endl;

    {
      cout << "\nsqConstrProblem cptNormal" << std::endl;
      TOptSqConstrProblem sqConstrProblem(2, { -10, -10 }, { 10, 10 }, { 0, 0 }, 0,
        cptNormal, 1.0);
      cout << "OptimumValue = " << sqConstrProblem.GetOptimumValue() << std::endl;
      cout << "OptimumPoint = " << sqConstrProblem.GetOptimumPoint()[0] << std::endl;
      cout << "Dimension = " << sqConstrProblem.GetDimension() << std::endl;
      cout << "ConstraintsNumber = " << sqConstrProblem.GetConstraintsNumber() << std::endl;
      cout << "ComputeFunction [0, 0] = " << sqConstrProblem.ComputeFunction({ 0, 0 }) << std::endl;
      int index = 0;
      cout << "cctAllConstraints ConstraintsNumber [1, 1] = " <<
        sqConstrProblem.ComputeConstraints({ 1, 1 }, cctAllConstraints,
          index)[0] << std::endl;
      cout << "ComputeConstraint 0 = " << sqConstrProblem.ComputeConstraint(index, { 0.5, 0.5 })
        << std::endl;
    }
    {
      cout << "\nsqConstrProblem cptInFeasibleDomain" << std::endl;
      TOptSqConstrProblem sqConstrProblem(2, { -10, -10 }, { 10, 10 }, { 0, 0 }, 0,
        cptInFeasibleDomain, 0.5);
      cout << "OptimumValue = " << sqConstrProblem.GetOptimumValue() << std::endl;
      cout << "OptimumPoint = " << sqConstrProblem.GetOptimumPoint()[0] << std::endl;
      cout << "Dimension = " << sqConstrProblem.GetDimension() << std::endl;
      cout << "ConstraintsNumber = " << sqConstrProblem.GetConstraintsNumber() << std::endl;
      cout << "ComputeFunction [0, 0] = " << sqConstrProblem.ComputeFunction({ 0, 0 }) << std::endl;
      int index = 0;
      cout << "cctAllConstraints ConstraintsNumber [1, 1] = " <<
        sqConstrProblem.ComputeConstraints({ 1, 1 }, cctAllConstraints,
          index)[0] << std::endl;
      cout << "ComputeConstraint 0 = " << sqConstrProblem.ComputeConstraint(index, { 0.5, 0.5 })
        << std::endl;
    }
    {
      cout << "\nsqConstrProblem cptOutFeasibleDomain" << std::endl;
      TOptSqConstrProblem sqConstrProblem(2, { -10, -10 }, { 10, 10 }, { 0, 0 }, 0,
        cptOutFeasibleDomain, 0.5);
      cout << "OptimumValue = " << sqConstrProblem.GetOptimumValue() << std::endl;
      cout << "OptimumPoint = " << sqConstrProblem.GetOptimumPoint()[0] << std::endl;
      cout << "Dimension = " << sqConstrProblem.GetDimension() << std::endl;
      cout << "ConstraintsNumber = " << sqConstrProblem.GetConstraintsNumber() << std::endl;
      cout << "ComputeFunction [0, 0] = " << sqConstrProblem.ComputeFunction({ 0, 0 }) << std::endl;
      int index = 0;
      cout << "cctAllConstraints ConstraintsNumber [1, 1] = " <<
        sqConstrProblem.ComputeConstraints({ 1, 1 }, cctAllConstraints,
          index)[0] << std::endl;
      cout << "ComputeConstraint 0 = " << sqConstrProblem.ComputeConstraint(index, { 0.5, 0.5 })
        << std::endl;
    }
    {
      cout << "\nsqConstrProblem cptOnFeasibleBorder" << std::endl;
      TOptSqConstrProblem sqConstrProblem(2, { -10, -10 }, { 10, 10 }, { 0, 0 }, 0,
        cptOnFeasibleBorder, 0.5);
      cout << "OptimumValue = " << sqConstrProblem.GetOptimumValue() << std::endl;
      cout << "OptimumPoint = " << sqConstrProblem.GetOptimumPoint()[0] << std::endl;
      cout << "Dimension = " << sqConstrProblem.GetDimension() << std::endl;
      cout << "ConstraintsNumber = " << sqConstrProblem.GetConstraintsNumber() << std::endl;
      cout << "ComputeFunction [0, 0] = " << sqConstrProblem.ComputeFunction({ 0, 0 }) << std::endl;
      int index = 0;
      cout << "cctAllConstraints ConstraintsNumber [1, 1] = " <<
        sqConstrProblem.ComputeConstraints({ 1, 1 }, cctAllConstraints,
          index)[0] << std::endl;
      cout << "ComputeConstraint 0 = " << sqConstrProblem.ComputeConstraint(index, { 0.5, 0.5 })
        << std::endl;
    }
  }

  catch (string s)
  {
    cout << s << std::endl;
  }
  system("pause");

  return 0;
}
