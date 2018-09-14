#pragma once

#include "IGeneralOptProblem.hpp"


// Methods of forming a problem with constraints
enum EConstrainedProblemType {
  cptNormal, cptInFeasibleDomain, cptOutFeasibleDomain,
  cptOnFeasibleBorder
};

// Schemes for calculating constraints
enum EConstraintComputationType { cctAllConstraints, cctIndexScheme };

class IConstrainedOptProblem : public IGeneralOptProblem
{
protected:
  // Type of constrained problem
  EConstrainedProblemType mProblemType;
  // Feasible domain fraction
  double mFeasibleDomainFraction;
  // Number of active constraints  (only for OnFeasibleBorder)
  // By default = 0, i.e. parameter is not set
  int mActiveConstraintNumber;

  /// Shift for the constraints (RHS)
  std::vector<double> mQ;
  /// Scaling factors for constraints
  std::vector<double> mZoomRatios;
  /// Coordinate shift of constraints
  std::vector<std::vector<double> > mShift;
  /// Coordinate shift of the global minimizer to the boundary point
  std::vector<double> mBoundaryShift;
  /// Improvement coefficients
  std::vector<double> mImprovementCoefficients;
  /// The accuracy of the search for the nearest point on the domain boundary (power of 0.5)
  int mBoundarySearchPrecision;


  /// Is there a task scaling
  bool mIsZoom;
  /// Is there a constaint shifting
  bool mIsShift;
  /// Is there a optimizer shifting 
  bool mIsBoundaryShift;
  /// Is there a common fraction of feasible domain
  bool mIsTotalDelta;
  /// Is there an objective funstion scaling
  bool mIsImprovementOfTheObjective;

  void InitProblem();

  // Perform the transformation of the objective function value
  double TransformValue(double val, vector<double> point, int index = 0) const;
  // Perform the coordinate transformation of the objective function points
  vector<double> TransformPoint(vector<double> point, int index = 0) const;
  // Perform the transformation of the constraint number index
  double TransformConstraintValue(double val, int index) const;
  // Perform the transformation of the constraint number index point
  vector<double> TransformConstraintPoint(vector<double> point, int index) const;

  /// Calculate maximum of thу generalized constraint 
  double MaxFunctionCalculate(std::vector<double> y);

  /** Computs the shift of constraints by the fraction of the feasible domain
  \param[in] delta - fraction of feasible domain
  \param[in] m - number of grid nodes in computing the feasible domain fraction 
  \param[in] Epsilon - accuraсy of the feasible domain fraction computing
  \param[in] maxM - maximum number of trials
  \return shift values for the constraints
  */
  double CalculateRHS(double delta, int m = 100, double Epsilon = 0.01, int maxM = 10000000);

  /// Scaling factors for the constraints
  virtual void SetZoom();
  /// Shift to the global minimizer for the constraints
  virtual void SetShift();
  /// Shift to the global minimizer of the objective function to the boundary point 
  virtual vector<double>  SetBoundaryShift();

public:
  IConstrainedOptProblem(int dim, vector<double> loBound, vector<double> upBound,
    int probIndex = -1, EConstrainedProblemType problemType = cptNormal,
    double fraction = 1, int activeConstrNum = 0);

  /// Get global minimizer 
  vector<double> GetOptimumPoint() const;
  /// Get global minimum value
  double GetOptimumValue() const;

  /// What is specified for the objective function
  bool GetStatus(enum EOptFunctionParameter param) const;
  /// Get global maximizer 
  vector<double> GetMaxPoint() const;
  /// Get global maximum value
  double GetMaxValue() const;
  /// Get a Lipschitz constant for the objective function
  double GetLipschitzConstant() const;

  /// Compute the objective function value at the point y
  double ComputeFunction(const vector<double>& y) const;
  /// Compute the objective function derivatives value at the point y
  vector<double> ComputeFunctionDerivatives(const vector<double>& y) const;

  /// Get the constraints number
  int GetConstraintsNumber() const;

  /// What is specified for the constraint number index
  bool GetConstraintStatus(int index, enum EOptFunctionParameter param) const;
  /// Return a Lipschitz constant for the constraint number index
  double GetConstraintLipschitzConstant(int index) const;

  /// Compute the value of the constraint number index at the point y
  double ComputeConstraint(int index, const vector<double>& y) const;

  /// Compute the value of the constraints at the point y, index is the number of the last computed constraint
  vector<double> ComputeConstraints(const vector<double>& y,
    EConstraintComputationType t, int &index) const;

  /// Compute the derivatives of the constraint number index at the point y
  vector<double> ComputeConstraintDerivatives(int index, const vector<double>& y) const;

  virtual ~IConstrainedOptProblem();
};