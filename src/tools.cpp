#include "tools.h"

using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
   int num_measurements = estimations.size();
   int num_variables = estimations[0].size();
   VectorXd residuals = VectorXd::Zero(num_variables);

   for (int i = 0; i < num_variables; i++) {
     for (int j = 0; j < num_measurements; j++) {
       residuals[i] += pow(estimations[j][i] - ground_truth[j][i], 2);
     }

     residuals[i] = sqrt(residuals[i] / num_measurements);
   }

   return residuals;
}
