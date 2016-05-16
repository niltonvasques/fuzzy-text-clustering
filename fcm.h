#ifndef __FCM_H__
#define __FCM_H__
#include "shared.h"

static inline void compute_prototypes() {
  double t[num_docs][num_clusters];
#pragma omp parallel for
  for (uint i = 0; i < num_docs; i++) {
    for (uint j = 0; j < num_clusters; j++) {
      double tij = pow(fabs(memberships[i][j]), fuzziness);
      t[i][j] = tij;
    }
  }
#pragma omp parallel for collapse(2)
  for (uint j = 0; j < num_clusters; j++) {
    for (uint k = 0; k < num_terms; k++) {
      double numerator = 0.0;
      double denominator = 0.0;
      for (uint i = 0; i < num_docs; i++) {
        numerator += t[i][j] * docs[i][k];
        denominator += t[i][j];
      }
      prototypes[j][k] = numerator / denominator;
    }
  }
}

static inline double update_memberships() {
  double max_diff = 0.0;
  double diffs[num_docs][num_clusters];
#pragma omp parallel for
  for (uint j = 0; j < num_clusters; j++) {
    for (uint i = 0; i < num_docs; i++) {
      double new_uij = get_new_value(i, j);
      double diff = new_uij - memberships[i][j];
      memberships[i][j] = new_uij;
      diffs[i][j] = diff;

    }
  }
#pragma omp parallel for reduction(max:max_diff)
  for (uint j = 0; j < num_clusters; j++) {
    for (uint i = 0; i < num_docs; i++) {
      max_diff = MAX(diffs[i][j], max_diff);
    }
  }
  return max_diff;
}

static inline void fcm(){
  generate_memberships();
  init_prototypes();

  save_matrix("inital_memberships.matrix", memberships, num_docs);
  save_matrix("inital_prototypes.matrix", prototypes, num_clusters);

  double max_diff;
  do {
    compute_prototypes();
    max_diff = update_memberships();
    if(arguments.verbose) 
      cout << max_diff << endl;
  } while (max_diff > epsilon);
}


#endif
