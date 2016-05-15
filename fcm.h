#ifndef __FCM_H__
#define __FCM_H__
#include "shared.h"

static inline void compute_prototypes() {
    vector<double> t[num_docs];
    for (int i = 0; i < num_docs; i++) {
        for (int j = 0; j < num_clusters; j++) {
            t[i].pb(pow(memberships[i][j], fuzziness));
        }
    }
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < num_clusters; j++) {
        for (int k = 0; k < num_terms; k++) {
            double numerator = 0.0;
            double denominator = 0.0;
            for (int i = 0; i < num_docs; i++) {
                numerator += t[i][j] * docs[i][k];
                denominator += t[i][j];
            }
            prototypes[j][k] = numerator / denominator;
        }
    }
}

static inline double update_memberships() {
    double max_diff = 0.0;
    for (int j = 0; j < num_clusters; j++) {
      for (int i = 0; i < num_docs; i++) {
        double new_uij = get_new_value(i, j);
        double diff = new_uij - memberships[i][j];
        memberships[i][j] = new_uij;

        max_diff = MAX(diff, max_diff);
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
