#ifndef __FCM_H__
#define __FCM_H__
#include "shared.h"

static inline void compute_prototypes() {
    int i, j, k;
    double numerator, denominator;
    vector<double> t[num_docs];
    for (i = 0; i < num_docs; i++) {
        for (j = 0; j < num_clusters; j++) {
            t[i].pb(pow(memberships[i][j], fuzziness));
        }
    }
    for (j = 0; j < num_clusters; j++) {
        for (k = 0; k < num_terms; k++) {
            numerator = 0.0;
            denominator = 0.0;
            for (i = 0; i < num_docs; i++) {
                numerator += t[i][j] * docs[i][k];
                denominator += t[i][j];
            }
            prototypes[j][k] = numerator / denominator;
        }
    }
}

static inline double update_memberships() {
    int i, j;
    double new_uij;
    double max_diff = 0.0, diff;
    for (j = 0; j < num_clusters; j++) {
      for (i = 0; i < num_docs; i++) {
        new_uij = get_new_value(i, j);
        diff = new_uij - memberships[i][j];
        if (diff > max_diff)
          max_diff = diff;
        memberships[i][j] = new_uij;
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
    cout << max_diff << endl;
  } while (max_diff > epsilon);
}


#endif
