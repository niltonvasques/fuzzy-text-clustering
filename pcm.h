#ifndef __PCM_H__
#define __PCM_H__
#include "shared.h"

static inline void pcm_compute_prototypes() {
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

static inline void estimate_gamas() {
  double k = 1;
  int i, j;

  for (j = 0; j < num_clusters; j++) {
    double numerator = 0;
    double denominator = 0;
    for (i = 0; i < num_docs; i++) {
      double dij = get_norm(i, j); 
      double uij = memberships[i][j]; 
      uij = pow(uij, fuzziness);
      //printf("dij %f uij %f\n", dij, uij);
      numerator += uij * dij;
      denominator += uij;
    }
    //printf("numerator %f denominator %f\n", numerator, denominator);
    gamas[j] = k * (numerator / denominator);
    printf("gamas[%d]: %f\n", j, gamas[j]);
  }
}

static inline double tipicality(double distance, int j){
  double denominator = distance / gamas[j];
  double exp = 1.0 / (fuzziness - 1.0);
  denominator = 1.0 + pow(denominator, exp);
  return 1.0 / denominator;
}

static inline double update_tipicalities() {
  int i, j;
  double new_uij;
  double tik;
  double sum_ic = 0;
  double sum_kn = 0;
  double sum_lc = 0;
  double sum_jn[num_clusters];
  double distance;

  for (j = 0; j < num_clusters; j++) {
    sum_jn[j] = 0;
  }

  for (i = 0; i < num_docs; i++) {
    sum_ic = 0;
    for (j = 0; j < num_clusters; j++) {
      distance = get_norm(i, j);
      new_uij = tipicality(distance, j);
      tik = pow(new_uij, fuzziness);
      sum_ic += tik * distance;
      tik = pow(1 - new_uij, fuzziness);
      sum_jn[j] += tik;
      memberships[i][j] = new_uij;
    }
    sum_kn += sum_ic;
  }

  for (j = 0; j < num_clusters; j++) {
    sum_lc += sum_jn[j] * gamas[j];
  }

  return sum_kn + sum_lc;
}

int pcm() {
  double max_diff;
  double curr_j = 0, old_j = 0;

  fcm();
  gamas.clear();
  times(i, num_clusters){
    gamas.pb(0);
  }
  estimate_gamas();
  do {
    pcm_compute_prototypes();
    curr_j = update_tipicalities();
    max_diff = abs(curr_j - old_j);
    old_j = curr_j;
    printf("max_diff: %f\n", max_diff);
  } while (max_diff > epsilon);
  return 0;
}


#endif
