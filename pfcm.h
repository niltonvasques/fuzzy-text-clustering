#ifndef __PFCM_H__
#define __PFCM_H__
#include "shared.h"
#include "pcm.h"
#include "fcm.h"

static inline void pfcm_compute_prototypes() {
  int i, j, k;
  double numerator, denominator;
  vector<double> t[num_docs];
  for (i = 0; i < num_docs; i++) {
    for (j = 0; j < num_clusters; j++) {
      double result = a * pow(memberships[i][j], fuzziness_m) +
        b * pow(tipicalities[i][j], fuzziness_n);
      t[i].pb(result);
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

static inline double pfcm_tipicality(double distance, int j){
  double denominator = distance / gamas[j];
  double exp = 1.0 / (fuzziness_n - 1.0);
  denominator = 1.0 + pow(denominator, exp);
  return 1.0 / denominator;
}

static inline double pfcm_update_memberships_and_tipicalities() {
  int i, j;
  double new_uij;
  double new_tik;
  double sum_kn = 0;
  double sum_lc = 0;
  double sum_jn[num_clusters];
  double distance;

  for (j = 0; j < num_clusters; j++) {
    sum_jn[j] = 0;
  }

  for (i = 0; i < num_docs; i++) {
    for (j = 0; j < num_clusters; j++) { 
      distance = get_norm(i, j); 
      new_uij = get_new_value(i, j);
      tipicalities[i][j] = pfcm_tipicality(distance, j);

      //degree_of_memb[i][j] = a * new_uij + b * new_tik;
      memberships[i][j] = new_uij;

      new_uij = pow(new_uij, fuzziness_m);

      sum_kn += (a * new_uij + b * pow(tipicalities[i][j], fuzziness_n)) * ( distance * distance );
      new_tik = pow(1 -tipicalities[i][j], fuzziness_n);
      sum_jn[j] += new_tik;
    }
  }

  for (j = 0; j < num_clusters; j++) {
    sum_lc += sum_jn[j] * gamas[j];
  }

  return sum_kn + sum_lc;
}

void init_tipicalities(){
  times(i, num_docs){
    tipicalities[i].clear();
    times(j, num_clusters){
      //double distance = get_norm(i, j); 
      //tipicalities[i].pb(pfcm_tipicality(distance, j));
      tipicalities[i].pb(0);
    }
  }
}

int pfcm() {
  double max_diff;
  double curr_j = 0, old_j = 0;

  fcm();
  gamas.clear();
  times(i, num_clusters){
    gamas.pb(0);
  }
  estimate_gamas();
  do {
    init_tipicalities();
    pfcm_compute_prototypes();
    curr_j = pfcm_update_memberships_and_tipicalities();
    max_diff = fabs(curr_j - old_j);
    old_j = curr_j;
    printf("max_diff: %f curr_j: %f\n", max_diff, curr_j);
  } while (max_diff > epsilon);
  return 0;
}


#endif
