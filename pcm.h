#ifndef __PCM_H__
#define __PCM_H__
#include "shared.h"

static inline void pcm_compute_prototypes() {
  double t[num_docs][num_clusters];
#pragma omp parallel for
  for (uint i = 0; i < num_docs; i++) {
    for (uint j = 0; j < num_clusters; j++) {
      // Its required to use fabs because pow of negative numbers sometimes
      // can result in NaN
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

static inline void estimate_gamas() {
  double k = 1;

#pragma omp parallel for
  for (uint j = 0; j < num_clusters; j++) {
    double numerator = 0;
    double denominator = 0;
    for (uint i = 0; i < num_docs; i++) {
      double dij = get_norm(i, j, docs, prototypes); 
      double uij = memberships[i][j]; 
      uij = pow(fabs(uij), fuzziness);
    //  printf("dij %f uij %f\n", dij, uij);
      numerator += uij * dij;
      denominator += uij;
    }
    gamas[j] = k * (numerator / denominator);
    if(isnan(gamas[j])) exit(1);
    if(arguments.verbose) 
      printf("gamas[%d]: %f\n", j, gamas[j]);
  }
}

static inline double tipicality(double distance, uint j){
  double denominator = pow(distance,2) / gamas[j];
  double exp = 1.0 / (fuzziness - 1.0);
  denominator = 1.0 + pow(denominator, exp);
  if(isnan(denominator)) return 0;
  return 1.0 / denominator;
}

static inline double update_tipicalities() {
  double sum_kn = 0;
  double sum_lc = 0;
  double sum_jn[num_clusters];

#pragma omp parallel for
  for (uint j = 0; j < num_clusters; j++) {
    sum_jn[j] = 0;
  }

#pragma omp parallel for reduction(+:sum_kn)
  for (uint j = 0; j < num_clusters; j++) {
    double sum_ic = 0;
    double sumjn = 0;
    for (uint i = 0; i < num_docs; i++) {
      double distance;
      double new_uij;
      double tik;
      distance = get_norm(i, j, docs, prototypes);
      new_uij = tipicality(distance, j);
      tik = pow(fabs(new_uij), fuzziness);
      sum_ic += tik * distance;
      tik = pow(fabs(1 - new_uij), fuzziness);
      sumjn += tik;
      memberships[i][j] = new_uij;
    }
    sum_kn += sum_ic;
    sum_jn[j] = sumjn;
  }

#pragma omp parallel for reduction(+:sum_lc)
  for (uint j = 0; j < num_clusters; j++) {
    sum_lc += sum_jn[j] * gamas[j];
  }

  return sum_kn + sum_lc;
}

static inline void default_prototypes(){
  prototypes[0].pb(-3.34);
  prototypes[0].pb(1.67);
  prototypes[1].pb(1.67);
  prototypes[1].pb(0.00);
}

static inline void init_gamas(){
  gamas.clear();
  times(i, num_clusters){
    gamas.pb(0);
  }
}

int pcm() {
  double max_diff;
  double curr_j = 0, old_j = 0;

  init_gamas();
  
  if(arguments.random){
    //default_prototypes();
    generate_memberships();
    init_prototypes();
    pcm_compute_prototypes();
  }else{
    fcm();
  }	
  save_matrix("initial_memberships.matrix", memberships, num_docs);
  save_matrix("initial_prototypes.matrix", prototypes, num_clusters);
  estimate_gamas();
  old_j = update_tipicalities();
  do {
    pcm_compute_prototypes();
    curr_j = update_tipicalities();
    max_diff = abs(curr_j - old_j);
    old_j = curr_j;
    if(arguments.verbose) 
      printf("max_diff: %f\n", max_diff);
  } while (max_diff > epsilon);
  return 0;
}


#endif
