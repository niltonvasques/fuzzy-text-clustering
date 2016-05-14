#include <bits/stdc++.h>
#include "shared.h"


using namespace std;


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

static inline void save_ranking(string fname, priority_queue<pdi> &ranking){
  FILE *f;
  if ((f = fopen(fname.c_str(), "w")) == NULL) {
    printf("Cannot create output file.\n");
    exit(1);
  }

  for (uint i = 0; i < ranking.size(); i++) {
    pdi item = ranking.top();
    ranking.pop();
    fprintf(f, "%s\t%lf\t", terms[item.second].c_str(), item.first);
    fprintf(f, "\n");
  }
  fclose(f);
}

static inline void soft_fdcl(){
  vector<bool> candidates(num_terms, true);

  double threshold = 1.0 / num_clusters;

  times(i, num_clusters){
    priority_queue<pdi> ranking;
    times(j, num_terms){
      double relevants = 0;
      double recovereds = 0;
      double correct_recovereds = 0;
      double precison = 0;
      double recall = 0;
      double f1 = 0;
      times(k, num_docs){
        if(docs[k][j] > 0.0) recovereds++;

        double degree = memberships[k][i];
        if(degree >= threshold){
          relevants++;
          if(docs[k][j] > 0.0) correct_recovereds++;
        }
      }
      if(recovereds == 0) precison = 1;
      else precison = correct_recovereds / recovereds;

      if(relevants == 0 ) recall = 1;
      else recall = correct_recovereds / relevants;
      
      if(precison == 0 && recall == 0) f1 = 0;
      else f1 = (2 * precison * recall)/(precison + recall);

      if(isnan(recovereds))         printf("%lf recovered is nan", recovereds);
      if(isnan(correct_recovereds)) printf("%lf correct_recovereds is nan", correct_recovereds);
      if(isnan(relevants))          printf("%lf relevants is nan", relevants);
      if(isnan(precison))           printf("%lf precison is nan", precison);
      if(isnan(recall))             printf("%lf recall is nan", recall);
      if(isnan(f1)){
        printf("  %lf recovered is nan", recovereds);
        printf("  %lf correct_recovereds is nan", correct_recovereds);
        printf("  %lf relevants is nan", relevants);
        printf("  %lf precison is nan", precison);
        printf("  %lf recall is nan", recall);
        printf("  %lf f1 is nan\n", f1);
      } 
      ranking.push(mp(f1,j)); 
    }
    ostringstream oss;
    oss << "ranking.cluster" << i;
    save_ranking(oss.str(), ranking);
  }
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

  save_matrix("memberships.matrix", memberships, num_docs);
  save_matrix("prototypes.matrix", prototypes, num_clusters);
}

int main(){

  read_data();

  fcm();
  soft_fdcl();

  return 0;
}
