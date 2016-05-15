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

        double degree = final_memberships[k][i];
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
}

static inline double aswc(){
  double fs = 0;
  double s = 0;
  times(i, num_docs){
    double max1_degree = 0;
    double max2_degree = 0;
    times(j, num_clusters){
      if(j == 0 || memberships[i][j] > max1_degree){
        max1_degree = memberships[i][j];
        pii item = crisp[i];
        crisp[i] = mp(j,item.first);
        max2_degree = memberships[i][item.first];
      }else if(memberships[i][j] > max2_degree){
        max2_degree = memberships[i][j];
        crisp[i].second = j;
      }
    }
  }
  double sum_up = 0, sum_down = 0;
  times(i, num_docs){
    int g = crisp[i].first;
    vector<double> alphas(num_clusters, 0);
    vector<int> alphas_count(num_clusters, 0);
    times(j, num_docs){
      int grupo = crisp[j].first;
      alphas[grupo] += norm_doc2doc(i, j);
      alphas_count[grupo]++;
    }
    double alpha_g = alphas[g] / alphas_count[g];
    double beta = 999;
    times(h, num_clusters){
      if(h != g){
        beta = MIN(beta, alphas[h]/alphas_count[h]);
      }
    }
    s = (beta - alpha_g)/MAX(beta, alpha_g);
    double u1 = memberships[i][crisp[i].first];
    double u2 = memberships[i][crisp[i].second];
    sum_up += s * (u1 - u2);
    sum_down += (u1 - u2);
    //printf("d%d: %lf beta %lf alpha %lf u1 %lf u2 %lf sumup %lf sumdown\n", i, beta,
     //   alpha_g, u1, u2, sum_up, sum_down);
  }
  fs = sum_up / sum_down;
  return fs;
}

int main(){

  read_data();

  double max_fs = -1;
  double fs;
  int max_groups = 2;
  for(int i = 3; i < 3; i++){
    num_clusters = i;
    fcm();
    fs = aswc();
    cout << "FS: " << fs << endl;
    printf("%lf aswc %d groups\n", fs, i);
    if(fs > max_fs){
      max_fs = fs;
      max_groups = i;
      copy_memberships();
      save_matrix("", final_memberships, num_docs);
      printf("max fs found: %lf aswc %d groups\n", fs, i);
    }
  }
  printf("final fs: %lf aswc %d groups\n", max_fs, max_groups);
  save_matrix("memberships.matrix", final_memberships, num_docs);
  save_matrix("prototypes.matrix", prototypes, num_clusters);
  soft_fdcl();

  return 0;
}
