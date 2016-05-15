#ifndef __DESCRIPTORS_H__
#define __DESCRIPTORS_H__
#include "shared.h"

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

#endif
