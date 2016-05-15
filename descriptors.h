#ifndef __DESCRIPTORS_H__
#define __DESCRIPTORS_H__
#include "shared.h"

struct SCORE {
  double score;
  uint term;
  uint cluster;
  SCORE() {} 
  SCORE(double _score, uint _term, uint _cluster) : score(_score), term(_term), cluster(_cluster) {} 
};

class SCOREComparisson
{
  bool reverse;
  public:
  SCOREComparisson(const bool& revparam=false)
  {reverse=revparam;}
  bool operator() (const SCORE& lhs, const SCORE&rhs) const
  {
    if (reverse) return (lhs.score>rhs.score);
    else return (lhs.score<rhs.score);
  }
};
// using mycomparison:
typedef std::priority_queue<SCORE,std::vector<SCORE>,SCOREComparisson> pq_score;

#define pdpii pair<double, pair<uint,uint> >


static inline void save_ranking(string fname, vector<SCORE> &ranking){
  FILE *f;
  if ((f = fopen(fname.c_str(), "w")) == NULL) {
    printf("Cannot create output file.\n");
    exit(1);
  }

  for (uint i = 0; i < ranking.size(); i++) {
    SCORE item = ranking[i];
    fprintf(f, "%s\t%lf\t", terms[item.term].c_str(), item.score);
    fprintf(f, "\n");
  }
  fclose(f);
}

static inline void soft_fdcl(){
  vector<bool> candidates(num_terms, true);

  double threshold = 1.0 / num_clusters;
  pq_score ranking;

  times(i, num_clusters){
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
      ranking.push(SCORE(f1,j,i)); 
    }
  }

  vector<SCORE> descriptors[num_clusters];
  vector<SCORE> ranking_by_cluster[num_clusters];
  vector<bool> choosed(num_terms, false);
  uint selected = 0;
  while(selected < num_descriptors*num_clusters && !ranking.empty()){
    SCORE item = ranking.top(); ranking.pop();
    ranking_by_cluster[item.cluster].pb(item);
    if(!choosed[item.term] && descriptors[item.cluster].size() < num_descriptors){
      choosed[item.term] = true;
      descriptors[item.cluster].pb(item);
      selected++;
    }
  }
  times(i, num_clusters){
    ostringstream oss, oss2;
    oss << "descriptors.cluster" << i;
    save_ranking(oss.str(), descriptors[i]);
    oss2 << "ranking.cluster" << i;
    save_ranking(oss2.str(), ranking_by_cluster[i]);
  }
}

#endif