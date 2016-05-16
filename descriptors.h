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

static inline void save_arff(string fname, vector<SCORE> *descriptors){
  FILE *f;
  if ((f = fopen(fname.c_str(), "w")) == NULL) {
    printf("Cannot create output file.\n");
    exit(1);
  }
  fprintf(f, "@RELATION relation_name\n\n");
  times(i, max_groups){
    times(j, descriptors[i].size()){
      SCORE item = descriptors[i][j];
      fprintf(f, "@ATTRIBUTE %s\tNUMERIC\n",terms[item.term].c_str());
    }
  }
  fprintf(f, "@ATTRIBUTE class \t{");
  times(i, max_groups){
    if(i > 0)
      fprintf(f, ",");
    fprintf(f, "cluster_%d", i);
  }
  fprintf(f, "}\n\n");
  fprintf(f, "@DATA\n");
  times(k, num_docs){
    double max_degree = 0;
    uint g = 0;
    times(i, max_groups){
      if(memberships[k][i] > max_degree){
        g = i;
        max_degree = memberships[k][i];
      }
      times(j, descriptors[i].size()){
        SCORE item = descriptors[i][j];
        fprintf(f, "%lf,",docs[k][item.term]);
      }
    }
    fprintf(f, "cluster_%d\n",g);
  }
  fclose(f);
}

static inline void save_crisp(string fname, vector<SCORE> *descriptors){
  FILE *f;
  if ((f = fopen(fname.c_str(), "w")) == NULL) {
    printf("Cannot create output file.\n");
    exit(1);
  }
  times(k, num_docs){
    double max_degree = 0;
    uint g = 0;
    times(i, max_groups){
      if(memberships[k][i] > max_degree){
        g = i;
        max_degree = memberships[k][i];
      }
    }
    fprintf(f, "%d\n",g);
  }
  fclose(f);
}

static inline void soft_fdcl(){
  vector<bool> candidates(num_terms, true);

  double threshold = 1.0 / max_groups;
  pq_score ranking;

  times(i, max_groups){
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

  vector<SCORE> descriptors[max_groups];
  vector<SCORE> ranking_by_cluster[max_groups];
  vector<bool> choosed(num_terms, false);
  uint selected = 0;
  while(!ranking.empty()){
    SCORE item = ranking.top(); ranking.pop();
    ranking_by_cluster[item.cluster].pb(item);

    if(selected < num_descriptors*max_groups && !choosed[item.term] &&
        descriptors[item.cluster].size() < num_descriptors){
      choosed[item.term] = true;
      descriptors[item.cluster].pb(item);
      selected++;
    }
  }
  times(i, max_groups){
    ostringstream oss, oss2;
    oss << arguments.path << "descriptors.soft-fdcl.cluster" << i;
    save_ranking(oss.str(), descriptors[i]);
    oss2 << arguments.path << "ranking.soft-fdcl.cluster" << i;
    save_ranking(oss2.str(), ranking_by_cluster[i]);
  }
  ostringstream oss, oss2;
  oss << arguments.path << "clusters.soft-fdcl.arff";
  save_arff(oss.str(), descriptors);

  oss2 << arguments.path << "soft-fdcl.clusters";
  save_crisp(oss2.str(), descriptors);
}


static inline void pdcl(){
  double threshold = 1.0 / max_groups;
  pq_score ranking;

  vector<double> normalized[num_docs];
  times(i, num_docs){
    double sum = 0;
    times(j, max_groups){
      sum += final_memberships[i][j];
    }
    times(j, max_groups){
      normalized[i].pb(final_memberships[i][j]/sum); 
    }
  }

  times(i, max_groups){
    times(j, num_terms){
      double relevants = 0;
      double recovereds = 0;
      double correct_recovereds = 0;
      double precison = 0;
      double recall = 0;
      double f1 = 0;
      times(k, num_docs){

        double degree = final_memberships[k][i];
        double degree_norm = normalized[k][i];

        if(docs[k][j] > 0.0) recovereds += degree;
        if(degree_norm >= threshold){
          relevants += degree;
          if(docs[k][j] > 0.0) correct_recovereds += degree;
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

  vector<SCORE> descriptors[max_groups];
  vector<SCORE> ranking_by_cluster[max_groups];
  vector<bool> choosed(num_terms, false);
  uint selected = 0;
  while(!ranking.empty()){
    SCORE item = ranking.top(); ranking.pop();
    ranking_by_cluster[item.cluster].pb(item);

    if(selected < num_descriptors*max_groups && !choosed[item.term] &&
        descriptors[item.cluster].size() < num_descriptors){
      choosed[item.term] = true;
      descriptors[item.cluster].pb(item);
      selected++;
    }
  }
  times(i, max_groups){
    ostringstream oss, oss2;
    oss << arguments.path << "descriptors.pdcl.cluster" << i;
    save_ranking(oss.str(), descriptors[i]);
    oss2 << arguments.path << "ranking.pdcl.cluster" << i;
    save_ranking(oss2.str(), ranking_by_cluster[i]);
  }
  ostringstream oss,oss2;
  oss << arguments.path << "clusters.pdcl.arff";
  save_arff(oss.str(), descriptors);
  oss2 << arguments.path << "pdcl.clusters";
  save_crisp(oss2.str(), descriptors);
}

static inline void mixed_pdcl(){
  double threshold = 1.0 / max_groups;
  pq_score ranking;

  vector<double> merged[num_docs];
  times(i, num_docs){
    double sum = 0;
    times(j, max_groups){
      sum += final_tipicalities[i][j];
    }
    times(j, max_groups){
      double t_degree = final_tipicalities[i][j]/sum;
      double m_degree = final_memberships[i][j]/sum;
      merged[i].pb(( a * m_degree + b * t_degree )/(a+b)); 
    }
  }

  times(i, max_groups){
    times(j, num_terms){
      double relevants = 0;
      double recovereds = 0;
      double correct_recovereds = 0;
      double precison = 0;
      double recall = 0;
      double f1 = 0;
      times(k, num_docs){

        double degree = final_tipicalities[k][i];
        double merged_degree = merged[k][i];

        if(docs[k][j] > 0.0) recovereds += degree;
        if(merged_degree >= threshold){
          relevants += degree;
          if(docs[k][j] > 0.0) correct_recovereds += degree;
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

  vector<SCORE> descriptors[max_groups];
  vector<SCORE> ranking_by_cluster[max_groups];
  vector<bool> choosed(num_terms, false);
  uint selected = 0;
  while(!ranking.empty()){
    SCORE item = ranking.top(); ranking.pop();
    ranking_by_cluster[item.cluster].pb(item);

    if(selected < num_descriptors*max_groups && !choosed[item.term] &&
        descriptors[item.cluster].size() < num_descriptors){
      choosed[item.term] = true;
      descriptors[item.cluster].pb(item);
      selected++;
    }
  }
  times(i, max_groups){
    ostringstream oss, oss2;
    oss << arguments.path << "descriptors.mixed-pdcl.cluster" << i;
    save_ranking(oss.str(), descriptors[i]);
    oss2 << arguments.path << "ranking.mixed-pdcl.cluster" << i;
    save_ranking(oss2.str(), ranking_by_cluster[i]);
  }
  ostringstream oss,oss2;
  oss << arguments.path << "clusters.mixed-pdcl.arff";
  save_arff(oss.str(), descriptors);
  oss2 << arguments.path << "mixed-pdcl.clusters";
  save_crisp(oss2.str(), descriptors);
}

#endif
