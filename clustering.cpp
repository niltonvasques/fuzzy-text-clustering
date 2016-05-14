#include <bits/stdc++.h>

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define ll long long
#define pb push_back
#define mp make_pair
#define pii pair<int,int>
#define pdd pair<double,double>
#define pid pair<int,double>
#define pdi pair<double,int>
#define pld pair<ll,double>
#define pll pair<ll,ll>
#define DEBUG true
#define DEBUG_LOG(A) if(DEBUG) cout << (A) << endl
#define repeat(SIZE) for(int index = 0; index < (SIZE); index++)
#define times(VAR,SIZE) for(int VAR = 0; VAR < (SIZE); VAR++)
#define MAX_DOCS 20000
#define MAX_CLUSTERS 100 
#define uint unsigned int
#define ull unsigned long long 

using namespace std;

vector<string> terms;
vector<double> docs[MAX_DOCS];
vector<double> prototypes[MAX_CLUSTERS];
vector<double> memberships[MAX_DOCS];
int num_terms;
int num_docs;
int num_clusters = 3;
double fuzziness = 1.2; 
double epsilon = 0.001;

void read_data(){

  string line;

  cin >> num_terms >> num_docs;

  cin >> line; // First discover.names line is useless

  repeat(num_terms-1){
    cin >> line;
    int token_index = line.find(":");
    line.erase(line.begin()+token_index, line.end());
    line = line.substr(1,line.size()-2);
    terms.pb(line);
    
    //DEBUG_LOG(line);
  }

  times(i, num_docs){
    getline(cin, line, ','); // first value is just the filename
    //DEBUG_LOG("TITLE");
    //DEBUG_LOG(line);
    double frequency; 
    times(j, num_terms-1){
      if(j == (num_terms-2)){ //Last field don't have a comma
        cin >> frequency; 
      }else{
        getline(cin, line, ',');
        frequency = atof(line.c_str());
      }
      docs[i].pb(frequency);
    }
    //DEBUG_LOG("endline");
    //DEBUG_LOG(line);
    //DEBUG_LOG(frequency);
  }
}

static inline void save_matrix(string fname, vector<double> *matrix, uint size) {
    uint i, j;
    FILE *f;
    if ((f = fopen(fname.c_str(), "w")) == NULL) {
        printf("Cannot create output file.\n");
        exit(1);
    }
    for (i = 0; i < size; i++) {
        for (j = 0; j < matrix[i].size(); j++) {
            fprintf(f, "%lf\t", matrix[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

static inline double get_norm(int i, int j) {
    int k;
    double sum = 0.0;
    for (k = 0; k < num_terms; k++) {
        sum += pow(docs[i][k] - prototypes[j][k], 2);
    }
    return sqrt(sum);
}

static inline double get_new_value(int i, int j) {
    int k;
    double t, p, sum;
    sum = 0.0;
    p = 2 / (fuzziness - 1);
    for (k = 0; k < num_clusters; k++) {
        t = get_norm(i, j) / get_norm(i, k);
        t = pow(t, p);
        sum += t;
    }
    return 1.0 / sum;
}

void debug_memberships(){
  times(i, num_docs){
    times(j, num_clusters){
      printf("%.8f\t", memberships[i][j]);
    }
    printf("\n");
  }
}

void generate_memberships(){
  double s,rval;
  int r;
  for (int i = 0; i < num_docs; i++) {
    s = 0.0;
    r = 100000;
    for (int j = 0; j < num_clusters-1; j++) {
      rval = rand() % (r + 1);
      r -= rval;
      memberships[i].pb(rval / 100000.0);
      s += memberships[i][j];
    }
    memberships[i].pb(1.0 - s);
  }
}

void init_prototypes(){
  times(i, num_clusters){
    times(j, num_terms){
      prototypes[i].pb(0);
    }
  }
}

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

void save_ranking(string fname, priority_queue<pdi> &ranking){
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

void soft_fdcl(){
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

void fcm(){
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
