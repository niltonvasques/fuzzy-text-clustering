#ifndef __SHARED_H__
#define __SHARED_H__

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

#ifndef __MAIN_FILE__

using namespace std;

vector<string> terms;
vector<double> docs[MAX_DOCS];
vector<double> prototypes[MAX_CLUSTERS];
vector<int>    clusters[MAX_CLUSTERS];
vector<double> memberships[MAX_DOCS];
vector<double> tipicalities[MAX_DOCS];
vector<double> gamas(MAX_CLUSTERS, 0);
vector<double> final_memberships[MAX_DOCS];
vector<double> final_tipicalities[MAX_DOCS];
vector<pii > crisp(MAX_DOCS, mp(0,1));
int num_terms;
int num_docs;
int num_clusters = 3;
double fuzziness = 1.2; 
double epsilon = 0.01;
double a = 1;
double b = 2;
double fuzziness_n = 1.2; 
double fuzziness_m = 1.2; 

#else

extern std::vector<string> terms;
extern std::vector<double> docs[MAX_DOCS];
extern std::vector<double> prototypes[MAX_CLUSTERS];
extern std::vector<double> tipicalities[MAX_DOCS];
extern std::vector<int>  clusters[MAX_CLUSTERS];
extern std::vector<double> memberships[MAX_DOCS];
extern std::vector<double> final_memberships[MAX_DOCS];
extern vector<double> gamas;
extern std::vector<pii> crisp;
extern int num_terms;
extern int num_docs;
extern int num_clusters;
extern double fuzziness; 
extern double epsilon;
extern double a;
extern double b;
extern double fuzziness_n = 1.2; 
extern double fuzziness_m = 1.4; 

#endif

static inline void read_data(){

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
    if(fname == ""){
      f = stdout;
    }else if ((f = fopen(fname.c_str(), "w")) == NULL) {
        printf("Cannot create output file.\n");
        exit(1);
    }
    for (i = 0; i < size; i++) {
        for (j = 0; j < matrix[i].size(); j++) {
            fprintf(f, "%lf\t", matrix[i][j]);
        }
        fprintf(f, "\n");
    }
    if(fname != "")
      fclose(f);
}

static inline double norm_doc2doc(int i, int j) {
    int k;
    double sum = 0.0;
    for (k = 0; k < num_terms; k++) {
        sum += pow(docs[i][k] - docs[j][k], 2);
    }
    return sqrt(sum);
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

static inline void debug_memberships(){
  times(i, num_docs){
    times(j, num_clusters){
      printf("%.8f\t", memberships[i][j]);
    }
    printf("\n");
  }
}

static inline void generate_memberships(){
  double s,rval;
  int r;
  crisp.clear();
  for (int i = 0; i < num_docs; i++) {
    s = 0.0;
    r = 100000;
    crisp.pb(mp(0,1));
    memberships[i].clear();
    for (int j = 0; j < num_clusters-1; j++) {
      rval = rand() % (r + 1);
      r -= rval;
      memberships[i].pb(rval / 100000.0);
      s += memberships[i][j];
    }
    memberships[i].pb(1.0 - s);
  }
}

static inline void init_prototypes(){
  times(i, num_clusters){
    prototypes[i].clear();
    times(j, num_terms){
      prototypes[i].pb(0);
    }
  }
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
      if(h != g && alphas_count[h] > 0){
        beta = MIN(beta, alphas[h]/alphas_count[h]);
      }
    }
    s = (beta - alpha_g)/MAX(beta, alpha_g);
    double u1 = memberships[i][crisp[i].first];
    double u2 = memberships[i][crisp[i].second];
    sum_up += s * (u1 - u2);
    sum_down += (u1 - u2);
    printf("d%d: %lf beta %lf alpha %lf u1 %lf u2 %lf sumup %lf sumdown\n", i, beta,
        alpha_g, u1, u2, sum_up, sum_down);
  }
  fs = sum_up / sum_down;
  return fs;
}

static inline void store_final_memberships(){
  times(i, num_docs){
    final_memberships[i].swap(memberships[i]);
  }
}
static inline void store_final_tipicalities(){
  times(i, num_docs){
    final_tipicalities[i].swap(tipicalities[i]);
  }
}

#endif
