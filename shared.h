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

std::vector<string> terms;
std::vector<double> docs[MAX_DOCS];
std::vector<double> prototypes[MAX_CLUSTERS];
std::vector<double> memberships[MAX_DOCS];
int num_terms;
int num_docs;
int num_clusters = 3;
double fuzziness = 1.2; 
double epsilon = 0.01;

#else

extern std::vector<string> terms;
extern std::vector<double> docs[MAX_DOCS];
extern std::vector<double> prototypes[MAX_CLUSTERS];
extern std::vector<double> memberships[MAX_DOCS];
extern int num_terms;
extern int num_docs;
extern int num_clusters;
extern double fuzziness; 
extern double epsilon;

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

static inline void init_prototypes(){
  times(i, num_clusters){
    times(j, num_terms){
      prototypes[i].pb(0);
    }
  }
}

#endif
