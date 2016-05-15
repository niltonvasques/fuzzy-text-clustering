#include <bits/stdc++.h>
#include "shared.h"
#include "fcm.h"
#include "pcm.h"
#include "pfcm.h"
#include "descriptors.h"

using namespace std;

#include <argp.h>
#include <stdbool.h>

const char *argp_program_version = "Fuzzy Clustering";
const char *argp_program_bug_address = "<niltonvasques@dcc.ufba.br>";
static char doc[] = "Implementing some fuzzy clustering algorithms for flexible document organization.";
static char args_doc[] = "[FILENAME]...";
static struct argp_option options[] = { 
  { "fcm", 'f', 0, 0, "Run fuzzy c means."},
  { "pcm", 'p', 0, 0, "Run possibilistc c means."},
  { "pfcm", 'x', 0, 0, "Run possibilistc fuzzy c means."},
  { "input", 'i', "INPUT", 0, "Path for input data."},
  { "verbose", 'v', 0, 0, "More verbose execution."},
  { "a", 'a', "VALUE", 0, "PFCM a param. DEFAULT=1"},
  { "b", 'b', "VALUE", 0, "PFCM b param. DEFAULT=2"},
  { "m", 'm', "VALUE", 0, "Set m fuzziness param. DEFAULT=1.2"},
  { "n", 'n', "VALUE", 0, "Set n fuzziness param. DEFAULT=1.2"},
  { "g", 'c', "VALUE", 0, "Set max clusters. DEFAULT=3"},
  { "k", 'k', 0, 0, "Set cosine norm. DEFAULT=euclidian"},
  { "d", 'd', 0, 0, "Set euclidian norm. DEFAULT=euclidian"},
  { "r", 'r', "VALUE", 0, "Set random initials. DEFAULT=5"},
  { 0 } 
};


static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  struct arguments *arguments = (struct arguments *)state->input;
  switch (key) {
    case 'f': arguments->mode = FCM; break;
    case 'p': arguments->mode = PCM; break;
    case 'x': arguments->mode = PFCM; break;
    case 'i': arguments->input = arg; break;
    case 'v': arguments->verbose = true; break;
    case 'a': arguments->a = atof(arg); break;
    case 'b': arguments->b = atof(arg); break;
    case 'm': arguments->m = atof(arg); break;
    case 'n': arguments->n = atof(arg); break;
    case 'c': arguments->c = atof(arg); break;
    case 'r': arguments->r = atof(arg); break;
    case 'd': arguments->norm = EUCLIDIAN; break;
    case 'k': arguments->norm = COSINE; break;
    case ARGP_KEY_ARG: return 0;
    default: return ARGP_ERR_UNKNOWN;
  }   
  return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };

int main(int argc, char *argv[]){

  arguments.input = NULL;

  arguments.mode = FCM;
  arguments.verbose = false;
  arguments.a = 1;
  arguments.b = 2;
  arguments.m = 1.2;
  arguments.n = 1.2;
  arguments.c = 3;
  arguments.r = 5;
  arguments.norm = EUCLIDIAN;

  argp_parse(&argp, argc, argv, 0, 0, &arguments);

  a = arguments.a;
  b = arguments.a;
  fuzziness = fuzziness_m = arguments.m;
  fuzziness_n = arguments.n;
  uint max_clusters = arguments.c;
  uint random_initials = arguments.r;

  //if(!arguments.input) {
  //  printf("--input FILE is required.\n");
  //  return 0;
  //}

  printf ("------------------------------------------------------------------------\n");
  if(arguments.mode == PCM){
    printf ("clustering data with pcm (possibilistic c means) method\n");
  }else if(arguments.mode == FCM){
    printf ("clustering data with fcm (fuzzy c means) method\n");
  }else if(arguments.mode == PFCM){
    printf ("clustering data with pfcm (possibilistic fuzzy c means) method\n");
  }
  string norm = "euclidian";
  string method = "fcm";
  if(arguments.norm == COSINE) norm = "cosine";
  if(arguments.mode == PCM) method = "pcm";
  if(arguments.mode == PFCM) method = "pfcm";
  printf ("parameters:\n");
  printf ("method\tm\tn\ta\tb\t\tnorm\n");
  printf ("%s\t%.2lf\t%.2lf\t%.2lf\t\t%.2lf\t%s\n", method.c_str(), fuzziness,
      fuzziness_n, a, b, norm.c_str());
  printf ("------------------------------------------------------------------------\n");

  printf("reading data\n");
  read_data();

  double max_fs = -2;
  double fs;
  uint max_groups = 2;
  printf("find optimal cluster number from 2 to %d\n", max_clusters);
  for(uint i = 2; i <= max_clusters; i++){
    printf("computing clustering with %d groups\n", i);
    num_clusters = i;
    times(j, random_initials){
      if(arguments.mode == PCM){
        pcm();
      }else if(arguments.mode == FCM){
        fcm();
      }else if(arguments.mode == PFCM){
        pfcm();
      }
      fs = aswc();
      if(arguments.verbose) cout << "FS: " << fs << endl;
      if(arguments.verbose) printf("%lf aswc %d groups\n", fs, i);
      if(fs > max_fs){
        max_fs = fs;
        max_groups = i;
        store_final_memberships();
        if(arguments.mode == PFCM){
          store_final_tipicalities();
        }
        //save_matrix("", final_memberships, num_docs);
        //save_matrix("", final_tipicalities, num_docs);
        if(arguments.verbose) printf("max fs found: %lf aswc %d groups\n", fs, i);
      }
    }
  }
  printf ("------------------------------------------------------------------------\n");
  printf ("saving matrix data...\n");
  if(arguments.mode == PCM){
    save_matrix("tipicalities.matrix", final_memberships, num_docs);
  }else if(arguments.mode == FCM){
    save_matrix("memberships.matrix", final_memberships, num_docs);
  }else if(arguments.mode == PFCM){
    save_matrix("tipicalties.matrix", final_tipicalities, num_docs);
    save_matrix("memberships.matrix", final_memberships, num_docs);
  }
  save_matrix("prototypes.matrix", prototypes, num_clusters);

  printf ("extracting descriptors...\n");
  soft_fdcl();

  printf ("the clustering process has finished\n");
  printf ("------------------------------------------------------------------------\n");
  printf ("max(aswc)\tgroups\n");
  printf ("%lf\t%d\n", max_fs, max_groups);
  printf ("------------------------------------------------------------------------\n");

  return 0;
}