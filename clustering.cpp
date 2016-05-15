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

  argp_parse(&argp, argc, argv, 0, 0, &arguments);

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
  printf ("------------------------------------------------------------------------\n");

  read_data();

  double max_fs = -1;
  double fs;
  int max_groups = 2;
  for(int i = 3; i <= 3; i++){
    printf("computing clustering with %d groups\n", i);
    num_clusters = i;
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
