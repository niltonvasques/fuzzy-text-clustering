#include <bits/stdc++.h>
#include "shared.h"
#include "fcm.h"
#include "pcm.h"
#include "pfcm.h"
#include "descriptors.h"

using namespace std;

int main(){

  read_data();

  double max_fs = -1;
  double fs;
  int max_groups = 2;
  for(int i = 3; i <= 3; i++){
    num_clusters = i;
    pfcm();
    fs = aswc();
    cout << "FS: " << fs << endl;
    printf("%lf aswc %d groups\n", fs, i);
    if(fs > max_fs){
      max_fs = fs;
      max_groups = i;
      store_final_memberships();
      store_final_tipicalities();
      //save_matrix("", final_memberships, num_docs);
      //save_matrix("", final_tipicalities, num_docs);
      printf("max fs found: %lf aswc %d groups\n", fs, i);
    }
  }
  printf("final fs: %lf aswc %d groups\n", max_fs, max_groups);
  save_matrix("tipicalties.matrix", final_tipicalities, num_docs);
  save_matrix("memberships.matrix", final_memberships, num_docs);
  save_matrix("prototypes.matrix", prototypes, num_clusters);
  soft_fdcl();

  return 0;
}
