#include <Rcpp.h>
#include <stdlib.h>
#include <math.h>
#define ROOT 0
#define HITS 0
#define BINS 1
#define MAX_SOL 60000
using namespace Rcpp;

/* Type Definitions */
typedef struct node {
  int index;
  int shift;
  int up, down;
} Node;
typedef struct column {
  int index;
  int prev, next;
  int first;
} Column;

/* Function Declarations */
void solve(Column* col_list, Node** row_list, int& ni, int& nb,
           int** choice, int level, int& solN,
           int* total_shift, int* total_energy,
           bool& con, int& min_shift, int& min_energy);
int next_col(Column* col_list, Node** row_list, int& ni);
void cover(Column* col_list, Node** row_list, int& cind);
void uncover(Column* col_list, Node** row_list, int& rind);
void print_row_list (Node** row_list);
void print_col_list (Column* col_list);

// [[Rcpp::export]]
List exact_cover(IntegerVector ind, IntegerVector bin, IntegerVector shift,
                 bool con) {
  int i = 0;
  int j = 0;
  int ni = unique(ind).size();
  int nb = unique(bin).size();
  int nr = ind.size();
  int ncl = ni + nb;
  int ri = ncl + 1;
  /* Initialize Column List */
  Column* col_list = new Column[ncl+1];
  col_list[ROOT].index = 0;
  col_list[ROOT].prev = ncl;
  col_list[ROOT].next = 1;
  col_list[ncl].index = ncl;
  col_list[ncl].prev = ncl-1;
  col_list[ncl].next = 0;
  for (i = 1; i < ncl; i++) {
    col_list[i].index = i;
    col_list[i].prev = i-1;
    col_list[i].next = i+1;
  }
  /* Initialize Row Header */
  Node** row_list = new Node* [ncl+nr+1];
    for (i = 0; i <= ncl+nr; i++)
      row_list[i] = new Node [2];
  for (i = 0; i <= ncl; i++) {
    row_list[i][HITS].index = i;
    row_list[i][BINS].index = i;
    row_list[i][HITS].up = 0;
    row_list[i][BINS].up = 0;
    row_list[i][HITS].down = 0;
    row_list[i][BINS].down = 0;
  }
  /* Initialize Row List */
  int last_node[ncl+1];
  for (i = 0; i <= ncl; i++)
    last_node[i] = i;
  for (i = 0, ri = ncl + 1; i < nr; i++, ri++) {
    row_list[ri][HITS].index = ind[i];
    row_list[ri][BINS].index = bin[i];
    row_list[ri][HITS].shift = shift[i];
    row_list[ri][BINS].shift = shift[i];
    for(j = 0; j < 2; j++) {
      row_list[ri][j].up = last_node[row_list[ri][j].index];
      row_list[row_list[ri][j].up][j].down = ri;
      last_node[row_list[ri][j].index] = ri;
    }
  }
  /* Close the Loops */
  for (i = 0; i <= ncl; i++) {
    if (i <= ni) {
      row_list[i][HITS].up = last_node[i];
      row_list[last_node[i]][HITS].down = i;
      col_list[i].first = row_list[i][HITS].down;
    }
    else {
      row_list[i][BINS].up = last_node[i];
      row_list[last_node[i]][BINS].down = i;
      col_list[i].first = row_list[i][BINS].down;
    }
  }
  int** choice = new int* [MAX_SOL];
  int level = 0;
  int solN = 0;
  choice[solN] = new int [nb];
  //float mshift = -114 + 25.7 * log(nr);
  int* total_shift = new int [MAX_SOL];
  int min_shift = INT_MAX - 50;
  if (min_shift < 0) {
    min_shift = nb;
  }
  int* total_energy = new int [MAX_SOL];
  int min_energy = INT_MAX - 50;
  for(i = 0; i < MAX_SOL; i++) {
    total_shift[i] = 0;
    total_energy[i] = 0;
  }
  //print_row_list(row_list); // print-out
  solve(col_list, row_list, ni, nb, choice, level, solN,
        total_shift, total_energy, con, min_shift, min_energy);
  IntegerMatrix cover_mat(solN, nb);
  IntegerVector energy(solN);
  IntegerVector shift_vec(solN);
  for (i = 0; i < solN; i++) {
    shift_vec(i) = total_shift[i];
    energy(i) = total_energy[i];
    for (j = 0; j < nb; j++) {
      cover_mat(i,j) = choice[i][j] - ncl;
    }
  }
  List ans;
  ans["cover"] = cover_mat;
  ans["energy"] = energy;
  ans["shift"] = shift_vec;
  delete col_list;
  delete row_list;
  delete choice;
  delete total_shift;
  return ans;
}

/* Subroutines */
void solve(Column* col_list, Node** row_list, int& ni, int& nb,
           int** choice, int level, int& solN,
           int* total_shift, int* total_energy,
           bool& con, int& min_shift, int& min_energy) {
  int rind = next_col(col_list, row_list, ni);
  int i = 0;
  bool limit = FALSE;
  // print_col_list(col_list); // print-out
  if (rind == 0)
    return;
  int cind = row_list[rind][HITS].index;
  while (rind != cind) {
    //Rprintf("%d:",level); // print-out
    choice[solN][level] = rind;
    total_shift[solN] += row_list[rind][HITS].shift;
    if (level > 0)
      total_energy[solN] += row_list[choice[solN][level]][HITS].index -
                            row_list[choice[solN][level-1]][HITS].index - 1;
    if (con) {
      if (total_shift[solN] > min_shift ||
          total_energy[solN] > min_energy)
        limit = TRUE;
    }
    if (limit) {
      total_shift[solN] -= row_list[rind][HITS].shift;
      if (level > 0)
        total_energy[solN] -= row_list[choice[solN][level]][HITS].index -
          row_list[choice[solN][level-1]][HITS].index - 1;
      return;
    }
    //Rprintf("%d\n",rind); // print-out
    //Rprintf("Covering:\n"); // print-out
    cover(col_list, row_list, rind);
    level++;
    if (level == nb) {
      //Rprintf("Uncovering:\n"); // print-out
      uncover(col_list, row_list, rind);
      //print_col_list(col_list); // print-out
      solN++;
      //Rprintf("-----\n"); // print-out
      if (solN == MAX_SOL)
        return;
      choice[solN] = new int [nb];
      total_shift[solN] = total_shift[solN-1];
      total_energy[solN] = total_energy[solN-1];
      if (con) {
        if (total_energy[solN] < min_energy) {
          min_energy = total_energy[solN];
          //Rprintf("New min energy = %d\n", min_energy); // print-out
        }
        if (total_shift[solN] < min_shift) {
          min_shift = total_shift[solN];
        }
      }
      for (i = 0; i < nb; i++)
        choice[solN][i] = choice[solN-1][i];
      total_shift[solN] -= row_list[rind][HITS].shift;
      total_energy[solN] -= row_list[choice[solN][level-1]][HITS].index -
                            row_list[choice[solN][level-2]][HITS].index - 1;
      return;
    }
    solve(col_list, row_list, ni, nb, choice, level, solN,
          total_shift, total_energy,
          con, min_shift, min_energy);
    //Rprintf("Uncovering:\n"); // print-out
    uncover(col_list, row_list, rind);
    // print_col_list(col_list); // print-out
    level--;
    if (solN == MAX_SOL)
      return;
    total_shift[solN] -= row_list[rind][HITS].shift;
    if (level > 0)
      total_energy[solN] -= row_list[choice[solN][level]][HITS].index -
                            row_list[choice[solN][level-1]][HITS].index - 1;
    next_col(col_list, row_list, ni);
    rind = row_list[rind][HITS].down;
  }
  return;
}
int next_col(Column* col_list, Node** row_list, int& ni) {
  int cind, rhead, row;
  do {
    cind = col_list[ROOT].next;
    rhead = col_list[cind].index;
    row = row_list[rhead][HITS].down;
    col_list[col_list[cind].prev].next = col_list[cind].next;
    col_list[col_list[cind].next].prev = col_list[cind].prev;
    if (row != rhead)
      // Rprintf("Column %d\n",cind); // print-out
      return row;
  } while (cind <= ni);
  return 0;
}
void cover(Column* col_list, Node** row_list, int& rind) {
  int i, j;
  int ind;
  int cind[2] = {row_list[rind][HITS].index, row_list[rind][BINS].index};
  for (i = 0; i < 2; i++) {
    ind = rind;
    while (ind != cind[i]) {
      for (j = 0; j < 2; j++) {
        row_list[row_list[ind][j].up][j].down = row_list[ind][j].down;
        row_list[row_list[ind][j].down][j].up = row_list[ind][j].up;
      }
      ind = row_list[ind][i].down;
    }
    ind = row_list[rind][i].up;
    while (ind != cind[i]) {
      for (j = 0; j < 2; j++) {
        row_list[row_list[ind][j].up][j].down = row_list[ind][j].down;
        row_list[row_list[ind][j].down][j].up = row_list[ind][j].up;
      }
      ind = row_list[ind][i].up;
    }
    //print_row_list(row_list); // print-out
  }
}
void uncover(Column* col_list, Node** row_list, int& rind) {
  int i, j;
  int ind;
  int cind[2] = {row_list[rind][HITS].index, row_list[rind][BINS].index};
  for (i = 1; i >= 0; i--) {
    ind = rind;
    //Rprintf("%d:",rind); // print-out
    while (ind != cind[i]) {
      //Rprintf("%d,",ind); // print-out
      for (j = 1; j >= 0; j--) {
          row_list[row_list[ind][j].down][j].up = ind;
          row_list[row_list[ind][j].up][j].down = ind;
      }
      ind = row_list[ind][i].up;
    }
    ind = row_list[rind][i].down;
    while (ind != cind[i]) {
      //Rprintf("%d,",ind); // print-out
      for (j = 1; j >= 0; j--) {
          row_list[row_list[ind][j].down][j].up = ind;
          row_list[row_list[ind][j].up][j].down = ind;
      }
      ind = row_list[ind][i].down;
    }
    //Rprintf("\n");  // print-out
    //print_row_list(row_list); // print-out
  }
  col_list[col_list[cind[HITS]].next].prev = cind[HITS];
  col_list[col_list[cind[HITS]].prev].next = cind[HITS];
}
void print_row_list (Node** row_list) {
  for (int i = 0; i <= 15; i++) { // hardcoded... :(
    Rprintf("%d:\t",i);
    for (int j = 0; j < 2; j++) {
      Rprintf("index=%d,up=%d,down=%d\t",
              row_list[i][j].index,row_list[i][j].up,row_list[i][j].down);
    }
    Rprintf("\n");
  }
  Rprintf("\n");
}
void print_col_list (Column* col_list) {
  int col = ROOT;
  do {
    Rprintf("index=%d,prev=%d,next=%d\n",
            col,col_list[col].prev,col_list[col].next);
    col = col_list[col].next;
  } while (col != 0);
  Rprintf("\n");
}
