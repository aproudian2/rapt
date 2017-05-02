#include <Rcpp.h>
#include <stdlib.h>
#define ROOT 0
#define HITS 0
#define BINS 1
using namespace Rcpp;

/* Type Definitions */
typedef struct node {
  int index;
  int up, down;
} Node;

typedef struct column {
  int index;
  int prev, next;
} Column;

/* Function Declarations */
int next_col(Column* col_list, Node** row_list, int& ni);
void cover(Node** row_list, int& cind, int* sol_row, int& ni);
int uncover(Node** row_list, int& cind);

// [[Rcpp::export]]
IntegerVector exact_cover(IntegerVector ind, IntegerVector bin) {
  int i = 0;
  int j = 0;
  int s = 0;
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
  for (i=1; i<ncl; i++) {
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
    for(j = 0; j < 2; j++) {
      row_list[ri][j].up = last_node[row_list[ri][j].index];
      row_list[row_list[ri][j].up][j].down = ri;
      last_node[row_list[ri][j].index] = ri;
    }
  }
  for (i = 0; i <= ncl; i++) {
    for (j = 0; j < 2; j++) {
      row_list[i][j].up = last_node[i];
    }
    if (i > ni)
      row_list[last_node[i]][BINS].down = i;
    else
      row_list[last_node[i]][HITS].down = i;
  }
  int cind[nb];
  int** sol = new int* [nb];
  sol[s] = new int [2];
  for (i = 0; i < nb; i++)
    sol[s][i] = 0;
  for(i = 0; i < nb; i++) {
    cind[i] = next_col(col_list, row_list, ni);
    cover(row_list, cind[i], sol[s], ni);
  }
  //IntegerVector ans = IntegerVector(cind, cind + nb);
  IntegerVector ans = IntegerVector(sol[s], sol[s] + nb);
  return ans;
}

/* Subroutines */
int next_col(Column* col_list, Node** row_list, int &ni) {
  int cind, rhead, row;
  do {
    cind = col_list[ROOT].next;
    rhead = col_list[cind].index;
    row = row_list[rhead][HITS].down;
    col_list[ROOT].next = col_list[cind].next;
    if (row != rhead)
      return cind;
  } while (cind <= ni);
  return 0;
}

void cover(Node** row_list, int& cind, int* sol_row, int& ni) {
  int i,j;
  int ind[2];
  int rind = row_list[cind][HITS].down;
  int bind = row_list[rind][BINS].index;
  sol_row[bind-(ni+1)] = cind;
  do {
    ind[HITS] = row_list[cind][HITS].down;
    ind[BINS] = row_list[bind][BINS].down;
    for(i = 0; i < 2; i++) {
      for(j = 0; j < 2; j++) {
        row_list[row_list[ind[i]][j].down][j].up =
          row_list[ind[i]][j].up;
        row_list[row_list[ind[i]][j].up][j].down =
          row_list[ind[i]][j].down;
      }
    }
  } while ((ind[HITS] != cind) || (ind[BINS] != bind));
}

//int uncover(Node ** row_list, int& cind) {

//};
