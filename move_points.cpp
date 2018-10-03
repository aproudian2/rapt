#include <Rcpp.h>
#include <stdlib.h>
using namespace std;

// [[Rcpp::export]]

Rcpp::DataFrame move_points(Rcpp::IntegerMatrix match,
                      Rcpp::ComplexVector source, Rcpp::ComplexVector sink,
                      Rcpp::List domain, Rcpp::DataFrame points) {
  Rcpp::DataFrame ans = Rcpp::clone(points);
  int i, j, n;
  n = match.nrow();
  int ss_link[n][2];
  for(i = 0; i < n; i++) {
    for(j = 0; j < 2; j++) {
      ss_link[i][j] = match(i,j) - 1;
    }
  }
  vector< complex<double> > so = Rcpp::as< vector <complex<double> > >(source);
  vector< complex<double> > si = Rcpp::as< vector <complex<double> > >(sink);
  Rcpp::List::iterator li;
  vector< vector< complex<double> > > dom_list(domain.size());
  vector< vector< complex<double> > >::iterator di;
  for(li = domain.begin(), di = dom_list.begin();
      li != domain.end(); ++li, ++di)
    *di = Rcpp::as<vector< complex<double> > >(*li);
  vector< complex<double> >::iterator soi, sii, bi;
  vector< complex<double> > bin =
    Rcpp::as< vector <complex<double> > >(points["bin"]);
  for(i = 0; i < n; i++) {
    // Rcpp::Rcout << so[ss_link[i][0]] << "," << si[ss_link[i][1]] << "\n\t";
    for(di = dom_list.begin(); di != dom_list.end(); ++di) {
      // for(bi = di->begin(); bi != di->end(); ++bi)
        // Rcpp::Rcout << *bi;
      // Rcpp::Rcout << "\n\t";
      soi = find(di->begin(), di->end(), so[ss_link[i][0]]);
      sii = find(di->begin(), di->end(), si[ss_link[i][1]]);
      if(soi != di->end() || sii != di->end()) {
        // Rcpp::Rcout << "\r match " << i << ". ";
        if(soi == di->end() || sii == di->end()) {
          bi = find(bin.begin(), bin.end(), so[ss_link[i][0]]);
          *bi = si[ss_link[i][1]];
          // Rcpp::Rcout << "move!";
        }
        // Rcpp::Rcout << "\n\n";
        break;
      }
    }
  }
  ans["bin"] = Rcpp::wrap(bin);
  return(ans);
}
