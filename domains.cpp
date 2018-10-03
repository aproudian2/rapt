#include <Rcpp.h>
#include <stdlib.h>
using namespace std;

bool sortComplex(complex<double> a, complex<double> b);
bool sortVec(vector< complex<double> > a, vector< complex<double> > b);

// [[Rcpp::export]]
Rcpp::List domains(Rcpp::List rect, int bridge, int max_dom) {
  vector< vector< complex<double> > > dom_list(rect.length());
  Rcpp::List::iterator li;
  Rcpp::List ans;
  vector< vector< complex<double> > >::iterator i,j;
  for(li = rect.begin(), i = dom_list.begin(); li != rect.end(); ++li, ++i) {
    *i = Rcpp::as<vector< complex<double> > >(*li);
    sort(i->begin(), i->end(), sortComplex);
  }
  int n,m;
  vector< complex<double> > cv;
  vector< complex<double> >::iterator ci;
  for(i = dom_list.begin(); i != dom_list.end(); ++i) {
    for(j = i + 1; j != dom_list.end(); ++j) {
      n = i->size();
      m = j->size();
      (n>m) ? cv.resize(n) : cv.resize(m);
      ci = set_intersection(i->begin(), i->end(),
                            j->begin(), j->end(),
                            cv.begin(), sortComplex);
      cv.resize(ci - cv.begin());
      if(cv.size() > bridge) {
        cv.resize(i->size() + j->size());
        ci = set_union(i->begin(), i->end(),
                       j->begin(), j->end(),
                       cv.begin(), sortComplex);
        cv.resize(ci - cv.begin());
        *i = cv;
        sort(i->begin(), i->end(), sortComplex);
        dom_list.erase(j);
        sort(dom_list.begin(), dom_list.end(), sortVec);
        i = dom_list.begin() - 1;
        break;
      }
    }
  }
  sort(dom_list.begin(), dom_list.end(), sortVec);
  if(bridge > 0) {
    for(i = dom_list.begin(); i != dom_list.end(); ++i) {
      for(j = i + 1; j != dom_list.end(); ++j) {
        n = i->size();
        m = j->size();
        (n>m) ? cv.resize(n) : cv.resize(m);
        ci = set_intersection(i->begin(), i->end(),
                              j->begin(), j->end(),
                              cv.begin(), sortComplex);
        cv.resize(ci - cv.begin());
        if(!cv.empty()) {
          if (n + m - cv.size() > max_dom) {
            if (n > m) {
              ci = find(i->begin(), i->end(), cv[0]);
              i->erase(ci);
            } else {
              ci = find(j->begin(), j->end(), cv[0]);
              j->erase(ci);
            }
          } else {
            cv.resize(i->size() + j->size());
            ci = set_union(i->begin(), i->end(),
                           j->begin(), j->end(),
                           cv.begin(), sortComplex);
            cv.resize(ci - cv.begin());
            *i = cv;
            sort(i->begin(), i->end(), sortComplex);
            dom_list.erase(j);
            sort(dom_list.begin(), dom_list.end(), sortVec);
            i = dom_list.begin() - 1;
            break;
          }
        }
      }
    }
  }
  for(i = dom_list.begin(); i != dom_list.end(); ++i)
    ans.push_back(Rcpp::wrap(*i));
  return(ans);
}

bool sortComplex(complex<double> a, complex<double> b)
{
  if (real(a) == real(b))
    return imag(a) < imag(b);
  return real(a) < real(b);
}

bool sortVec(vector< complex<double> > a, vector< complex<double> > b) {
  return a.size() < b.size();
}
