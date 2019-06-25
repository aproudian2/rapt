#ifndef READER_H
#define READER_H

#include <cstdlib>
#include <tuple>
#include <vector>
#include <fstream>
using namespace std;

// [[Rcpp::plugins(cpp11)]]

// Defined based on v1.0 TAPSim spec in howto.pdf included with distribution.
class Results;

class Result {
    int index, id;
    unsigned int number;
    float voltage;
    tuple<float,float,float> start;
    tuple<float,float,float> stop;
    float tof;
    float probability;
    float potentialBefore;
    tuple<float,float,float> fieldBefore;
    float potentialAfter;
    tuple<float,float,float> fieldAfter;
    tuple<float,float,float> normal;
    tuple<float,float,float> apex;
  public:
    Result (int index,
            int id,
            unsigned int number,
            float voltage,
            tuple<float,float,float> start,
            tuple<float,float,float> stop,
            float tof,
            float probability,
            float potentialBefore,
            tuple<float,float,float> fieldBefore,
            float potentialAfter,
            tuple<float,float,float> fieldAfter,
            tuple<float,float,float> normal,
            tuple<float,float,float> apex);
    Result (Results, unsigned int);
    int get_index() {return index;}
    int get_id() {return id;}
    unsigned int get_number() {return number;}
    float get_voltage() {return voltage;}
    tuple<float,float,float> get_start() {return start;}
    float get_startX() {return get<0>(start);}
    float get_startY() {return get<1>(start);}
    float get_startZ() {return get<2>(start);}
    tuple<float,float,float> get_stop() {return stop;}
    float get_stopX() {return get<0>(stop);}
    float get_stopY() {return get<1>(stop);}
    float get_stopZ() {return get<2>(stop);}
    float get_tof() {return tof;}
    float get_probability() {return probability;}
    float get_potentialBefore() {return potentialBefore;}
    tuple<float,float,float> get_fieldBefore() {return fieldBefore;}
    float get_fieldBeforeX() {return get<0>(fieldBefore);}
    float get_fieldBeforeY() {return get<1>(fieldBefore);}
    float get_fieldBeforeZ() {return get<2>(fieldBefore);}
    float get_potentialAfter() {return potentialAfter;}
    tuple<float,float,float> get_fieldAfter() {return fieldAfter;}
    float get_fieldAfterX() {return get<0>(fieldAfter);}
    float get_fieldAfterY() {return get<1>(fieldAfter);}
    float get_fieldAfterZ() {return get<2>(fieldAfter);}
    tuple<float,float,float> get_normal() {return normal;}
    float get_normalX() {return get<0>(normal);}
    float get_normalY() {return get<1>(normal);}
    float get_normalZ() {return get<2>(normal);}
    tuple<float,float,float> get_apex() {return apex;}
    float get_apexX() {return get<0>(apex);}
    float get_apexY() {return get<1>(apex);}
    float get_apexZ() {return get<2>(apex);}
    void print();
};

class Tapsim {
    char* filepath;
    bool binary;
    unsigned int head_len;
    unsigned int head_lines;
    unsigned int entries;
  public:
    Tapsim (char* fp);
    char* get_filepath() {return filepath;}
    bool get_binary() {return binary;}
    unsigned int get_head_len() {return head_len;}
    unsigned int get_head_lines() {return head_lines;}
    unsigned int get_entries() {return entries;}
    void print();
};

class Results : public Tapsim {
    vector<Result> results;
  public:
    Results (char* fp);
    Result get_result(unsigned int n) {return results[n];}
    void print_result(unsigned int n);
};

#endif
