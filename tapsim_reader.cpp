#include <fstream>
#include <Rcpp.h>
#include "tapsim_reader.h"
#define R_BYTES 104

using namespace std;
// [[Rcpp::plugins(cpp11)]]

Result::Result (Results results, unsigned int n) {
  char* fp = results.get_filepath();
  if (results.get_binary()) {
    printf("Binary file.\n");
    ifstream file(fp, ios::binary);
    unsigned int fwd = results.get_head_len() + n * R_BYTES;
    file.seekg(fwd);

    file.read(reinterpret_cast<char*>(&index), sizeof(int));
    file.read(reinterpret_cast<char*>(&id), sizeof(int));
    file.read(reinterpret_cast<char*>(&number), sizeof(int));
    file.read(reinterpret_cast<char*>(&voltage), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<0>(start)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<1>(start)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<2>(start)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<0>(stop)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<1>(stop)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<2>(stop)), sizeof(float));
    file.read(reinterpret_cast<char*>(&tof), sizeof(float));
    file.read(reinterpret_cast<char*>(&probability), sizeof(float));
    file.read(reinterpret_cast<char*>(&potentialBefore), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<0>(fieldBefore)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<1>(fieldBefore)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<2>(fieldBefore)), sizeof(float));
    file.read(reinterpret_cast<char*>(&potentialAfter), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<0>(fieldAfter)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<1>(fieldAfter)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<2>(fieldAfter)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<0>(normal)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<1>(normal)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<2>(normal)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<0>(apex)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<1>(apex)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<2>(apex)), sizeof(float));

    file.close();
  }
  else {
    ifstream file(fp);
    string line;
    unsigned int beg = 0;
    unsigned int end = 0;
    unsigned int fwd = results.get_head_lines();
    for (int i = 0; i <= fwd + n; i++) {
      getline(file, line);
    }

    getline(file, line);

    end = line.find("\t", beg);
    index = stoi(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    id = stoi(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    number = stoul(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    voltage = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<0>(start) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<1>(start) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<2>(start) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<0>(stop) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<1>(stop) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<2>(stop) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    tof = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    probability = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    potentialBefore = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<0>(fieldBefore) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<1>(fieldBefore) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<2>(fieldBefore) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    potentialAfter = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<0>(fieldAfter) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<1>(fieldAfter) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<2>(fieldAfter) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<0>(normal) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<1>(normal) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<2>(normal) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<0>(apex) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<1>(apex) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<2>(apex) = stof(line.substr(beg, end));
    beg = end + 1;

  file.close();
  }
}
void Result::print(void) {
  printf("Index\tID\tNumber\n");
  printf("%i\t%i\t%i\n", index, id, number);
  printf("Voltage\n");
  printf("%e\n", voltage);
  printf("Start X\t\tStart Y\t\tStart Z\n");
  printf("%e\t%e\t%e\n", get<0>(start), get<1>(start), get<2>(start));
  printf("Stop X\t\tStop Y\t\tStop Z\n");
  printf("%e\t%e\t%e\n", get<0>(stop), get<1>(stop), get<2>(stop));
  printf("TOF\t\tProbability\tPotential Before\n");
  printf("%e\t%e\t%e\n", tof, probability, potentialBefore);
  printf("Field Before X\tField Before Y\tField Before Z\n");
  printf("%e\t%e\t%e\n", get<0>(fieldBefore), get<1>(fieldBefore),
         get<2>(fieldBefore));
  printf("Potential After\n");
  printf("%e\n", potentialAfter);
  printf("Field After X\tField After Y\tField After Z\n");
  printf("%e\t%e\t%e\n", get<0>(fieldAfter), get<1>(fieldAfter),
         get<2>(fieldAfter));
  printf("Normal X\tNormal Y\tNormal Z\n");
  printf("%e\t%e\t%e\n", get<0>(normal), get<1>(normal), get<2>(normal));
  printf("Apex X\t\tApex Y\t\tApex Z\n");
  printf("%e\t%e\t%e\n", get<0>(apex), get<1>(apex), get<2>(apex));
}

Tapsim::Tapsim (char* fp) {
  filepath = fp;
  ifstream file(fp);
  head_len = 0;
  head_lines = 0;
  string buff;
  bool head = true;
  int i = 0;
  for (i = 0; head && i < 50; i++) {
    getline(file, buff);
    if (buff.find("ASCII") != string::npos) {
      head = false;
      binary = false;
      break;
    }
    else if (buff.find("BINARY") != string::npos) {
      head = false;
      binary = true;
      break;
    }
    else {
      head_len += buff.length();
    }
  }
  head_lines = i;
  if (i >= 50) {
    printf("Error reading header.\n");
    printf("Head Length Read: %u\n", head_len);
    return;
  }
  if (binary == true) {
    string s_entries = buff.substr(7);
    entries = stoi(s_entries);
  }
  else {
    for (i = 0; !file.eof(); i++) {
      getline(file, buff);
    }
    entries = i - 1;
  }
}
void Tapsim::print() {
  printf("Filepath: %s\n", filepath);
  printf("Header Length: %u\n", head_len);
  printf("Header Lines: %u\n", head_lines);
  printf("Binary: %s\n", binary ? "true" : "false");
  printf("Entries: %u\n", entries);
}

Results::Results (char* fp) : Tapsim (fp) {
  for (int i = 0; i < this->get_entries(); i++) {
    results.push_back(Result(*this, i));
  }
}
void Results::print_result(unsigned int n) {
  Result result = this->get_result(n);
  result.print();
}

Rcpp::DataFrame readResults(string fp) {
  char* c_fp = new char [fp.length()+1];
  strcpy(c_fp, fp.c_str());
  Results c_data(c_fp);
  Rcpp::IntegerVector index_vec;
  Rcpp::IntegerVector id_vec;
  Rcpp::IntegerVector number_vec;
  Rcpp::NumericVector voltage_vec;
  Rcpp::NumericVector startX_vec;
  Rcpp::NumericVector startY_vec;
  Rcpp::NumericVector startZ_vec;
  Rcpp::NumericVector stopX_vec;
  Rcpp::NumericVector stopY_vec;
  Rcpp::NumericVector stopZ_vec;
  Rcpp::NumericVector tof_vec;
  Rcpp::NumericVector prob_vec;
  Rcpp::NumericVector potentialBefore_vec;
  Rcpp::NumericVector fieldBeforeX_vec;
  Rcpp::NumericVector fieldBeforeY_vec;
  Rcpp::NumericVector fieldBeforeZ_vec;
  Rcpp::NumericVector potentialAfter_vec;
  Rcpp::NumericVector fieldAfterX_vec;
  Rcpp::NumericVector fieldAfterY_vec;
  Rcpp::NumericVector fieldAfterZ_vec;
  Rcpp::NumericVector normalX_vec;
  Rcpp::NumericVector normalY_vec;
  Rcpp::NumericVector normalZ_vec;
  Rcpp::NumericVector apexX_vec;
  Rcpp::NumericVector apexY_vec;
  Rcpp::NumericVector apexZ_vec;

  unsigned int entries = c_data.get_entries();
  for(int i = 0; i < entries; i++) {
    Result result = c_data.get_result(i);
    index_vec.push_back(result.get_index());
    id_vec.push_back(result.get_id());
    number_vec.push_back(result.get_number());
    voltage_vec.push_back(result.get_voltage());
    startX_vec.push_back(result.get_startX());
    startY_vec.push_back(result.get_startY());
    startZ_vec.push_back(result.get_startZ());
    stopX_vec.push_back(result.get_stopX());
    stopY_vec.push_back(result.get_stopY());
    stopZ_vec.push_back(result.get_stopZ());
    tof_vec.push_back(result.get_tof());
    prob_vec.push_back(result.get_probability());
    potentialBefore_vec.push_back(result.get_potentialBefore());
    fieldBeforeX_vec.push_back(result.get_fieldBeforeX());
    fieldBeforeY_vec.push_back(result.get_fieldBeforeY());
    fieldBeforeZ_vec.push_back(result.get_fieldBeforeZ());
    potentialAfter_vec.push_back(result.get_potentialAfter());
    fieldAfterX_vec.push_back(result.get_fieldAfterX());
    fieldAfterY_vec.push_back(result.get_fieldAfterY());
    fieldAfterZ_vec.push_back(result.get_fieldAfterZ());
    normalX_vec.push_back(result.get_normalX());
    normalY_vec.push_back(result.get_normalY());
    normalZ_vec.push_back(result.get_normalZ());
    apexX_vec.push_back(result.get_apexX());
    apexY_vec.push_back(result.get_apexY());
    apexZ_vec.push_back(result.get_apexZ());
  }

  Rcpp::StringVector col_names;
  col_names.push_back("Index");
  col_names.push_back("ID");
  col_names.push_back("Number");
  col_names.push_back("Voltage");
  col_names.push_back("StartX");
  col_names.push_back("StartY");
  col_names.push_back("StartZ");
  col_names.push_back("StopX");
  col_names.push_back("StopY");
  col_names.push_back("StopZ");
  col_names.push_back("TOF");
  col_names.push_back("Prob");
  col_names.push_back("PotBef");
  col_names.push_back("FieldBefX");
  col_names.push_back("FieldBefY");
  col_names.push_back("FieldBefZ");
  col_names.push_back("PotAft");
  col_names.push_back("FieldAftX");
  col_names.push_back("FieldAftY");
  col_names.push_back("FieldAftZ");
  col_names.push_back("NormalX");
  col_names.push_back("NormalY");
  col_names.push_back("NormalZ");
  col_names.push_back("ApexX");
  col_names.push_back("ApexY");
  col_names.push_back("ApexZ");

  Rcpp::List r_data(26);
  r_data(0) = index_vec;
  r_data(1) =  id_vec;
  r_data(2) = number_vec;
  r_data(3) = voltage_vec;
  r_data(4) = startX_vec;
  r_data(5) = startY_vec;
  r_data(6) = startZ_vec;
  r_data(7) = stopX_vec;
  r_data(8) = stopY_vec;
  r_data(9) = stopZ_vec;
  r_data(10) = tof_vec;
  r_data(11) = prob_vec;
  r_data(12) = potentialBefore_vec;
  r_data(13) = fieldBeforeX_vec;
  r_data(14) = fieldBeforeY_vec;
  r_data(15) = fieldBeforeZ_vec;
  r_data(16) = potentialAfter_vec;
  r_data(17) = fieldAfterX_vec;
  r_data(18) = fieldAfterY_vec;
  r_data(19) = fieldAfterZ_vec;
  r_data(20) = normalX_vec;
  r_data(21) = normalY_vec;
  r_data(22) = normalZ_vec;
  r_data(23) = apexX_vec;
  r_data(24) = apexY_vec;
  r_data(25) = apexZ_vec;
  r_data.attr("names") = col_names;
  r_data.attr("row.names") =
    Rcpp::IntegerVector::create(NA_INTEGER, entries);
  r_data.attr("class") = "data.frame";

  return r_data;
}

// [[Rcpp::export]]
Rcpp::DataFrame read_TAPSim(string fp, string type) {
  if (type == "results") {
    Rcpp::DataFrame results = readResults(fp);
    return results;
  }
  else {
    Rprintf("Not yet implemented.");
    return NULL;
  }
}
