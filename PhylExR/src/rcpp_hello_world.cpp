// [[Rcpp::plugins("cpp17")]]

#include <Rcpp.h>
#include <deque>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "node.hpp"

using namespace Rcpp;
using namespace std;

double log_add(double x, double y)
{
  // make x the max
  if (y > x) {
    double temp = x;
    x = y;
    y = temp;
  }
  // now x is bigger
  if (x == R_NegInf) {
    return x;
  }
  double negDiff = y - x;
  if (negDiff < -20) {
    return x;
  }
  return x + log(1.0 + exp(negDiff));
}

//' @export
// [[Rcpp::export]]
String ReadParentVector(String file_path, size_t mutation_count) {
  //Rcout << "Reading...\n";

  ifstream dat_file (file_path);

  NumericVector vec(mutation_count);
  string line;
  // Optimal parent vector is on line 2 * mutation_count + 4.
  size_t line_no = 2 * mutation_count + 4;
  for (size_t i = 0; i < line_no; i++) {
    getline(dat_file, line);
  }
  //Rcout << line_no << ": " << line << "\n";
  dat_file.close();
  return line;
}

//' @export
// [[Rcpp::export]]
NumericVector GetChains(NumericVector parent_vector) {
  Rcout << "Identifying chains from mutation tree...\n";
  size_t node_count = parent_vector.size();
  vector<Node *> nodes;
  NumericVector chain_labels(node_count);
  for (size_t i = 0; i < node_count; i++) {
    nodes.push_back(new Node(i));
  }
  Node *root = new Node(node_count);
  nodes.push_back(root);
  for (size_t i = 0; i < node_count; i++) {
    size_t parent_idx = parent_vector[i];
    nodes[parent_idx]->add_child(nodes[i]);
    //Rcout << "node: " << nodes[i]->get_id() << ", parent: " << parent_idx << "\n";
  }

  size_t chain_id = 0;
  deque<Node *> queue;
  for (auto child : root->get_children()) {
      queue.push_back(child);
  }
  while (!queue.empty()) {
    auto node = queue.back();
    //Rcout << node->get_id() << "\n";
    queue.pop_back();
    while (true) {
        chain_labels(node->get_id()) = chain_id;
        //Rcout << "num children: " << node->get_children().size() << "\n";
        if (node->get_children().size() == 1) {
            node = node->get_children()[0];
        } else {
            break;
        }
    }
    chain_labels(node->get_id()) = chain_id;
    chain_id++;

    for (auto child_node : node->get_children()) {
        queue.push_back(child_node);
    }
  }

  return chain_labels;
}

void split_string(std::string str, std::vector<std::string> &vec, char sep)
{
  vec.clear();
  std::istringstream str_stream(str);
  std::string token;
  while (getline(str_stream, token, sep)) {
    vec.push_back(token);
  }
}

// [[Rcpp::export]]
std::string get_parent_name(std::string node) {
  std::vector<std::string> split_vec;
  split_string(node, split_vec, '_');
  std::string parent_name = "";
  size_t parent_name_len = split_vec.size() - 1;
  for (size_t i = 0; i < parent_name_len; i++) {
    parent_name += split_vec[i];
    if (i < parent_name_len - 1) {
      parent_name += "_";
    }
  }

  return parent_name;
}

std::unordered_map<std::string, std::unordered_set<std::string>> ProcessDatum2Node(DataFrame datum2node)
{
  std::unordered_map<std::string, std::unordered_set<std::string>> node2data;
  size_t n_rows = datum2node.nrows();
  CharacterVector snvs = datum2node[0];
  CharacterVector nodes = datum2node[1];
  for (size_t i = 0; i < n_rows; i++) {
    std::string snv = as<std::string>(snvs[i]);
    std::string node = as<std::string>(nodes[i]);
    node2data[node].insert(snv);
  }
  return node2data;
}

//' @export
// [[Rcpp::export]]
NumericMatrix GetConfigMatrix(DataFrame datum2node, std::vector<std::string> ordered_nodes)
{
  size_t n_clones = ordered_nodes.size();
  size_t n_snvs = datum2node.nrows();

  // Get SNVs for each node.
  std::unordered_map<std::string, std::unordered_set<std::string>> node2data = ProcessDatum2Node(datum2node);
  for (std::string node : ordered_nodes) {
    std::string parent_node = get_parent_name(node);
    if (parent_node == "0") {
      continue;
    }
    node2data[node].insert(node2data[parent_node].begin(), node2data[parent_node].end());
  }

  // Config matrix.
  NumericMatrix config(n_snvs, n_clones);

  CharacterVector snvs = datum2node[0];
  // Map from SNV name to idx.
  unordered_map<string, size_t> snv2idx;
  for (size_t i = 0; i < n_snvs; i++) {
    snv2idx[as<std::string>(snvs[i])] = i;
  }

  for (size_t col_idx = 0; col_idx < n_clones; col_idx++) {
    for (string snv : node2data[ordered_nodes[col_idx]]) {
      size_t row_idx = snv2idx[snv];
      config(row_idx, col_idx) = 1;
    }
  }

  return config;
}

// [[Rcpp::export]]
double LogAdd(double x, double y)
{
  // make x the max
  if (y > x) {
    double temp = x;
    x = y;
    y = temp;
  }
  // now x is bigger
  if (x == R_NegInf) {
    return x;
  }
  double negDiff = y - x;
  if (negDiff < -20) {
    return x;
  }
  return x + log(1.0 + exp(negDiff));
}

// [[Rcpp::export]]
double LogSumExp(NumericVector x)
{
  double max = R_NegInf;
  double maxIndex = 0;
  for (unsigned int i = 0; i < x.size(); i++)
  {
    if (x[i] > max) {
      max = x[i];
      maxIndex = i;
    }
  }
  if (max == R_NegInf) return R_NegInf;
  // compute the negative difference
  double threshold = max - 20;
  double sumNegativeDifferences = 0.0;
  for (unsigned int i = 0; i < x.size(); i++) {
    if (i != maxIndex && x[i] > threshold) {
      sumNegativeDifferences += exp(x[i] - max);
    }
  }
  if (sumNegativeDifferences > 0.0) {
    return max + log(1.0 + sumNegativeDifferences);
  } else {
    return max;
  }
}

// [[Rcpp::export]]
NumericMatrix IdentifyNodeMutationStatus(DataFrame datum2node,
                                         const std::vector<std::string> &ordered_nodes,
                                         const std::vector<std::string> &ordered_mutations)
{
  size_t node_count = ordered_nodes.size();
  size_t snv_count = ordered_mutations.size();
  NumericMatrix node2snv(node_count, snv_count);
  
  std::unordered_set<std::string> empty_set;
  std::unordered_map<std::string, std::unordered_set<std::string>> node2data = ProcessDatum2Node(datum2node);
  node2data["0"] = empty_set;
  for (std::string node : ordered_nodes) {
    std::string parent_node = get_parent_name(node);
    if (parent_node == "0") {
      continue;
    }
    node2data[node].insert(node2data[parent_node].begin(), node2data[parent_node].end());
  }
  
  for (size_t i = 0; i < ordered_nodes.size(); i++) {
    for (size_t j = 0; j < ordered_mutations.size(); j++) {
      if (node2data[ordered_nodes.at(i)].count(ordered_mutations.at(j))) {
        node2snv(i,j) = 1;
      } else {
        node2snv(i,j) = 0;
      }
    }
  }
  return node2snv;
}
