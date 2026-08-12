#include <Rcpp.h>

// [[Rcpp::plugins("cpp11")]]
using namespace Rcpp;

void partitionStudents(
  const std::vector<int>& students,
  const std::vector<int>& n_proposals_made,
  const IntegerMatrix& s_prefs,
  std::vector<int>& temp_singles,
  std::vector<int>& final_singles
);

void removeInterruptingPair(
  IntegerMatrix& s_prefs,
  int student,
  int college
);

// [[Rcpp::export]]
List eadam_cpp(
  IntegerMatrix s_prefs,
  IntegerMatrix c_prefs,
  IntegerVector n_slots,
  String acceptance,
  LogicalVector consent,
  int bound_ea_iters
) {
  int n_students = s_prefs.ncol();
  int n_colleges = c_prefs.ncol();
  
  // create vector of maps storing the college preferences more efficiently
  // (key -> value) <-> (student -> rank of student in college ranking)
  // this allows for O(1) access of the rank of a student in a college's ranking
  std::vector<std::unordered_map<int,int>> c_prefs_by_student(n_colleges);
  for(int c = 0; c < n_colleges; ++c) {
    for(int rank = 0; rank < n_students; ++rank) {
      if(c_prefs(rank, c) == NA_INTEGER) {
        break;
      }
      c_prefs_by_student[c][c_prefs(rank, c)] = rank;
    }
  }
  
  // vector of students that are temporarily single / definitely stay single
  std::vector<int> temp_singles, final_singles;
  
  // vector containing all students (used for partitioning later)
  std::vector<int> all_students(s_prefs.ncol());
  std::iota (std::begin(all_students), std::end(all_students), 0);
  
  // vector of proposals made by each student (zero)
  std::vector<int> n_props(n_students);
  
  // vector of assignment maps
  // for each college, a map: (key -> value)
  // <-> (rank of temporarily assigned student in college's ranking -> student)
  // this allows for fast decisions on which student to admit in the deferred
  // acceptance algorithm
  std::vector<std::map<int,int>> assignments(n_colleges);
  
  // 
  std::map<int,int> interrupting_pairs; 
  
  // number of efficiency-adjustments
  int ea_iters = -1;
  
  // number of inner algorithm (immediate/deferred acceptance) iterations
  int iter;
  
  do {
    // update number of efficiency-adjustments
    ++ea_iters;
    
    iter = 0;
    
    // vector of number of proposals that each student made so far (zero)
    std::vector<int> n_props(n_students);
    
    // reset vector of assignment maps for next efficiency adjustment iteration 
    std::vector<std::map<int,int>>(n_colleges).swap(assignments);
    
    // remove interrupters from previous round from student preferences
    for(const auto& interrupting_pair: interrupting_pairs) {
      int s = interrupting_pair.first, c = interrupting_pair.second;
      removeInterruptingPair(s_prefs, s, c);
    }
    
    // reset map of interrupting pairs
    interrupting_pairs.clear();
   
    // partition students into vector of students that are temporarily single
    // and vector of students that definitely stay single
    temp_singles.clear();
    final_singles.clear();   
    partitionStudents(all_students, n_props, s_prefs, temp_singles, final_singles);
    
    // for each college whether a student was rejected in any previous round
    std::vector<bool> is_interrupting(n_colleges, false);
    
    while(!temp_singles.empty()) {
      // update number of algorithm iterations
      ++iter;
      
      // let all temporarily single singles make an offer
      std::map<int, std::vector<int>> offers;
      for(auto s : temp_singles) {
        int next_college = s_prefs(n_props[s], s);
        offers[next_college].push_back(s);
        n_props[s]++;
      }
      
      // vector of students rejected in the current round
      std::vector<int> rejected;
      
      // map of interrupting pairs due to this deferred-acceptance round
      // interrupting student -> interrupting college
      std::map<int,int> current_interrupting_pairs;
      
      // iterate over approached colleges and decide which of the proposing
      // students to accept
      // (either using deferred-acceptance or immediate-acceptance)
      for(auto const& offer : offers) {
        int approached = offer.first;
        std::vector<int> proposers = offer.second;
     
        if(acceptance == "deferred") {
          // map of assignments that were created this DA iteration
          std::map<int,int> new_assignments; 
          
          // temporarily accept all student that are ranked by the college; the
          // remaining students are rejected
          for(auto const& p: proposers) {
            if(c_prefs_by_student[approached].count(p)) {
              int rank = c_prefs_by_student[approached][p];
              assignments[approached][rank] = p;
              new_assignments[rank] = p;
            } else {
              rejected.push_back(p);
            }
          }
        
          // will rejected students that were not newly assigned create an
          // interrupting pair?
          bool interrupting_flag = is_interrupting[approached];
         
          // as long as there are more students assigned to the college as it has
          // capacity, reject highest ranked student assigned
          while(assignments[approached].size() > n_slots[approached]) {
            std::map<int,int>::iterator it = prev(assignments[approached].end());
            rejected.push_back(it->second);
           
            is_interrupting[approached] = true;
            if(interrupting_flag && !new_assignments.count(it->first) && consent[it->second]) {
              current_interrupting_pairs[it->second] = approached;
            }
            
            assignments[approached].erase(it);
          }
        } else {
          std::map<int,int> possible_assignments;
          
          // students that are ranked by the college are considered; the remaining
          // students are immediately rejected
          for(auto const& p: proposers) {
            if(c_prefs_by_student[approached].count(p)) {
              int rank = c_prefs_by_student[approached][p];
              possible_assignments[rank] = p;
            } else {
              rejected.push_back(p);
            }
          }
         
          // transfer students from possible assignment to final assignment 
          // until there is no capacity left; 
          // the remaining students are rejected
          for(auto const& assignment: possible_assignments) {
            if(assignments[approached].size() >= n_slots[approached]) {
              rejected.push_back(assignment.second);
            } else {
              assignments[approached][assignment.first] = assignment.second;
            }
          }
        }
      }
      
      // update interrupting pairs if necessary
      if(current_interrupting_pairs.size()) {
        interrupting_pairs = current_interrupting_pairs;
      }
      
      // parition rejected students into finally / temporarily rejected students
      temp_singles.clear();
      partitionStudents(rejected, n_props, s_prefs, temp_singles, final_singles);
    }
    
  } while(interrupting_pairs.size() && ea_iters < bound_ea_iters);
  
  
  // initialize, fill, and return results list
  List results;

  // in R, results["matches"] will be a list of named double vectors
  // list index <-> college
  // entries of double vector <-> matched students
  // names of double vector <-> ranks of matched students in college ranking
  results["matchings"] = assignments;
  results["singles"] = final_singles;
  results["iters"] = iter;
  results["ea_iters"] = ea_iters;
  results["s_prefs"] = s_prefs;
  results["interrupting_pairs"] = interrupting_pairs;
  return results;
}

void partitionStudents(
    const std::vector<int>& students,
    const std::vector<int>& n_props,
    const IntegerMatrix& s_prefs,
    std::vector<int>& temp_singles,
    std::vector<int>& final_singles
) {
  int n_colleges = s_prefs.nrow();
  for(const auto& s : students) {
    if(n_props[s] < n_colleges && s_prefs(n_props[s], s) != NA_INTEGER) {
      temp_singles.push_back(s);
    } else {
      final_singles.push_back(s);
    }
  }
}

void removeInterruptingPair(
  IntegerMatrix& s_prefs, int student, int college
) {
  bool found = false;
  for(int rank = 0; rank < s_prefs.nrow(); ++rank) {
    if(s_prefs(rank, student) == college) {
      s_prefs(rank, student) = NA_INTEGER;
      found = true;
    } else {
      s_prefs(rank-found, student) = s_prefs(rank, student);
      if(s_prefs(rank, student) == NA_INTEGER) break;
    }
  }
  if(found) s_prefs(s_prefs.nrow() - 1, student) = NA_INTEGER;
}
