#ifndef cPopulation_h
#define cPopulation_h

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "cSubpopulation.h"

//class cSubpopulation;

#include "lineageTree.h"


// Boost
#include <boost/program_options.hpp>

using namespace boost::program_options;

class cPopulation {

private:

  long double m_total_pop_size;
  long double m_new_pop_size;
  long double m_ratio;

  int m_total_mutations;
  int m_total_subpopulations_lost;
  int m_transfers;
  int m_number_of_subpopulations;  
  int m_verbose;
  int m_replicates;
  int m_minimum_printed;
  int m_total_transfers;
  int m_transfer_interval_to_print;
  int m_lineage;

  std::vector<cSubpopulation> m_populations;
  std::vector<int> m_divided_lineages;
  std::vector< std::vector<double> > m_runs;
  std::vector<double> m_this_run;

  double m_divisions_until_mutation;
  double m_desired_divisions;
  double m_completed_divisions;
  double m_update_time;
  double m_current_cells;
  double m_whole_cells;
  double m_this_time_to_next_whole_cell;
  double m_time_to_next_whole_cell;
  double m_max_w;
  double m_pop_size_before_dilution;
  double m_dilution_factor;
  double m_transfer_binomial_sampling_p;
  double m_lambda;
  double m_by_color[2];
  double m_max_divergence_factor;
  double m_binomial_sampling_threshold;

  bool m_keep_transferring;

  char m_beneficial_mutation_distribution;

  // Simulation parameters that should be arguments
  uint64_t m_initial_population_size;
  uint64_t m_pop_size_after_dilution;             // N sub 0 --int is to get rid of warning
  double m_mutation_rate_per_division;            // mu
  double m_average_mutation_s;                    // s
  double m_growth_phase_generations;
  
public:

  //CONSTRUCTOR  
  cPopulation()
  {
     m_transfers = 0;
     m_verbose = 0;
     m_max_w = 1;
     m_number_of_subpopulations = 0;
     m_total_mutations = 0;
     m_total_subpopulations_lost = 0;
  }
  
  //DESTRUCTOR
  virtual ~cPopulation() { ; };
  
  //GETTERS
  const long double GetTotalPopSize() { return m_total_pop_size; }
  const long double GetRatio() { return m_ratio; }
  
  const int GetTotalMutations() { return m_total_mutations; }
  const int GetTotalSubpopulationsLost() { return m_total_subpopulations_lost; }
  const int GetTransfers() { return m_transfers; }
  const int GetNumberOfSubpopulations() { return m_number_of_subpopulations; }
  const int GetNewPopSize() { return m_new_pop_size; }
  const int GetVerbose() { return m_verbose; }
  const int GetTransferIntervalToPrint() {return m_transfer_interval_to_print; }
  const int GetTotalTransfers() {return m_total_transfers; }
  const int GetReplicates() {return m_replicates; }
  const int GetMinimumPrinted() {return m_minimum_printed; }
  const int GetLineageTree() {return m_lineage; }
  //std::vector< std::vector<double> > GetRuns() { return &m_runs; }
  //std::vector<double> GetThisRun() { return &m_this_run; }  

  const double GetDivisionsUntilMutation() { return m_divisions_until_mutation; }
  const double GetDesiredDivisions() { return m_desired_divisions; }
  const double GetCompletedDivisions() { return m_completed_divisions; }
  const double GetUpdateTime() { return m_update_time; }
  const double GetCurrentCells() { return m_current_cells; }
  const double GetWholeCells() { return m_whole_cells; }
  const double GetThisTimeToNextWholeCell() { return m_this_time_to_next_whole_cell; }
  const double GetTimeToNextWholeCell() { return m_time_to_next_whole_cell; }
  const double GetMaxW() { return m_max_w; }
  const double GetPopSizeBeforeDilution() { return m_pop_size_before_dilution; }
  const double GetDilutionFactor() { return m_dilution_factor; }
  const double GetTransferBinomialSamplingP() { return m_transfer_binomial_sampling_p; }
  const double GetLambda() {return m_lambda; }
  const double GetMaxDivergenceFactor() {return m_max_divergence_factor; }
  const double GetBinomialSamplingThreshold() { return m_binomial_sampling_threshold; }
  
  const bool GetKeepTransferring() {return m_keep_transferring; }
  const char GetBeneficialMutationDistribution() {return m_beneficial_mutation_distribution; }  

  const uint64_t GetInitialPopulationSize() {return m_initial_population_size; }
  const uint64_t GetPopSizeAfterDilution() {return m_pop_size_after_dilution; }
  const double GetMutationRatePerDivision() {return m_mutation_rate_per_division; }
  const double GetAverageMutationS() {return m_average_mutation_s; }
  const double GetGrowthPhaseGenerations() { return m_growth_phase_generations; }

  std::vector<cSubpopulation> GetPopulation() { return m_populations; }


  enum e_colors {
    RED=0,
    WHITE=1,
  };

  //Overloaded operators
  
  cSubpopulation& operator [] (int index) { return m_populations[index]; }
  const cSubpopulation& operator [] (int index) const { return m_populations[index]; }
  
  //  virtual MutationList& GetMutations(const char in_marker) {};

  //SETTERS
  void SetTotalPopSize(long double in_total_pop_size) { m_total_pop_size = in_total_pop_size; }
  void SetTotalMutations(int in_total_mutations) { m_total_mutations = in_total_mutations; }
  void SetTotalSubpopulationsLost(int in_total_subpopulations_lost) { m_total_subpopulations_lost=in_total_subpopulations_lost; }
  void SetTransfers(int in_transfers) { m_transfers = in_transfers; }
  void SetDivisionsUntilMutation(double in_divisions_until_mutation){ m_divisions_until_mutation = in_divisions_until_mutation; }  
  void SetDesiredDivisions(double in_desired_divisions){ m_desired_divisions = in_desired_divisions; }
  void SetNumberOfSubpopulations(int in_number_of_subpopulations){ m_number_of_subpopulations = in_number_of_subpopulations; }
  void SetNewPopSize(int in_new_pop_size) {m_new_pop_size = in_new_pop_size; }
  void SetCompletedDivisions(int in_completed_divisions) {m_completed_divisions = in_completed_divisions; }
  void SetUpdateTime(double in_update_time) {m_update_time = in_update_time; }
  void SetCurrentCells(double in_current_cells) { m_current_cells = in_current_cells; }
  void SetWholeCells(double in_whole_cells) { m_whole_cells = in_whole_cells; }
  void SetThisTimeToNextWholeCell(double in_this_time_to_next_whole_cell) { m_this_time_to_next_whole_cell = in_this_time_to_next_whole_cell; }
  void SetTimeToNextWholeCell(double in_time_to_next_whole_cell) { m_time_to_next_whole_cell = in_time_to_next_whole_cell; }
  void SetMaxW(double in_max_w) { m_max_w = in_max_w; }
  void SetVerbose(int in_verbose) { m_verbose = in_verbose; }
  void SetPopSizeBeforeDilution(double in_pop_size_before_dilution) { m_pop_size_before_dilution = in_pop_size_before_dilution; }
  void SetDilutionFactor(double in_dilution_factor) { m_dilution_factor = in_dilution_factor; }
  void SetTransferBinomialSamplingP(double in_transfer_binomial_sampling_p) { m_transfer_binomial_sampling_p = in_transfer_binomial_sampling_p; }
  void SetLambda(double in_lambda) { m_lambda = in_lambda; }
  void SetKeepTransferring(bool in_keep_transferring) { m_keep_transferring = in_keep_transferring; }
  void SetRatio(long double in_ratio) { m_ratio = in_ratio; }
  void SetTransferIntervalToPrint(int in_transfer_interval_to_print) { m_transfer_interval_to_print = in_transfer_interval_to_print; }
  void SetTotalTransfers(int in_total_transfers) { m_total_transfers = in_total_transfers; }
  void SetMaxDivergenceFactor( double in_max_divergence_factor) { m_max_divergence_factor = in_max_divergence_factor; }
  void SetReplicates ( int in_replicates) { m_replicates = in_replicates; }
  void SetMinimumPrinted (int in_minimum_printed) { m_minimum_printed = in_minimum_printed; }
  void SetBinomialSamplingThreshold (double in_binomial_sampling_threshold) { m_binomial_sampling_threshold = in_binomial_sampling_threshold; }
  void SetInitialPopulationSize(uint64_t in_initial_population_size) {m_initial_population_size =in_initial_population_size; }
  void SetPopSizeAfterDilution(uint64_t in_pop_size_after_dilution) {m_pop_size_after_dilution = in_pop_size_after_dilution; }
  void SetMutationRatePerDivision(double in_mutation_rate_per_division) {m_mutation_rate_per_division = in_mutation_rate_per_division; }
  void SetAverageMutationS(double in_average_mutation_s) {m_average_mutation_s = in_average_mutation_s; }
  void SetGrowthPhaseGenerations(double in_growth_phase_generations) { m_growth_phase_generations= in_growth_phase_generations; }
  void SetBeneficialMutationDistribution(char in_beneficial_mutation_distribution) { m_beneficial_mutation_distribution = in_beneficial_mutation_distribution; }
  void SetLineageTree(int in_lineage) { m_lineage = in_lineage; }

  //METHODS
  void AddSubpopulation(cSubpopulation& subpop);
  void UpdateLineages();
  void DetermineDivisionTime();
  void Mutate(gsl_rng * randomgenerator, lineageTree& red, lineageTree& white);
  void Resample(gsl_rng * randomgenerator);
  void PushBackRuns();
  void PrintOut(const std::string& output_file_name);
  void ClearRuns(lineageTree& red, lineageTree& white);
  void RunSummary();
  void ResetRunStats();
  void SetParameters(const variables_map &options);
  void DisplayParameters();
  void CalculateDivisions();
  void SeedSubpopulations(lineageTree& red, lineageTree& white);
};


#endif
