#ifndef cPopulation_h
#define cPopulation_h

#include "cSubpopulation.h"

using namespace boost::program_options;

//@agm This is a remnant from cLineageTree.h
class cLineageTree : public tree<cGenotype> {};

class cPopulation {

private:

  double m_ratio;
  
  // we calculate the population size on demand
  uint32_t m_population_size;
  uint32_t m_transfers;
  uint32_t m_total_mutations;
  uint32_t m_total_subpopulations_lost;
  uint32_t m_number_of_subpopulations;
  uint32_t m_total_transfers;
    
  bool m_population_size_stale;

  uint16_t m_verbose;
  uint16_t m_replicates;
  uint16_t m_minimum_printed;
  uint16_t m_transfer_interval_to_print;
  uint16_t m_lineage;

  std::vector<cSubpopulation> m_populations;
  std::vector<uint32_t> m_divided_lineages;
  std::vector< std::vector<double> > m_runs;
  std::vector<double> m_this_run;

  // @JEB: An uint16_t rather than a uint because this can go negative by a few cells
  //       when cells (usually the ancestors) divide simultaneously.
  int64_t m_divisions_until_mutation; 
  
  double m_desired_divisions;
  double m_completed_divisions;
  double m_update_time;
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
  char m_red_white_only;

  // Simulation parameters that should be arguments
  uint32_t m_initial_population_size;
  uint32_t m_pop_size_after_dilution;             // N sub 0 --int is to get rid of warning
  uint16_t m_seed;
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
     m_population_size_stale = true;
  }
	
  //DESTRUCTOR
  virtual ~cPopulation() { ; };
  
  //GETTERS
  const double GetRatio() { return m_ratio; }
  
  const uint32_t GetPopulationSize(); // calculated on demand
  const uint32_t GetTransfers() { return m_transfers; }
  const uint32_t GetTotalMutations() { return m_total_mutations; }
  const uint32_t GetTotalSubpopulationsLost() { return m_total_subpopulations_lost; }
  const uint32_t GetNumberOfSubpopulations() { return m_number_of_subpopulations; }
  const int32_t GetTotalTransfers() { return m_total_transfers; }
  
  const uint16_t GetVerbose() { return m_verbose; }
  const uint16_t GetTransferIntervalToPrint() { return m_transfer_interval_to_print; }
  const uint16_t GetReplicates() { return m_replicates; }
  const uint16_t GetMinimumPrinted() { return m_minimum_printed; }
  const uint16_t GetLineageTree() { return m_lineage; }

  const int64_t GetDivisionsUntilMutation() { return m_divisions_until_mutation; }   //!@JEB - keep
  const double GetCompletedDivisions() { return m_completed_divisions; }
  const double GetMaxW() { return m_max_w; }
  const double GetPopSizeBeforeDilution() { return m_pop_size_before_dilution; }
  const double GetDilutionFactor() { return m_dilution_factor; }
  const double GetTransferBinomialSamplingP() { return m_transfer_binomial_sampling_p; }
  const double GetLambda() {return m_lambda; }
  const double GetMaxDivergenceFactor() { return m_max_divergence_factor; }
  const double GetBinomialSamplingThreshold() { return m_binomial_sampling_threshold; }
  
  const bool GetKeepTransferring() { return m_keep_transferring; }
  const char GetBeneficialMutationDistribution() { return m_beneficial_mutation_distribution; } 
  const char GetRedWhiteOnly() { return m_red_white_only; }

  const uint32_t GetInitialPopulationSize() {return m_initial_population_size; }
  const uint32_t GetPopSizeAfterDilution() {return m_pop_size_after_dilution; }
  const double GetMutationRatePerDivision() {return m_mutation_rate_per_division; }
  const double GetAverageMutationS() {return m_average_mutation_s; }
  const double GetGrowthPhaseGenerations() { return m_growth_phase_generations; }
	
  const uint16_t GetSeed() { return m_seed; }

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
  void SetTotalMutations(uint32_t in_total_mutations) { m_total_mutations = in_total_mutations; }
  void SetTotalSubpopulationsLost(uint32_t in_total_subpopulations_lost) { m_total_subpopulations_lost=in_total_subpopulations_lost; }
  void SetTransfers(uint32_t in_transfers) { m_transfers = in_transfers; }
  void SetDivisionsUntilMutation(int64_t in_divisions_until_mutation){ m_divisions_until_mutation = in_divisions_until_mutation; }  //!@JEB - keep
  void SetNumberOfSubpopulations(uint32_t in_number_of_subpopulations){ m_number_of_subpopulations = in_number_of_subpopulations; }
  void SetCompletedDivisions(uint32_t in_completed_divisions) {m_completed_divisions = in_completed_divisions; }
  void SetMaxW(double in_max_w) { m_max_w = in_max_w; }
  void SetVerbose(int in_verbose) { m_verbose = in_verbose; }
  void SetPopSizeBeforeDilution(double in_pop_size_before_dilution) { m_pop_size_before_dilution = in_pop_size_before_dilution; }
  void SetDilutionFactor(double in_dilution_factor) { m_dilution_factor = in_dilution_factor; }
  void SetTransferBinomialSamplingP(double in_transfer_binomial_sampling_p) { m_transfer_binomial_sampling_p = in_transfer_binomial_sampling_p; }
  void SetLambda(double in_lambda) { m_lambda = in_lambda; }
  void SetKeepTransferring(bool in_keep_transferring) { m_keep_transferring = in_keep_transferring; }
  void SetRedWhiteOnly(char in_red_white_only) { m_red_white_only = in_red_white_only; }
  void SetRatio(double in_ratio) { m_ratio = in_ratio; }
  void SetTransferIntervalToPrint(int in_transfer_interval_to_print) { m_transfer_interval_to_print = in_transfer_interval_to_print; }
  void SetTotalTransfers(uint32_t in_total_transfers) { m_total_transfers = in_total_transfers; }
  void SetMaxDivergenceFactor( double in_max_divergence_factor) { m_max_divergence_factor = in_max_divergence_factor; }
  void SetReplicates ( uint16_t in_replicates) { m_replicates = in_replicates; }
  void SetMinimumPrinted (int in_minimum_printed) { m_minimum_printed = in_minimum_printed; }
  void SetBinomialSamplingThreshold (double in_binomial_sampling_threshold) { m_binomial_sampling_threshold = in_binomial_sampling_threshold; }
  void SetInitialPopulationSize(uint32_t in_initial_population_size) {m_initial_population_size =in_initial_population_size; }
  void SetPopSizeAfterDilution(uint32_t in_pop_size_after_dilution) {m_pop_size_after_dilution = in_pop_size_after_dilution; }
  void SetMutationRatePerDivision(double in_mutation_rate_per_division) {m_mutation_rate_per_division = in_mutation_rate_per_division; }
  void SetAverageMutationS(double in_average_mutation_s) {m_average_mutation_s = in_average_mutation_s; }
  void SetGrowthPhaseGenerations(double in_growth_phase_generations) { m_growth_phase_generations= in_growth_phase_generations; }
  void SetBeneficialMutationDistribution(char in_beneficial_mutation_distribution) { m_beneficial_mutation_distribution = in_beneficial_mutation_distribution; }
  void SetLineageTree(int in_lineage) { m_lineage = in_lineage; }
  void SetSeedParams(uint16_t seed_type) { m_seed = seed_type; }

  //METHODS
	//@agm To keep the lines of manageable length, if a method has multiple variables, each variable got a new line
    
  void SetParameters(const variables_map &options);
	
	//! Move time forward by this increment, growing all subpopulations
  void UpdateSubpopulations(double update_time);

  //! Calculate the time until the next subpopulation divides (passes a whole number of cells)
  double TimeToNextWholeCell();
  void FrequenciesPerTransferPerNode(tree<cGenotype> *newtree, 
                                     std::vector< std::vector<cGenotypeFrequency> >& frequencies);
  void FrequenciesOfSubpops(tree<cGenotype> newtree,
                            std::vector< std::vector<cGenotypeFrequency> >& freqs_for_muller);
  void Resample(gsl_rng * randomgenerator);
  void PushBackRuns();
  void PrintOut(const std::string& output_file_name,
                std::vector< std::vector<cGenotypeFrequency> > frequencies);
  void ClearRuns(cLineageTree* tree);
  void RunSummary();
  void ResetRunStats();
  void DisplayParameters();
  void CalculateDivisions();
  void SeedSubpopulationForRedWhite(cLineageTree* newtree, 
                                    uint32_t& node_id);
  void SeedPopulationWithOneColony(cLineageTree* newtree,
                                   uint32_t& node_id);
  void AddSubpopulation(cSubpopulation& subpop,
                        uint32_t& node_id);
  void Mutate(gsl_rng * randomgenerator, 
              cLineageTree* newtree, 
              uint32_t& node_id);
  void PrintFrequenciesToScreen(std::vector< std::vector<cGenotypeFrequency> > frequencies);
  float Logarithm(float mantissa);
  float fast_log(float val);
};

#endif
