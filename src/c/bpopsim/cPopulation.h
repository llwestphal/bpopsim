#ifndef cPopulation_h
#define cPopulation_h

#include "cSubpopulation.h"

using namespace std;

namespace bpopsim {

  //@agm This is a remnant from cLineageTree.h
  class cLineageTree : public tree<cGenotype> {};

  class cPopulation {

  private:
    
    // Simulation parameters that should be arguments
    // Initialized outside of constructor
    double m_initial_population_size;
    double m_pop_size_after_dilution;
    uint32_t m_max_transfers_from_cli;
    uint16_t m_coarse_graining;
    double m_mutation_rate_per_division;
    double m_max_w; 
    
    //All vector types initialized at constructor
    std::vector<cSubpopulation> m_current_subpopulations;
    std::vector<uint32_t> m_divided_lineages;
    std::vector< std::vector<double> > m_runs;
    std::vector<double> m_this_run;
    std::vector<double> m_total_cells;
    std::vector< vector<double> > m_all_subpopulations_at_all_times;
    std::vector<double> m_first_mutational_vals;
    std::vector<double> m_average_fitness;
    std::vector<cSubpopulation>  m_mutations_since_last_transfer;
    std::vector< std::vector<cGenotypeFrequency> > m_frequencies;
    
    //All of the following should be initialized in the constructor
    //@agm DO NOT MAKE A CLASS VARIABLE THAT IS NOT INITIALIZED AT THE CONSTRUCTOR!!!
    uint32_t m_num_completed_transfers;   //Number of transfers that have been completed thus far
    uint32_t m_total_mutations;           //Number of mutations so far
    uint32_t m_total_subpopulations_lost; //Number of subpopulations that no longer exist
    uint16_t m_seed;                      //Random number generator
    uint64_t m_genotype_count;            //used to assign node ids in tree, should equal number of nodes
    uint64_t m_muller_rez;                //Vertical Resolution of muller plot... set at command line (default: 2500)
    
    // @JEB: An uint16_t rather than a uint because this can go negative by a few cells
    //       when cells (usually the ancestors) divide simultaneously.
    int64_t m_divisions_until_mutation; 
    
    double m_population_size;           //Current population size
    double m_average_mutation_s;           //Not sure what this is... it is a returned value.
    double m_growth_phase_generations;     //Again not sure
    double m_pop_size_before_dilution;     //Self explanatory
    double m_dilution_factor;              //The fraction by which to be diluted
    double m_transfer_binomial_sampling_p;
    double m_lambda;
    double m_by_color[2]; 
    double m_max_divergence_factor; 
    double m_binomial_sampling_threshold; 
    double m_completed_divisions;
    double m_ratio;
    double m_initial_fitness;
    double m_mutation_timer;
    double m_elapsed_time;
    double m_cell_equivalents;
    
    //Other sporadic types
    bool m_keep_transferring;

    char m_beneficial_mutation_distribution;
    
    gsl_rng* m_rng;
    
    cLineageTree m_tree;
    
  public:

    //CONSTRUCTOR  
    cPopulation() : 
    m_population_size(0),
    m_num_completed_transfers(0),
    m_total_mutations(0),
    m_total_subpopulations_lost(0),
    m_divisions_until_mutation(0),
    m_average_mutation_s(0),
    m_growth_phase_generations(0),
    m_pop_size_before_dilution(0),
    m_dilution_factor(0),
    m_transfer_binomial_sampling_p(0),
    m_lambda(0),
    m_max_divergence_factor(0),
    m_binomial_sampling_threshold(0),
    m_completed_divisions(0),
    m_ratio(0),
    m_seed(0),
    m_max_w(1),
    m_genotype_count(0),
    m_muller_rez(0),
    m_divided_lineages(0),
    m_current_subpopulations(0),
    m_runs(0),
    m_this_run(0),
    m_total_cells(0),
    m_all_subpopulations_at_all_times(0),
    m_first_mutational_vals(0),
    m_average_fitness(0),
    m_keep_transferring(true),
    m_mutations_since_last_transfer(0),
    m_frequencies(0),
    m_initial_fitness(1.0),
    m_mutation_timer(0),
    m_elapsed_time(0),
    m_cell_equivalents(0) { };
    
    //DESTRUCTOR
    virtual ~cPopulation() { };
    
    //GETTERS
    const double GetRatio() { return m_ratio; }
    
    const double GetPopulationSize() { return m_population_size; }; // don't delete @JEB
    const double   GetInitialFitness() { return m_initial_fitness; };
    // "Calculate" the population size by iterating through subpops - slow, but to check if m_population_size is correct!
    const double CalculatePopulationSize(); // don't delete this one @JEB
    
    const uint32_t GetTransfers() { return m_num_completed_transfers; }
    const uint32_t GetTotalMutations() { return m_total_mutations; }
    const uint32_t GetTotalSubpopulationsLost() { return m_total_subpopulations_lost; }
    const int32_t GetTotalTransfers() { return m_max_transfers_from_cli; }

    const int64_t GetDivisionsUntilMutation() { return m_divisions_until_mutation; }   //don't delete @JEB
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

    const double GetInitialPopulationSize() {return m_initial_population_size; }
    const double GetPopSizeAfterDilution() {return m_pop_size_after_dilution; }
    const double GetMutationRatePerDivision() {return m_mutation_rate_per_division; }
    const double GetAverageMutationS() {return m_average_mutation_s; }
    const double GetGrowthPhaseGenerations() { return m_growth_phase_generations; }
    
    const uint16_t GetSeed() { return m_seed; }

    std::vector<cSubpopulation> GetPopulation() { return m_current_subpopulations; }


    enum e_colors {
      RED=0,
      WHITE=1,
    };

    //Overloaded operators
    
    cSubpopulation& operator [] (int index) { return m_current_subpopulations[index]; }
    const cSubpopulation& operator [] (int index) const { return m_current_subpopulations[index]; }
    
    //  virtual MutationList& GetMutations(const char in_marker) {};

    //SETTERS
    void SetTotalMutations(uint32_t in_total_mutations) { m_total_mutations = in_total_mutations; }
    void SetTotalSubpopulationsLost(uint32_t in_total_subpopulations_lost) { m_total_subpopulations_lost=in_total_subpopulations_lost; }
    void SetTransfers(uint32_t in_transfers) { m_num_completed_transfers = in_transfers; }
    void SetDivisionsUntilMutation(int64_t in_divisions_until_mutation){ m_divisions_until_mutation = in_divisions_until_mutation; }  //!@JEB - keep
    void SetCompletedDivisions(uint32_t in_completed_divisions) {m_completed_divisions = in_completed_divisions; }
    void SetMaxW(double in_max_w) { m_max_w = in_max_w; }
    void SetPopSizeBeforeDilution(double in_pop_size_before_dilution) { m_pop_size_before_dilution = in_pop_size_before_dilution; }
    void SetDilutionFactor(double in_dilution_factor) { m_dilution_factor = in_dilution_factor; }
    void SetTransferBinomialSamplingP(double in_transfer_binomial_sampling_p) { m_transfer_binomial_sampling_p = in_transfer_binomial_sampling_p; }
    void SetLambda(double in_lambda) { m_lambda = in_lambda; }
    void SetRatio(double in_ratio) { m_ratio = in_ratio; }
    void SetTotalTransfers(uint32_t in_total_transfers) { m_max_transfers_from_cli = in_total_transfers; }
    void SetMaxDivergenceFactor( double in_max_divergence_factor) { m_max_divergence_factor = in_max_divergence_factor; }
    void SetBinomialSamplingThreshold (double in_binomial_sampling_threshold) { m_binomial_sampling_threshold = in_binomial_sampling_threshold; }
    void SetInitialPopulationSize(double in_initial_population_size) {m_initial_population_size = in_initial_population_size; }
    void SetPopSizeAfterDilution(double in_pop_size_after_dilution) {m_pop_size_after_dilution = in_pop_size_after_dilution; }
    void SetMutationRatePerDivision(double in_mutation_rate_per_division) {m_mutation_rate_per_division = in_mutation_rate_per_division; }
    void SetGrowthPhaseGenerations(double in_growth_phase_generations) { m_growth_phase_generations= in_growth_phase_generations; }
    void SetBeneficialMutationDistribution(char in_beneficial_mutation_distribution) { m_beneficial_mutation_distribution = in_beneficial_mutation_distribution; }
    void SetSeedParams(uint16_t in_seed_type) { m_seed = in_seed_type; }
    void SetInitialMutVals(std::vector<double> in_mut_vals) { m_first_mutational_vals = in_mut_vals; }
    void SetCoarseGraining(uint16_t in_coarse) { m_coarse_graining = in_coarse; }
    void SetInitialFitness(double in_initial_fitness ) { m_initial_fitness = in_initial_fitness; }
    void SetMullerRez(uint64_t in_muller_rez ) { m_muller_rez = in_muller_rez; }
    
    void SetRNG(gsl_rng * in_rng) { m_rng = in_rng; } //@JEB

    //METHODS
    //@agm To keep the lines of manageable length, if a method has multiple variables, each variable got a new line
      
    void   SetParameters(AnyOption & options);
    //! Move time forward by this increment, growing all subpopulations
    void   UpdateSubpopulations(double update_time);
    //! Calculate the time until the next subpopulation divides (passes a whole number of cells)
    double TimeToNextWholeCell();
    void   FrequenciesPerTransferPerNode(); 
    void   ConvertExternalData(const string &input_file);
    void   PrintFreqsQuick();
    std::vector<cGenotypeFrequency>::iterator Find_Node_in_Freq(std::vector<cGenotypeFrequency> &frequencies, 
                                                                tree<cGenotype>::sibling_iterator this_node);
    double Find_Node_in_Freq_By_NodeID(std::vector<cGenotypeFrequency> &frequencies,
                                       uint32_t this_node);
    cSubpopulation* Find_Node_in_Populations_By_NodeID(uint32_t this_node);
    
    double AssignChildFreq(tree<cGenotype>::sibling_iterator child_node,
                           double parent_low,
                           double parent_high,
                           std::vector<cFrequencySlice> * child_freqs, 
                           std::vector<cGenotypeFrequency> &frequencies,
                           int depth = 0);
    void   DrawMullerMatrix(std::string output_folder,
                          std::vector< std::vector<int> > muller_matrix);
    void   Resample();
    void   Deterministic_Resample();
    void   CullPopulations();
    void   PushBackRuns();
    void   RunSummary();
    void   DisplayParameters();
    void   CalculateDivisions();
    void   CalculateDivisionsNew();
    void   SeedSubpopulationForRedWhite();
    void   SeedPopulationWithOneColony();
    void   AddSubpopulation(cSubpopulation& subpop);
    void   Mutate();
    void   MutateNew();
    std::vector<uint16_t> CurrentUniqueGenotypes();
    void   PrintUniqueGenotypes(const std::string& output_folder,
                              std::vector< std::vector<uint16_t> > * number_of_unique_genotypes);
    void   PrintOut(const std::string& output_folder, uint32_t on_run);
    void   PrintExpectationValue(const std::string& output_folder);
    void   PrintOut_RedWhiteOnly(const std::string& output_folder,
                               std::vector< std::vector<double> > *red_white_ratios,
                               uint16_t transfer_interval_to_print);
    void   PrintFrequenciesToScreen(std::string output_folder);
    void   PrintFrequenciesToScreen_RedWhiteOnly(std::string output_folder);
    void   PrintSingleFitness(std::string out_folder);
    
    uint32_t MutationAboveThreshold_2(float threshold);
    std::vector<uint32_t> MutationAboveThreshold(float threshold);
    double CalculateSimilarity(std::string output_folder);
    double CountMutipleDivergedSubpops();
    void   TimeToSweep(std::string output_folder);
    void   CalculateAverageFitness();
    void   PrintFitness(std::string output_folder);
    float  Logarithm(float mantissa);
    // Prints out the tree using bracket notation.
    void   PrintTree();
  };

}

#endif
