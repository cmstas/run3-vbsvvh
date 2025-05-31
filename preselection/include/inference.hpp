// RInferenceUtils.hxx modified to use RVecs instead of just primitive types
#ifndef TMVA_RINFERENCEUTILS
#define TMVA_RINFERENCEUTILS
 
#include <utility> // std::forward, std::index_sequence
#include <vector>  // std::vector
#include <ROOT/RVec.hxx> // ROOT::RVec
 
namespace TMVA {
namespace Experimental {
 
namespace Internal {
 
/// Vector Compute helper
template <typename I, typename T, typename F>
class ComputeHelper;
 
template <std::size_t... N, typename T, typename F>
class ComputeHelper<std::index_sequence<N...>, T, F> {
   template <std::size_t Idx>
   using AlwaysVectorT = ROOT::RVec<T>;
   F fFunc;
 
public:
   ComputeHelper(F &&f) : fFunc(std::forward<F>(f)) {}
   
   // Vector version that processes each element and returns a vector of results
   ROOT::RVec<T> operator()(AlwaysVectorT<N>... args) {
      // Get the size of the first vector (assuming all vectors have the same size)
      const size_t size = std::get<0>(std::make_tuple(args...)).size();
      ROOT::RVec<T> results;
      
      // Process each set of elements
      for (size_t i = 0; i < size; ++i) {
         // Extract the ith element from each vector and compute result
         // Compute returns a vector, we need to extract first element
         auto result_vec = fFunc.Compute({(args[i])...});
         if (!result_vec.empty()) {
            results.push_back(result_vec[0]);
         }
      }
      
      return results;
   }
};
 
} // namespace Internal
 
/// Helper to pass TMVA model to RDataFrame.Define nodes
template <std::size_t N, typename T, typename F>
auto Compute(F &&f) -> Internal::ComputeHelper<std::make_index_sequence<N>, T, F>
{
   return Internal::ComputeHelper<std::make_index_sequence<N>, T, F>(std::forward<F>(f));
}
 
} // namespace Experimental
} // namespace TMVA
 
#endif // TMVA_RINFERENCEUTILS