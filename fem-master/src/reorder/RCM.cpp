// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "RCM.h"

#include "OmpInterface.h"

#include <algorithm>
#include <queue>

namespace gmshfem::reorder
{


  RCM::RCM() {}

  RCM::~RCM() {}

  struct SortClass {
    const unsigned long long *_degree;
    bool operator()(const unsigned long long i, const unsigned long long j)
    {
      return _degree[i] < _degree[j];
    }
  };

  void RCM::apply(std::vector< unsigned long long > &sorted, const unsigned long long *const row, const unsigned long long *const indices)
  {
    std::vector< unsigned long long > min; // long long for padding (64bits)
    std::vector< unsigned long long > Nmin;
    std::vector< unsigned long long > degree;
    std::vector< unsigned long long > newIndices;
#pragma omp parallel num_threads(omp::getNumThreads())
    {
      unsigned int numThreads = omp::getNumThreads();
      unsigned int myThreadID = omp::getThreadNum();
#pragma omp single
      {
        min.resize(numThreads, row[sorted.size()] - row[0]);
        Nmin.resize(numThreads, 0);
        degree.resize(sorted.size(), 0);
      }
#pragma omp for
      for(auto i = 0ULL; i < sorted.size(); ++i) {
        degree[i] = row[i + 1] - row[i];
        if(min[myThreadID] > degree[i]) {
          min[myThreadID] = degree[i];
          Nmin[myThreadID] = i;
        }
      }
#pragma omp single
      {
        for(auto i = 0U; i < numThreads; ++i) {
          if(min[0] > min[i]) {
            min[0] = min[i];
            Nmin[0] = Nmin[i];
          }
        }
      }
      SortClass sortClass;
      sortClass._degree = degree.data();

#pragma omp single
      {
        newIndices.resize(row[sorted.size()], 0);
      }
#pragma omp for
      for(auto i = 0ULL; i < sorted.size(); ++i) {
        for(auto j = row[i]; j < row[i + 1]; ++j) {
          newIndices[j] = indices[j];
        }
        std::sort(&newIndices[row[i]], &newIndices[row[i + 1]], sortClass);
      }
    }

    sorted[Nmin[0]] = sorted.size();
    unsigned long long nbrNodesSorted = sorted.size() - 1;

    std::queue< unsigned long long > queue;
    for(auto i = row[Nmin[0]]; i < row[Nmin[0] + 1]; ++i) {
      queue.push(newIndices[i]);
    }

    while(nbrNodesSorted != 0) {
      if(queue.size() != 0) {
        if(sorted[queue.front()] == 0) {
          sorted[queue.front()] = nbrNodesSorted--;
          for(auto i = row[queue.front()]; i < row[queue.front() + 1]; ++i) {
            queue.push(newIndices[i]);
          }
        }
        queue.pop();
      }
      else {
        unsigned long long min = row[sorted.size()] - row[0];
        unsigned long long Nmin = 0;
        for(auto i = 0ULL; i < sorted.size(); ++i) {
          if(min > degree[i] && sorted[i] == 0) {
            min = degree[i];
            Nmin = i;
          }
        }
        queue.push(Nmin);
      }
    }
  }


} // namespace gmshfem::reorder
