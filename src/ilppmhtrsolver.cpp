/*
 * ilppmhtrsolver.cpp
 *
 *  Created on: 20-aug-2017
 *      Author: M. El-Kebir
 */

#include "ilppmhtrsolver.h"
#include <lemon/time_measure.h>
#include <queue>

IlpPmhTrSolver::IlpPmhTrSolver(const CloneTree& T,
                               const std::string& primary,
                               MigrationGraph::Pattern pattern,
                               const std::string& gurobiLogFilename,
                               const StringPairList& forcedComigrations)
  : IlpPmhSolver(T, primary, pattern, gurobiLogFilename, forcedComigrations)
  , _pNodeToStateSet(NULL)
  , _pNodeToRootState(NULL)
  , _pTprime(NULL)
  , _zz()
  , _r()
  , _pCallback(NULL)
{
}

IlpPmhTrSolver::~IlpPmhTrSolver()
{
  if (_pNodeToStateSet != NULL){
    for (int i=0;i<anrSolutions; i++){
      delete _pNodeToStateSet[i];
    }
    delete _pNodeToStateSet;
  }
  const int anrSolutions = _model.get(GRB_IntAttr_SolCount);
  if (_pNodeToRootState != NULL){
    for (int i=0;i<anrSolutions; i++){
      delete _pNodeToRootState[i];
    }
    delete _pNodeToRootState;
  }
  delete _pTprime;
  delete _pCallback;
}

void IlpPmhTrSolver::initMultiSourceSeedingConstraints()
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    const int i = (*_pNodeToIndex)[v_i];
    for (NodeIt v_j(getTree()); v_j != lemon::INVALID; ++v_j)
    {
      const int j = (*_pNodeToIndex)[v_j];
      if (v_i == v_j) continue;
      if (isAncestor(v_i, v_j))
      {
        for (int s = 0; s < nrAnatomicalSites; ++s)
        {
          const int size_L_s = _L[s].size();
          for (int c = 0; c < size_L_s; ++c)
          {
            for (int d = 0; d < size_L_s; ++d)
            {
              if (c == d) continue;
              _model.addConstr(_r[i][s][c] + _x[j][s][d] <= 1);
            }
          }
        }
      }
    }
  }
}

void IlpPmhTrSolver::processSolution(bool post_processing)
{
  const int anrSolutions = _model.get(GRB_IntAttr_SolCount);
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  const int nrNodes = _indexToNode.size();

  // std::cout << "............." << std::endl;

  //   for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i){
  //   if (lemon::countOutArcs(getTree(), v_i) > 0){
  //     const int i = (*_pNodeToIndex)[v_i];
  //     for (int s = 0; s < nrAnatomicalSites; ++s){
  //       const int size_L_s = _L[s].size();
  //       for (int c = 0; c < size_L_s; ++c){
  //         for (int t = 0; t < nrAnatomicalSites; ++t){
  //           const int size_L_t = _L[t].size();
  //           for (int d = 0; d < size_L_t; ++d){
  //             std::cout<< _zz[i][s][c][t][d].get(GRB_StringAttr_VarName) << _zz[i][s][c][t][d].get(GRB_DoubleAttr_Xn) << " -> ";
  //             for(OutArcIt m(getTree(), v_i); m != lemon::INVALID; ++m){
  //               const int ij = (*_pArcToIndex)[m];
  //               std::cout << _xx[ij][s][c][t][d].get(GRB_StringAttr_VarName) << _xx[ij][s][c][t][d].get(GRB_DoubleAttr_Xn) <<" ";
  //             }
  //             std::cout << std::endl;
  //           }
  //         }
  //       }
  //     }
  //   }
  // }
  
  /**_model.set(GRB_IntParam_SolutionNumber, 2);
  for(int s = 0; s < nrAnatomicalSites; ++s){
    int size_L_s = _L[s].size();
    for (int c = 0; c < size_L_s; ++c){
      for(int t = 0; t < nrAnatomicalSites; ++t){
        int size_L_t = _L[t].size();
        for (int d = 0; d < size_L_t; ++d){
          std::cout << s<<' '<<c<<' '<<t<<' '<<d << ' ' << _z[s][c][t][d].get(GRB_DoubleAttr_Xn)<<'\n';
        }
      }
    }
  }

  _model.set(GRB_IntParam_SolutionNumber, 3);
  for(int s = 0; s < nrAnatomicalSites; ++s){
    int size_L_s = _L[s].size();
    for (int c = 0; c < size_L_s; ++c){
      for(int t = 0; t < nrAnatomicalSites; ++t){
        int size_L_t = _L[t].size();
        for (int d = 0; d < size_L_t; ++d){
          std::cout << s<<' '<<c<<' '<<t<<' '<<d << ' ' << _z[s][c][t][d].get(GRB_DoubleAttr_Xn)<<'\n';
        }
      }
    }
  }
  std::cout << nrAnatomicalSites<<' '<<nrNodes<<std::endl;
  for (int i=0;i<nrNodes;i++){
    for(int s=0;s<nrAnatomicalSites;s++){
      int size_L_s = _L[s].size();
      for (int c=0;c<size_L_s;c++){
        _model.set(GRB_IntParam_SolutionNumber, 2);
        int p= _x[i][s][c].get(GRB_DoubleAttr_Xn);
        int r= _r[i][s][c].get(GRB_DoubleAttr_Xn);
        _model.set(GRB_IntParam_SolutionNumber, 3);
        //if (p!= _x[i][s][c].get(GRB_DoubleAttr_Xn) || r!= _r[i][s][c].get(GRB_DoubleAttr_Xn)){
          std::cout <<i<< isLeaf(_indexToNode[i]) << ' '<<s<< ' '<<c<< ' '<<p<<' '<< _x[i][s][c].get(GRB_DoubleAttr_Xn)<<
          ' '<<r<<' '<< _r[i][s][c].get(GRB_DoubleAttr_Xn)<<std::endl;
        //}
      }
    }
  }**/
  /**for (int i=0;i<nrNodes;i++){
    std::cout<<isLeaf(_indexToNode[i])<<'\n';
  }for (int s=0;s<nrAnatomicalSites;s++){
    std::cout<<s<<'s';
    for (Node c:_L[s]){
        std::cout << (*_pNodeToIndex)[c]<<' ';
      }std::cout <<'\n';
  }**/

//  const int nrArcs = _indexToArc.size();

//  for (int i = 0; i < nrNodes; ++i)
//  {
//    for (int s = 0; s < nrAnatomicalSites; ++s)
//    {
//      int size_L_s = _L[s].size();
//      for (int c = 0; c < size_L_s; ++c)
//      {
//        if (_r[i][s][c].get(GRB_DoubleAttr_X) >= 0.4)
//        {
//          std::cout << _r[i][s][c].get(GRB_StringAttr_VarName)
//          << " = " << _r[i][s][c].get(GRB_DoubleAttr_X) << std::endl;
//        }
//      }
//    }
//  }
//  
//  for (int i = 0; i < nrNodes; ++i)
//  {
//    Node v_i = _indexToNode[i];
//    for (int s = 0; s < nrAnatomicalSites; ++s)
//    {
//      int size_L_s = _L[s].size();
//      for (int c = 0; c < size_L_s; ++c)
//      {
//        if (_x[i][s][c].get(GRB_DoubleAttr_X) >= 0.4)
//        {
//          std::cout << _x[i][s][c].get(GRB_StringAttr_VarName)
//          << " = " << _x[i][s][c].get(GRB_DoubleAttr_X) << std::endl;
//        }
//      }
//    }
//  }
//  
//  for (int i = 0; i < nrNodes; ++i)
//  {
//    for (int s = 0; s < nrAnatomicalSites; ++s)
//    {
//      int size_L_s = _L[s].size();
//      for (int c = 0; c < size_L_s; ++c)
//      {
//        for (int t = 0; t < nrAnatomicalSites; ++t)
//        {
//          int size_L_t = _L[t].size();
//          for (int d = 0; d < size_L_t; ++d)
//          {
//            if (_zz[i][s][c][t][d].get(GRB_DoubleAttr_X) >= 0.4)
//            {
//              std::cout << _zz[i][s][c][t][d].get(GRB_StringAttr_VarName)
//              << " = " << _zz[i][s][c][t][d].get(GRB_DoubleAttr_X) << std::endl;
//            }
//          }
//        }
//      }
//    }
//  }
  
  //    for (int s = 0; s < nrAnatomicalSites; ++s)
  //    {
  //      int size_L_s = _L[s].size();
  //      for (int c = 0; c < size_L_s; ++c)
  //      {
  //        for (int t = 0; t < nrAnatomicalSites; ++t)
  //        {
  //          int size_L_t = _L[t].size();
  //          for (int d = 0; d < size_L_t; ++d)
  //          {
  //            if (_zz[i][s][c][t][d].get(GRB_DoubleAttr_X) >= 0.4)
  //            {
  //              std::cout << _zz[i][s][c][t][d].get(GRB_StringAttr_VarName)
  //              << " = " << _zz[i][s][c][t][d].get(GRB_DoubleAttr_X) << std::endl;
  //            }
  //          }
  //        }
  //      }
  //    }


  if (_pNodeToStateSet != NULL){
    for (int i=0;i<anrSolutions; i++){
      delete _pNodeToStateSet[i];
    }
    delete _pNodeToStateSet;
  }
  if (_pNodeToRootState != NULL){
    for (int i=0;i<anrSolutions; i++){
      delete _pNodeToRootState[i];
    }
    delete _pNodeToRootState;
  }

  _pNodeToStateSet = new IntPairSetNodeMap*[anrSolutions];
  _pNodeToRootState = new IntPairNodeMap*[anrSolutions];
  _pTprime = new CloneTree*[anrSolutions];
  _pLPlus = new StringNodeMap*[anrSolutions];

  for(int nSoln=0; nSoln<anrSolutions; nSoln++){
    _model.set(GRB_IntParam_SolutionNumber, nSoln);
    constructGraph();
    _pNodeToStateSet[nSoln] = new IntPairSetNodeMap(getTree());
    _pNodeToRootState[nSoln] = new IntPairNodeMap(getTree());
//       for (int i = 0; i < nrNodes; ++i)
//  {
//    Node v_i = _indexToNode[i];
//    for (int s = 0; s < nrAnatomicalSites; ++s)
//    {
//      int size_L_s = _L[s].size();
//      for (int c = 0; c < size_L_s; ++c)
//      {
//        if (_x[i][s][c].get(GRB_DoubleAttr_Xn) >= 0.4)
//        {
//          std::cout << _x[i][s][c].get(GRB_StringAttr_VarName)
//          << " = " << _x[i][s][c].get(GRB_DoubleAttr_X) << std::endl;
//        }
//      }
//    }
//  }
// std::cout <<"..........."<<std::endl;
//  for (int s = 0; s < nrAnatomicalSites; ++s)
//  {
//    for (int t = 0; t < nrAnatomicalSites; ++t)
//    {
//      int size_L_t = _L[t].size();
//      for (int d = 0; d < size_L_t; ++d)
//      {
//        for (int e = 0; e < size_L_t; ++e)
//        {
//          if (_w[s][t][d][e].get(GRB_DoubleAttr_X) >= 0.4)
//          {
//            std::cout << _w[s][t][d][e].get(GRB_StringAttr_VarName)
//            << " = " << _w[s][t][d][e].get(GRB_DoubleAttr_X) << std::endl;
//          }
//        }
//      }
//    }
//  }
  
    for (int i = 0; i < nrNodes; ++i)
    {
      Node v_i = _indexToNode[i];
      for (int s = 0; s < nrAnatomicalSites; ++s)
      {
        int size_L_s = _L[s].size();
        for (int c = 0; c < size_L_s; ++c)
        {
          if (_r[i][s][c].get(GRB_DoubleAttr_Xn) >= 0.4)
          {
            //if(nSoln !=1) std::cout << nSoln<<" " << i << ' ' << s<<' '<<c<<' '
            //          << " = " << _r[i][s][c].get(GRB_DoubleAttr_Xn) << std::endl;
            (*_pNodeToRootState[nSoln])[v_i] = std::make_pair(s, c);
          }
        }
      }
    }
    
    for (int i = 0; i < nrNodes; ++i)
    {
      Node v_i = _indexToNode[i];
      for (int s = 0; s < nrAnatomicalSites; ++s)
      {
        int size_L_s = _L[s].size();
        for (int c = 0; c < size_L_s; ++c)
        {
          if (_x[i][s][c].get(GRB_DoubleAttr_Xn) >= 0.4)
          {
  //          std::cout << _x[i][s][c].get(GRB_StringAttr_VarName)
  //                    << " = " << _x[i][s][c].get(GRB_DoubleAttr_X) << std::endl;
            (*_pNodeToStateSet[nSoln])[v_i].insert(std::make_pair(s, c));
          }
        }
      }
    }
  //  const int nrArcs = _indexToArc.size();
  //  for (int ij = 0; ij < nrArcs; ++ij)
  //  {
  //    Arc a_ij = _indexToArc[ij];
  
  //    for (int s = 0; s < nrAnatomicalSites; ++s)
  //    {
  //      int size_L_s = _L[s].size();
  //      for (int c = 0; c < size_L_s; ++c)
  //      {
  //        for (int t = 0; t < nrAnatomicalSites; ++t)
  //        {
  //          int size_L_t = _L[t].size();
  //          for (int d = 0; d < size_L_t; ++d)
  //          {
  //            if (_xx[ij][s][c][t][d].get(GRB_DoubleAttr_X) >= 0.4)
  //            {
  //              std::cout << _xx[ij][s][c][t][d].get(GRB_StringAttr_VarName)
  //              << " = " << _xx[ij][s][c][t][d].get(GRB_DoubleAttr_X) << std::endl;
  //            }
  //          }
  //        }
  //      }
  //    }
  //  }
    
    BoolNodeMap leafPresence(getTree(), true);
    StringToStringMap toMutLabel;
    refine(leafPresence, toMutLabel, nSoln, post_processing);
  }
}

void IlpPmhTrSolver::initVariables()
{
  IlpPmhSolver::initVariables();
  
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  const int nrNodes = _indexToNode.size();
  
  char buf[1024];
  
  _zz = Var5Matrix(nrNodes);
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    const int i = (*_pNodeToIndex)[v_i];
    _zz[i] = Var4Matrix(nrAnatomicalSites);
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      _zz[i][s] = Var3Matrix(size_L_s);
      for (int c = 0; c < size_L_s; ++c)
      {
        _zz[i][s][c] = VarMatrix(nrAnatomicalSites);
        for (int t = 0; t < nrAnatomicalSites; ++t)
        {
          const int size_L_t = _L[t].size();
          _zz[i][s][c][t] = VarArray(size_L_t);
          for (int d = 0; d < size_L_t; ++d)
          {
            snprintf(buf, 1024, "zz;%s;%s;%d;%s;%d",
                     getLabel(v_i).c_str(),
                     _indexToAnatomicalSite[s].c_str(), c,
                     _indexToAnatomicalSite[t].c_str(), d);
            _zz[i][s][c][t][d] = strlen(buf) > 255 ? _model.addVar(0, 1, 0, GRB_BINARY) :  _model.addVar(0, 1, 0, GRB_BINARY, buf);
          }
        }
      }
    }
  }
  
  _r = Var3Matrix(nrNodes);
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (isLeaf(v_i)) continue;
    const int i = (*_pNodeToIndex)[v_i];
    _r[i] = VarMatrix(nrAnatomicalSites);
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      _r[i][s] = VarArray(size_L_s);
      for (int c = 0; c < size_L_s; ++c)
      {
        snprintf(buf, 1024, "r;%s;%s;%d",
                 getLabel(v_i).c_str(),
                 _indexToAnatomicalSite[s].c_str(), c);
        _r[i][s][c] = strlen(buf) > 255 ? _model.addVar(0, 1, 0, GRB_BINARY) : _model.addVar(0, 1, 0, GRB_BINARY, buf);
      }
    }
  }

  _p = Var3Matrix(nrNodes);
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    //if (isLeaf(v_i)) continue;
    const int i = (*_pNodeToIndex)[v_i];
    _p[i] = VarMatrix(nrAnatomicalSites);
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      _p[i][s] = VarArray(size_L_s);
      for (int c = 0; c < size_L_s; ++c)
      {
        snprintf(buf, 1024, "p;%s;%s;%d",
                 getLabel(v_i).c_str(),
                 _indexToAnatomicalSite[s].c_str(), c);
        _p[i][s][c] = strlen(buf) > 255 ? _model.addVar(0, 1, 0, GRB_BINARY) : _model.addVar(0, 1, 0, GRB_BINARY, buf);
      }
    }
  }
  
  _model.update();
}

void IlpPmhTrSolver::initConstraintsNonEdgesG()
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  
  GRBLinExpr sum;
  
  for (ArcIt a_ij(getTree()); a_ij != lemon::INVALID; ++a_ij)
  {
    const int ij = (*_pArcToIndex)[a_ij];
    Node v_i = getTree().source(a_ij);
    Node v_j = getTree().target(a_ij);
    const int i = (*_pNodeToIndex)[v_i];
    const int j = (*_pNodeToIndex)[v_j];
    
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        for (int t = 0; t < nrAnatomicalSites; ++t)
        {
          const int size_L_t = _L[t].size();
          for (int d = 0; d < size_L_t; ++d)
          {
            _model.addConstr(_xx[ij][s][c][t][d] <= _x[i][s][c]);
            _model.addConstr(_xx[ij][s][c][t][d] <= _r[j][t][d]);
            _model.addConstr(_xx[ij][s][c][t][d] >= _x[i][s][c] + _r[j][t][d] - 1);
          }
        }
      }
    }
  }
  
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const int size_L_s = _L[s].size();
    for (int c = 0; c < size_L_s; ++c)
    {
      for (int t = 0; t < nrAnatomicalSites; ++t)
      {
        const int size_L_t = _L[t].size();
        for (int d = 0; d < size_L_t; ++d)
        {
          if (s == t && c == d) continue;
          
          for (ArcIt a_ij(getTree()); a_ij != lemon::INVALID; ++a_ij)
          {
            const int ij = (*_pArcToIndex)[a_ij];
            sum += _xx[ij][s][c][t][d];
          }
          
          for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
          {
            const int i = (*_pNodeToIndex)[v_i];
            sum += _zz[i][s][c][t][d];
          }
          
          _model.addConstr(_z[s][c][t][d] <= sum);
          sum.clear();
        }
      }
    }
  }

  // cycle inequality constraints
  
  // only consider subset of size at most maximum out-degree
//  BoolVector subset(nrAnatomicalSites, false);
//  while (nextCombinationAnatomicalSites(subset))
//  {
//    IntPairVector selectedStates;
//    for (int i = 0; i < nrAnatomicalSites; ++i)
//    {
//      if (subset[i])
//        selectedStates.push_back(std::make_pair(i, 0));
//    }
//    
//    if (selectedStates.size() < 2) continue;
//    
//    do
//    {
//      do
//      {
//        for (int i = 1; i < selectedStates.size(); ++i)
//        {
//          int s = selectedStates[i-1].first;
//          int c = selectedStates[i-1].second;
//          int t = selectedStates[i].first;
//          int d = selectedStates[i].second;
//          assert(0 <= s && s < nrAnatomicalSites);
//          assert(0 <= t && t < nrAnatomicalSites);
//          assert(0 <= c && c < _L[s].size());
//          assert(0 <= d && d < _L[t].size());
//          sum += _z[s][c][t][d];
//        }
//        sum += _z[selectedStates.back().first][selectedStates.back().second][selectedStates.front().first][selectedStates.front().second];
//        _model.addConstr(sum <= selectedStates.size() - 1);
////        std::cout << "SUCCESS" << std::endl;
//        sum.clear();
//      } while (std::next_permutation(selectedStates.begin(),
//                                     selectedStates.end()));
//    }
//    while (nextCombinationStates(selectedStates));
//  }

  _model.update();
}

bool IlpPmhTrSolver::nextCombinationStates(IntPairVector& states) const
{
  int n = states.size();
  
  int i = 0;
  for (; i < n; ++i)
  {
    int s_i = states[i].first;
    if (states[i].second < _L[s_i].size() - 1)
    {
      break;
    }
  }
  
  if (i == n)
    return false;
  
  states[i].second += 1;
  
  for (int j = 0; j < i; ++j)
  {
    states[j].second = 0;
  }
  
  return true;
}

void IlpPmhTrSolver::initLeafVariables()
{
  IlpPmhSolver::initLeafVariables();
  
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  
  char buf[1024];
  
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (!isLeaf(v_i)) continue;
    const int i = (*_pNodeToIndex)[v_i];
    _r[i] = VarMatrix(nrAnatomicalSites);
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      _r[i][s] = VarArray(size_L_s);
      for (int c = 0; c < size_L_s; ++c)
      {
        snprintf(buf, 1024, "r;%s;%s;%d",
                 getLabel(v_i).c_str(),
                 _indexToAnatomicalSite[s].c_str(), c);
        _r[i][s][c] = strlen(buf) > 255 ? _model.addVar(0, 1, 0, GRB_BINARY) : _model.addVar(0, 1, 0, GRB_BINARY, buf);
      }
    }
  }
  
  _model.update();
}

void IlpPmhTrSolver::initLeafConstraints()
{
  IlpPmhSolver::initLeafConstraints();
  
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  GRBLinExpr sum, sum2;
  
  // one root color
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (!isLeaf(v_i)) continue;
    const int i = (*_pNodeToIndex)[v_i];
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        _model.addConstr(_r[i][s][c] <= _x[i][s][c]);
        sum += _r[i][s][c];
      }
    }
    _model.addConstr(sum == 1);
    sum.clear();
  }
  
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (!isLeaf(v_i)) continue;
    const int i = (*_pNodeToIndex)[v_i];
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        sum2 += _x[i][s][c];
        for (int t = 0; t < nrAnatomicalSites; ++t)
        {
          const int size_L_t = _L[t].size();
          for (int d = 0; d < size_L_t; ++d)
          {
            sum += _zz[i][s][c][t][d];
            _model.addConstr(_zz[i][s][c][t][d] <= _x[i][s][c]);
            _model.addConstr(_zz[i][s][c][t][d] <= _x[i][t][d]);
          }
        }
      }
    }
    _model.addConstr(sum2 - 1 == sum);
    sum.clear();
    sum2.clear();
  }
  
  _model.update();
}

void IlpPmhTrSolver::initConstraintsG()
{
  IlpPmhSolver::initConstraintsG();
  
  GRBLinExpr sum, sum2;
  
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  
  // require _zz[i] to induce a subtree
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (isLeaf(v_i)) continue;
    const int i = (*_pNodeToIndex)[v_i];
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        sum2 += _x[i][s][c];
        for (int t = 0; t < nrAnatomicalSites; ++t)
        {
          const int size_L_t = _L[t].size();
          for (int d = 0; d < size_L_t; ++d)
          {
            sum += _zz[i][s][c][t][d];
            _model.addConstr(_zz[i][s][c][t][d] <= _x[i][s][c]);
            _model.addConstr(_zz[i][s][c][t][d] <= _x[i][t][d]);
          }
        }
      }
    }
    _model.addConstr(sum2 - 1 == sum);
    sum.clear();
    sum2.clear();
  }
  
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (isLeaf(v_i)) continue;
    const int i = (*_pNodeToIndex)[v_i];
    
    for (int t = 0; t < nrAnatomicalSites; ++t)
    {
      const int size_L_t = _L[t].size();
      for (int d = 0; d < size_L_t; ++d)
      {
        for (int s = 0; s < nrAnatomicalSites; ++s)
        {
          const int size_L_s = _L[s].size();
          for (int c = 0; c < size_L_s; ++c)
          {
            sum += _zz[i][s][c][t][d];
          }
        }
        _model.addConstr(sum == _x[i][t][d] - _r[i][t][d]);
        sum.clear();
      }
    }
  }
  
  // one arc in G
  for (ArcIt a_ij(getTree()); a_ij != lemon::INVALID; ++a_ij)
  {
    Node v_j = getTree().target(a_ij);
    if (isLeaf(v_j)) continue;
    const int ij = (*_pArcToIndex)[a_ij];
    for (int t = 0; t < nrAnatomicalSites; ++t)
    {
      const int size_L_t = _L[t].size();
      for (int d = 0; d < size_L_t; ++d)
      {
        for (int s = 0; s < nrAnatomicalSites; ++s)
        {
          const int size_L_s = _L[s].size();
          for (int c = 0; c < size_L_s; ++c)
          {
            sum += _xx[ij][s][c][t][d];
          }
        }
        
//        _model.addConstr(sum <= 1);
        sum.clear();
      }
    }
  }
  
  // relate _zz to _z
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    const int i = (*_pNodeToIndex)[v_i];
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        for (int t = 0; t < nrAnatomicalSites; ++t)
        {
          const int size_L_t = _L[t].size();
          for (int d = 0; d < size_L_t; ++d)
          {
            _model.addConstr(_z[s][c][t][d] >= _zz[i][s][c][t][d]);
          }
        }
      }
    }
  }
  
  // no incoming edges to (P,0), root of _G
  for (int t = 0; t < nrAnatomicalSites; ++t)
  {
    const int size_L_t = _L[t].size();
    for (int d = 0; d < size_L_t; ++d)
    {
      sum += _z[t][d][_primaryIndex][0];
    }
  }
  _model.addConstr(sum == 0);
  sum.clear();
  
  _model.update();
}

void IlpPmhTrSolver::initConstraints()
{
  IlpPmhSolver::initConstraints();
  
  GRBLinExpr sum;
  
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  
  for (ArcIt a_ij(getTree()); a_ij != lemon::INVALID; ++a_ij)
  {
    Node v_i = getTree().source(a_ij);
    Node v_j = getTree().target(a_ij);
    
    int i = (*_pNodeToIndex)[v_i];
    int j = (*_pNodeToIndex)[v_j];
    
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        for (int t = s + 1; t < nrAnatomicalSites; ++t)
        {
          const int size_L_t = _L[t].size();
          for (int d = 0; d < size_L_t; ++d)
          {
            _model.addConstr(_x[i][s][c] + _x[j][s][c] + _x[i][t][d] + _x[j][t][d] <= 3);
          }
        }
      }
    }
  }
  
  // one root color
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    if (isLeaf(v_i)) continue;
    const int i = (*_pNodeToIndex)[v_i];
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        _model.addConstr(_r[i][s][c] <= _x[i][s][c]);
        sum += _r[i][s][c];
      }
    }
    _model.addConstr(sum == 1);
    sum.clear();
  }
  _model.addConstr(_r[(*_pNodeToIndex)[getRoot()]][_primaryIndex][0] == 1);
  
  _model.update();
}

GRBLinExpr IlpPmhTrSolver::initObjective(const IntTriple& bounds)
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
//  const int nrNodes = _indexToNode.size();
  
  // migration number
  GRBLinExpr migrationNumber;
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    const int size_L_s = _L[s].size();
    for (int c = 0; c < size_L_s; ++c)
    {
      if (s == _primaryIndex && c == 0) continue;
      migrationNumber += _y[s][c];
    }
  }
  
  if (bounds.first != -1)
  {
    _model.addConstr(migrationNumber <= bounds.first);
  }

  // comigration number
  GRBLinExpr comigrationNumber;
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    for (int t = 0; t < nrAnatomicalSites; ++t)
    {
      comigrationNumber += _gamma[s][t];
    }
  }
  
  if (bounds.second.first != -1)
  {
    _model.addConstr(comigrationNumber <= bounds.second.first);
  }
  
  GRBLinExpr obj;
  obj += migrationNumber
      + (1. / (nrAnatomicalSites * nrAnatomicalSites)) * comigrationNumber;
  
//  double factor = (1. / (nrAnatomicalSites * nrAnatomicalSites))
//                * (1. / (nrNodes * nrNodes));
//
//  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
//  {
//    const int i = (*_pNodeToIndex)[v_i];
//    for (int s = 0; s < nrAnatomicalSites; ++s)
//    {
//      const int size_L_s = _L[s].size();
//      for (int c = 0; c < size_L_s; ++c)
//      {
//        for (int t = 0; t < nrAnatomicalSites; ++t)
//        {
//          const int size_L_t = _L[t].size();
//          for (int d = 0; d < size_L_t; ++d)
//          {
//            obj += factor * _x[i][s][c];
//          }
//        }
//      }
//    }
//  }
  
  GRBLinExpr seedingSiteNumber;
  for (int s = 0; s < nrAnatomicalSites; ++s)
  {
    seedingSiteNumber += _sigma[s];
  }
  
  if (bounds.second.second != -1)
  {
    _model.addConstr(seedingSiteNumber <= bounds.second.second);
  }
  
  double factor = (1. / (nrAnatomicalSites * nrAnatomicalSites))
//         * (1. / (nrNodes * nrNodes))
                * (1. / (nrAnatomicalSites + 1));
  
  obj += factor * seedingSiteNumber;
  
  obj = 1000 * obj;
  
  _model.setObjective(obj, GRB_MINIMIZE);
  _model.update();
  
  return obj;
}

void IlpPmhTrSolver::initVertexLabelingConstraints()
{
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  GRBLinExpr sum, sum2;
  
  // Every vertex is labeled by at least one (s,c)
  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i)
  {
    const int i = (*_pNodeToIndex)[v_i];
    if (isLeaf(v_i)) continue;
    for (int s = 0; s < nrAnatomicalSites; ++s)
    {
      const int size_L_s = _L[s].size();
      for (int c = 0; c < size_L_s; ++c)
      {
        sum += _x[i][s][c];
        sum2 += _x[i][s][c];
      }
      _model.addConstr(sum2 <= 1);
      sum2.clear();
    }
    
    if (lemon::countOutArcs(getTree(), v_i) > 2)
    {
      _model.addConstr(sum >= 1);
    }
    else
    {
      _model.addConstr(sum == 1);
    }
    sum.clear();
  }
}

void IlpPmhTrSolver::symmetryBreakingConstraints(){

  //p related constraints
  _model.addConstr(_p[(*_pNodeToIndex)[getRoot()]][_primaryIndex][0] == 1);
  const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  const int nrNodes = _indexToNode.size();
  GRBLinExpr sum;
  for (int t = 0; t < nrAnatomicalSites; t++){
    const int size_L_t = _L[t].size();
    for (int d=0; d<size_L_t; d++){
      for(int i=0; i<nrNodes; i++){
        sum += _p[i][t][d];
        Node v_i = _indexToNode[i];
        for (int s = 0; s < nrAnatomicalSites; s++){
          const int size_L_s = _L[s].size();
          for (int c=0; c<size_L_s; c++){
            _model.addConstr(_p[i][t][d] >= _zz[i][s][c][t][d]);
            const int p_i = (*_pNodeToIndex)[getParent(v_i)];
            Arc a_pii = InArcIt(getTree(), v_i);
            const int pii = (*_pArcToIndex)[a_pii];
            _model.addConstr(_p[i][t][d] >= _xx[pii][s][c][t][d] - _x[p_i][t][d]);
          }
        }
      }
      //_model.addConstr(sum == _y[t][d]);
      sum.clear();
    }
  }

  for (int t = 0; t < nrAnatomicalSites; t++){
    for(int i=0; i<nrNodes; i++){
      for(int j=i+1; j<nrNodes; j++){
        const int size_L_t = _L[t].size();
        for (int d1=0; d1<size_L_t; d1++){
          for (int d2=d1+1; d2<size_L_t;d2++){
            _model.addConstr(_p[i][t][d2]+_p[j][t][d1] <= 1);
          }
        }
      }
    }
  }




  // const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  // const int nrNodes = _indexToNode.size();
  // GRBLinExpr sum;
  // for (int s = 0; s < nrAnatomicalSites; s++){
  //   const int size_L_s = _L[s].size();
  //   for (int c1=0; c1<size_L_s; c1++){
  //     for (int c2=c1+1; c2<size_L_s;c2++){
  //       for(int i=0; i<nrNodes; i++){
  //         sum += std::pow(2, nrNodes - i - 1) * (_x[i][s][c1] - _x[i][s][c2]);
  //       }
  //       _model.addConstr(sum >= 0);
  //       sum.clear();
  //     }
  //   }
  // }

  // const int nrAnatomicalSites = _anatomicalSiteToIndex.size();
  // const int nrNodes = _indexToNode.size();
  // GRBLinExpr sum;
  // for (int t = 0; t < nrAnatomicalSites; t++){
  //   for(int i=0; i<nrNodes; i++){
  //     for(int j=i+1; j<nrNodes; j++){
  //       const int size_L_t = _L[t].size();
  //       for (int d1=0; d1<size_L_t; d1++){
  //         for (int d2=d1+1; d2<size_L_t;d2++){
  //           for (int s1 = 0; s1 < nrAnatomicalSites; s1++){
  //             const int size_L_s1 = _L[s1].size();
  //             for (int c1=0; c1<size_L_s1; c1++){
  //               for (int s2 = 0; s2 < nrAnatomicalSites; s2++){
  //                 const int size_L_s2 = _L[s2].size();
  //                 for (int c2=0; c2<size_L_s2; c2++){
  //                   Node v_i = _indexToNode[i];
  //                   Node v_j = _indexToNode[j];
  //                   const int p_i = (*_pNodeToIndex)[getParent(v_i)];
  //                   const int p_j = (*_pNodeToIndex)[getParent(v_j)];
  //                   if (v_j == getRoot()) continue;
  //                   Arc a_pii = InArcIt(getTree(), v_i);
  //                   Arc a_pjj = InArcIt(getTree(), v_j);
  //                   const int pii = (*_pArcToIndex)[a_pii];
  //                   const int pjj = (*_pArcToIndex)[a_pjj];
  //                    _model.addConstr(_xx[pjj][s1][c1][t][d1] - _x[p_j][t][d1] + _zz[i][s2][c2][t][d2] <= 1);
  //                    _model.addConstr(_zz[j][s1][c1][t][d1] + _zz[i][s2][c2][t][d2] <= 1);
  //                    _model.addConstr(_xx[pjj][s1][c1][t][d1] - _x[p_j][t][d1] + _xx[pii][s2][c2][t][d2] - _x[p_i][t][d2] <= 1);
  //                    _model.addConstr(_zz[j][s1][c1][t][d1] + _xx[pii][s2][c2][t][d2] - _x[p_i][t][d2] <= 1);
  //                 }
  //               }
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  // }


  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i){
    if (lemon::countOutArcs(getTree(), v_i) > 2){
      const int i = (*_pNodeToIndex)[v_i];
      for (int s = 0; s < nrAnatomicalSites; ++s){
        const int size_L_s = _L[s].size();
        for (int c = 0; c < size_L_s; ++c){
          for(OutArcIt m(getTree(), v_i); m != lemon::INVALID; ++m){
            const int j = (*_pNodeToIndex)[getTree().target(m)];
            sum += _r[j][s][c];
          }
          _model.addConstr(sum >= 2 * (_x[i][s][c] - _r[i][s][c]));
          sum.clear();
        }
      }
    }
  }

  for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i){
    if (lemon::countOutArcs(getTree(), v_i) > 2){
      const int i = (*_pNodeToIndex)[v_i];
      for (int s = 0; s < nrAnatomicalSites; ++s){
        const int size_L_s = _L[s].size();
        for (int c = 0; c < size_L_s; ++c){
          for (int t = 0; t < nrAnatomicalSites; ++t){
            const int size_L_t = _L[t].size();
            for (int d = 0; d < size_L_t; ++d){
              sum += _zz[i][s][c][t][d];
            }
          }
          for(OutArcIt m(getTree(), v_i); m != lemon::INVALID; ++m){
            const int j = (*_pNodeToIndex)[getTree().target(m)];
            sum += _r[j][s][c];
          }
          for (int t = 0; t < nrAnatomicalSites; ++t){
            const int size_L_t = _L[t].size();
            for (int d = 0; d < size_L_t; ++d){
              for(OutArcIt m(getTree(), v_i); m != lemon::INVALID; ++m){
                const int j = (*_pNodeToIndex)[getTree().target(m)];
                const int ij = (*_pArcToIndex)[m];
                sum += (_r[j][t][d] - _xx[ij][t][d][t][d]);
              }
            }
          }
          _model.addConstr(sum >= 2 * _r[i][s][c]);
          sum.clear();
        }
      }
    }
  }

  // for (NodeIt v_i(getTree()); v_i != lemon::INVALID; ++v_i){
  //   if (lemon::countOutArcs(getTree(), v_i) > 0){
  //     const int i = (*_pNodeToIndex)[v_i];
  //     for (int s = 0; s < nrAnatomicalSites; ++s){
  //       const int size_L_s = _L[s].size();
  //       for (int c = 0; c < size_L_s; ++c){
  //         for (int t = 0; t < nrAnatomicalSites; ++t){
  //           const int size_L_t = _L[t].size();
  //           for (int d = 0; d < size_L_t; ++d){
  //             sum += _zz[i][s][c][t][d];
  //             for(OutArcIt m(getTree(), v_i); m != lemon::INVALID; ++m){
  //               const int ij = (*_pArcToIndex)[m];
  //               sum += _xx[ij][s][c][t][d];
  //             }
  //           }
  //         }
  //       }
  //     }
  //     _model.addConstr(sum >= 2 * _r[i][s][c]);
  //     sum.clear();
  //   }
  //}

  
  for (int s = 0; s < nrAnatomicalSites; ++s){
    const int size_L_s = _L[s].size();
    for (int c = 1; c < size_L_s; ++c){
      _model.addConstr(_y[s][c] <= _y[s][c - 1]);
    }
  }
}

IntTriple IlpPmhTrSolver::run(const CloneTree& T,
                              const std::string& primary,
                              const std::string& outputDirectory,
                              const std::string& outputPrefix,
                              const StringToIntMap& colorMap,
                              MigrationGraph::Pattern pattern,
                              int nrThreads,
                              bool outputILP,
                              bool outputSearchGraph,
                              int timeLimit,
                              const IntTriple& bounds,
                              const StringPairList& forcedComigrations,
                              int nrSolutions,
                              bool post_processing,
                              bool count_mode)
{
  std::string filenameGurobiLog;
  if (!outputDirectory.empty())
  {
    char buf[1024];
    snprintf(buf, 1024, "%s/%slog-%s-%s.txt",
             outputDirectory.c_str(),
             outputPrefix.c_str(),
             primary.c_str(),
             MigrationGraph::getPatternString(pattern).c_str());
    
    filenameGurobiLog = buf;
  }
  
  IlpPmhTrSolver solver(T,
                        primary,
                        pattern,
                        filenameGurobiLog,
                        forcedComigrations);
  
  return IlpPmhSolver::run(solver,
                           T,
                           primary,
                           outputDirectory,
                           outputPrefix,
                           colorMap,
                           pattern,
                           nrThreads,
                           outputILP,
                           outputSearchGraph,
                           timeLimit,
                           bounds,
                           forcedComigrations,
                           nrSolutions,
                           post_processing,
                           count_mode);
}

void IlpPmhTrSolver::refine(const BoolNodeMap& leafPresence,
                            StringToStringMap& toMutLabel,
                            int solIdx,
                            bool post_processing)
{
  // 1. refine
  Digraph Tprime;
  Node root_Tprime = Tprime.addNode();
  StringNodeMap label(Tprime);
  StringNodeMap lPlus(Tprime);
  label[root_Tprime] = getLabel(getRoot());
  toMutLabel[label[root_Tprime]] = getLabel(getRoot());
  lPlus[root_Tprime] = _indexToAnatomicalSite[_primaryIndex];
  
  refine(leafPresence, toMutLabel, getRoot(), Tprime, root_Tprime, label, lPlus, solIdx);
  
  if(post_processing){
    std::queue<Node> q;
    q.push(root_Tprime);
    while(!q.empty()){
      Node curr = q.front();
      q.pop();
      for(OutArcIt p(Tprime, curr); p != lemon::INVALID; ++p){
        Node child = Tprime.target(p);
        if (label[child].find('^') != std::string::npos && lemon::countOutArcs(Tprime, child) == 1){
          for(OutArcIt m(Tprime, child); m != lemon::INVALID; ++m){
            Node gchild = Tprime.target(m);
            Tprime.addArc(curr, gchild);
            q.push(gchild);
          }
          Tprime.erase(child);
        } else {
          q.push(child);
        }
      }
    }
  }
  
  
  // 2. construct clone tree
  _pTprime[solIdx] = new CloneTree(Tprime, root_Tprime, label, lPlus);
  _pLPlus[solIdx] = new StringNodeMap(_pTprime[solIdx]->tree());
  
  for (NodeIt v(Tprime); v != lemon::INVALID; ++v)
  {
    const std::string& label_v = label[v];
    Node vv = _pTprime[solIdx]->getNodeByLabel(label_v);
    _pLPlus[solIdx]->set(vv, lPlus[v]);
  }
}

bool IlpPmhTrSolver::isConnected(const IntPairSet& Sigma_u,
                                 IntPair& root_sc) const
{
  // Check whether Sigma_u induces a connected subgraph of S
  BoolNodeMap filter(_G, false);
  for (IntPair sc : Sigma_u)
  {
    assert(_subLabelToNodeG[sc.first][sc.second] != lemon::INVALID);
    filter[_subLabelToNodeG[sc.first][sc.second]] = true;
  }
  
  root_sc = std::make_pair(-1, -1);
  for (IntPair sc : Sigma_u)
  {
    assert(_subLabelToNodeG[sc.first][sc.second] != lemon::INVALID);
    Node v_sc = _subLabelToNodeG[sc.first][sc.second];
    Arc a = InArcIt(_G, v_sc);
    if (a == lemon::INVALID)
    {
      root_sc = sc;
    }
    else if (!filter[_G.source(a)])
    {
      if (root_sc == std::make_pair(-1, -1))
      {
        root_sc = sc;
      }
      else
      {
        root_sc = std::make_pair(-1, -1);
        return false;
      }
    }
  }
  assert (root_sc != std::make_pair(-1, -1));
  return true;
}

bool IlpPmhTrSolver::areSiblings(const IntPairSet& Sigma_u,
                                 IntPair& parent_sc) const
{
  parent_sc = std::make_pair(-1, -1);
  for (IntPair sc : Sigma_u)
  {
    assert(_subLabelToNodeG[sc.first][sc.second] != lemon::INVALID);
    Node v_sc = _subLabelToNodeG[sc.first][sc.second];
    
    if (v_sc == _rootG)
    {
      return false;
    }
    
    Node v_dt = _G.source(InArcIt(_G, v_sc));
    if (parent_sc == std::make_pair(-1, -1))
    {
      parent_sc = _nodeGToSubLabel[v_dt];
    }
    else if (parent_sc != _nodeGToSubLabel[v_dt])
    {
      parent_sc = std::make_pair(-1, -1);
      return false;
    }
  }
  
  return true;
}

bool IlpPmhTrSolver::nextCombinationAnatomicalSites(BoolVector& subset) const
{
  int n = subset.size();
  for (int i = 0; i < n; ++i)
  {
    if (subset[i] == false)
    {
      subset[i] = true;
      for (int j = 0; j < i; ++j)
      {
        subset[j] = false;
      }
      return true;
    }
  }
  
  return false;
}

void IlpPmhTrSolver::refine(const BoolNodeMap& leafPresence,
                            StringToStringMap& toMutLabel,
                            Node u, //Node v_inT,
                            Digraph& Tprime,
                            Node uu, //Node v_inTprime,
                            StringNodeMap& label,
                            StringNodeMap& lPlus,
                            int solIdx)
{
  const IntPairSet& Sigma_u = (*_pNodeToStateSet[solIdx])[u];
  if (Sigma_u.size() <= 1)
  {
    // leave vertex u unperturbed
    if (isLeaf(u))
    {
      lPlus[uu] = getLeafAnatomicalSiteLabel(u);
    }
    else
    {
      if (u != getRoot())
      {
        if (Sigma_u.empty())
        {
          Node pi_uu = Tprime.source(InArcIt(Tprime, uu));
          lPlus[uu] = lPlus[pi_uu];std::cout <<"does it happen?";
        }
        else
        {
          IntPair sc = *(Sigma_u.begin());
          lPlus[uu] = _indexToAnatomicalSite[sc.first];
        }
      }
      for (OutArcIt a(getTree(), u); a != lemon::INVALID; ++a)
      {
        Node v = getTree().target(a);
        if (!leafPresence[v])
          continue;
        Node vv = Tprime.addNode();
        Tprime.addArc(uu, vv);
        
        label[vv] = getLabel(v);
        assert(toMutLabel.count(label[vv]) == 0);
        toMutLabel[label[vv]] = getLabel(v);
        refine(leafPresence, toMutLabel, v, Tprime, vv, label, lPlus, solIdx);
      }
    }
  }
  else
  {
    // Check whether Sigma_u induces a connected subgraph of S
    IntPair pi_sc;
    IntPairToNodeMap toBackBone;
    if (isConnected(Sigma_u, pi_sc))
    {
      // add vertices
      for (IntPair sc : Sigma_u)
      {
        if (sc != pi_sc)
        {
          Node v_sc = Tprime.addNode();
          toBackBone[sc] = v_sc;
          label[v_sc] = getLabel(u) + "^" + _indexToAnatomicalSite[sc.first];
          lPlus[v_sc] = _indexToAnatomicalSite[sc.first];
          
          assert(toMutLabel.count(label[v_sc]) == 0);
          toMutLabel[label[v_sc]] = getLabel(u);
        }
        else
        {
          toBackBone[sc] = uu;
          lPlus[uu] = _indexToAnatomicalSite[sc.first];
        }
      }
      
      // add edges
      for (IntPair sc : Sigma_u)
      {
        if (sc != pi_sc)
        {
          Node v_sc = toBackBone[sc];
          
          IntPair td = _nodeGToSubLabel[_G.source(InArcIt(_G, _subLabelToNodeG[sc.first][sc.second]))];
          Node v_td = toBackBone[td];
          Tprime.addArc(v_td, v_sc);
        }
      }
    }
    else if (areSiblings(Sigma_u, pi_sc))
    {
      lPlus[uu] = _indexToAnatomicalSite[pi_sc.first];
      
      // add vertices
      for (IntPair sc : Sigma_u)
      {
        Node v_sc = Tprime.addNode();
        toBackBone[sc] = v_sc;
        label[v_sc] = getLabel(u) + "^" + _indexToAnatomicalSite[sc.first];
        lPlus[v_sc] = _indexToAnatomicalSite[sc.first];
        Tprime.addArc(uu, v_sc);
      }
    }
    else
    {
      /*writeDOT(std::cout);
      writeDOT(std::cout, S);
      for (IntPair sc : Sigma_u)
      {
        std::cout << _indexToAnatomicalSite[sc.first] << " , " << sc.second << std::endl;
      }*/
      assert(false);
    }
    
    // add original children
    for (OutArcIt a(getTree(), u); a != lemon::INVALID; ++a)
    {
      Node v = getTree().target(a);
      const IntPair& td = (*_pNodeToRootState[solIdx])[v];
      
      if (!leafPresence[v])
        continue;
      
      Node vv = Tprime.addNode();
      label[vv] = getLabel(v);

      if (toBackBone.count(td) == 1)
      {
        Tprime.addArc(toBackBone[td], vv);
      }
      else
      {
        Node v_td = _subLabelToNodeG[td.first][td.second];
        assert(v_td != _rootG);
        Node parent_v_td = _G.source(InArcIt(_G, v_td));
        const IntPair& parent_td = _nodeGToSubLabel[parent_v_td];

        assert(toBackBone.count(parent_td) == 1);
        Tprime.addArc(toBackBone[parent_td], vv);
      }
      
      assert(toMutLabel.count(label[vv]) == 0);
      toMutLabel[label[vv]] = getLabel(v);
      refine(leafPresence, toMutLabel, v, Tprime, vv, label, lPlus, solIdx);
    }
  }
}

void IlpPmhTrSolver::initCallbacks()
{
  _model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
  _pCallback = new IlpPmhPrSolverCycleElimination(_indexToAnatomicalSite, _primaryIndex, _y, _z);
  _model.setCallback(_pCallback);
}
