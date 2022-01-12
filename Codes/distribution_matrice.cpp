#include <vector>

std::vector<std::vector<int>> distrubtion_matrix(std::vector<std::vector<int>> v){
  int n = v.size(), m = v.front().size();
  assert(n == m);
  std::vector<std::vector<int>> ret(n + 1, std::vector<int> (m + 1, 0));
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < m; ++j){
      for(int x = i; x < n; ++x){
        for(int y = 0; y <= j; ++y){
          ret[i][j+1] += v[x][y];
        }
      }
    }
  }
  return ret;
}

