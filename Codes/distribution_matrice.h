#include <vector>

std::vector<std::vector<int>> distrubtion_matrix(std::vector<std::vector<int>> v){
  int n = v.size(), m = v.front().size();
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

// perm[i] = row my column.
std::vector<std::vector<int>> semi_local_lcs(std::vector<int> perm, int m){
	int n = perm.size();
	std::vector<std::vector<int>> perm_matrix(n,std::vector<int> (n,0));
	for(int i = 0; i < n; ++i){
		perm_matrix[perm[i]][i] = 1;
	}
	auto dist_matrix = distrubtion_matrix(perm_matrix);
	std::vector<std::vector<int>> H(n + 1, std::vector<int> (n + 1));
	for(int i = 0; i <= n; ++i){
		for(int j = 0; j <= n; ++j){
			H[i][j] = j - i - dist_matrix[i][j] + m;
			std::cout << H[i][j] << " ";
		}
		std::cout << "\n";
	}
	return H;
}