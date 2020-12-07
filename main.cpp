#include <cmath>
#include <iostream>
#include <random>
#include <vector>

void vec_noisy_fun_value(std::vector<double> &noisy_fun_value) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(-0.25, 0.25);
  for (auto k = 0; k < 100; ++k) {
    double x = k * M_PI / 100;
    noisy_fun_value.push_back(sin(x) + 0.5 + dist(gen));
  }
}

void vector_vec_alpha(std::vector<std::vector<double>> &v_v_alpha) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(0, 1);
  for (auto k = 0; k < 100; ++k) {
    double arg = dist(gen);
    v_v_alpha.push_back(std::vector<double>{arg, (1 - arg) / 2, (1 - arg) / 2});
  }
}

void vec_filtered_fun_value(const std::vector<double> &noisy_fun_value,
                            const std::vector<std::vector<double>> &v_v_alpha,
                            std::vector<double> &filtered_fun_value,
                            const int &r) {
  int m = (r - 1) / 2;
  for (int k = m; k < 100 - m; ++k) {
    double sum = 0;
    for (int i = k - m; i < k + m; ++i) {
      sum += v_v_alpha[k][i + m + 1 - k] / noisy_fun_value[i];
    }
    filtered_fun_value.push_back(pow(sum, -1));
  }
}

int main() {
  const int r = 3;
  std::vector<double> noisy_fun_value = {};
  std::vector<std::vector<double>> v_v_alpha = {};
  std::vector<double> filtered_fun_value = {};
  vec_noisy_fun_value(noisy_fun_value);
  vector_vec_alpha(v_v_alpha);
  vec_filtered_fun_value(noisy_fun_value, v_v_alpha, filtered_fun_value, r);
  for (auto &it : filtered_fun_value) {
    std::cout << it << std::endl;
  }
  return 0;
}
