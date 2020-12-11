#include <cmath>
#include <iomanip>
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

void vector_vec_alpha(std::vector<double> &alpha, const int &r,
                      std::vector<double> rand_values, const int &iter) {
  alpha.push_back(rand_values[iter]);
  for (auto i = 0; i < r - 1; ++i) {
    alpha.push_back((1 - rand_values[iter]) / (r - 1));
  }
}

void vec_filtered_fun_value(const std::vector<double> &noisy_fun_value,
                            const std::vector<double> &alpha,
                            std::vector<double> &filtered_fun_value,
                            const int &r) {
  int m = (r - 1) / 2;
  for (int k = m; k < 100 - m; ++k) {
    double sum = 0;
    for (int i = k - m; i <= k + m; ++i) {
      sum += alpha[i + m + 1 - k - 1] / noisy_fun_value[i];
    }
    filtered_fun_value.push_back(pow(sum, -1));
  }
}

double search_noise_criteria(const std::vector<double> &filtered_fun_value) {
  double sum = 0;
  for (auto i = 1; i < 98; ++i) {
    sum += pow(filtered_fun_value[i] - filtered_fun_value[i - 1], 2);
  }
  return pow(sum, 0.5);
}

double
search_difference_criterion(const std::vector<double> &filtered_fun_value,
                            const std::vector<double> &noisy_fun_value) {
  double sum = 0;
  for (auto i = 0; i < 98; ++i) {
    sum += pow(filtered_fun_value[i] - noisy_fun_value[i], 2);
  }
  return pow(sum * 0.01, 0.5);
}

double search_functionality(const double &noise_criteria,
                            const double &difference_criterion) {
  double n = 10;
  double min = 1000000;
  double lymbda = 0;
  for (auto i = 0; i <= n; ++i) {
    double arg =
        (i / n) * noise_criteria + (1 - (i / n)) * difference_criterion;
    if (arg < min) {
      min = arg;
      lymbda = i / n;
    }
  }
  return min;
}

double search_dist(const double &noise_criteria,
                   const double &difference_criterion) {
  return pow(noise_criteria * noise_criteria +
                 difference_criterion * difference_criterion,
             0.5);
}

void vec_rand_value(std::vector<double> &rand_values) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(0, 1);
  for (auto i = 0; i < 10; ++i) {
    double arg = dist(gen);
    rand_values.push_back(arg);
  }
}

void print_head(const int &r) {
  std::cout << "+-----+---------+";
  for (auto i = 0; i < r * 8; ++i) {
    std::cout << "-";
  }
  std::cout << "+--------+--------+" << std::endl;
  std::cout << "|  h  |"
            << "   dis   |";
  for (auto i = 0; i < r * 4 - 2; ++i) {
    std::cout << " ";
  }
  std::cout << "alph";
  for (auto i = 0; i < r * 4 - 2; ++i) {
    std::cout << " ";
  }
  std::cout << "|   w    |"
            << "    d   |" << std::endl;
  std::cout << "+-----+---------+";
  for (auto i = 0; i < r * 8; ++i) {
    std::cout << "-";
  }
  std::cout << "+--------+--------+" << std::endl;
}

void print_data(const std::vector<double> &alpha, const int &iter,
                const double &dist, const double &noise_criteria,
                const double &difference_criterion) {
  std::cout << std::right << "|" << std::setw(5) << iter * 0.1
            << std::setprecision(4) << "|" << std::setw(9) << dist << "|";
  for (auto &it : alpha) {
    std::cout << std::setprecision(3) << std::setw(8) << it;
  }
  std::cout << "|" << std::setprecision(4) << std::right << std::setw(8)
            << noise_criteria << "|" << std::setw(8) << difference_criterion
            << "|" << std::endl;
}

void print_end(const int &min_iter, const double &min_functionality,
               const double &min_noise_criteria,
               const double &min_difference_criterion) {
  std::cout << std::endl << "+-----+--------+--------+--------+" << std::endl;
  std::cout << "|  h* |   J    |   w    |    d   |" << std::endl;
  std::cout << "+-----+--------+--------+--------+" << std::endl;
  std::cout << "|" << std::setw(5) << min_iter*0.1 << "|" << std::setw(8)
            << min_functionality << "|" << std::setw(8)
            << min_noise_criteria << "|" << std::setw(8)
            << min_difference_criterion <<"|"<< std::endl;
  std::cout << "+-----+--------+--------+--------+" << std::endl;
}

void experiment(const int &r) {
  double noise_criteria = 0;
  double difference_criterion = 0;
  double functionality = 0;
  double dist = 0;
  double min_dist = 10000;
  double min_functionality = 0;
  double min_noise_criteria = 0;
  double min_difference_criterion = 0;
  int min_iteration = 0;
  std::vector<double> noisy_fun_value = {};
  vec_noisy_fun_value(noisy_fun_value);
  print_head(r);
  for (auto i = 0; i < 11; ++i) {
    std::vector<double> filtered_fun_value = {};
    std::vector<double> rand_values = {};
    vec_rand_value(rand_values);
    std::vector<double> alpha = {};
    vector_vec_alpha(alpha, r, rand_values, i);
    std::cout << std::endl;
    vec_filtered_fun_value(noisy_fun_value, alpha, filtered_fun_value, r);
    noise_criteria = search_noise_criteria(filtered_fun_value);
    difference_criterion =
        search_difference_criterion(filtered_fun_value, noisy_fun_value);
    dist = search_dist(noise_criteria, difference_criterion);
    functionality = search_functionality(noise_criteria, difference_criterion);
    print_data(alpha, i, dist, noise_criteria, difference_criterion);
    if (min_dist > dist) {
      min_dist = dist;
      min_functionality = functionality;
      min_noise_criteria = noise_criteria;
      min_difference_criterion = difference_criterion;
      min_iteration = i;
    }
  }
  print_end(min_iteration, min_functionality, min_noise_criteria,
            min_difference_criterion);
}

int main() {
  int r = 3;
  experiment(r);
  r = 5;
  experiment(r);

  return 0;
}
