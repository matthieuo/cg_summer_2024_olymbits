#pragma GCC optimize("-Ofast")
#pragma GCC optimize("inline")
#pragma GCC optimize("omit-frame-pointer")
#pragma GCC optimize("unroll-loops")
#pragma GCC optimize("tree-vectorize")
#pragma GCC optimize("fast-math")
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using namespace std;
typedef tuple<std::string, int, int, int, int, int, int, int> input_tuple;

enum class Action { UP, LEFT, DOWN, RIGHT };

std::string to_string(Action action) {
  switch (action) {
  case Action::UP:
    return "UP";
  case Action::LEFT:
    return "LEFT";
  case Action::DOWN:
    return "DOWN";
  case Action::RIGHT:
    return "RIGHT";
  default:
    return "Unknown";
  }
}

std::ostream &operator<<(std::ostream &os, Action action) {
  return os << to_string(action);
}

constexpr std::size_t NP = 3;

// ***************************************************
std::vector<std::array<Action, 3>> generate_combinations() {
  std::vector<std::array<Action, 3>> combinations;
  for (auto a1 : {Action::UP, Action::DOWN, Action::LEFT, Action::RIGHT}) {
    for (auto a2 : {Action::UP, Action::DOWN, Action::LEFT, Action::RIGHT}) {
      for (auto a3 : {Action::UP, Action::DOWN, Action::LEFT, Action::RIGHT}) {
        combinations.push_back({a1, a2, a3});
      }
    }
  }
  return combinations;
}

std::vector<int> string_to_vector(const std::string &str) {
  std::vector<int> result;

  for (char ch : str) {
    if (std::isdigit(ch)) {
      result.push_back(ch - '0');
    }
  }

  return result;
}

template <std::size_t N, typename T>
std::array<int, N> rank_numbers(const std::array<T, N> &numbers) {
  std::array<std::pair<T, int>, N> number_indices;
  for (size_t i = 0; i < N; ++i) {
    number_indices[i] = std::make_pair(numbers[i], i);
  }

  std::sort(number_indices.begin(), number_indices.end(), std::greater<>());
  std::array<int, N> ranks;

  int current_rank = 1;
  ranks[number_indices[0].second] = current_rank;
  for (size_t i = 1; i < N; ++i) {
    if (number_indices[i].first == number_indices[i - 1].first) {
      ranks[number_indices[i].second] = ranks[number_indices[i - 1].second];
    } else {
      current_rank = i + 1;
      ranks[number_indices[i].second] = current_rank;
    }
  }

  return ranks;
}

//-------------- COURSE ----------------

template <std::size_t N> class Course {
public:
  Course(const std::string &gpu_, const std::array<int, N> &pos,
         const std::array<int, N> &diz)
      : pos(pos), diz(diz), medals({4, 4, 4}) {

    is_finish = (gpu_ == "GAME_OVER");
    if (is_finish) {
      auto sorted_array = rank_numbers(pos);
      for (int i = 0; i < sorted_array.size(); ++i) {
        medals[i] = sorted_array[i];
      }
    } else {
      cerr << gpu_ << gpu_.length() << endl;
      assert(gpu_.length() == 30);

      for (int i = 0; i < gpu_.length(); ++i) {
        gpu[i] = gpu_[i];
      }
    }
  }

  Course()
      : gpu(""), pos({0, 0, 0}), diz({0, 0, 0}), medals({4, 4, 4}),
        is_finish(false) {}

  bool game_finish() const { return is_finish; }

  std::array<int, N> get_medals() const { return medals; }

  Course action(const std::array<Action, N> &actions) const {

    if (is_finish) {
      return *this;
    }

    std::array<int, N> new_pos = pos;
    std::array<int, N> new_diz = diz;

    bool is_game_finish = false;

    for (int i = 0; i < 3; ++i) {
      auto [x, y] = update_agent_pos_diz(pos[i], diz[i], actions[i]);
      new_pos[i] = x;
      new_diz[i] = y;
      if (x == 30 - 1) { // gpu.size() - 1) {
        is_game_finish = true;
      }
    }

    auto c = *this;

    c.pos = new_pos;
    c.diz = new_diz;

    if (is_game_finish) {
      auto sorted_array = rank_numbers(new_pos);
      for (int i = 0; i < sorted_array.size(); ++i) {
        c.medals[i] = sorted_array[i];
        c.is_finish = true;
      }
    }

    return c;
  }

  void print_game_state() const {
    if (is_finish) {
      std::cerr << "FINISH  " << std::endl;
    }
    std::cerr << "Piste : " << gpu << std::endl;
    for (int i = 0; i < N; ++i) {
      std::cerr << "Joueur " << i + 1 << " : Position " << pos[i]
                << ", Étourdissement " << diz[i] << "M" << medals[i]
                << std::endl;
    }
  }

  // private:
  std::tuple<int, int> update_agent_pos_diz(int init_pos, int init_diz,
                                            Action action) const {

    if (init_diz > 0) { // diz -> not move
      return {init_pos, init_diz - 1};
    }
    int target_pos = init_pos;
    assert(init_diz == 0);
    int target_diz = 0;

    switch (action) {
    case Action::UP:
      target_pos += 2;
      break;
    case Action::LEFT:
      target_pos += 1;
      break;
    case Action::DOWN:
      target_pos += 2;
      break;
    case Action::RIGHT:
      target_pos += 3;
      break;
    }

    target_pos = std::min(target_pos, 30 - 1);

    for (int i = init_pos + 1; i <= target_pos; ++i) {
      if ((target_pos - i == 1) && action == Action::UP) {
        // don't check colision at the first if we UP
        continue;
      }
      if (gpu[i] == '#') {
        // colision
        target_diz = 2;
        target_pos = i;
        break;
      }
    }
    return {target_pos, target_diz};
  }

  char gpu[30];
  std::array<int, N> pos;
  std::array<int, N> diz;

  std::array<int, N> medals; // Medals
  bool is_finish;
};

//------------------------------------ ARCHERY ------

template <std::size_t N> class Archery {
public:
  Archery(const std::string gpu,
          const std::array<std::array<int, 2>, N> cursors)
      : cursors(cursors), medals({4, 4, 4}) {

    if (gpu == "GAME_OVER") {
      is_finish = true;
      // compute medals
      std::array<double, NP> v;
      for (int i = 0; i < NP; ++i) {
        v[i] = -(std::pow(cursors[i][0], 2) + std::pow(cursors[i][1], 2));
      }
      auto sorted_array = rank_numbers(v);
      for (int i = 0; i < sorted_array.size(); ++i) {
        medals[i] = sorted_array[i];
      }
      initial_size = 0;

    } else {
      is_finish = false;
      auto in_v = string_to_vector(gpu);
      wind = in_v;
      initial_size = wind.size();
    }
  }

  Archery()
      : initial_size(0), wind({}), cursors({}), medals({4, 4, 4}),
        is_finish(false) {}

  Archery action(const std::array<Action, N> &actions) const;

  bool game_finish() const { return is_finish; }

  void print_game_state() const {
    if (is_finish) {
      std::cerr << "FINISH  " << std::endl;
    }
    std::cerr << "Wind " << std::endl;
    std::copy(wind.begin(), wind.end(),
              std::ostream_iterator<int>(std::cerr, " "));

    for (int i = 0; i < N; ++i) {
      std::cerr << "Joueur " << i + 1 << " : Position ";
      std::copy(cursors[i].begin(), cursors[i].end(),
                std::ostream_iterator<int>(std::cerr, " "));
      cerr << "med" << medals[i];
      cerr << std::endl;
    }
  }

  int initial_size;
  std::array<int, N> get_medals() const { return medals; }
  // private:
  std::vector<int> wind;
  std::array<std::array<int, 2>, N> cursors;
  std::array<int, N> medals; // Medals
  bool is_finish;
};

template <std::size_t N>
Archery<N> Archery<N>::action(const std::array<Action, N> &actions) const {

  if (is_finish) {
    return *this;
  }

  auto new_cursors = cursors;
  auto new_wind = wind;
  for (std::size_t i = 0; i < actions.size(); ++i) {
    Action a = actions[i];

    int offset = new_wind[0];

    int dx = 0;
    int dy = 0;
    if (a == Action::DOWN) {
      dy = offset;
    } else if (a == Action::LEFT) {
      dx = -offset;
    } else if (a == Action::RIGHT) {
      dx = offset;
    } else {
      dy = -offset;
    }
    std::array<int, 2> &cursor = new_cursors[i];
    cursor[0] += dx;
    cursor[1] += dy;
    int max_dist = 20;
    if (cursor[0] > max_dist) {
      cursor[0] = max_dist;
    }
    if (cursor[1] > max_dist) {
      cursor[1] = max_dist;
    }
    if (cursor[0] < -max_dist) {
      cursor[0] = -max_dist;
    }
    if (cursor[1] < -max_dist) {
      cursor[1] = -max_dist;
    }
  }

  new_wind.erase(new_wind.begin());

  auto ret_val = *this;
  ret_val.wind = new_wind;
  ret_val.cursors = new_cursors;

  if (new_wind.empty()) {
    // game finish
    ret_val.is_finish = true;

    std::array<double, NP> v;
    for (int i = 0; i < NP; ++i) {
      v[i] = -(std::pow(new_cursors[i][0], 2) + std::pow(new_cursors[i][1], 2));
    }

    auto sorted_array = rank_numbers(v);
    for (int i = 0; i < sorted_array.size(); ++i) {
      ret_val.medals[i] = sorted_array[i];
    }
  }

  return ret_val;
}

//------------------------------------ DIVING ------

template <std::size_t N> class Diving {
public:
  Diving(const std::string &gpu_str, const std::array<int, N> &points,
         const std::array<int, N> &combos)
      : points(points), combos(combos), medals({4, 4, 4}) {

    if (gpu_str == "GAME_OVER") {
      is_finish = true;
      auto sorted_array = rank_numbers(points);
      for (int i = 0; i < sorted_array.size(); ++i) {
        medals[i] = sorted_array[i];
      }
    } else {
      is_finish = false;
      for (char c : gpu_str) {
        switch (c) {
        case 'U':
          gpu.push_back(Action::UP);
          break;
        case 'L':
          gpu.push_back(Action::LEFT);
          break;
        case 'D':
          gpu.push_back(Action::DOWN);
          break;
        case 'R':
          gpu.push_back(Action::RIGHT);
          break;
        default:
          throw std::invalid_argument("Invalid character in GPU string");
        }
      }
    }
  }

  Diving()
      : gpu({}), points({}), combos({}), medals({4, 4, 4}), is_finish(false) {}

  bool game_finish() const { return is_finish; }
  std::array<int, N> get_medals() const { return medals; }
  Diving action(const std::array<Action, N> &actions) const {
    if (is_finish) {
      return *this;
    }

    Diving new_state = *this;
    std::vector<Action> new_gpu = new_state.gpu;
    std::array<int, N> new_points = new_state.points;
    std::array<int, N> new_combos = new_state.combos;

    for (size_t i = 0; i < actions.size(); ++i) {
      Action action = actions[i];
      if (new_gpu[0] == action) {
        new_combos[i]++;
        new_points[i] += new_combos[i];
      } else {
        new_combos[i] = 0;
      }
    }

    new_gpu.erase(new_gpu.begin());

    new_state.points = new_points;
    new_state.combos = new_combos;
    new_state.gpu = new_gpu;

    if (new_gpu.empty()) {
      auto sorted_array = rank_numbers(new_points);
      for (int i = 0; i < sorted_array.size(); ++i) {
        new_state.medals[i] = sorted_array[i];
      }
      new_state.is_finish = true;
    }

    return new_state;
  }

  void print_game_state() const {
    std::cerr << "GPU: ";
    for (const auto &action : gpu) {
      switch (action) {
      case Action::UP:
        std::cerr << "UP ";
        break;
      case Action::LEFT:
        std::cerr << "LEFT ";
        break;
      case Action::DOWN:
        std::cerr << "DOWN ";
        break;
      case Action::RIGHT:
        std::cerr << "RIGHT ";
        break;
      }
    }
    std::cerr << std::endl;
    std::cerr << "Points: ";
    for (const auto &point : points) {
      std::cerr << point << " ";
    }
    std::cerr << std::endl;
    std::cerr << "Combos: ";
    for (const auto &combo : combos) {
      std::cerr << combo << " ";
    }
    std::cerr << std::endl;
    std::cerr << "Medals: ";
    for (const auto &medal : medals) {
      std::cerr << medal << " ";
    }
    std::cerr << std::endl;
    std::cerr << "Is Finish: " << std::boolalpha << is_finish << std::endl;
  }

  // private:
  std::vector<Action> gpu;
  std::array<int, N> points;
  std::array<int, N> combos;
  std::array<int, N> medals; // Medals
  bool is_finish;
};

//------------------------------------ ROLLER ------

Action char_to_action(char c) {
  switch (c) {
  case 'U':
    return Action::UP;
  case 'L':
    return Action::LEFT;
  case 'D':
    return Action::DOWN;
  case 'R':
    return Action::RIGHT;
  default:
    throw std::invalid_argument("Invalid character for Action");
  }
}

template <std::size_t N> class Roller {
public:
  Roller(const std::string &gpu, const std::array<int, N> &positions,
         const std::array<int, N> &risk, int timer)
      : positions(positions), risk(risk), medals({4, 4, 4}), length(10),
        timer(timer) {

    if (gpu == "GAME_OVER") {
      is_finish = true;
      auto sorted_array = rank_numbers(positions);
      for (int i = 0; i < sorted_array.size(); ++i) {
        medals[i] = sorted_array[i];
      }
    } else {
      is_finish = false;
      if (gpu.size() != 4) {
        throw std::invalid_argument(
            "GPU string must contain exactly 4 characters.");
      }
      for (size_t i = 0; i < 4; ++i) {
        directions[i] = char_to_action(gpu[i]);
      }
    }
  }

  Roller()
      : positions({}), risk({}), medals({4, 4, 4}), length(0), timer(0),
        is_finish(false) {
    directions = {Action::RIGHT, Action::DOWN, Action::LEFT, Action::UP};
  }

  bool game_finish() const { return is_finish; }
  std::array<int, N> get_medals() const { return medals; }
  Roller action(const std::array<Action, N> &actions) const {
    if (is_finish) {
      return *this;
    }

    Roller new_state = *this;
    std::array<int, N> new_positions = new_state.positions;
    std::array<int, N> new_risk = new_state.risk;

    for (size_t i = 0; i < actions.size(); ++i) {
      Action action = actions[i];

      if (new_risk[i] < 0) {
        new_risk[i]++;
        continue;
      }

      int idx = std::find(new_state.directions.begin(),
                          new_state.directions.end(), action) -
                new_state.directions.begin();
      int dx = idx == 0 ? 1 : (idx == 3 ? 3 : 2);

      new_positions[i] += dx;
      int riskValue = -1 + idx;
      new_risk[i] = std::max(0, new_risk[i] + riskValue);
    }

    for (size_t i = 0; i < new_positions.size(); ++i) {
      if (new_risk[i] < 0) {
        continue;
      }

      bool clash = false;
      for (size_t k = 0; k < new_positions.size(); ++k) {
        if (k == i) {
          continue;
        }
        if (new_positions[k] % new_state.length ==
            new_positions[i] % new_state.length) {
          clash = true;
          break;
        }
      }
      if (clash) {
        new_risk[i] += 2;
      }

      if (new_risk[i] >= 5) {
        new_risk[i] = -2; // stun
      }
    }

    new_state.timer--;

    // Mettre à jour l'état
    new_state.positions = new_positions;
    new_state.risk = new_risk;

    static std::random_device rd;
    static std::mt19937 g(rd());
    std::shuffle(new_state.directions.begin(), new_state.directions.end(), g);

    if (new_state.timer == 0) {
      auto sorted_array = rank_numbers(new_state.positions);
      for (int i = 0; i < sorted_array.size(); ++i) {
        new_state.medals[i] = sorted_array[i];
      }
      new_state.is_finish = true;
    }
    return new_state;
  }

  void print_game_state() const {
    std::cerr << "Positions: ";
    for (const auto &pos : positions) {
      std::cerr << pos << " ";
    }
    std::cerr << std::endl;
    std::cerr << "Risk: ";
    for (const auto &r : risk) {
      std::cerr << r << " ";
    }
    std::cerr << std::endl;
    std::cerr << "Timer: " << timer << std::endl;
    std::cerr << "Is Finish: " << std::boolalpha << is_finish << std::endl;
    std::cerr << "Directions: ";
    for (const auto &dir : directions) {
      switch (dir) {
      case Action::UP:
        std::cerr << "UP ";
        break;
      case Action::LEFT:
        std::cerr << "LEFT ";
        break;
      case Action::DOWN:
        std::cerr << "DOWN ";
        break;
      case Action::RIGHT:
        std::cerr << "RIGHT ";
        break;
      }
    }

    std::cerr << "Medals: ";
    for (const auto &medal : medals) {
      std::cerr << medal << " ";
    }
    std::cerr << std::endl;
  }

  // private:
  std::array<int, N> positions;
  std::array<int, N> risk;
  std::array<Action, 4> directions;
  std::array<int, N> medals;
  int length;
  int timer;
  bool is_finish;
};

//*********

template <std::size_t N> class FullGame {
public:
  FullGame(input_tuple cou, input_tuple arr, input_tuple rol, input_tuple div,
           const std::vector<std::vector<int>> &real_s,
           const std::vector<std::array<Action, NP>> &poss_actions)
      : poss_actions(poss_actions) {
    {
      auto [gpu, reg_0, reg_1, reg_2, reg_3, reg_4, reg_5, reg_6] = cou;
      co = Course<N>(gpu, {reg_0, reg_1, reg_2}, {reg_3, reg_4, reg_5});
    }
    {
      auto [gpu, reg_0, reg_1, reg_2, reg_3, reg_4, reg_5, reg_6] = arr;
      std::array<int, 2> a1 = {reg_0, reg_1};
      std::array<int, 2> a2 = {reg_2, reg_3};
      std::array<int, 2> a3 = {reg_4, reg_5};
      ar = Archery<N>(gpu, {a1, a2, a3});
    }
    {
      auto [gpu, reg_0, reg_1, reg_2, reg_3, reg_4, reg_5, reg_6] = rol;
      ro = Roller<N>(gpu, {reg_0, reg_1, reg_2}, {reg_3, reg_4, reg_5}, reg_6);
    }
    {
      auto [gpu, reg_0, reg_1, reg_2, reg_3, reg_4, reg_5, reg_6] = div;
      di = Diving<N>(gpu, {reg_0, reg_1, reg_2}, {reg_3, reg_4, reg_5});
    }

    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < 13; ++j) {
        real_scores[i][j] = real_s[i][j];
      }
    }
  }

  FullGame action(const std::array<Action, N> &actions) const;
  FullGame random_action() const;
  FullGame random_play_until_end() const;
  const std::array<float, N> players_stats() const;

  bool is_all_finish() const {
    return co.game_finish() && ar.game_finish() && ro.game_finish() &&
           di.game_finish();
  }

  void set_all_finish() {
    co.is_finish = true;
    ar.is_finish = true;
    ro.is_finish = true;
    di.is_finish = true;
  }

  void print_game_state() const {
    co.print_game_state();
    ar.print_game_state();
    ro.print_game_state();
    di.print_game_state();
  }
  std::vector<std::array<Action, N>> poss_actions;

private:
  std::array<Action, 3> select_random_combination() const;

  std::array<std::array<int, 13>, N> real_scores;

  Course<N> co;
  Archery<N> ar;
  Roller<N> ro;
  Diving<N> di;
};

template <std::size_t N>
const std::array<float, N> FullGame<N>::players_stats() const {
  std::array<std::array<float, N>, 4> ret_ar = {};
  std::array<float, N> ret_ar_final = {};
  auto mc = co.get_medals();
  auto ma = ar.get_medals();
  auto mr = ro.get_medals();
  auto md = di.get_medals();

  for (int i = 0; i < N; ++i) {
    ret_ar[0][i] = real_scores[i][1] * 3 + real_scores[i][2];
    ret_ar[1][i] = real_scores[i][4] * 3 + real_scores[i][5];
    ret_ar[2][i] = real_scores[i][7] * 3 + real_scores[i][8];
    ret_ar[3][i] = real_scores[i][10] * 3 + real_scores[i][11];
  }

  for (int i = 0; i < N; ++i) {
    if (mc[i] == 1) {
      ret_ar[0][i] += 3;
    } else if (mc[i] == 2) {
      ret_ar[0][i] += 1;
    } else if (mc[i] == 3) {
      ret_ar[0][i] += 0;
    } else {
      cerr << "HEE mc " << mc[i] << endl;
    }

    if (ma[i] == 1) {
      ret_ar[1][i] += 3;
    } else if (ma[i] == 2) {
      ret_ar[1][i] += 1;
    } else if (ma[i] == 3) {
      ret_ar[1][i] += 0;
    } else {
      cerr << "HEE ma " << ma[i] << endl;
    }

    if (mr[i] == 1) {
      ret_ar[2][i] += 2.8;
    } else if (mr[i] == 2) {
      ret_ar[2][i] += 0.8;
    } else if (mr[i] == 3) {
      ret_ar[2][i] += 0.1;
    } else {
      cerr << "HEE mr " << mr[i] << endl;
    }

    if (md[i] == 1) {
      ret_ar[3][i] += 3;
    } else if (md[i] == 2) {
      ret_ar[3][i] += 1;
    } else if (md[i] == 3) {
      ret_ar[3][i] += 0;
    } else {
      cerr << "HEE md " << md[i] << endl;
    }
  }

  for (int j = 0; j < N; ++j) {
    float mult = 1.0;
    for (int i = 0; i < 4; ++i) {
      mult *= ret_ar[i][j];
    }
    ret_ar_final[j] = mult;
  }

  auto ordered_a = rank_numbers(ret_ar_final);

  std::array<float, N> ret_ar_final_final = {};

  for (int i = 0; i < N; ++i) {
    if (ordered_a[i] == 1) {
      ret_ar_final_final[i] = 1;
    } else if (ordered_a[i] == 2) {
      ret_ar_final_final[i] = -1;
    } else if (ordered_a[i] == 3) {
      ret_ar_final_final[i] = -1;
    } else {
      cerr << "HEE mc 22 " << mc[i] << endl;
    }
  }
  return ret_ar_final_final;
}

template <std::size_t N>
FullGame<N> FullGame<N>::random_play_until_end() const {
  FullGame<N> curr_game = *this;

  while (!curr_game.is_all_finish()) {
    curr_game = curr_game.random_action();
  }

  return curr_game;
}

template <std::size_t N>
std::array<Action, 3> FullGame<N>::select_random_combination() const {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, poss_actions.size() - 1);

  int random_index = dis(gen);
  return poss_actions[random_index];
}

template <std::size_t N>
FullGame<N> FullGame<N>::action(const std::array<Action, N> &actions) const {
  auto nco = co.action(actions);
  auto nar = ar.action(actions);
  auto nro = ro.action(actions);
  auto ndi = di.action(actions);

  FullGame fg = *this;
  fg.co = nco;
  fg.ar = nar;
  fg.ro = nro;
  fg.di = ndi;

  return fg;
}

template <std::size_t N> FullGame<N> FullGame<N>::random_action() const {
  auto sel_a = select_random_combination();
  return action(sel_a);
}

// =================== MCTS ==================
template <std::size_t N> struct Node {
  std::array<float, N> rewards;
  std::array<int, N> n_sims;
  FullGame<N> game;
  const int n_moves = 4 * 4 * 4;
  Node *parent;
  std::vector<std::unique_ptr<Node>> children;

  Node(const Node &) = delete;
  Node(const FullGame<N> &game);
  Node(const FullGame<N> &game, Node *parent,
       const std::array<Action, N> &actions);
};

template <std::size_t N>
Node<N>::Node(const FullGame<N> &game)
    : rewards({}), n_sims({}), game(game), parent(nullptr) {
  children.reserve(n_moves);
}

template <std::size_t N>
Node<N>::Node(const FullGame<N> &game, Node *parent,
              const std::array<Action, N> &actions)
    : rewards({}), n_sims({}), game(game.action(actions)), parent(parent) {
  // this->game = game.action(actions);
  children.reserve(n_moves);
}

// MCTS move
template <std::size_t N> Node<N> *select_and_expand(Node<N> *root) {
  Node<N> *n = root;
  while (true) {
    // return node if game terminated
    if (n->game.is_all_finish())
      return n;
    // expand if new child found
    const int k = n->children.size();
    if (k < n->n_moves) {
      n->children.push_back(
          std::make_unique<Node<N>>(n->game, n, n->game.poss_actions[k]));
      return n->children.back().get();
    }
    // select child node using UCB
    n = select_ucb(n);
  }
}

template <std::size_t N> void backpropagate(Node<N> *node) {
  const auto eval_val = node->game.random_play_until_end().players_stats();
  Node<N> *n = node;
  while (n) {
    for (size_t curr_play = 0; curr_play < N; ++curr_play) {
      n->rewards[curr_play] += eval_val[curr_play];
      n->n_sims[curr_play] += 1;
    }
    n = n->parent;
  }
}

inline float ucb1(float reward, float c_nsims, int p_nsims) {
  const float exploitation = 0; // reward / c_nsims;
  const float exploration = std::sqrt(std::log(1.0 + p_nsims) / c_nsims);

  return exploitation + exploration;
}

template <std::size_t N> Node<N> *select_ucb(const Node<N> *n) {
  std::array<Action, N> ac_play;
  for (int play_id = 0; play_id < N; ++play_id) {
    int best_i = -1;
    float best_score = -std::numeric_limits<float>::infinity();
    for (int i = 0; i < n->n_moves; ++i) {
      const auto &c = n->children[i];
      const float s =
          ucb1(c->rewards[play_id], c->n_sims[play_id], n->n_sims[play_id]);
      if (s > best_score) {
        best_score = s;
        best_i = i;
      }
    }
    ac_play[play_id] = n->game.poss_actions[best_i][play_id];
  }
  int found_i = std::distance(n->game.poss_actions.begin(),
                              std::find(n->game.poss_actions.begin(),
                                        n->game.poss_actions.end(), ac_play));

  return n->children[found_i].get();
}

template <std::size_t N> Action best_node(const Node<N> &root, int curr_play) {
  cerr << "here0 " << endl;
  const auto &cs = root.children;

  std::map<Action, float> m{{Action::DOWN, 0},
                            {Action::LEFT, 0},
                            {Action::RIGHT, 0},
                            {Action::UP, 0}};

  for (int i = 0; i < cs.size(); ++i) {
    m[root.game.poss_actions[i][curr_play]] +=
        cs[i]->rewards[curr_play] / cs[i]->n_sims[curr_play];
  }

  auto x = *std::max_element(
      m.begin(), m.end(),
      [](const pair<Action, float> &p1, const pair<Action, float> &p2) {
        return p1.second < p2.second;
      });

  return x.first;
}

template <std::size_t N>
Action get_best_move(const FullGame<N> &game, int player_idx) {
  Node root(game);
  auto duration = std::chrono::milliseconds(42);

  auto start = std::chrono::steady_clock::now();

  int nb_it = 0;
  while (true) {
    Node<N> *node = select_and_expand(&root);
    backpropagate(node);

    ++nb_it;
    auto now = std::chrono::steady_clock::now();
    auto elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
    if (elapsed >= duration) {
      cerr << "Break at " << elapsed.count() << " " << nb_it << endl;
      break;
    }
  }
  return best_node(root, player_idx);
}

int main() {
  const std::vector<std::array<Action, NP>> all_actions_players =
      generate_combinations();

  int player_idx;
  cin >> player_idx;
  cin.ignore();
  int nb_games;
  cin >> nb_games;
  cin.ignore();

  int step = 0;

  while (1) {
    std::vector<std::vector<int>> scores;
    for (int i = 0; i < 3; i++) {
      string score_info;
      getline(cin, score_info);

      std::vector<int> numbers;
      std::istringstream iss(score_info);
      int number;

      while (iss >> number) {
        numbers.push_back(number);
      }
      assert(numbers.size() == 13);
      scores.push_back(numbers);
    }

    cerr << scores[1][0] << endl;

    vector<input_tuple> vit;
    for (int i = 0; i < nb_games; i++) {
      string gpu;
      int reg_0;
      int reg_1;
      int reg_2;
      int reg_3;
      int reg_4;
      int reg_5;
      int reg_6;
      cin >> gpu >> reg_0 >> reg_1 >> reg_2 >> reg_3 >> reg_4 >> reg_5 >> reg_6;
      cin.ignore();

      cerr << gpu << " " << reg_0 << " " << reg_1 << " " << reg_2 << " "
           << reg_3 << " " << reg_4 << " " << reg_5 << " " << reg_6 << endl;
      vit.push_back(
          make_tuple(gpu, reg_0, reg_1, reg_2, reg_3, reg_4, reg_5, reg_6));
    }
    FullGame<NP> fg(vit[0], vit[1], vit[2], vit[3], scores,
                    all_actions_players);
    auto bm = get_best_move<NP>(fg, player_idx);
    cout << bm << endl;

    ++step;
  }
}
