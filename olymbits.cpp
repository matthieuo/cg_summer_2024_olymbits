#include <iostream>
#include <array>
#include <cassert>
#include <ostream>
#include <algorithm>
#include <vector>
#include <iterator>
#include <cmath>


using namespace std;

enum class Action { UP, LEFT, DOWN, RIGHT };

std::string to_string(Action action) {
    switch (action) {
        case Action::UP: return "UP";
        case Action::LEFT: return "LEFT";
        case Action::DOWN: return "DOWN";
        case Action::RIGHT: return "RIGHT";
        default: return "Unknown";
    }
}

// Surcharge de l'opérateur << pour imprimer les valeurs de l'enum class
std::ostream& operator<<(std::ostream& os, Action action) {
    return os << to_string(action);
}


//enum Medal { GOLD = 1, SILVER, BRONZE, NONE};

constexpr std::size_t NP = 3;


template <std::size_t N, typename T>
std::array<int, N> rank_numbers(const std::array<T, N>& numbers) {
    // Create an array of pairs to store numbers and their original indices
    std::array<std::pair<T, int>, N> number_indices;
    for (size_t i = 0; i < N; ++i) {
        number_indices[i] = std::make_pair(numbers[i], i);
    }

    // Sort the pairs based on the values of the numbers in descending order
    std::sort(number_indices.begin(), number_indices.end(), std::greater<>());

    // Initialize the ranks array
    std::array<int, N> ranks;

    // Assign ranks based on the sorted numbers
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


template <std::size_t N>
class Course {
public:
  // Constructeur initialisant l'état du jeu
  Course(const std::string& gpu, const std::array<int, N>& pos, const std::array<int, N>& diz)
    : gpu(gpu), pos(pos), diz(diz), medals({4,4,4}), is_finish(false) {}

  Course()
    : gpu(""), pos({0,0,0}), diz({0,0,0}), medals({4,4,4}), is_finish(false) {}

  // Fonction pour effectuer une action et mettre à jour l'état du jeu
  Course action(const std::array<Action, 3>& actions) const {

    if (is_finish) {
      return *this;
    }
    
    // Mise à jour de l'état du jeu pour chaque joueur
    std::array<int, N> new_pos = pos;
    std::array<int, N> new_diz = diz;

    bool is_game_finish = false;
    
    for (int i = 0; i < 3; ++i) {
      auto [x,y] = update_agent_pos_diz(pos[i], diz[i], actions[i]);
      new_pos[i] = x;
      new_diz[i] = y;
      if (x == gpu.size() - 1) {
	is_game_finish = true;
      }
      //std::cout << x << "," << y << std::endl;
    }

    auto c =  Course(gpu, new_pos, new_diz);
    //std::cout << new_pos << std::endl;
    if (is_game_finish) {
      auto sorted_array = rank_numbers(new_pos);
      for (int i = 0; i< sorted_array.size(); ++i) {
	c.medals[i] = sorted_array[i];
	c.is_finish = true;
	//std::cout << c.medals[i] << std::endl;
      }
    }

    
    return c;
  }


  // Fonction pour afficher l'état du jeu
  void print_game_state() const {
    if (is_finish){
      std::cerr << "FINISH  "<< std::endl;
    }
    std::cerr << "Piste : " << gpu << std::endl;
    for (int i = 0; i < N; ++i) {
      std::cerr << "Joueur " << i + 1 << " : Position " << pos[i] << ", Étourdissement " << diz[i] << "M" << medals[i] << std::endl;
    }
  }

private:
  // Fonction pour mettre à jour la position d'un agent apres un move
  std::tuple<int,int> update_agent_pos_diz(int init_pos, int init_diz, Action action) const {

    if (init_diz > 0) { // diz -> not move
      return {init_pos, init_diz - 1};
    }
    int target_pos = init_pos;
    assert(init_diz == 0);
    int target_diz = 0;
    // Déterminer le mouvement en fonction de l'action

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
    case  Action::RIGHT:
      target_pos += 3;
      break;
    }

    target_pos = std::min(target_pos, (int)gpu.size()-1);

    //std::cout << init_pos << " " << target_pos << std::endl;

    for (int i = init_pos+1; i<=target_pos; ++i) {
      if ((target_pos - i == 1) && action == Action::UP) {
	//don't check colision at the first if we UP
	
	continue;
      }
      if (gpu[i] == '#') {
	//colision
	target_diz = 2;
	target_pos = i;
	break;
      }
    }

    //std::cout << "*" << init_pos << " " << target_pos << std::endl;
    return {target_pos, target_diz};
  }


  // Variables représentant l'état du jeu
  std::string gpu;
  std::array<int, N> pos; // Positions des agents
  std::array<int, N> diz; // Décomptes d'étourdissement des agents


  std::array<int, N> medals; // Medals
  bool is_finish;
};

//------------------------------------ ARCHERY ------
//------------------------------------ ARCHERY ------
//------------------------------------ ARCHERY ------

template <std::size_t N>
class Archery {
public:
    Archery(const std::vector<int> gpu, const std::array<std::array<int, 2>, N> cursors)
      : wind(gpu), cursors(cursors), is_finish(false) ,medals({4,4,4}){}

  Archery()
    : wind({}), cursors({}), is_finish(false),medals({4,4,4}) {}
  
    Archery action(const std::array<Action, N>& actions) const;

  void print_game_state() const {
    if (is_finish){
      std::cerr << "FINISH  "<< std::endl;
    }
    std::cerr << "Wind " << std::endl;
    std::copy(wind.begin(), wind.end(), std::ostream_iterator<int>(std::cerr, " "));
    
    for (int i = 0; i < N; ++i) {
      std::cerr << "Joueur " << i + 1 << " : Position ";
      std::copy(cursors[i].begin(), cursors[i].end(), std::ostream_iterator<int>(std::cerr, " "));
      cerr << "med" << medals[i];
      cerr  << std::endl;
    }
  }


private:
  std::vector<int> wind; // Supposition que wind est un vector d'entiers
  std::array<std::array<int, 2>, N> cursors;
  std::array<int, N> medals; // Medals
  bool is_finish;

};

template <std::size_t N>
Archery<N> Archery<N>::action(const std::array<Action, N>& actions) const {

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
    std::array<int, 2>& cursor = new_cursors[i];
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

  auto ret_val = Archery<N>(new_wind, new_cursors);
  
  if (new_wind.empty()) {
    //game finish
    ret_val.is_finish = true;
  
    std::array<double, NP> v;
    for (int i=0;i<NP;++i) {
       v[i] = -(std::pow(new_cursors[i][0],2) + std::pow(new_cursors[i][1],2));
  
    }
  
    auto sorted_array = rank_numbers(v);
    for (int i=0; i< sorted_array.size(); ++i) {
	ret_val.medals[i] = sorted_array[i];
 
	//std::cout << c.medals[i] << std::endl;
      }
  }
  
  return ret_val;
}



// ***************************************************
std::vector<int> string_to_vector(const std::string& str) {
    std::vector<int> result;

    for (char ch : str) {
        if (std::isdigit(ch)) {
            result.push_back(ch - '0');
        }
    }
    
    return result;
}


int main()
{
    int player_idx;
    cin >> player_idx; cin.ignore();
    int nb_games;
    cin >> nb_games; cin.ignore();

    // game loop

    int step = 0;
    Course<NP> c;
    Archery<NP> ar;
    
    while (1) {
        for (int i = 0; i < 3; i++) {
            string score_info;
            getline(cin, score_info);
        }


	
        for (int i = 0; i < nb_games; i++) {
            string gpu;
            int reg_0;
            int reg_1;
            int reg_2;
            int reg_3;
            int reg_4;
            int reg_5;
            int reg_6;
            cin >> gpu >> reg_0 >> reg_1 >> reg_2 >> reg_3 >> reg_4 >> reg_5 >> reg_6; cin.ignore();
	    if (step == 0 && i == 0) {
	      c = Course<NP>(gpu, {reg_0,reg_1,reg_2},{reg_3,reg_4,reg_5});
	    }
	    if (step == 0 && i == 1) {
	      std::array<int, 2> a1 = {reg_0,reg_1};
	      std::array<int, 2> a2 = {reg_2,reg_3};
	      std::array<int, 2> a3 = {reg_4,reg_5};
	      auto in_v=string_to_vector(gpu);
	      //std::array<std::array<int, 2>, NP> array = { a1, a2,a3 };
	      ar = Archery<NP>(in_v,{ a1, a2,a3 } );
	    }
	    if (i == 1) {
	      std::cerr << gpu << " " << reg_0 << " " << reg_1 << endl;
	    }
        }
	//c.print_game_state();
	//c = c.action({Action::UP,Action::LEFT,Action::LEFT});
	ar.print_game_state();
	ar = ar.action({Action::UP,Action::LEFT,Action::LEFT});
        // Write an action using cout. DON'T FORGET THE "<< endl"
        // To debug: cerr << "Debug messages..." << endl;

        cout << "UP" << endl;

	++step;
    }
}




