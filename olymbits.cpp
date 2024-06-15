#include <iostream>
#include <array>
#include <cassert>
#include <ostream>
#include <algorithm>
#include <string>
#include <tuple>
#include <vector>
#include <iterator>
#include <cmath>
#include <stdexcept>
#include <optional>
#include <random>


using namespace std;
typedef tuple<std::string, int, int, int, int, int, int, int> input_tuple;

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

/*template <std::size_t N>
class MiniGame {
public:
  virtual MiniGame action(const std::array<Action, N>& actions) const = 0;
  virtual void print_game_state() const =0;

  bool game_finish() const {return is_finish;}
  
private:

  bool is_finish;
  };*/



//-------------- COURSE ----------------
//-------------- COURSE ----------------
//-------------- COURSE ----------------

template <std::size_t N>
class Course {
public:
  // Constructeur initialisant l'état du jeu
  Course(const std::string& gpu, const std::array<int, N>& pos, const std::array<int, N>& diz)
    : gpu(gpu), pos(pos), diz(diz), medals({4,4,4}) {
    is_finish = (gpu == "GAME_OVER");
    if (is_finish) {
      auto sorted_array = rank_numbers(pos);
      for (int i = 0; i< sorted_array.size(); ++i) {
	medals[i] = sorted_array[i];
      }
    }
    }
    

  Course()
    : gpu(""), pos({0,0,0}), diz({0,0,0}), medals({4,4,4}), is_finish(false) {}

  // Fonction pour effectuer une action et mettre à jour l'état du jeu
   Course action(const std::array<Action, N>& actions) const {

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
    Archery(const std::string gpu, const std::array<std::array<int, 2>, N> cursors)
      : cursors(cursors) ,medals({4,4,4}){

      if (gpu == "GAME_OVER") {
	is_finish = true;
	//compute medals
	std::array<double, NP> v;
	for (int i=0;i<NP;++i) {
	  v[i] = -(std::pow(cursors[i][0],2) + std::pow(cursors[i][1],2));
  
	}
	auto sorted_array = rank_numbers(v);
	for (int i=0; i< sorted_array.size(); ++i) {
	  medals[i] = sorted_array[i];
      }
	
      } else {
	is_finish = false;
	auto in_v=string_to_vector(gpu);
	wind = in_v;
      }

    }

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
  std::vector<int> wind; 
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


//------------------------------------ DIVING ------
//------------------------------------ DIVING ------
//------------------------------------ DIVING ------

template <std::size_t N>
class Diving {
public:

    Diving(const std::string& gpu_str, const std::array<int, N>& points, const std::array<int, N>& combos)
        : points(points), combos(combos), medals({4, 4, 4}) {


      if (gpu_str == "GAME_OVER") {
	is_finish = true;
	auto sorted_array = rank_numbers(points);
	for (int i = 0; i< sorted_array.size(); ++i) {
	  medals[i] = sorted_array[i];
	}
      } else {
	is_finish = false;
        // Conversion de la chaîne de caractères en std::vector<Action>
        for (char c : gpu_str) {
            switch (c) {
                case 'U': gpu.push_back(Action::UP); break;
                case 'L': gpu.push_back(Action::LEFT); break;
                case 'D': gpu.push_back(Action::DOWN); break;
                case 'R': gpu.push_back(Action::RIGHT); break;
                default: throw std::invalid_argument("Invalid character in GPU string");
            }
        }
      }
    }

    Diving()
        : gpu({}), points({}), combos({}), medals({4, 4, 4}), is_finish(false) {}

    // Fonction pour effectuer une action et mettre à jour l'état du jeu
    Diving action(const std::array<Action, N>& actions) const {
        if (is_finish) {
            return *this;
        }

        // Créer une copie de l'état actuel
        Diving new_state = *this;
        std::vector<Action> new_gpu = new_state.gpu;
        std::array<int, N> new_points = new_state.points;
        std::array<int, N> new_combos = new_state.combos;

        // Effectuer les actions
        for (size_t i = 0; i < actions.size(); ++i) {
            Action action = actions[i];
            if (!new_gpu.empty() && new_gpu[0] == action) {
                new_combos[i]++;
                new_points[i] += new_combos[i];
            } else {
                new_combos[i] = 0;
            }
        }

	new_gpu.erase(new_gpu.begin()); // Suppression du premier élément de gpu
        

        new_state.points = new_points;
        new_state.combos = new_combos;
        new_state.gpu = new_gpu;

        if (new_gpu.empty()) {
	  auto sorted_array = rank_numbers(new_points);
	  for (int i = 0; i< sorted_array.size(); ++i) {
	    new_state.medals[i] = sorted_array[i];
	  }
	  new_state.is_finish = true;
        }

        return new_state;
    }

    // Fonction pour afficher l'état du jeu
    void print_game_state() const {
        std::cerr << "GPU: ";
        for (const auto& action : gpu) {
            switch (action) {
                case Action::UP: std::cerr << "UP "; break;
                case Action::LEFT: std::cerr << "LEFT "; break;
                case Action::DOWN: std::cerr << "DOWN "; break;
                case Action::RIGHT: std::cerr << "RIGHT "; break;
            }
        }
        std::cerr << std::endl;
        std::cerr << "Points: ";
        for (const auto& point : points) {
            std::cerr << point << " ";
        }
        std::cerr << std::endl;
        std::cerr << "Combos: ";
        for (const auto& combo : combos) {
            std::cerr << combo << " ";
        }
        std::cerr << std::endl;
        std::cerr << "Medals: ";
        for (const auto& medal : medals) {
            std::cerr << medal << " ";
        }
        std::cerr << std::endl;
        std::cerr << "Is Finish: " << std::boolalpha << is_finish << std::endl;
    }

private:
    // Variables représentant l'état du jeu
    std::vector<Action> gpu;
    std::array<int, N> points;
    std::array<int, N> combos;
    std::array<int, N> medals; // Medals
    bool is_finish;
};


//------------------------------------ ROLLER ------
//------------------------------------ ROLLER ------
//------------------------------------ ROLLER ------



Action char_to_action(char c) {
    switch (c) {
        case 'U': return Action::UP;
        case 'L': return Action::LEFT;
        case 'D': return Action::DOWN;
        case 'R': return Action::RIGHT;
        default: throw std::invalid_argument("Invalid character for Action");
    }
}

template <std::size_t N>
class Roller {
public:
    // Constructeur initialisant l'état du jeu
  Roller( const std::string& gpu,const std::array<int, N>& positions, const std::array<int, N>& risk, int timer)
    : positions(positions), risk(risk), medals({4,4,4}),length(10), timer(timer) {


    if (gpu == "GAME_OVER") {
      is_finish = true;
      auto sorted_array = rank_numbers(positions);
      for (int i = 0; i< sorted_array.size(); ++i) {
	medals[i] = sorted_array[i];
      }
    } else {
      is_finish = false;
      if (gpu.size() != 4) {
	throw std::invalid_argument("GPU string must contain exactly 4 characters.");
      }
      for (size_t i = 0; i < 4; ++i) {
	directions[i] = char_to_action(gpu[i]);
      }
    }
  }

    Roller()
        : positions({}), risk({}), medals({4,4,4}), length(0), timer(0), is_finish(false) {
        directions = {Action::RIGHT, Action::DOWN, Action::LEFT, Action::UP};
    }

    // Fonction pour effectuer une action et mettre à jour l'état du jeu
    Roller action(const std::array<Action, N>& actions) const {
        if (is_finish) {
            return *this;
        }

        // Créer une copie de l'état actuel
        Roller new_state = *this;
        std::array<int, N> new_positions = new_state.positions;
        std::array<int, N> new_risk = new_state.risk;

        // Effectuer les actions
        for (size_t i = 0; i < actions.size(); ++i) {
            Action action = actions[i];

            if (new_risk[i] < 0) {
                new_risk[i]++;
                continue;
            }

            int idx = std::find(new_state.directions.begin(), new_state.directions.end(), action) - new_state.directions.begin();
            int dx = idx == 0 ? 1 : (idx == 3 ? 3 : 2);

            new_positions[i] += dx;
            int riskValue = -1 + idx;
            new_risk[i] = std::max(0, new_risk[i] + riskValue);
        }

        // Vérifier les collisions et ajuster le risque
        for (size_t i = 0; i < new_positions.size(); ++i) {
            if (new_risk[i] < 0) {
                continue;
            }

            bool clash = false;
            for (size_t k = 0; k < new_positions.size(); ++k) {
                if (k == i) {
                    continue;
                }
                if (new_positions[k] % new_state.length == new_positions[i] % new_state.length) {
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

        // Décrémenter le timer
        new_state.timer--;

        // Mettre à jour l'état
        new_state.positions = new_positions;
        new_state.risk = new_risk;

        if (new_state.timer == 0) {
	  auto sorted_array = rank_numbers(new_state.positions);
	  for (int i = 0; i< sorted_array.size(); ++i) {
	    new_state.medals[i] = sorted_array[i];
	  }
            new_state.is_finish = true;
        }

        return new_state;
    }

    // Fonction pour mélanger les directions
    void shuffle_directions(const std::optional<std::string>& new_gpu = std::nullopt) {
        if (!new_gpu) {
            // Mélanger aléatoirement
            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(directions.begin(), directions.end(), g);
        } else {
            // Mélanger en fonction de new_gpu
            if (new_gpu->size() != 4) {
                throw std::invalid_argument("GPU string must contain exactly 4 characters.");
            }
            for (size_t i = 0; i < 4; ++i) {
                directions[i] = char_to_action(new_gpu->at(i));
            }
        }
    }

    // Fonction pour afficher l'état du jeu
    void print_game_state() const {
        std::cerr << "Positions: ";
        for (const auto& pos : positions) {
            std::cerr << pos << " ";
        }
        std::cerr << std::endl;
        std::cerr << "Risk: ";
        for (const auto& r : risk) {
            std::cerr << r << " ";
        }
        std::cerr << std::endl;
        std::cerr << "Timer: " << timer << std::endl;
        std::cerr << "Is Finish: " << std::boolalpha << is_finish << std::endl;
        std::cerr << "Directions: ";
        for (const auto& dir : directions) {
            switch (dir) {
                case Action::UP: std::cerr << "UP "; break;
                case Action::LEFT: std::cerr << "LEFT "; break;
                case Action::DOWN: std::cerr << "DOWN "; break;
                case Action::RIGHT: std::cerr << "RIGHT "; break;
            }
        }

	std::cerr << "Medals: ";
        for (const auto& medal : medals) {
	  std::cerr << medal << " ";
        }
        std::cerr << std::endl;
    }

private:
  // Variables représentant l'état du jeu
  std::array<int, N> positions;
  std::array<int, N> risk;
  std::array<Action, 4> directions;
  std::array<int, N> medals; // Medals
  int length;
  int timer;
  bool is_finish;
};







//*********


    
template <std::size_t N>
class FullGame {
public:

  
  
  FullGame(input_tuple cou, input_tuple arr, input_tuple rol, input_tuple div ) {
    {
      auto [gpu, reg_0,reg_1,reg_2,reg_3,reg_4,reg_5,reg_6] = cou;
      co = Course<N>(gpu, {reg_0,reg_1,reg_2},{reg_3,reg_4,reg_5});
    }
    {
      auto [gpu, reg_0,reg_1,reg_2,reg_3,reg_4,reg_5,reg_6] = arr;
      std::array<int, 2> a1 = {reg_0,reg_1};
      std::array<int, 2> a2 = {reg_2,reg_3};
      std::array<int, 2> a3 = {reg_4,reg_5};
      ar = Archery<N>(gpu,{ a1, a2,a3 } );
    }
    {
      auto [gpu, reg_0,reg_1,reg_2,reg_3,reg_4,reg_5, reg_6] = rol;
      ro = Roller<N>(gpu, {reg_0,reg_1,reg_2},{reg_3,reg_4,reg_5},reg_6);
    }
    {
      auto [gpu, reg_0,reg_1,reg_2,reg_3,reg_4,reg_5,reg_6] = div;
      di = Diving<N>(gpu, {reg_0,reg_1,reg_2},{reg_3,reg_4,reg_5});
    }
  }
  FullGame action(const std::array<Action, N>& actions) const;

  void print_game_state() const {
    co.print_game_state();
    ar.print_game_state();
    ro.print_game_state();
    di.print_game_state();
  }
private:
  Course<N> co;
  Archery<N> ar;
  Roller<N> ro;
  Diving<N> di;
  
};

template <std::size_t N>
FullGame<N> FullGame<N>::action(const std::array<Action, N>& actions) const {
  auto nco = co.action(actions);
  auto nar = ar.action(actions);
  auto nro = ro.action(actions);
  auto ndi = di.action(actions);

  FullGame fg;
  fg.co = nco;
  fg.ar = nar;
  fg.ro = nro;
  fg.di = ndi;

  return fg;
  
}
  
int main()
{
    int player_idx;
    cin >> player_idx; cin.ignore();
    int nb_games;
    cin >> nb_games; cin.ignore();

    // game loop

    int step = 0;
    /*Course<NP> c;
    Archery<NP> ar;
    Diving<NP> div;
    Roller<NP> rol;*/
    
    while (1) {
        for (int i = 0; i < 3; i++) {
            string score_info;
            getline(cin, score_info);
        }


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
            cin >> gpu >> reg_0 >> reg_1 >> reg_2 >> reg_3 >> reg_4 >> reg_5 >> reg_6; cin.ignore();
	    vit.push_back(make_tuple(gpu,reg_0,reg_1,reg_2,reg_3,reg_4,reg_5,reg_6));
	    /*    if (step == 0 && i == 0) {
	      c = Course<NP>(gpu, {reg_0,reg_1,reg_2},{reg_3,reg_4,reg_5});
	    }
	    if (step == 0 && i == 1) {
	      std::array<int, 2> a1 = {reg_0,reg_1};
	      std::array<int, 2> a2 = {reg_2,reg_3};
	      std::array<int, 2> a3 = {reg_4,reg_5};
	      //auto in_v=string_to_vector(gpu);
	      //std::array<std::array<int, 2>, NP> array = { a1, a2,a3 };
	      ar = Archery<NP>(gpu,{ a1, a2,a3 } );
	    }

	    if (step == 0 && i == 3) {
	      div = Diving<NP>(gpu, {reg_0,reg_1,reg_2},{reg_3,reg_4,reg_5});
	    }
	    if (step == 0 && i == 2) {
	      rol = Roller<NP>(gpu, {reg_0,reg_1,reg_2},{reg_3,reg_4,reg_5},reg_6);
	    }
	    if (i == 2) {
	      std::cerr << gpu << " " << reg_0 << " " << reg_3 << " " << reg_6 << endl;
	      if (gpu != "GAME_OVER") {
	      rol.shuffle_directions(gpu);
	      }
	      }*/
        }
	FullGame<NP> fg(vit[0],vit[1],vit[2],vit[3]);
	fg.print_game_state();
	//c.print_game_state();
	//c = c.action({Action::UP,Action::LEFT,Action::LEFT});
	//ar.print_game_state();
	//div.print_game_state();
	//div = div.action({Action::UP,Action::LEFT,Action::LEFT});
	//rol.print_game_state();
	
	//rol = rol.action({Action::UP,Action::UP,Action::UP});

	
        // Write an action using cout. DON'T FORGET THE "<< endl"
        // To debug: cerr << "Debug messages..." << endl;

        cout << "UP" << endl;

	++step;
    }
}




